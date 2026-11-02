#!/usr/bin/env bash
#SBATCH --job-name=FiberTools_MergeHeatmaps
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=85G
#SBATCH --time=24:00:00
#SBATCH --output=FiberTools_MergeHeatmaps.%j.out
#SBATCH --error=FiberTools_MergeHeatmaps.%j.err

# =============================================================================
# Merge nucs.bam files across replicates + both runs, then generate
# fractional TSS heatmaps and multi-sample comparison plots.
#
# Bedgraph column layout from ft pileup --m6a --cpg:
#   col1=chrom col2=start col3=end col4=coverage col5=fire_coverage
#   col6=score col7=nuc_count col8=msp_count col9=m6a_count col10=cpg_count
#
# Signal is RAW COUNTS — divided by col4 (total coverage) to get 0-1 fraction.
# NOTE: set -e is intentionally NOT used so one group failure doesn't kill all.
# =============================================================================

set -uo pipefail

module load Miniforge3/24.11.3-0 ucsc/443
source activate /home/ry00555/fibertools

FT_RESULTS="/lustre2/scratch/ry00555/ONTRun9_10Combined/fibertools_results"
OUT_DIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/MergedForHeatmaps"
GENOME="/home/ry00555/Research/Genomes/GenBankNcrassachromsizes.txt"
TSS_BED="/home/ry00555/Research/Genomes/neurospora_TSS_stranded.bed"
GROUPS_FILE="/home/ry00555/Research/FungiCutAndRun/fiber_groups.txt"

mkdir -p "$OUT_DIR"

if [ ! -f "$GROUPS_FILE" ]; then
    echo "❌  Groups file not found: $GROUPS_FILE"; exit 1
fi
if [ ! -f "$TSS_BED" ]; then
    echo "❌  TSS BED not found: $TSS_BED"; exit 1
fi

# =============================================================================
# Helper: extract fractional signal (signal_col / coverage) → bedgraph + BigWig
# =============================================================================
make_fraction_tracks() {
    local BG="$1"
    local BG4="$2"
    local BW="$3"
    local COL="$4"

    if [ ! -f "$BG4" ]; then
        echo "      Extracting fraction (col$COL / col4) → $(basename $BG4)"
        awk -v col="$COL" '
            NR>1 && $4 > 0 {
                frac = $col / $4
                print $1"\t"$2"\t"$3"\t"frac
            }' "$BG" \
            | sort -k1,1 -k2,2n \
            | awk -v genome="$GENOME" '
                BEGIN { while ((getline line < genome) > 0) {
                    split(line,a,"\t"); chromlen[a[1]]=a[2] } }
                { if ($2<0) $2=0
                  if (chromlen[$1] && $3>chromlen[$1]) $3=chromlen[$1]
                  if ($2<$3) print }
            ' > "$BG4" || { echo "      ❌  awk failed"; return 1; }

        if [ ! -s "$BG4" ]; then
            echo "      ❌  Empty output"; rm -f "$BG4"; return 1
        fi
        echo "      ✅  $(basename $BG4)"
    else
        echo "      ⏭️  $(basename $BG4) exists"
    fi

    if [ ! -f "$BW" ]; then
        echo "      Converting → $(basename $BW)"
        bedGraphToBigWig "$BG4" "$GENOME" "$BW" \
            && echo "      ✅  $(basename $BW)" \
            || echo "      ❌  bedGraphToBigWig failed (non-fatal)"
    else
        echo "      ⏭️  $(basename $BW) exists"
    fi
    return 0
}

# =============================================================================
# Main: process all groups — merge, pileup, fractional tracks
# =============================================================================
echo "============================================================"
echo " Merging replicates + runs for heatmap generation"
echo " Groups file : $GROUPS_FILE"
echo " TSS BED     : $TSS_BED"
echo " Output      : $OUT_DIR"
echo "============================================================"

PROCESSED_GROUPS=""

while read -r LINE; do
    [[ -z "$LINE" || "$LINE" == \#* ]] && continue

    GROUP=$(echo "$LINE" | awk '{print $1}')
    BARCODES=$(echo "$LINE" | awk '{$1=""; print $0}' | xargs)

    echo ""
    echo "╔══════════════════════════════════════════════════════╗"
    echo "  Group: $GROUP  (barcodes: $BARCODES)"
    echo "╚══════════════════════════════════════════════════════╝"

    GROUP_DIR="$OUT_DIR/${GROUP}"
    MERGED_BAM="$OUT_DIR/${GROUP}_merged.nucs.bam"
    mkdir -p "$GROUP_DIR"

    # ── Collect nucs.bam files ───────────────────────────────────
    BAMS_TO_MERGE=""
    for RUN in ONTRun9 ONTRun10; do
        for BARCODE in $BARCODES; do
            for SAMPLE_DIR in "$FT_RESULTS/$RUN"/*${BARCODE}*/; do
                [ -d "$SAMPLE_DIR" ] || continue
                NUCS_BAM=$(find "$SAMPLE_DIR" -maxdepth 1 -name "*.nucs.bam" | head -1)
                if [ -n "$NUCS_BAM" ] && [ -f "$NUCS_BAM" ]; then
                    echo "  Found: $(basename $NUCS_BAM)"
                    BAMS_TO_MERGE="$BAMS_TO_MERGE $NUCS_BAM"
                else
                    echo "  ⚠️  No nucs.bam in $SAMPLE_DIR"
                fi
            done
        done
    done

    BAMS_TO_MERGE="${BAMS_TO_MERGE# }"
    if [ -z "$BAMS_TO_MERGE" ]; then
        echo "  ❌  No nucs.bam files found — skipping"
        continue
    fi

    echo "  Total: $(echo $BAMS_TO_MERGE | wc -w) nucs.bam(s)"

    # ── Merge ────────────────────────────────────────────────────
    if [ ! -f "$MERGED_BAM" ]; then
        BAM_COUNT=$(echo "$BAMS_TO_MERGE" | wc -w)
        if [ "$BAM_COUNT" -eq 1 ]; then
            cp $BAMS_TO_MERGE "$MERGED_BAM" && samtools index "$MERGED_BAM"
        else
            samtools merge -f -@ 8 "$MERGED_BAM" $BAMS_TO_MERGE && samtools index "$MERGED_BAM"
        fi
        echo "  ✅  Merged: $(basename $MERGED_BAM)"
    else
        echo "  ⏭️  Merged BAM already exists"
    fi

    # ── Single pileup ────────────────────────────────────────────
    PILEUP_BG="$GROUP_DIR/${GROUP}.allpileup.bedgraph"
    if [ ! -f "$PILEUP_BG" ]; then
        echo "  Running pileup..."
        ft pileup --m6a --cpg --per-base \
            --out "$PILEUP_BG" "$MERGED_BAM" \
            || { echo "  ❌ pileup failed — skipping group"; continue; }
        echo "  ✅  Pileup done"
    else
        echo "  ⏭️  Pileup exists"
    fi

    echo "  Columns: $(head -1 $PILEUP_BG)"

    # ── Fractional bedgraphs + BigWigs ───────────────────────────
    echo "  Making fractional tracks..."
    make_fraction_tracks "$PILEUP_BG" "$GROUP_DIR/${GROUP}.nuc.frac.bedgraph" "$GROUP_DIR/${GROUP}.nuc.bw" 7
    make_fraction_tracks "$PILEUP_BG" "$GROUP_DIR/${GROUP}.msp.frac.bedgraph" "$GROUP_DIR/${GROUP}.msp.bw" 8
    make_fraction_tracks "$PILEUP_BG" "$GROUP_DIR/${GROUP}.m6a.frac.bedgraph" "$GROUP_DIR/${GROUP}.m6A.bw" 9
    make_fraction_tracks "$PILEUP_BG" "$GROUP_DIR/${GROUP}.5mC.frac.bedgraph" "$GROUP_DIR/${GROUP}.5mC.bw" 10

    echo "  Sample m6A fractions (first 3 non-zero):"
    awk '$4 > 0' "$GROUP_DIR/${GROUP}.m6a.frac.bedgraph" 2>/dev/null | head -3 | sed 's/^/    /' || true

    PROCESSED_GROUPS="$PROCESSED_GROUPS $GROUP"
    echo "  ✅  $GROUP complete"

done < "$GROUPS_FILE"

echo ""
echo "============================================================"
echo " All groups processed: $PROCESSED_GROUPS"
echo " Generating heatmaps..."
echo "============================================================"

# =============================================================================
# Per-group heatmaps
# =============================================================================
for GROUP in $PROCESSED_GROUPS; do
    GROUP_DIR="$OUT_DIR/${GROUP}"
    M6A_BW="$GROUP_DIR/${GROUP}.m6A.bw"
    MSP_BW="$GROUP_DIR/${GROUP}.msp.bw"
    NUC_BW="$GROUP_DIR/${GROUP}.nuc.bw"
    CPG_BW="$GROUP_DIR/${GROUP}.5mC.bw"

    [ -f "$M6A_BW" ] || { echo "  ⚠️  $GROUP BigWigs missing — skipping heatmaps"; continue; }

    echo ""
    echo "  ── $GROUP ──"

    MATRIX="$GROUP_DIR/${GROUP}.TSS.gz"
    if [ ! -f "$MATRIX" ]; then
        computeMatrix reference-point \
            --referencePoint TSS \
            -b 2000 -a 2000 \
            -R "$TSS_BED" \
            -S "$M6A_BW" "$MSP_BW" "$NUC_BW" "$CPG_BW" \
            -o "$MATRIX" \
            --outFileNameMatrix "$GROUP_DIR/${GROUP}.TSS.tab" \
            --missingDataAsZero \
            --samplesLabel "m6A" "MSP" "nucleosome" "5mC" \
            -p 12 \
            && echo "    ✅  matrix" \
            || { echo "    ❌  computeMatrix failed"; continue; }
    else
        echo "    ⏭️  matrix exists"
    fi

    plotHeatmap \
        -m "$MATRIX" \
        -out "$GROUP_DIR/${GROUP}.TSS_heatmap.png" \
        --sortUsingSamples 1 \
        --sortRegions descend \
        --colorMap Greens Greens Reds Blues \
        --whatToShow 'heatmap and colorbar' \
        --plotTitle "$GROUP" \
        --heatmapHeight 12 --heatmapWidth 3 \
        --zMin 0 0 0 0 --zMax 1 1 1 1 \
        && echo "    ✅  heatmap" || echo "    ❌  plotHeatmap failed"

    plotProfile \
        -m "$MATRIX" \
        -out "$GROUP_DIR/${GROUP}.TSS_profile.png" \
        --plotTitle "$GROUP — fraction signal around TSS" \
        --samplesLabel "m6A" "MSP (accessible)" "nucleosome" "5mC" \
        --colors darkgreen green red blue \
        --yMin 0 --yMax 1 \
        --plotHeight 6 --plotWidth 8 \
        && echo "    ✅  profile" || echo "    ❌  plotProfile failed"

done

# =============================================================================
# Multi-sample comparison
# =============================================================================
echo ""
echo "============================================================"
echo " Multi-sample comparison"
echo "============================================================"

COMPARE_GROUPS="WT_Eddie WT_Rochelle cac-1 cac-2 rtt109 rtt109FLAG"

for SIGNAL in m6A msp nuc 5mC; do

    BW_LIST=""
    LABELS=""
    for GROUP in $COMPARE_GROUPS; do
        case $SIGNAL in
            m6A) F="$OUT_DIR/${GROUP}/${GROUP}.m6A.bw" ;;
            msp) F="$OUT_DIR/${GROUP}/${GROUP}.msp.bw" ;;
            nuc) F="$OUT_DIR/${GROUP}/${GROUP}.nuc.bw" ;;
            5mC) F="$OUT_DIR/${GROUP}/${GROUP}.5mC.bw" ;;
        esac
        [ -f "$F" ] && BW_LIST="$BW_LIST $F" && LABELS="$LABELS $GROUP"
    done

    [ -z "$BW_LIST" ] && { echo "  ⚠️  No BigWigs for $SIGNAL — skipping"; continue; }

    echo ""
    echo "  ── All strains: $SIGNAL ──"

    MATRIX="$OUT_DIR/allgroups.${SIGNAL}.TSS.gz"
    if [ ! -f "$MATRIX" ]; then
        computeMatrix reference-point \
            --referencePoint TSS \
            -b 2000 -a 2000 \
            -R "$TSS_BED" \
            -S $BW_LIST \
            -o "$MATRIX" \
            --missingDataAsZero \
            --samplesLabel $LABELS \
            -p 12 \
            && echo "    ✅  matrix" \
            || { echo "    ❌  computeMatrix failed"; continue; }
    else
        echo "    ⏭️  matrix exists"
    fi

    case $SIGNAL in
        m6A|msp) CMAP="Greens" ;;
        nuc)     CMAP="Reds"   ;;
        5mC)     CMAP="Blues"  ;;
    esac

    plotProfile \
        -m "$MATRIX" \
        -out "$OUT_DIR/allgroups.${SIGNAL}.profile.png" \
        --plotTitle "${SIGNAL} fraction around TSS — all strains" \
        --yMin 0 \
        --plotHeight 6 --plotWidth 10 \
        --perGroup \
        && echo "    ✅  comparison profile" || echo "    ❌  plotProfile failed"

    plotHeatmap \
        -m "$MATRIX" \
        -out "$OUT_DIR/allgroups.${SIGNAL}.heatmap.png" \
        --sortUsingSamples 1 --sortRegions descend \
        --colorMap $CMAP \
        --whatToShow 'heatmap and colorbar' \
        --plotTitle "${SIGNAL} fraction — all strains" \
        --heatmapHeight 12 --heatmapWidth 2 \
        --zMin 0 --zMax 1 \
        && echo "    ✅  comparison heatmap" || echo "    ❌  plotHeatmap failed"

done

echo ""
echo "============================================================"
echo " All done."
echo "============================================================"
find "$OUT_DIR" -name "*.png" | sort | while read -r PNG; do
    echo "  $(basename $PNG)"
done
