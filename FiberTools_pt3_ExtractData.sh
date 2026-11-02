#!/bin/bash
#SBATCH --job-name=FiberTools_Pt3_Extract
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=85G
#SBATCH --time=24:00:00
#SBATCH --output=FiberTools_Pt3_Extract.%j.out
#SBATCH --error=FiberTools_Pt3_Extract.%j.err

# =============================================================================
# FiberTools Part 3 — BigWig conversion + ft extract (per-barcode)
# =============================================================================
# This script does NOT generate heatmaps. Heatmaps are made by
# FiberTools_MergeForHeatmaps.sh which runs AFTER this script using
# merged nucs.bam files for better coverage.
#
# This script does two things per sample:
#   1. bedgraph → BigWig  (for IGV and per-sample inspection)
#   2. ft extract         (per-fiber nucleosome/m6A/5mC BED files
#                          for single-fiber spacing analysis)
#
# Run order:
#   Part1 → Part2 → Part3 (this) → MergeForHeatmaps
# =============================================================================

set -euo pipefail

module load Miniforge3/24.11.3-0 ucsc/443
source activate /home/ry00555/fibertools

WORKDIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/fibertools_results"
GENOME="/home/ry00555/Research/Genomes/GenBankNcrassachromsizes.txt"

echo "============================================================"
echo " FiberTools Pt3: BigWig conversion + ft extract"
echo " (Heatmaps handled by FiberTools_MergeForHeatmaps.sh)"
echo " Working dir: $WORKDIR"
echo "============================================================"

# ── Helper: bedgraph → BigWig ────────────────────────────────────
make_bigwig() {
    local BG="$1"
    local BW="$2"

    if [ -f "$BW" ]; then
        echo "      ⏭️  $(basename $BW) already exists"
        return 0
    fi

    if [ ! -f "$BG" ]; then
        echo "      ⚠️  bedgraph not found: $(basename $BG) — skipping"
        return 0
    fi

    echo "      Converting $(basename $BG) → $(basename $BW)"

    local BG4="${BG%.bedgraph}.4col.bedgraph"

    sort -k1,1 -k2,2n "$BG" \
        | cut -f1-4 \
        | awk -v genome="$GENOME" '
            BEGIN {
                while ((getline line < genome) > 0) {
                    split(line, a, "\t"); chromlen[a[1]] = a[2]
                }
            }
            {
                if ($2 < 0) $2 = 0
                if (chromlen[$1] && $3 > chromlen[$1]) $3 = chromlen[$1]
                if ($2 < $3) print $1"\t"$2"\t"$3"\t"$4
            }
        ' > "$BG4"

    bedGraphToBigWig "$BG4" "$GENOME" "$BW" \
        && rm -f "$BG4" \
        || { echo "      ❌ bedGraphToBigWig failed for $(basename $BG)"; rm -f "$BG4"; return 1; }

    echo "      ✅  $(basename $BW)"
}

# ── Main loop ────────────────────────────────────────────────────
for RUN in ONTRun9 ONTRun10; do

    RUN_DIR="$WORKDIR/$RUN"
    if [ ! -d "$RUN_DIR" ]; then
        echo "⚠️  $RUN_DIR not found — skipping $RUN"
        continue
    fi

    echo ""
    echo "╔══════════════════════════════════════════════════════╗"
    echo "  Processing $RUN"
    echo "╚══════════════════════════════════════════════════════╝"

    for NUCS_BAM in "$RUN_DIR"/*/*.nucs.bam; do
        [ -f "$NUCS_BAM" ] || continue

        DIR=$(dirname "$NUCS_BAM")
        FILE=$(basename "$NUCS_BAM")
        SAMPLE="${FILE%.nucs.bam}"

        echo ""
        echo "  ── $SAMPLE ──"

        # ── [1] BigWig conversion ────────────────────────────────
        echo "    [1/2] BigWig conversion"

        NUC_BG="$DIR/${SAMPLE}.nucspileup.bedgraph"
        M6A_BG="$DIR/${SAMPLE}.m6Apileup.bedgraph"
        CPG_BG="$DIR/${SAMPLE}.5mcpileup.bedgraph"

        make_bigwig "$NUC_BG" "${NUC_BG%.bedgraph}.bw"
        make_bigwig "$M6A_BG" "${M6A_BG%.bedgraph}.bw"
        make_bigwig "$CPG_BG" "${CPG_BG%.bedgraph}.bw"

        # ── [2] ft extract (per-fiber BED files) ─────────────────
        echo "    [2/2] ft extract"

        NUC_BED="$DIR/${SAMPLE}_nucleosome.bed"
        M6A_BED="$DIR/${SAMPLE}_m6a.bed"
        CPG_BED="$DIR/${SAMPLE}_cpg.bed"

        if [ ! -f "$NUC_BED" ]; then
            ft extract "$NUCS_BAM" \
                --nuc "$NUC_BED" \
                --m6a "$M6A_BED" \
                --cpg "$CPG_BED" \
                --threads 8 \
                && echo "      ✅  ft extract done" \
                || echo "      ❌  ft extract failed for $SAMPLE"
        else
            echo "      ⏭️  ft extract already done"
        fi

    done
done

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "============================================================"
echo " Part 3 complete."
echo "============================================================"
for RUN in ONTRun9 ONTRun10; do
    echo ""
    echo "  $RUN BigWigs:"
    find "$WORKDIR/$RUN" -name "*.bw" | sort | while read -r BW; do
        SIZE=$(du -h "$BW" | cut -f1)
        echo "    $SIZE  $(basename $BW)"
    done
done

echo ""
echo "Next step: sbatch FiberTools_MergeForHeatmaps.sh"
echo "  This will merge nucs.bams by strain, generate pileups,"
echo "  convert to BigWig, and produce TSS heatmaps."
