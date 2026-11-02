#!/bin/bash
#SBATCH --job-name=FiberTools_Pt1_MergeBams
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=FiberTools_Pt1_MergeBams.%j.out
#SBATCH --error=FiberTools_Pt1_MergeBams.%j.err

# =============================================================================
# FiberTools Part 1 — Merge split BAMs within each run (runs kept separate)
# =============================================================================
# For each run (ONTRun9, ONTRun10) and each barcode, merges all the split
# BAM files (e.g. _0.bam, _1.bam ... _20.bam) into one merged BAM named
# by strain/sample from the metadata.
#
# ONTRun9 and ONTRun10 are kept SEPARATE — each produces its own merged BAM.
#
# Input layout:
#   InputBams/
#     ONTRun9/barcode01/FBE00208_pass_barcode01_*_0.bam  ...
#     ONTRun10/barcode01/...
#
# Output layout:
#   MergedBams/
#     ONTRun9/
#       ONTRun9_WT_Eddie_barcode01_merged.bam
#       ONTRun9_cac-1_Eddie_barcode03_merged.bam  ...
#     ONTRun10/
#       ONTRun10_WT_Eddie_barcode01_merged.bam  ...
# =============================================================================

set -euo pipefail

module load Miniforge3/24.11.3-0
source activate /home/ry00555/fibertools

BASE_DIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/InputBams"
OUT_DIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/MergedBams"
SORTED_TMP="${OUT_DIR}/sorted_q10_tmp"

mkdir -p "$OUT_DIR/ONTRun9"
mkdir -p "$OUT_DIR/ONTRun10"
mkdir -p "$SORTED_TMP"

# =============================================================================
# Barcode → sample name mapping (from FiberMeta.csv)
# Barcodes 23 & 24 are empty in metadata — skipped.
# =============================================================================
declare -A SAMPLE_MAP=(
    [barcode01]="WT_Eddie_barcode01"
    [barcode02]="WT_Eddie_barcode02"
    [barcode03]="cac-1_Eddie_barcode03"
    [barcode04]="cac-1_Eddie_barcode04"
    [barcode05]="cac-2_Eddie_barcode05"
    [barcode06]="cac-2_Eddie_barcode06"
    [barcode07]="WT_Rochelle_barcode07"
    [barcode08]="WT_Rochelle_barcode08"
    [barcode09]="rtt109_Rochelle_barcode09"
    [barcode10]="rtt109_Rochelle_barcode10"
    [barcode11]="rtt109FLAG_Rochelle_barcode11"
    [barcode12]="rtt109FLAG_Rochelle_barcode12"
    [barcode13]="WT_Eddie_barcode13"
    [barcode14]="cac-1_Eddie_barcode14"
    [barcode15]="cac-2_Eddie_barcode15"
    [barcode16]="gDNA_Hia5_ctrl_Eddie_barcode16"
    [barcode17]="WizardgDNA_HMW_Eddie_barcode17"
    [barcode18]="WT_Rochelle_barcode18"
    [barcode19]="rtt109_Rochelle_barcode19"
    [barcode20]="rtt109FLAG_Rochelle_barcode20"
    [barcode21]="gDNA_Hia5_ctrl_Rochelle_barcode21"
    [barcode22]="WizardgDNA_HMW_Rochelle_barcode22"
)

echo "============================================================"
echo " FiberTools Pt1: Merge split BAMs (runs kept separate)"
echo " Output dir: $OUT_DIR"
echo "============================================================"

for RUN in ONTRun9 ONTRun10; do

    echo ""
    echo "╔══════════════════════════════════════════════════════╗"
    echo "  Processing $RUN"
    echo "╚══════════════════════════════════════════════════════╝"

    for BARCODE in "${!SAMPLE_MAP[@]}"; do

        SAMPLE="${SAMPLE_MAP[$BARCODE]}"
        SRC_DIR="$BASE_DIR/$RUN/$BARCODE"
        FINAL_BAM="$OUT_DIR/$RUN/${RUN}_${SAMPLE}_merged.bam"

        # Skip if source barcode dir doesn't exist for this run
        if [ ! -d "$SRC_DIR" ]; then
            echo "  ⚠️  $RUN/$BARCODE — directory not found, skipping"
            continue
        fi

        # Skip if already done
        if [ -f "$FINAL_BAM" ] && [ -f "${FINAL_BAM}.bai" ]; then
            echo "  ⏭️  Skipping ${RUN}/${BARCODE} — already merged"
            continue
        fi

        BAM_COUNT=$(find "$SRC_DIR" -maxdepth 1 -name "*.bam" ! -name "*.bai" | wc -l)
        if [ "$BAM_COUNT" -eq 0 ]; then
            echo "  ⚠️  $RUN/$BARCODE — no BAMs found, skipping"
            continue
        fi

        echo ""
        echo "  ── $RUN / $BARCODE ($SAMPLE) ──"
        echo "     Found $BAM_COUNT split BAM(s)"

        SORTED_BAMS=()

        for BAM in "$SRC_DIR"/*.bam; do
            [ -f "$BAM" ] || continue

            BAM_BASE=$(basename "$BAM" .bam)
            SORTED_BAM="$SORTED_TMP/${RUN}_${SAMPLE}_${BAM_BASE}_q10_sorted.bam"

            if [ ! -f "$SORTED_BAM" ]; then
                samtools view -b -q 10 "$BAM" \
                    | samtools sort -@ 4 -o "$SORTED_BAM"
            fi
            SORTED_BAMS+=("$SORTED_BAM")
        done

        echo "     Merging ${#SORTED_BAMS[@]} sorted BAM(s)..."

        if [ ${#SORTED_BAMS[@]} -eq 1 ]; then
            cp "${SORTED_BAMS[0]}" "$FINAL_BAM"
        else
            samtools merge -f -@ 4 "$FINAL_BAM" "${SORTED_BAMS[@]}"
        fi

        samtools index "$FINAL_BAM"
        echo "     ✅  $(basename $FINAL_BAM)"

    done
done

# Clean up temp sorted BAMs
echo ""
echo "Cleaning up temp files..."
rm -rf "$SORTED_TMP"

# Summary
echo ""
echo "============================================================"
echo " All done. Merged BAMs:"
echo "============================================================"
for RUN in ONTRun9 ONTRun10; do
    echo ""
    echo "  $RUN:"
    find "$OUT_DIR/$RUN" -name "*_merged.bam" | sort | while read -r BAM; do
        SIZE=$(du -h "$BAM" | cut -f1)
        echo "    $SIZE  $(basename $BAM)"
    done
done

echo ""
echo "Next: update IN_DIR in FiberTools_Pt2_CallNucs.sh to:"
echo "  $OUT_DIR/ONTRun9   and/or   $OUT_DIR/ONTRun10"
