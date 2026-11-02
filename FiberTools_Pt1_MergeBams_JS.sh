#!/bin/bash
#SBATCH --job-name=FiberTools_Pt1_MergeBams
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jms71286@uga.edu
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
# Input layout:
#   InputBams/
#     ONTRun13/barcode01/FBE00208_pass_barcode01_*_0.bam  ...
#     ONTRun13/barcode01/...
#
# Output layout:
#   MergedBams/
#     ONTRun13/
#       ONTRun13_CEA10_conidia1_barcode07_merged.bam  ...
#       ONTRun13_cac-1_Eddie_barcode03_merged.bam  ...
#
# =============================================================================

#download fiber tools first in terminal: curl --proto '=https' --tlsv1.2 -LsSf https://github.com/fiberseq/fibertools-rs/releases/download/v0.13.1/fibertools-rs-installer.sh | sh

set -euo pipefail

module load Miniforge3/24.11.3-0
#curl --proto '=https' --tlsv1.2 -LsSf https://github.com/fiberseq/fibertools-rs/releases/download/v0.13.1/fibertools-rs-installer.sh | sh
#conda install bioconda::fibertools-rs
#conda create --name myenv python=3.10 ?
#source activate /home/jms71286/fibertools

BASE_DIR="/lustre2/scratch/jms71286/FIBER/ONTRun13/bam_pass"
OUT_DIR="/lustre2/scratch/jms71286/FIBER/ONTRun13/MergedBams"
SORTED_TMP="${OUT_DIR}/sorted_ONTRun13_tmp"

mkdir -p "$SORTED_TMP"

# =============================================================================
# Barcode → sample name mapping (from FiberMeta.csv)
# Barcodes 7-9.
# =============================================================================
declare -A SAMPLE_MAP=(
    [barcode07]="CEA10_conidia1_barcode07"
    [barcode08]="CEA10_conidia2_barcode08"
    [barcode09]="CEA10_conidia3_barcode09"
)

echo "============================================================"
echo " FiberTools Pt1: Merge split BAMs (runs kept separate)"
echo " Output dir: $OUT_DIR"
echo "============================================================"

for RUN in ONTRun13; do

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
            SORTED_BAM="$SORTED_TMP/${RUN}_${SAMPLE}_${BAM_BASE}_ONTRun13_sorted.bam"

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
for RUN in ONTRun13; do
    echo ""
    echo "  $RUN:"
    find "$OUT_DIR/$RUN" -name "*_merged.bam" | sort | while read -r BAM; do
        SIZE=$(du -h "$BAM" | cut -f1)
        echo "    $SIZE  $(basename $BAM)"
    done
done

echo ""
echo "Next: update IN_DIR in FiberTools_Pt2_CallNucs_JS.sh to:"
echo "  $OUT_DIR/ONTRun13"
