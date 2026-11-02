#!/bin/bash
#SBATCH --job-name=FiberTools_Pt2_CallNucs
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=FiberTools_Pt2_CallNucs.%j.out
#SBATCH --error=FiberTools_Pt2_CallNucs.%j.err

# =============================================================================
# FiberTools Part 2 — Call nucleosomes, pileups, track decorators, and QC
# =============================================================================
# For each merged BAM:
#   1. ft add-nucleosomes  → *.nucs.bam
#   2. ft pileup (m6A, 5mC, nucleosome coverage)
#   3. ft track-decorators → *.nuctracks.bed + decorated bed (for IGV)
#   4. ft qc               → *.qc.txt (QC metrics table)
# =============================================================================

set -euo pipefail

module load Miniforge3/24.11.3-0
source activate /home/ry00555/fibertools

MERGED_DIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/MergedBams"
OUT_DIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/fibertools_results"
THREADS=8

mkdir -p "$OUT_DIR"

echo "============================================================"
echo " FiberTools Pt2: Nucleosomes + pileups + decorators + QC"
echo " Input : $MERGED_DIR"
echo " Output: $OUT_DIR"
echo "============================================================"

for RUN in ONTRun9 ONTRun10; do

    RUN_IN="$MERGED_DIR/$RUN"
    RUN_OUT="$OUT_DIR/$RUN"
    mkdir -p "$RUN_OUT"

    if [ ! -d "$RUN_IN" ]; then
        echo "⚠️  $RUN_IN not found — skipping $RUN"
        continue
    fi

    echo ""
    echo "╔══════════════════════════════════════════════════════╗"
    echo "  Processing $RUN"
    echo "╚══════════════════════════════════════════════════════╝"

    for BAM in "$RUN_IN"/*_merged.bam; do
        [ -f "$BAM" ] || continue

        SAMPLE=$(basename "$BAM" .bam)
        SAMPLE_DIR="$RUN_OUT/$SAMPLE"
        mkdir -p "$SAMPLE_DIR"

        NUCS_BAM="$SAMPLE_DIR/${SAMPLE}.nucs.bam"

        echo ""
        echo "  ── $SAMPLE ──"

        # ── Step 1: add-nucleosomes ──────────────────────────────
        if [ ! -f "$NUCS_BAM" ]; then
            echo "    ft add-nucleosomes..."
            ft add-nucleosomes "$BAM" "$NUCS_BAM" \
                || { echo "    ❌ add-nucleosomes failed for $SAMPLE"; continue; }
            samtools index "$NUCS_BAM"
        else
            echo "    ⏭️  nucs.bam already exists — skipping"
        fi

        # ── Step 2: pileups ──────────────────────────────────────
        M6A_BG="$SAMPLE_DIR/${SAMPLE}.m6Apileup.bedgraph"
        CPG_BG="$SAMPLE_DIR/${SAMPLE}.5mcpileup.bedgraph"
        NUC_BG="$SAMPLE_DIR/${SAMPLE}.nucspileup.bedgraph"

        if [ ! -f "$M6A_BG" ]; then
            echo "    ft pileup — m6A..."
            ft pileup --m6a --per-base --fiber-coverage \
                --out "$M6A_BG" "$NUCS_BAM" \
                || echo "    ❌ m6A pileup failed for $SAMPLE"
        else
            echo "    ⏭️  m6A pileup exists — skipping"
        fi

        if [ ! -f "$CPG_BG" ]; then
            echo "    ft pileup — 5mC..."
            ft pileup --cpg --per-base --fiber-coverage \
                --out "$CPG_BG" "$NUCS_BAM" \
                || echo "    ❌ 5mC pileup failed for $SAMPLE"
        else
            echo "    ⏭️  5mC pileup exists — skipping"
        fi

        if [ ! -f "$NUC_BG" ]; then
            echo "    ft pileup — nucleosome coverage..."
            ft pileup --per-base --fiber-coverage \
                --out "$NUC_BG" "$NUCS_BAM" \
                || echo "    ❌ nuc pileup failed for $SAMPLE"
        else
            echo "    ⏭️  nuc pileup exists — skipping"
        fi

        # ── Step 3: track-decorators (for IGV) ───────────────────
        NUCTRACKS_BED="$SAMPLE_DIR/${SAMPLE}.nuctracks.bed"
        DECORATED_BED="$SAMPLE_DIR/decorated_${SAMPLE}.bed"

        if [ ! -f "$NUCTRACKS_BED" ]; then
            echo "    ft track-decorators..."
            ft track-decorators \
                --bed12 "$NUCTRACKS_BED" \
                --decorator "$DECORATED_BED" \
                "$NUCS_BAM" \
                || echo "    ❌ track-decorators failed for $SAMPLE"
            echo "    ✅  track-decorators done"
        else
            echo "    ⏭️  nuctracks.bed already exists — skipping"
        fi

        # ── Step 4: QC metrics ───────────────────────────────────
        # ft qc reports per-fiber metrics including:
        #   - total fibers, mean fiber length
        #   - m6A count and fraction per fiber
        #   - number of nucleosomes per fiber
        #   - MSP (methylated span) length distribution
        #   - % fibers with 0 m6A (indicates labeling failure if high)
        # --m6a-per-msp also reports m6A density within each
        # methylated span (linker/accessible element), useful for
        # checking Hia5 labeling efficiency
        QC_OUT="$SAMPLE_DIR/${SAMPLE}.qc.txt"

        if [ ! -f "$QC_OUT" ]; then
            echo "    ft qc..."
            ft qc \
                --m6a-per-msp \
                "$NUCS_BAM" \
                "$QC_OUT" \
                || echo "    ❌ qc failed for $SAMPLE"
            echo "    ✅  QC done → $(basename $QC_OUT)"
        else
            echo "    ⏭️  QC already exists — skipping"
        fi

        echo "    ✅  $SAMPLE complete"

    done
done

# =============================================================================
# Summary — print key QC metrics across all samples
# =============================================================================
echo ""
echo "============================================================"
echo " All done. QC summary:"
echo "============================================================"
echo ""

for RUN in ONTRun9 ONTRun10; do
    echo "  $RUN:"
    find "$OUT_DIR/$RUN" -name "*.qc.txt" | sort | while read -r QC; do
        SAMPLE=$(basename "$QC" .qc.txt)
        echo ""
        echo "    ── $SAMPLE ──"
        # Print the header + data rows, indented
        cat "$QC" | head -20 | sed 's/^/      /'
    done
    echo ""
done

echo ""
echo "IGV tip: load decorated_*.bed to see nucleosome footprints"
echo "per fiber. Color nucs.bam by: View → Color by → base"
echo "modification 2-color (6mA)."
echo ""
echo "QC tip: check % fibers with 0 m6A — should be <1% for a"
echo "successful Hia5 labeling. Mean m6A fraction should be ~5-6%."
echo ""
echo "Next: run FiberTools_Pt3_ExtractData.sh"
