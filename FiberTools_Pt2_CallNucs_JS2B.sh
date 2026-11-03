#!/bin/bash
#SBATCH --job-name=FiberTools_Pt2_CallNucs
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jms71286@uga.edu
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
#!/usr/bin/env bash
#!/usr/bin/env bash
#load miniforge and samtools
set -Eeuo pipefail
trap 'echo "ERROR at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

module load Miniforge3/24.11.3-0
module load SAMtools/1.23.1-GCC-13.3.0

RUN="ONTRun13"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
CONDA_ENV="fibertools-rs"

BASE_DIR="/lustre2/scratch/jms71286/FIBER/ONTRun13/bam_pass"
MERGED_BASE="/lustre2/scratch/jms71286/FIBER/ONTRun13/MergedBams"
MERGED_DIR="${MERGED_BASE}/${RUN}"
RESULTS_BASE="/lustre2/scratch/jms71286/FIBER/ONTRun13/fibertools_results"
RUN_RESULTS="${RESULTS_BASE}/${RUN}"
SORTED_TMP="${MERGED_BASE}/sorted_${RUN}_tmp"

# Use conda run rather than conda activate in a noninteractive job.
FT=(conda run --no-capture-output -n "$CONDA_ENV" ft)

shopt -s nullglob
BAMS=("${MERGED_DIR}"/*_merged.bam)

if (( ${#BAMS[@]} == 0 )); then
    echo "ERROR: no *_merged.bam files found in $MERGED_DIR" >&2
    exit 1
fi

for BAM in "${BAMS[@]}"; do
    SAMPLE="$(basename "$BAM" .bam)"
    SAMPLE_DIR="${RUN_RESULTS}/${SAMPLE}"

    NUCS_BAM="${SAMPLE_DIR}/${SAMPLE}.nucs.bam"
    M6A_PILEUP="${SAMPLE_DIR}/${SAMPLE}.m6a.pileup.tsv"
    CPG_PILEUP="${SAMPLE_DIR}/${SAMPLE}.5mc.pileup.tsv"
    NUC_PILEUP="${SAMPLE_DIR}/${SAMPLE}.nuc.pileup.tsv"
    NUCTRACKS_BED="${SAMPLE_DIR}/${SAMPLE}.nuctracks.bed"
    DECORATED_BED="${SAMPLE_DIR}/${SAMPLE}.decorators.bed"
    QC_OUT="${SAMPLE_DIR}/${SAMPLE}.qc.tsv"
    VALIDATE_OUT="${SAMPLE_DIR}/${SAMPLE}.validate.txt"

    mkdir -p "$SAMPLE_DIR"

    echo
    echo "------------------------------------------------------------"
    echo "Sample: $SAMPLE"
    echo "Input:  $BAM"
    echo "Output: $SAMPLE_DIR"
    echo "------------------------------------------------------------"

    samtools quickcheck -v "$BAM"

    if [[ ! -s "${BAM}.bai" ]]; then
        samtools index -@ "$THREADS" "$BAM"
    fi

    # Test whether the merged BAM contains the annotations FiberTools needs.
    echo "Validating m6A annotations in input BAM..."

    if ! "${FT[@]}" validate --reads 5000 "$BAM" \
        >"$VALIDATE_OUT" 2>&1; then
        cat "$VALIDATE_OUT" >&2
        echo "ERROR: FiberTools validation failed for $BAM" >&2
        echo "The BAM may not contain usable m6A MM/ML tags." >&2
        exit 1
    fi

    cat "$VALIDATE_OUT"

    # Add nucleosome and MSP annotations.
    if [[ ! -s "$NUCS_BAM" ]]; then
        echo "Adding nucleosome and MSP calls..."

        "${FT[@]}" \
            --threads "$THREADS" \
            add-nucleosomes \
            "$BAM" \
            "$NUCS_BAM"

        samtools quickcheck -v "$NUCS_BAM"
        samtools index -@ "$THREADS" "$NUCS_BAM"
    else
        echo "Nucleosome BAM already exists. Skipping."
    fi

    # Validate m6A plus nucleosome annotations.
    echo "Validating nucleosome BAM..."
    "${FT[@]}" validate --reads 5000 "$NUCS_BAM"

    # m6A-only pileup.
    if [[ ! -s "$M6A_PILEUP" ]]; then
        echo "Creating m6A pileup..."

        "${FT[@]}" \
            --threads "$THREADS" \
            pileup \
            --m6a \
            --no-msp \
            --no-nuc \
            --per-base \
            --fiber-coverage \
            --out "$M6A_PILEUP" \
            "$NUCS_BAM"
    else
        echo "m6A pileup already exists. Skipping."
    fi

    # 5mC-only pileup.
    # This can be empty or uninformative if the BAM does not contain 5mC calls.
    if [[ ! -s "$CPG_PILEUP" ]]; then
        echo "Creating 5mC pileup..."

        "${FT[@]}" \
            --threads "$THREADS" \
            pileup \
            --cpg \
            --no-msp \
            --no-nuc \
            --per-base \
            --fiber-coverage \
            --out "$CPG_PILEUP" \
            "$NUCS_BAM"
    else
        echo "5mC pileup already exists. Skipping."
    fi

    # Nucleosome pileup. MSPs are excluded.
    if [[ ! -s "$NUC_PILEUP" ]]; then
        echo "Creating nucleosome pileup..."

        "${FT[@]}" \
            --threads "$THREADS" \
            pileup \
            --no-msp \
            --per-base \
            --fiber-coverage \
            --out "$NUC_PILEUP" \
            "$NUCS_BAM"
    else
        echo "Nucleosome pileup already exists. Skipping."
    fi

    # Decorator tracks for IGV.
    if [[ ! -s "$NUCTRACKS_BED" || ! -s "$DECORATED_BED" ]]; then
        echo "Creating IGV decorator tracks..."

        "${FT[@]}" \
            --threads "$THREADS" \
            track-decorators \
            --bed12 "$NUCTRACKS_BED" \
            --decorator "$DECORATED_BED" \
            "$NUCS_BAM"
    else
        echo "Decorator tracks already exist. Skipping."
    fi

    # QC metrics.
    if [[ ! -s "$QC_OUT" ]]; then
        echo "Running FiberTools QC..."

        "${FT[@]}" \
            --threads "$THREADS" \
            qc \
            --m6a-per-msp \
            "$NUCS_BAM" \
            "$QC_OUT"
    else
        echo "QC output already exists. Skipping."
    fi

    echo "Completed: $SAMPLE"
done
