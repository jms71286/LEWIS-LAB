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

echo "============================================================"
echo "FiberTools pipeline"
echo "Run:        $RUN"
echo "Input:      $BASE_DIR"
echo "Merged:     $MERGED_DIR"
echo "Results:    $RUN_RESULTS"
echo "Threads:    $THREADS"
echo "============================================================"

# ---------------------------------------------------------------------------
# Environment checks
# ---------------------------------------------------------------------------

command -v samtools >/dev/null 2>&1 || {
    echo "ERROR: samtools was not found after loading the module." >&2
    exit 1
}

conda env list | awk '{print $1}' | grep -qx "$CONDA_ENV" || {
    echo "ERROR: Conda environment '$CONDA_ENV' does not exist." >&2
    echo "Create it once from the login node with:" >&2
    echo "  conda create -y -n $CONDA_ENV -c conda-forge -c bioconda fibertools-rs=0.13.1" >&2
    exit 1
}

"${FT[@]}" --version

# ---------------------------------------------------------------------------
# Part 1: filter, sort, merge, and index BAMs
# ---------------------------------------------------------------------------

echo
echo "============================================================"
echo "Part 1: Merge split BAMs"
echo "============================================================"

for BARCODE in barcode07 barcode08 barcode09; do
    SAMPLE="${SAMPLE_MAP[$BARCODE]}"
    SRC_DIR="${BASE_DIR}/${BARCODE}"
    FINAL_BAM="${MERGED_DIR}/${RUN}_${SAMPLE}_merged.bam"

    echo
    echo "Processing $BARCODE: $SAMPLE"

    if [[ ! -d "$SRC_DIR" ]]; then
        echo "WARNING: source directory not found: $SRC_DIR" >&2
        continue
    fi

    mapfile -d '' INPUT_BAMS < <(
        find "$SRC_DIR" -maxdepth 1 -type f -name "*.bam" -print0 | sort -z
    )

    if (( ${#INPUT_BAMS[@]} == 0 )); then
        echo "WARNING: no BAM files found in $SRC_DIR" >&2
        continue
    fi

    if [[ -s "$FINAL_BAM" && -s "${FINAL_BAM}.bai" ]]; then
        echo "Merged BAM and index already exist. Skipping merge."
    else
        SORTED_BAMS=()

        for BAM in "${INPUT_BAMS[@]}"; do
            BAM_BASE="$(basename "$BAM" .bam)"
            SORTED_BAM="${SORTED_TMP}/${RUN}_${SAMPLE}_${BAM_BASE}.q10.sorted.bam"

            if [[ ! -s "$SORTED_BAM" ]]; then
                echo "Filtering and sorting: $(basename "$BAM")"

                samtools view \
                    -@ "$THREADS" \
                    -b \
                    -q 10 \
                    -F 2308 \
                    "$BAM" |
                samtools sort \
                    -@ "$THREADS" \
                    -o "$SORTED_BAM" \
                    -
            fi

            samtools quickcheck -v "$SORTED_BAM"
            SORTED_BAMS+=("$SORTED_BAM")
        done

        if (( ${#SORTED_BAMS[@]} == 1 )); then
            cp -- "${SORTED_BAMS[0]}" "$FINAL_BAM"
        else
            samtools merge \
                -f \
                -@ "$THREADS" \
                "$FINAL_BAM" \
                "${SORTED_BAMS[@]}"
        fi

        samtools quickcheck -v "$FINAL_BAM"
        samtools index -@ "$THREADS" "$FINAL_BAM"
    fi

    echo "Created: $FINAL_BAM"
done
