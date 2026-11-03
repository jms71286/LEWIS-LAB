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

for QC in "${RUN_RESULTS}"/*/*.qc.tsv; do
    [[ -f "$QC" ]] || continue

    SAMPLE="$(basename "$QC" .qc.tsv)"

    echo
    echo "Sample: $SAMPLE"
    head -n 20 "$QC" | sed 's/^/  /'
done

echo
echo "Merged BAMs:"
find "$MERGED_DIR" -maxdepth 1 -type f -name "*_merged.bam" \
    -printf '%f\n' |
sort

echo
echo "FiberTools results:"
find "$RUN_RESULTS" -maxdepth 2 -type f \
    -printf '%P\n' |
sort

echo
echo "Cleaning temporary sorted BAMs..."
rm -rf -- "$SORTED_TMP"

echo "Done."
