#!/bin/bash
#SBATCH --job-name=MapCutAndRun_JS
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jms71286@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../MapCutAndRun.%j.out
#SBATCH --error=../MapCutAndRun.%j.err

cd /scratch/jms71286/cutandrun
source configJS.txt

OUTDIR="/scratch/jms71286/cutandrun/output"

#process reads using trimGalore
module load Trim_Galore/0.6.10-GCCcore-12.3.0
#!/bin/bash
#SBATCH --job-name=MapCutAndRun_JS
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jms71286@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --output=../MapCutAndRun.%j.out
#SBATCH --error=../MapCutAndRun.%j.err


cd /scratch/jms71286/cutandrun
source configJS.txt

# Required configuration:
# FASTQ="/path/to/raw_fastq_directory"
# GENOME="/path/to/bwa_index_prefix"
# BIN=10
# SMOOTH=25

: "${FASTQ:?Set FASTQ in configJS.txt}"
: "${GENOME:?Set GENOME in configJS.txt}"
: "${BIN:?Set BIN in configJS.txt}"
: "${SMOOTH:?Set SMOOTH in configJS.txt}"

THREADS="${SLURM_CPUS_PER_TASK:-24}"
OUTDIR="/scratch/jms71286/cutandrun/output"

# Make an unmatched filename pattern expand to an empty array.
shopt -s nullglob

# -------------------------
# Trim reads, one pair at a time
# -------------------------
module load Trim_Galore/0.6.10-GCCcore-12.3.0

raw_r1_files=( "$FASTQ"/*_R1_001.fastq.gz )

if (( ${#raw_r1_files[@]} == 0 )); then
    echo "ERROR: No *_R1_001.fastq.gz files found in $FASTQ" >&2
    exit 1
fi

for read1 in "${raw_r1_files[@]}"; do
    read2="${read1%_R1_001.fastq.gz}_R2_001.fastq.gz"

    if [[ ! -s "$read2" ]]; then
        echo "ERROR: Missing or empty R2 file: $read2" >&2
        exit 1
    fi

    trim_galore \
        --paired \
        --length 20 \
        --fastqc \
        --gzip \
        --output_dir "$OUTDIR/TrimmedReads" \
        "$read1" "$read2"
done

# -------------------------
# Load downstream tools
# -------------------------
module load SAMtools/1.23.1-GCC-13.3.0
module load BWA/0.7.17-GCCcore-11.3.0
module load deepTools/3.5.5-gfbf-2023a
module load MACS3/3.0.1-gfbf-2023a

# -------------------------
# Map, create tracks, call peaks
# -------------------------
trimmed_r1_files=(
    "$OUTDIR"/TrimmedReads/*_R1_001_val_1.fq.gz
)

if (( ${#trimmed_r1_files[@]} == 0 )); then
    echo "ERROR: No validated R1 files found after trimming" >&2
    exit 1
fi

for read1 in "${trimmed_r1_files[@]}"; do
    read2="${read1%_R1_001_val_1.fq.gz}_R2_001_val_2.fq.gz"

    if [[ ! -s "$read2" ]]; then
        echo "ERROR: Missing or empty trimmed R2 file: $read2" >&2
        exit 1
    fi

    file="${read1##*/}"
    name="${file%_R1_001_val_1.fq.gz}"

    bam="$OUTDIR/SortedBamFiles/${name}.bam"
    bigwig="$OUTDIR/BigWigs/${name}"

    echo "Processing: $name"

    # Use one additional sorting thread while BWA runs.
    bwa mem \
        -M \
        -v 3 \
        -t "$((THREADS - 2))" \
        "$GENOME" \
        "$read1" \
        "$read2" |
        samtools sort \
            -@ 1 \
            -T "$OUTDIR/SortedBamFiles/${name}.tmp.${SLURM_JOB_ID:-$$}" \
            -o "$bam" \
            -

    samtools index "$bam"

    bamCoverage \
        -p "$THREADS" \
        --binSize "$BIN" \
        --normalizeUsing BPM \
        --minMappingQuality 10 \
        --smoothLength "$SMOOTH" \
        --outFileFormat bigwig \
        --bam "$bam" \
        --outFileName \
        "${bigwig}.bin_${BIN}.smooth_${SMOOTH}.Bulk.bw"

    bamCoverage \
        -p "$THREADS" \
        --MNase \
        --binSize 1 \
        --normalizeUsing BPM \
        --minMappingQuality 10 \
        --smoothLength 25 \
        --outFileFormat bigwig \
        --bam "$bam" \
        --outFileName \
        "${bigwig}.bin_1.smooth_25.MNase.bw"

    macs3 callpeak \
        -t "$bam" \
        -f BAMPE \
        -n "$name" \
        --broad \
        -g 41037538 \
        --broad-cutoff 0.1 \
        --outdir "$OUTDIR/Peaks" \
        --min-length 800 \
        --max-gap 500 \
        --nolambda
done
