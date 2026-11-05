#!/bin/bash
#SBATCH --job-name=cutANDrun
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jms71286@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=./log.%j.out
#SBATCH --error=./log.%j.err

cd $SLURM_SUBMIT_DIR
source configJS.txt

GENOME="/home/zlewis/Genomes/AfumagatisGenome/Af_CEA10/Afum_A1163_edit6.fasta"
OUTDIR="/scratch/jms71286/cutandrun/output"
FASTQ="/scratch/jms71286/cutandrun/FastqFiles"
THREADS="12"
BIN="25"
SMOOTH="75"
#
module load Trim_Galore/0.6.10-GCCcore-12.3.0
module load Trim_Galore/0.6.10-GCCcore-12.3.0
module load SAMtools/1.23.1-GCC-13.3.0
module load BWA/0.7.17-GCCcore-11.3.0
module load deepTools/3.5.5-gfbf-2023a
module load MACS3/3.0.1-gfbf-2023a
#process reads using trimGalore
trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz
#
FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1.fq\.gz" #Don't forget the *
#
for f in "${OUTDIR}"/TrimmedReads/*_R1_001_val_1.fq.gz; do
    name=$(basename "$f" _R1_001_val_1.fq.gz)
    read2="${OUTDIR}/TrimmedReads/${name}_R2_001_val_2.fq.gz"
    bam="${OUTDIR}/SortedBamFiles/${name}.bam"
    bigwig="${OUTDIR}/BigWigs/${name}"

    bwa mem -t "$THREADS" "$GENOME" "$f" "$read2" |
        samtools view -bhSu - |
        samtools sort -@ "$THREADS" \
            -T "${OUTDIR}/SortedBamFiles/${name}.tmp" \
            -o "$bam" -

# 	# File Vars
# 	#use sed to get the name of the second read matching the input file
	read2=$(echo "$f" | sed 's/R1_001_val_1.fq\.gz/R2_001_val_2.fq\.gz/g')
	#variable for naming bam file
 	bam="${OUTDIR}/SortedBamFiles/${name}.bam"
	#variable name for bigwig output
	bigwig="${OUTDIR}/BigWigs/${name}"
	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
#
samtools index "$bam"
samtools flagstat "$bam" > "${OUTDIR}/QC/${name}.flagstat.txt"
samtools sort -@ "$THREADS" -T "${OUTDIR}/SortedBamFiles/${name}.tmp" -o "$bam" -
samtools flagstat "$bam" > "${OUTDIR}/QC/${name}.flagstat.txt"

#Plot all reads with deeptools
bamCoverage -p "$THREADS" -bs "$BIN" --normalizeUsing BPM --minMappingQuality 10 --smoothLength "$SMOOTH" -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"

#plot mononucleosomes
bamCoverage -p "$THREADS" --MNase -bs 1 --normalizeUsing BPM --minMappingQuality 10 --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"
#
macs3 callpeak -t "${bam}" -f BAMPE -n "${name}" --broad -g 41037538 --broad-cutoff 0.1 --outdir "${OUTDIR}/Peaks" --min-length 800 --max-gap 500 --nolambda
done
done
