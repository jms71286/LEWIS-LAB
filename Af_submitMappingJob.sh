#!/bin/bash
# Usage:  bash submitMappingJob.sh sampleList.txt
# Submits one MapData.sh job per SampleID in the list.

if [[ -z "$1" ]]; then
  echo "
  ################
  ERROR: You must include the path to your sample list in the command line call.
  eg. bash submitMappingJob.sh sampleList.txt
  ################
  "
  exit 1
fi

fastqPath="/scratch/ry00555/Afum/CEA10/FastQ"   # renamed <SampleID>_R1/_R2.fastq.gz files
outdir="/scratch/ry00555/Afum/CEA10"            # TrimmedFastQs, bamFiles, bigWig, Peaks, logs go here

mkdir -p "${outdir}/logs" "${outdir}/TrimmedFastQs" "${outdir}/bamFiles" \
         "${outdir}/bigWig" "${outdir}/Peaks"

while read -r line || [[ -n "$line" ]]; do
  line="${line//$'\r'/}"            # strip Windows/Excel line endings
  [[ -z "$line" || "$line" == \#* ]] && continue

  if [[ ! -s "${fastqPath}/${line}_R1.fastq.gz" ]]; then
    echo "SKIPPED $line: ${fastqPath}/${line}_R1.fastq.gz not found"
    continue
  fi

  sbatch \
    --job-name="${line}" \
    --output="${outdir}/logs/%x.out" \
    --error="${outdir}/logs/%x.err" \
    --export=ALL,accession="${line}",fastqPath="${fastqPath}",outdir="${outdir}" \
    MapAfData.sh
  echo "$line mapping job submitted"
  sleep 1
done < "$1"
