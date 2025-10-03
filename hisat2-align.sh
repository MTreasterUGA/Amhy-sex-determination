#!/bin/bash
#SBATCH --partition=batch       # other queues for >7 days, more memory, gpu
#SBATCH --job-name=hisat2-align
#SBATCH --ntasks=1
#SBATCH --time=24:00:00          # D-HH:MM:SS format
#SBATCH --mem=24G              # memory per node
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=mjt86304@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/mjt86304/slurm_out/%j-%x.out
#SBATCH --error=/scratch/mjt86304/slurm_out/%j-%x.err
#SBATCH --array=1-36


set -euo pipefail
date| tee /dev/stderr

cd /scratch/mjt86304/embryo_rna
mkdir -p aligned
mkdir -p summary

if [[ ! -f "samples.txt" ]]; then
    echo "Error: samples.txt file not found."
    exit 1
fi

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" samples.txt)
if [[ -z "${SAMPLE}" ]]; then
    echo "Error: No sample found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

INPUT="trimmed/${SAMPLE}_"
INPUT1="${INPUT}R1_PE_trimmed.fastq.gz"
INPUT2="${INPUT}R1_SE_trimmed.fastq.gz"
INPUT3="${INPUT}R2_PE_trimmed.fastq.gz"
INPUT4="${INPUT}R2_SE_trimmed.fastq.gz"

OUTPUT="aligned/${SAMPLE}_"
OUTPUT1="${OUTPUT}PE.bam"
OUTPUT2="${OUTPUT}SE1.bam"
OUTPUT3="${OUTPUT}SE2.bam"

if [[ ! -f "${INPUT1}" || \
    ! -f "${INPUT2}" || \
    ! -f "${INPUT3}" || \
    ! -f "${INPUT4}" ]]; then
    echo "Error: Input files for sample ${SAMPLE} not found."
    exit 1
fi

INDEX="/home/mjt86304/refgenome/hisat2_build/gac_uga_v5_tran"       #Hisat2 index path, include characters before .1.ht2
if [[ ! -e "${INDEX}.1.ht2" ]]; then
    echo "Error: HISAT2 index not found at ${INDEX}"
    exit 1
fi

ml HISAT2/3n-20201216-gompi-2022a
ml SAMtools/1.17-GCC-12.2.0

echo "start ${SAMPLE} paired" | tee /dev/stderr

hisat2 -q --threads 16 -x "${INDEX}" --rna-strandness RF \
    --summary-file "summary/${SAMPLE}_PE.txt" \
    -1 "$INPUT1" -2 "$INPUT3" | \
    samtools view -bS --threads 16 | samtools sort --threads 16 -o "$OUTPUT1"

echo "start ${SAMPLE} SE1" | tee /dev/stderr

hisat2 -q --threads 16 -x "${INDEX}" --rna-strandness R \
    --summary-file "summary/${SAMPLE}_SE1.txt" \
    -U "$INPUT2" | \
    samtools view -bS --threads 16 | samtools sort --threads 16 -o "$OUTPUT2"

echo "start ${SAMPLE} SE2" | tee /dev/stderr

hisat2 -q --threads 16 -x "${INDEX}" --rna-strandness F \
    --summary-file "summary/${SAMPLE}_SE2.txt" \
    -U "$INPUT4" | \
    samtools view -bS --threads 16 | samtools sort --threads 16 -o "$OUTPUT3"

for i in "${OUTPUT}"*
do
    echo "indexing $i" | tee /dev/stderr
    samtools index -@ 16 "$i"
done
