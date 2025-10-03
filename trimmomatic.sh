#!/bin/bash
#SBATCH --partition=batch       # other queues for >7 days, more memory, gpu
#SBATCH --job-name=trimmomatic
#SBATCH --ntasks=1
#SBATCH --time=48:00:00         # D-HH:MM:SS format
#SBATCH --mem=4G               # memory per node
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=mjt86304@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/mjt86304/slurm_out/%j-%x.out
#SBATCH --error=/scratch/mjt86304/slurm_out/%j-%x.err
#SBATCH --array=1-36

set -euo pipefail
date | tee /dev/stderr

cd /scratch/mjt86304/embryo_rna
mkdir -p trimmed

if [[ ! -f "samples.txt" ]]; then
    echo "Error: samples.txt file not found."
    exit 1
fi

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" samples.txt)
if [[ -z "${SAMPLE}" ]]; then
    echo "Error: No sample found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

INPUT="raw/${SAMPLE}_"
INPUT1="${INPUT}R1_001.fastq.gz"
INPUT2="${INPUT}R2_001.fastq.gz"

OUTPUT="trimmed/${SAMPLE}_"
OUTPUT1="${OUTPUT}R1_PE_trimmed.fastq.gz"
OUTPUT2="${OUTPUT}R1_SE_trimmed.fastq.gz"
OUTPUT3="${OUTPUT}R2_PE_trimmed.fastq.gz"
OUTPUT4="${OUTPUT}R2_SE_trimmed.fastq.gz"

if [[ ! -f "${INPUT1}" || ! -f "${INPUT2}" ]]; then
    echo "Error: Input files for sample ${SAMPLE} not found."
    exit 1
fi

echo $SAMPLE | tee /dev/stderr

ml Trimmomatic/0.39-Java-13

java -jar "${EBROOTTRIMMOMATIC}/trimmomatic-0.39.jar" PE -threads 4 \
    "${INPUT1}" "${INPUT2}" \
    "${OUTPUT1}" "${OUTPUT2}" "${OUTPUT3}" "${OUTPUT4}" \
    ILLUMINACLIP:"${EBROOTTRIMMOMATIC}/adapters/TruSeq3-PE.fa":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
