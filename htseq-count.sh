#!/bin/bash
#SBATCH --partition=batch       # other queues for >7 days, more memory, gpu
#SBATCH --job-name=htseq-count
#SBATCH --ntasks=1
#SBATCH --time=72:00:00         # D-HH:MM:SS format
#SBATCH --mem=20G              # memory per node
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=mjt86304@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/mjt86304/slurm_out/%j-%x.out
#SBATCH --error=/scratch/mjt86304/slurm_out/%j-%x.err
#SBATCH --array=1-36


set -euo pipefail
date | tee /dev/stderr

cd /scratch/mjt86304/embryo_rna
mkdir -p counts

if [[ ! -f "samples.txt" ]]; then
    echo "Error: samples.txt file not found."
    exit 1
fi

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" samples.txt)
if [[ -z "${SAMPLE}" ]]; then
    echo "Error: No sample found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

INPUT="aligned/${SAMPLE}_"
INPUT1="${INPUT}PE.bam"
INPUT2="${INPUT}SE1.bam"
INPUT3="${INPUT}SE2.bam"

OUTPUT="counts/${SAMPLE}_"
OUTPUT1="${OUTPUT}PE_counts.txt"
OUTPUT2="${OUTPUT}SE1_counts.txt"
OUTPUT3="${OUTPUT}SE2_counts.txt"

if [[ ! -f "${INPUT1}" || \
    ! -f "${INPUT2}" || \
    ! -f "${INPUT3}" ]]; then
    echo "Error: Input files for sample ${SAMPLE} not found."
    exit 1
fi

REFGENOME="/home/mjt86304/refgenome/GCF_016920845.1_GAculeatus_UGA_version5_genomic.gtf"

if [[ ! -f "${REFGENOME}" ]]; then
    echo "Error: Reference Genome file not found."
    exit 1
fi

echo $SAMPLE | tee /dev/stderr

ml HTSeq/2.0.2-foss-2022a

# for Truseq stranded prep kit, use --stranded=reverse for paired-end reads and R1 single-end reads
# use --stranded=yes for R2 single-end reads

htseq-count --stranded=reverse --format bam --order pos "${INPUT1}" "${REFGENOME}" > "${OUTPUT1}"

htseq-count --stranded=reverse --format bam "${INPUT2}" "${REFGENOME}" > "${OUTPUT2}"

htseq-count --stranded=yes --format bam "${INPUT3}" "${REFGENOME}" > "${OUTPUT3}"
