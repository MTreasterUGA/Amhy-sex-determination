#!/bin/bash
#SBATCH --partition=batch       # other queues for >7 days, more memory, gpu
#SBATCH --job-name=htseq-parallel
#SBATCH --ntasks=1
#SBATCH --time=72:00:00         # D-HH:MM:SS format
#SBATCH --mem=20G              # memory per node

#SBATCH --cpus-per-task=1
#SBATCH --mail-user=mjt86304@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/mjt86304/slurm_out/%j-%x.out
#SBATCH --error=/scratch/mjt86304/slurm_out/%j-%x.out

date

cd /scratch/mjt86304/hisat2

mkdir /scratch/mjt86304/htseq

ml HTSeq/2.0.2-foss-2022a

for f in *sorted.bam
do
    htseq-count --stranded=no --format bam $f -r pos /scratch/mjt86304/refgenome/gac_uga_v5.gtf > /scratch/mjt86304/htseq/${f}.txt
done