#!/bin/bash
#SBATCH --partition=batch       # other queues for >7 days, more memory, gpu
#SBATCH --job-name=hisat2-unpaired
#SBATCH --ntasks=1
#SBATCH --time=24:00:00          # D-HH:MM:SS format
#SBATCH --mem=100G              # memory per node

#SBATCH --cpus-per-task=24
#SBATCH --mail-user=mjt86304@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/mjt86304/slurm_out/%j-%x.out
#SBATCH --error=/scratch/mjt86304/slurm_out/%j-%x.out

date

cd /scratch/mjt86304/trimm

mkdir hisat2

ml HISAT2/3n-20201216-gompi-2022a

for u in *_R1.trimmed.unpaired.fastq.gz
do
    n=${u%_R1.trimmed.unpaired.fastq.gz}
    v=${n}_R2.trimmed.unpaired.fastq.gz
    echo start_$n
    hisat2 -q --threads 24 -x /scratch/mjt86304/refgenome/gac_uga_v5_tran -U $u,$v -S hisat2/${n}_unpaired.sam
    echo end_$n
done