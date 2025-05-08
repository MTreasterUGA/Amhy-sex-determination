#!/bin/bash
#SBATCH --partition=batch       # other queues for >7 days, more memory, gpu
#SBATCH --job-name=trimmomatic
#SBATCH --ntasks=1
#SBATCH --time=48:00:00         # D-HH:MM:SS format
#SBATCH --mem=15G               # memory per node

#SBATCH --cpus-per-task=4
#SBATCH --mail-user=mjt86304@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/mjt86304/slurm_out/%j-%x.out
#SBATCH --error=/scratch/mjt86304/slurm_out/%j-%x.out

date

cd /scratch/mjt86304/embryoseq

mkdir trimm

ml Trimmomatic/0.39-Java-13

for f in *1_001.fastq.gz 
do
    g=${f%1_001.fastq.gz} 
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 4 $f ${g}2_001.fastq.gz trimm/${g}1.trimmed.paired.fastq.gz trimm/${g}1.trimmed.unpaired.fastq.gz trimm/${g}2.trimmed.paired.fastq.gz trimm/${g}2.trimmed.unpaired.fastq.gz ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done