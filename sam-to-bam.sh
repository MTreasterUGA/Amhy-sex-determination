#!/bin/bash
#SBATCH --partition=batch       # other queues for >7 days, more memory, gpu
#SBATCH --job-name=sam-to-bam
#SBATCH --ntasks=2
#SBATCH --time=24:00:00         # D-HH:MM:SS format
#SBATCH --mem=10G               # memory per node

#SBATCH --cpus-per-task=4
#SBATCH --mail-user=mjt86304@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/mjt86304/slurm_out/%j-%x.out
#SBATCH --error=/scratch/mjt86304/slurm_out/%j-%x.out

date

cd /scratch/mjt86304/trimm/hisat2

ml SAMtools/1.17-GCC-12.2.0

for f in *
do
    samtools view -bS --threads 4 $f | samtools sort --threads 4 -o ${f%.sam}_sorted.bam
done

for g in *.bam
do
    samtools index -b -@ 4 $g
done