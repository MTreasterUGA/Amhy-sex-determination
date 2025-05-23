#!/bin/bash
#SBATCH --partition=batch       # other queues for >7 days, more memory, gpu
#SBATCH --job-name=hisat2-build
#SBATCH --ntasks=1
#SBATCH --time=8:00:00          # D-HH:MM:SS format
#SBATCH --mem=200G              # memory per node

#SBATCH --cpus-per-task=4
#SBATCH --mail-user=mjt86304@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/mjt86304/slurm_out/%j-%x.out
#SBATCH --error=/scratch/mjt86304/slurm_out/%j-%x.out

date

cd /scratch/mjt86304/

ml HISAT2/3n-20201216-gompi-2022a

hisat2_extract_splice_sites.py gac_uga_v5.gtf > gac_uga_v5.ss
hisat2_extract_exons.py gac_uga_v5.gtf > gac_uga_v5.exon
hisat2-build --threads 16 --exon gac_uga_v5.exon --ss gac_uga_v5.ss gac_uga_v5.fna gac_uga_v5_tran