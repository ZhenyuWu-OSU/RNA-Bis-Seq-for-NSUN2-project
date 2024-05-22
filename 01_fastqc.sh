#!/bin/bash
#SBATCH --nodes=1 --ntasks=2
#SBATCH --mem=16GB
#SBATCH --time=4:00:00
#SBATCH --job-name=Prealignment_QC
#SBATCH --mail-type=ALL
#SBATCH --account=PAS1412

cd /fs/ess/PAS1412/ZY_BS_seq2/30-707841147/    # go to current working directory

module load fastqc

mkdir fastqc_out
fastqc  -o  fastqc_out /fs/ess/PAS1412/ZY_BS_seq2/30-707841147/00_fastq/*.fastq.gz
