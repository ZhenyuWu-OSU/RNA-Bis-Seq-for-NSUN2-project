#!/bin/bash
#SBATCH --nodes=1 --ntasks=2
#SBATCH --mem=32GB
#SBATCH --time=20:00:00
#SBATCH --job-name=trim
#SBATCH --account=PAS1412

cd /fs/ess/PAS1412/ZY_BS_seq2/30-707841147/00_fastq    # go to current working directory

module load java
module load trimmomatic

mkdir trim

ls *.fastq.gz |cut -d "_" -f 1| sort -u| while read id; do \
    echo ${id}_R1_001.fastq; \
    java -jar /users/PAS1475/zhenyuwu/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
        -threads 32 ${id}_R1_001.fastq.gz ${id}_R2_001.fastq.gz \
        -baseout ./trim/${id}_clean.fastq \
        ILLUMINACLIP:/users/PAS1475/zhenyuwu/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:8 MINLEN:36; done
