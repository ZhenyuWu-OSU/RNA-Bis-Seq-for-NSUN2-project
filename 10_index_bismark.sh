#!/usr/bin/bash
#SBATCH --job-name makeindex
#SBATCH --account PAS1412
#SBATCH --time=8:00:00
#SBATCH --nodes=1 
#SBATCH --mem=64GB

cd /fs/ess/PCON0022/zhenyuwu/2022ZY_BS_RNA
export PATH="/fs/ess/PCON0022/zhenyuwu/2022ZY_BS_RNA/Bismark-master:$PATH"

ml python/3.9-2022.05 
ml bowtie2
ml hisat2
cd Bismark-master
chmod u+x bismark_genome_preparation

cd /fs/scratch/PAS1412/Zhenyu/

bismark_genome_preparation --hisat2 /fs/scratch/PAS1412/Zhenyu/11_bismark_HISAT


bismark_genome_preparation --bowtie2 /fs/scratch/PAS1412/Zhenyu/11_bismark_Bowtie
