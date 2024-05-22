#!/usr/bin/bash
#SBATCH --job-name makeindex
#SBATCH --account PAS1412
#SBATCH --time=6:00:00
#SBATCH --nodes=1 
#SBATCH --mem=64GB

cd /fs/ess/PCON0022/zhenyuwu/2022ZY_BS_RNA
export PATH="/fs/ess/PCON0022/zhenyuwu/2022ZY_BS_RNA/meRanTK-1.3.0:$PATH"

ml python/3.9-2022.05 
source activate merantk
cd ./meRanTK-1.3.0
chmod u+x meRanGh

cd /fs/scratch/PAS1412/Zhenyu

meRanGh mkbsidx \
 -t 16 \
 -fa /fs/scratch/PAS1412/Zhenyu/10_human_KSHV_Dhfr_genome/human_KSHV_Dhfr.fa \
 -id /fs/scratch/PAS1412/Zhenyu/11_ref_meRanGh   