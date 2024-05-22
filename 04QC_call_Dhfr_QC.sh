#!/usr/bin/bash
#SBATCH --job-name call_Gh
#SBATCH --account PAS1412
#SBATCH --time=40:00:00
#SBATCH --nodes=1 
#SBATCH --mem=64GB

export PATH="/fs/ess/PAS1412/Zy_transfer_PCON0022/2022ZY_BS_RNA/meRanTK-1.3.0:$PATH"
ml python/3.9-2022.05 
source activate merantk
cd /fs/ess/PAS1412/Zy_transfer_PCON0022/2022ZY_BS_RNA/meRanTK-1.3.0
chmod u+x meRanGh
chmod u+x meRanGs
chmod u+x meRanCall
chmod u+x meRanCompare

DATAPATH=/fs/ess/PAS1412/ZY_BS_seq2/30-707841147/00_fastq/trim

cd $DATAPATH

ls *.fastq |cut -d "_" -f 1| sort -u| while read id; do \
    echo ${id}_Dhfr_QC ; \
    cd /fs/ess/PAS1412/ZY_BS_seq2/30-707841147/3_call_meranGh/; \
    meRanCall \
         -s /fs/ess/PAS1412/ZY_BS_seq2/30-707841147/2_align_meranGh/meRanGhResult/${id}_sorted.bam \
         -f /fs/ess/PAS1412/Zhenyu/10_human_KSHV_Dhfr_genome/human_KSHV_Dhfr.fa \
         -p 32 \
         -fs5 6 \
         -fs3 0 \
         -rs5 0 \
         -rs3 0 \
         -rl 150 \
         -ccr \
         -tref \
         --cSeqID "Dhfr" ; done


### call ended