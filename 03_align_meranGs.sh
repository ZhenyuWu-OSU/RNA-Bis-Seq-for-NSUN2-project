#!/usr/bin/bash
#SBATCH --job-name align
#SBATCH --account PAS1412
#SBATCH --time=80:00:00
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
    echo ${id}_clean_1P.fastq; \
    cd /fs/ess/PAS1412/ZY_BS_seq2/30-707841147/2_align_meranGs/; \
    meRanGs align \
         -o ./meRanGhResult \
         -f /fs/ess/PAS1412/ZY_BS_seq2/30-707841147/00_fastq/trim/${id}_clean_1P.fastq \
         -r /fs/ess/PAS1412/ZY_BS_seq2/30-707841147/00_fastq/trim/${id}_clean_2P.fastq \
         -t 12 \
         -S ${id}.sam \
         -ud ./meRanGhUnaligned \
         -un \
         -MM \
         -id  /fs/scratch/PAS1412/Zhenyu/11_ref_meRanGs \
         -star_outFilterMultimapNmax 20 \
         -sjO 149 \
         -GTF /fs/ess/PAS1412/Zhenyu/10_human_KSHV_Dhfr_genome/human_KSHV_Dhfr.gtf \
         -bg \
         -mbgc 10 \
         -mbp ; done
