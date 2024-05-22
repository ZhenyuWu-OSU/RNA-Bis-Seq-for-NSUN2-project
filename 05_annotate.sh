#!/usr/bin/bash
#SBATCH --job-name annoate
#SBATCH --account PAS1412
#SBATCH --time=10:00:00
#SBATCH --nodes=1 
#SBATCH --mem=64GB

export PATH="/fs/ess/PAS1412/Zy_transfer_PCON0022/2022ZY_BS_RNA/meRanTK-1.3.0:$PATH"
ml python/3.9-2022.05 
source activate merantk
cd /fs/ess/PAS1412/Zy_transfer_PCON0022/2022ZY_BS_RNA/meRanTK-1.3.0
chmod u+x meRanT
chmod u+x meRanGh
chmod u+x meRanCall
chmod u+x meRanCompare
chmod u+x meRanAnnotate

DATAPATH=/fs/ess/PAS1412/ZY_BS_seq2/30-707841147/00_fastq/trim

cd $DATAPATH

ls *.fastq |cut -d "_" -f 1| sort -u| while read id; do \
    echo ${id}_clean_1P.fastq; \
    cd /fs/ess/PAS1412/ZY_BS_seq2/30-707841147/4_annotation/; \
    meRanAnnotate \
         -p 8 \
         -o ./annotate_${id}.txt \
         -t /fs/ess/PAS1412/ZY_BS_seq2/30-707841147/3_call_meranGh/meRanCallResult_${id}_sorted.bam_FDR_0.01.txt  \
         -g /fs/ess/PAS1412/Zhenyu/10_human_KSHV_Dhfr_genome/human_KSHV_Dhfr.gff3 \
         -f 'gene' \
         -rd -er ; done
         



# cd /fs/ess/PAS1412/ZY_BS_seq2/30-707841147/4_annotation


# cd /fs/scratch/PAS1412/Zhenyu/4_5_merantK_compare_annoate

# meRanAnnotate \
# -p 8 \
# -t /fs/scratch/PAS1412/Zhenyu/3_call_meranGh/meRanCallResult_818T_sorted.bam_FDR_0.01.txt \
# -f 'gene' \
# -g /fs/scratch/PAS1412/Zhenyu/10_human_KSHV_Dhfr_genome/human_KSHV_Dhfr.gff3 \
# -o ./annoated_818_T.txt \
# -rd -er


# meRanAnnotate \
# -p 8 \
# -t /fs/scratch/PAS1412/Zhenyu/3_call_meranGh/meRanCallResult_RESUB-8-18TDP1_sorted_FDR_0.01.txt \
# -f 'gene' \
# -g /fs/scratch/PAS1412/Zhenyu/10_human_KSHV_Dhfr_genome/human_KSHV_Dhfr.gff3 \
# -o ./annoated_818TD.txt \
# -rd -er


# meRanAnnotate \
# -p 8 \
# -t /fs/scratch/PAS1412/Zhenyu/3_call_meranGh/meRanCallResult_8-1_TD_sorted.bam_FDR_0.01.txt \
# -f 'gene' \
# -g /fs/scratch/PAS1412/Zhenyu/10_human_KSHV_Dhfr_genome/human_KSHV_Dhfr.gff3 \
# -o ./annoated_81_TD.txt \
# -rd -er


# meRanAnnotate \
# -p 8 \
# -t /fs/scratch/PAS1412/Zhenyu/3_call_meranGh/meRanCallResult_RESUB-8-1_T_sorted.bam_FDR_0.01.txt \
# -f 'gene' \
# -g /fs/scratch/PAS1412/Zhenyu/10_human_KSHV_Dhfr_genome/human_KSHV_Dhfr.gff3 \
# -o ./annoated_81_T.txt \
# -rd -er


# meRanAnnotate \
# -p 8 \
# -t /fs/scratch/PAS1412/Zhenyu/3_call_meranGh/meRanCallResult_8-26-Bac16-1_sorted.bam_FDR_0.01.txt \
# -f 'gene' \
# -g /fs/scratch/PAS1412/Zhenyu/10_human_KSHV_Dhfr_genome/human_KSHV_Dhfr.gff3 \
# -o ./annoated_826Bac1.txt \
# -rd -er


# meRanAnnotate \
# -p 8 \
# -t /fs/scratch/PAS1412/Zhenyu/3_call_meranGh/meRanCallResult_8-26-Bac16-3_sorted.bam_FDR_0.01.txt \
# -f 'gene' \
# -g /fs/scratch/PAS1412/Zhenyu/10_human_KSHV_Dhfr_genome/human_KSHV_Dhfr.gff3 \
# -o ./annoated_826Bac3.txt \
# -rd -er


# meRanAnnotate \
# -p 8 \
# -t /fs/scratch/PAS1412/Zhenyu/3_call_meranGh/meRanCallResult_8-26-M-1_sorted.bam_FDR_0.01.txt \
# -f 'gene' \
# -g /fs/scratch/PAS1412/Zhenyu/10_human_KSHV_Dhfr_genome/human_KSHV_Dhfr.gff3 \
# -o ./annoated_M1.txt \
# -rd -er


# meRanAnnotate \
# -p 8 \
# -t /fs/scratch/PAS1412/Zhenyu/3_call_meranGh/meRanCallResult_8-26-M-3_sorted.bam_FDR_0.01.txt \
# -f 'gene' \
# -g /fs/scratch/PAS1412/Zhenyu/10_human_KSHV_Dhfr_genome/human_KSHV_Dhfr.gff3 \
# -o ./annoated_M3.txt \
# -rd -er



# meRanAnnotate \
# -p 8 \
# -t /fs/scratch/PAS1412/Zhenyu/3_call_meranGh/meRanCallResult_8-26-M-3_sorted.bam_FDR_0.01.txt \
# -f 'gene|mRNA|transcript|ncRNA' \
# -g /fs/scratch/PAS1412/Zhenyu/10_human_KSHV_Dhfr_genome/human_KSHV_Dhfr.gff3 \
# -o ./annoated_M3.txt \
# -rd -er

# meRanAnnotate \
# -p 8 \
# -t meRanCallResult_818T_sorted.bam_FDR_0.01.txt \
# -g /fs/ess/scratch/PCON0022/zhenyu/00_ref_genome_raw/gencode.v41.primary_assembly.annotation.gtf \
# -o /fs/ess/scratch/PCON0022/zhenyu/3_annotate_meranGh/meRanCallResult_818T_sorted2.txt \
# -rd -er

 