library("DESeq2")
library(ggplot2)
library(stringr)


##### 非常奇怪， 和论文上传的文件得到的结果差很多！！ 难道是因为我没有align EBV genome 然后这个组里面EBV 的转录组很多么？ 
####  我还是相信文章的作者，所以直接用他们的分析好了

getwd()
setwd("/Users/zhenyu/Desktop/Bioinfo_analysis/2024.02.19 EBV_RNA_SEQ/GSE240008/")

## 
## PRO-seq not sure 
GSE240008=read.delim("./out.txt",comment.char = "#",row.names = 1,check.names = F)
countmatrix=GSE240008[,-c(1:5)]
meta=read.csv("./SraRunTable.txt",)

colnames(countmatrix)
substr(colnames(countmatrix),47, 57)
colnames(countmatrix)=substr(colnames(countmatrix),47, 57)
write.csv(countmatrix,"GSE240008_count.csv")


tem=colnames(countmatrix)
tem=colnames(countmatrix)[order(colnames(countmatrix))]
tem
write.csv(tem,"name.csv",quote = F,row.names = F)
info=data.frame(meta$Run)
info=read.csv("./name.csv",header = F)
info$SampleID=meta$Run
info$condition=paste(meta$cell_line, meta$treatment, sep='_')
condition=paste(meta$cell_line, meta$treatment, sep='_')
info$V2=condition
info=info[info$V1 %in% colnames(countmatrix),]
countmatrix=countmatrix[,info$SampleID]


####

colnames(info)=c("SampleID","condition")


dds <- DESeqDataSetFromMatrix(countData = countmatrix,
                              colData = info,
                              design = ~ condition)

#### DEG analysis 
keep <- rowSums(counts(dds)) >=5
dds <- dds[keep,]
dds=DESeq(dds)

# View(counts(dds)) ## rawcounts
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts,"normalized_counts.csv")
vsd <- vst(dds, blind=FALSE)
vsd_matrix=assay(vsd)
# rld <- rlog(dds, blind=FALSE)
# rld_matrx=assay(rld)
write.csv(vsd_matrix,"vsd_matrix.csv")


pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_classic()
ggsave("PCA.png", device = "png")



dds$condition
colnames(dds)
dds@colData

DEG <- results(dds,contrast=c("condition","Akata_anti-IgG","Akata_(-)"))


DEG <- results(dds,contrast=c("condition","Mutu_anti-IgM","Mutu_(-)"))
DEG <- results(dds,contrast=c("condition","Mutu_24h Zta","Mutu_24h Ctl"))
DEG <- results(dds,contrast=c("condition","Mutu_12h Zta","Mutu_12h Ctl"))
DEG <- results(dds,contrast=c("condition","Mutu_6h Zta","Mutu_6h Ctl"))



res_symbol=as.data.frame(DEG)
res_symbol=cbind(res_symbol,normalized_counts)
resp005=subset(res_symbol,res_symbol$padj<0.05)

#
write.csv(resp005,"AKATAHDEG_005.csv")
write.csv(res_symbol,"AKATADEG_all.csv")


write.csv(resp005,"MUTUIGMHDEG_005.csv")
write.csv(res_symbol,"MUTUIGMHDEG_all.csv")


write.csv(resp005,"MUTU_12HDEG_005.csv")
write.csv(res_symbol,"MUTU_12HDEG_all.csv")


write.csv(resp005,"MUTU_6HDEG_005.csv")
write.csv(res_symbol,"MUTU_6HDEG_all.csv")

intersest=subset(res_symbol,rownames(res_symbol) %in% c("NOP2","NSUN2","MYC","TRIM25","TRIM65","TRIM66") )
write.csv(intersest,"MUTUIGMInterest_gene.csv")


write.csv(intersest,"AKATAInterest_gene.csv")
write.csv(intersest,"MUTU_12HInterest_gene.csv")
write.csv(intersest,"MUTU_6HInterest_gene.csv")

write.csv(intersest,"Interest_gene.csv")
write.csv(intersest,"Interest_gene.csv")
write.csv(intersest,"Interest_gene.csv")


### 


AKATA_Data=countmatrix[,15:20]
AKATA_info=info[info$SampleID %in% colnames(AKATA_Data),]


dds <- DESeqDataSetFromMatrix(countData = AKATA_Data,
                              colData = AKATA_info,
                              design = ~ condition)

#### DEG analysis 
keep <- rowSums(counts(dds)) >=5
dds <- dds[keep,]
dds=DESeq(dds)

# View(counts(dds)) ## rawcounts
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts,"normalized_counts.csv")
vsd <- vst(dds, blind=FALSE)
vsd_matrix=assay(vsd)
# rld <- rlog(dds, blind=FALSE)
# rld_matrx=assay(rld)
write.csv(vsd_matrix,"vsd_matrix.csv")


pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_classic()
ggsave("PCA.png", device = "png")







