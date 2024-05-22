
setwd("~/Desktop/Bioinfo_analysis/2023.04.09 Public KD NSUN2 m5c data/")

## 433 first example 
library("DESeq2")
library(ggplot2)
library(stringr)

countmatrix=read.delim("./GSE133621_reads-count-cell-line.txt",header = T,row.names = 1)


info=colnames(countmatrix)
treatment=c("Control","KD_NSUN2","Control","KD_NSUN2")
str_link=stringr::str_c(info,treatment)
info2=data.frame(info,treatment,str_link)


colnames(info2)=c("SampleID","condition","extra")

######## GSE743333
dds <- DESeqDataSetFromMatrix(countData = countmatrix,
                              colData = info2,
                              design = ~condition)

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
# res <- results(dds)
# # res <- results(dds, name="condition_treated_vs_untreated")
# head(res)
# sum(res$padj < 0.05, na.rm=TRUE)

dds$condition
colnames(dds)
dds@colData
# dds$condition <- relevel(dds$condition, ref = "0hpi")
levels(dds$condition)
resultsNames(dds)
# res=results(dds)
KDNSUN2 <- results(dds,contrast=c("condition","KD_NSUN2","Control"))



res_symbol=as.data.frame(KDNSUN2)
normalized_counts=countmatrix[rownames(res_symbol),]

# res_symbol["GeneID"]=rownames(res_symbol)
res_symbol=cbind(res_symbol,normalized_counts)
# res_symbol$log2FoldChange=0-res_symbol$log2FoldChange
write.csv(res_symbol,"GSE133621_DEG_all.csv")

resp005=subset(res_symbol,res_symbol$padj<0.05)
head(resp005)
write.csv(test_pes005,"GSE74333_DEG0.05.csv")

write.csv(test_pes005,"GSE93750_DEG0.05.csv")

## annoate 
library(biomaRt)
trans_gene_id <- rownames(resp005)
ensembl <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
result_gene_entrez <- getBM(values = trans_gene_id,
                            filters = 'ensembl_gene_id', 
                            mart = ensembl,
                            attributes = c("ensembl_gene_id",'hgnc_symbol') )#'entrezgene_id'))

rownames(result_gene_entrez)=result_gene_entrez[,1]
test_pes005=merge.data.frame(resp005,result_gene_entrez,by=0,all.x = T)
head(test_pes005)
write.csv(test_pes005,"GSE93750_DEG0.05.csv")
