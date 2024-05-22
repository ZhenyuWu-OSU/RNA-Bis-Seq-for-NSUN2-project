
setwd("~/Desktop/Bioinfo_analysis/2023.04.09 Public KD NSUN2 m5c data/")

GSE93750=readxl::read_xls("./NSUN2 as the methyltransferase and ALYREF /GSE93750_all-samples-reads-count.xls")
GSE93750_rawdata=GSE93750[,1:5]
write.csv(GSE93750_rawdata,"GSE93750_rawdata.csv")
GSE74333_c1=read.delim("./unpublished—GSE74333_RAW/GSM1917422_siCTRL_HeLa_rep1_union.txt",row.names = 1,header = F)
GSE74333_c2=read.delim("./unpublished—GSE74333_RAW/GSM1917423_siCTRL_HeLa_rep2_union.txt",row.names = 1,header = F)
GSE74333_KD1=read.delim("./unpublished—GSE74333_RAW/GSM1917424_siNSUN2_HeLa_rep1_union.txt",row.names = 1,header = F)
GSE74333_KD2=read.delim("./unpublished—GSE74333_RAW/GSM1917425_siNSUN2_HeLa_rep2_union.txt",row.names = 1,header = F)

GSE74333_raw=cbind(GSE74333_c1,GSE74333_c2,GSE74333_KD1,GSE74333_KD2)
colnames(GSE74333_raw)=c("c1","c2","KD1","KD2")
write.csv(GSE74333_raw,"GSE74333_raw.csv")



## 433 first example 
library("DESeq2")
library(ggplot2)
library(stringr)

info=meta[,c("Run","source_name")]
countmatrix=countmatrix[,info$Run]
colnames(info)=c("SampleID","condition")

CELLINE=c(rep("BC3",4),rep("BCBL1",4),rep("HEL6",4))
treatment=rep(c("Control","PEP005","JQ1","JQ1+PEP005"),3)
str_link=str_c(CELLINE,treatment)
info2=data.frame(info$SampleID,CELLINE,treatment,str_link)
#### 
info=colnames(GSE74333_raw)
treatment=c("Control","Control","KD_NSUN2","KD_NSUN2")
str_link=stringr::str_c(info,treatment)
info2=data.frame(info,treatment,str_link)

countmatrix=GSE74333_raw
######## GSE743333
dds <- DESeqDataSetFromMatrix(countData = countmatrix,
                              colData = info2,
                              design = ~treatment)


####### GSE93750_rawdata

countmatrix=GSE93750_rawdata[,2:5]
rownames(GSE93750_rawdata)=GSE93750_rawdata[,1]
GSE93750_rawdata=as.data.frame(GSE93750_rawdata)
rownames(countmatrix)=GSE93750_rawdata[,1]


info=colnames(countmatrix)
treatment=c("Control","Control","KD_NSUN2","KD_NSUN2")
str_link=stringr::str_c(info,treatment)
info2=data.frame(info,treatment,str_link)

dds <- DESeqDataSetFromMatrix(countData = countmatrix,
                              colData = info2,
                              design = ~treatment)


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


pcaData <- plotPCA(vsd, intgroup=c("treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=treatment)) +
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
KDNSUN2 <- results(dds,contrast=c("treatment","KD_NSUN2","Control"))
KDNSUN2



res_symbol=as.data.frame(KDNSUN2)
normalized_counts=countmatrix[rownames(res_symbol),]

# res_symbol["GeneID"]=rownames(res_symbol)
res_symbol=cbind(res_symbol,normalized_counts)
# res_symbol$log2FoldChange=0-res_symbol$log2FoldChange
write.csv(res_symbol,"res96_0_DEG_all.csv")

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


##
GSE93750_DEG0.05=read.csv("./GSE93750_DEG0.05.csv")
GSE74333_DEG0.05=read.csv("./GSE74333_DEG0.05.csv")



