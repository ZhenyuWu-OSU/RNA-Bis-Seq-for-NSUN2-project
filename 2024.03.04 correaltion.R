## TrY 1 https://malouche.github.io/DataVisuWithR/ggpubr.html
library(RTCGA)
installTCGA("RTCGA.mRNA")
library(RTCGA.mRNA)
load("data_mRNA.RData")

installTCGA("RTCGA.rnaseq.20160128") 
install_github("RTCGA/RTCGA", build_vignettes = TRUE)

?expressionsTCGA
library(ExperimentHub)
eh <- ExperimentHub()
myfiles <- query(eh,  "RTCGA.rnaseq.20160128" )
myfiles[[1]]  ## load the first resource in the list
?RTCGA::readTCGA() 

expr <- expressionsTCGA(DLBC.mRNA,extract.cols = c("NOP2","NSUN2","MYC","GATA1","E2F1","ELF1","SP1","SPI1"))


##  TRY2 https://waldronlab.io/cBioPortalData/articles/cBioPortalData.html
## https://www.karissawhiting.com/cbioportalR/articles/overview-of-workflow.html
library(cBioPortalData)
library(AnVIL)
cbio <- cBioPortal()
studies <- getStudies(cbio, buildReport = TRUE)
head(studies)                        

acc <- cBioPortalData(api = cbio, by = "hugoGeneSymbol", studyId = "dlbc_tcga",
                      genePanelId = "IMPACT341",
                     )
)


meta=metadata(acc)

## 
?cBioDataPack
cbio <- cBioPortal()
(getStudies(cbio)[["studyId"]])
laml <- cBioDataPack("dlbc_tcga", ask = FALSE,names.field = "Hugo_Symbol")
laml
laml@ExperimentList$mrna_seq_v2_rsem
dim(laml@ExperimentList$mrna_seq_v2_rsem_zscores_ref_all_samples)


data=(laml@ExperimentList$mrna_seq_v2_rsem_zscores_ref_all_samples)
### 
do.call(rbind.data.frame, lapply(x, function(y) {
  tmp <- getNodeCenter(y)
  tmp <- getPoints(tmp)
  tmp <- setNames(tmp, c("x", "y"))
  as.list(tmp)
}))


####  https://stackoverflow.com/questions/57317958/general-way-to-transform-s4-object-to-dataframe-r

S4_to_dataframe <- function(s4obj) {
  nms <- slotNames(s4obj)
  
  lst <- lapply(nms, function(nm) slot(s4obj, nm))
  as.data.frame(setNames(lst, nms))
}


data2=S4_to_dataframe(data)

data@assays@data@listData
data2=data.frame(data@assays@data@listData)
data3=as.data.frame(t(data2))
###

install.packages("ggstatsplot")
library(ggstatsplot)
ggscatterstats(
  data  = data3,
  x     = NSUN2,
  y     = MYC,
  xlab  = "NSUN2",
  ylab  = "C-MYC",
  title = "Correlation",
)

ggplot(data3, aes(x=NSUN2, y=MYC))  + geom_point()+geom_smooth(method=lm)
  
cor.test()
corrplot::corrplot()
test=data.frame(rownames(data2))



##

corrplot::corrplot(data3,)

?corrplot
ggplotsca

## http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/78-perfect-scatter-plots-with-correlation-and-marginal-histograms/

## https://cran.r-project.org/web/packages/TCGAretriever/vignettes/getting_started_with_TCGAretriever.html

library(ggplot2)
library(ggpubr)

expr <- expressionsTCGA(DLBC.mRNA,extract.cols = c("NOP2","NSUN2","MYC","GATA1","E2F1","ELF1","SP1","SPI1"))


ggscatter(data3, x = "MYC", y = c("NOP2","NSUN2"),
          combine = TRUE, ylab = "Expression",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE, 
          repel = T,
          #cor.coef = TRUE,# Add confidence interval
          xlab = "Expression of C-MYC",
          add.params = list(color = "blue",
                            fill = "lightgray")
)+ stat_cor(method = "spearman")  # Add correlation coefficient 



ggscatter(data3, x = "SPI1", y = c("NOP2","NSUN2"),
          combine = TRUE, ylab = "Expression",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE, 
          repel = T,
          #cor.coef = TRUE,# Add confidence interval
          xlab = "Expression of SPI1",
          add.params = list(color = "blue",
                            fill = "lightgray")
)+ stat_cor(method = "spearman")  # Add correlation coefficient 



?stat_cor




##
cor(xRle,yRle)
library(survival)
library(survminer)

