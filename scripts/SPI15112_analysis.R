library(BiocManager)
library(tidyverse)
library(tximport)
library(tidyverse)
library("vsn")
library("pheatmap")
library("apeglm")
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DESeq2)

# Female samples 
SPI15112_bamF <- "/g/data/xl04/eh8642/RNAseq_run/SPI15112_run/Female/workdirectory/star_salmon"
SPI15112_quantPracF <- "/g/data/xl04/eh8642/RNAseq_run/SPI15112_run/Female/workdirectory/star_salmon"
SPI15112_contrastsF <- read.csv("/g/data/xl04/eh8642/RNAseq_run/metadata/SPI15112_contrasts.csv") %>%
  {. <- .[.$gender=="Female",];.}

for(i in 1:length(SPI15112_contrastsF$group)){
  if(SPI15112_contrastsF[i,3] == "Familial" || SPI15112_contrastsF[i,3] =="Sporadic"){
    SPI15112_contrastsF[i,3] <- "Affected"
  }
}


SPI15112_tx2geneF <- "/g/data/xl04/eh8642/RNAseq_run/SPI15112_run/Female/workdirectory/star_salmon/tx2gene.tsv"

# Contrast sheet 
SPI15112_female_inputnames <- list.files("/g/data/xl04/eh8642/RNAseq_run/SPI15112_run/Female/workdirectory/star_salmon", pattern="*.bai") %>%
  {. <- gsub(".markdup.sorted.bam.bai","",.);.}


SPI15112_female_readIn <- deSeq2Input(SPI15112_bamF, 
                        contrastInput=SPI15112_contrastsF, 
                        SPI15112_quantPracF, 
                        SPI15112_tx2geneF, 
                        factorKeep="ALL")

SPI15112_female_analysis <- deSeq2Analysis(SPI15112_female_readIn, 
                            geneNames = T, 
                            basicSummary = T, 
                            varInterest = "group", 
                            csvOut = "/g/data/xl04/eh8642/RNAseq_run/deSeq_results/", 
                            nameOut = "SPI15112")







































































