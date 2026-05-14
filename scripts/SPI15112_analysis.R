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
library(biomaRt)


.libPaths("/home/272/eh8642/R/x86_64-pc-linux-gnu-library/4.4")

#Female samples 
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


plotMA(SPI15112_female_analysis$results, ylim=c(-2,2))
plot_MA <- plotMA(SPI15112_female_analysis$results, ylim=c(-2,2))

plotDispEsts(SPI15112_female_analysis$dds)
install.packages("targets")


# Male samples 
SPI15112_bamM <- "/g/data/xl04/eh8642/RNAseq_run/SPI15112_run/Male/workdirectory/star_salmon"
SPI15112_quantPracM <- "/g/data/xl04/eh8642/RNAseq_run/SPI15112_run/Male/workdirectory/star_salmon"
SPI15112_contrastsM <- read.csv("/g/data/xl04/eh8642/RNAseq_run/metadata/SPI15112_contrasts.csv") %>%
  {. <- .[.$gender=="Male",];.}

SPI15112_tx2geneM <- "/g/data/xl04/eh8642/RNAseq_run/SPI15112_run/Male/workdirectory/star_salmon/tx2gene.tsv"


for(i in 1:length(SPI15112_contrastsM$group)){
  if(SPI15112_contrastsM[i,3] == "Familial" || SPI15112_contrastsM[i,3] =="Sporadic"){
    SPI15112_contrastsM[i,3] <- "Affected"
  }
}


SPI15112_male_readIn <- deSeq2Input(SPI15112_bamM, 
                                      contrastInput=SPI15112_contrastsM, 
                                      SPI15112_quantPracM, 
                                      SPI15112_tx2geneM, 
                                      factorKeep="ALL")

SPI15112_male_analysis <- deSeq2Analysis(SPI15112_male_readIn, 
                                           geneNames = T, 
                                           basicSummary = T, 
                                           varInterest = "group", 
                                           csvOut = "/g/data/xl04/eh8642/RNAseq_run/deSeq_results/", 
                                           nameOut = "SPI15112")



# Transcript mapping 
SPI15112_female_readIn$tx2gene
SPI15112_female_readIn$tx2gene
ensemblIDs <- SPI15112_female_readIn$tx2gene$transcript_id
ensemblIDs
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
ensembl
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
results <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand'),
                 filters='ensembl_gene_id',
                 values=c("ENSG00000139618"),
                 mart=ensembl)
library(GenomicFeatures)
BiocManager::install("txdbmaker")
txdb <- makeTxDbFromGFF(
  "/g/data/xl04/hrp561/adrna/reference/Homo_sapiens-GCA_009914755.4-2022_07-genes.chrnames.gtf"
)

# Gene coordinates
genes_gr <- genes(txdb)

genes_gr["ENSG05220038328"]

# Transcript coordinates
tx_gr <- transcripts(txdb)

tx_gr["ENST05220147413"]



















































