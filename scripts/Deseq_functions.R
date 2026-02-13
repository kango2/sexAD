library(tidyverse)
library(tximport)
library(tidyverse)
library("vsn")
library("pheatmap")
library("apeglm")
library(org.Hs.eg.db)
library(AnnotationDbi)

make_fastq_sample_sheet <- function(filepaths,
                                    output_csv = "samples.csv",
                                    strandedness = "auto") {
  
  # Ensure character vector
  filepaths <- as.character(filepaths)
  
  # Keep only FASTQ files
  filepaths <- filepaths[grepl("\\.fastq\\.gz$", filepaths)]
  
  # Extract filename only
  filenames <- basename(filepaths)
  
  # Identify read direction
  read_type <- ifelse(grepl("_R1_", filenames), "R1",
                      ifelse(grepl("_R2_", filenames), "R2", NA))
  
  if (any(is.na(read_type))) {
    stop("Some files do not contain _R1_ or _R2_ in the name.")
  }
  
  # Extract sample + lane (everything before _R1_001.fastq.gz etc.)
  sample_name <- sub("(_R[12]_001\\.fastq\\.gz)$", "", filenames)
  
  # Create data frame
  df <- data.frame(
    sample = sample_name,
    filepath = filepaths,
    read = read_type,
    stringsAsFactors = FALSE
  )
  
  # Reshape to wide format (R1 and R2 columns)
  df_wide <- reshape(
    df,
    idvar = "sample",
    timevar = "read",
    direction = "wide"
  )
  
  # Clean column names
  colnames(df_wide) <- gsub("filepath\\.", "fastq_", colnames(df_wide))
  
  # Add strandedness column
  df_wide$strandedness <- strandedness
  
  # Reorder columns
  df_wide <- df_wide[, c("sample", "fastq_R1", "fastq_R2", "strandedness")]
  colnames(df_wide)[2:3] <- c("fastq_1", "fastq_2")
  
  # Write CSV
  write.csv(df_wide, output_csv, row.names = FALSE, quote = FALSE)
  
  return(df_wide)
}


deSeq2Input <- function(bamInput, 
                        contrastInput, 
                        quantsfInput, 
                        factorKeep="ALL", 
                        tx2geneLoc){
  
  # Error checks 
  if(!(typeof(bamInput) == "character" | type(bamInput) == "list")){
      stop("bamInput wrong format")
  }

  if(!(typeof(quantsfInput) == "character" | type(quantsfInput) == "list")){
      stop("bamInput wrong format")
  }

  finalOutputList <- NULL

  # Bam + contrast table 
  bamInputs <- c()
  contrastTable <- NULL
  if(typeof(bamInput) == "list"){
      for(locs in bamInput){
          readIn <- list.files(bamInput,pattern="*.bai") %>%
              {. <- gsub(".markdup.sorted.bam.bai","",.);.}
          bamInputs <- merge(bamInput,readIn)
      }
  }else{
      bamInputs <- list.files(bamInput,pattern="*.bai") %>%
              {. <- gsub(".markdup.sorted.bam.bai","",.);.}
  }
  
  if("ALL" %in% factorKeep){
    contrastTable <- as.data.frame(contrastInput) 
    for(i in 2:ncol(contrastTable)){
      contrastTable[,i] <- as.factor(contrastTable[,i])
    }
  }else{
    contrastTable <- as.data.frame(contrastInput[,1]) %>%
      {.[,factorKeep] <- contrastInput[,factorKeep];.} 
    for(i in 2:ncol(contrastTable)){
      contrastTable[,i] <- as.factor(contrastTable[,i])
    }
  }
  
  finalOutputList[["contrastTable"]] <- contrastTable

  # quantSF files and tx2gene list 
  quantInput <- c()
  if(type(quantsfInput) == "list"){
      for(locs in quantsfInput){
          readIn <- list.files(quantsfInput,pattern="*.bai", full.names = T) %>%
              {. <- gsub(".markdup.sorted.bam.bai","/quant.sf",.);.}
          quantInput <- merge(quantInput,readIn)
      }
  }else{
      quantInput <- list.files(quantsfInput,pattern="*.bai", full.names = T) %>%
              {. <- gsub(".markdup.sorted.bam.bai","/quant.sf",.);.}
  }

  finalOutputList[["quantInput"]] <- quantInput

  tx2gene <- NULL
  if(type(tx2geneLoc) == "list"){
      for(locs in tx2geneLoc){
          readIn <- read_tsv(locs,
              col_names = TRUE)
          tx2gene <- merge(tx2gene,readIn) %>% 
              unique()
      }
  }else{
      tx2gene <- read_tsv(tx2geneLoc,
              col_names = TRUE)
  }

  finalOutputList[["tx2gene"]] <- tx2gene

  return(finalOutputList)

}

bam_prac <- "/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out_run2/workdirectory/star_salmon"
quantPrac <- "/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out_run2/workdirectory/star_salmon"
Femalecontrast_2
tx2gene <- "/g/data/xl04/eh8642/RNAseq_run/Female_RNAseq_out_run2/workdirectory/star_salmon/tx2gene.tsv"

practice <- deSeq2Input(bam_prac, 
                        contrastInput=Femalecontrast_2, 
                        quantPrac, 
                        tx2gene, 
                        factorKeep="ALL")

deSeq2Analysis <- function(deSeq2Input, 
                           geneNames=FALSE, 
                           varInterest=NA,
                           basicSummary=TRUE, 
                           sigValue=0.1){
  # Error checks 
  if(!(type(deSeq2Input)=="list")){
      stop("deSeq2 input incorrect format")
  }else{
    #inputValues <- c("contrastTable","tx2gene","quantInput")
    #  if(!(inputValues %in% names(deSeq2Input))){
    #      stop("deSe2 input doesn't have all correct inputs")
    #  }
  }
  
  deSeq2AnalysisOutput<-NULL
  contrastTable <- deSeq2Input[["contrastTable"]]
  contrastTable$group <- as.factor(contrastTable$group)
  quantsfInput <- deSeq2Input[["quantInput"]]
  tx2gene <- deSeq2Input[["tx2gene"]]
  design_formula <- as.formula(paste("~", varInterest))
  
  if(is.na(varInterest)){
    varInterest <- colnames(contrastTable)[2]
  }
  
  Mergedtxi <- tximport(
      quantsfInput,
      type = "salmon",
      txOut = FALSE,
      tx2gene = tx2gene,
      countsFromAbundance = "no"  # <-- critical
  )   
  
  dds <- DESeqDataSetFromTximport(
      Mergedtxi,
      colData = contrastTable,
      design = design_formula
  )
 
  if(geneNames){
      gene_names <- mergedtx2gene |>
          dplyr::select(gene_id, gene_name) |>
          dplyr::distinct()
      rowData(dds)$gene_name <- gene_names$gene_name[
          match(rownames(dds), gene_names$gene_id)
      ]
      rowData(dds)$gene_name <- ifelse(
          is.na(rowData(dds)$gene_name),
          rownames(dds),
          rowData(dds)$gene_name
      )
      rowData(dds)$gene_name <- make.unique(rowData(dds)$gene_name)
  }

  smallestGroup <- 0
  for(levs in levels(contrastTable$group)){
      if(smallestGroup == 0){
          smallestGroup <- length(which(contrastTable$group %in% levs))
      }else if((length(which(contrastTable$group %in% levs)) < smallestGroup)){ 
          smallestGroup <- length(which(contrastTable$group %in% levs))
      }
  }

  keep <- rowSums(counts(dds) >= 10) >= smallestGroup
  dds <- dds[keep,]
  dds$group <- factor(dds$group, levels = levels(contrastTable$group))
    
  dds <- DESeq(dds)
  res <- results(dds)
  resNames<-resultsNames(dds)

  resLFC <- lfcShrink(dds,coef=resNames[2],type="apeglm")

  resOrdered <- res[order(res$pvalue),]
  resSig <- subset(resOrdered, padj < sigValue)
  
  deSeq2AnalysisOutput[["dds"]] <- dds
  
  if(geneNames){
    resOrdered <- as.data.frame(resOrdered) |>
      tibble::rownames_to_column("gene_id") |>
      dplyr::mutate(
        gene_name = rowData(dds)$gene_name[
          match(gene_id, rownames(dds))
        ]
      ) |>
      dplyr::relocate(gene_name, .before = gene_id)
    
    resSig <- subset(resOrdered, padj < sigValue)
    resSig <- as.data.frame(resSig) |>
      tibble::rownames_to_column("gene_id") |>
      dplyr::mutate(
        gene_name = rowData(dds)$gene_name[
          match(gene_id, rownames(dds))
        ]
      ) |>
      dplyr::relocate(gene_name, .before = gene_id)

  }
  
  deSeq2AnalysisOutput[["results"]] <- res
  deSeq2AnalysisOutput[["resOrdered"]] <- resOrdered
  deSeq2AnalysisOutput[["resSig"]] <- resSig
  
  return(deSeq2AnalysisOutput)
}

practice2 <- deSeq2Analysis(practice, 
                            geneNames = F, 
                            basicSummary = T, 
                            varInterest = "group")



deseq2Plots <- function(dds,
                        results,
                        intGroup=NA,
                        plotMA=T,
                        plotCounts=T, 
                        plotDisp=F, 
                        plotHeat=F){
  
  plots <- NULL
  
  if(plotMA){
    plots[["plotMA"]] <- plotMA(results, ylim=c(-2,2))
  }
  
  if(plotCounts){
    plots[["plotCounts"]] <- plotCounts(dds, 
                                        gene=which.min(results$padj), 
                                        intgroup=intGroup)
  }
  
  if(plotDisp){
    plots[["plotDisp"]] <- plotDispEsts(dds)
   }
  
  if(plotHeat){

    ntd <- normTransform(dds)
    #colnames(ntd) <- contrastmerged$file
    meanSdPlot(assay(ntd))
    select <- order(rowMeans(counts(dds,normalized=TRUE)),
                    decreasing=TRUE)[1:20]
    df <- as.data.frame(colData(dds)[,"group"])
    colnames(df) <- intGroup
    mat <- assay(vst(dds))[select, ]
    rownames(mat) <- rowData(dds)$gene_name[select]
    # ensure uniqueness
    rownames(mat) <- make.unique(rownames(mat))
    
    plots[["heatMap"]] <- pheatmap(mat, cluster_rows=FALSE, show_rownames=TRUE,
             cluster_cols=FALSE, annotation_col=df, show_colnames = TRUE)
    
  }
  
  return(plots)
  
}

plotPrac <- deseq2Plots(practice2$dds, 
                        practice2$results, 
                        intGroup="group")
plotPrac$plotCounts





asdf <- plotMA(practice2$results, ylim=c(-2,2))
class(asdf)












