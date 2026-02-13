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
