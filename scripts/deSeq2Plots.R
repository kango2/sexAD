deSeq2Plots <- function(dds,
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
