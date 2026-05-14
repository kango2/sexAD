volcanoplot <- function(DESeq2Results, 
                        log2FoldUp=0.5, 
                        pvalSig=0.05){
  
  de <- DESeq2Results
  de$diffexpressed <- "NO"
  de$diffexpressed[de$log2FoldChange > log2FoldUp & de$pvalue < pvalSig] <- "UP"
  de$diffexpressed[de$log2FoldChange < -log2FoldUp & de$pvalue < pvalSig] <- "DOWN"
  
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("DOWN", "UP", "NO")
  
  de$gene_symbol <- rownames(de)
  de$delabel <- NA
  de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]
  
  plotOut <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
    geom_point() + 
    theme_minimal() +
    geom_text()
  
  return(plotOut)
  
}