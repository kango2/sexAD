pairwisePlot <- function(ds2results,
                         yAxis,
                         annotation.gtf,
                         specChr){
  
  # Get gene co-ordinates 
  txdb <- makeTxDbFromGFF(annotation.gtf) %>%
    genes() 
  
  txdb <- txdb[seqnames(genes_gr) %in% specChr]
  
  
  overLap <- rownames(ds2results) %>%
    {. <- .[which(. %in% txdb$gene_id)];.}
  
  txdb <- txdb[txdb$gene_id%in%overLap]
  
  rangesSelect <- ranges(txdb) %>%
    as.data.frame() %>%
    {.$mean <- rep(0, nrow(.));.}
  
  for(i in 1:nrow(rangesSelect)){
    rangesSelect[i,5] <- (rangesSelect[i,1]+rangesSelect[i,2])/2
  }
  
  resultsMapping <- ds2results %>%
    {. <- .[rownames(.)%in%overLap,]} %>%
    as.data.frame() %>%
    {.$location <- rangesSelect$mean;.} %>%
    {.$location <- .$location/1000000;.} %>%
    {.$sig <- rep("UnSig",nrow(.));.}
  
  # Plot 
  resultsMapPlot <- ggplot(data = resultsMapping,mapping = aes(location,resultsMapping[,yAxis]))+
    geom_point()+
    ylim(min(resultsMapping[,yAxis]),max(resultsMapping[,yAxis])) + 
    labs(
      x = "X-chrom Location (Mbp)",
      y = paste0("F/M ratio ",yAxis),
      color = "Legend Title" # Changes the title for the 'color' aesthetic
    )
  
  
  return(resultsMapPlot)
  
}