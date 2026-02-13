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