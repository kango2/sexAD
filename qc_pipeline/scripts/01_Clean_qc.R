# Cleanup txt files 
inputs <- commandArgs(trailingOnly = TRUE)
metadata <- inputs[1] 
workdir <- inputs[2]
keypair <- inputs[3]
data_cols <- inputs[4] |>
    strsplit(" ") |>
    unlist()

metadata_txt <- read.csv(metadata, header = TRUE, stringsAsFactors = FALSE)
individuals <- unique(metadata_txt$sample)
keypair <- read.csv(keypair, header = TRUE, stringsAsFactors = FALSE)

file_list <- list.files(path=paste0(workdir,"/01_Inputs"), 
                        pattern="\\.txt$", 
                        full.names=TRUE)

for(files in file_list){
  cur_file <- read.delim(
    files,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) 
  if (!"Sample" %in% colnames(cur_file)) {
    stop(paste("No Sample column found in", files))
  }
  cur_file <- cur_file[cur_file$Sample %in% individuals,]
  col_select <- c("Sample",colnames(cur_file)[which(colnames(cur_file) %in% data_cols)])
  cur_file <- cur_file[,col_select]
  
  cleaned_path <- paste0(files, ".cleaned")
  if(file.exists(cleaned_path)){
    next
  }else{
    write.table(cur_file, cleaned_path, row.names = FALSE, quote= FALSE, sep = "\t")
  }

}

write.table(individuals, paste0(workdir, "/01_Inputs/individuals.csv"), row.names = FALSE, quote= FALSE,col.names = FALSE)
