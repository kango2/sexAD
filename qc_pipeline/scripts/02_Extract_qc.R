# Read in metadata table and sort out mess 
inputs <- commandArgs(trailingOnly = TRUE)
workdir <- inputs[1] 
data_cols <- inputs[2] |>
    strsplit(" ") |>
    unlist()


file_list <- list.files(path=paste0(workdir,"/01_Inputs"), 
                        pattern="\\.txt\\.cleaned$", 
                        full.names=TRUE)

tables <- list()

for(files in file_list){
  cur_file <- read.delim(
    files,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) 
  col_select <- c("Sample",colnames(cur_file)[which(colnames(cur_file) %in% data_cols)])  
  tables[[length(tables) + 1]] <- cur_file[,col_select, drop = FALSE]
}

if (length(tables) == 0) {
  stop("No cleaned QC tables were found in 01_Inputs")
}

final_qc <- Reduce(function(x, y) merge(x, y, by = "Sample", all = TRUE), tables)

write.table(final_qc, paste0(workdir, "/02_Metadata/combined_qc.tsv"), 
    row.names = FALSE, 
    quote= FALSE, 
    sep = "\t")
