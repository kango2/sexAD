# Cleanup txt files 


inputs <- commandArgs(trailingOnly = TRUE)
metadata <- "/g/data/xl04/eh8642/RNAseq_run/samplesheets/SPI5112_samplesheet_malesBatch2.csv"
workdir <- "/g/data/xl04/eh8642/RNAseq_run/qc_pipeline_practice"
keypair <- "/g/data/xl04/eh8642/RNAseq_run/qc_pipeline/inputs/metric_pair_no_spaces.csv"
data_cols <- c("samtools_stats-reads_mapped",
               "star-total_reads",
               "samtools_stats-reads_mapped_percent",
               "total_passed",
               "total_failed")

metadata_txt <- read.csv(metadata, header = TRUE, stringsAsFactors = FALSE)
individuals <- unique(metadata_txt$sample)
keypair <- read.csv(keypair, header = TRUE, stringsAsFactors = FALSE)

file_list <- list.files(path=paste0(workdir,"/01_Inputs"), 
                        pattern="*.txt", 
                        full.names=TRUE)
for(files in file_list){
  cur_file <- read.delim(
    files,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) 
  cur_file <- cur_file[cur_file$Sample %in% individuals,]
  col_select <- c("Sample",colnames(cur_file)[which(colnames(cur_file) %in% data_cols)])
  
  cur_file <- cur_file[,col_select]
  file_base <- basename(files)
  metric <- keypair$metric[match(file_base, keypair$file)]
  
  
  write.table(cur_file, paste0(files,".cleaned"), row.names = FALSE, quote= FALSE, sep = "\t")
}

write.table(individuals, paste0(workdir, "/01_Inputs/individuals.csv"), row.names = FALSE, quote= FALSE,col.names = FALSE)

