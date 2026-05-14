# Find the required metadata txt and move them 
inputs <- commandArgs(trailingOnly = TRUE)
keypair <- "/g/data/xl04/eh8642/RNAseq_run/qc_pipeline/inputs/metric_pair_no_spaces.csv"
data_cols <- c("samtools_stats-reads_mapped",
               "star-total_reads",
               "samtools_stats-reads_mapped_percent",
               "total_passed",
               "total_failed")
nextflow_out <- "/g/data/xl04/eh8642/RNAseq_run/SPI15112_run/Male/batch2"
cur_dir <- "/g/data/xl04/eh8642/RNAseq_run/qc_pipeline_practice"
metadata <- "/g/data/xl04/eh8642/RNAseq_run/samplesheets/SPI5112_samplesheet_malesBatch2.csv"

keypair <- read.csv(keypair, header = TRUE, stringsAsFactors = FALSE)
keypair <- keypair[match(data_cols, keypair$metric),]
metadata_txt <- read.csv(metadata, header = TRUE, stringsAsFactors = FALSE)
individuals <- unique(metadata_txt$sample)

for (i in 1:nrow(keypair)) {
    metric <- keypair$metric[i]
    filename <- keypair$file[i]

    # Find the file in the current directory
    file_path <- list.files(paste0(nextflow_out,"/workdirectory/multiqc/star_salmon/multiqc_report_data"),
                            pattern = filename, recursive = TRUE, full.names = TRUE)
    
    if (length(file_path) == 0) {
        warning(paste("File not found for metric:", metric))
    } else {
        # Move the file to the 01_Inputs directory
        dest_path <- paste0(cur_dir, "/01_Inputs/", filename)

        if(file.exists(dest_path)){ 
            next 
        }else{
            file.copy(file_path, dest_path, overwrite = TRUE)
            message(paste("Moved", filename, "to", dest_path))
        }
    }
}







list.files("/g/data/xl04/eh8642/RNAseq_run/SPI15112_run/Male/batch2/workdirectory/multiqc/star_salmon/multiqc_report_data", 
           pattern = "multiqc_general_stats.txt")









