# Find the required metadata txt and move them 
inputs <- commandArgs(trailingOnly = TRUE)
keypair <- inputs[1]
data_cols <- inputs[2] |>
    strsplit(" ") |>
    unlist()
nextflow_out <- inputs[3]
cur_dir <- inputs[4]
metadata <- inputs[5] 
keypair <- read.csv(keypair, header = TRUE, stringsAsFactors = FALSE)
keypair <- keypair[match(data_cols, keypair$metric),]
metadata_txt <- read.csv(metadata, header = TRUE, stringsAsFactors = FALSE)
individuals <- unique(metadata_txt$sample)


for (i in 1:nrow(keypair)) {
    metric <- keypair$metric[i]
    filename <- keypair$file[i] 
    
    report_dir <- paste0(nextflow_out, "/workdirectory/multiqc/star_salmon/multiqc_report_data")
    file_candidates <- list.files(
        report_dir,
        recursive = TRUE,
        full.names = TRUE
    )
    file_path <- file_candidates[basename(file_candidates) == filename]
   
    if (length(file_path) == 0) {
        warning(paste("File not found for metric:", metric))
    } else if (length(file_path) > 1) {
        warning(paste("Multiple files found for metric:", metric, "- using the first match"))
        file_path <- file_path[1]
    } else {
        dest_path <- paste0(cur_dir, "/01_Inputs/", filename)
        if (!file.exists(dest_path)) {
            file.copy(file_path, dest_path, overwrite = TRUE)
            message(paste("Copied", filename, "to", dest_path))
        }
    }
}
