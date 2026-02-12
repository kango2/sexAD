library(BiocManager)

rDNAOne <- readDNAStringSet(con)
con <- gzfile("/g/data/xl04/hrp561/adrna/data/SPI13459-Run1/Control_F_1_SPI13459A1_222V2FLT1_S1_L001_R1_001.fastq.gz", "rt")

fasta_lines <- readLines(con)



