library(tidyverse)

##FileDate=20190503
##ENA_FILE_PATH=path to ENA file on ENA ftp site
##MD5=md5sum of file
##RUN_ID=SRA/ERA run accession
##STUDY_ID=SRA/ERA study accession
##STUDY_NAME=Name of study
##CENTER_NAME=Submission centre name
##SUBMISSION_ID=SRA/ERA submission accession
##SUBMISSION_DATE=Date sequence submitted, YYYY-MM-DD
##SAMPLE_ID=SRA/ERA sample accession
##SAMPLE_NAME=Sample name
##POPULATION=Sample population. Further information may be available with the data collection.
##EXPERIMENT_ID=Experiment accession
##INSTRUMENT_PLATFORM=Type of sequencing machine
##INSTRUMENT_MODEL=Model of sequencing machine
##LIBRARY_NAME=Library name
##RUN_NAME=Name of machine run
##INSERT_SIZE=Submitter specifed insert size/paired nominal length
##LIBRARY_LAYOUT=Library layout, this can be either PAIRED or SINGLE
##PAIRED_FASTQ=Name of mate pair file if exists (Runs with failed mates will have a library layout of PAIRED but no paired fastq file)
##READ_COUNT=Read count for the file
##BASE_COUNT=Basepair count for the file
##ANALYSIS_GROUP=Analysis group is used to identify groups, or sets, of data. Further information may be available with the data collection.

setwd("/g/data/xl04/eh8642")

sexData <- "data/igsr_samples.tsv"
oneT <- "data/1000G_2504_high_coverage.sequence.index.txt"
filereportRun <- "data/filereport_read_run_ERA1783081_tsv.txt"
filereportRead <- read_tsv(filereportRun, col_names = T)
sexdataRead <- read_tsv(sexData, col_names = T) %>%
  na.omit()
oneTGenomes <- read_tsv(oneT, col_names = T)
unique(sexdataRead$`Biosample ID`)


nContain <- c()
for(i in 1:length(sexdataRead$`Biosample ID`)){
  if(str_split(sexdataRead$`Biosample ID`[1], "")[[1]][[4]] == "E"){
    nContain <- c(nContain, i)
  }
}

eContain <- c()
for(i in 1:length(sexdataRead$`Biosample ID`)){
  if(str_split(sexdataRead$`Biosample ID`[1], "")[[1]][[4]] == "E"){
    nContain <- c(nContain, i)
  }
}

subpopSamples <- NULL


for(pops in unique(sexdataRead$`Population code`)){
  popSelect <- sexdataRead[sexdataRead$`Population code` == pops,]
  males <- popSelect[popSelect$Sex == "male",] %>%
    {. <- .[sample(1:nrow(.),1),];.}
  females <- popSelect[popSelect$Sex == "female",] %>%
    {. <- .[sample(1:nrow(.),1),];.}
  subpopSamples[[pops]] <- rbind(males, females)
  
}

sexdataReadJoin <- sexdataRead %>%
  {colnames(.)[1] <- "SAMPLE_NAME";.}

joinOne <- left_join(oneTGenomes, sexdataReadJoin, by="SAMPLE_NAME")

filereportReadJoin <- filereportRead %>%
  {colnames(.)[1] <- "RUN_ID";.}

FinalJoin <- left_join(joinOne,filereportReadJoin, by="RUN_ID")

cat(paste("Sex", "SubPop", "ENA file path", sep = " | "), "\n")
subpopSamplesJoin <- NULL
for(pops in unique(FinalJoin$`Population code`)){
  popSelect <- FinalJoin[FinalJoin$`Population code` == pops,]
  males <- popSelect[popSelect$Sex == "male",] %>%
    {. <- .[sample(1:nrow(.),1),];.}
  females <- popSelect[popSelect$Sex == "female",] %>%
    {. <- .[sample(1:nrow(.),1),];.}
  subpopSamplesJoin[[pops]] <- rbind(males, females)
  cat(paste("male", pops, males$`#ENA_FILE_PATH`, sep=" | "), "\n")
  cat(paste("female", pops, females$`#ENA_FILE_PATH`, sep=" | "),"\n")
}

subpopSamplesJoinDf <- do.call(rbind,subpopSamplesJoin) %>%
  as.data.frame()
write_csv(subpopSamplesJoinDf, "data/subpopMetadata.csv")

#saveRDS(subpopSamplesJoinDf, file = "data/subpopMetadata.rds")




















































