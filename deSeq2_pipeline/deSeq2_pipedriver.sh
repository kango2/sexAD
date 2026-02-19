#!/bin/bash
#PBS -P xl04
#PBS -N deSeq2_pipedriver
#PBS -q normal
#PBS -l walltime=16:00:00
#PBS -l ncpus=1
#PBS -l mem=1GB
#PBS -l jobfs=200GB
#PBS -l storage=scratch/xl04+gdata/xl04+gdata/if89
#PBS -l wd
#PBS -o /g/data/xl04/eh8642/dart_Taro/09_pipelineOut/taropipeDriver.o
#PBS -e /g/data/xl04/eh8642/dart_Taro/09_pipelineOut/taroPipeDriver.e

# Taro pipeline - driver script 
export TMPDIR=/jobfs/$PBS_JOBID
mkdir -p "$TMPDIR"

# Read in input script 

while IFS="," read -r csv_location input_loc output_loc ref_loc databaseName
do 
    echo "Reading in-$csv_location"
    echo "Reading in-$input_loc"
    echo "Reading in-$output_loc"
    echo "Reading in $ref_loc"
    echo "Reading in $databaseName"

done < /g/data/xl04/eh8642/dart_Taro/taroPipeline/inputScript.csv


if [ ! -n $csv_location ]; then 
    echo "CSV location is not a string"
    set -e 
fi 

if [ ! -d $csv_location ]; then 
    echo "CSV location does not exist" 
    set -e 
fi 

# Manage conda installation

if [ ! -f "$csv_location/bin/activate" ]; then 
    echo "Activate file doesn't exist" 
    set -e 
else 
    source "$csv_location/bin/activate" 
    if [ $CONDA_DEFAULT_ENV == "base" ]; then 
        :
    else
        echo "Issue with conda activation" 
        set -e
    fi
fi 

# * other things to include (later) - check that r_env what other name is
# called is installed in conda environment etc... 

# Create directories
echo "Creating directories" 
if [ ! -d $output_loc ]; then 
    echo "Output directory doesn't exist" 
    set -e 
else 
    cd $output_loc
    mkdir 00_logs 01_Repmod 02_Repmask 03_Database 04_RepRemove 05_Filtered
fi
echo "done" 

# 1. Remove repeat regions 
# 1.1 RepeatModeler 

cd 01_Repmod

Repmodloc="${output_loc}/01_Repmod"

#RepModelJob=$(qsub -v databasename=$databaseName,refGenome=$ref_loc,01_Repmod_loc=$Repmodloc /g/data/xl04/eh8642/dart_Taro/taroPipeline/scripts/repeatmodeler_copy.sh)

RepModelJob=$(qsub \
  -v databasename="$databaseName",refGenome="$ref_loc",Repmodloc="$Repmodloc" \
  /g/data/xl04/eh8642/dart_Taro/taroPipeline/scripts/repeatmodeler_copy.sh \
  | tr -d '[:space:]'
)


# 1.2 RepeatMasker

#RepMaskerJob=$(qsub -N gen_validation -W depend=afterok:${RepModelJob})

Repmaskloc="${output_loc}/02_Repmask"

RepeatMaskJob=$(qsub \
  -W depend=afterok:$RepModelJob \
  -v databasename="$databaseName",refGenome="$ref_loc",Repmodloc="$Repmodloc",Repmaskloc="$Repmaskloc" \
  /g/data/xl04/eh8642/dart_Taro/taroPipeline/scripts/repeatmasker_copy.sh \
)

echo "RepeatModeler done" 


# 1.3 Remove repeats 
bed_file="${output_loc}/03_Database"
RepeatsRemove="${output_loc}/04_RepRemove"

RemoveRepeatsJob=$(qsub \
  -W depend=afterok:$RepeatMaskJob \
  -v bed_file="$bed_file",bam_files="$bam_files",RepeatsRemove="$RepeatsRemove" \
  /g/data/xl04/eh8642/dart_Taro/taroPipeline/scripts/04_remove_repeats_copy.sh \
)

echo "RepeatMasker done"

# 1.4 Filter low quality
FilterLQ="${output_loc}/05_Filtered"

FilterLowQ=$(qsub \
  -W depend=afterok:$RemoveRepeatsJob \
  -v indir="$RepeatsRemove",outdir="$FilterLQ" \
  /g/data/xl04/eh8642/dart_Taro/taroPipeline/scripts/05_filter_lowQ_copy.sh \
)

echo "Repeats removed"


































