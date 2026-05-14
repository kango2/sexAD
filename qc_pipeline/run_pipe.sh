#!/bin/bash
#PBS -P xl04
#PBS -N taroPipeDriver
#PBS -q normal
#PBS -l walltime=16:00:00
#PBS -l ncpus=1
#PBS -l mem=1GB
#PBS -l jobfs=200GB
#PBS -l storage=scratch/xl04+gdata/xl04+gdata/if89
#PBS -l wd
#PBS -W umask=022
#PBS -M jose.mijangosaraujo@sydney.edu.au
#PBS -o /g/data/xl04/eh8642/dart_Taro/09_pipelineOut/taropipeDriver.o
#PBS -e /g/data/xl04/eh8642/dart_Taro/09_pipelineOut/taroPipeDriver.e

module purge 
module load samtools 
module load singularity
module load R 

cur_dir=$(pwd)
YQ="${cur_dir}/bin/yq"
KEYPAIR=inputs/metric_pair_no_spaces.csv 

# nf-core QC pipeline driver script
export TMPDIR=/jobfs/$PBS_JOBID
mkdir -p "$TMPDIR"

# Load config file
mapfile -t my_list < <(grep -E '^[a-zA-Z0-9_]+=' config.sh | cut -d= -f1)

source config.sh

# Load yaml file 
mapfile -t MULTIQC_ARGS < <($YQ -r '.multiqc_args[]' config.yaml)
mapfile -t PLOT_METRICS < <($YQ -r '.plot_metrics[]' config.yaml)
mapfile -t PLOT_BAR < <($YQ -r '.plot_panels.bar[]' config.yaml)
mapfile -t PLOT_HEATMAP < <($YQ -r '.plot_panels.heatmap[]' config.yaml)
PLOT_SCATTER=$($YQ -r '.plot_panels.scatter[] | [.x, .y, .label] | @tsv' config.yaml | paste -sd'|' -)

WORKDIR=$($YQ -r '.workdir' config.yaml)
NEXTFLOW_OUT=$($YQ -r '.nextflow_out' config.yaml)
METADATA=$($YQ -r '.metadata' config.yaml)


# Create directories
echo "Creating directories" 
if [ ! -d $WORKDIR ]; then 
    mkdir -p "$WORKDIR"
fi
mkdir -p "$WORKDIR/01_Inputs" "$WORKDIR/02_Metadata" "$WORKDIR/03_MultiQC" "$WORKDIR/04_Plots"
echo "done" 

cd $cur_dir

# Find and move txt files required 

Rscript scripts/00_Find_qc.R "$KEYPAIR" "${MULTIQC_ARGS[*]}" "$NEXTFLOW_OUT" "$WORKDIR" "$METADATA"

# Cleanup txt files 

Rscript scripts/01_Clean_qc.R "$METADATA" "$WORKDIR" "$KEYPAIR" "${MULTIQC_ARGS[*]}"

# Extract_qc

Rscript scripts/02_Extract_qc.R "$WORKDIR" "${MULTIQC_ARGS[*]}"

# Plot_qc

Rscript scripts/03_Plot_qc.R "$WORKDIR" "${PLOT_METRICS[*]}" "${PLOT_BAR[*]}" "${PLOT_HEATMAP[*]}" "$PLOT_SCATTER"

# Create different qc tables for general stats and individuals


















