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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
YQ="${SCRIPT_DIR}/bin/yq"
CONFIG_SH="${SCRIPT_DIR}/config.sh"
CONFIG_YAML="${SCRIPT_DIR}/config.yaml"
KEYPAIR="${SCRIPT_DIR}/inputs/metric_pair_no_spaces.csv"

# nf-core QC pipeline driver script
export TMPDIR=/jobfs/$PBS_JOBID
mkdir -p "$TMPDIR"

# Load config file
mapfile -t my_list < <(grep -E '^[a-zA-Z0-9_]+=' "$CONFIG_SH" | cut -d= -f1)

source "$CONFIG_SH"

# Load yaml file 
mapfile -t MULTIQC_ARGS < <($YQ -r '.multiqc_args[]' "$CONFIG_YAML")
mapfile -t PLOT_METRICS < <($YQ -r '.plot_metrics[]' "$CONFIG_YAML")
mapfile -t PLOT_BAR < <($YQ -r '.plot_panels.bar[]' "$CONFIG_YAML")
mapfile -t PLOT_HEATMAP < <($YQ -r '.plot_panels.heatmap[]' "$CONFIG_YAML")
PLOT_SCATTER=$($YQ -r '.plot_panels.scatter[] | [.x, .y, .label] | @tsv' "$CONFIG_YAML" | paste -sd'|' -)
mapfile -t METADATA_VARS < <($YQ -r '.metadata_vars[]' "$CONFIG_YAML")
mapfile -t METADATA_PLOT_METRICS < <($YQ -r '.metadata_plot_metrics[]' "$CONFIG_YAML")

WORKDIR=$($YQ -r '.workdir' "$CONFIG_YAML")
NEXTFLOW_OUT=$($YQ -r '.nextflow_out' "$CONFIG_YAML")
METADATA=$($YQ -r '.metadata' "$CONFIG_YAML")


# Create directories
echo "Creating directories" 
if [ ! -d $WORKDIR ]; then 
    mkdir -p "$WORKDIR"
fi
mkdir -p "$WORKDIR/01_Inputs" "$WORKDIR/02_Metadata" "$WORKDIR/03_Plots"
echo "done" 

if [ -f "$WORKDIR/metadata_input.csv" ]; then
    echo "Using existing run-local metadata_input.csv"
elif [ -f "$SCRIPT_DIR/metadata_input.csv" ]; then
    cp "$SCRIPT_DIR/metadata_input.csv" "$WORKDIR/metadata_input.csv"
    echo "Copied metadata_input.csv into the run folder"
else
    echo "metadata_input.csv not found in either $WORKDIR or $SCRIPT_DIR"
fi

cd "$SCRIPT_DIR"

# Find and move txt files required 

Rscript "$SCRIPT_DIR/scripts/00_Find_qc.R" "$KEYPAIR" "${MULTIQC_ARGS[*]}" "$NEXTFLOW_OUT" "$WORKDIR" "$METADATA"

# Cleanup txt files 

Rscript "$SCRIPT_DIR/scripts/01_Clean_qc.R" "$METADATA" "$WORKDIR" "$KEYPAIR" "${MULTIQC_ARGS[*]}"

# Extract_qc

Rscript "$SCRIPT_DIR/scripts/02_Extract_qc.R" "$WORKDIR" "${MULTIQC_ARGS[*]}"

# Plot_qc

Rscript "$SCRIPT_DIR/scripts/03_Plot_qc.R" "$WORKDIR" "${METADATA_VARS[*]}" "${METADATA_PLOT_METRICS[*]}" "${PLOT_METRICS[*]}" "${PLOT_BAR[*]}" "${PLOT_HEATMAP[*]}" "$PLOT_SCATTER"

# Create different qc tables for general stats and individuals












