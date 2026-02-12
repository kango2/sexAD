#!/bin/bash
#PBS -l ncpus=1
#PBS -l mem=32GB
#PBS -l jobfs=32GB
#PBS -q normal
#PBS -lother=mdss
#PBS -P xl04
#PBS -l walltime=10:00:00
#PBS -l storage=gdata/xl04+gdata/if89
#PBS -l wd

module purge
module load nextflow 
module load singularity

nextflow run nf-core/rnaseq \
    -profile singularity,nci_gadi \
    --input samplesheet.csv \
    --outdir RNAseq_out \
    --genome CHM13

















