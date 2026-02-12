#!/bin/bash
#PBS -N deepvariant
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l storage=gdata/if89+gdata/te53+gdata/xl04
#PBS -l wd
#PBS -j oe 

module load singularity 

singularity exec \
  /g/data/if89/singularityimg/deepvariant_b350512297c12.sif \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type ONT_R104 \
  --ref /g/data/xl04/hrp561/adrna/data/SPI13459-Run1/Control_F_1_SPI13459A1_222V2FLT1_S1_L001_R1_001.fastq.gz \
  --reads /g/data/xl04/hrp561/adrna/analyses/alnbam/Control_F_1.bam \
  --output_vcf /g/data/xl04/eh8642/XY-walkthrough/output/Control_F_1_SPI13459A1_222V2FLT1_S1_L001_R1_001.output.vcf.gz \
  --output_gvcf /g/data/xl04/eh8642/XY-walkthrough/output/Control_F_1_SPI13459A1_222V2FLT1_S1_L001_R1_001.output.g.vcf.gz \
  --num_shards ${PBS_NCPUS} \
  --regions chrX chrY \
  --intermediate_results_dir /g/data/xl04/eh8642/XY-walkthrough/output/intermediate_results







