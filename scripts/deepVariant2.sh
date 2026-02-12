#!/bin/bash
#PBS -N deepvariant
#PBS -P xl04
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l storage=gdata/if89+gdata/xl04
#PBS -l wd
#PBS -j oe 
#PBS -l jobfs=50GB

mkdir -p output_3
mkdir -p output_3/intermediate_results_dir

module load singularity 


singularity exec \
    /g/data/if89/singularityimg/deepvariant_b350512297c12.sif \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type WES \
    --customized_model=model/model.ckpt \
    --ref /g/data/xl04/hrp561/adrna/reference/chm13-t2t.ebv.phix.chrQ.xy.fa \
    --reads /g/data/xl04/hrp561/adrna/analyses/alnbam/Control_F_3.bam \
    --output_vcf /g/data/xl04/eh8642/DeepVariantRun/output_3/Control_F_3.output.vcf.gz \
    --output_gvcf /g/data/xl04/eh8642/DeepVariantRun/output_3/Control_F_3.output.g.vcf.gz \
    --num_shards ${PBS_NCPUS} \
    --regions /g/data/xl04/eh8642/DeepVariantRun/Xregions.bed \
    --intermediate_results_dir /g/data/xl04/eh8642/DeepVariantRun/output_3/intermediate_results_dir
