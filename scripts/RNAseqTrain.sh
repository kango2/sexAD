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

module load singularity 

OUTPUT_DIR="/g/data/xl04/eh8642/RNAseqTraining"
LOG_DIR="/g/data/xl04/eh8642/RNAseqTraining/log_dir"

REF="/g/data/xl04/hrp561/adrna/reference/chm13-t2t.ebv.phix.chrQ.xy.fa"
BAM_CONTROl1="/g/data/xl04/hrp561/adrna/analyses/alnbam/Control_F_1.bam"
BAM_CONTROL2="/g/data/xl04/hrp561/adrna/analyses/alnbam/Control_F_2.bam"
BAM_CONTROL3="/g/data/xl04/hrp561/adrna/analyses/alnbam/Control_F_3.bam"
TRUTH_VCF="/g/data/xl04/eh8642/RNAseqTraining"
TRUTH_BED="/g/data/xl04/eh8642/DeepVariantRun/Xregions.bed"

N_SHARDS=${PBS_NCPUS}


singularity exec \
    /g/data/if89/singularityimg/deepvariant_b350512297c12.sif \
    /opt/deepvariant/bin/make_examples \
      --mode=training \
      --ref "${REF}" \
      --reads "${BAM_CONTROl1}" \
      --examples "${OUTPUT_DIR}/training_set.with_label.tfrecord@${N_SHARDS}.gz" \
      --truth_variants "${TRUTH_VCF}" \
      --confident_regions "${TRUTH_BED}" \
      --task {} \
      --regions /g/data/xl04/eh8642/DeepVariantRun/Xregions.bed \
      --channel_list "BASE_CHANNELS,insert_size" \
  