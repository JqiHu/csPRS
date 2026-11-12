#!/bin/bash
#SBATCH --job-name LDSC_paritioned_heritability
#SBATCH -c 2 --mem-per-cpu 80g
#SBATCH -t 24:00:00
#SBATCH -p day,REPLACE_WITH_PARTITION,scavenge
#SBATCH --mail-type=ALL --mail-user=user@example.org
#SBATCH --requeue

module load miniconda
conda activate ldsc_r
cd /path/to/csPRS/heritability_par/res

# a1 is the effective allele
# function to calculate 
getH (){
  prefix=$1

  ldsc='/path/to/software/ldsc/ldsc.py'

  #  ───> note the trailing dot on each of these:
  gwas="/path/to/resources/prs/sumstats/munged_sumstats/munged_MTAG_cad_cardiogramc4a.txt.sumstats.gz"
#  annot_prefix="../annot_specificity_top10pct/${prefix}/Annot_${prefix}_chr"
#  annot_prefix="../annot_DEG/${prefix}/Annot_${prefix}_chr"
  annot_prefix="../annot_specificity_top10pct_merge251127/Annot_${prefix}_chr"
  baseline_prefix='/path/to/software/ldsc/ref/baseline/baselineLD.'
  weight_prefix='/path/to/software/ldsc/ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.'
  frq_prefix='/path/to/software/ldsc/ref/1000G_Phase3_frq/1000G.EUR.QC.'

  python ${ldsc} \
    --h2 ${gwas} \
    --ref-ld-chr ${annot_prefix}@,${baseline_prefix} \
    --w-ld-chr  ${weight_prefix} \
    --frqfile-chr ${frq_prefix} \
    --overlap-annot \
    --out "CAD_${prefix}"
}
export -f getH

getH merge251127
#module load parallel
#parallel getH {} ::: cell_state cell_type
