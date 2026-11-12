module load miniconda
conda activate ldsc_r

cd "/path/to/csPRS/heritability_par/annot_specificity_top10pct_merge251127" || exit 1
ldsc='/path/to/software/ldsc/ldsc.py'
ref_bfile='/path/to/software/ldsc/ref/1000G_EUR_Phase3_plink'

prefix=$1
chr=$2


python ${ldsc} \
        --l2 \
        --bfile ${ref_bfile}/1000G.EUR.QC.${chr} \
        --ld-wind-cm 1 \
        --annot "Annot_${prefix}_chr${chr}.annot.gz" \
        --out Annot_${prefix}_chr${chr}

