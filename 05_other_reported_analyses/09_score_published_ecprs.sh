#!/bin/bash
#SBATCH --job-name ecPRS
#SBATCH -c 1 --mem-per-cpu 100g
#SBATCH -t 24:00:00
#SBATCH -p day,REPLACE_WITH_PARTITION,scavenge
#SBATCH --mail-type=ALL --mail-user=user@example.org
#SBATCH --requeue 

module load PLINK/1.9b_6.21-x86_64
cd /path/to/csPRS/score_split_analysis/analysis_compare/ec_prs
# Calculate PRS
geno='/path/to/data/genotypes/ukb_genotypes'
plink --bfile ${geno} \
	--score weight_ecPRS.txt 1 6 7 \
	--out ukb_ecprs

