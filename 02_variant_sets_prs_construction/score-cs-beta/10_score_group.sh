#!/bin/bash
# Usage: bash 10_score_group.sh main|complement GROUP
set -e
if [ "$#" -ne 2 ]; then
    echo "Usage: bash 10_score_group.sh main|complement GROUP" >&2
    exit 2
fi
branch=$1
cluster=$2
project_root="${CSPRS_PROJECT_ROOT:-/path/to/csPRS}"
work_root="${project_root}/score_cs_beta"
case "${branch}" in
    main)
        weight_dir="${work_root}/beta"
        out_path="${work_root}/score"
        ;;
    complement)
        weight_dir="${work_root}/beta_remain"
        out_path="${work_root}/score_remain"
        ;;
    *)
        echo "Unknown branch: ${branch}; use main or complement" >&2
        exit 2
        ;;
esac
weight="${weight_dir}/${cluster}_weight.txt"
if [ ! -f "${weight}" ]; then
    echo "Weight file not found: ${weight}" >&2
    exit 1
fi
mkdir -p "${out_path}"
module load PLINK/2_avx2_20221024
geno='/path/to/data/genotypes/ukb_genotypes'
plink2 --bfile "${geno}" \
    --score "${weight}" 1 2 \
    --score-col-nums 3,4 \
    --out "${out_path}/ukb_${cluster}"
