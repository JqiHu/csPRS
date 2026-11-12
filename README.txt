csPRS analysis code

Overview
Scripts for cell-specific gene selection, cell-group clustering, construction
of coronary artery disease polygenic risk scores, and downstream association
and prediction analyses.

Script order
Two-digit prefixes indicate the order within each analysis module. Independent
analysis branches do not need to be run consecutively. The optional hg19
workflow is separate from the hg38 SNP-mapping workflow. Alternative interaction
models in candidates_pending_selection/ are not consecutive pipeline steps.

01_cell_gene_sets_clustering
    01_average_expression.R
        Average expression by cell type and state.
    02_select_specificity_genes.R
        Top-decile expression-specificity gene sets.
    03_calculate_cell_gene_overlap.R
        Pairwise gene-set overlap.
    04_cluster_cells_and_plot_jaccard.R
        Hierarchical clustering, variance explained, and heatmap.
    05_merge_cluster_gene_sets.R
        Merge genes using the curated ten-group mapping.
    06_calculate_cluster_gene_overlap.R
        Overlap between merged gene sets.
    07_pathway_enrichment.R
        Pathway enrichment.

02_variant_sets_prs_construction
    01_prepare_hg38_gene_lists.R
        hg38 gene coordinates and per-cell expressed-gene lists.
    02_map_snps_to_genes.R
        Coordinate conversion, SNP/gene mapping, per-cell weights.
    03_select_specificity_variant_weights.R
        Select weights for specificity genes.
    04_prepare_cell_gwas_and_sdpr_jobs.R
        Per-cell GWAS inputs and optional per-cell SDPR tasks.
    score-cs-beta/05_prepare_cluster_gwas.R
        Combine per-cell GWAS inputs by final group.
    score-cs-beta/06_prepare_complement_gwas.R
        Complementary SNP sets.
    score-cs-beta/07_generate_sdpr_jobs.R
        SDPR tasks for the selected PRS branch.
    score-cs-beta/08_merge_sdpr_weights.R
        Merge chromosome outputs and genome-wide weights.
    score-cs-beta/09_generate_score_jobs.R
        PLINK tasks for the selected PRS branch.
    score-cs-beta/10_score_group.sh
        Score one group with PLINK 2.
    11_aggregate_prs_scores.R
        Aggregate main and complement participant scores.
    optional/01_prepare_hg19_gene_coordinates.R
        Separate hg19 coordinate workflow; not a main-pipeline step.

03_analysis_data_preparation
    01_prepare_cohort_and_splits.R
        Cohort preparation and analysis splits.
    02_define_testing_prs_tertiles.R
        Testing-sample PRS tertiles.
    03_define_validation_prs_tertiles.R
        Validation-sample PRS tertiles.

04_association_interaction_validation
    01_overall_cad_associations.R
        Overall CAD associations.
    02_overall_prs_discrimination.R
        Overall PRS discrimination.
    03_testing_cad_associations.R
        Testing-sample CAD associations.
    04_testing_risk_factor_associations.R
        Risk-factor associations.
    05_testing_continuous_interactions.R
        Continuous PRS/lipid interactions.
    06_testing_tertile_binary_lipid_interactions.R
        PRS-tertile/binary-lipid interactions.
    07_testing_stratified_lipid_associations.R
        Testing-sample stratified associations.
    08_testing_absolute_risk.R
        Testing-sample absolute-risk summaries.
    09_validation_stratified_lipid_associations.R
        Validation-sample stratified associations.
    10_validation_absolute_risk.R
        Validation-sample absolute-risk summaries.
    candidates_pending_selection/01_testing_tertile_continuous_lipid_interactions.R
        Alternative testing interaction model.
    candidates_pending_selection/02_validation_continuous_interactions.R
        Alternative validation interaction model.
    candidates_pending_selection/03_validation_tertile_continuous_lipid_interactions.R
        Alternative validation interaction model.

05_other_reported_analyses
    01_prs_correlations.R
        PRS correlations.
    02_fit_ridge_model.R
        Ridge prediction model.
    03_fit_elastic_net_model.R
        Elastic-net prediction model.
    04_fit_xgboost_model.R
        XGBoost prediction model.
    05_compare_prediction_models.R
        Compare model predictions after fitting all three models.
    06_analyze_disease_embeddings.R
        Disease embeddings in case subgroups.
    07_compare_published_pdprs.R
        Comparison with published pdPRS.
    08_prepare_published_ecprs_variants.R
        Prepare published EC PRS variants.
    09_score_published_ecprs.sh
        Score the published EC PRS.
    10_compare_published_ecprs.R
        Comparison with published EC PRS.
    ldsc/01_prepare_cluster_annotations.R
        Ten-group SNP annotations.
    ldsc/02_calculate_ld_scores.sh
        Annotation-specific LD scores for one chromosome.
    ldsc/03_partition_heritability.sh
        Partitioned heritability with baseline annotations.
    ldsc/04_collect_heritability_results.R
        Extract cell-group results.

Input and output connections
The group-mapping workbook in 01_cell_gene_sets_clustering/ maps 77 annotations
to ten groups: aCM, vCM, Adip, Immune1, Immune2, EC, EC_Meso, FB, NC, and PC_SMC.
Place it in the configured merge_251127 working directory for the gene-set,
GWAS, and LDSC scripts that read it.

Gene overlap writes gene_list/shared_gene/specitificy_top10.csv. The clustering
script reads that file and uses the mapping workbook in
 gene_list/genes_specitificy_top10/merge_251127/.
The merged gene sets are written as top10_group_10.rds.

Module 02 uses hg38 gene_info.csv and the per-cell expressed-gene lists from
01_prepare_hg38_gene_lists.R. The optional hg19 script produces gene_info_2.csv
and does not replace these inputs. Upstream mapping and cell-type export sections
in 02_map_snps_to_genes.R are commented out; its active code requires an existing
snp_to_gene.rds and exports cell-state weights. Per-cell GWAS inputs feed grouped
GWAS preparation, followed by complement preparation and the branch-specific
training/scoring steps below. Both score branches are inputs to
11_aggregate_prs_scores.R.

After PRS aggregation, prepare the cohort/splits and tertiles in module 03 before
the corresponding module 04 analyses. In module 05, fit the three prediction
models before 05_compare_prediction_models.R. For the published EC PRS, run
08_prepare_published_ecprs_variants.R, 09_score_published_ecprs.sh, and then
10_compare_published_ecprs.R. Other module 05 analyses are separate branches.

Main and complement PRSs
Scripts 07, 08, and 09 in 02_variant_sets_prs_construction/score-cs-beta/ require
one explicit argument: main or complement. The same branch must be used across
training, weight merging, and scoring.

    Branch       GWAS input        Weights       Score output
    main         gwas_cs/          beta/         score/
    complement   gwas_cs_remain/   beta_remain/  score_remain/

These directories are under /path/to/csPRS/score_cs_beta/. The branch scripts
also accept CSPRS_PROJECT_ROOT as an environment override for /path/to/csPRS.
It should refer to the analysis workspace, not this code collection. Other
scripts retain explicit /path/to/... placeholders that must be configured.

Example, from the score-cs-beta code directory:
    export CSPRS_PROJECT_ROOT=/path/to/csPRS
    Rscript 07_generate_sdpr_jobs.R main
    Rscript 07_generate_sdpr_jobs.R complement

Execute the generated jobs/07_sdpr_main.job and jobs/07_sdpr_complement.job
commands using the local scheduler, and wait for all chromosome outputs.
Then run:
    Rscript 08_merge_sdpr_weights.R main
    Rscript 08_merge_sdpr_weights.R complement
    Rscript 09_generate_score_jobs.R main
    Rscript 09_generate_score_jobs.R complement

Execute jobs/09_score_main.job and jobs/09_score_complement.job to calculate
scores. The job generator records the absolute path of its sibling
10_score_group.sh; regenerate job files if the code directory is moved. Export
CSPRS_PROJECT_ROOT in the job environment when using a non-default workspace.
A single task can also be run as:
    bash 10_score_group.sh main Adip
    bash 10_score_group.sh complement Adip

LDSC
Within 05_other_reported_analyses/ldsc/, generate annotations, then run
02_calculate_ld_scores.sh merge251127 CHROMOSOME for chromosomes 1 through 22.
After those jobs complete, run 03_partition_heritability.sh and
04_collect_heritability_results.R. Annotation files and LD-score files share
heritability_par/annot_specificity_top10pct_merge251127/; heritability output
and collection use CAD_merge251127.results under heritability_par/res/.

Software and data
R dependencies are specified through library() and namespace calls, including
Seurat, SeuratDisk, Matrix, data.table, tidyverse packages, biomaRt,
Bioconductor annotation/enrichment packages, survival, glmnet, and xgboost.
External tools include SDPR, PLINK/PLINK 2, liftOver, and LDSC. Versions are not
pinned by a lock file. Adapt HPC module and scheduler commands locally.

Only the final group-mapping workbook is bundled. Other inputs include heart
single-cell expression data, GWAS summary statistics, genetic reference panels,
participant genotypes and phenotypes, and external score weights or embeddings.
Participant-level UK Biobank data require appropriate access.

Path configuration
Replace placeholders with local locations:
    /path/to/csPRS          analysis workspace
    /path/to/resources      shared inputs and reference resources
    /path/to/software       external tools
    /path/to/data/genotypes/ukb_genotypes
                           PLINK genotype-file prefix
Configure scheduler partition and email placeholders. The generic covariate
filename is ukbb_kinship_eur_PCs.csv. The numbered code directories do not mirror
the working data directories referenced by setwd() and relative paths. Prepare
those directories and inputs before execution.

Execution limitations
This collection is not an automated end-to-end pipeline. GWAS preprocessing
and genome-wide PRS training/scoring are not supplied as a complete workflow.
The cell gene-overlap script still uses table operations on list inputs; the
cluster gene-merge script requires attention when the excluded doublets label
is absent; SNP mapping and LDSC annotation currently deduplicate SNPs before
later gene filtering. These implementation limitations require review before
using the affected outputs. Numeric prefixes document order, not successful
validation of every scientific analysis.
