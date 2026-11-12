# Generate PLINK tasks for main or complement PRSs.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1]] %in% c("main", "complement")) {
  stop("Provide one branch argument: main or complement")
}
branch <- args[[1]]
project_root <- Sys.getenv("CSPRS_PROJECT_ROOT", unset = "/path/to/csPRS")
work_root <- file.path(project_root, "score_cs_beta")
gwas_dir <- file.path(work_root, if (branch == "main") "gwas_cs" else "gwas_cs_remain")
weight_dir <- file.path(work_root, if (branch == "main") "beta" else "beta_remain")

if (!dir.exists(weight_dir)) stop("Weight directory not found: ", weight_dir)
files <- list.files(weight_dir, pattern = "_weight\\.txt$")
if (!length(files)) stop("No group weight files found in: ", weight_dir)
groups <- sub("_weight\\.txt$", "", files)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript")
script_file <- sub("^--file=", "", script_arg)
# Some Rscript builds encode spaces as ~+~ in the --file argument.
if (!file.exists(script_file)) {
  script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
}
script_dir <- dirname(normalizePath(script_file, mustWork = TRUE))
score_script <- file.path(script_dir, "10_score_group.sh")
if (!file.exists(score_script)) stop("Scoring script not found: ", score_script)
jobs <- vapply(groups, function(group) {
  paste("bash", shQuote(score_script), shQuote(branch), shQuote(group))
}, character(1))
job_dir <- file.path(work_root, "jobs")
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
job_file <- file.path(job_dir, paste0("09_score_", branch, ".job"))
writeLines(jobs, job_file)
message("Generated ", length(jobs), " commands: ", job_file)
