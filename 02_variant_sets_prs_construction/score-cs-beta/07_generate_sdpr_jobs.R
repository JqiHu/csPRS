# Generate chromosome-specific SDPR jobs for main or complement PRSs.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1]] %in% c("main", "complement")) {
  stop("Provide one branch argument: main or complement")
}
branch <- args[[1]]
project_root <- Sys.getenv("CSPRS_PROJECT_ROOT", unset = "/path/to/csPRS")
work_root <- file.path(project_root, "score_cs_beta")
gwas_dir <- file.path(work_root, if (branch == "main") "gwas_cs" else "gwas_cs_remain")
weight_dir <- file.path(work_root, if (branch == "main") "beta" else "beta_remain")

job_dir <- file.path(work_root, "jobs")
output_dir <- file.path(weight_dir, "by_chr")
if (!dir.exists(gwas_dir)) stop("GWAS directory not found: ", gwas_dir)
gwas_files <- list.files(gwas_dir, pattern = "\\.txt$", full.names = TRUE)
if (!length(gwas_files)) stop("No GWAS input files found in: ", gwas_dir)
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sample_size <- 184305
sdpr_dir <- "/path/to/software/SDPR"
reference_dir <- "/path/to/resources/sdpr_reference/"
jobs <- character()
for (gwas_file in gwas_files) {
  group <- sub("\\.txt$", "", basename(gwas_file))
  for (chr in 1:22) {
    output <- file.path(output_dir, paste0(group, "_chr", chr))
    jobs <- c(jobs, paste(
      "cd", shQuote(sdpr_dir), "&& ./SDPR -mcmc -ref_dir", shQuote(reference_dir),
      "-ss", shQuote(gwas_file), "-N", sample_size, "-chr", chr,
      "-out", shQuote(output), "-n_threads 3"
    ))
  }
}
job_file <- file.path(job_dir, paste0("07_sdpr_", branch, ".job"))
writeLines(jobs, job_file)
message("Generated ", length(jobs), " commands: ", job_file)
