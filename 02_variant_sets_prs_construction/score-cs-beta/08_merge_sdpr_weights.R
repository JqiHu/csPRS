# Organize weights for selected SNPs
# From two sources:
## 1. SDPR_GW
## 2. SDPR_CS
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1]] %in% c("main", "complement")) {
  stop("Provide one branch argument: main or complement")
}
branch <- args[[1]]
project_root <- Sys.getenv("CSPRS_PROJECT_ROOT", unset = "/path/to/csPRS")
work_root <- file.path(project_root, "score_cs_beta")
gwas_dir <- file.path(work_root, if (branch == "main") "gwas_cs" else "gwas_cs_remain")
weight_dir <- file.path(work_root, if (branch == "main") "beta" else "beta_remain")

if (!dir.exists(file.path(weight_dir, "by_chr"))) {
  stop("SDPR chromosome outputs not found under: ", weight_dir)
}
setwd(weight_dir)
library(data.table)
library(stringr)

# ----------
# Type 1
# ----------
# GWAS summary statistics with position and weights from SDPR
weight <- fread('/path/to/resources/prs/beta/SDPR_cad_cardiogramc4a.txt',data.table=F)

# ----------
# Type 2
# ----------
groups <- list.files(gwas_dir, pattern = "\\.txt$")
groups <- sub("\\.txt$", "", groups)
if (!length(groups)) stop("No GWAS groups found in: ", gwas_dir)

### Read one by one and add SDPR_GW
for(g in groups){
  ## Type 2 by chr
  dat <- c()
  for(chr in 1:22){
    add <- fread(paste0('./by_chr/',g,'_chr',chr),data.table=F)
    dat <- rbind(dat,add)
  }
  tmp <- merge(dat,weight,by=c('SNP','A1'),all.x=T)
  colnames(tmp)[3:4] <- c('beta_cs','beta')
  print(sum(is.na(tmp$beta))) # need to be 0!!!!

  write.table(tmp,paste0(g,'_weight.txt'),
	      row.names=F,sep='\t',quote=F)
}
