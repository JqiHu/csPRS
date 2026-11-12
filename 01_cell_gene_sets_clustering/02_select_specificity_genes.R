# Generate gene list based on expression data downloaded 
# !!! R 4.4
# require large memory to read the data (>30G)
setwd('/path/to/csPRS/gene_list/genes_specitificy_top10')
library(Seurat)
library(SeuratDisk)
library(Matrix)

# 1) Load the Seurat object made from freeze.h5ad
#Convert("hca_heart_global_ctl200723_freeze.h5ad", dest = "h5seurat", overwrite = TRUE)
obj <- LoadH5Seurat("../data/hca_heart_global_ctl200723_freeze.h5seurat")
DefaultAssay(obj) <- "RNA"
#### It seems that the data has not been normalized

# 2) Pick the identity as your cell-type column
idcol <- 'cell_type'
#idcol <- 'cell_states' # more refined definition of cell subpopulations
Idents(obj) <- obj[[idcol, drop = TRUE]]

# 3) Build a linear expression matrix (CP10k) from raw counts
cts <- GetAssayData(obj, assay = "RNA", slot = "counts")  # genes x cells (dgCMatrix)
libsize <- Matrix::colSums(cts)

# Keep only cells with valid identity and nonzero library size
ct_vec <- as.character(Idents(obj))
keep_cells <- !is.na(ct_vec) & libsize > 0
cts      <- cts[, keep_cells, drop = FALSE]
libsize  <- libsize[keep_cells]
ct_vec   <- ct_vec[keep_cells]

# Normalize counts per cell to CP10k (linear)
cp10k <- t( t(cts) / libsize ) * 1e4   # genes x cells (dgCMatrix)

# QC: keep genes expressed in >1% of cells (any type)
#    "expressed" = counts > 0 in a cell
detected_prop <- Matrix::rowMeans(cts > 0)  # proportion of cells with nonzero counts
gene_keep <- detected_prop > 0.01           # >1% 
sum(gene_keep==F) # remove 19,334 genes
cp10k <- cp10k[gene_keep, , drop = FALSE]   # filter early for speed

# 4) Average expression per gene in each cell type (on linear scale)
M <- sparse.model.matrix(~ 0 + factor(ct_vec))   # columns = cell types
colnames(M) <- sub("^factor\\(ct_vec\\)", "", colnames(M))
celltype_counts <- Matrix::colSums(M)            # number of cells per type
# Sum per cell type, then divide by counts to get MEANS (genes x celltypes)
sum_by_ct <- cp10k %*% M                         # (genes x celltypes), sums
avg_by_ct <- sum_by_ct %*% Diagonal(x = 1 / as.numeric(celltype_counts))
# avg_by_ct: mean CP10k expression for each gene in each cell type (linear)

# 5) Specificity S_g: proportion of each gene's expression accounted for by each cell type
row_totals <- Matrix::rowSums(avg_by_ct)
row_totals[row_totals == 0] <- .Machine$double.eps
Sg <- avg_by_ct / row_totals    
colnames(Sg) <- colnames(M)

write.table(Sg,'./cell_type/Sg.tsv',row.names=T,sep='\t',quote=F)

# 6) For each cell, select the top decile genes 
get_top_Sg_genes <- function(Sg, pct = 0.10, min_per_ct = 1, na.rm = TRUE) {
  stopifnot(is.matrix(Sg) || inherits(Sg, "Matrix"))
  ct_names <- colnames(Sg)
  res <- vector("list", length(ct_names))
  names(res) <- ct_names
  genes <- rownames(Sg)

  for (j in seq_along(ct_names)) {
    v <- Sg[, j]
    #### remove non-expressed genes
    v2 <- v[v>0]
    ## cut-off value
    cutoff <- quantile(v2,0.9)

    add <- genes[v>=cutoff] 
    res[[j]] <- add
  }
  return(res)
}

# Run it (top 10%)
top10pct_Sg_genes <- get_top_Sg_genes(Sg, pct = 0.10)

## Count number of genes
count <- lengths(top10pct_Sg_genes)

saveRDS(top10pct_Sg_genes,'./cell_type/top10_gene.rds')
write.csv(count,'./cell_type/top10_gene_cell_count.csv')

# ------------------------------------------------------------------------------
# More refined cell types
# ------------------------------------------------------------------------------
idcol <- 'cell_states' # more refined definition of cell subpopulations
Idents(obj) <- obj[[idcol, drop = TRUE]]

# 3) Build a linear expression matrix (CP10k) from raw counts
cts <- GetAssayData(obj, assay = "RNA", slot = "counts")  # genes x cells (dgCMatrix)
libsize <- Matrix::colSums(cts)

# Keep only cells with valid identity and nonzero library size
ct_vec <- as.character(Idents(obj))
keep_cells <- !is.na(ct_vec) & libsize > 0
cts      <- cts[, keep_cells, drop = FALSE]
libsize  <- libsize[keep_cells]
ct_vec   <- ct_vec[keep_cells]

# Normalize counts per cell to CP10k (linear)
cp10k <- t( t(cts) / libsize ) * 1e4   # genes x cells (dgCMatrix)

# QC: keep genes expressed in >1% of cells (any type)
#    "expressed" = counts > 0 in a cell
detected_prop <- Matrix::rowMeans(cts > 0)  # proportion of cells with nonzero counts
gene_keep <- detected_prop > 0.01           # >1% 
(sum(gene_keep==F)) # 19,334 genes removed
cp10k <- cp10k[gene_keep, , drop = FALSE]   # filter early for speed

# 4) Average expression per gene in each cell type (on linear scale)
M <- sparse.model.matrix(~ 0 + factor(ct_vec))   # columns = cell types
colnames(M) <- sub("^factor\\(ct_vec\\)", "", colnames(M))
celltype_counts <- Matrix::colSums(M)            # number of cells per type
# Sum per cell type, then divide by counts to get MEANS (genes x celltypes)
sum_by_ct <- cp10k %*% M                         # (genes x celltypes), sums
avg_by_ct <- sum_by_ct %*% Diagonal(x = 1 / as.numeric(celltype_counts))
# avg_by_ct: mean CP10k expression for each gene in each cell type (linear)

# 5) Specificity S_g: proportion of each gene's expression accounted for by each cell type
row_totals <- Matrix::rowSums(avg_by_ct)
row_totals[row_totals == 0] <- .Machine$double.eps
Sg <- avg_by_ct / row_totals
colnames(Sg) <- colnames(M)

write.table(Sg,'./cell_state/Sg.tsv',row.names=T,sep='\t',quote=F)

# 6) For each cell, select the top decile genes
get_top_Sg_genes <- function(Sg, pct = 0.10, min_per_ct = 1, na.rm = TRUE) {
  stopifnot(is.matrix(Sg) || inherits(Sg, "Matrix"))
  ct_names <- colnames(Sg)
  res <- vector("list", length(ct_names))
  names(res) <- ct_names
  genes <- rownames(Sg)

  for (j in seq_along(ct_names)) {
    v <- Sg[, j]
    #### remove non-expressed genes
    v2 <- v[v>0]
    ## cut-off value
    cutoff <- quantile(v2,0.9)

    add <- genes[v>=cutoff]
    res[[j]] <- add
  }
  return(res)
}

# Run it (top 10%)
top10pct_Sg_genes <- get_top_Sg_genes(Sg, pct = 0.10)

## Count number of genes
count <- lengths(top10pct_Sg_genes)

saveRDS(top10pct_Sg_genes,'./cell_state/top10_gene.rds')
write.csv(count,'./cell_state/top10_gene_cell_count.csv')

