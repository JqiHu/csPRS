# Generate gene list based on expression data downloaded 
# !!! R 4.4
setwd('/path/to/csPRS/gene_list')
library(Seurat)
library(SeuratDisk)

# 1) Load the Seurat object made from freeze.h5ad
#Convert("hca_heart_global_ctl200723_freeze.h5ad", dest = "h5seurat", overwrite = TRUE)
obj <- LoadH5Seurat("hca_heart_global_ctl200723_freeze.h5seurat")

### Some explorations
# How many cells and genes?
dim(obj)              # 33,538 cells x 486,134 genes

# Cell-level metadata
head(obj@meta.data)

# What assays are present?
Assays(obj) # RNA

# Default assay
DefaultAssay(obj)

# Which embeddings are available?
Reductions(obj) # pca and umap

# View counts and normalized expression for a few genes
head(GetAssayData(obj, slot = "counts")[, 1:5])
head(GetAssayData(obj, slot = "data")[, 1:5])

## Check whether the data has been lognormalized
data   <- GetAssayData(obj, slot = "data")
counts <- GetAssayData(obj, slot = "counts")

summary(as.vector(counts[1:100, 1:10]))
summary(as.vector(data[1:100, 1:10]))
range(counts@x, na.rm=TRUE)
range(data@x, na.rm = TRUE)
obj@commands
#### It seems that the data has not been normalized

# Make sure the data slot has log1p-normalized values
obj.norm <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)
data   <- GetAssayData(obj.norm, slot = "data")
range(data@x, na.rm = TRUE) # normalized range!

# 2) Pick the identity as your cell-type column
idcol <- 'cell_type'
#idcol <- 'cell_states' # more refined definition of cell subpopulations
Idents(obj.norm) <- obj.norm[[idcol, drop = TRUE]]

# 3) Get average expression (log-normalized 'data' slot by default)
avg <- AverageExpression(obj.norm, assays = "RNA", slot = "data")$RNA
# avg is a genes x cell_types matrix

# 4) Save as wide CSV (genes as rows, cell types as columns)
avg_df <- as.data.frame(avg)
dim(avg_df)
head(avg_df)
avg_df$gene <- rownames(avg_df)
avg_df <- avg_df[, c("gene", setdiff(colnames(avg_df), "gene"))]

write.csv(avg_df, file = "average_expression_by_celltype.csv", row.names = FALSE)

### for more refined cell types
idcol <- 'cell_states' # more refined definition of cell subpopulations
Idents(obj.norm) <- obj.norm[[idcol, drop = TRUE]]

# 3) Get average expression (log-normalized 'data' slot by default)
avg <- AverageExpression(obj.norm, assays = "RNA", slot = "data")$RNA
# avg is a genes x cell_types matrix

# 4) Save as wide CSV (genes as rows, cell types as columns)
avg_df <- as.data.frame(avg)
dim(avg_df)
head(avg_df)
avg_df$gene <- rownames(avg_df)
avg_df <- avg_df[, c("gene", setdiff(colnames(avg_df), "gene"))]

write.csv(avg_df, file = "average_expression_by_cellstate.csv", row.names = FALSE)
