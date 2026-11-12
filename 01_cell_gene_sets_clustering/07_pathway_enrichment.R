# Over-representation 
setwd('/path/to/csPRS/gene_list/genes_specitificy_top10/merge_251127')
library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)

## 1. Your data ------------------------------------
# Genes selected for each cluster
gene_sets <- readRDS('top10_group_10.rds')
# Background genes [all genes passing QC]
info <- fread('../cell_type/Sg.tsv',data.table=F)
genes_bg <- info$V1

# Map gene symbols to ENTREZID
all_symbols <- unique(c(unlist(gene_sets), genes_bg))
id_map <- bitr(all_symbols,
               fromType = "SYMBOL",
               toType   = "ENTREZID",
               OrgDb    = org.Hs.eg.db)
entrez_map <- setNames(id_map$ENTREZID, id_map$SYMBOL)
# background in Entrez
bg_entrez <- na.omit(entrez_map[genes_bg])

## 2. Get pathways from MSigDB ----------------------
###### Hallmark gene sets
msig_h <- msigdbr(species = "Homo sapiens", category = "H")
term2gene_h <- msig_h %>% dplyr::select(gs_name, gene_symbol)

hallmark_list <- lapply(names(gene_sets), function(nm) {
  sel <- intersect(gene_sets[[nm]], genes_bg)      # enforce universe

  res <- enricher(gene      = sel,
                  TERM2GENE = term2gene_h,
                  universe  = genes_bg)

  if (is.null(res)) return(NULL)

  df <- as.data.frame(res)
  df$gene_set <- nm
  df
})

hallmark_res <- do.call(rbind, hallmark_list)
head(hallmark_res)

###### KEGG pathway sets
kegg_list <- lapply(names(gene_sets), function(nm) {
  sel <- intersect(gene_sets[[nm]], genes_bg)
  sel_entrez <- na.omit(entrez_map[sel])
  if (length(sel_entrez) == 0) return(NULL)

  res <- enrichKEGG(gene         = sel_entrez,
                    universe     = bg_entrez,
                    organism     = "hsa",
                    pAdjustMethod = "BH",
		    keyType      = "ncbi-geneid")
  if (is.null(res)) return(NULL)

  res <- setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  df <- as.data.frame(res)
  df$gene_set <- nm
  df
})

kegg_res <- do.call(rbind, kegg_list)
head(kegg_res)

###### GO terms 
go_list <- lapply(names(gene_sets), function(nm) {
  sel <- intersect(gene_sets[[nm]], genes_bg)
  sel_entrez <- na.omit(entrez_map[sel])
  if (length(sel_entrez) == 0) return(NULL)

  res <- enrichGO(gene          = sel_entrez,
                  universe      = bg_entrez,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  readable      = TRUE)
  if (is.null(res)) return(NULL)

  df <- as.data.frame(res)
  df$gene_set <- nm
  df
})

go_res <- do.call(rbind, go_list)
head(go_res)

out <- rbind(hallmark_res,kegg_res,go_res)

write.csv(out,'pathway_over_represent.csv',row.names=F)
