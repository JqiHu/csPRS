# Get reference list of genes with positions!
# !!!! R 4.4
setwd('/path/to/csPRS/gene_list')
library(data.table)
library(biomaRt) # for gene position

## ---------- Install (run once) ----------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("HGNChelper","org.Hs.eg.db","AnnotationDbi",
          "TxDb.Hsapiens.UCSC.hg19.knownGene","GenomicRanges","dplyr","tibble","stringr")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(to_install)) BiocManager::install(to_install, ask=FALSE, update=FALSE)

## ---------- Library ----------
library(HGNChelper)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(dplyr)
library(tibble)
library(stringr)

# First, read average expression per cell
dat1 <- fread('average_expression_by_celltype.csv',data.table=F)
dat2 <- fread('average_expression_by_cellstate.csv',data.table=F)
## Get gene list
genes <- union(dat1$gene,dat2$gene)
length(genes) # 33,538 

# Then map!
## ---------- Main function ----------
map_genes_hg19 <- function(genes, keep_non_std_chr = FALSE) {
  stopifnot(is.character(genes))
  genes_in <- trimws(genes)

  # 1) Normalize symbols (fix outdated ones)
  hgnc_fix <- HGNChelper::checkGeneSymbols(genes_in)
  # Use Suggested.Symbol when available; otherwise keep input as-is
  norm_symbol <- ifelse(is.na(hgnc_fix$Suggested.Symbol), hgnc_fix$x, hgnc_fix$Suggested.Symbol)

  # 2) Map via org.Hs.eg.db using both SYMBOL and ALIAS
  # SYMBOL path
  map_sym <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = unique(norm_symbol),
                                   keytype = "SYMBOL",
                                   columns = c("SYMBOL","ENTREZID","ENSEMBL")) %>%
    as_tibble() %>%
    rename(matched_symbol = SYMBOL)

  # ALIAS path (try original strings too, not only normalized)
  map_alias <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = unique(c(genes_in, norm_symbol)),
                                     keytype = "ALIAS",
                                     columns = c("ALIAS","SYMBOL","ENTREZID","ENSEMBL")) %>%
    as_tibble() %>%
    rename(input_key = ALIAS, matched_symbol = SYMBOL)

  # Build a long “mapping” table keyed by original input
  # Join back by either exact normalized symbol or alias match
  df_in <- tibble(input = genes_in, norm_symbol = norm_symbol)

  # SYMBOL joins
  m1 <- df_in %>%
    left_join(map_sym %>% distinct(matched_symbol, ENTREZID, ENSEMBL),
              by = c("norm_symbol" = "matched_symbol"), keep = TRUE) %>%
    mutate(mapping_source = ifelse(!is.na(ENTREZID) | !is.na(ENSEMBL), "SYMBOL", NA_character_))

  # For entries still unmapped, try ALIAS
  need_alias <- m1 %>% filter(is.na(ENTREZID) & is.na(ENSEMBL)) %>% pull(input) %>% unique()
  m2 <- tibble()
  if (length(need_alias)) {
    m2 <- df_in %>%
      filter(input %in% need_alias) %>%
      left_join(map_alias %>% distinct(input_key, matched_symbol, ENTREZID, ENSEMBL),
                by = c("input" = "input_key")) %>%
      mutate(mapping_source = ifelse(!is.na(ENTREZID) | !is.na(ENSEMBL), "ALIAS", NA_character_))
  }

  map_all <- bind_rows(
    m1 %>% filter(!(input %in% need_alias) | !is.na(mapping_source)),
    m2
  ) %>%
    distinct(input, norm_symbol, matched_symbol, ENTREZID, ENSEMBL, mapping_source)

  # 3) Get hg19 gene positions from TxDb (Entrez-based)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  gr_genes <- genes(txdb)  # GRanges, mcols has "gene_id" (Entrez IDs)
  # Prepare lookup table
  kg_df <- as_tibble(gr_genes) %>%
    transmute(
      ENTREZID = as.character(mcols(gr_genes)$gene_id),
      seqnames = as.character(seqnames),
      start = start,
      end = end,
      strand = as.character(strand)
    )

  # 4) Join positions
  res <- map_all %>%
    left_join(kg_df, by = "ENTREZID") %>%
    mutate(
      chromosome = seqnames,
      # Helper flags
      chr = str_remove(chromosome %||% "", "^chr"),
      chr_std = chr %in% c(as.character(1:22), "X", "Y", "M", "MT")
    )

  if (!keep_non_std_chr) {
    res <- res %>% filter(is.na(chromosome) | chr_std)
  }

  # 5) Build a GRanges (when positions available)
  gr <- res %>%
    filter(!is.na(chromosome) & !is.na(start) & !is.na(end)) %>%
    mutate(seqnames = chromosome) %>%
    { GRanges(seqnames = .$seqnames,
              ranges   = IRanges(start = .$start, end = .$end),
              strand   = ifelse(is.na(.$strand), "*", .$strand),
              input    = .$input,
              norm_symbol = .$norm_symbol,
              matched_symbol = .$matched_symbol,
              ENTREZID = .$ENTREZID,
              ENSEMBL  = .$ENSEMBL) }

  list(
    table   = res %>% select(input, norm_symbol, matched_symbol, ENTREZID, ENSEMBL,
                             chromosome, start, end, strand, mapping_source) %>% arrange(input),
    granges = gr
  )
}

# Run for all genes
out <- map_genes_hg19(genes)
out.df <- as.data.frame(out$table)

out.df$chromosome <- as.numeric(str_remove_all(out.df$chromosome,'chr'))
out.df <- out.df[out.df$chromosome %in% 1:22,]
dim(out.df) # 20,907

write.csv(out.df,'gene_info_2.csv',row.names=F)

missing <- setdiff(genes,out.df$input)

write.table(missing,'gene_missing2.txt',
	  row.names=F,sep='\t',quote=F)
