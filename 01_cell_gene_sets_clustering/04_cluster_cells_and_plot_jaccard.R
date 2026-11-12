# Plot heatmaps for jaccard similarity
setwd('/path/to/csPRS/gene_list/genes_specitificy_top10/merge_251127')
library(data.table)
library(stringr)
library(ggplot2)
library(ggpubr)
rm(list=ls())

## read data
dat <- fread('/path/to/csPRS/gene_list/shared_gene/specitificy_top10.csv', data.table=FALSE)

# --- Ordering by clustering on covariance profiles ---------------------------
# build wide matrix 
cov_mat <- xtabs(Jaccard ~ Var1 + Var2, data = dat)
# similarity = correlation between traits' covariance profiles
# use pairwise complete obs to handle NAs (diagonals or missing pairs)
S2 <- cor(cov_mat, use = "pairwise.complete.obs")
hc <- hclust(as.dist(1 - S2))  # distance = 1 - correlation
ord <- rownames(S2)[hc$order]

dat$Var1 <- factor(dat$Var1,levels=ord)
dat$Var2 <- factor(dat$Var2,levels=ord)



#### Decide the number of clusters
# --- Prepare matrix of similarity profiles ------------------------------
# order rows/cols according to the dendrogram
X  <- S2[ord, ord]

# center columns (no scaling) so total variance = sum of squares
Xc <- scale(X, center = TRUE, scale = FALSE)

# --- Function: variance explained by k clusters -------------------------
var_explained_k <- function(hc, Xc, k) {
  cl <- cutree(hc, k = k)
  cl <- cl[rownames(Xc)]          # make sure order matches Xc
  
  # total sum of squares
  TSS <- sum(Xc^2, na.rm = TRUE)
  
  # between-cluster sum of squares
  BSS <- 0
  for (g in unique(cl)) {
    idx   <- which(cl == g)
    mu_g  <- colMeans(Xc[idx, , drop = FALSE], na.rm = TRUE)
    BSS   <- BSS + length(idx) * sum(mu_g^2, na.rm = TRUE)
  }
  
  BSS / TSS  # proportion of variance explained
}

# --- Try different numbers of clusters ----------------------------------
ks <- 2:20
var_df <- data.frame(
  k          = ks,
  prop_var   = sapply(ks, function(k) var_explained_k(hc, Xc, k))
)

# decide elbow point 
# by observation 
elbow_k=9

# quick elbow plot
pdf('jaccard_sim_cluster_elbow.pdf')
plot(var_df$k, var_df$prop_var, type = "b",
     xlab = "k", ylab = "Proportion of variance explained")
points(elbow_k, var_df$prop_var[var_df$k==elbow_k], pch = 19, cex = 1.5)
dev.off()

# define clusters
clusters <- cutree(hc, k = elbow_k)

# long to wide
dat2 <- dat[order(dat$Var1,dat$Var2),]
library(tidyr)

df_wide <- dat2 %>%
  pivot_wider(
    names_from = Var2,
    values_from = Jaccard
  )
clusters2 <- clusters[match(df_wide$Var1,names(clusters))]
df_wide <- cbind.data.frame(clusters2,df_wide)

writexl::write_xlsx(df_wide,'jaccard_similarity_specitificy_wide.xlsx')

### Edit the cluster manually
cluster <- readxl::read_xlsx('jaccard_similarity_group_edit.xlsx')
cluster <- cluster[order(cluster$Cluster2),]

dat$Var1 <- factor(dat$Var1,levels=cluster$Var1)
dat$Var2 <- factor(dat$Var2,levels=cluster$Var1)
# for each cluster, find the start and end index in the ordered list
cluster_blocks <- lapply(unique(cluster$Cluster2), function(cl) {
  idx <- which(cluster$Cluster2== cl)
  c(start = min(idx), end = max(idx))
})
names(cluster_blocks) <- paste0("cl", unique(cluster$Cluster2))

p <- ggplot(dat, aes(Var1, Var2, fill = Jaccard)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red",
                      name = "Jaccard Similarity") +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  )

# add rectangles for diagonal clusters
for (i in seq_along(cluster_blocks)) {
  blk <- cluster_blocks[[i]]
  p <- p + 
    annotate("rect",
             xmin = blk["start"] - 0.5,
             xmax = blk["end"] + 0.5,
             ymin = blk["start"] - 0.5,
             ymax = blk["end"] + 0.5,
             fill = NA, color = "black", size = 0.7)
}

pdf("jaccard_sim.pdf", width = 14, height = 12)
print(p)
dev.off()

