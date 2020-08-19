suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pacman))
p_load(reticulate)
p_load(pheatmap)


# `synasper` not working on local computer, so using `reticulate` to
# log in to Synapse instead
use_condaenv("synapse-r")
synapseclient <- reticulate::import('synapseclient')
syn <- synapseclient$Synapse()
syn$login(silent = T)


## 
# drug vs drug
##

# Goldstandard AUCs
auc_file <- syn$get("syn21212720")$path
aucs <- suppressMessages(read_csv(auc_file))

# Create pseudo-drug (mean response), GRD
grd <- aucs %>% 
  group_by(lab_id) %>%
  summarize(auc = mean(auc))  %>%
  add_column(inhibitor = "GRD", .after = "lab_id")

# Add GRD to the tibble then transform it to a matrix (lab_id x inhibitor)
drug.mat <- aucs %>%
  bind_rows(grd) %>% 
  spread(inhibitor, auc) %>%
  column_to_rownames("lab_id")

# Plot the Pearson correlations as a heatmap
pheatmap(
  cor(drug.mat, method = "pearson", use = "complete.obs"), 
  main = "Drug vs. Drug Correlations",
  fontsize = 5, 
  clustering_distance_rows = "correlation", 
  clustering_distance_cols = "correlation",
  file = "drug_v_drug.png",
  width = 10, height = 10
)


## 
# gene vs drug
##

# RNAseq data
rnaseq_file <- syn$get("syn21212718")$path
rnaseq <- suppressMessages(read_csv(rnaseq_file)) %>%
  select(-Symbol) %>%
  column_to_rownames("Gene")

# Since Symbol is not unique, Gene must be used as rownames instead. To
# help change the row labels on the heatmap later, create a mapping of
# the Ensembl ID to symbol
gene_symbols <- suppressMessages(read_csv(rnaseq_file)) %>%
  select(Gene, Symbol)

# Check to make sure expression data is within the range we want
hist(rnaseq[,2])

# Peak around -5 so normalize data so that peak is at 0
rnaseq <- rnaseq - min(rnaseq) + .001
hist(rnaseq[,2])

# For every gene, find the mean and standard deviation, then perform a 
# smoothScatter to help filter the data
# mu <- apply(rnaseq, 1, mean)
# SD <- apply(rnaseq, 1, sd)
# smoothScatter(mu, SD)
# filtered.rnaseq <- rnaseq[mu > 5 & SD > 3.5,] %>%
#   rownames_to_column("Gene") %>%
#   left_join(gene_symbols, by = "Gene") %>%
#   select(-Gene) %>%
#   column_to_rownames("Symbol")
  

# Plot the Pearson correlations as a heatmap
# pheatmap(
#   cor(t(filtered.rnaseq), drug.mat, method = "pearson", use = "complete.obs"),
#   main = "Gene vs. Drug Correlations",
#   fontsize = 8,
#   clustering_distance_rows = "correlation", 
#   clustering_distance_cols = "correlation"
# )


# Compute the matrix of correlations and p-values
#   (reference: http://www.sthda.com/upload/rquery_cormat.r)
cor.mats <- function(mat1, mat2) {
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2)
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  p.mat <- c.mat <- matrix(NA, n1, n2)
  for (i in 1:n1) {
    for (j in 1:n2) {
      tmp <- cor.test(mat1[, i], mat2[, j], method = "pearson", use = "complete.obs")
      p.mat[i, j] <- tmp$p.value
      c.mat[i, j] <- tmp$estimate
    }
  }
  rownames(p.mat) <- rownames(c.mat) <- colnames(mat1)
  colnames(p.mat) <- colnames(c.mat) <- colnames(mat2)
  list(r = c.mat, p = p.mat)
}

# corrs <- cor.mats(drug.mat, t(rnaseq))
# p.mat <- corrs$p
# c.mat <- corrs$r
# 
# Export matrices to RDS for subsequent faster loading times.
# saveRDS(p.mat, file = "p_values.rds")
# saveRDS(c.mat, file = "correlations.rds")

# Load RDS for correlation and p-value matrices
p.mat <- readRDS("p_values.rds")
c.mat <- readRDS("correlations.rds")

# Remap Ensembl ID (colnames) back to gene symbols
colnames(p.mat) <- gene_symbols$Symbol[which(gene_symbols$Gene %in% colnames(p.mat))]
colnames(c.mat) <- gene_symbols$Symbol[which(gene_symbols$Gene %in% colnames(c.mat))]
