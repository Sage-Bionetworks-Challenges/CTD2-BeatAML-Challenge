suppressPackageStartupMessages(library(tidyverse))
library(reticulate)
library(pheatmap)


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
  fontsize = 8,
  clustering_distance_rows = "correlation", 
  clustering_distance_cols = "correlation"
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
mu <- apply(rnaseq, 1, mean)
SD <- apply(rnaseq, 1, sd)
smoothScatter(mu, SD)
filtered.rnaseq <- rnaseq[mu > 5 & SD > 3.5,] %>%
  rownames_to_column("Gene") %>%
  left_join(gene_symbols, by = "Gene") %>%
  select(-Gene) %>%
  column_to_rownames("Symbol")
  

# Plot the Pearson correlations as a heatmap
pheatmap(
  cor(t(filtered.rnaseq), drug.mat, method = "pearson", use = "complete.obs"),
  main = "Gene vs. Drug Correlations",
  fontsize = 8,
  clustering_distance_rows = "correlation", 
  clustering_distance_cols = "correlation"
)


# Compute the matrix of correlation p-values
#   (source: http://www.sthda.com/upload/rquery_cormat.r)
cor.pmat <- function(x, ...) {
  mat <- as.matrix(x)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
