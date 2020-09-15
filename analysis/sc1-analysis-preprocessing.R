suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(tidyverse))
use.reticulate <- FALSE
if(use.reticulate) {
    p_load(reticulate)
} else {
    p_load(synapser)    
}
p_load(pheatmap)
suppressPackageStartupMessages(p_load(reshape2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))


# `synasper` not working on local computer, so using `reticulate` to
# log in to Synapse instead
if(use.reticulate) {
    use_condaenv("synapse-r")
    synapseclient <- reticulate::import('synapseclient')
    syn <- synapseclient$Synapse()
    syn$login(silent = T)
} else {
    synLogin()
}

suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))
suppressPackageStartupMessages(p_load("ComplexHeatmap"))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

auc.synIds <- list("training" = "syn21212913",
                   "leaderboard" = "syn21671784",
                   "validation" = "syn21212720")

expr.synIds <- list("training" = "syn21212911",
                    "leaderboard" = "syn21671825",
                    "validation" = "syn21212718")

phases <- names(auc.synIds)
names(phases) <- phases

## Read in the drug response AUCs
raw.aucs <-
    llply(auc.synIds,
          .fun = function(synId) {
              if(use.reticulate) {
                  auc_file <- syn$get(synId)$path
              } else {
                  auc_file <- synGet(synId)$path
              }
              aucs <- suppressMessages(read_csv(auc_file))
              aucs
          })

## Read in the expression matrices
raw.exprs <-
    llply(expr.synIds,
          .fun = function(synId) {
              # RNAseq data
              if(use.reticulate) {
                  rnaseq_file <- syn$get(synId)$path
              } else {
                  rnaseq_file <- synGet(synId)$path
              }
              rnaseq <- suppressMessages(read_csv(rnaseq_file))
              rnaseq
          })

## Append pseudo-drug (mean response), GRD, to AUC matrices
aucs <-
    llply(raw.aucs,
          .fun = function(auc.df) {
              grd <- ddply(auc.df, .variables = c("lab_id"),
                           .fun = function(df) {
                               data.frame(inhibitor = "GRD",
                                          auc = mean(df$auc, na.rm = TRUE),
                                          stringsAsFactors = FALSE)
                           })
              grd <- grd[, colnames(auc.df)]
              ## Add GRD to the tibble then transform it to a matrix (lab_id x inhibitor)
              drug.mat <- auc.df %>%
                  bind_rows(grd) %>% 
                  spread(inhibitor, auc) %>%
                  column_to_rownames("lab_id")
              drug.mat
          })

## Create symbol and ensembl versions of the expression matrices
ensg.exprs <-
    llply(raw.exprs,
          .fun = function(df) {
              df %>%
                  select(-Symbol) %>%
                  column_to_rownames("Gene")
          })

symbol.exprs <-
    llply(raw.exprs,
          .fun = function(df) {
              rnaseq_symbol <- df %>% select(-Gene)
              rnaseq_symbol <- as.data.frame(rnaseq_symbol)
              tmp <-
                  ddply(rnaseq_symbol, .variables = c("Symbol"),
                        .fun = function(df) {
                            colMeans(df[, !(colnames(df) == "Symbol")])
                        })
              rnaseq_symbol <- tmp[, !(colnames(tmp) == "Symbol")]
              rownames(rnaseq_symbol) <- tmp$Symbol
              rnaseq_symbol
          })

## Create long versions of the AUC and (symboled-based) expression matrices
long.aucs <-
    llply(aucs,
          .fun = function(df) {
              drug.long <- melt(as.matrix(df))
              colnames(drug.long) <- c("sample", "drug", "auc")
              drug.long
          })

long.symbol.exprs <-
    llply(symbol.exprs,
          .fun = function(df) {
              gene.long <- melt(as.matrix(df))
              colnames(gene.long) <- c("gene", "sample", "expr")
              gene.long
          })

cat("Computing gene vs drug correlations\n")
gene.symbol.expr.cors <-
    llply(phases,
          .fun = function(phase) {
              drug.long <- long.aucs[[phase]]              
              ddply(drug.long, .variables = c("drug"),
                    .parallel = TRUE,
                    .fun = function(drug.df) {
                        ddply(long.symbol.exprs[[phase]],
                              .variables = c("gene"),
                              .fun = function(gene.df) {
                                  df <- merge(drug.df, gene.df)
                                  pearson.ct <-
                                      cor.test(df$auc, df$expr,
                                               method="pearson",
                                               use="pairwise.complete.obs")
                                  spearman.ct <-
                                      cor.test(df$auc, df$expr,
                                               method="spearman",
                                               use="pairwise.complete.obs")
                                  ret <-
                                      data.frame(pearson.pval = as.numeric(pearson.ct$p.value),
                                                 pearson.cor = as.numeric(pearson.ct$estimate),
                                                 spearman.pval = as.numeric(spearman.ct$p.value),
                                                 spearman.cor = as.numeric(spearman.ct$estimate))
                              })
                    })
          })

## Compute GSVA of expression data (against KEGG)

suppressPackageStartupMessages(p_load("msigdbr"))
kegg <- msigdbr(species = "Homo sapiens", subcategory = "CP:KEGG")
## mp <- symbols.to.ensg.mapping(unique(kegg$gene_symbol))

## kegg <- merge(as.data.frame(kegg), mp, by.x = c("gene_symbol"), by.y = c("symbol"))
## kegg.gene.sets <- dlply(kegg, .variables = c("gs_name"), .fun = function(df) as.character(df$ensg))
kegg.gene.sets <-
    dlply(kegg, .variables = c("gs_name"), .fun = function(df) as.character(df$gene_symbol))

suppressPackageStartupMessages(library("GSVA"))

cat("Computing GSVA for KEGG pathways\n")
long.kegg.eses <-
    llply(symbol.exprs,
          .fun = function(rnaseq_symbol) {
              es <- gsva(as.matrix(rnaseq_symbol),
                         kegg.gene.sets, parallel.sz=num.processes, verbose=TRUE)
              es.long <- melt(as.matrix(es))
              colnames(es.long) <- c("pathway", "sample", "es")
              es.long
          })

cat("Computing correlations of KEGG enrichment score vs drug\n")
kegg.es.cors <- 
    llply(phases,
          .fun = function(phase) {
              drug.long <- long.aucs[[phase]]
              ddply(drug.long, .variables = c("drug"),
                    .parallel = TRUE,
                    .fun = function(drug.df) {
                        ddply(long.kegg.eses[[phase]],
                              .variables = c("pathway"),
                              .fun = function(gene.df) {
                                  df <- merge(drug.df, gene.df)
                                  pearson.ct <-
                                      cor.test(df$auc, df$es,
                                               method="pearson",
                                               use="pairwise.complete.obs")
                                  spearman.ct <-
                                      cor.test(df$auc, df$es,
                                               method="spearman",
                                               use="pairwise.complete.obs")
                                  ret <-
                                      data.frame(pearson.pval = as.numeric(pearson.ct$p.value),
                                                 pearson.cor = as.numeric(pearson.ct$estimate),
                                                 spearman.pval = as.numeric(spearman.ct$p.value),
                                                 spearman.cor = as.numeric(spearman.ct$estimate))
                              })
                    })
          })


## Store everything to Synapse
analysis.folder.synId <- "syn22397014"

for(phase in phases) {
    file <- paste0(phase, "_", "aucs.csv")
    cat(paste0("Writing file ", file, "\n"))
    write.table(file = file, aucs[[phase]], sep = ",", 
                row.names = TRUE, col.names = TRUE, quote = FALSE)
    f <- File(file, parentId = analysis.folder.synId, synapseStore = TRUE)
    cat(paste0("Storing file ", file , " to Synapse\n"))
    synStore(f)
    file.remove(file)

    file <- paste0(phase, "_", "ensg_exprs.csv")
    cat(paste0("Writing file ", file, "\n"))
    write.table(file = file, ensg.exprs[[phase]], sep = ",", 
                row.names = TRUE, col.names = TRUE, quote = FALSE)
    f <- File(file, parentId = analysis.folder.synId, synapseStore = TRUE)
    cat(paste0("Storing file ", file , " to Synapse\n"))
    synStore(f)
    file.remove(file)

    file <- paste0(phase, "_", "symbol_exprs.csv")
    cat(paste0("Writing file ", file, "\n"))
    write.table(file = file, symbol.exprs[[phase]], sep = ",", 
                row.names = TRUE, col.names = TRUE, quote = FALSE)
    f <- File(file, parentId = analysis.folder.synId, synapseStore = TRUE)
    cat(paste0("Storing file ", file , " to Synapse\n"))
    synStore(f)
    file.remove(file)

    file <- paste0(phase, "_", "long_aucs.csv")
    cat(paste0("Writing file ", file, "\n"))
    write.table(file = file, long.aucs[[phase]], sep = ",", 
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    f <- File(file, parentId = analysis.folder.synId, synapseStore = TRUE)
    cat(paste0("Storing file ", file , " to Synapse\n"))
    synStore(f)
    file.remove(file)

    file <- paste0(phase, "_", "long_symbol_exprs.csv")
    cat(paste0("Writing file ", file, "\n"))
    write.table(file = file, long.symbol.exprs[[phase]], sep = ",", 
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    f <- File(file, parentId = analysis.folder.synId, synapseStore = TRUE)
    cat(paste0("Storing file ", file , " to Synapse\n"))
    synStore(f)
    file.remove(file)

    file <- paste0(phase, "_", "gene_symbol_expr_cors.csv")
    cat(paste0("Writing file ", file, "\n"))
    write.table(file = file, gene.symbol.expr.cors[[phase]], sep = ",", 
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    f <- File(file, parentId = analysis.folder.synId, synapseStore = TRUE)
    cat(paste0("Storing file ", file , " to Synapse\n"))
    synStore(f)
    file.remove(file)

    file <- paste0(phase, "_", "long_kegg_eses.csv")
    cat(paste0("Writing file ", file, "\n"))
    write.table(file = file, long.kegg.eses[[phase]], sep = ",", 
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    f <- File(file, parentId = analysis.folder.synId, synapseStore = TRUE)
    cat(paste0("Storing file ", file , " to Synapse\n"))
    synStore(f)
    file.remove(file)

    file <- paste0(phase, "_", "kegg_es_cors.csv")
    cat(paste0("Writing file ", file, "\n"))
    write.table(file = file, kegg.es.cors[[phase]], sep = ",", 
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    f <- File(file, parentId = analysis.folder.synId, synapseStore = TRUE)
    cat(paste0("Storing file ", file , " to Synapse\n"))
    synStore(f)
    file.remove(file)
}

cat("Exiting successfully\n")
q(status = 0)
              
    
## 
# drug vs drug
##

## Compute drug/drug correlations
drug.cor.mats <-
    llply(aucs,
          .fun = function(drug.mat) cor(drug.mat, use = "pairwise.complete.obs"))





if(FALSE) {
cors <-
    ddply(drug.long, .variables = c("drug"),
          .parallel = TRUE,
          .fun = function(drug.df) {
              ddply(gene.long, .variables = c("gene"),
                    .fun = function(gene.df) {
                        df <- merge(drug.df, gene.df)
                        ct <- cor.test(df$auc, df$expr, method="pearson", use="pairwise.complete.obs")
                        ret <- data.frame(pval = as.numeric(ct$p.value), cor = as.numeric(ct$estimate))
                    })
              })
}
load(".Rdata.corr")

mu <- apply(rnaseq, 1, mean)
SD <- apply(rnaseq, 1, sd)

expressed.genes <- names(mu)[mu > 0]

# Since Symbol is not unique, Gene must be used as rownames instead. To
# help change the row labels on the heatmap later, create a mapping of
# the Ensembl ID to symbol
gene_symbols <- suppressMessages(read_csv(rnaseq_file)) %>%
  select(Gene, Symbol)

cors.bh <-
    ddply(subset(cors, gene %in% expressed.genes), .variables = c("drug"),
          .fun = function(df) {
              data.frame(gene = df$gene, cor=df$cor, pval = df$pval, padj = p.adjust(df$pval))
          })

sig.cors <- subset(cors.bh, padj < 0.1)


sig.cors <- merge(sig.cors, gene_symbols, by.x = "gene", by.y = "Gene")

o <- order(abs(sig.cors$cor), decreasing=TRUE)
sig.cors <- sig.cors[o,]

sig.cor.mat <- acast(data = sig.cors[, c("drug", "Symbol", "cor")], formula = Symbol ~ drug)

num.non.na.cols <- colSums(!is.na(sig.cor.mat))
o <- order(num.non.na.cols)
sig.cor.mat <- sig.cor.mat[, o]

num.non.na.rows <- rowSums(!is.na(sig.cor.mat))
o <- order(num.non.na.rows)
sig.cor.mat <- sig.cor.mat[o, ]

sig.cor.genes <- as.character(unique(subset(cors.bh, padj < 0.1 & !(drug %in% c("Venetoclax")))$gene))
sig.cors.all <- subset(cors.bh, gene %in% sig.cor.genes)
sig.cors.all <- merge(sig.cors.all, gene_symbols, by.x = "gene", by.y = "Gene")
sig.cor.mat.all <- acast(data = sig.cors.all[, c("drug", "Symbol", "cor")], formula = Symbol ~ drug)


drug.families <- read.table("41586_2018_623_MOESM3_ESM-s11.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
drug.families$value <- "Yes"
drug.family.mat <- acast(data = drug.families, formula = family ~ inhibitor, fill = NA)

common.drugs <- intersect(colnames(drug.family.mat), colnames(sig.cor.mat))
drug.family.mat <- drug.family.mat[, common.drugs]
sub <- sig.cor.mat[, common.drugs]
num.non.na.cols <- colSums(!is.na(sub))
o <- order(num.non.na.cols)
num.non.na.cols <- num.non.na.cols[o]
drug.family.mat <- drug.family.mat[, o]

synLogin()

synId <- "syn22267537"
tbl <- synTableQuery(paste0("SELECT * FROM ", synId))
tbl <- as.data.frame(tbl)

## Exclude "gold" (which must mean gold-standard because its pearson and spearman = 1)
tbl <- subset(tbl, team != "gold")
tbl$pearson <- as.numeric(tbl$pearson)

drug.mean.tbl <-
    ddply(tbl, .variables = c("inhibitor"),
          .fun = function(df) data.frame(pearson = mean(df$pearson, na.rm=TRUE)))

if(!(all(common.drugs %in% drug.mean.tbl$inhibitor))) { stop("Missing drugs\n") }
sub.drug.mean.tbl <- subset(drug.mean.tbl, inhibitor %in% common.drugs)
rownames(sub.drug.mean.tbl) <- sub.drug.mean.tbl$inhibitor
sub.drug.mean.tbl <- sub.drug.mean.tbl[names(num.non.na.cols), ]
drug.means <- sub.drug.mean.tbl$pearson
names(drug.means) <- sub.drug.mean.tbl$inhibitor

make.plot <- function(mat, sig.vec, perf.vec, cor.mat, drug.gene.cor.mat) {
    breaks = c(1, 10, 100, 1000)
    col_ha = HeatmapAnnotation(sig = anno_barplot(log10(as.numeric(sig.vec)),
                                                  axis_param = list(at = log10(breaks), labels =breaks)),
                               perf = anno_barplot(perf.vec))
    names(col_ha) = c("# Sig Genes", "Mean Pearson")
    row.flag <- rowSums(is.na(mat)) != ncol(mat)
    h <- Heatmap(mat[row.flag,], cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = col_ha,
                 column_names_rot = 45, na_col = "white", show_heatmap_legend = FALSE,
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8))
    h2 <- Heatmap(cor.mat, cluster_rows = FALSE, cluster_columns = FALSE)
    h3 <- Heatmap(drug.gene.cor.mat, cluster_rows = FALSE, cluster_columns = FALSE)    
    h %v% h2 %v% h3
}

drug.cor.mat <- cor(drug.mat, use = "pairwise.complete.obs")

symbols.to.ensg.mapping <- function(symbols) {
  suppressPackageStartupMessages(p_load(mygene))
  dummy <- data.frame(query = symbols, ensembl = NA)
  bm <- tryCatch({queryMany(symbols, scopes="symbol", fields=c("ensembl.gene"), species="human")}, error = function(e) { return(dummy) })
  flag <- grepl(pattern="ensembl", colnames(bm))
  if(length(which(flag)) != 1) {
    stop(paste0("Could not find ensembl col in: ", paste(colnames(bm), collapse=" "), "\n"))
  }
  ensg.col <- colnames(bm)[flag]
  bm <- bm[, c("query", ensg.col)]
  lst <- bm[,ensg.col]
  names(lst) <- bm$query
  bm <- ldply(lst, .fun = function(comp) data.frame(unlist(comp)))
  colnames(bm) <- c("symbol", "ensg")
  bm <- bm[!(bm$ensg %in% c("")),,drop=F]
  bm <- bm[!is.na(bm$ensg) & !is.na(bm$symbol),,drop=F]
  bm
}



es.cors.bh <-
    ddply(subset(es.cors, !(drug %in% c("Venetoclax"))), .variables = c("drug"),
          .fun = function(df) {
              data.frame(pathway = df$pathway, cor=df$cor, pval = df$pval, padj = p.adjust(df$pval))
          })


sig.cor.pathways <- as.character(unique(subset(es.cors.bh, padj < 0.1)$pathway))
sig.cors.pathways.all <- subset(es.cors.bh, pathway %in% sig.cor.pathways)
sig.cor.pathway.mat.all <- acast(data = sig.cors.pathways.all[, c("drug", "pathway", "cor")], formula = pathway ~ drug)
sig.cor.pathways <- subset(es.cors.bh, padj < 0.1)
sig.cor.pathway.mat.all <- as.data.frame(acast(data = sig.cor.pathways[, c("drug", "pathway", "cor")], formula = pathway ~ drug))


png("all-sig-genes-vs-drug-family-perf-ordered.png")
o <- order(drug.means)
drugs <- names(drug.means)[o]
drugs <- drugs[!(drugs %in% c("Venetoclax"))]
drugs <- drugs[hclust(dist(drug.cor.mat[drugs,drugs], method="euclidean"), method="complete")$order]
for(d in drugs[!(drugs %in% colnames(sig.cor.pathway.mat.all))]) { sig.cor.pathway.mat.all[,d] <- NA }
## make.plot(drug.family.mat[, drugs], num.non.na.cols[drugs], drug.means[drugs], drug.cor.mat[drugs, drugs], sig.cor.mat.all[, drugs])
make.plot(drug.family.mat[, drugs], num.non.na.cols[drugs], drug.means[drugs], drug.cor.mat[drugs, drugs], sig.cor.pathway.mat.all[, drugs])
d <- dev.off()

stop("stop")

png("multiple-sig-genes-vs-drug-family-perf-ordered.png")
sub.flag <- as.numeric(num.non.na.cols) > 1
o <- order(drug.means[sub.flag])
make.plot(drug.family.mat[, sub.flag][, o], num.non.na.cols[sub.flag][o], drug.means[sub.flag][o])
d <- dev.off()


df <- data.frame(drug = names(num.non.na.cols), sig.genes = as.numeric(num.non.na.cols))
df <- subset(df, sig.genes > 2)
df <- df[order(df$sig.genes),]
df$drug <- factor(df$drug, levels = df$drug)
g <- ggplot(data = df, aes(x = drug, y = sig.genes)) + geom_col()
g <- g + theme(axis.text.x = element_text(angle = 45, hjust=1))
g <- g + ylab("Number of Significantly Correlated Genes")
png("num-sig-genes-vs-drugs.png")
print(g)
d <- dev.off()

heatmap(sig.cor.mat, Colv = NA, Rowv = NA, scale = "none")

num.non.na.rows <- rowSums(!is.na(sig.cor.mat))

heatmap(sig.cor.mat[num.non.na.rows > 1, ], Colv = NA, Rowv = NA, scale = "none")

heatmap(sig.cor.mat[num.non.na.rows == 1, ], Colv = NA, Rowv = NA, scale = "none")



stop("stop")

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
