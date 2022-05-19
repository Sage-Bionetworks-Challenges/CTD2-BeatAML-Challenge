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
suppressPackageStartupMessages(p_load(pvclust))
suppressPackageStartupMessages(p_load(data.table))
suppressPackageStartupMessages(p_load(cowplot))

plot.dir <- "plots"
dir.create(plot.dir)

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

suppressPackageStartupMessages(p_load("fgsea"))
suppressPackageStartupMessages(p_load("msigdbr"))

suppressPackageStartupMessages(p_load(biomaRt))
suppressPackageStartupMessages(p_load(AnnotationHub))
suppressPackageStartupMessages(p_load(ensembldb))
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(openxlsx))
suppressPackageStartupMessages(p_load(stringi))
suppressPackageStartupMessages(p_load(rcdk))
suppressPackageStartupMessages(p_load(rpubchem))
suppressPackageStartupMessages(p_load(fingerprint))
suppressPackageStartupMessages(p_load(reshape2))
suppressPackageStartupMessages(p_load(ggpubr))
suppressPackageStartupMessages(p_load(GGally))
suppressPackageStartupMessages(p_load(ggbeeswarm))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1


phases <- c("training", "leaderboard", "validation")
names(phases) <- phases

limit.genes.to.protein.coding <- function(genes, use.symbols = TRUE) {

    use.biomaRt <- TRUE
    exclude.mt <- FALSE
    if(use.biomaRt) {
        ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        pc.tbl <- getBM(attributes = c("ensembl_gene_id","external_gene_name","description","chromosome_name"),
                        filters = 'biotype', values = c('protein_coding'), mart = ensembl)
        if(exclude.mt) {
            pc.tbl <- subset(pc.tbl, chrosome_name != "MT")
        }
        keys <- as.character(pc.tbl$external_gene_name)
        if(!use.symbols) { keys <- as.character(pc.tbl$ensembl_gene_id) }
        genes <- genes[genes %in% keys]
        return(genes)
    }
    if(exclude.mt) { stop("Have not implemented exclude.mt for annotationHub -- but is easy\n") }

    ah <- AnnotationHub()
    flag <- (ah$species == "Homo sapiens") & (ah$genome == "GRCh38") & (ah$dataprovider == "Ensembl") & (ah$rdataclass == "EnsDb")
    ah2 <- ah[flag, ]
    ## as.data.frame(mcols(ah2))[1:10,c("title"),drop=FALSE]
    edb <- ah2[["AH73881"]]

    ## keytypes(edb)
    ## columns(edb)
    keys <- keys(edb, "GENENAME")
    columns <- c("GENEID", "ENTREZID", "GENEBIOTYPE")
    tbl <- ensembldb::select(edb, keys, columns, keytype = "GENENAME")
    pc.tbl <- subset(tbl, GENEBIOTYPE == "protein_coding")

    keys <- as.character(pc.tbl$GENENAME)
    if(!use.symbols) { keys <- as.character(pc.tbl$GENEID) }
    genes <- genes[genes %in% keys]
    return(genes)
}

read.csv.from.synapse <- function(synId, ...) {
    if(use.reticulate) {
        file <- syn$get(synId)$path
    } else {
        file <- synGet(synId)$path
    }
    cat(paste0("Reading ", file, " from synapse\n"))
    # ret <- suppressMessages(read.table(file, ...))
    ret <- suppressMessages(fread(file))
    cat(paste0("Done reading ", file, " from synapse\n"))
    return(ret)
}

read.csv.table.from.synapse <- function(synId, sep=",") {
    if(use.reticulate) {
        file <- syn$get(synId)$path
    } else {
        file <- synGet(synId)$path
    }
    ret <- suppressMessages(read.table(file, sep=sep, header=TRUE))
}


## Read in the prediction results from the validation phase
synId <- "syn22267537"
validation.result.tbl <- synTableQuery(paste0("SELECT * FROM ", synId))
## validation.result.tbl <- as.data.frame(validation.result.tbl)
validation.result.tbl <- validation.result.tbl$asDataFrame()

## Exclude "gold" (which must mean gold-standard because its pearson and spearman = 1)
## Signal and ymemari are the same (empirically, they have identical scores)
teams.to.drop <- c("gold", "Signal")
validation.result.tbl <- subset(validation.result.tbl, !(team %in% teams.to.drop))
validation.result.tbl$pearson <- as.numeric(validation.result.tbl$pearson)
validation.result.tbl$spearman <- as.numeric(validation.result.tbl$spearman)

cat(paste0(length(unique(validation.result.tbl$team)), " teams submitted to SC1\n")) 

cat(paste0("Teams: ", paste(sort(unique(validation.result.tbl$team)), sep=", "), "\n"))

drug.mean.validation.result.tbl <-
    ddply(validation.result.tbl, .variables = c("inhibitor"),
          .fun = function(df) data.frame(pearson = mean(df$pearson, na.rm=TRUE), spearman = mean(df$spearman, na.rm=TRUE)))

o <- order(drug.mean.validation.result.tbl$spearman, decreasing=FALSE)
drug.mean.validation.result.tbl <- drug.mean.validation.result.tbl[o, ]

summarized.validation.result.tbl <-
  ddply(validation.result.tbl, .variables = c("inhibitor"),
        .fun = function(df) {
                 data.frame(spearman = median(as.numeric(df$spearman), na.rm=TRUE), pearson = median(as.numeric(df$pearson), na.rm=TRUE))
               })

o <- order(summarized.validation.result.tbl$spearman, decreasing=FALSE)
summarized.validation.result.tbl <- summarized.validation.result.tbl[o, ]


write.table(file=paste0(plot.dir, "/mean-drug-perf.tsv"), drug.mean.validation.result.tbl, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

write.table(file=paste0(plot.dir, "/median-drug-perf.tsv"), summarized.validation.result.tbl, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

validation.result.tbl$inhibitor <- factor(validation.result.tbl$inhibitor, levels = summarized.validation.result.tbl$inhibitor)

validation.result.tbl$color <- "black"
validation.result.tbl[validation.result.tbl$team == "AAUH","color"] <- "red"

is_low_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x))
}

is_high_outlier <- function(x) {
  return(x > quantile(x, 0.75) + 1.5 * IQR(x))
}

is_outlier <- function(x) {
  return(is_low_outlier(x) | is_high_outlier(x))
}

validation.result.tbl <- 
    validation.result.tbl %>% 
    group_by(inhibitor) %>% 
    mutate(high_outlier = is_high_outlier(spearman)) %>%
    mutate(low_outlier = is_low_outlier(spearman)) 

write.table(file=paste0(plot.dir, "/validation-results-with-outliers.tsv"), validation.result.tbl, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

g <- ggplot()
g <- g + geom_boxplot(data = validation.result.tbl, aes(x = inhibitor, y = spearman), outlier.shape=NA)
g <- g + geom_point(data = subset(validation.result.tbl, high_outlier == TRUE), aes(x = inhibitor, y = spearman, colour = team))
g <- g + geom_point(data = subset(validation.result.tbl, low_outlier == TRUE), aes(x = inhibitor, y = spearman))
g <- g + xlab("Inhibitor") + ylab("Spearman Correlation\n(Observed vs Predicted)")
g <- g + theme(axis.text.y = element_text(size = 4))
g <- g + labs(colour='Team')
g <- g + coord_flip()
g.drug <- g

png(paste0(plot.dir, "/beat-aml-sc1-per-drug-score-distributions.png"))
print(g.drug)
d <- dev.off()
pdf(paste0(plot.dir, "/beat-aml-sc1-per-drug-score-distributions.pdf"))
print(g.drug)
d <- dev.off()

cat(paste0("Number of high outliers: ", nrow(subset(validation.result.tbl, high_outlier == TRUE)), "\n"))
cat(paste0("Number of teams w/ high outliers: ", length(unique(subset(validation.result.tbl, high_outlier == TRUE)$team)), "\n"))
high.outlier.teams <- as.data.frame(table(subset(validation.result.tbl, high_outlier == TRUE)$team))
colnames(high.outlier.teams) <- c("team", "num.high.outliers")
cat(paste0("Number of teams w/ multiple high-outliers: ", nrow(subset(high.outlier.teams, num.high.outliers > 1)), "\n"))

cat(paste0("Number of low outliers: ", nrow(subset(validation.result.tbl, low_outlier == TRUE)), "\n"))
cat(paste0("Number of teams w/ low outliers: ", length(unique(subset(validation.result.tbl, low_outlier == TRUE)$team)), "\n"))
low.outlier.teams <- as.data.frame(table(subset(validation.result.tbl, low_outlier == TRUE)$team))
colnames(low.outlier.teams) <- c("team", "num.low.outliers")
cat(paste0("Number of teams w/ multiple low-outliers: ", nrow(subset(low.outlier.teams, num.low.outliers > 1)), "\n"))
outlier.teams <- merge(high.outlier.teams, low.outlier.teams, all=TRUE)
print(outlier.teams)

write.table(file=paste0(plot.dir, "/outlier-drug-perf.tsv"), outlier.teams, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

indices <- seq(from=nrow(summarized.validation.result.tbl),to=2)
names(indices) <- summarized.validation.result.tbl[indices,"inhibitor"]
consecutive.diffs <- 
   ldply(indices,
         .fun = function(i) {
                  data.frame(spearman.diff = summarized.validation.result.tbl[i,"spearman"] - summarized.validation.result.tbl[i-1,"spearman"],
                             pearson.diff = summarized.validation.result.tbl[i,"pearson"] - summarized.validation.result.tbl[i-1,"pearson"])
         })
colnames(consecutive.diffs)[1] <- "inhibitor"
# Venetoclax has the larger improvement relative to other drugs
print(consecutive.diffs[which.max(abs(consecutive.diffs$spearman.diff)),])
o <- order(consecutive.diffs$spearman.diff, decreasing=TRUE)
consecutive.diffs <- consecutive.diffs[o,]
print(consecutive.diffs)

write.table(file=paste0(plot.dir, "/drug-perf-differences.tsv"), consecutive.diffs, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")



## We only use training -- don't read the others, as it takes a lot of time.
expr.synIds <- list("training" = "syn21212911")
## expr.synIds[["leaderboard"]] <- "syn21671825"
## expr.synIds[[""validation"]] <- "syn21212718"

## Read in the expression matrices
raw.exprs <-
    llply(expr.synIds,
          .fun = function(synId) {
              # RNAseq data
	      cat("Reading expr ", synId, "\n")
              if(use.reticulate) {
                  rnaseq_file <- syn$get(synId)$path
              } else {
                  rnaseq_file <- synGet(synId)$path
              }
	      cat("Reading expr file ", rnaseq_file, "\n")
              # rnaseq <- suppressMessages(read_csv(rnaseq_file))
              rnaseq <- suppressMessages(fread(rnaseq_file))
	      cat("Done reading expr file ", rnaseq_file, "\n")
              rnaseq
          })

## Read in the intermediate results generated by sc1-analysis-preprocessing.R

## We are only using training data, so only read in that

gene.symbol.expr.cor.synIds <- list()
gene.symbol.expr.cor.synIds[["training"]] <- "syn22397862"
##gene.symbol.expr.cor.synIds[["leaderboard"]] <- "syn22398045"
##gene.symbol.expr.cor.synIds[["validation"]] <-  "syn22398135"

gene.symbol.expr.cors <-
    llply(gene.symbol.expr.cor.synIds,
          .fun = function(synId) read.csv.from.synapse(synId, sep=",", header=TRUE))

symbol.expr.synIds <-
    list("training" = "syn22397592",
         "leaderboard" = "syn22397931",
         "validation" = "syn22398088")

symbol.exprs <-
    llply(symbol.expr.synIds,
          .fun = function(synId) read.csv.table.from.synapse(synId))

auc.synIds <-
    list("training" = "syn22397379",
         "leaderboard" = "syn22397867",
         "validation" = "syn22398053")

aucs <-
    llply(auc.synIds,
          .fun = function(synId) read.csv.table.from.synapse(synId))

dnaseq.synIds <-
    list("training" = "syn21297970",
         "leaderboard" = "syn21297972",
         "validation" = "syn21297966")

snvs <-
    llply(dnaseq.synIds,
          .fun = function(synId) read.csv.table.from.synapse(synId))

snv.mats <-
    llply(snvs,
          .fun = function(df) {
              df[, "mut"] <- 1
              acast(data = df, formula = Hugo_Symbol ~ lab_id, fun.aggregate = function(x) ifelse(length(x) > 0, 1, 0),
                    value.var = "mut", fill = 0)
          })

## Read in the highly variable and expressed genes from our FIMM/OHSU Beat AML manuscript
synId <- "syn21899141"
synId <- "syn20628997"
aml.expressed.variable.genes.tbl <- read.csv.table.from.synapse(synId, sep="\t")
gene_symbols <- as.data.frame(unique(raw.exprs[["training"]][, c("Gene", "Symbol")]))
aml.expressed.variable.symbols <- subset(gene_symbols, Gene %in% as.character(aml.expressed.variable.genes.tbl$gene))$Symbol

## filter.genes <- limit.genes.to.protein.coding(aml.expressed.variable.symbols)
filter.genes <- aml.expressed.variable.symbols

## Read in the drug family annotations
drug.families <- read.table("41586_2018_623_MOESM3_ESM-s11.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
drug.families$value <- "Yes"
drug.family.mat <- acast(data = drug.families, formula = family ~ inhibitor, fill = NA)
drug.family.mat <- as.data.frame(drug.family.mat)


## Read in drug annotations
synId <- "syn23002655"
drug.anno <- read.csv.table.from.synapse(synId, sep="\t")

## Look up the fingerpints
smiles <- as.character(drug.anno$smiles)
names(smiles) <- drug.anno$inhibitor

fingerprints <-
    llply(smiles,
          .fun = function(sm) {
              if(is.null(sm)) { return(NA) }
              if(is.na(sm)) { return(NA) }
              parsed <- parse.smiles(sm)
              if(is.null(parsed)) { return(NA) }
              if(is.null(parsed[[1]])) { return(NA) }              
              get.fingerprint(parsed[[1]])
          })

## Define expressed genes from the distribution of the mean expression
## of the _training_ data
mu <- apply(symbol.exprs[["training"]], 1, mean)
SD <- apply(symbol.exprs[["training"]], 1, sd)

cutoff <- 0

plot(density(mu))
abline(v=cutoff)

## filter.genes <- names(mu)[mu > cutoff]




gene.symbol.expr.cors.bh <-
    llply(gene.symbol.expr.cors,
          .fun = function(cors) {
              ddply(subset(cors, gene %in% filter.genes), .variables = c("drug"),
                    .fun = function(df) {
                        data.frame(gene = df$gene, cor=df$pearson.cor,
                                   pval = df$pearson.pval,
                                   padj = p.adjust(df$pearson.pval))
                    })
          })

gene.symbol.expr.cors.sig <-
    llply(gene.symbol.expr.cors.bh,
          .fun = function(df) subset(df, padj < 0.05))

gene.symbol.expr.cors.sig.mats <-
    llply(gene.symbol.expr.cors.sig,
          .fun = function(sig.cors) {
              sig.cor.mat <-
                  acast(data = sig.cors[, c("drug", "gene", "cor")],
                        formula = gene ~ drug)
          })

sig.cor.mat <- gene.symbol.expr.cors.sig.mats[["training"]]

tmp.mat <- gene.symbol.expr.cors.sig[["training"]]
drugs.to.exclude <- c("Venetoclax")
tmp.mat <- subset(tmp.mat, !(drug %in% drugs.to.exclude))
sig.cor.genes <- as.character(unique(tmp.mat$gene))
sig.cors.all <- subset(gene.symbol.expr.cors[["training"]], gene %in% sig.cor.genes)
sig.cor.mat.all <- acast(data = sig.cors.all[, c("drug", "gene", "pearson.cor")],
                         formula = gene ~ drug)

l <- list(drug.family.mat, sig.cor.mat, sig.cor.mat.all)
all.drugs <- unique(unlist(lapply(l, colnames)))
all.drugs <- unique(c(all.drugs, gene.symbol.expr.cors$training$drug))

rename.drug.map <- data.frame(drug = unique(all.drugs),
                              safe.name = make.names(unique(all.drugs)),
                              stringsAsFactors = FALSE)

rownames(rename.drug.map) <- rename.drug.map$safe.name

renamed.aucs <- aucs
for(nm in names(renamed.aucs)) {
    missing <- colnames(renamed.aucs[[nm]])[!(colnames(renamed.aucs[[nm]]) %in% rownames(rename.drug.map))]
    if(length(missing) > 0) {
        cat(paste0("Missing some drugs: ", paste(missing, sep=", "), "\n"))
    }
    renamed.aucs[[nm]] <- renamed.aucs[[nm]][, colnames(renamed.aucs[[nm]]) %in% rownames(rename.drug.map)]
    colnames(renamed.aucs[[nm]]) <- rename.drug.map[colnames(renamed.aucs[[nm]]), "drug"]
}

auc.ranges <-
    llply(renamed.aucs,
          .fun = function(aucs) {
              apply(aucs, 2, mad, na.rm=TRUE)
          })

l <- list(drug.family.mat, sig.cor.mat, sig.cor.mat.all, renamed.aucs[["training"]])
all.drugs <- unique(unlist(lapply(l, colnames)))

rownames(drug.mean.validation.result.tbl) <- drug.mean.validation.result.tbl$inhibitor
missing.predictions <- all.drugs[!(all.drugs %in% rownames(drug.mean.validation.result.tbl))]

prediction.drugs <- rownames(drug.mean.validation.result.tbl)
all.drugs <- intersect(all.drugs, prediction.drugs)
# mean.drug.predictions <- drug.mean.validation.result.tbl$pearson
# Spearman was the primary metric in the challenge
mean.drug.predictions <- drug.mean.validation.result.tbl$spearman
names(mean.drug.predictions) <- rownames(drug.mean.validation.result.tbl)

median.drug.predictions <- summarized.validation.result.tbl$spearman
names(median.drug.predictions) <- rownames(summarized.validation.result.tbl)

all.drugs <- unique(c(all.drugs, "GRD"))
all.drugs <- all.drugs[!(all.drugs %in% drugs.to.exclude)]

drug.cor.mat <- cor(renamed.aucs[["training"]], use = "pairwise.complete.obs", method="spearman")

sig.cor.mat <- as.data.frame(sig.cor.mat)
sig.cor.mat.all <- as.data.frame(sig.cor.mat.all)
drug.cor.mat <- as.data.frame(drug.cor.mat)

training.auc.ranges <- auc.ranges[["training"]]
validation.auc.ranges <- auc.ranges[["validation"]]

structural.similarity.mat <-
    ldply(fingerprints,
          .fun = function(f1) {
              ret <- ldply(fingerprints,
                           .fun = function(f2) {
                               if( ( class(f1) == "logical" ) && is.na(f1) ) {
                                   return(data.frame(sim = NA))
                               }
                               if( ( class(f2) == "logical" ) && is.na(f2) ) {
                                   return(data.frame(sim = NA))
                               }
                               data.frame(sim = distance(f1, f2, method = "tanimoto"))
                           })
              colnames(ret)[1] <- "drug2"
              ret
          })
colnames(structural.similarity.mat)[1] <- "drug1"
structural.similarity.mat <- acast(structural.similarity.mat, formula = drug1 ~ drug2)
structural.similarity.mat <- as.data.frame(structural.similarity.mat)


for(d in all.drugs) {
    if(!(d %in% colnames(drug.family.mat))) { drug.family.mat[, d] <- NA }
    if(!(d %in% colnames(sig.cor.mat))) { sig.cor.mat[, d] <- NA }
    if(!(d %in% rownames(sig.cor.mat))) { sig.cor.mat[d, ] <- NA }    
    if(!(d %in% colnames(drug.cor.mat))) { drug.cor.mat[, d] <- NA }    
    if(!(d %in% colnames(sig.cor.mat.all))) { sig.cor.mat.all[, d] <- NA }
    if(!(d %in% rownames(sig.cor.mat.all))) { sig.cor.mat.all[d, ] <- NA }    
    if(!(d %in% colnames(structural.similarity.mat))) { structural.similarity.mat[, d] <- NA }
    if(!(d %in% rownames(structural.similarity.mat))) { structural.similarity.mat[d, ] <- NA }    
    if(!(d %in% names(mean.drug.predictions))) { mean.drug.predictions[d] <- NA }
    if(!(d %in% names(training.auc.ranges))) { training.auc.ranges[d] <- NA }    
    if(!(d %in% names(validation.auc.ranges))) { validation.auc.ranges[d] <- NA }    
    
}    
drug.family.mat <- drug.family.mat[, all.drugs]
sig.cor.mat <- sig.cor.mat[, all.drugs]
drug.cor.mat <- drug.cor.mat[all.drugs, all.drugs]
sig.cor.mat.all <- sig.cor.mat.all[all.drugs, all.drugs]
mean.drug.predictions <- mean.drug.predictions[all.drugs]
training.auc.ranges <- training.auc.ranges[all.drugs]
validation.auc.ranges <- validation.auc.ranges[all.drugs]
num.sig.genes <- colSums(!is.na(sig.cor.mat))
structural.similarity.mat <- structural.similarity.mat[all.drugs, all.drugs]

source("plot-dendrogram.R")

drugs <- all.drugs

drug.mat <- drug.cor.mat[drugs, drugs]

if(FALSE) {
pvc <- pvclust(drug.mat, method.hclust = "complete",
               method.dist = "euclidean", nboot=1000,
               iseed = 1234, parallel = TRUE)
}

column.hc <- hclust(dist(drug.mat, method="euclidean"), method="complete")
column_dend <- as.dendrogram(column.hc)
row_dend <-
    as.dendrogram(hclust(dist(t(drug.mat), method="euclidean"),
                         method="complete"))

num.clusters <- 7
sv <- structural.similarity.mat

structural.similarity.mat <- sv
structural.similarity.mat <- as.matrix(structural.similarity.mat )
diag(structural.similarity.mat) <- NA

#file <- "multi_uni_plot_data.csv"
#plotData <- read.table(file, sep=",", header=TRUE)
# Top performer/Rasmus model: multi-drug vs uni-drug on _validation_ data.
file <- "rasmus-multi-vs-uni-validation.tsv"
plotData <- read.table(file, sep="\t", header=TRUE)
#multi_vs_uni_variate <- plotData$pearson.cor.multi - plotData$pearson.cor.uni
multi_vs_uni_variate <- plotData$spearman.cor.multi - plotData$spearman.cor.uni
names(multi_vs_uni_variate) <- plotData$inhibitor

##cutoff <- 1
## flag <- !is.na(structural.similarity.mat) &
##     structural.similarity.mat > cutoff
## structural.similarity.mat[flag] <- cutoff
png(paste0(plot.dir, "/sc1-training-data-vs-validation-performance.png"))
pdf(paste0(plot.dir, "/sc1-training-data-vs-validation-performance.pdf"))
ht <- make.plot.dendro(drug.family.mat[, drugs],
                       num.sig.genes[drugs], # training
		       multi_vs_uni_variate[drugs], # validation
                       mean.drug.predictions[drugs], # validation scores
                       training.auc.ranges[drugs], # training drug AUC MADs
                       drug.cor.mat[drugs, drugs], # training
                       structural.similarity.mat[drugs, drugs],
                       sig.cor.mat.all[, drugs], 
                       column_dend,
                       row_dend,
                       column_split = num.clusters,
                       row_split = num.clusters,                 
                       column_title_gp = gpar(fill = c("red", "blue", "green", "orange", "purple", "yellow", "brown"), font = 1:7))

d <- dev.off()

##ct <- cutree(column.hc, k=num.clusters)
##drug.clusters <- data.frame(inhibitor = names(ct), cluster = unname(ct))

## Get the column dendrogram from the plot
col.dend <- column_dend(ht)
clusters <- 1:length(col.dend)
names(clusters) <- clusters

get.children.labels <- function(node) {
    if(is.leaf(node)) { return(attr(node, "label")) }
    unlist(llply(node, .fun = function(kid) get.children.labels(kid)))
}

drug.clusters <-
    ldply(clusters,
          .fun = function(cluster) {
              data.frame(inhibitor =
                             unlist(llply(col.dend[[cluster]],
                                          .fun = function(node) get.children.labels(node))))
          })
colnames(drug.clusters) <- c("cluster", "inhibitor")
drug.clusters <- drug.clusters[, c("inhibitor", "cluster")]

write.table(file="drug-clusters.tsv", drug.clusters, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

drug.range.tbl <- data.frame(inhibitor = names(training.auc.ranges), training.range = unname(training.auc.ranges))
drug.range.tbl <- merge(drug.range.tbl, data.frame(inhibitor = names(validation.auc.ranges), validation.range = unname(validation.auc.ranges)))
mean.perf.tbl <- data.frame(inhibitor = names(mean.drug.predictions), mean.correlation = unname(mean.drug.predictions))
median.perf.tbl <- data.frame(inhibitor = names(median.drug.predictions), median.correlation = unname(median.drug.predictions))

multi.vs.uni.df <- data.frame(inhibitor = names(multi_vs_uni_variate), multi.vs.uni = as.vector(multi_vs_uni_variate))

all.tbl <- merge(drug.clusters, drug.range.tbl, all.x = TRUE)
all.tbl <- merge(all.tbl, mean.perf.tbl, all.x = TRUE)
all.tbl <- merge(all.tbl, multi.vs.uni.df, all.x = TRUE)
all.tbl <- all.tbl[!(all.tbl$inhibitor == "GRD"),]
all.tbl$cluster <- factor(all.tbl$cluster, levels = sort(unique(all.tbl$cluster)))
all.tbl$inhibitor.name <- make.names(all.tbl$inhibitor)


# Plot multi vs uni for each cluster

g1 <- ggplot(data = all.tbl)
g1 <- g1 + geom_boxplot(aes(x = cluster, y = mean.correlation))
g1 <- g1 + ylab("Mean Pearson Correlation\n(over Teams)")

g2 <- ggplot(data = all.tbl)
g2 <- g2 + geom_boxplot(aes(x = cluster, y = training.range))
g2 <- g2 + ylab("Dynamic Range\n(AUC MAD over Samples)")

g3 <- ggplot(data = all.tbl)
g3 <- g3 + geom_boxplot(aes(x = cluster, y = multi.vs.uni))
g3 <- g3 + ylab("Multi- vs Uni-variate Performance")


png(paste0(plot.dir, "/sc1-training-clusters-vs-validation-performance.png"))
grid.arrange(g1,g2,g3)
d <- dev.off()

pdf(paste0(plot.dir, "/sc1-training-clusters-vs-validation-performance.pdf"))
grid.arrange(g1,g2,g3)
d <- dev.off()

png(paste0(plot.dir, "/sc1-training-clusters-vs-multi-vs-uni-performance.png"))
print(g3)
d <- dev.off()

pdf(paste0(plot.dir, "/sc1-training-clusters-vs-multi-vs-uni-performance.pdf"))
print(g3)
d <- dev.off()

all.tbl <- merge(drug.clusters, drug.range.tbl, all.x = TRUE)
all.tbl <- merge(all.tbl, mean.perf.tbl, all.x = TRUE)
all.tbl <- all.tbl[!(all.tbl$inhibitor == "GRD"),]
ohsu.drug.targets <- read.table("ohsu-drug-targets.tsv", sep="\t", header=TRUE)
# Ensure a minimum drug binding
cutoff <- -log10(100*10^-9)
ohsu.drug.targets$log.num.targets <- log10(ohsu.drug.targets$num.targets)
ohsu.drug.targets[ohsu.drug.targets$max.mean.pchembl < cutoff, "num.targets"] <- NA
ohsu.drug.targets[ohsu.drug.targets$max.mean.pchembl < cutoff, "log.num.targets"] <- NA
all.tbl <- merge(all.tbl, ohsu.drug.targets[, c("inhibitor", "num.targets", "log.num.targets", "max.mean.pchembl")], all.x=TRUE)
all.tbl <- merge(all.tbl, multi.vs.uni.df, all.x=TRUE)
all.tbl$inhibitor.name <- make.names(all.tbl$inhibitor)

cluster.names <- unique(all.tbl$cluster)
my_comparisons <-
    llply(cluster.names[cluster.names != "4"], .fun = function(cl) c("4", cl))

all.tbl$cluster <- factor(all.tbl$cluster, levels = sort(unique(all.tbl$cluster)))
g1 <- ggplot(data = all.tbl)
g1 <- g1 + geom_boxplot(aes(x = cluster, y = mean.correlation))
g1 <- g1 + ylab("Mean Pearson Correlation\n(over Teams)")

g2 <- ggplot(data = all.tbl)
g2 <- g2 + geom_boxplot(aes(x = cluster, y = training.range))
g2 <- g2 + ylab("Dynamic Range\n(AUC MAD over Samples)")

png(paste0(plot.dir, "/sc1-training-clusters-vs-validation-performance.png"))
grid.arrange(g1,g2)
d <- dev.off()

pdf(paste0(plot.dir, "/sc1-training-clusters-vs-validation-performance.pdf"))
grid.arrange(g1,g2)
d <- dev.off()








g1 <- ggboxplot(data = all.tbl, x = "cluster", y = "mean.correlation")
g1 <- g1 + ylab("Mean Pearson Correlation\n(over Teams)")
g1 <- g1 + stat_compare_means(comparison = my_comparisons)

g2 <- ggboxplot(data = all.tbl, x = "cluster", y = "training.range")
g2 <- g2 + ylab("Dynamic Range\n(AUC MAD over Samples)")
g2 <- g2 + stat_compare_means(comparison = my_comparisons)

g3 <- ggboxplot(data = all.tbl, x = "cluster", y = "log.num.targets")
g3 <- g3 + ylab("Num Targets\n(Log_10)")
g3 <- g3 + stat_compare_means(comparison = my_comparisons)

g4 <- ggboxplot(data = all.tbl, x = "cluster", y = "multi.vs.uni")
g4 <- g4 + ylab("Multi - Uni")
g4 <- g4 + stat_compare_means(comparison = my_comparisons)

# grid.arrange(g1,g2,g3,g4)

klaeger.reconciled.synId <- "syn26320868"
klaeger.targets.synId <- "syn26320930"

klaeger.reconciled <- read.csv.from.synapse(klaeger.reconciled.synId, sep="\t", comment.char="", header=TRUE)
klaeger.targets <- read.csv.from.synapse(klaeger.targets.synId, sep="\t", header=TRUE)


ki.tbl <- merge(all.tbl, klaeger.reconciled, by.x = c("inhibitor", "inhibitor.name"), by.y = c("inhibitor", "inhibitor.name"))

grd.cors <-
  llply(aucs,
        .fun = function(aucs.df) {
                 grd <- rowMeans(aucs.df, na.rm=TRUE)
                 drugs <- colnames(aucs.df)
                 names(drugs) <- drugs
                 method = "spearman"
                 method = "pearson"
                 ret.grd.cors <- ldply(drugs,
                                       .fun = function(drug) data.frame(grd.cor = cor(aucs.df[,drug], grd, method=method, use='pairwise.complete.obs')))
                 colnames(ret.grd.cors)[1] <- "inhibitor"
                 ret.grd.cors
               })

median.grd.cors <-
  llply(aucs,
        .fun = function(aucs.df) {
                 grd <- apply(aucs.df, 1, function(x) median(x, na.rm=TRUE))
                 drugs <- colnames(aucs.df)
                 names(drugs) <- drugs
                 method = "spearman"
                 method = "pearson"
                 ret.grd.cors <- ldply(drugs,
                                       .fun = function(drug) data.frame(grd.cor = cor(aucs.df[,drug], grd, method=method, use='pairwise.complete.obs')))
                 colnames(ret.grd.cors)[1] <- "inhibitor"
                 ret.grd.cors
               })

all.grd.cors <- list("mean" = grd.cors, "median" = median.grd.cors)

all.grd.plts <- list()

for(met in c("mean", "median")) {
  all.grd.plts[[met]] <- list()
for(ds in c("training", "validation")) {

  all.tbl.with.grd <- merge(all.tbl, all.grd.cors[[met]][[ds]], by.x = c("inhibitor.name"), by.y = c("inhibitor"))
  all.tbl.with.grd$multi.higher <- "FALSE"
  flag <- !is.na(all.tbl.with.grd$multi.vs.uni) & (all.tbl.with.grd$multi.vs.uni > 0)
  all.tbl.with.grd[flag,"multi.higher"] <- "TRUE"
  all.tbl.with.grd$multi.higher <- factor(all.tbl.with.grd$multi.higher, levels = c("TRUE", "FALSE"))

  g <- ggplot()
  g <- g + geom_violin(data=all.tbl.with.grd, aes(x = multi.higher, y = grd.cor))
  g <- g + geom_boxplot(data=all.tbl.with.grd, aes(x = multi.higher, y = grd.cor), width = 0.25)
  g <- g + geom_beeswarm(data=all.tbl.with.grd, aes(x = multi.higher, y = grd.cor))
  g <- g + xlab("Joint Modeling Performance > Independent Modeling Performance")
  g <- g + ylab("Drug Response vs Mean Drug Response\n(Pearson Correlation)")
  all.grd.plts[[met]][[ds]] <- g
  print(g)
  ggsave(paste0(plot.dir, "/multi-minus-uni-vs-", ds, "-", met, "-grd-box.png"))

  print(g)
  ggsave(paste0(plot.dir, "/multi-minus-uni-vs-", ds, "-", met, "-grd-box.pdf"))

  wt <- wilcox.test(grd.cor ~ multi.higher, data=all.tbl.with.grd, alternative="greater")

  cat(paste0(ds, "; ", met, ": Comparing correlation with GRD between those with high vs lower performance in joint vs univariate: wilcox one-sided p-value: ", wt$p.value, "\n"))
  higher <- subset(all.tbl.with.grd, multi.higher=="TRUE")
  qs <- as.numeric(quantile(higher$grd.cor,probs=c(0.25,0.5,0.75)))
  cat(paste0(ds, "; ", met, ": Joint > Uni (n=", nrow(higher), ") correlation with GRD quantiles: 25% = ", qs[1], " 50% = ", qs[2], " 75% = ", qs[3], "\n"))
  lower <- subset(all.tbl.with.grd, multi.higher=="FALSE")
  qs <- as.numeric(quantile(lower$grd.cor,probs=c(0.25,0.5,0.75)))
  cat(paste0(ds, "; ", met, ": Joint > Uni correlation with GRD quantiles: 25% = ", qs[1], " 50% = ", qs[2], " 75% = ", qs[3], "\n"))
  print(wt)
}
}



ki.tbl <- merge(ki.tbl, grd.cors[["training"]], by.x = c("inhibitor.name"), by.y = c("inhibitor"))
colnames(ki.tbl)[colnames(ki.tbl) == "grd.cor"] <- "training.grd.cor"
ki.tbl <- merge(ki.tbl, grd.cors[["validation"]], by.x = c("inhibitor.name"), by.y = c("inhibitor"))
colnames(ki.tbl)[colnames(ki.tbl) == "grd.cor"] <- "validation.grd.cor"

cat(paste0("Num unique kinase inhibitors: ", length(unique(ki.tbl$inhibitor.name)), "\n"))
cat(paste0("Len of kinase inhibitor table: ", nrow(ki.tbl), "\n"))

cols <- c("mean.correlation", "Number.of.targets", "training.range", "validation.range", "training.grd.cor", "validation.grd.cor")
#upper = list(continuous = wrap("cor", method = "spearman"))
# Report pearson correlation p-values (which should be the default for ggpairs)
upper = list(continuous = wrap("cor", method = "pearson"))

p_load(Hmisc)
cor.res <- rcorr(as.matrix(ki.tbl[,cols]),type="pearson")
# write.table(file=paste0(plot.dir, "/perf-targets-range-grd.tsv"), cor(ki.tbl[,cols]), row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
write.table(file=paste0(plot.dir, "/perf-targets-range-grd.tsv"), cor.res$r, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
write.table(file=paste0(plot.dir, "/perf-targets-range-grd-values.tsv"), cor.res$P, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

g <- ggpairs(ki.tbl[,cols], columnLabels = c("Performance", "Num Targets", "MAD\n(Training)", "MAD\n(Validation)", "GRD\n(Training)", "GRD\n(Validation)"), upper = upper, lower = list(continuous = "smooth"))
png(paste0(plot.dir, "/perf-targets-range-grd.png"))
print(g)
d <- dev.off()

pdf(paste0(plot.dir, "/perf-targets-range-grd.pdf"))
print(g)
d <- dev.off()


cat(paste0("Multivarate test of mean correlation association\n"))
lmf <- lm(formula("mean.correlation ~ 0 + Number.of.targets + training.range + validation.range + training.grd.cor + validation.grd.cor"), data = ki.tbl)
print(summary(lmf))

png(paste0(plot.dir, "/sc1-individ-drug-perf-and-correlates.png"), width = 2 * 480)
plot_grid(g.drug, ggmatrix_gtable(g), labels="AUTO")
d <- dev.off()

pdf(paste0(plot.dir, "/sc1-individ-drug-perf-and-correlates.pdf"), width = 2 * 7, onefile=FALSE)
plot_grid(g.drug, ggmatrix_gtable(g), labels="AUTO")
d <- dev.off()



all.cluster.tbl <- ki.tbl %>%
  group_by(cluster) %>%
  dplyr::summarise(across(c(Number.of.targets, training.range, training.grd.cor, mean.correlation), mean, na.rm= TRUE))

cols <- c("mean.correlation", "training.range", "training.grd.cor")
g <- ggpairs(all.cluster.tbl[,cols], columnLabels = c("Performance", "MAD", "GRD"), lower = list(continuous = "smooth"))
png(paste0(plot.dir, "/perf-targets-range-grd-cluster.png"))
print(g)
d <- dev.off()

pdf(paste0(plot.dir, "/perf-targets-range-grd-cluster.pdf"))
print(g)
d <- dev.off()

lm_corr_eqn <- function(df, method = "pearson", display.r2 = FALSE, display.pval = FALSE){
    m <- lm(y ~ x, df);
    ct <- cor.test(df$x, df$y, method = method)
    estimate <- as.numeric(ct$estimate)
    if(display.r2 == TRUE) { estimate <- estimate*estimate }
    pval <- ct$p.value
    cat(paste0("method = ", method, " estimate = ", estimate, " pval = ", pval, "\n"))
    eq <- NULL
    digits <- 2
    if((method == "pearson") && (display.r2 == TRUE)) { 
      if(display.pval) { 
        eq <- substitute(italic(r)^2~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=digits, scientific=0),
                              pval = format(pval, digits=digits, scientific=0)))
      } else {
        eq <- substitute(italic(r)^2~"="~est, 
                         list(est = format(estimate, digits=digits, scientific=0)))

      }
    } else if((method == "pearson") && (display.r2 == FALSE)) {
      if(display.pval) {
        eq <- substitute(italic(r)~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=digits, scientific=0),
                              pval = my_format(pval)))
                              # pval = format(pval, digits=digits, scientific=0)))
      } else {
        eq <- substitute(italic(r)~"="~est, 
                         list(est = format(estimate, digits=digits, scientific=0)))

      }
    } else if((method == "spearman") && (display.r2 == FALSE)) {
      if(display.pval) { 
        eq <- substitute(rho~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=digits, scientific=0),
                              pval = format(pval, digits=digits, scientific=0)))
      } else {
        eq <- substitute(rho~"="~est, 
                         list(est = format(estimate, digits=digits, scientific=0)))

      }
    } else {
      stop(paste("lm_corr_eqn does not know how to handle method = ", method,  " display.r2 = ", display.r2, "\n"))
    }
    as.character(as.expression(eq));                 
}

plot.correlation <- function(x, y, labels = NULL, colors = NULL, display.r2 = FALSE, method = "pearson", display.pval = FALSE, xoffset = 0.5, yoffset = 0.8, 
                             sz = 25, geom.text.size = 5, ...) {
  df <- data.frame(x = x, y = y)
  if(!is.null(labels)) {
    df$labels <- labels
  }
  g <- NULL
  if(is.null(labels)) {
    g <- ggplot(df, aes(x = x, y = y))
  } else {
    g <- ggplot(df, aes(x = x, y = y, label = labels))
  }
  if(!is.null(colors)) {
    g <- g + geom_point(aes(colour = colors))
  } else {
    g <- g + geom_point()
  }
  if(!is.null(labels)) {
    g <- g + geom_text(vjust = "inward", hjust = "inward")
  }
  g <- g + geom_smooth(data = df, aes(x = x, y = y), method='lm')
  x.min <- min(df$x, na.rm=TRUE)
  x.max <- max(df$x, na.rm=TRUE)
  y.min <- min(df$y, na.rm=TRUE)
  y.max <- max(df$y, na.rm=TRUE)

  ylimits <- NULL

  xlimits <- ggplot_build(g)$layout$panel_params[[1]]$x.range
  ylimits <- ggplot_build(g)$layout$panel_params[[1]]$y.range

  g <- g + geom_text(size = geom.text.size, x = xlimits[1] + xoffset * (xlimits[2] - xlimits[1]), y = ylimits[1] + yoffset * (ylimits[2] - ylimits[1]), 
                     label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
  g <- g +theme(text = element_text(size = sz),
             axis.text.x = element_text(size=sz),
             axis.text.y = element_text(size=sz),
             axis.title.x = element_text(size=sz),
             axis.title.y = element_text(size=sz),
             title = element_text(size=sz),
             plot.title = element_text(hjust = 0.5, size=sz))
  g
}

g <- plot.correlation(all.cluster.tbl$mean.correlation, all.cluster.tbl$training.range)
g <- g + xlab("Performance") + ylab("Dynamic Range")
g

png(paste0(plot.dir, "/perf-vs-range-cluster.png"))
print(g)
d <- dev.off()

pdf(paste0(plot.dir, "/perf-vs-range-cluster.pdf"))
print(g)
d <- dev.off()

labels <- ki.tbl$inhibitor
flag <- ki.tbl$Number.of.targets < 25 & ki.tbl$training.range > 50
labels[!flag] <- "" 

g <- plot.correlation(ki.tbl$Number.of.targets, ki.tbl$training.range, labels=labels)
g <- g + xlab("Number of Targets") + ylab("Dynamic Range")
g

png(paste0(plot.dir, "/num-targets-vs-range.png"))
print(g)
d <- dev.off()

pdf(paste0(plot.dir, "/num-targets-vs-range.pdf"))
print(g)
d <- dev.off()

png(paste0(plot.dir, "/num-targets-vs-range-outliers.png"))
grid.table(subset(drug.families, inhibitor %in% ki.tbl$inhibitor[flag])[, c("inhibitor", "family")], rows = NULL)
d <- dev.off()

pdf(paste0(plot.dir, "/num-targets-vs-range-outliers.pdf"))
grid.table(subset(drug.families, inhibitor %in% ki.tbl$inhibitor[flag])[, c("inhibitor", "family")], rows = NULL)
d <- dev.off()

lmf <- lm(formula("mean.correlation ~ 0 + training.range + Number.of.targets"), data = ki.tbl)
summary(lmf)

train.models.using.targets <- function(auc.matrix, expr.matrix, drug.targets) {
    colnames(auc.matrix) <- make.names(colnames(auc.matrix))
    auc.matrix <- t(auc.matrix)
    expr.matrix <- t(expr.matrix)
    common.samples <- intersect(rownames(expr.matrix), rownames(auc.matrix))
    expr.matrix <- expr.matrix[common.samples,]
    auc.matrix <- auc.matrix[common.samples,]
    drugs <- names(drug.targets)
    drugs <- drugs[drugs %in% colnames(auc.matrix)]
    names(drugs) <- drugs
    fits <- llply(drugs, .parallel = FALSE,
                 .fun = function(drug) {
                     y.var <- drug
                     xs <- drug.targets[[drug]]
                     xs <- xs[xs %in% colnames(expr.matrix)]
                     data.mat <- cbind(auc.matrix[,drug,drop=F], expr.matrix[,xs,drop=FALSE])
                     fit <- lm(formula(paste0(y.var, ' ~ 0 + ', paste0(xs, collapse=" + "))), data=as.data.frame(data.mat))
                     fit
                 })
    fits
}

test.models <- function(fits, auc.matrix, expr.matrix) {
    colnames(auc.matrix) <- make.names(colnames(auc.matrix))
    auc.matrix <- t(auc.matrix)
    expr.matrix <- t(expr.matrix)
    common.samples <- intersect(rownames(expr.matrix), rownames(auc.matrix))
    expr.matrix <- expr.matrix[common.samples,]
    auc.matrix <- auc.matrix[common.samples,]
    drugs <- names(fits)
    drugs <- drugs[drugs %in% colnames(auc.matrix)]
    names(drugs) <- drugs
    res <- llply(drugs, .parallel = TRUE,
                 .fun = function(drug) {
                     y.var <- drug
                     #xs <- symbol.drug.targets[[drug]]
                     #xs <- xs[xs %in% colnames(expr.matrix)]
                     pred <- predict(fits[[drug]], newdata = as.data.frame(expr.matrix))
                     data.frame(pred = pred, obs = as.numeric(auc.matrix[,drug,drop=T]), sample = rownames(auc.matrix))
                })
    res
}

cat("Stopping here\n")
stop("stop")

fits <- train.models.using.targets(t(aucs[['training']]), symbol.exprs[['training']], symbol.drug.targets)

res <- test.models(fits, t(aucs[['validation']]), symbol.exprs[['validation']])
cors <- ldply(res, .fun = function(df) {
    data.frame(cor.p = cor(df$pred, df$obs, use='pairwise.complete.obs'),
               cor.s = cor(df$pred, df$obs, use='pairwise.complete.obs', method='spearman'))
    })
colnames(cors)[1] <- "inhibitor"

sims.long <- melt(structural.similarity.mat)
colnames(sims.long) <- c("inhibitor1", "inhibitor2", "similarity")
sims.long <- merge(sims.long, drug.clusters, by.x = c("inhibitor1"),
                   by.y = c("inhibitor"))
sims.long <- merge(sims.long, drug.clusters, by.x = c("inhibitor2"),
                   by.y = c("inhibitor"))
sims.long$samecluster <- sims.long$cluster.x == sims.long$cluster.y
sims.long$intracluster <-
    unlist(apply(sims.long[, c("cluster.x", "cluster.y")],
                 1, function(row) {
                     if(row[1] == row[2]) {
                         return(paste0("intracluster ", row[1]))
                     }
                     return("different\nclusters")
                 }))


g <- ggplot()
g <- g + geom_boxplot(data = sims.long, aes(x = intracluster, y = similarity))
g <- g + xlab("")

png(paste0(plot.dir, "/sc1-training-clusters-vs-similarity.png"))
print(g)
d <- dev.off()

pdf(paste0(plot.dir, "/sc1-training-clusters-vs-similarity.pdf"))
print(g)
d <- dev.off()

sims.long <- melt(structural.similarity.mat)
colnames(sims.long) <- c("inhibitor1", "inhibitor2", "similarity")
sims.long <- merge(sims.long, drug.clusters, by.x = c("inhibitor1"),
                   by.y = c("inhibitor"))
sims.long <- merge(sims.long, drug.clusters, by.x = c("inhibitor2"),
                   by.y = c("inhibitor"))
sims.long$samecluster <- sims.long$cluster.x == sims.long$cluster.y
sims.long$intracluster <-
    unlist(apply(sims.long[, c("cluster.x", "cluster.y")],
                 1, function(row) {
                     if(row[1] == row[2]) {
                         return(paste0("intracluster ", row[1]))
                     }
                     return("different\nclusters")
                 }))


g <- ggplot()
g <- g + geom_boxplot(data = sims.long, aes(x = intracluster, y = similarity))
g <- g + xlab("")

png(paste0(plot.dir, "/sc1-training-clusters-vs-similarity.png"))
print(g)
d <- dev.off()

pdf(paste0(plot.dir, "/sc1-training-clusters-vs-similarity.pdf"))
print(g)
d <- dev.off()

drugs.by.cluster <- dlply(drug.clusters, .variables = c("cluster"), .fun = function(df) as.character(df$inhibitor))

## Take the mean drug/gene correlation across all drugs in a cluster
cluster.gene.correlation.profiles <-
    llply(drugs.by.cluster,
          .fun = function(drugs.in.cluster) {
              rowMeans(sig.cor.mat.all[, drugs.in.cluster])
          })

## Perform GSEA on mean drug/gene correlation profile of each cluster

kegg <- as.data.frame(msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG"))
kegg.pathways <- dlply(kegg, .variables = c("gs_name"), .fun = function(df) as.character(df$human_gene_symbol))

go.bp <- as.data.frame(msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"))
go.bp.terms <- dlply(go.bp, .variables = c("gs_name"), .fun = function(df) as.character(df$human_gene_symbol))

pathway.dbs <- list("kegg" = kegg.pathways, "go.bp" = go.bp.terms)

gsea.res <-
    llply(cluster.gene.correlation.profiles,
          .fun = function(profile) {
              stats <- profile
              llply(pathway.dbs,
                    .fun = function(pathway.db) {
                        msz <- 10
                        res <- fgsea(pathways = pathway.db,
                                     stats = stats,
                                     minSize = msz,
                                     maxSize = 500,
                                     nperm = 10000)
                        o <- order(res$padj)
                        res[o, ]
                    })
          })


drugs.by.cluster <- lapply(drugs.by.cluster, sort)

res <- as.data.frame(stri_list2matrix(drugs.by.cluster))
colnames(res) <- paste0("Cluster ", 1:ncol(res))

file <- "cluster-pathway-enrichments.xlsx"
sheet <- "clusters"

wb <- createWorkbook()

## write.xlsx(res, file = file, sheetName = sheet, append = FALSE)
addWorksheet(wb, sheetName = sheet)
writeData(wb, sheet = sheet, res)

nms <- names(gsea.res)
names(nms) <- nms
for(nm in nms) {
    print(nm)
    sheet <- paste0("C", nm, " KEGG")
    tbl <- as.data.frame(gsea.res[[nm]][["kegg"]])
    addWorksheet(wb, sheetName = sheet)
    writeData(wb, sheet = sheet, tbl)
##    write.xlsx(tbl, file = file, sheetName = sheet, append = FALSE)    
    sheet <- paste0("C", nm, " GO BP")
    tbl <- as.data.frame(gsea.res[[nm]][["go.bp"]])
##    write.xlsx(tbl, file = file, sheetName = sheet, append = FALSE)        
    addWorksheet(wb, sheetName = sheet)
    writeData(wb, sheet = sheet, tbl)
}

saveWorkbook(wb, file = file, overwrite = TRUE)


save.image(".Rdata")



