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

read.csv.from.synapse <- function(synId) {
    if(use.reticulate) {
        file <- syn$get(synId)$path
    } else {
        file <- synGet(synId)$path
    }
    ret <- suppressMessages(read_csv(file))
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
validation.result.tbl <- subset(validation.result.tbl, team != "gold")
validation.result.tbl$pearson <- as.numeric(validation.result.tbl$pearson)

cat(paste0(length(unique(validation.result.tbl$team)), " teams submitted to SC1\n")) 

drug.mean.validation.result.tbl <-
    ddply(validation.result.tbl, .variables = c("inhibitor"),
          .fun = function(df) data.frame(pearson = mean(df$pearson, na.rm=TRUE)))


expr.synIds <- list("training" = "syn21212911",
                    "leaderboard" = "syn21671825",
                    "validation" = "syn21212718")

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
              rnaseq <- suppressMessages(read_csv(rnaseq_file))
              rnaseq
          })

## Read in the intermediate results generated by sc1-analysis-preprocessing.R

gene.symbol.expr.cor.synIds <-
    list("training" = "syn22397862",
         "leaderboard" = "syn22398045",
         "validation" = "syn22398135")

gene.symbol.expr.cors <-
    llply(gene.symbol.expr.cor.synIds,
          .fun = function(synId) read.csv.from.synapse(synId))

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
mean.drug.predictions <- drug.mean.validation.result.tbl$pearson
names(mean.drug.predictions) <- rownames(drug.mean.validation.result.tbl)

all.drugs <- unique(c(all.drugs, "GRD"))
all.drugs <- all.drugs[!(all.drugs %in% drugs.to.exclude)]

drug.cor.mat <- cor(renamed.aucs[["training"]], use = "pairwise.complete.obs")

sig.cor.mat <- as.data.frame(sig.cor.mat)
sig.cor.mat.all <- as.data.frame(sig.cor.mat.all)
drug.cor.mat <- as.data.frame(drug.cor.mat)

training.auc.ranges <- auc.ranges[["training"]]

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
    
}    
drug.family.mat <- drug.family.mat[, all.drugs]
sig.cor.mat <- sig.cor.mat[, all.drugs]
drug.cor.mat <- drug.cor.mat[all.drugs, all.drugs]
sig.cor.mat.all <- sig.cor.mat.all[all.drugs, all.drugs]
mean.drug.predictions <- mean.drug.predictions[all.drugs]
training.auc.ranges <- training.auc.ranges[all.drugs]
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
png("sc1-training-data-vs-validation-performance.png")
sv <- structural.similarity.mat

structural.similarity.mat <- sv
structural.similarity.mat <- as.matrix(structural.similarity.mat )
diag(structural.similarity.mat) <- NA

file <- "multi_uni_plot_data.csv"
plotData <- read.table(file, sep=",", header=TRUE)
multi_vs_uni_variate <- plotData$multi - plotData$uni
names(multi_vs_uni_variate) <- plotData$drug

##cutoff <- 1
## flag <- !is.na(structural.similarity.mat) &
##     structural.similarity.mat > cutoff
## structural.similarity.mat[flag] <- cutoff
ht <- make.plot.dendro(drug.family.mat[, drugs],
                       num.sig.genes[drugs],
		       multi_vs_uni_variate[drugs],
                       mean.drug.predictions[drugs],
                       training.auc.ranges[drugs],
                       drug.cor.mat[drugs, drugs],
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

drug.range.tbl <- data.frame(inhibitor = names(training.auc.ranges), range = unname(training.auc.ranges))
mean.perf.tbl <- data.frame(inhibitor = names(mean.drug.predictions), mean.correlation = unname(mean.drug.predictions))

multi.vs.uni.df <- data.frame(inhibitor = names(multi_vs_uni_variate), multi.vs.uni = as.vector(multi_vs_uni_variate))

all.tbl <- merge(drug.clusters, drug.range.tbl, all.x = TRUE)
all.tbl <- merge(all.tbl, mean.perf.tbl, all.x = TRUE)
all.tbl <- merge(all.tbl, multi.vs.uni.df, all.x = TRUE)
all.tbl <- all.tbl[!(all.tbl$inhibitor == "GRD"),]
all.tbl$cluster <- factor(all.tbl$cluster, levels = sort(unique(all.tbl$cluster)))

# Plot multi vs uni for each cluster

g1 <- ggplot(data = all.tbl)
g1 <- g1 + geom_boxplot(aes(x = cluster, y = mean.correlation))
g1 <- g1 + ylab("Mean Pearson Correlation\n(over Teams)")

g2 <- ggplot(data = all.tbl)
g2 <- g2 + geom_boxplot(aes(x = cluster, y = range))
g2 <- g2 + ylab("Dynamic Range\n(AUC MAD over Samples)")

g3 <- ggplot(data = all.tbl)
g3 <- g3 + geom_boxplot(aes(x = cluster, y = multi.vs.uni))
g3 <- g3 + ylab("Multi- vs Uni-variate Performance")


png("sc1-training-clusters-vs-validation-performance.png")
grid.arrange(g1,g2,g3)
d <- dev.off()

png("sc1-training-clusters-vs-multi-vs-uni-performance.png")
print(g3)
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

png("sc1-training-clusters-vs-similarity.png")
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

png("sc1-training-clusters-vs-similarity.png")
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



