library(pacman)
library(synapser)
p_load(openxlsx)

synLogin()

## Download Robert Allaway's Drug-Target Explorer (DTE) links
compound.structure.synId <- "syn25388383"
compound.target.assocations.synId <- "syn25388377"
compound.links.synId <- "syn25388390"
compound.synonyms.synId <- "syn25388380"

compound.structure.tbl <- read.table(synGet(compound.structure.synId)$path, header=TRUE, sep=",")
compound.target.associations.tbl <- read.table(synGet(compound.target.assocations.synId)$path, header=TRUE, sep=",")
compound.links.tbl <- read.table(synGet(compound.links.synId)$path, header=TRUE, sep=",", comment.char="", quote="\"")
chembl.links.tbl <- subset(compound.links.tbl, database == "chembl")
compound.synonyms.tbl <- read.table(synGet(compound.synonyms.synId)$path, header=TRUE, sep=",", comment.char="", quote="\"")
chembl.synonyms.tbl <- subset(compound.synonyms.tbl, grepl(synonym, pattern="CHEMBL"))[, c("inchikey", "synonym")]
my.dup <- function(x) duplicated(x, fromLast = TRUE) | duplicated(x, fromLast = FALSE)
chembl.synonyms.tbl <- chembl.synonyms.tbl[!my.dup(chembl.synonyms.tbl),]

ohsu.drug.anno <- read.table("ohsu-dtex-translation.txt", sep="\t", header=TRUE, comment.char="")
## ohsu.drug.anno <- ohsu.drug.anno[, c("molecule.name", "inhibitor", "smiles", "link", "external_id", "database", "chembl.id")]
ohsu.drug.anno <- ohsu.drug.anno[, c("molecule.name", "inhibitor", "smiles", "database", "chembl.id")]
ohsu.drug.anno <- subset(ohsu.drug.anno, database == "chembl")
ohsu.drug.anno <- na.omit(ohsu.drug.anno)

# Add in the InChiKey
ohsu.drug.anno <- merge(ohsu.drug.anno, chembl.synonyms.tbl, all.x=TRUE, by.x="chembl.id", by.y = "synonym")
ohsu.drug.associations <- merge(ohsu.drug.anno, compound.target.associations.tbl)
ohsu.drug.associations <- subset(ohsu.drug.associations, !is.na(mean_pchembl))

## Assume that the highest pChembl (e.g., corresponding to the least K_d) corresponds to the drug target.
## Note that pChembl is in log10 space
## Ignore anything with 50x less binding affinity.
## Count the number of targets within this range
off.target.offset <- log10(50)

ohsu.drug.targets <-
    ddply(ohsu.drug.associations, .variables = c("inhibitor"),
          .fun = function(df) {
                   max.mean.pchembl <- max(df$mean_pchembl)
                   data.frame(num.genes.assayed = length(unique(df$hugo_gene)), max.mean.pchembl = max.mean.pchembl,
                              num.targets = length(which(df$mean_pchembl >= (max.mean.pchembl - off.target.offset))))
              })

write.table(ohsu.drug.targets, file="ohsu-drug-targets.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

klaeger.anno <- read.xlsx("aan4368_Table_S1.xlsx", sheet=2)

source("structure-mapping.R")

fp.ohsu <- parseInputFingerprint(unique(as.character(ohsu.drug.anno$smiles)))
klaeger.smiles <- na.omit(as.character(klaeger.anno$SMILES.code))
names(klaeger.smiles) <- klaeger.smiles
fp.klaeger <- parseInputFingerprint(unique(klaeger.smiles))
# fp.klaeger <- llply(klaeger.smiles, parseSingleInputFingerprint)

struct.mapping <- map.drug.structures.based.on.fingerprints(fp.ohsu, fp.klaeger, "ohsu.smiles", "klaeger.smiles", tanimoto.cutoff = 0.9)

# Read in the Klaeger CATDS data
klaeger.catds <- read.xlsx("EMS82808-supplement-Table_S6.xlsx", sheet=2)

klaeger.catds.matrix <- read.xlsx("EMS82808-supplement-Table_S6.xlsx", sheet=3)

# There are multiple CATDS per drug (with the protein in the target column used as
# reference / numerator). Select the max CATDS value (i.e., defined
# relative to the highest binding target)
klaeger.catds <- 
  ddply(klaeger.catds, .variables = c("Drug"),
        .fun = function(df) { data.frame(CATDS = max(df$CATDS)) } )

# Take statistics calculated using the "most potent" protein as the reference
klaeger.stats.tabs <- list("most" = 4, "least" = 5, "most.global" = 6)
klaeger.stats <- 
  llply(klaeger.stats.tabs,
        .fun = function(tab) {
                      klaeger.stats <- read.xlsx("EMS82808-supplement-Table_S6.xlsx", sheet=tab, startRow = 5)
                      stat.cols <- c("Number.of.targets", "Gini", "Selectivity", "Partition.index.single", "CATDS.single", "Entropy")
                      stat.cols <- c(stat.cols, colnames(klaeger.stats)[grepl(x=colnames(klaeger.stats), pattern="Gini.")])
                      stat.cols <- c(stat.cols, colnames(klaeger.stats)[grepl(x=colnames(klaeger.stats), pattern="Selectivity.")])
                      stat.cols <- unique(stat.cols)
                      klaeger.stats <- klaeger.stats[, c("Drug", stat.cols)]
                      klaeger.stats
                    })

analysis.folder.synId <- "syn22397014"

reconciled.tbls <- list()
for(nm in names(klaeger.stats)) {
    reconciled <- merge(ohsu.drug.anno[, c("inhibitor", "smiles")], struct.mapping, by.x = "smiles", by.y = "ohsu.smiles")
    reconciled <- merge(reconciled, klaeger.anno[, c("Drug", "SMILES.code")], by.x = "klaeger.smiles", "SMILES.code")
    reconciled <- merge(reconciled, klaeger.catds)
    reconciled <- merge(reconciled, klaeger.stats[[nm]])
    # reconciled <- merge(reconciled, all.tbl)
    reconciled$inhibitor.name <- unlist(lapply(reconciled$inhibitor, make.names))
    reconciled.tbls[[nm]] <- reconciled
    file = paste0("klaeger-drug-harmonization-table-", nm, ".tsv")
    write.table(file=file, reconciled, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

    f <- File(file, parentId = analysis.folder.synId, synapseStore = TRUE)
    cat(paste0("Storing file ", file , " to Synapse\n"))
    synStore(f)
    file.remove(file)

}


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

symbol.ensg.map <- symbols.to.ensg.mapping(colnames(klaeger.catds.matrix))
symbol.ensg.map$symbol <- as.character(symbol.ensg.map$symbol)
symbol.ensg.map$ensg <- as.character(symbol.ensg.map$ensg)

klaeger.catds.targets <- read.xlsx("EMS82808-supplement-Table_S6.xlsx", sheet=2)
klaeger.catds.targets <- merge(klaeger.catds.targets, reconciled.tbls[["most.global"]][, c("Drug", "inhibitor.name", "inhibitor")], by.x="Drug", by.y="Drug")

file = paste0("klaeger-targets.tsv")
write.table(file=file, klaeger.catds.targets, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

f <- File(file, parentId = analysis.folder.synId, synapseStore = TRUE)
cat(paste0("Storing file ", file , " to Synapse\n"))
synStore(f)
file.remove(file)






tbl <- read.xlsx("Supplementary Table 3 Drug Matrices.xlsx", sheet=2)
clusters <- read.table("drug-clusters.tsv", sep="\t", header=TRUE)


## Get the OHSU drug class annotations
#file <- "~/Downloads/41586_2018_623_MOESM3_ESM.xlsx"
#ohsu.drug.metadata <- openxlsx:::read.xlsx(file, sheet = 11)

# file <- "~/Downloads/ctd2/fimm-ohsu-drug-map.tsv"
file <- "fimm-ohsu-drug-map.tsv"
drug.map <- read.table(file, sep="\t", header=TRUE)

drug.map.kinases <- subset(drug.map, grepl(Class.explained, pattern="Kinase"))
inhibs <- clusters$inhibitor
potential.kinases <-
    inhibs[(inhibs %in% drug.map.kinases$OHSU_DRUG_NAME) ||
           !(inhibs %in% drug.map$OHSU_DRUG_NAME)]
potential.kinases <- as.character(potential.kinases)

missing <- potential.kinases[!(potential.kinases %in% colnames(tbl))]

ohsu.drug.anno$chembl.id %in% klaeger.anno$ChEMBL.ID
ohsu.drug.anno$smiles %in% klaeger.anno$SMILES.code

tmp1 <- ohsu.drug.anno[, c("inhibitor", "chembl.id")]
colnames(tmp1) <- c("ohsu.inhibitor", "ChEMBL.ID")
tmp1 <- na.omit(tmp1)
tmp2 <- klaeger.anno[, c("Drug", "ChEMBL.ID")]
colnames(tmp2) <- c("klaeger.Drug", "ChEMBL.ID")
tmp2 <- na.omit(tmp2)
m <- merge(tmp1, tmp2, all=FALSE)

missing <- potential.kinases[!(potential.kinases %in% m$ohsu.inhibitor)]

missing.smiles <- subset(ohsu.drug.anno[, c("inhibitor", "smiles")], inhibitor %in% missing)
missing.smiles <- unique(missing.smiles)
missing.smiles <- na.omit(missing.smiles)

cat("Calculating drug fingerprints for OHSU\n")
fp.ohsu <- parseInputFingerprint(as.character(missing.smiles$smiles))

cat("Calculating drug fingerprints for Kraeger\n")
                                        # fp.klaeger <- parseSingleInputFingerprint(as.character(klaeger.anno$SMILES.code))
klaeger.smiles <- na.omit(as.character(klaeger.anno$SMILES.code))
names(klaeger.smiles) <- klaeger.smiles
fp.klaeger <- parseInputFingerprint(klaeger.smiles)
# fp.klaeger <- llply(klaeger.smiles, parseSingleInputFingerprint)

struct.mapping <- map.drug.structures.based.on.fingerprints(fp.ohsu, fp.klaeger, "ohsu.smiles", "klaeger.smiles", tanimoto.cutoff = 0.9)

if(!is.null(struct.mapping)) {
    reconciled <- merge(ohsu.drug.anno[, c("inhibitor", "smiles")], struct.mapping, by.x = "smiles", by.y = "ohsu.smiles")
    reconciled <- merge(reconciled, klaeger.anno[, c("Drug", "SMILES.code")], by.x = "klaeger.smiles", "SMILES.code")

}


