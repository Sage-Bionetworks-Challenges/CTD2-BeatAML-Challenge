suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(tidyverse))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(rcdk))
suppressPackageStartupMessages(p_load(rpubchem))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))

synLogin()

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

aucs <-
    llply(auc.synIds,
          .fun = function(synId) read.csv.table.from.synapse(synId))

all.drugs <-
    sort(unique(unlist(lapply(aucs, function(df) unique(df$inhibitor)))))
all.drugs <- as.character(all.drugs)

## Separate the string "drug1 (drug2)" into "drug1, drug2"
extract.drugs.as.csv <- function(str) {
    if(!grepl(str, pattern="\\(")) {
        return(str)
    }
    pat <- "(\\S+)\\s+\\(([^\\)]+)\\)"
    res <- regexpr(pattern = pat, str, perl=TRUE)
    if(res != 1) {
        stop(paste0("Did not get a match from ", str, "\n"))
    }
    cap.start <- attr(res, "capture.start")
    cap.length <- attr(res, "capture.length")
    indices <- 1:length(cap.start)
    strs <-
        lapply(indices, function(idx)
            substr(str, cap.start[idx], cap.start[idx] + cap.length[idx]-1))
    paste0(strs, collapse=", ")
}

    
drug.anno <- data.frame(inhibitor = all.drugs, stringsAsFactors = FALSE)
drug.anno$synonyms <- unlist(lapply(drug.anno$inhibitor,
                                    extract.drugs.as.csv))

compound.name.to.smiles <- function(compound) {
    df <- get.synonyms(compound)
    if(nrow(df) == 0) {
        warnings(paste0("Could not find synonyms of ", compound, "\n"))
        return(NA)
    }
    df <- df[, c("Name", "CID", "Synonym")]
    cids <- as.numeric(unique(df$CID))
    smiles <-
        unlist(llply(cids,
                     .fun = function(cid) {
                         unique(as.character(get.cid(cid)[, "Canonical.SMILES"]))
                     }))
    smiles
}

reconcile.smiles.across.synonyms <- function(synonyms) {
    compounds <- unlist(strsplit(synonyms, split=",[ ]*"))
    smiles <- unlist(lapply(compounds, compound.name.to.smiles))
    smiles <- smiles[!is.na(smiles)]
    if(length(smiles) == 0) { return(NA) }
    if(any(is.na(smiles))) { return(NA) }
    if(is.logical(smiles) && is.na(smiles)) { return(NA) }
    if(!all(smiles == smiles[1])) {
        print(smiles)
        warnings(paste0("Got multiple SMILES for ", synonyms, "\n"))
        return(NA)
    }
    smiles[1]
}

## drug.anno$smiles <- unlist(lapply(drug.anno$inhibitor, compound.name.to.smiles))
drug.anno$smiles <-
    unlist(llply(drug.anno$synonyms, .fun = reconcile.smiles.across.synonyms))

## Store everything to Synapse
analysis.folder.synId <- "syn22397014"

file <- "compound-annotations.tsv"
write.table(file = file, drug.anno, row.names = FALSE, col.names = TRUE,
            quote = FALSE, sep = "\t")
f <- File(file, parentId = analysis.folder.synId, synapseStore = TRUE)
cat(paste0("Storing file ", file , " to Synapse\n"))
synStore(f)

