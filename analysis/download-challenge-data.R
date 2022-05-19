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

suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))

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

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

raw.data.dir <- "../Data/"
dir.create(raw.data.dir)

phases <- c("training", "leaderboard", "validation")
names(phases) <- phases

for(phase in phases) { dir.create(paste0(raw.data.dir, phase)) }

read.csv.from.synapse <- function(synId, ...) {
    if(use.reticulate) {
        file <- syn$get(synId)$path
    } else {
        file <- synGet(synId)$path
    }
    ret <- suppressMessages(read.table(file, ...))
}

read.csv.table.from.synapse <- function(synId, sep=",") {
    if(use.reticulate) {
        file <- syn$get(synId)$path
    } else {
        file <- synGet(synId)$path
    }
    ret <- suppressMessages(read.table(file, sep=sep, header=TRUE))
}

response.synIds <- list("training" = "syn21297969")
response.synIds[["leaderboard"]] <- "syn21671789"
response.synIds[["validation"]] <- "syn21297965"

## Download response data
resp <-
    llply(phases,
          .fun = function(phase) {
              synId <- response.synIds[[phase]]
	      cat("Reading response file ", synId, "\n")
              if(use.reticulate) {
                  resp_file <- syn$get(synId)$path
              } else {
                  resp_file <- synGet(synId, downloadLocation = paste0(raw.data.dir, phase))$path
              }
          })

expr.synIds <- list("training" = "syn21212911")
expr.synIds[["leaderboard"]] <- "syn21671825"
expr.synIds[["validation"]] <- "syn21212718"

## Download expression data
raw.exprs <-
    llply(phases,
          .fun = function(phase) {
              synId <- expr.synIds[[phase]]
	      cat("Reading expr ", synId, "\n")
              if(use.reticulate) {
                  rnaseq_file <- syn$get(synId)$path
              } else {
                  rnaseq_file <- synGet(synId, downloadLocation = paste0(raw.data.dir, phase))$path
              }
          })

auc.synIds <-
    list("training" = "syn21212913", 
         "leaderboard" = "syn21671784",
         "validation" = "syn21212720"
        )

## Download auc data
aucs <-
    llply(phases,
          .fun = function(phase) {
              synId <- auc.synIds[[phase]]
	      cat("Reading auc ", synId, "\n")
              if(use.reticulate) {
                  file <- syn$get(synId)$path
              } else {
                  file <- synGet(synId, downloadLocation = paste0(raw.data.dir, phase))$path
              }
          })


dnaseq.synIds <-
    list("training" = "syn21297970",
         "leaderboard" = "syn21671788",
         "validation" = "syn21297966"
        )

## Download dna data
dna <-
    llply(phases,
          .fun = function(phase) {
              synId <- dnaseq.synIds[[phase]]
	      cat("Reading dna ", synId, "\n")
              if(use.reticulate) {
                  file <- syn$get(synId)$path
              } else {
                  file <- synGet(synId, downloadLocation = paste0(raw.data.dir, phase))$path
              }
          })

clinical.cat.synIds <-
    list("training" = "syn21340776",
         "leaderboard" = "syn21671786",
         "validation" = "syn21340794"
        )

## Download clinical cat data
clinical.cat <-
    llply(phases,
          .fun = function(phase) {
              synId <- clinical.cat.synIds[[phase]]
	      cat("Reading clinical cat ", synId, "\n")
              if(use.reticulate) {
                  file <- syn$get(synId)$path
              } else {
                  file <- synGet(synId, downloadLocation = paste0(raw.data.dir, phase))$path
              }
          })

clinical.cat.legend.synIds <-
    list("training" = "syn21340738",
         "leaderboard" = "syn21671785",
         "validation" = "syn21340795"
        )

## Download clinical cat legend data
clinical.cat.legends <-
    llply(phases,
          .fun = function(phase) {
              synId <- clinical.cat.legend.synIds[[phase]]
	      cat("Reading clinical cat legend ", synId, "\n")
              if(use.reticulate) {
                  file <- syn$get(synId)$path
              } else {
                  file <- synGet(synId, downloadLocation = paste0(raw.data.dir, phase))$path
              }
          })

clinical.numerical.synIds <-
    list("training" = "syn21340778",
         "leaderboard" = "syn21671787",
         "validation" = "syn21340796"
        )


## Download clinical numerical data
clinical.numerical <-
    llply(phases,
          .fun = function(phase) {
              synId <- clinical.numerical.synIds[[phase]]
	      cat("Reading clinical numerical ", synId, "\n")
              if(use.reticulate) {
                  file <- syn$get(synId)$path
              } else {
                  file <- synGet(synId, downloadLocation = paste0(raw.data.dir, phase))$path
              }
          })

