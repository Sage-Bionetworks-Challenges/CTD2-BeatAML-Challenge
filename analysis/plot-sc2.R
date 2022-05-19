# library(tidyverse)
library(reticulate)
library(challengescoring)
library(challengerutils)
library(reshape2)
library(pacman)
library(ggplot2)

# Synapse setup -- using reticulate because synasper has problems.
# Follow directions for installed challengerutils here: https://github.com/Sage-Bionetworks/challengerutils
# conda create  -c sebp -c conda-forge -n challenge python=3.7 lifelines scikit-survival scipy
## # conda create -n challenge python=3.7 lifelines scikit-survival scipy
# pip install challengeutils
## conda install -c conda-forge lifelines
use_condaenv("challenge")
# remotes::install_github("Sage-Bionetworks/challengerutils")
synapseclient <- reticulate::import('synapseclient')
utils <- reticulate::import('challengeutils.utils')

syn <- synapseclient$Synapse()
syn$login()

# challengerutils requires a login inside its package, so logging in again
challengerutils::syn_login()

## sc2 <- reticulate::import_from_path('sc2_utils', path = "../Docker")
source_python("../Docker/sc2_utils.py")

gold.sc2 <- read.csv(syn$get("syn21297965")$path)
head(gold.sc2)

submissions2 <- evaluation_queue_query(
  "select submitterId, prediction_fileid from evaluation_9614422 where round == 'final_sub' ORDER BY concordance_index DESC"
)

# Get all predictions and put into a single dataframe.
sub.mat <- lapply(submissions2$prediction_fileid, function(sub) {
  pred <- read.csv(syn$get(sub)$path)
  
  # order the predictions to match the goldstandard's order
  temp <- pred[match(gold.sc2$lab_id, pred$lab_id), 2]
  return(temp)
}) %>%
  as.data.frame()

# Add several models based on yasin's approach
# model.dir <- "/Users/whitebr/work/sage/beataml-challenge/yasin/beatAMLdreamChallenge/output/"
model.dir <- "../yasin/beatAMLdreamChallenge/output/"

# Ensure data are downloaded by sourcing ../../analysis/download-challenge-data.R"
# data.dir <- "/Users/whitebr/work/sage/beataml-challenge/Data/training/"
data.dir <- "../Data/"
# data.dir <- "/Users/whitebr/work/sage/beataml-challenge/Data/"
validation.response.file <- paste0(data.dir, "/validation/response.csv")
leaderboard.response.file <- paste0(data.dir, "/leaderboard/response.csv")
training.response.file <- paste0(data.dir, "/training/response.csv")

models <-
  list(
       "Base: age" = "coxph-fit-age",
       "ymemari (uncor)" = "coxph-fit-uncor")

models <-
  list(
       "Base: age + mean AUC" = "coxph-fit-age-grd",
       "Base: age" = "coxph-fit-age",
       "ymemari + mean AUC" = "coxph-fit-mean-auc",
# ymemari-rerun is identical to the submitted ymemari results
#       "ymemari-rerun" = "coxph-fit",
       "ymemari - PC5" = "coxph-fit-no-PC5",
       "Base: mean AUC" = "coxph-fit-mean-auc-only",  
       "Base: LSC17" = "coxph-fit-sig-LSC17",     
       "ymemari (uncor)" = "coxph-fit-uncor",
       "ymemari (uncor) + PC5" = "coxph-fit-uncor-with-PC5",
       "ymemari (uncor) + mean AUC" = "coxph-fit-uncor-grd")



for(mdl in models) {
  model.output <- paste0(model.dir, mdl, "_predictions_response.csv")
  pred <- read.csv(model.output)
  
  # order the predictions to match the goldstandard's order
  temp <- pred[match(gold.sc2$lab_id, pred$lab_id), 2]
  sub.mat <- cbind(sub.mat, as.numeric(temp))  
}


# Rename the columns with the team/participant names
sc2.names <- sapply(submissions2$submitterId, function(sub) {
  name <- tryCatch({
    syn$getUserProfile(sub)$userName
  }, error = function(err) {
    syn$getTeam(sub)$name
  })
  return(name)
})

flag <- sc2.names == "vchung"
sc2.names[flag] <- "Base: all data"

participant.methods <- sc2.names[!flag]
participant.methods <- as.vector(participant.methods)

sc2.names <- c(as.character(sc2.names), names(models))


colnames(sub.mat) <- sc2.names

auc <- evaluation_queue_query(
  "select auc from evaluation_9614422 where round == 'final_sub' ORDER BY concordance_index DESC"
)$auc %>%
  as.numeric()

names(auc) <- sc2.names[1:length(auc)]

validation_survival <- read.table(validation.response.file, sep=",", header=TRUE)
leaderboard_survival <- read.table(leaderboard.response.file, sep=",", header=TRUE)
train_survival <- read.table(training.response.file, sep=",", header=TRUE)

# preds is a matrix whose rows are samples and columns are predictions
ensemble.pred <- function(preds) {
  rnks <- apply(preds, 2, rank)
  rowMeans(rnks)
}


my.auc <-
  unlist(lapply(1:ncol(sub.mat), function(i) scoreSC2_auc_with_r(train_survival, validation_survival, data.frame(survival = sub.mat[,i]))))
names(my.auc) <- colnames(sub.mat)  

ens.pred <- ensemble.pred(sub.mat[,participant.methods])
ens.auc <- scoreSC2_auc_with_r(train_survival, validation_survival, data.frame(survival = ens.pred))
my.auc["Ensemble"] <- ens.auc

comp <- data.frame(lab_id = names(auc), auc = as.numeric(auc))
comp <- merge(comp, data.frame(lab_id = names(my.auc), my.auc = as.numeric(my.auc)))

set.seed(8)

# Bootstrap CI 1000 times with random sampling, then compute Bayes.
N <- 1000
bs_indices.sc2 <- matrix(1:nrow(sub.mat), nrow(sub.mat), N) %>%
  apply(2, sample, replace = T)

boot.sc2 <- apply(bs_indices.sc2, 2, function(ind) {
  tmp.gold <- gold.sc2[ind,]
  apply(sub.mat[ind,], 2, function(pred) {
    ## ci <- sc2$scoreSC2_with_r(pred, tmp.gold)
    ci <- scoreSC2_with_r(pred, tmp.gold)
    return(ci)
  })
}) %>%
  t()


bayes.sc2 <- computeBayesFactor(boot.sc2, 1, T) %>%
  as.data.frame()
bayes.sc2

# Apply an ensemble method that is just the mean rank
# NB: intentionally do this _after_ calculating bayes factor

ens.boot.sc2 <- apply(bs_indices.sc2, 2, function(ind) {
  tmp.gold <- gold.sc2[ind,]
  ens.pred <- ensemble.pred(sub.mat[ind,participant.methods])
  ci <- scoreSC2_with_r(ens.pred, tmp.gold)
}) 
boot.sc2 <- cbind(boot.sc2, ens.boot.sc2)
colnames(boot.sc2)[ncol(boot.sc2)] <- "Ensemble"



## vchung is in here as the baseline
flag <- colnames(boot.sc2) == "vchung"
colnames(boot.sc2)[flag] <- "Base: all data"

# lvls <- rev(colnames(boot.sc2))

df <- data.frame(team = colnames(boot.sc2), score = as.numeric(colMeans(boot.sc2)))
o <- order(df$score, decreasing=FALSE)
df <- df[o,]
lvls <- df$team
print(sort(lvls))

scores <- melt(boot.sc2)
colnames(scores) <- c("boot", "team", "ci")
scores$team <- factor(scores$team, levels = lvls)
scores$fill <- "#56B4E9"
ties <- c(rownames(bayes.sc2)[1], rownames(bayes.sc2)[bayes.sc2[,1] < 3], "Ensemble")
scores[scores$team %in% ties,"fill"] <- "#9FE600"

source("geom_boxplotMod.R")

## Plot boxplot
g1 <- ggplot(data = scores, aes_string(x = "team", y = "ci", fill = "fill"))
# g1 <- g1 + geom_boxplotMod(fill = "#56B4E9")
g1 <- g1 + geom_boxplotMod()
g1 <- g1 + coord_flip()
g1 <- g1 + xlab("Method")
g1 <- g1 + ylab("Concordance Index")
g1 <- g1 + theme(text = element_text(size=18), title = element_text(size = 20),
                 axis.text.x = element_text(angle = 45, hjust = 1),
		 axis.title.y = element_blank(),
		 legend.position = "none")

# auc.df <- data.frame(team = names(auc), auc = as.numeric(auc))
auc.df <- data.frame(team = names(my.auc), auc = as.numeric(my.auc))
auc.df <- merge(scores, auc.df, all.x = TRUE)
auc.df$team <- factor(auc.df$team, levels = lvls)
auc.df$fill <- "#E69F00"

write.table(file="scores.tsv", auc.df, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

## Plot barplot
g2 <- ggplot(data = auc.df)
g2 <- g2 + geom_col(aes_string(x = "team", y = "auc"), fill = "#E69F00")
# g2 <- g2 + geom_col(aes_string(x = "team", y = "auc", fill = "fill"))
g2 <- g2 + coord_flip()
g2 <- g2 + xlab("Method")
g2 <- g2 + ylab("AUC")
g2 <- g2 + theme(text = element_text(size=18))    
g2 <- g2 + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1))

suppressPackageStartupMessages(p_load(cowplot)) ## for plot_grid

pg <- plot_grid(g1, g2, nrow=1, align="h", rel_widths = c(3,0.6))

png("sc2-scores.png", width = 2 * 480)
print(pg)
d <- dev.off()