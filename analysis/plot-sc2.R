library(tidyverse)
library(reticulate)
library(challengescoring)
library(challengerutils)
library(reshape2)
library(pacman)

# Synapse setup -- using reticulate because synasper has problems.
## use_condaenv("py37")
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

# Rename the columns with the team/participant names
sc2.names <- sapply(submissions2$submitterId, function(sub) {
  name <- tryCatch({
    syn$getUserProfile(sub)$userName
  }, error = function(err) {
    syn$getTeam(sub)$name
  })
  return(name)
})
colnames(sub.mat) <- sc2.names

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

auc <- evaluation_queue_query(
  "select auc from evaluation_9614422 where round == 'final_sub' ORDER BY concordance_index DESC"
)$auc %>%
  as.numeric()

## vchung is in here as the baseline
flag <- colnames(boot.sc2) == "vchung"
colnames(boot.sc2)[flag] <- "Baseline"
names(auc) <- colnames(boot.sc2)

auc.df <- data.frame(team = names(auc), auc = as.numeric(auc))

lvls <- rev(colnames(boot.sc2))

scores <- melt(boot.sc2)
colnames(scores) <- c("boot", "team", "ci")
scores$team <- factor(scores$team, levels = lvls)
auc.df$team <- factor(auc.df$team, levels = lvls)

source("geom_boxplotMod.R")

## Plot boxplot
g1 <- ggplot(data = scores, aes_string(x = "team", y = "ci"))
g1 <- g1 + geom_boxplotMod(fill = "#56B4E9")
g1 <- g1 + coord_flip()
g1 <- g1 + xlab("Method")
g1 <- g1 + ylab("Concordance Index")
g1 <- g1 + theme(text = element_text(size=18), title = element_text(size = 20),
                 axis.text.x = element_text(angle = 45, hjust = 1),
		 axis.title.y = element_blank())

## Plot barplot
g2 <- ggplot(data = auc.df)
g2 <- g2 + geom_col(aes_string(x = "team", y = "auc"), fill = "#E69F00")
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