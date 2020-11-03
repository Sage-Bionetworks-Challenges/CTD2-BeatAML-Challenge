library(tidyverse)
library(reticulate)
library(challengescoring)
library(challengerutils)

# Synapse setup -- using reticulate because synasper has problems.
## use_condaenv("py37")
synapseclient <- reticulate::import('synapseclient')
utils <- reticulate::import('challengeutils.utils')

syn <- synapseclient$Synapse()
syn$login()

# challengerutils requires a login inside its package, so logging in again
challengerutils::syn_login()

sc2 <- reticulate::import_from_path('sc2_utils', path = "../Docker")

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
    ci <- sc2$scoreSC2_with_r(pred, tmp.gold)
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

