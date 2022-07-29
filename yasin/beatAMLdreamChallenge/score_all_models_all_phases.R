suppressPackageStartupMessages(require(pacman))
suppressPackageStartupMessages(p_load(survivalROC))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(data.table))

input_dir = "/input/"
output_dir = "output/"
rds_dir = "/rds/"

rds_dir = "./"

models <- c(
  "coxph-fit",
  "coxph-fit-sig-LSC17",
  "coxph-fit-sig-LSC17-recompute",
  "coxph-fit-sig-Gentles-recompute",
  "coxph-fit-sig-Wagner",
  "coxph-fit-sig-Wagner-recompute",
  "coxph-fit-sig-Wang",
  "coxph-fit-sig-Wang-recompute",
  "coxph-fit-sig-Elsayed",
  "coxph-fit-sig-Elsayed-recompute",
  "coxph-fit-no-PC5",
  "coxph-fit-age-grd",
  "coxph-fit-age",
  "coxph-fit-mean-auc",
  "coxph-fit-mean-auc-only",  
  "coxph-fit-uncor",
  "coxph-fit-uncor-with-PC5",
  "coxph-fit-uncor-grd")

names(models) <- models
gt.file <- "../../Data/concatentated/response.csv"
gt <- as.data.frame(fread(gt.file))
#colnames(gt)[colnames(gt) == "vitalStatus"] <- "truthVitalStatus"
#colnames(gt)[colnames(gt) == "overallSurvival"] <- "truthOverallSurvival"


score.sc2 <- function(tbl, pred.col, status.truth.col, os.truth.col, predict.time = 365) {
  score <- survivalROC(Stime=tbl[,os.truth.col], status=tbl[,status.truth.col], marker=tbl[,pred.col], predict.time = predict.time, method = "KM")
  score$AUC
}

datasets <- c("training", "validation", "leaderboard", "concatentated")
names(datasets) <- datasets

output_dir <- "output/"

phase <- "training"
model <- models[1]
res.file <- paste0(output_dir, "/", model, "_predictions_response_", phase, ".csv")

flag <- gt$vitalStatus == "Dead"
gt[flag,"vitalStatus"] <- 1
flag <- gt$vitalStatus == "Alive"
gt[flag,"vitalStatus"] <- 0
gt$vitalStatus <- as.numeric(gt$vitalStatus)

# score <- score.sc2(m, "survival", "vitalStatus", "overallSurvival")

aucs <-
  ldply(datasets,
        .fun = function(phase) {
	         df <- ldply(models,
		             .fun = function(model) {
                                      res.file <- paste0(output_dir, "/", model, "_predictions_response_", phase, ".csv")
				      res <- as.data.frame(fread(res.file))
				      m <- merge(res, gt, by = "lab_id")
                                      data.frame(auc = score.sc2(m, "survival", "vitalStatus", "overallSurvival"))
                                    })
                 colnames(df)[1] <- "model"
		 df

               })
colnames(aucs)[1] <- "phase"	       

orig.aucs <- aucs

flag <- aucs$model == "coxph-fit"
aucs[flag,"model"] <- "yasin"
flag <- aucs$model == "coxph-fit-no-PC5"
aucs[flag,"model"] <- "yasin-no-PC5"
flag <- aucs$model == "coxph-fit-uncor-no-PC5"
aucs[flag,"model"] <- "yasin-uncor-no-PC5"
flag <- aucs$model == "coxph-fit-uncor"
aucs[flag,"model"] <- "yasin-uncor"

aucs$model <- unlist(lapply(aucs$model, function(str) gsub(str, pattern="coxph-fit-", replacement="")))

flag <- aucs$phase == "concatentated"
aucs[flag,"phase"] <- "concatenated"

aucs$phase <- factor(aucs$phase, levels = c("training", "leaderboard", "validation", "concatenated"))

o <- order(aucs$auc, decreasing = TRUE)
aucs <- aucs[o,]
head(aucs)

valid.aucs <- subset(aucs, phase=="validation")
o <- order(valid.aucs$auc, decreasing=FALSE)
lvls <- valid.aucs[o, "model"]
aucs$model <- factor(aucs$model, levels = lvls)

library(ggplot2)
g <- ggplot(data = aucs)
g <- g + geom_col(aes(x = model, y = auc))
g <- g + ylab("survivalROC AUC")
g <- g + facet_wrap(~phase, nrow=1)
g <- g + coord_flip()

png("survivalROC-over-models-and-phases.png", width = 2 * 480)
print(g)
d <- dev.off()