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

score <- score.sc2(m, "survival", "vitalStatus", "overallSurvival")

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

o <- order(aucs$auc, decreasing = TRUE)
aucs <- aucs[o,]
head(aucs)


