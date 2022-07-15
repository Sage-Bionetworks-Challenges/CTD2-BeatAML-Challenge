suppressPackageStartupMessages(require(pacman))
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(survival))
suppressPackageStartupMessages(p_load(survminer))

# model.dir <- "/Users/whitebr/work/sage/beataml-challenge/yasin/beatAMLdreamChallenge/output/"
model.dir <- "../yasin/beatAMLdreamChallenge/output/"

# Ensure data are downloaded by sourcing ../../analysis/download-challenge-data.R
data.dir <- "../Data/"

datasets <- c("training", "leaderboard", "validation")
names(datasets) <- datasets

responses <-
  llply(datasets,
        .fun = function(dataset) {
                 read.table(paste0(data.dir, dataset, "/response.csv"), header=TRUE, sep=",")
               })

response.df <- ldply(responses)
colnames(response.df)[1] <- "phase"

plot.kms.across.datasets <- function(response.df, censor.time.days) {

  flag <- response.df$vitalStatus == "Dead"
  response.df[flag,"vitalStatus"] <- 1
  flag <- response.df$vitalStatus == "Alive"
  response.df[flag,"vitalStatus"] <- 0
  response.df$vitalStatus <- as.numeric(response.df$vitalStatus)

  response.df$censoredVitalStatus <- response.df$vitalStatus
  response.df$censoredOverallSurvival <- response.df$overallSurvival
  flag <- response.df$overallSurvival > censor.time.days
  response.df[flag, "censoredOverallSurvival"] <- censor.time.days
  response.df[flag, "censoredVitalStatus"] <- 0

  # os <- survfit( Surv(overallSurvival, vitalStatus) ~ phase, data = response.df)
  os <- survfit( Surv(censoredOverallSurvival, censoredVitalStatus) ~ phase, data = response.df)
  # os <- survfit( Surv(os.time, os.status) ~ 1, data = d.demo)

  sd <- survdiff( Surv(censoredOverallSurvival, censoredVitalStatus) ~ phase, data = response.df)
  sd.pval <- round(1 - pchisq(sd$chisq, length(sd$n) - 1), digits=2)
  print(sd)

  g <- ggsurvplot(os, data = response.df, combine = TRUE, # Combine curves
             risk.table = TRUE,                  # Add risk table
             conf.int = FALSE,                    # Add confidence interval
             conf.int.style = "step",            # CI style, use "step" or "ribbon"
             censor = TRUE,                     # Show censor points
             tables.theme = theme_cleantable(),  # Clean risk table
             palette = "jco")
  g <- g + ggtitle(paste0("BeatAML censored at ", censor.time.days, " p = ", sd.pval))
  g
}

# In this challenge, we used dynamic/cumulative AUC
# (https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.metrics.cumulative_dynamic_auc.html)
# evaluated at 1 year.
# So, let's censor to 1 year.
for(num.years in c(1, 2, 3)) {
  censor.time.days <- 365.25 * num.years
  g <- plot.kms.across.datasets(response.df, censor.time.days)
  png(paste0("plots/", "km-across-phases-censor-", num.years, "-yrs.png"))
  print(g)
  d <- dev.off()
}

