suppressPackageStartupMessages(require(pacman))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(data.table))
suppressPackageStartupMessages(p_load(xlsx))

# Ensure data are downloaded by sourcing ../../analysis/download-challenge-data.R
data.dir <- "../Data/"

# # Ensure concatentated data are created by sourcing ../../analysis/combine-datasets-across-phases.R
# concat.dir <- paste0(data.dir, "concatentated")
# challenge.survival <- as.data.frame(fread(paste0(concat.dir, "/", "response.csv")))

datasets <- c("training", "leaderboard", "validation")
names(datasets) <- datasets

responses <-
  llply(datasets,
        .fun = function(dataset) {
                 read.table(paste0(data.dir, dataset, "/response.csv"), header=TRUE, sep=",")
               })

response.df <- ldply(responses)
colnames(response.df)[1] <- "phase"

wave1234.metadata <- read.xlsx("Table_S1.xlsx", sheetIndex = 1)
wave1234.key <- read.xlsx("lab_dbgap_key.xlsx", sheetIndex = 1)

m <- merge(wave1234.key, wave1234.metadata, by=c("dbgap_dnaseq_sample", "dbgap_rnaseq_sample", "dbgap_subject_id"))
stopifnot(nrow(m) == nrow(wave1234.key))
stopifnot(nrow(m) == nrow(wave1234.metadata))

colnames(response.df)[colnames(response.df)=="vitalStatus"] <- "challengeVitalStatus"
colnames(response.df)[colnames(response.df)=="overallSurvival"] <- "challengeOverallSurvival"

response.df.with.metadata <- merge(response.df, m, by.x = c("lab_id"), by.y = c("labId"))
stopifnot(nrow(response.df) == nrow(response.df.with.metadata))

fwrite(response.df.with.metadata, "all-responses-with-metadata.csv")

