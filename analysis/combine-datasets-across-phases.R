suppressPackageStartupMessages(require(pacman))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(data.table))

# Ensure data are downloaded by sourcing ../../analysis/download-challenge-data.R
data.dir <- "../Data/"

datasets <- c("training", "leaderboard", "validation")
names(datasets) <- datasets

dest.dir <- paste0(data.dir, "concatentated")

concat.files <- c("aucs.csv", "clinical_categorical.csv", "clinical_numerical.csv", "response.csv", "dnaseq.csv")

dir.create(dest.dir)
file.copy(paste0(data.dir, "training/clinical_categorical_legend.csv"), paste0(dest.dir, "/clinical_categorical_legend.csv"))

expr.mats <-
  llply(datasets,
        .fun = function(dataset) as.data.frame(fread(paste0(data.dir, dataset, "/rnaseq.csv"))))


all.cols <- llply(expr.mats, .fun = function(mat) colnames(mat)[!(colnames(mat) %in% c("Gene", "Symbol"))])
all.cols <- Reduce("c", all.cols)
stopif(any(duplicated(all.cols)))

merged.expr.mat <- Reduce(function(...) merge(..., all=T), expr.mats)

merged.expr.mat <- merged.expr.mat[, c("Gene", "Symbol", all.cols)]

fwrite(merged.expr.mat, file = paste0(dest.dir, "/rnaseq.csv"))

for(file in concat.files) {
  dfs <- llply(datasets,
               .fun = function(dataset) as.data.frame(fread(paste0(data.dir, dataset, "/", file))))
  cols <- Reduce(intersect, lapply(dfs, colnames))
  dfs <- llply(dfs, .fun = function(df) df[, cols])
  df <- ldply(dfs, data.frame)
  df <- df[, make.names(cols)]
  colnames(df) <- cols
  fwrite(df, file = paste0(dest.dir, "/", file))
}