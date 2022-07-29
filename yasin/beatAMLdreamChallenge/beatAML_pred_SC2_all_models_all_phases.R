library(readr)
library(dplyr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(limma)
library(survival)

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

source("signature-utils.R")



read.sig <- function(signature.file) {
  sig <- read.table(signature.file, sep="\t", header=TRUE)

  # We assume that the signature file has columns coef (coefficient) and id (gene symbol)
  signature <- sig$coef
  names(signature) <- sig$id
  signature
}

lsc17.signature.name <- "LSC17"
lsc17.signature.file <- paste0(lsc17.signature.name, "-genelist.tsv")
lsc17.signature <- read.sig(lsc17.signature.file)

gentles.signature.name <- "Gentles"
gentles.signature.file <- paste0(gentles.signature.name, "-genelist.tsv")
gentles.signature <- read.sig(gentles.signature.file)

wagner.signature.name <- "Wagner"
wagner.signature.file <- paste0(wagner.signature.name, "-genelist.tsv")
wagner.signature <- read.sig(wagner.signature.file)

wang.signature.name <- "Wang"
wang.signature.file <- paste0(wang.signature.name, "-genelist.tsv")
wang.signature <- read.sig(wang.signature.file)

elsayed.signature.name <- "Elsayed"
elsayed.signature.file <- paste0(elsayed.signature.name, "-genelist.tsv")
elsayed.signature <- read.sig(elsayed.signature.file)

library(survminer)



apply.model <- function(model.name, phase, response_data_t) {
  # load the model fit
  my_coxph <- readRDS(file = paste0(rds_dir,model.name,".rds"))

  #coxph.pred <- survfit(my_coxph, newdata = response_data_t)
  coxph.pred <- predict(my_coxph, response_data_t)

  output_df <- data.frame(
    lab_id = names(coxph.pred), #names(summary(coxph.pred)$table[,"median"]),
    survival = -1 * coxph.pred, #summary(coxph.pred)$table[,"median"], 
    row.names=NULL)
  output_df$survival[which(is.na(output_df$survival))]=0

  out_file=paste0(output_dir,model.name,"_predictions_response_", phase, ".csv")

  write_csv(output_df, out_file)

}



apply.models.to.phase <- function(phase) {

input_dir <- paste0("../../Data/", phase, "/")

dnaseq_t <- read_csv(paste0(input_dir,"dnaseq.csv"))
rnaseq_t <- read_csv(paste0(input_dir,"rnaseq.csv"))
aucs <- read_csv(paste0(input_dir, "aucs.csv"))
clinical_categorical_t <- read_csv(paste0(input_dir,"clinical_categorical.csv"))
clinical_numerical_t <- read_csv(paste0(input_dir,"clinical_numerical.csv"))
clinical_categorical_legend_t <- read_csv(paste0(input_dir,"clinical_categorical_legend.csv"))
colnames(clinical_categorical_t) <- make.names(colnames(clinical_categorical_t))
response <- read_csv(paste0(input_dir, "response.csv"))

response$vitalStatus[which(response$vitalStatus=="Alive")] <- 0
response$vitalStatus[which(response$vitalStatus=="Dead")] <- 1
response$vitalStatus <- as.numeric(response$vitalStatus)

# drop columns with missing values
clinical_numerical_t <- clinical_numerical_t %>% select(-c("%.Blasts.in.PB","WBC.Count"))

# prepare rnaseq data
rna_log2counts_t <- as.matrix(rnaseq_t[,3:dim(rnaseq_t)[2]])
rownames(rna_log2counts_t) <- rnaseq_t$Gene
rna_counts_t <- round(2^rna_log2counts_t)
rownames(rna_counts_t) <- rownames(rna_log2counts_t)

sample.cols <- colnames(rna_log2counts_t)
gene.col <- "Symbol"

expr.df <- as.data.frame(rna_counts_t)
expr.df[, gene.col] <- rnaseq_t[, gene.col]

# rnaseq cpm sum of total check
print(colSums(2^rna_log2counts_t))

# variance stabilizing transformation
#vstcounts_t <- vst(rna_counts_t, blind = FALSE)
#rownames(vstcounts_t) <- rownames(rna_log2counts_t)
#rna_log2counts_t <- vstcounts_t

# prepare input for predict
countdata_t <- as.data.frame(rna_counts_t)
rnaseq_voom_t = voom(countdata_t)$E
topGenes_rnaseq_t = data.frame(Gene = rownames(rnaseq_voom_t[
  order(apply(rnaseq_voom_t,1,mad), decreasing = T)[1:500],]))
sigDE_log2counts_t <- rna_log2counts_t[which(rownames(rna_log2counts_t) 
                                         %in% topGenes_rnaseq_t$Gene),]
rnaseq_pcDat_t <- prcomp(t(sigDE_log2counts_t))

response_data_t <- response %>%
  inner_join(clinical_categorical_t, by = "lab_id") %>% 
  inner_join(clinical_numerical_t, by = "lab_id") %>% 
  inner_join(as.data.frame(rnaseq_pcDat_t[["x"]]) %>% 
               select(1:min(5,dim(rnaseq_pcDat_t[["x"]])[2])) %>% 
               mutate(lab_id=rownames(rnaseq_pcDat_t[["x"]])), by = "lab_id")
library(plyr)
aucs.z <- melt(t(scale(t(acast(aucs, inhibitor ~ lab_id)), center=TRUE, scale=TRUE)))
colnames(aucs.z) <- c("inhibitor", "lab_id", "auc")
mean.aucs <- ddply(aucs.z, .variables = c("lab_id"), .fun = function(df) data.frame(mean_auc = mean(df$auc, na.rm=TRUE)))
## mean.aucs <- ddply(aucs, .variables = c("lab_id"), .fun = function(df) data.frame(mean_auc = mean(df$auc, na.rm=TRUE)))

tmp <- as.data.frame(t(rna_log2counts_t)) %>% mutate(lab_id = rownames(t(rna_log2counts_t)))
response_data_t <- response_data_t %>% inner_join(tmp, by = "lab_id")


response_data_t <- merge(response_data_t, mean.aucs)

sample.cols <- unique(response_data_t$lab_id)
expr.df <- expr.df[, c(gene.col, sample.cols)]

lsc17.score.df <- compute.lsc17.score(expr.df, lsc17.signature, lsc17.signature.name, gene.col, sample.cols)
#gentles.score.df <- compute.gentles.score(expr.df, gentles.signature, gentles.signature.name, gene.col, sample.cols)
wagner.score.df <- compute.score.in.min.max.cpm.space(expr.df, wagner.signature, wagner.signature.name, gene.col, sample.cols)
wang.score.df <- compute.score.in.z.cpm.space(expr.df, wang.signature, wang.signature.name, gene.col, sample.cols)
elsayed.score.df <- compute.score.in.log2.cpm.space(expr.df, elsayed.signature, elsayed.signature.name, gene.col, sample.cols)


response_data_t <- response_data_t %>%
  inner_join(lsc17.score.df, by = "lab_id") %>% 
#  inner_join(gentles.score.df, by = "lab_id") %>% 
  inner_join(wagner.score.df, by = "lab_id") %>% 
  inner_join(wang.score.df, by = "lab_id") %>%
  inner_join(elsayed.score.df, by = "lab_id") 

if(nrow(response_data_t) != nrow(lsc17.score.df)) { stop("Wrong number of patients in LSC17 model\n") }
if(nrow(response_data_t) != nrow(wang.score.df)) { stop("Wrong number of patients in wang model\n") }
if(nrow(response_data_t) != nrow(wagner.score.df)) { stop("Wrong number of patients in wagner model\n") }
if(nrow(response_data_t) != nrow(elsayed.score.df)) { stop("Wrong number of patients in elsayed model\n") }

rownames(response_data_t)=response_data_t$lab_id
response_data_t <- response_data_t %>% select(-c("lab_id"))

pred.dfs <- list()
pred.dfs[[paste0(lsc17.signature.name, "-score")]] <- lsc17.score.df
pred.dfs[[paste0(wagner.signature.name, "-score")]] <- wagner.score.df
pred.dfs[[paste0(wang.signature.name, "-score")]] <- wang.score.df
pred.dfs[[paste0(elsayed.signature.name, "-score")]] <- elsayed.score.df

for(nm in names(pred.dfs)) {
  df <- pred.dfs[[nm]]
  colnames(df) <- c("survival", "lab_id")
  df <- df[, c("lab_id", "survival")]
  # These are _risk_ scores, so flip the sign -> such that high risk, means low/negative survival 
  # (as above for the coxph models)
  # You can confirm this in the coxph models based on the scores that the signs of the score coefficient
  # are positive -- _except_ for Wang
  factor <- -1
  if(nm == paste0(wang.signature.name, "-score")) { factor <- 1 }
  df$survival <- factor * df$survival
  out_file=paste0(output_dir,nm,"_predictions_response_", phase, ".csv")  
  write_csv(df, out_file)
}

for(model in models) {
  print(model)
  apply.model(model, phase, response_data_t)
}



}

print("Training")
apply.models.to.phase("training")
print("Validation")
apply.models.to.phase("validation")
print("Leaderboard")
apply.models.to.phase("leaderboard")
print("Concat")
apply.models.to.phase("concatentated")
print("Done")
