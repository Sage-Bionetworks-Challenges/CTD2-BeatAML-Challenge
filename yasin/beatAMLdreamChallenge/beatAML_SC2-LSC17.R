suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(readr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(reshape2))
suppressPackageStartupMessages(p_load(S4Vectors))
suppressPackageStartupMessages(p_load(tibble))
suppressPackageStartupMessages(p_load(limma))
suppressPackageStartupMessages(p_load(survival))
suppressPackageStartupMessages(p_load(survminer))

# Ensure data are downloaded by sourcing ../../analysis/download-challenge-data.R"
# data.dir <- "/Users/whitebr/work/sage/beataml-challenge/Data/training/"
data.dir <- "../../Data/training/"

source("signature-utils.R")

signature.file <- "lsc17-genelist.tsv"
signature.name <- "LSC17"

sig <- read.table(signature.file, sep="\t", header=TRUE)

# We assume that the signature file has columns coef (coefficient) and id (gene symbol)
signature <- sig$coef
names(signature) <- sig$id


# dnaseq <- read_csv(paste0(data.dir, "dnaseq.csv"))
rnaseq <- read_csv(paste0(data.dir, "rnaseq.csv"))
aucs <- read_csv(paste0(data.dir, "aucs.csv"))
clinical_categorical <- read_csv(paste0(data.dir, "clinical_categorical.csv"))
clinical_numerical <- read_csv(paste0(data.dir, "clinical_numerical.csv"))
clinical_categorical_legend <- read_csv(paste0(data.dir, "clinical_categorical_legend.csv"))
colnames(clinical_categorical) <- make.names(colnames(clinical_categorical))
response <- read_csv(paste0(data.dir, "response.csv"))


# drop columns with missing values 
clinical_numerical <- clinical_numerical %>% select(-c("%.Blasts.in.PB","WBC.Count"))

# prepare rnaseq data
rna_log2counts <- as.matrix(rnaseq[,3:dim(rnaseq)[2]])
rownames(rna_log2counts) <- rnaseq$Gene
rna_counts <- round(2^rna_log2counts)
rownames(rna_counts) <- rownames(rna_log2counts)

sample.cols <- colnames(rna_log2counts)
gene.col <- "Symbol"

adj.signature <- 
  adjust.scores.to.reflect.multimappers(signature, as.data.frame(rnaseq), gene.col) 
score <- calculate.score(adj.signature, as.data.frame(rnaseq), gene.col, sample.cols) 
score.df <- data.frame(signature = as.vector(score), lab_id = names(score))
colnames(score.df)[1] <- signature.name

countdata <- as.data.frame(rna_counts)
rnaseq_voom = voom(countdata)$E
topGenes_rnaseq = data.frame(Gene = rownames(rnaseq_voom[
  order(apply(rnaseq_voom,1,mad), decreasing = T)[1:500],]))
sigDE_log2counts <- rna_log2counts[which(rownames(rna_log2counts) 
                                         %in% topGenes_rnaseq$Gene),]
rnaseq_pcDat <- prcomp(t(sigDE_log2counts))

response$vitalStatus[which(response$vitalStatus=="Alive")] <- 0
response$vitalStatus[which(response$vitalStatus=="Dead")] <- 1
response$vitalStatus <- as.numeric(response$vitalStatus)
kmsurvival <- survfit(Surv(overallSurvival, vitalStatus) ~ 1, data=response)
plot(kmsurvival,xlab="time",ylab="Survival probability")

response_data <- response %>% 
  inner_join(clinical_categorical, by = "lab_id") %>% 
  inner_join(clinical_numerical, by = "lab_id") %>% 
  inner_join(score.df, by = "lab_id") %>% 
  inner_join(aucs %>% group_by(lab_id) %>% 
               summarise(mean_auc = mean(auc)), by = "lab_id") %>% 
  inner_join(as.data.frame(rnaseq_pcDat[["x"]]) %>% 
               select(1:min(5,dim(rnaseq_pcDat[["x"]])[2])) %>% 
               mutate(lab_id=rownames(rnaseq_pcDat[["x"]])), by = "lab_id")

rownames(response_data)=response_data$lab_id
response_data <- response_data %>% select(-c("lab_id"))

covars <- signature.name
multiv_formula <- as.formula(paste('Surv(overallSurvival, vitalStatus) ~', paste(covars, collapse=" + ")))

print(multiv_formula)
coxph.fit <- coxph(multiv_formula, data = response_data)
survConcordance(Surv(overallSurvival, vitalStatus) ~ predict(coxph.fit, response_data), response_data)

saveRDS(coxph.fit, file=paste0("coxph-fit-sig-", signature.name, ".rds"))

ggforest(coxph.fit, data=response_data, fontsize = 1.05)
ggsave(paste0("sc2-forest-sig-", signature.name, ".png"), width = 14)