suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(readr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(reshape2))
suppressPackageStartupMessages(p_load(S4Vectors))
suppressPackageStartupMessages(p_load(tibble))
suppressPackageStartupMessages(p_load(limma))
suppressPackageStartupMessages(p_load(survival))
suppressPackageStartupMessages(p_load(survminer))

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--optimize-coeffs"), action="store_true",
                default=FALSE,
                help="(Re-)optimize the published coefficients in a signature")
)

descr <- "\
Train a model
"

parser <- OptionParser(usage = "%prog [options]", option_list=option_list, description=descr)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
recompute.coefs <- opt$`optimize-coeffs`

if(recompute.coefs) {
  cat("Recomputing coefficients\n")
} else {
  cat("Using published coefficients\n")
}

# Ensure data are downloaded by sourcing ../../analysis/download-challenge-data.R"
# data.dir <- "/Users/whitebr/work/sage/beataml-challenge/Data/training/"
data.dir <- "../../Data/training/"

source("signature-utils.R")

signature.name <- "Wagner"
signature.file <- paste0(signature.name, "-genelist.tsv")

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

expr.df <- as.data.frame(rna_counts)
expr.df[, gene.col] <- rnaseq[, gene.col]

# See "Statistical analyses" of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6482359/
# which says that
# Gene-expression values were normalized using the min-max approach. 
# (though not immediately clear what the original, non min-max'ed, scale was).
score.df <- compute.score.in.min.max.cpm.space(expr.df, signature, signature.name, gene.col, sample.cols)

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

covars <- signature.name
model.name <- signature.name
if(recompute.coefs) {
  model.name <- paste0(model.name, "-recompute")
  covars <- as.data.frame(rnaseq)[rnaseq$Symbol %in% names(signature), "Gene"]
  tmp <- as.data.frame(t(rna_log2counts)) %>% mutate(lab_id = rownames(t(rna_log2counts)))
  tmp <- tmp[, c("lab_id", covars)]
  response_data <- response_data %>% inner_join(tmp, by = "lab_id")
}

rownames(response_data)=response_data$lab_id
response_data <- response_data %>% select(-c("lab_id"))


multiv_formula <- as.formula(paste('Surv(overallSurvival, vitalStatus) ~', paste(covars, collapse=" + ")))

print(multiv_formula)
coxph.fit <- coxph(multiv_formula, data = response_data)
survConcordance(Surv(overallSurvival, vitalStatus) ~ predict(coxph.fit, response_data), response_data)

saveRDS(coxph.fit, file=paste0("coxph-fit-sig-", model.name, ".rds"))

ggforest(coxph.fit, data=response_data, fontsize = 1.05)
ggsave(paste0("sc2-forest-sig-", model.name, ".png"), width = 14)