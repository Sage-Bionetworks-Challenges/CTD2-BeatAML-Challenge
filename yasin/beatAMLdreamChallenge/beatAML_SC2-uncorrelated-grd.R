library(readr)
library(dplyr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(limma)
library(survival)
library(survminer)
library(plyr)

data.dir <- "/Users/whitebr/work/sage/beataml-challenge/Data/training/"

dnaseq <- read_csv(paste0(data.dir, "dnaseq.csv"))
rnaseq <- read_csv(paste0(data.dir, "rnaseq.csv"))
aucs <- read_csv(paste0(data.dir, "aucs.csv"))
clinical_categorical <- read_csv(paste0(data.dir, "clinical_categorical.csv"))
clinical_numerical <- read_csv(paste0(data.dir, "clinical_numerical.csv"))
clinical_categorical_legend <- read_csv(paste0(data.dir, "clinical_categorical_legend.csv"))
colnames(clinical_categorical) <- make.names(colnames(clinical_categorical))
response <- read_csv(paste0(data.dir, "response.csv"))

# mean.aucs <- ddply(aucs, .variables = c("lab_id"), .fun = function(df) data.frame(meanAUC = mean(df$auc)))
aucs.z <- melt(t(scale(t(acast(aucs, inhibitor ~ lab_id)), center=TRUE, scale=TRUE)))
colnames(aucs.z) <- c("inhibitor", "lab_id", "auc")
mean.aucs <- ddply(aucs.z, .variables = c("lab_id"), .fun = function(df) data.frame(mean_auc = mean(df$auc, na.rm=TRUE)))

# drop columns with missing values 
clinical_numerical <- clinical_numerical %>% select(-c("%.Blasts.in.PB","WBC.Count"))

# prepare rnaseq data
rna_log2counts <- as.matrix(rnaseq[,3:dim(rnaseq)[2]])
rownames(rna_log2counts) <- rnaseq$Gene
rna_counts <- round(2^rna_log2counts)
rownames(rna_counts) <- rownames(rna_log2counts)

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
  inner_join(mean.aucs, by = "lab_id") %>%
#  inner_join(aucs %>% group_by(lab_id) %>% 
#               summarise(mean_auc = mean(auc)), by = "lab_id") %>% 
  inner_join(as.data.frame(rnaseq_pcDat[["x"]]) %>% 
               select(1:min(5,dim(rnaseq_pcDat[["x"]])[2])) %>% 
               mutate(lab_id=rownames(rnaseq_pcDat[["x"]])), by = "lab_id")

rownames(response_data)=response_data$lab_id
response_data <- response_data %>% select(-c("lab_id"))
#coxph(Surv(overallSurvival, vitalStatus) ~ . , data = response_data)  
covariates <- names(response_data)[-c(1,2)]
#below copied from r-bloggers
univ_formulas <- sapply(covariates, function(x) 
  as.formula(paste('Surv(overallSurvival, vitalStatus) ~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = response_data)})
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

cols <- rownames(res[which(as.numeric(res[,4])<0.05),])
if(TRUE) {
if(("ageAtDiagnosis" %in% cols) && ("ageAtSpecimenAcquisition" %in% cols)) {
  cols <- cols[(cols != "ageAtSpecimenAcquisition")]
}
}

multiv_formula <- paste('Surv(overallSurvival, vitalStatus) ~', paste(cols,collapse=" + "))
multiv_formula <- paste0(multiv_formula, " + ", "mean_auc")
multiv_formula <- as.formula(multiv_formula)

coxph.fit <- coxph(multiv_formula, data = response_data)
survConcordance(Surv(overallSurvival, vitalStatus) ~ predict(coxph.fit, response_data), response_data)

## Use stepwise regression to remove correlation
coxph.fit <- step(coxph.fit, direction="backward")

cat("After stepwise regression\n")
print(coxph.fit)

cols <- names(coefficients(coxph.fit))
cols <- unique(c(cols, "mean_auc"))

multiv_formula <- paste('Surv(overallSurvival, vitalStatus) ~', paste(cols,collapse=" + "))
multiv_formula <- paste0(multiv_formula, " + ", "mean_auc")
multiv_formula <- as.formula(multiv_formula)

coxph.fit <- coxph(multiv_formula, data = response_data)
survConcordance(Surv(overallSurvival, vitalStatus) ~ predict(coxph.fit, response_data), response_data)

saveRDS(coxph.fit, file="coxph-fit-uncor-grd.rds")

ggforest(coxph.fit, data=response_data, fontsize = 1.05)
ggsave("sc2-forest-uncor-grd.png", width = 14)

if(FALSE) {
df <- as.data.frame(response_data)
for(i in 1:4) {
for(j in (i+1):5) {
coli <- cols[i]
colj <- cols[j]
cat(paste0(coli, " vs ", colj, "\n"))
print(table(df[,coli], df[,colj]))
}}
}