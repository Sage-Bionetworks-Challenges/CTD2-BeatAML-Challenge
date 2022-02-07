suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(MESS))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(ggrepel))
suppressPackageStartupMessages(p_load(reshape2))

suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1


file <- "~/Downloads/inhibitor_data_points_10_04_2018.txt"

raw.drug.resp <- read.table(file, sep="\t", header=TRUE)

print(head(raw.drug.resp))

## Only keep those with replication -- nah, but note the replicates
raw.drug.resp$replicant <- as.numeric(as.character(raw.drug.resp$replicant))
replicated <- unique(subset(raw.drug.resp, replicant %in% c(2, 3))[, c("lab_id", "inhibitor", "time_of_read")])
print(dim(raw.drug.resp))
raw.drug.resp <- merge(raw.drug.resp, replicated, all.x=FALSE)
print(dim(raw.drug.resp))

## Only keep the most commonly occurring time of read
tbl <- as.data.frame(table(raw.drug.resp$time_of_read))
tor <- as.numeric(as.character(tbl[which.max(tbl$Freq),1]))
raw.drug.resp <- subset(raw.drug.resp, time_of_read == tor)

print(table(raw.drug.resp$run_type))
raw.drug.resp <- subset(raw.drug.resp, run_type == "Standard")

aucs <-
  ddply(raw.drug.resp,
        .variables = c("lab_id", "run_type", "inhibitor", "replicant", "time_of_read"),
	.fun = function(df) {
	         x <- log10(as.numeric(df$well_concentration))
		 y <- as.numeric(df$normalized_viability)
		 if(length(x) != 7) {
		   return(NULL)
		   print(df)
		   stop(paste0("Len(x) = ", length(x)))
		 }
		 y[y < 0] <- 0
		 y[y > 100] <- 100
		 data <- data.frame(x = x, y = y / 100)
		 # mdl <- glm(y ~ 1 + x, data = data, family = binomial(link = "probit"))
		 mdl <- lm(y ~ 1 + x, data = data)
		 auc.probit <- integrate(f = function(x) predict(mdl, newdata = data.frame(x = x)),
		                         lower = min(data$x), upper = max(data$x))
		 data.frame(auc = auc(x, y), auc.probit = auc.probit[["value"]]*100)
	       })

replicated.aucs <- merge(aucs, replicated)

cols <- c("lab_id", "inhibitor", "replicant", "auc.probit")

tmp <- replicated.aucs[, cols]
tmp <- subset(tmp, replicant %in% c(1,2))
tmp2 <- ddply(tmp, .variables = c("lab_id", "inhibitor"),
                   .fun = function(df) {
		    	    if(nrow(df) != 2) { return() }
		            auc1 = subset(df, replicant == 1)$auc.probit
		            auc2 = subset(df, replicant == 2)$auc.probit
			    data.frame(auc1 = auc1, auc2 = auc2)
			  })

g <- ggplot(data = tmp2, aes(x = auc1, y = auc2))
g <- g + xlab("AUC (Rep 1)") + ylab("AUC (Rep 2)")
g <- g + geom_point()
g <- g + geom_smooth(method = "lm")
g <- g + facet_wrap(~ inhibitor, labeller = labeller(inhibitor = label_wrap_gen(10)))
## g + geom_density_2d()

g.cors <- ddply(tmp2, .(inhibitor), summarise, cor = round(cor(auc1, auc2), 2))
g <- g + geom_text(data=g.cors, aes(label=paste("r=", cor, sep="")), x=200, y=20)
ggsave("per-drug-replicate-correlations.pdf")

ggsave("per-drug-replicate-correlations.png")

rmses <- ddply(tmp2, .variables = c("inhibitor"),
               .fun = function(df) {
	       	        rmse <- sqrt(sum((df$auc1 - df$auc2)^2)/nrow(df))
	                data.frame(rmse = rmse, rel.err = rmse / mean(mean(df$auc1), mean(df$auc2)))
		      })

all.errs <- merge(g.cors, rmses)
o <- order(all.errs$cor)
all.errs$inhibitor <- factor(all.errs$inhibitor, levels = all.errs$inhibitor[o])
all.errs.long <- melt(all.errs)

my_labels <- c("rmse" = "RMSD", "rel.err" = "Normalized RMSD", "cor" = "Pearson Correlation")

g <- ggplot(data = all.errs.long, aes(x = inhibitor, y = value))
g <- g + geom_point()
g <- g + facet_wrap(~ variable, nrow = 3, scales = "free_y", labeller = as_labeller(my_labels))
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g <- g + xlab("") + ylab("")
g <- g + theme(plot.margin = theme_get()$plot.margin + unit(c(0,0,0,25),"pt"))
ggsave("per-drug-replicate-errors.pdf")

ggsave("per-drug-replicate-errors.png")

write.table(file = "drug-replicate-aucs.tsv", tmp2, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

# Take mean over replicates
mean.aucs <-
  ddply(aucs,
        .variables = c("lab_id", "run_type", "inhibitor", "time_of_read"),
	.fun = function(df) {
	         data.frame(mean.auc = mean(df$auc), mean.auc.probit = mean(df$auc.probit))
	       })
	       


drug.cors <-
  ddply(subset(replicated.aucs, replicant %in% c(1,2)),
        .variables = c("run_type", "inhibitor", "time_of_read"),
	.fun = function(df) {
                 df1 <- subset(df, replicant == 1)
                 df2 <- subset(df, replicant == 2)
		 m <- merge(df1, df2, by = "lab_id")
		 flag <- !is.na(m$auc.probit.x) & !is.na(m$auc.probit.y)
		 if(length(which(flag)) < 3) { return(NULL) }
		 data.frame(cor.p = cor(m$auc.probit.x, m$auc.probit.y, use = "pairwise.complete.obs", method = "pearson"),
	 	            cor.s = cor(m$auc.probit.x, m$auc.probit.y, use = "pairwise.complete.obs", method = "spearman"))
               })

# Download performance of each group vs each drug

synLogin()

synId <- "syn22267537"
tbl <- synTableQuery(paste0("SELECT * FROM ", synId))
tbl <- as.data.frame(tbl)

## Exclude "gold" (which must mean gold-standard because its pearson and spearman = 1)
tbl <- subset(tbl, team != "gold")
tbl$pearson <- as.numeric(tbl$pearson)
tbl$spearman <- as.numeric(tbl$spearman)

drug.mean.tbl <-
    ddply(tbl, .variables = c("inhibitor"),
          .fun = function(df) data.frame(pearson = mean(df$pearson, na.rm=TRUE), spearman = mean(df$spearman, na.rm=TRUE)))

tbl$inhibitor <- factor(tbl$inhibitor)

drug.cors.vs.perf <- merge(drug.cors, drug.mean.tbl)

source("plotting_utils.R")

g1 <- plot.correlation(drug.cors.vs.perf$cor.p, drug.cors.vs.perf$pearson, labels = drug.cors.vs.perf$inhibitor, display.pval = TRUE, display.r2 = FALSE)
if(FALSE) {
g1 <- ggplot(data = drug.cors.vs.perf, aes(x = cor.p, y = pearson, label = inhibitor))
g1 <- g1 + geom_point()
g1 <- g1 + geom_text_repel()
g1 <- g1 + geom_smooth(method='lm', formula=y~x)
}
g1 <- g1 + xlab("Replicate Correlation") + ylab("Mean Pearson Correlation\n(Observed vs Predicted Response)")
file <- "performance-vs-replication-correlation.png"
png(file)
print(g1)
d <- dev.off()

file <- "performance-vs-replication-correlation.pdf"
pdf(file)
print(g1)
d <- dev.off()

## Plot dynamic range vs response prediction performance (and label "good" and "poor" drugs)
## these are validation aucs
synId <- "syn21212720"
file <- synGet(synId, downloadFile = TRUE)$path
gt.aucs <- read.table(file, sep=",", header=TRUE, stringsAsFactors = FALSE)

gt.auc.synIds <- list("training" = "syn21212913",
                   "leaderboard" = "syn21671784",
                   "validation" = "syn21212720")
if(FALSE) {

gt.auc.synIds <-
    list("training" = "syn22397379",
         "leaderboard" = "syn22397867",
         "validation" = "syn22398053")
}

read.csv.table.from.synapse <- function(synId, sep=",") {
        file <- synGet(synId)$path
    ret <- suppressMessages(read.table(file, sep=sep, header=TRUE))
}

all.gt.aucs <-
    llply(gt.auc.synIds,
          .fun = function(synId) read.csv.table.from.synapse(synId))

all.gt.aucs <- ldply(all.gt.aucs)

# Compare the results used in the challenge
m <- merge(all.gt.aucs, mean.aucs)
m.replicated <- merge(m, replicated)

p_load(GGally)
# p_load("VineCopula")
# p_load("spatstat")

g <- ggpairs(m.replicated[, c("auc", "mean.auc", "mean.auc.probit")])
g
ggsave("auc-correlations.pdf")

ggsave("auc-correlations.png")


auc.range <-
    ddply(gt.aucs, .variables = c("inhibitor"),
          .fun = function(df) {
              data.frame(range = mad(df$auc))
          })

## Plot AUC dynamic range vs mean prediction performance
m <- merge(auc.range, drug.mean.tbl)

g2 <- plot.correlation(m$range, m$pearson, labels = m$inhibitor, display.pval = TRUE, display.r2 = FALSE)
if(FALSE) {
g2 <- ggplot(data = m, aes(x = range, y = pearson, label = inhibitor))
g2 <- g2 + geom_point()
g2 <- g2 + geom_text_repel()
g2 <- g2 + geom_smooth(method='lm', formula=y~x)
}
g2 <- g2 + xlab("Dynamic Range (MAD)") + ylab("Mean Pearson Correlation\n(Observed vs Predicted Response)")
file <- "performance-vs-dynamic-range.png"
png(file)
print(g2)
d <- dev.off()

file <- "performance-vs-dynamic-range.pdf"
pdf(file)
print(g2)
d <- dev.off()

m.all <- merge(m, drug.cors, all = TRUE)

file <- "inhibitor-stats.tsv"
write.table(file=file, m.all, sep="\t", row.names=FALSE, col.names=TRUE)

g3 <- plot.correlation(m.all$range, m.all$cor.p, labels = m.all$inhibitor, display.pval = TRUE, display.r2 = FALSE)
if(FALSE) {
g3 <- ggplot(data = m.all, aes(x = range, y = cor.p, label = inhibitor))
g3 <- g3 + geom_point()
g3 <- g3 + geom_text_repel()
g3 <- g3 + geom_smooth(method='lm', formula=y~x)
}
g3 <- g3 + xlab("Dynamic Range (MAD)") + ylab("Replicate Correlation")
file <- "replicate-correlation-vs-dynamic-range.png"
png(file)
print(g3)
d <- dev.off()


file <- "replicate-correlation-vs-dynamic-range.pdf"
pdf(file)
print(g3)
d <- dev.off()

suppressPackageStartupMessages(p_load(cowplot))
pg <- plot_grid(plotlist = list(g1, g2, g3), nrow=1, labels = "AUTO")
pg
ggsave("drug-replicate-perf-mad-correlation.png", width = 14.5)
ggsave("drug-replicate-perf-mad-correlation.pdf", width = 14.5)

tmp <- m.all[, c("inhibitor", "range", "pearson", "cor.p", "cor.s")]
colnames(tmp) <- c("inhibitor", "MAD", "performance", "replicate.p", "replicate.s")
if(FALSE) {
for(col in c("MAD", "performance", "replicate.p")) {
  tmp[,col] <- tmp[,col] / max(tmp[, col], na.rm=TRUE)
}
}
melted <- melt(tmp)
g <- ggplot(data = melted)
# g <- g + geom_density(aes(colour = variable, x=value))
g <- g + geom_density(aes(x=value))
g <- g + facet_wrap(~ variable, scales = "free")
g
ggsave("drug-replicate-perf-mad-correlation-density.png", width = 14.5)
ggsave("drug-replicate-perf-mad-correlation-density.pdf", width = 14.5)
