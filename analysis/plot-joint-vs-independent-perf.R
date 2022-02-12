library(pacman)
p_load(ggplot2)
p_load(ggrepel)
p_load(plyr)
suppressPackageStartupMessages(p_load(gridExtra))
suppressPackageStartupMessages(p_load(ggbeeswarm))

plot.dir <- "plots"
dir.create(plot.dir)

# multi_uni_plot_data.csv are Rasmus' original results using the training data
#file <- "multi_uni_plot_data.csv"
#plotData <- read.table(file, sep=",", header=TRUE)
# These are results on the validation data using his code
file <- "rasmus-multi-vs-uni-validation.tsv"
plotData <- read.table(file, sep="\t", header=TRUE, as.is=TRUE)



use.validation.results <- TRUE

if(use.validation.results) {
  plotData$multi <- plotData$spearman.cor.multi
  plotData$uni <- plotData$spearman.cor.uni
  plotData$drug <- plotData$inhibitor
} 
plotData$label <- plotData$drug
plotData$multi.vs.uni <- plotData$multi - plotData$uni
plotData <- plotData[!(plotData$inhibitor == "GRD"),]
#flag <- !is.na(plotData$multi) & !is.na(plotData$uni) & ( (plotData$multi - plotData$uni) > 0.25 )
#plotData[!flag,"label"] <- ""

st <- boxplot.stats(plotData$multi.vs.uni)
all_out_index <- which(plotData$multi.vs.uni %in% st$out)
others <- which(!(plotData$multi.vs.uni %in% st$out))
plotData[others, "label"] <- ""
flag <- plotData$multi.vs.uni %in% st$out[st$out > 0]
multi.drugs <- plotData[flag, "inhibitor"]
flag <- plotData$multi.vs.uni %in% st$out

cat(paste0(length(which(!is.na(plotData$multi.vs.uni) & (plotData$multi.vs.uni > 0))), " of ", length(which(!is.na(plotData$multi.vs.uni))), " drugs (with non-NA) have multi > uni\n"))

drug.families <- read.table("41586_2018_623_MOESM3_ESM-s11.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

png(paste0(plot.dir, "/multi-vs-uni-outliers.png"))
grid.table(subset(drug.families, inhibitor %in% multi.drugs)[, c("inhibitor", "family")], rows = NULL)
d <- dev.off()

png(paste0(plot.dir, "/multi-vs-uni-outliers.pdf"))
grid.table(subset(drug.families, inhibitor %in% multi.drugs)[, c("inhibitor", "family")], rows = NULL)
d <- dev.off()


if(FALSE) {
# , col = mean_drug_cor_spearman)) 
g <- ggplot()
g <- g + geom_point(data=plotData, aes(x =1:nrow(plotData) , y = (multi-uni)))
g <- g + xlab("")
g <- g + ylab("Performance Difference (Joint - Independent)")
g <- g + geom_text_repel(data=plotData[flag,],
			 aes(x =1:nrow(plotData), y = (multi-uni), label = label))
}

sz <- 12

set.text.size <- function(g, sz) {
g <- g +theme(text = element_text(size = sz),
             axis.text.x = element_text(size=sz),
             axis.text.y = element_text(size=sz),
             axis.title.x = element_text(size=sz),
             axis.title.y = element_text(size=sz),
             title = element_text(size=sz))
g
}

g1 <- ggplot()
g1 <- g1 + geom_point(data=plotData, aes(x = uni, y = multi))
g1 <- g1 + geom_text_repel(data=plotData[flag,],
			 aes(x = uni, y = multi, label = label))
g1 <- g1 + xlab("Spearman Correlation\n(Observed vs Predicted; Validation; Independent)")
g1 <- g1 + ylab("Spearman Correlation\n(Observed vs Predicted; Validation; Joint)")
g1 <- g1 + geom_abline(slope = 1, intercept = 0, linetype="dashed")
g1 <- set.text.size(g1, sz)
	     
# Does not appear to be a correlation between better multi - uni performance
# and mean_drug_cor (i.e., correlation between a given drug and all other drugs)
if(TRUE) {
g <- ggplot()
if(use.validation.results) {

  use.reticulate <- FALSE
  if(use.reticulate) {
      p_load(reticulate)
  } else {
      p_load(synapser)    
  }
  
  # `synasper` not working on local computer, so using `reticulate` to
  # log in to Synapse instead
  if(use.reticulate) {
      use_condaenv("synapse-r")
      synapseclient <- reticulate::import('synapseclient')
      syn <- synapseclient$Synapse()
      syn$login(silent = T)
  } else {
      synLogin()
  }


  read.csv.table.from.synapse <- function(synId, sep=",") {
      if(use.reticulate) {
          file <- syn$get(synId)$path
      } else {
          file <- synGet(synId)$path
      }
      ret <- suppressMessages(read.table(file, sep=sep, header=TRUE))
  }

  auc.synIds <-
      list("training" = "syn22397379",
           "leaderboard" = "syn22397867",
           "validation" = "syn22398053")

  aucs <-
      llply(auc.synIds,
            .fun = function(synId) read.csv.table.from.synapse(synId))

  auc.ranges <-
      llply(aucs,
            .fun = function(aucs.df) {
                apply(aucs.df, 2, mad, na.rm=TRUE)
            })


grd.cors <-
  llply(aucs,
        .fun = function(aucs.df) {
                 grd <- rowMeans(aucs.df, na.rm=TRUE)
                 drugs <- colnames(aucs.df)
                 names(drugs) <- drugs
                 method = "spearman"
                 method = "pearson"
                 ret.grd.cors <- ldply(drugs,
                                       .fun = function(drug) data.frame(grd.cor = cor(aucs.df[,drug], grd, method=method, use='pairwise.complete.obs')))
                 colnames(ret.grd.cors)[1] <- "inhibitor"
                 ret.grd.cors
               })

  plotData$inhibitor.name <- make.names(plotData$inhibitor)
  plotData <- merge(plotData, grd.cors[["validation"]], by.x = c("inhibitor.name"), by.y = c("inhibitor"))
  g <- g + geom_point(data=plotData, aes(x = multi - uni, y = grd.cor))
  g <- g + geom_text_repel(data=plotData[flag,],
			   aes(x = multi - uni, y = grd.cor, label = label))
  g <- g + ylab("Drug Response vs Mean Drug Response\n(Pearson Correlation)")
} else {
  g <- g + geom_point(data=plotData, aes(x = multi - uni, y = mean_drug_cor_spearman))
  g <- g + geom_text_repel(data=plotData[flag,],
			   aes(x = multi - uni, y = mean_drug_cor_spearman, label = label))
  g <- g + ylab("Mean Drug Spearman Correlation")
}
g <- g + xlab("Joint - Independent\n(Spearman Correlation; Observed vs Predicted)")
print(g)
ggsave(paste0(plot.dir, "/multi-minus-uni-vs-grd.png"))

print(g)
ggsave(paste0(plot.dir, "/multi-minus-uni-vs-grd.pdf"))
}

# I believe these are validation mads in inhibitor-stats.tsv; see calculate-replication-variance.R
# but instead, recalculate it here.
#file <- "inhibitor-stats.tsv"
#drug.stats <- read.table(file = file, sep="\t", header=TRUE)

if(!all(plotData$inhibitor.name %in% names(auc.ranges[["validation"]]))) {
  stop("Drug name mismatch\n")
}


val.drug.stats <- data.frame(inhibitor = names(auc.ranges[["validation"]]), validation.range = as.numeric(auc.ranges[["validation"]]))
train.drug.stats <- data.frame(inhibitor = names(auc.ranges[["training"]]), training.range = as.numeric(auc.ranges[["training"]]))
drug.stats <- merge(val.drug.stats, train.drug.stats)

df <- merge(plotData, drug.stats, by.x = c("inhibitor.name"), by.y = c("inhibitor"), all = TRUE)
# NB: exclude GRD, but unlike sc1-analysis.R, we do not exclude venetoclax
df <- df[!(df$inhibitor == "GRD"),]


df$multi.vs.uni <- df$multi - df$uni
df$multi.higher <- "FALSE"
flag <- !is.na(df$multi.vs.uni) & (df$multi.vs.uni > 0)
df[flag,"multi.higher"] <- "TRUE"
df$multi.higher <- factor(df$multi.higher, levels = c("TRUE", "FALSE"))

g <- ggplot()
g <- g + geom_violin(data=df, aes(x = multi.higher, y = validation.range))
g <- g + geom_boxplot(data=df, aes(x = multi.higher, y = validation.range), width = 0.25)
g <- g + geom_beeswarm(data=df, aes(x = multi.higher, y = validation.range))
g <- g + xlab("Joint Modeling > Independent Modeling Performance")
g <- g + ylab("Validation MAD")
print(g)
ggsave(paste0(plot.dir, "/multi-minus-uni-vs-valid-mad-box.png"))

print(g)
ggsave(paste0(plot.dir, "/multi-minus-uni-vs-valid-mad-box.pdf"))


g <- ggplot()
g <- g + geom_violin(data=df, aes(x = multi.higher, y = training.range))
g <- g + geom_boxplot(data=df, aes(x = multi.higher, y = validation.range), width = 0.25)
g <- g + geom_beeswarm(data=df, aes(x = multi.higher, y = training.range))
g <- g + xlab("Joint Modeling > Independent Modeling Performance")
g <- g + ylab("Training MAD")
print(g)
ggsave(paste0(plot.dir, "/multi-minus-uni-vs-train-mad-box.png"))

print(g)
ggsave(paste0(plot.dir, "/multi-minus-uni-vs-train-mad-box.pdf"))


col <- "range"
col <- "pearson"
col <- "cor"
# plot(df$multi - df$uni, df[,col])

flag <- df$label != ""
g2 <- ggplot(data = df, aes(y = multi - uni, x = validation.range))
g2 <- g2 + geom_point()
g2 <- g2 + ylab("Joint - Independent\n(Spearman Correlation; Observed vs Predicted)")
g2 <- g2 + xlab("Dynamic Range (MAD)")
g2 <- g2 + geom_text_repel(data=df[flag,],
			 aes(y = multi - uni, x = validation.range, label = label))
g2 <- g2 +theme(text = element_text(size = sz),
             axis.text.x = element_text(size=sz),
             axis.text.y = element_text(size=sz),
             axis.title.x = element_text(size=sz),
             axis.title.y = element_text(size=sz),
             title = element_text(size=sz))


grd.plts <- list()

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}




for(ds in c("training", "validation")) {

  all.tbl.with.grd <- merge(df[, !(colnames(df) == "grd.cor")], grd.cors[[ds]], by.x = c("inhibitor.name"), by.y = c("inhibitor"))

st <- boxplot.stats(all.tbl.with.grd$multi.vs.uni)
all.tbl.with.grd$label <- NA
flag <- !is.na(all.tbl.with.grd$multi.vs.uni) & (all.tbl.with.grd$multi.vs.uni > 0)
all.tbl.with.grd[flag, "label"] <- "Joint"
flag <- !is.na(all.tbl.with.grd$multi.vs.uni) & (all.tbl.with.grd$multi.vs.uni < 0)
all.tbl.with.grd[flag, "label"] <- "Uni"

all_out_index <- which(all.tbl.with.grd$multi.vs.uni %in% st$out)
pos_out_index <- which(all.tbl.with.grd$multi.vs.uni %in% st$out[st$out > 0])
neg_out_index <- which(all.tbl.with.grd$multi.vs.uni %in% st$out[st$out < 0])
all.tbl.with.grd[pos_out_index, "label"] <- "Joint++"
all.tbl.with.grd[neg_out_index, "label"] <- "Joint--"

#summary(lm(training.range ~ label, data=subset(all.tbl.with.grd, label %in% c("Uni","Joint++"))))
#summary(lm(validation.range ~ label, data=subset(all.tbl.with.grd, label %in% c("Uni","Joint++"))))

  g <- ggplot()
  g <- g + geom_violin(data=all.tbl.with.grd, aes_string(x = "label", y = "grd.cor"))
  g <- g + geom_boxplot(data=all.tbl.with.grd, aes_string(x = "label", y = "grd.cor"), width = 0.25)
  g <- g + geom_beeswarm(data=all.tbl.with.grd, aes_string(x = "label", y = "grd.cor"))
  g <- g + xlab("Joint Modeling > Independent Modeling Performance\n(Validation)")
  g <- set.text.size(g, sz)
print(g)
ggsave(paste0(plot.dir, "/joint-ind-label-vs-", ds, "-grd.png"), width = 14)

  g <- ggplot()
  g <- g + geom_violin(data=all.tbl.with.grd, aes_string(x = "label", y = paste0(ds, ".range")))
  g <- g + geom_boxplot(data=all.tbl.with.grd, aes_string(x = "label", y = paste0(ds, ".range")), width = 0.25)
  g <- g + geom_beeswarm(data=all.tbl.with.grd, aes_string(x = "label", y = paste0(ds, ".range")))
  g <- g + xlab("Joint Modeling > Independent Modeling Performance\n(Validation)")
  g <- set.text.size(g, sz)
print(g)
ggsave(paste0(plot.dir, "/joint-ind-label-vs-", ds, "-range.png"), width = 14)

  g <- ggplot()
  g <- g + geom_violin(data=all.tbl.with.grd, aes(x = multi.higher, y = grd.cor))
  g <- g + geom_boxplot(data=all.tbl.with.grd, aes(x = multi.higher, y = grd.cor), width = 0.25)
  g <- g + geom_beeswarm(data=all.tbl.with.grd, aes(x = multi.higher, y = grd.cor))
  g <- g + xlab("Joint Modeling > Independent Modeling Performance\n(Validation)")
  g <- set.text.size(g, sz)
  
  g <- g + ylab(paste0("Drug Response vs Mean Drug Response\n(Pearson Correlation; ", firstup(ds), ")"))
  grd.plts[[ds]] <- g

  flag <- all.tbl.with.grd$inhibitor %in% multi.drugs

  wt <- wilcox.test(grd.cor ~ multi.higher, data=all.tbl.with.grd, alternative="greater")
  met <- "mean"
  cat(paste0(ds, "; ", met, ": Comparing correlation with GRD between those with high vs lower performance in joint vs univariate: wilcox one-sided p-value: ", wt$p.value, "\n"))
  higher <- subset(all.tbl.with.grd, multi.higher=="TRUE")
  qs <- as.numeric(quantile(higher$grd.cor,probs=c(0.25,0.5,0.75)))
  cat(paste0(ds, "; ", met, ": Joint > Uni (n=", nrow(higher), ") correlation with GRD quantiles: 25% = ", qs[1], " 50% = ", qs[2], " 75% = ", qs[3], "\n"))
  lower <- subset(all.tbl.with.grd, multi.higher=="FALSE")
  qs <- as.numeric(quantile(lower$grd.cor,probs=c(0.25,0.5,0.75)))
  cat(paste0(ds, "; ", met, ": Joint > Uni correlation with GRD quantiles: 25% = ", qs[1], " 50% = ", qs[2], " 75% = ", qs[3], "\n"))
  print(wt)

  print(all.tbl.with.grd[flag, c("inhibitor.name", "grd.cor")])
  cat(paste0("mean grd of multi drugs = ", mean(all.tbl.with.grd[flag,"grd.cor"]), " vs rest = ", mean(all.tbl.with.grd[!flag,"grd.cor"]), "\n"))

  col <- paste0(ds, ".range")
  print(all.tbl.with.grd[flag, c("inhibitor.name", col)])
  cat(paste0("mean ", col, " of multi drugs = ", mean(all.tbl.with.grd[flag,col]), " vs rest = ", mean(all.tbl.with.grd[!flag,col]), "\n"))
}


library(cowplot)
source("plotting_utils.R")

g3 <- plot.correlation(df$validation.range, df$spearman.cor.uni, labels = df$label, display.pval = TRUE, display.r2 = FALSE, sz = sz)
g3 <- g3 + ylab("Spearman Correlation\n(Observed vs Predicted; Independent)")
g3 <- g3 + xlab("Dynamic Range (MAD)")
if(FALSE) {
g3 <- ggplot(data = df, aes(y = spearman, x = range))
g3 <- g3 + geom_point()
g3 <- g3 + geom_text_repel(data=df[flag,], aes(y = spearman, x = range, label = label))
g3 <- g3 +theme(text = element_text(size = sz),
             axis.text.x = element_text(size=sz),
             axis.text.y = element_text(size=sz),
             axis.title.x = element_text(size=sz),
             axis.title.y = element_text(size=sz),
             title = element_text(size=sz))
}

pg <- plot_grid(g1, g2, g3, nrow=1, labels="AUTO")
ggsave(paste0(plot.dir, "/joint-vs-ind-perf.png"), width = 14)

ggsave(paste0(plot.dir, "/joint-vs-ind-perf.pdf"), width = 14)


print(g1)
ggsave(paste0(plot.dir, "/joint-vs-ind-perf-only.png"), width = 14)

ggsave(paste0(plot.dir, "/joint-vs-ind-perf-only.pdf"), width = 14)


pg <- plot_grid(g1, grd.plts[["training"]], grd.plts[["validation"]], nrow=1, labels="AUTO")
ggsave(paste0(plot.dir, "/joint-vs-ind-perf-vs-grd.png"), width = 14)

ggsave(paste0(plot.dir, "/joint-vs-ind-perf-vs-grd.pdf"), width = 14)

pg <- plot_grid(g1, g2, nrow=1, labels="AUTO")
ggsave(paste0(plot.dir, "/joint-vs-ind-perf-no-mad.png"), width = 14)

ggsave(paste0(plot.dir, "/joint-vs-ind-perf-no-mad.pdf"), width = 14)
