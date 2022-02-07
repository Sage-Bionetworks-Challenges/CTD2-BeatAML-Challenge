library(pacman)
p_load(ggplot2)
p_load(ggrepel)

file <- "multi_uni_plot_data.csv"
plotData <- read.table(file, sep=",", header=TRUE)

mean_diff <- mean(plotData$multi - plotData$uni, na.rm = T)

plotData$label <- plotData$drug
flag <- (plotData$multi - plotData$uni) > 0.25
plotData[!flag,"label"] <- ""

if(FALSE) {
# , col = mean_drug_cor_spearman)) 
g <- ggplot()
g <- g + geom_point(data=plotData, aes(x =1:nrow(plotData) , y = (multi-uni)))
g <- g + xlab("")
g <- g + ylab("Performance Difference (Joint - Independent)")
g <- g + geom_text_repel(data=plotData[flag,],
			 aes(x =1:nrow(plotData), y = (multi-uni), label = label))
}

sz <- 15

g1 <- ggplot()
g1 <- g1 + geom_point(data=plotData, aes(x = uni, y = multi))
g1 <- g1 + xlab("Spearman Correlation\n(Observed vs Predicted; Independent)")
g1 <- g1 + ylab("Spearman Correlation\n(Observed vs Predicted; Joint)")
g1 <- g1 + geom_text_repel(data=plotData[flag,],
			 aes(x = uni, y = multi, label = label))
g1 <- g1 + geom_abline(slope = 1, intercept = 0, linetype="dashed")
g1 <- g1 +theme(text = element_text(size = sz),
             axis.text.x = element_text(size=sz),
             axis.text.y = element_text(size=sz),
             axis.title.x = element_text(size=sz),
             axis.title.y = element_text(size=sz),
             title = element_text(size=sz))
	     
# Does not appear to be a correlation between better multi - uni performance
# and mean_drug_cor (i.e., correlation between a given drug and all other drugs)
if(TRUE) {
g <- ggplot()
g <- g + geom_point(data=plotData, aes(x = multi - uni, y = mean_drug_cor_spearman))
g <- g + geom_text_repel(data=plotData[flag,],
			 aes(x = multi - uni, y = mean_drug_cor_spearman, label = label))
g <- g + xlab("Joint - Independent\n(Spearman Correlation; Observed vs Predicted)")
g <- g + ylab("Mean Drug Spearman Correlation")
print(g)
ggsave("multi-minus-uni-vs-grd.png")
}

file <- "inhibitor-stats.tsv"
drug.stats <- read.table(file = file, sep="\t", header=TRUE)

df <- merge(plotData, drug.stats, by.x = c("drug"), by.y = c("inhibitor"))
col <- "range"
col <- "pearson"
col <- "cor"
# plot(df$multi - df$uni, df[,col])

flag <- df$label != ""
g2 <- ggplot(data = df, aes(y = multi - uni, x = range))
g2 <- g2 + geom_point()
g2 <- g2 + ylab("Joint - Independent\n(Spearman Correlation; Observed vs Predicted)")
g2 <- g2 + xlab("Dynamic Range (MAD)")
g2 <- g2 + geom_text_repel(data=df[flag,],
			 aes(y = multi - uni, x = range, label = label))
g2 <- g2 +theme(text = element_text(size = sz),
             axis.text.x = element_text(size=sz),
             axis.text.y = element_text(size=sz),
             axis.title.x = element_text(size=sz),
             axis.title.y = element_text(size=sz),
             title = element_text(size=sz))

library(cowplot)
source("plotting_utils.R")

g3 <- plot.correlation(df$range, df$spearman, labels = df$label, display.pval = TRUE, display.r2 = FALSE, sz = sz)
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
ggsave("joint-vs-ind-perf.png", width = 14)

pg <- plot_grid(g1, g2, nrow=1, labels="AUTO")
ggsave("joint-vs-ind-perf-no-mad.png", width = 14)
