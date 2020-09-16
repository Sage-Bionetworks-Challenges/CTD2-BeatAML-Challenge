suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(synapser))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(ggrepel))
suppressPackageStartupMessages(p_load(reshape2))
suppressPackageStartupMessages(p_load(ggcorrplot))

synLogin()

synId <- "syn22267537"
tbl <- synTableQuery(paste0("SELECT * FROM ", synId))
tbl <- as.data.frame(tbl)

## Exclude "gold" (which must mean gold-standard because its pearson and spearman = 1)
tbl <- subset(tbl, team != "gold")
tbl$pearson <- as.numeric(tbl$pearson)

drug.mean.tbl <-
    ddply(tbl, .variables = c("inhibitor"),
          .fun = function(df) data.frame(pearson = mean(df$pearson, na.rm=TRUE)))

tbl$inhibitor <- factor(tbl$inhibitor)

o <- order(drug.mean.tbl$pearson, decreasing = TRUE)
drug.mean.tbl <- drug.mean.tbl[o, ]
drug.mean.tbl$inhibitor <- factor(drug.mean.tbl$inhibitor, levels = drug.mean.tbl$inhibitor)
drug.mean.tbl$label <- ""
drug.mean.tbl[1:10, "label"] <- as.character(drug.mean.tbl[1:10, "inhibitor"])
drug.mean.tbl[(nrow(drug.mean.tbl)-9):nrow(drug.mean.tbl), "label"] <-
    as.character(drug.mean.tbl[(nrow(drug.mean.tbl)-9):nrow(drug.mean.tbl), "inhibitor"])

poor.performing.drugs <-
    as.character(drug.mean.tbl[(nrow(drug.mean.tbl)-9):nrow(drug.mean.tbl), "inhibitor"])

good.performing.drugs <-
    as.character(drug.mean.tbl[1:10, "inhibitor"])

## Taken from:
## https://stackoverflow.com/questions/4357031/qqnorm-and-qqline-in-ggplot2
gg_qq <- function(x, distribution = "norm", ..., line.estimate = NULL, conf = 0.95,
                  labels = names(x)){
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))
  x <- na.omit(x)
  ord <- order(x)
  n <- length(x)
  P <- ppoints(length(x))
  df <- data.frame(ord.x = x[ord], z = q.function(P, ...))

  if(is.null(line.estimate)){
    Q.x <- quantile(df$ord.x, c(0.25, 0.75))
    Q.z <- q.function(c(0.25, 0.75), ...)
    b <- diff(Q.x)/diff(Q.z)
    coef <- c(Q.x[1] - b * Q.z[1], b)
  } else {
    coef <- coef(line.estimate(ord.x ~ z))
  }

  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE

  if(!is.null(labels)){ 
    df$label <- ifelse(df$ord.x > df$upper | df$ord.x < df$lower, labels[ord],"")
    }

  p <- ggplot(df, aes(x=z, y=ord.x)) +
    geom_point() + 
    geom_abline(intercept = coef[1], slope = coef[2]) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2) 
  if(!is.null(labels)) p <- p + geom_text_repel( aes(label = label))
  p
}


g <- ggplot()
g <- g + geom_col(data = drug.mean.tbl, aes(x = inhibitor, y = pearson))
g <- g + coord_flip()
g <- g + scale_x_discrete(labels = drug.mean.tbl$label)
g <- g + theme(text = element_text(size = 5))

file <- "sc1-mean-drug-pearson-barplot.png"
file <- "sc1-mean-drug-pearson-barplot.pdf"
pdf(file)
print(g)
d <- dev.off()

## qqnorm(drug.mean.tbl$pearson)
## qqline(drug.mean.tbl$pearson)

g <- gg_qq(drug.mean.tbl$pearson, labels = as.character(drug.mean.tbl$inhibitor))
file <- "sc1-mean-drug-pearson-qqplot.png"
png(file)
print(g)
d <- dev.off()

poor.tbl <- subset(tbl, inhibitor %in% poor.performing.drugs)
good.tbl <- subset(tbl, inhibitor %in% good.performing.drugs)

g <- ggplot(data = poor.tbl)
g <- g + facet_wrap(~ inhibitor)
g <- g + geom_histogram(aes(x = pearson))
g <- g + ggtitle("Poor-performing drugs")
file <- "sc1-poor-performing-histograms.png"
png(file)
print(g)
d <- dev.off()


g <- ggplot(data = good.tbl)
g <- g + facet_wrap(~ inhibitor)
g <- g + geom_histogram(aes(x = pearson))
g <- g + ggtitle("Good-performing drugs")
file <- "sc1-good-performing-histograms.png"
png(file)
print(g)
d <- dev.off()

## Plot dynamic range vs response prediction performance (and label "good" and "poor" drugs)
synId <- "syn21212720"
file <- synGet(synId, downloadFile = TRUE)$path
aucs <- read.table(file, sep=",", header=TRUE, stringsAsFactors = FALSE)

auc.range <-
    ddply(aucs, .variables = c("inhibitor"),
          .fun = function(df) {
              data.frame(range = mad(df$auc))
          })

## Plot AUC dynamic range vs mean prediction performance
m <- merge(auc.range, drug.mean.tbl)

g <- ggplot(data = m, aes(x = range, y = pearson, label = label))
g <- g + geom_point()
g <- g + xlab("Dynamic Range (MAD)") + ylab("Mean Pearson Correlation")
g <- g + geom_text_repel()
g <- g + geom_smooth(method='lm', formula=y~x)
file <- "performance-vs-dynamic-range.png"
png(file)
print(g)
d <- dev.off()

auc.matrix <- acast(aucs, inhibitor ~ lab_id)
corr <- cor(t(auc.matrix[good.performing.drugs,]), use = "pairwise.complete.obs")
ggcorrplot(corr, title = "Good performing drugs")
file <- 

corr <- cor(t(auc.matrix[poor.performing.drugs,]), use = "pairwise.complete.obs")
ggcorrplot(corr, title = "Poor performing drugs")
