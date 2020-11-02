suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(ggplot2))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))

suppressPackageStartupMessages(p_load("foreach"))
suppressPackageStartupMessages(p_load("parallel"))

suppressPackageStartupMessages(p_load("reshape2"))

suppressPackageStartupMessages(p_load("synapser"))
synLogin()

num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}
num.processes <- num.cores - 1

## See https://stackoverflow.com/questions/27803710/ggplot2-divide-legend-into-two-columns-each-with-its-own-title
plot.anno.heatmap.with.multiple.legends <-
    function(df, id.col, heatmap.columns, legend.columns, legend.pals) {

        suppressPackageStartupMessages(p_load("RColorBrewer"))
        all.columns <- unique(c(legend.columns, heatmap.columns))
        df <- df[, c(id.col, all.columns)]

        ## Assume annotations are characters
        for(col in c(all.columns)) {
            df[, col] <- as.character(df[, col])
        }

        columns <- 1:length(legend.columns)
        names(columns) <- legend.columns

        color.vecs <-
            llply(columns,
                  .fun = function(idx) {
                      anno.col <- legend.columns[idx]
                      vec <- unique(df[, anno.col])
                      len <- length(vec)
                      colors <- brewer.pal(len, legend.pals[idx])
                      names(colors) <- vec
                      colors
                  })

        all.colors <- Reduce("c", color.vecs)
        names(all.colors) <- Reduce("c", unlist(lapply(color.vecs, names)))

        names(heatmap.columns) <- heatmap.columns
        names(all.columns) <- all.columns
        anno.df <- ldply(all.columns,
                     .fun = function(anno.col) {
                         data.frame(val = df[, anno.col], id = df[, id.col])
                     })
        colnames(anno.df)[1] <- "type"

        tmp <- subset(anno.df, type %in% heatmap.columns)
        tmp$type <- factor(tmp$type, levels = heatmap.columns)
        print(tmp$id)
        print(levels(tmp$id))
        full.plot <-
            ggplot(tmp, aes(y = id, x = type, fill = val)) + geom_tile() +
            scale_fill_manual(values = all.colors) +
            theme(legend.position="none")

        full.plot <- full.plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                                       axis.ticks.y = element_blank(), text = element_text(size = 18),
                                       axis.text.x = element_text(angle = 45, hjust = 1),
                                       axis.title.x = element_blank())
        

        names(legend.columns) <- legend.columns
        legends <-
            llply(legend.columns,
                  .fun = function(anno.col) {
                      flag <- anno.df$type == anno.col
                      g <- ggplot(anno.df[flag, ], aes_string(x = "id", y = "type", fill = "val"))
                      g <- g + geom_tile()
                      g <- g + scale_fill_manual(values = all.colors, name = anno.col)
                  })

        return(list("full.plot" = full.plot, "legends" = legends))
    }

## Modified slightly from
## https://stackoverflow.com/questions/53170465/how-to-make-a-base-r-style-boxplot-using-ggplot2
geom_boxplotMod <- function(mapping = NULL, data = NULL, stat = "boxplot", 
    position = "dodge2", ..., outlier.colour = NULL, outlier.color = NULL, 
    outlier.fill = NULL, outlier.shape = 1, outlier.size = 1.5, 
    outlier.stroke = 0.5, outlier.alpha = NULL, notch = FALSE, notchwidth = 0.5,
    varwidth = FALSE, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
    linetype = "dashed") # to know how these come here use: args(geom_boxplot)
    {
    list(geom_boxplot(
            mapping = mapping, data = data, stat = stat, position = position,
            outlier.colour = outlier.colour, outlier.color = outlier.color, 
            outlier.fill = outlier.fill, outlier.shape = outlier.shape, 
            outlier.size = outlier.size, outlier.stroke = outlier.stroke, 
            outlier.alpha = outlier.alpha, notch = notch, 
            notchwidth = notchwidth, varwidth = varwidth, na.rm = na.rm, 
            show.legend = show.legend, inherit.aes = inherit.aes, linetype = 
            linetype, ...),
        stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.25),
        #the width of the error-bar heads are decreased
        stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.25),
        stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), ...)
        )
    }


## skhurana just submitted the baseline
teams.to.drop <- c("skhurana", "gold")

## Read in the annotations; keep only the columns of interest; encode them as characters (not levels); and translate the names
suppressPackageStartupMessages(p_load(xlsx))
anno <- read.xlsx("Beat-AML-participant-summary.xlsx", sheetIndex = 1)
anno.df <- anno[, c("team", "SC1.method.class", "SC1.expr.covariates", "SC1.genomic.covariates", "SC1.clinical.covariates")]
for(col in colnames(anno.df)) { anno.df[, col] <- as.character(anno.df[, col]) }
anno.cols <- c("Method", "Expr", "Mut", "Clin")
colnames(anno.df) <- c("team", anno.cols)
## Set NAs to unknown
anno.df[is.na(anno.df)] <- "Unk"
rename.list <- list("yes" = "Yes", "no" = "No", "independent" = "Ind", "joint" = "Joint")
for(col in anno.cols) {
    for(nm in names(rename.list)) {
        flag <- anno.df[, col] == nm
        anno.df[flag,col] <- rename.list[[nm]]
    }
}
## Add a row for baseline (ridge) and MKL
anno.df <- rbind(anno.df, data.frame("team" = "ridge", "Method" = "Ind", "Expr" = "Yes", "Mut" = "Yes", "Clin" = "Yes"))
anno.df <- rbind(anno.df, data.frame("team" = "MKL", "Method" = "Joint", "Expr" = "Yes", "Mut" = "Yes", "Clin" = "Yes"))
flag <- anno.df$team == "Bioinformatics_Class_Challenge"
anno.df[flag, "team"] <- "BCC"
anno.df <- subset(anno.df, !(team %in% teams.to.drop))

## Begin generate dummy data
teams <- unique(anno.df$team)
names(teams) <- teams
n.boot <- 1000
use.dummy <- FALSE

if(use.dummy) {
    scores <- ldply(teams,
                    .fun = function(team) {
                        rand.mean <- runif(1, min = 0.5, max = 0.75)
                        data.frame(boot = 1:n.boot, pearson = rnorm(n.boot, rand.mean, sd = 0.05),
                                   spearman = rnorm(n.boot, rand.mean, sd = 0.05))
                    })
    colnames(scores)[1] <- "team"
    mean.pearson.scores <- ddply(scores, .variables = c("team"), .fun = function(df) data.frame(pearson = mean(df$pearson)))
} else {
    ## Read in the raw prediction results
    synId <- "syn22156457"
    pred.sc1 <- synTableQuery(paste0("SELECT * FROM ", synId))
    pred.sc1 <- as.data.frame(pred.sc1)
    
    pred.mat.sc1 <- unique(pred.sc1[, c("lab_id", "inhibitor", "auc_gold")])
    colnames(pred.mat.sc1)[3] <- "gold"
    for(team in unique(pred.sc1$team_name)) {
        df <- subset(pred.sc1, team_name == team)
        df <- df[, c("lab_id", "inhibitor", "auc_pred")]
        ## Rename team with long name
        if(team == "Bioinformatics_Class_Challenge") { team <- "BCC" }
        colnames(df)[3] <- team
        pred.mat.sc1 <- merge(pred.mat.sc1, df, all.x = TRUE)
    }

    ## Add in the baseline/ridge results
    ridge.res <- read.table("baseline_scores_by_drug.csv", sep=",", header=TRUE)
    ridge.res <- ridge.res[, c("lab_id", "inhibitor", "baseline")]
    colnames(ridge.res) <- c("lab_id", "inhibitor", "ridge")
    pred.mat.sc1 <- merge(pred.mat.sc1, ridge.res, all.x = TRUE)

    ## Add in Anna's MKL results
    mkl.res <- read.table("predictions__baseline_model__K_expr_mut__17102020.csv", sep=",", header=TRUE)
    colnames(mkl.res) <- c("lab_id", "inhibitor", "MKL")
    pred.mat.sc1 <- merge(pred.mat.sc1, mkl.res, all.x = TRUE)    
    
    pred.mat.sc1 <- subset(pred.mat.sc1, !is.na(gold))
    
    ## FUNCTION:
    ## Calculate average correlation per inhibitor of each column.
    sc1_metric <- function(df, cor_method) {
        cors <- plyr::daply(df, .variables = "inhibitor", .fun = function(x) {
            
                                        # Ignore the first column (inhibitor)
            apply(x[,-1], 2, function(y) {
                cor(y, x[,2], method = cor_method)
            })
        })
        cors[is.na(cors)] <- 0
        apply(cors, 2, mean)
    }
    
    set.seed(10)
    
    N <- 1000
    bs_indices.sc1 <- matrix(1:nrow(pred.mat.sc1), nrow(pred.mat.sc1), N) %>%
        apply(2, sample, replace = T)
    
    boot.sc1 <- apply(bs_indices.sc1, 2, function(ind) {
        sc1_metric(pred.mat.sc1[ind, -1], "spearman")
    }) %>% t()

    boot.sc1 <- adply(bs_indices.sc1, .margins = 2, .parallel = TRUE, .fun = function(ind) {
        tmp <- t(as.data.frame(sc1_metric(pred.mat.sc1[ind, -1], "spearman")))
        rownames(tmp) <- NULL
        tmp
    }) 
    boot.sc1 <- boot.sc1[, -1]
    scores <- melt(boot.sc1)
    colnames(scores) <- c("team", "spearman")
    scores$team <- as.character(scores$team)
    ## scores <- subset(scores, team != "gold")
    scores <- subset(scores, !(team %in% teams.to.drop))
    
    tmp <- sc1_metric(pred.mat.sc1[,-1], "pearson")
    mean.pearson.scores <- data.frame(team = names(tmp), pearson = as.numeric(tmp))
    ## mean.pearson.scores <- subset(mean.pearson.scores, team != "gold")
    mean.pearson.scores <- subset(mean.pearson.scores, !(team %in% teams.to.drop))
}

## End generate dummy data

## BEGIN EDIT HERE
## Can remove code between "Begin generate dummy data" and "End generate dummy data"
## Should add here 2 data.frames:
## (1) scores with columns "team" (team name), "boot" (1:number of bootstraps), "spearman" (pearson correlation)
## (2) mean.pearman.scores with columns "team" and "pearson"
## Note that mean.spearman.scores will be defined from scores
## END EDIT HERE

mean.spearman.scores <- ddply(scores, .variables = c("team"), .fun = function(df) data.frame(spearman = mean(df$spearman)))

## Order teams by the mean (including the annotations!)
o <- order(mean.spearman.scores$spearman, decreasing = FALSE)
lvls <- mean.spearman.scores[o, "team"]
scores$team <- factor(scores$team, levels = lvls)
mean.spearman.scores$team <- factor(mean.spearman.scores$team, levels = lvls)
mean.pearson.scores$team <- factor(mean.pearson.scores$team, levels = lvls)

anno.df$team <- factor(anno.df$team, levels = lvls)

## Plot boxplot
g1 <- ggplot(data = scores, aes_string(x = "team", y = "spearman"))
g1 <- g1 + geom_boxplotMod(fill = "#56B4E9")
g1 <- g1 + coord_flip()
g1 <- g1 + xlab("Method")
## g1 <- g1 + ylab("Pearson Correlation")
g1 <- g1 + ylab("Pearson")
g1 <- g1 + theme(text = element_text(size=18), title = element_text(size = 20),
                 axis.text.x = element_text(angle = 45, hjust = 1))

## Plot barplot
g2 <- ggplot(data = mean.pearson.scores)
g2 <- g2 + geom_col(aes_string(x = "team", y = "pearson"), fill = "#E69F00")
g2 <- g2 + coord_flip()
g2 <- g2 + xlab("Method")
g2 <- g2 + ylab("Spearman")
g2 <- g2 + theme(text = element_text(size=18))    
g2 <- g2 + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1))

suppressPackageStartupMessages(p_load(cowplot)) ## for plot_grid

## Create a dummy column that has all of the entries for Mut/Expr/Clin (i.e., yes, no, unk)
anno.df$Dummy <- rep(c("Yes", "No", "Unk"), len = nrow(anno.df))

ret <- plot.anno.heatmap.with.multiple.legends(anno.df, "team", heatmap.columns = c("Method", "Mut", "Expr", "Clin"),
                                               legend.columns = c("Method", "Dummy"),
                                               legend.pals = c("Paired", "Set2"))
    
full.plot <- ret[["full.plot"]]
for.first.legend <- ret[["legends"]][["Method"]]
for.first.legend <- for.first.legend + theme(text = element_text(size = 18))
for.second.legend <- ret[["legends"]][["Dummy"]]
for.second.legend <- for.second.legend + guides(fill=guide_legend(title = "Mut/Expr/Clin"))
for.second.legend <- for.second.legend + theme(text = element_text(size = 18))

legs <- plot_grid(get_legend(for.first.legend), get_legend(for.second.legend), nrow = 2, align = "v", rel_heights = c(1,1))
pg <- plot_grid(g1, g2, full.plot, legs, nrow=1, align="h", rel_widths = c(3,0.6,0.4,0.5))

png("sc1-scores.png", width = 2 * 480)
print(pg)
d <- dev.off()
