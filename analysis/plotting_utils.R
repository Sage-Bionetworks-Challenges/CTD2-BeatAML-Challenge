
my_format <- function(x) ifelse(x < 0.001, formatC(x, format="e",digits=2), signif(x,digits=2))

lm_corr_eqn <- function(df, method = "pearson", display.r2 = FALSE, display.pval = FALSE){
    m <- lm(y ~ x, df);
    ct <- cor.test(df$x, df$y, method = method)
    estimate <- as.numeric(ct$estimate)
    if(display.r2 == TRUE) { estimate <- estimate*estimate }
    pval <- ct$p.value
    cat(paste0("method = ", method, " estimate = ", estimate, " pval = ", pval, "\n"))
    eq <- NULL
    digits <- 2
    if((method == "pearson") && (display.r2 == TRUE)) { 
      if(display.pval) { 
        eq <- substitute(italic(r)^2~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=digits, scientific=0),
                              pval = my_format(pval)))
      } else {
        eq <- substitute(italic(r)^2~"="~est, 
                         list(est = format(estimate, digits=digits, scientific=0)))

      }
    } else if((method == "pearson") && (display.r2 == FALSE)) {
      if(display.pval) { 
        eq <- substitute(italic(r)~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=digits, scientific=0),
                              pval = my_format(pval)))			 
      } else {
        eq <- substitute(italic(r)~"="~est, 
                         list(est = format(estimate, digits=digits, scientific=0)))

      }
    } else if((method == "spearman") && (display.r2 == FALSE)) {
      if(display.pval) { 
        eq <- substitute(rho~"="~est*","~~italic(p)~"="~pval, 
                         list(est = format(estimate, digits=digits, scientific=0),
                              pval = my_format(pval)))			 			 
      } else {
        eq <- substitute(rho~"="~est, 
                         list(est = format(estimate, digits=digits, scientific=0)))

      }
    } else {
      stop(paste("lm_corr_eqn does not know how to handle method = ", method,  " display.r2 = ", display.r2, "\n"))
    }
    as.character(as.expression(eq));                 
}


geom.text.to.theme.text.size.ratio <- 5/14

plot.correlation <- function(x, y, labels = NULL, colors = NULL, display.r2 = FALSE, method = "pearson", display.pval = FALSE, xoffset = 0.5, yoffset = 0.9, sz = 25, geom.text.size = 5, ...) {
  df <- data.frame(x = x, y = y)
  if(!is.null(labels)) {
    df$labels <- labels
  }
  g <- NULL
  if(is.null(labels)) {
    g <- ggplot(df, aes(x = x, y = y))
  } else {
    g <- ggplot(df, aes(x = x, y = y, label = labels))
  }
  if(!is.null(colors)) {
    g <- g + geom_point(aes(colour = colors))
  } else {
    g <- g + geom_point()
  }
  if(!is.null(labels)) {
##    g <- g + geom_text(vjust = "inward", hjust = "inward")
    suppressPackageStartupMessages(p_load(ggrepel))
    g <- g + geom_text_repel(point.padding = NA, box.padding = 1)
  }
##  g <- g + theme(legend.position="none")
  g <- g + geom_smooth(data = df, aes(x = x, y = y), method='lm')
  x.min <- min(df$x, na.rm=TRUE)
  x.max <- max(df$x, na.rm=TRUE)
  y.min <- min(df$y, na.rm=TRUE)
  y.max <- max(df$y, na.rm=TRUE)

  ylimits <- NULL
if(FALSE) {
  use.ggplot.2.2.1.limit.code <- TRUE
  if(use.ggplot.2.2.1.limit.code) {
    ylimits <- ggplot_build(g)$layout$panel_ranges[[1]]$y.range
    xlimits <- ggplot_build(g)$layout$panel_ranges[[1]]$x.range
  } else {
    ylimits <- ggplot_build(g)$layout$panel_params[[1]]$y.range
    xlimits <- ggplot_build(g)$layout$panel_params[[1]]$x.range
  }
}
  xlimits <- ggplot_build(g)$layout$panel_params[[1]]$x.range
  ylimits <- ggplot_build(g)$layout$panel_params[[1]]$y.range

## to see why geom_text(size = sz) sz is different than in theme see: ratio of 14/5
## https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control/25062509 
##  g <- g + geom_text(x = x.min + 0.5 * (x.max - x.min), y = y.min + 1 * (y.max - y.min), label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
##  g <- g + geom_text(x = x.min + 0.5 * (x.max - x.min), y = y.min + 0.8 * (y.max - y.min), label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
##  g <- g + geom_text(x = x.min + 0.5 * (x.max - x.min), y = 0.8 * ylimits[2], label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
  g <- g + geom_text(size = geom.text.size, x = xlimits[1] + xoffset * (xlimits[2] - xlimits[1]), y = ylimits[1] + yoffset * (ylimits[2] - ylimits[1]), label = lm_corr_eqn(df, method = method, display.r2 = display.r2, display.pval = display.pval), parse=TRUE, ...)
  g <- g +theme(text = element_text(size = sz),
             axis.text.x = element_text(size=sz),
             axis.text.y = element_text(size=sz),
             axis.title.x = element_text(size=sz),
             axis.title.y = element_text(size=sz),
             title = element_text(size=sz),
             plot.title = element_text(hjust = 0.5, size=sz))
  g
}
