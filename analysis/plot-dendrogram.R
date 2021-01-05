
make.plot.dendro <- function(mat, sig.vec, perf.vec, range.vec, cor.mat, structural.mat, drug.gene.cor.mat, column_dend, row_dend,
                             column_split, row_split, ...) {
    breaks = c(1, 10, 100, 1000)
    breaks <- breaks[breaks < max(sig.vec)]
##    if(max(sig.vec) > max(breaks)) { breaks <- c(breaks, max(sig.vec)) }
    log.breaks <- c(0, log10(breaks))
    label.breaks <- c(0, breaks)
    vals <- log10(as.numeric(sig.vec))
    vals[is.infinite(vals)] <- 0
    n.targets <- colSums(!is.na(mat))    
    col_ha = HeatmapAnnotation(perf = anno_barplot(perf.vec),
                               range = anno_barplot(range.vec),
                               sig = anno_barplot(vals,
                                                  axis_param = list(at = log.breaks, labels = label.breaks)),
                               targets = anno_barplot(n.targets))
                               
    names(col_ha) = c("Mean Pearson", "Drug Range (MAD)", "# Sig Genes", "# Targets")
    mat <- as.matrix(mat)
    cor.mat <- as.matrix(cor.mat)
    structural.mat <- as.matrix(structural.mat)    
    drug.gene.cor.mat <- as.matrix(drug.gene.cor.mat)
    row.flag <- rowSums(is.na(mat)) != ncol(mat)
    colors <- structure(1, names = c("Yes"))
    h <- Heatmap(mat[row.flag,], cluster_rows = FALSE, cluster_columns = column_dend, top_annotation = col_ha,
                 column_names_rot = 45, na_col = "white", show_heatmap_legend = FALSE, 
                 row_names_gp = gpar(fontsize = 4),
                 column_names_gp = gpar(fontsize = 4),
                 show_column_dend = TRUE, column_split = column_split,
                 height = 1, col = colors,
                 ...)
    h2 <- Heatmap(cor.mat, cluster_rows = row_dend, cluster_columns = column_dend, show_row_dend = FALSE,
                  name = "Drug/drug\ncorrelation",
                  column_split = column_split, row_split = row_split,
                  row_names_gp = gpar(fontsize = 4),
                  show_row_names = FALSE,
                  row_title = "Drug",
                  height = 1,
                  ...)
    h3 <- Heatmap(structural.mat, cluster_rows = row_dend, cluster_columns = column_dend, show_row_dend = FALSE,
                  name = "Drug/drug\nstructure similarity",
                  column_split = column_split, row_split = row_split,
                  row_names_gp = gpar(fontsize = 4),
                  show_row_names = FALSE,
                  row_title = "Drug",
                  height = 1,
                  ...)
    h4 <- Heatmap(drug.gene.cor.mat, cluster_rows = TRUE, cluster_columns = column_dend,
                  show_row_dend = FALSE,
                  name = "Drug/gene\ncorrelation",
                  column_split = column_split,
                  column_names_rot = 45,
                  column_names_gp = gpar(fontsize = 4),
                  row_names_gp = gpar(fontsize = 4),
                  show_row_names = FALSE,
                  row_title = "Gene",
                  column_title = "Drug",
                  column_title_side = "bottom",
                  height = 1,
                  ...)
    ht_list = h %v% h2 %v% h3 %v% h4
    ht_list = h %v% h2 %v% h3
    draw(ht_list, column_title = "Drug", column_title_side = "bottom")
}
