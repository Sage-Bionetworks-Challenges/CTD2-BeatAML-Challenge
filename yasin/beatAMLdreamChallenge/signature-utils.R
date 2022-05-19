# Code modified from:
# https://github.com/Sage-Bionetworks-Challenges/Celgene-Multiple-Myeloma-Challenge/blob/master/publishedClassifiers/common/classifier-base.R:run.classifier.eset

# Calculate a score defined by signature based on the expression matrix in expr.df.
# signature is a named vector of coefficients, with the names corresponding to 
# entries in the gene.col column of expr.df.
# The score is calculated over the sample.cols of the expr.df.
calculate.score <- function(signature, expr.df, gene.col, sample.cols) {
    common.genes <- intersect(names(signature), expr.df[, gene.col])
    missing.genes <- names(signature)[!(names(signature) %in% common.genes)]
    if(length(missing.genes) > 0) {
      warning(paste('The following genes are missing from the data: ',paste(missing.genes,collapse=', ')))
    }

    available <- match(names(signature),expr.df[, gene.col])
    if(all(is.na(available))) {
      warning("No genes from signature in data\n")
      return(NULL)
    }

    signature.df <- data.frame(gene = names(signature), coefficient = as.vector(signature))
    colnames(signature.df)[1] <- gene.col

    expr.df <- merge(expr.df, signature.df, all = FALSE)

    prod <- as.matrix(t(expr.df[, sample.cols])) %*% as.vector(expr.df[, "coefficient", drop = TRUE])
    raw.score <- prod[,1]
    names(raw.score) <- rownames(prod)
    return(raw.score)
}

# Code modified from: https://github.com/Sage-Bionetworks-Challenges/Celgene-Multiple-Myeloma-Challenge/blob/master/publishedClassifiers/common/classifier-base.R:adjust.scores.to.reflect.multimappers
## Adjust the classifier coefficients to account for a probe that maps
## to multiple genes.  If n genes map to a probe, take the average of the
## n gene expression values.  i.e., if the probe-based score is beta,
## make the gene-based score beta/n for each of the n genes.
## signature is a named vector of coefficients,
## where the names are assumed to match the gene.col column of 
## the expr.df (expression) matrix.
## Though the coefficients will be revised to account for redundantly
## named genes in the expression matrix (e.g., multiple gene symbols mapping to the
## same ensembl id), other columns of the matrix will not be accessed here.
## This returns a vector of coefficients named according to their respective gene.
adjust.scores.to.reflect.multimappers <- function(signature, expr.df, gene.col) {
  
  ## Drop any genes from the coefficients if they are not in our data.
  signature <- signature[names(signature) %in% expr.df[, gene.col]]

  ## Now adjust the score to account for multi-mappers
  
  ## Multiplicative factor to account for multi-mappers
  expr.df <- expr.df[expr.df[,gene.col] %in% names(signature),]
  probe.tbl <- as.data.frame(table(expr.df[, gene.col]))
  probe.tbl$Freq <- 1.0 / probe.tbl$Freq
  names(probe.tbl) <- c("gene", "factor")
  
  signature.df <- data.frame(gene = names(signature), coefficient = as.vector(signature))
  signature.df <- merge(signature.df, probe.tbl)
  
  signature.df$factor <- as.numeric(signature.df$factor)
  signature.df$coefficient <- as.numeric(signature.df$coefficient)    
  signature.df$adj.coefficient <- apply(signature.df[,c("factor", "coefficient")], 1, function(row) row[1] * row[2])
  adj.signature <- signature.df$adj.coefficient
  names(adj.signature) <- signature.df$gene
  return(adj.signature)
}