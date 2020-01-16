##################
## DE functions
##################
## design_mat should have the intercept term in all of these

get_voom <- function(countdat, design_mat) {
  vout <- limma::voom(counts = countdat, design = design_mat)
  lout <- limma::lmFit(vout)
  eout <- limma::eBayes(lout)
  retlist <- list(bhat = coefficients(eout)[, 2],
                  pval = eout$p.value[, 2])
  return(retlist)
}

get_DESeq2 <- function(countdat, design_mat) {
  colnames(design_mat) <- c("I", "D")
  design_mat <- as.data.frame(design_mat)
  design_mat$D <- factor(design_mat$D)
  trash <- capture.output({
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdat, colData = design_mat, design = ~D)
    dds <- DESeq2::DESeq(object = dds, quiet = TRUE)
    res <- DESeq2::results(dds)
  })
  retlist <- list(bhat = res$log2FoldChange, pval = res$pvalue)
  return(retlist)
}


get_edgeR <- function(countdat, design_mat) {
  trash <- capture.output({
    dge <- edgeR::DGEList(counts = countdat, group = design_mat[, 2])
    dge <- edgeR::estimateDisp(y = dge)
    efit <- edgeR::exactTest(dge)
  })
  retlist <- list(bhat = efit$table$logFC,
                  pval = efit$table$PValue)
}
