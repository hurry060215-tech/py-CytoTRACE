args <- commandArgs(trailingOnly=TRUE)
upstream <- normalizePath(args[[1]], winslash='/')
out <- args[[2]]
source(file.path(upstream, 'R/CytoTRACE.R'))
load(file.path(upstream, 'data/marrow_10x_expr.rda'))
dir.create(out, recursive=TRUE, showWarnings=FALSE)
# Deterministic subset of real bone-marrow cells; retain the full gene axis.
mat <- as.matrix(marrow_10x_expr[, seq_len(120)])
write.csv(mat, file.path(out, 'counts.csv'))
set.seed(17)
res <- CytoTRACE(mat, enableFast=FALSE, ncores=1)
write.csv(data.frame(cell=colnames(mat), score=res$CytoTRACE, rank=res$CytoTRACErank,
                     counts=res$Counts, gcs=res$GCS), file.path(out, 'r_scores.csv'), row.names=FALSE)
write.csv(res$exprMatrix, file.path(out, 'r_expression.csv'))
writeLines(capture.output(sessionInfo()), file.path(out, 'R-sessionInfo.txt'))
# Deliberate boundary controls on the same real counts: an ERCC row and
# a cell with only five expressed genes. R restores that cell as NA.
edge <- rbind(mat, 'ERCC-control'=seq_len(ncol(mat)))
poor <- rep(0, nrow(edge)); poor[seq_len(5)] <- 1
edge <- cbind(edge, 'poor-control'=poor)
write.csv(edge, file.path(out, 'edge_counts.csv'))
edge_res <- CytoTRACE(edge, enableFast=FALSE, ncores=1)
write.csv(data.frame(cell=colnames(edge), score=edge_res$CytoTRACE,
                     rank=edge_res$CytoTRACErank, counts=edge_res$Counts,
                     gcs=edge_res$GCS), file.path(out, 'edge_r_scores.csv'), row.names=FALSE)
