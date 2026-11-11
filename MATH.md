# CytoTRACE core: mathematical scope and numerical evidence

The input matrix is genes by cells. The supported core has no batch correction
or fast subsampling. See RECONSTRUCTION_REPORT.md for the pinned R reference.

## Processing order

1. Remove zero-variance genes. Interpret input according to explicit counts,
   log2, or the optional R-compatible maximum-below-50 heuristic.
2. Remove ERCC rows and cells with at most ten detected genes.
3. Scale each cell to one million total counts and apply log2(x + 1).
4. Census normalization rescales inverse-log values to each cell's number of
   detected genes, then reapplies log2(x + 1).
5. Keep genes expressed in at least 5% of cells; select dispersion var/mean
   at or above the top-1000 cutoff, retaining ties. Remove cells expressing
   none of those genes.
6. Form a cell-cell Pearson correlation matrix from variable genes. Compute
   its mean before zeroing the diagonal; zero negative and below-cutoff entries,
   then row-normalize, keeping zero rows zero.
7. GCS is the mean normalized expression of the 200 genes most correlated with
   detected-gene counts. Fit nonnegative coefficients b minimizing ||S b-GCS||².
8. Diffuse the fitted score with v(t+1)=0.9 S v(t)+0.1 v(0), stopping when the
   mean absolute update is below 1e-6 or after 10,000 iterations.
9. Rank and rescale scores to [0,1]; restore filtered cells as missing values in
   original order. The Python helper returns zeros for a constant score range;
   this is not an assertion of R behavior in every degenerate input case.

## What is and is not established

R uses HiClimR/ccaPP for correlations and nnls for the nonnegative fit; Python
uses NumPy/SciPy. Pearson correlation is not a robust-correlation estimator.
Sharing an optimization objective or an algorithm family does not establish
bitwise equivalence between implementations or BLAS libraries.

The executed real-marrow comparisons give maximum score error 5.00e-16 and
identical ranks on the tested 120 cells, including the separately described
ERCC/poor-cell controls. These measurements replace the historical claim that
R could not be run. They do not establish bit equivalence or complete agreement
on every dataset, degeneracy, preprocessing mode or platform. Tests, fixtures,
session versions and measured errors are stored in tests/ and validation/.
