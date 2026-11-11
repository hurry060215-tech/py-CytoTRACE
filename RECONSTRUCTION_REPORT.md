# CytoTRACE reconstruction and measured validation

Validation date: 2026-11-10. This report supersedes the earlier synthetic-only parity claim. Correlation with detected-gene count is a sanity check, not an R comparison. The former test_exact_match.py is now test_synthetic_signal.py; actual R-output tests are in test_r_reference.py.

## Supported scope

The port implements the single-batch, non-subsampled core. ComBat, fast-mode subsampling, and iCytoTRACE remain unimplemented. The default accepts explicit genes-by-cells raw counts. input_scale="auto" reproduces R's maximum-below-50 inverse-log heuristic; input_scale="log2" accepts explicit log2 input. Explicit counts avoid misclassifying low-depth integer counts.

The repair retains the exact prevalence threshold, includes tied variable-gene dispersions at the cutoff, removes ERCC rows, filters cells with at most ten detected genes, and restores filtered cells as NaN in their original positions. Optional plotting no longer prevents importing the numerical core when ggplot2-python is absent.

## Direct R comparison

Reference: gunsagargulati/CytoTRACE@f4ff6c0ab93cc704eef44e17adcd5c6baf09703a. The unmodified CytoTRACE function ran with enableFast=FALSE, one core, and no batch labels. HiClimR and ccaPP were installed in an isolated validation library; reference functions were not replaced with Python.

| Input | Score maximum error | Rank difference | Count difference | GCS maximum error |
|---|---:|---:|---:|---:|
| Real marrow: 13,526 genes x 120 cells | 5.00e-16 | 0 | 0 | 7.78e-15 |
| Same counts + ERCC row + poor-quality cell | 5.00e-16 | 0 | 0 | 7.78e-15 |

Spearman agreement is 1 within numerical precision; the missing-score position also agrees. Controls are artificial, underlying marrow counts are real. No speedup or all-dataset equivalence is claimed.

## Reproduction

- examples/real_marrow_R_parity.ipynb: executed comparison and figures.
- tests/data/*.npz: compact inputs and actual R outputs used by regression tests.
- validation/provenance.json: hashes; metrics.json: errors; R-sessionInfo.txt: versions.
- validation/run_r_reference.R: original R execution with an upstream checkout.

```sh
Rscript validation/run_r_reference.R /path/to/CytoTRACE /path/to/output
python -m pytest -q
```

The fixture retains all genes in the first 120 columns of upstream marrow_10x_expr. The corrected official OmicVerse pancreas tutorial executes a documented 300-cell subset with all 27,998 genes and raw counts. Supplied PCA/UMAP coordinates are visualization only. Gene names are forwarded by the integration for ERCC filtering. This tutorial is not a full-cohort benchmark.
