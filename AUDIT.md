# Scoped R function mapping

| R function | Python entry point | Validation scope |
|---|---|---|
| CytoTRACE | cytotrace_run | Direct real-marrow comparison, no batch correction or fast subsampling |
| plotCytoTRACE | plot_cytotrace | Python plotting API; not pixel equivalence |
| plotCytoGenes | plot_cyto_genes | Python plotting API; not pixel equivalence |
| iCytoTRACE | Unimplemented | No equivalence claim |

Name matching alone missed the core entry point in the historical audit. This mapping describes supported APIs; it is not a percentage of algorithmic equivalence. See RECONSTRUCTION_REPORT.md for measured R errors and remaining limits.
