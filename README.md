# py-CytoTRACE

A **Python port of [CytoTRACE](https://github.com/gunsagargulati/CytoTRACE)** (Gulati et al. *Science* 2020) — single-cell developmental potential / stemness prediction from gene counts.

- Pure NumPy / SciPy implementation, no Rcpp dependency
- Core algorithm: detected-gene count + most-variable-genes + GCS + NNLS + diffusion smoothing → rank-normalised score in [0, 1]
- ggplot2-python visualisation: ``plot_cytotrace`` (scatter on embedding) + ``plot_cyto_genes`` (top-correlated genes bar)

## Install

```bash
pip install pycytotrace-bio
```
(module name is `pycytotrace`.)

For the optional plotting example, install `pycytotrace-bio[plot]`. The numerical core does not require the plotting extra.

## Quick-start

```python
import pycytotrace
res = pycytotrace.cytotrace_run(counts)        # counts: (n_genes × n_cells)
res.cytotrace                                  # 1.0 = stem-like, 0.0 = differentiated
res.cyto_genes.head()                          # top genes correlated with score
from ggplot2_py import ggsave
ggsave("out.png",
       plot=pycytotrace.plot_cytotrace(res, embedding=umap_coords),
       width=6, height=4, dpi=120)
```

## Function map

| Python | R | Status |
|---|---|---|
| `cytotrace_run` | `CytoTRACE` | ✅ |
| `CytoTRACE` dataclass | `CytoTRACE()` return list | ✅ |
| `plot_cytotrace` | `plotCytoTRACE` | ✅ |
| `plot_cyto_genes` | `plotCytoGenes` | ✅ |
| `iCytoTRACE` (multi-batch integration) | `iCytoTRACE` | ⏳ v0.2 |

## Known limitations (v0.1)

1. **Batch correction (ComBat) skipped** — single-batch fit only. v0.2 will add a `batch=` argument that delegates to scanpy's `sc.pp.combat`.
2. **Fast-mode (subsampling for >3000 cells) not yet implemented** — large datasets compute the full N×N similarity matrix, which is O(N²) memory.
3. **iCytoTRACE deferred** to v0.2.

## Citation

> Gulati, G. S. et al. *Single-cell transcriptional diversity is a hallmark of developmental potential.* Science 367, 405–411 (2020).

## License

MIT.

## Updated real-data validation

See [measured R comparisons and limitations](RECONSTRUCTION_REPORT.md) and the executed real-data notebook in `examples/`. Synthetic signal tests are not evidence of R equivalence.
