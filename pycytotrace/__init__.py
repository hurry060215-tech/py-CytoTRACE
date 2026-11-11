"""pycytotrace — pure-Python port of CytoTRACE (Gulati et al. *Science* 2020).

Predicts single-cell differentiation state from scRNA-seq counts using:
1. Sequencing-depth + Census normalisation
2. Gene Counts Signature (GCS) — mean of top 200 genes correlated with gene counts
3. NNLS regression of GCS onto cell-cell similarity matrix
4. Iterative diffusion smoothing
5. Rank-normalised → CytoTRACE score in [0,1] (1 = most stem-like, 0 = most differentiated)

v0.1 covers the core algorithm + ggplot2-python visualisation.
Batch correction (ComBat) and ``iCytoTRACE`` integration deferred to v0.2.

Citation: Gulati, G. S. et al. *Single-cell transcriptional diversity is a
hallmark of developmental potential.* Science 367, 405–411 (2020).
"""

from __future__ import annotations

__version__ = "0.1.0"

from .core import CytoTRACE, cytotrace_run
def __getattr__(name):
    if name in {"plot_cytotrace", "plot_cyto_genes"}:
        from . import plotting
        return getattr(plotting, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    "CytoTRACE",
    "cytotrace_run",
    "plot_cytotrace",
    "plot_cyto_genes",
    "__version__",
]
