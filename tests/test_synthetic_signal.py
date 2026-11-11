"""Synthetic signal sanity checks; direct R parity is tested separately."""
import sys
from pathlib import Path
import numpy as np
import pytest

_PORT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_PORT))
import pycytotrace


def test_recovers_gene_count_signal():
    """CytoTRACE's principal signal is detected-gene count; the algorithm's
    smoothed score should correlate strongly with raw gene counts. This is a synthetic sanity check, not an R/Python equivalence test."""
    rng = np.random.RandomState(42)
    n_cells, n_genes = 100, 200
    counts = rng.poisson(2.0, (n_genes, n_cells)).astype(float)
    for c in range(n_cells):
        extra = rng.choice(n_genes, size=max(1, n_cells - c), replace=False)
        counts[extra, c] += rng.poisson(5.0, len(extra))
    res = pycytotrace.cytotrace_run(counts)

    from scipy.stats import spearmanr
    rho_counts = spearmanr(res.counts, res.cytotrace)[0]
    assert rho_counts >= 0.80, (
        f"cytotrace vs detected-gene-count Spearman={rho_counts:.3f} below 0.80"
        f" — this is the algorithm's principal signal; failure suggests regression."
    )
