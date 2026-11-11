from pathlib import Path
import numpy as np
import pytest
from pycytotrace import cytotrace_run


@pytest.mark.parametrize('prefix', ['', 'edge_'])
def test_real_marrow_matches_executed_r_reference(prefix):
    fixture = np.load(Path(__file__).parent/'data'/f'{prefix}marrow_reference.npz')
    result = cytotrace_run(fixture['counts'], gene_names=fixture['genes'].tolist(),
                          cell_names=fixture['cells'].tolist())
    for field, key in [('cytotrace','scores'),('cytotrace_rank','ranks'),('counts','detected'),('gcs','gcs')]:
        np.testing.assert_allclose(getattr(result,field),fixture[key],rtol=1e-9,atol=1e-10,equal_nan=True)


def test_variable_gene_prevalence_is_not_rounded_down():
    from pycytotrace.core import _most_variable_genes
    values = np.zeros((2, 21))
    values[0, 0] = 10
    values[1, :2] = [1, 2]
    np.testing.assert_array_equal(_most_variable_genes(values), [1])


def test_variable_gene_dispersion_ties_are_retained():
    from pycytotrace.core import _most_variable_genes
    values = np.tile([1., 2., 3.], (5, 1))
    np.testing.assert_array_equal(_most_variable_genes(values, n_top=2), np.arange(5))


def test_r_log2_input_heuristic_matches_real_count_reference():
    fixture = np.load(Path(__file__).parent/'data/marrow_reference.npz')
    result = cytotrace_run(np.log2(fixture['counts'] + 1),
                          gene_names=fixture['genes'].tolist(), input_scale='auto')
    np.testing.assert_allclose(result.cytotrace, fixture['scores'], atol=1e-10, rtol=1e-9)
