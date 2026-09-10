from fractions import Fraction
from pathlib import Path
import numpy as np
from certify_circuit_normalization import gram_bound, rational_grid, certify


def test_exact_integer_gram():
    re = np.array([[3], [0]], dtype=object)
    im = np.array([[0], [4]], dtype=object)
    assert gram_bound(re, im, 5) == 0
    assert gram_bound(re, im, 4) == Fraction(9, 16)
    re, im = rational_grid(np.array([[0.5+0.25j]]), 4)
    assert re[0, 0] == 2 and im[0, 0] == 1


def test_saved_cumulative_bound():
    path = Path(__file__).parent/'qfi_results/completed_circuit_n3.npz'
    result, _ = certify(path, 3)
    assert Fraction(result['output_trace_norm_bound']) < Fraction(1, 10**10)
    assert len(result['layer_bounds']) == 4
