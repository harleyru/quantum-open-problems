from fractions import Fraction
from pathlib import Path
import numpy as np
from mpmath import iv
from certify_qcqc_lower import certify, hypograph
from repair_qcqc_primal import repair, as_interval
from qcqc_qfi import ad_ensemble, gram_coefficients


def test_interval_gram_matches_direct_contraction():
    iv.dps = 50
    s = np.eye(4)/2
    interval = as_interval(np.array([[Fraction(str(x)) for x in row] for row in s]),
                           np.zeros((4, 4), dtype=int))
    g = hypograph(interval, 1, Fraction(0))
    actual = np.array([[float(z.real.mid) for z in row] for row in g])
    c, dc = ad_ensemble(1)
    expected = (gram_coefficients(c, dc)@s.ravel()).real.reshape(5, 5)
    scales = np.diag([1, 1, 1, np.sqrt(2), np.sqrt(2)])
    assert np.allclose(actual, scales@expected@scales)


def test_repaired_certificate_and_excessive_bound():
    source = Path(__file__).parent/'qfi_results/qfi_n3.npz'
    assert certify(source, 3, '4.7368')['verified']
    assert not certify(source, 3, '4.74')['verified']


def test_repair_keeps_fraction_entries():
    source = Path(__file__).parent/'qfi_results/qfi_n2.npz'
    process, parts = repair(source, 2)
    assert all(isinstance(x, Fraction) for a in process for x in a.ravel())
    assert all(isinstance(x, Fraction) for nodes in parts for a in nodes.values() for x in a.ravel())
