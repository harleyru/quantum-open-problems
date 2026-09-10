import numpy as np
import pytest
from qcqc_qfi import ad_ensemble, gram_coefficients, hermitian_basis, solve


def test_full_complex_gauge_and_transpose():
    c, dc = ad_ensemble(2)
    rng = np.random.default_rng(1909)
    a = rng.normal(size=(16, 16))+1j*rng.normal(size=(16, 16))
    s = a@a.conj().T
    basis = hermitian_basis(4)
    assert len(basis) == 16
    assert np.allclose([[np.trace(a@b) for b in basis] for a in basis], np.eye(16))
    x = rng.normal(size=16)
    h = sum(v*b for v, b in zip(x, basis))
    d = dc-1j*c@h
    direct = 4*np.trace((d@d.conj().T).T@s).real
    g = (gram_coefficients(c, dc)@s.ravel()).real.reshape(17, 17)
    y = np.r_[1, x]
    assert np.allclose(y@g@y, direct)


@pytest.mark.parametrize('u', [0.1, 0.5, 0.9])
def test_single_use_analytic(u):
    result = solve(1, u)
    assert result['accepted'], result
    expected = 4*(1-u)/(1+np.sqrt(1-u))**2
    assert abs(result['objective']-expected) < 2e-5


def test_two_use_exceeds_parallel():
    result = solve(2)
    assert result['accepted'], result
    assert result['objective'] > 1.79457
    assert result['objective'] <= 4+2e-5
