import numpy as np
import pytest
import sympy as sp
from channel_first_jet import first_jet_factor, surrogate_kraus, choi


def test_complex_channel_surrogate():
    rng = np.random.default_rng(351)
    v, _ = np.linalg.qr(rng.normal(size=(4, 2))+1j*rng.normal(size=(4, 2)))
    h = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    dv = -1j*(h+h.conj().T)@v
    c = np.column_stack([k.T.ravel() for k in v.reshape(2, 2, 2)])
    dc = np.column_stack([k.T.ravel() for k in dv.reshape(2, 2, 2)])
    j, dj = c@c.conj().T, dc@c.conj().T+c@dc.conj().T
    a, da = first_jet_factor(j, dj)
    assert np.allclose(a@a.conj().T, j)
    assert np.allclose(da@a.conj().T+a@da.conj().T, dj)
    for theta in [-0.01, 0, 0.01]:
        ks = surrogate_kraus(a, da, 2, 2, theta)
        assert np.allclose(sum(k.conj().T@k for k in ks), np.eye(2))
        assert np.linalg.matrix_rank(choi(ks), tol=1e-10) == 2
    eps = 1e-6
    fd = (choi(surrogate_kraus(a, da, 2, 2, eps))-
          choi(surrogate_kraus(a, da, 2, 2, -eps)))/(2*eps)
    assert np.allclose(fd, dj, atol=1e-8)


def test_rank_birth_counterexample_has_same_first_jet():
    t = sp.Symbol('t', real=True)
    psi = sp.Matrix([sp.cos(t/2), sp.sin(t/2)])
    pure = psi*psi.T
    original = sp.kronecker_product(pure, sp.diag(1-t*t, t*t))
    surrogate = sp.kronecker_product(pure, sp.diag(1, 0))
    assert original.subs(t, 0) == surrogate.subs(t, 0)
    assert sp.diff(original, t).subs(t, 0) == sp.diff(surrogate, t).subs(t, 0)
    assert sp.diff(surrogate, t).subs(t, 0) != sp.zeros(4)


def test_reject_one_sided_rank_birth_tangent():
    with pytest.raises(ValueError, match='two-sided'):
        first_jet_factor(np.diag([1, 0]), np.diag([-1, 1]))
