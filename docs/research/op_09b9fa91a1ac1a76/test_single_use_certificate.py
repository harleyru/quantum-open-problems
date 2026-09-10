from fractions import Fraction
import numpy as np
import pytest
from single_use_certificate import half_noise_interval, symbolic_certificate


def test_exact_symbolic_certificate():
    assert all(value == 0 for value in symbolic_certificate().values())


def test_integer_only_enclosure():
    lo, hi = half_noise_interval()
    assert hi-lo == Fraction(8, 10**30)
    # F=12-8 sqrt(2): verify both endpoints by exact rational squaring.
    assert ((12-hi)/8)**2 < 2 < ((12-lo)/8)**2


@pytest.mark.parametrize('u', [0.1, 0.5, 0.9])
def test_explicit_probe_and_measurement(u):
    s = np.sqrt(1-u)
    p = s/(1+s)
    # Prepare R_y(2 arccos sqrt(p)) on target, then CNOT target -> memory.
    ry = np.array([[np.sqrt(p), -np.sqrt(1-p)],
                   [np.sqrt(1-p), np.sqrt(p)]])
    cnot = np.eye(4)[[0, 1, 3, 2]]
    state = cnot@np.kron(ry, np.eye(2))@np.array([1, 0, 0, 0])
    assert np.allclose(state, [np.sqrt(p), 0, 0, np.sqrt(1-p)])
    ks = [np.diag([1, s]), np.array([[0, np.sqrt(u)], [0, 0]])]
    generator = -0.5j*np.diag([1, -1])
    rho = np.zeros((4, 4), complex)
    drho = np.zeros_like(rho)
    for k in ks:
        v = np.kron(k, np.eye(2))@state
        dv = np.kron(k@generator, np.eye(2))@state
        rho += np.outer(v, v.conj())
        drho += np.outer(dv, v.conj())+np.outer(v, dv.conj())
    basis = [np.array([1, 0, 0, 1j])/np.sqrt(2),
             np.array([1, 0, 0, -1j])/np.sqrt(2),
             np.array([0, 1, 0, 0]), np.array([0, 0, 1, 0])]
    assert np.allclose(sum(np.outer(v, v.conj()) for v in basis), np.eye(4))
    probs = np.array([np.vdot(v, rho@v).real for v in basis])
    derivs = np.array([np.vdot(v, drho@v).real for v in basis])
    assert np.allclose(probs.sum(), 1)
    support = probs > 1e-12
    fi = np.sum(derivs[support]**2/probs[support])
    assert np.allclose(fi, 4*s*s/(1+s)**2)
