from pathlib import Path
import numpy as np
from extract_qcqc_circuit import extract, load_repaired, root, simulate


def circuit():
    path = Path(__file__).parent/'qfi_results/qfi_n3.npz'
    s, nodes = load_repaired(path, 3)
    ops, diagnostics = extract(s, nodes, 3)
    return s, ops, diagnostics


def test_effective_isometry_and_different_complex_channels():
    s, ops, diagnostics = circuit()
    assert max(d['isometry_residual'] for d in diagnostics.values()) < 1e-9
    rng = np.random.default_rng(92)
    for _ in range(3):
        channels = []
        c = np.ones(1, complex)
        for k in range(3):
            u, _ = np.linalg.qr(rng.normal(size=(2, 2))+1j*rng.normal(size=(2, 2)))
            channels.append(([u], [np.zeros((2, 2))]))
            c = np.kron(c, u.T.ravel())
        rho, drho = simulate(ops, 3, channels=channels)
        v = root(s).T@c
        assert np.allclose(rho, np.outer(v, v.conj()), atol=1e-10)
        assert abs(np.trace(rho)-1) < 1e-10
        assert np.max(np.abs(drho)) == 0


def test_channel_derivative_finite_difference_and_corrupted_map():
    _, ops, _ = circuit()
    rho, drho = simulate(ops, 3, theta=0.37)
    plus, _ = simulate(ops, 3, theta=0.370001)
    minus, _ = simulate(ops, 3, theta=0.369999)
    assert np.allclose(drho, (plus-minus)/2e-6, atol=1e-8)
    assert abs(np.trace(rho)-1) < 1e-10
    ops[frozenset()]['map'] *= 0.9
    bad, _ = simulate(ops, 3)
    assert abs(np.trace(bad)-1) > 0.1
