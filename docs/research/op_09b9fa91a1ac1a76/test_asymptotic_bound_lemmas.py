from itertools import product
from pathlib import Path
import numpy as np
from extract_qcqc_circuit import load_repaired


def ensemble(slots):
    columns = []
    for choices in product(*slots):
        column = np.ones(1, complex)
        for k in choices:
            column = np.kron(column, k.T.ravel())
        columns.append(column)
    return np.column_stack(columns)


def test_complex_tangent_domination_and_cross_cancellation():
    path = Path(__file__).parent/'qfi_results/qfi_n2.npz'
    s, _ = load_repaired(path, 2)
    rng = np.random.default_rng(164)
    for _ in range(4):
        v, _ = np.linalg.qr(rng.normal(size=(4, 2))+1j*rng.normal(size=(4, 2)))
        h = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
        h = h+h.conj().T
        dv = -1j*h@v
        beta = dv.conj().T@v
        p = v@beta.conj().T
        r = dv-p
        ks = list(v.reshape(2, 2, 2))
        ps = list(p.reshape(2, 2, 2))
        rs = list(r.reshape(2, 2, 2))
        def inner(a, b):
            return np.trace(a.conj().T@s.T@b)
        cc = ensemble([ks, ks])
        assert np.allclose(inner(cc, cc), 1)
        rr = [ensemble([rs, ks]), ensemble([ks, rs])]
        pp = [ensemble([ps, ks]), ensemble([ks, ps])]
        assert abs(inner(rr[0], rr[1])) < 1e-10
        residual_bound = np.linalg.norm(r.conj().T@r, 2)
        for a, b in zip(rr, pp):
            assert inner(a, a).real <= residual_bound+1e-10
            assert inner(b, b).real <= np.linalg.norm(beta, 2)**2+1e-10
        total = sum(rr+pp)
        bound = (2*np.linalg.norm(beta, 2)+np.sqrt(2*residual_bound))**2
        assert inner(total, total).real <= bound+1e-10
