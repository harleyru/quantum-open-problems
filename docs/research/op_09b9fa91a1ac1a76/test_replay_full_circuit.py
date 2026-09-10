from pathlib import Path
import numpy as np
from complete_qcqc_circuit import load_operations
from replay_full_circuit import replay
from extract_qcqc_circuit import simulate


def test_nonzero_completion_branch_is_retained():
    empty, full = frozenset(), frozenset([0])
    ops = {empty: {'map': np.array([[0], [1]], complex), 'incoming': [], 'outgoing': [(0, 2)]},
           full: {'map': np.array([[1, 0]], complex), 'incoming': [(0, 2)], 'outgoing': []}}
    nulls = {empty: np.empty((0, 1)), full: np.array([[0, 1]], complex)}
    rho, _ = replay(ops, nulls, 1)
    principal, _ = simulate(ops, 1)
    assert np.allclose(rho, [[1]])
    assert np.allclose(principal, [[0.5]])


def test_exported_full_replay_derivative():
    folder = Path(__file__).parent/'qfi_results'
    ops, nulls = load_operations(folder, 3, completed=True)
    rho, drho = replay(ops, nulls, 3, theta=0.31)
    plus, _ = replay(ops, nulls, 3, theta=0.310001)
    minus, _ = replay(ops, nulls, 3, theta=0.309999)
    assert abs(np.trace(rho)-1) < 1e-12
    assert np.allclose(drho, (plus-minus)/2e-6, atol=1e-8)
