import numpy as np
from qcqc_feasibility import QCQCFeasibility, subsets
from qcqc_support_dual import slacks, solve_support


def test_identity_support():
    result = solve_support(2, np.eye(16))
    assert abs(result['objective']-4) < 1e-6
    assert result['min_slack_eigenvalue'] > -1e-6
    assert not result['rigorous']


def test_adjoint_telescoping_arbitrary_matrices():
    model = QCQCFeasibility([2, 3], [3, 2])
    rng = np.random.default_rng(20)
    def hermitian(size):
        a = rng.normal(size=(size, size))+1j*rng.normal(size=(size, size))
        return a+a.conj().T
    y = {used: hermitian(model.size(model.axes(used))) for used in subsets(2)}
    nodes = {(used, k): hermitian(model.size(model.axes(used, k)))
             for used in subsets(2) for k in range(2) if k not in used}
    w = hermitian(model.size(model.full))
    equations = model.equations(w, nodes)
    lhs = sum(np.trace(y[used]@(a-b))
              for used, (a, b) in zip(subsets(2), equations))
    rhs = sum(np.trace(a@nodes[key]) for key, a in slacks(model, y).items())
    rhs += np.trace(y[frozenset([0, 1])]@w)-np.trace(y[frozenset()])
    assert np.allclose(lhs, rhs)
