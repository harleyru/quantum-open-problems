import numpy as np
import pytest
from qcqc_feasibility import QCQCFeasibility


@pytest.mark.parametrize('inputs,outputs,p,f', [([2],[3],2,2),
    ([2,3],[3,2],1,1), ([2,2,2],[2,2,2],1,1)])
def test_white_process(inputs, outputs, p, f):
    model = QCQCFeasibility(inputs, outputs, p, f)
    w = np.eye(model.size(model.full)) / (np.prod(inputs)*f)
    result = model.solve(w)
    assert result['accepted'], result
    assert len(result['nodes']) == len(inputs)*2**(len(inputs)-1)


def switch_process():
    # |0> -> A -> B -> F with control 0, and |0> -> B -> A -> F with control 1.
    v = np.zeros((2,2,2,2,2,2), complex)
    for a in range(2):
        for b in range(2):
            v[0,a,a,b,b,0] += 1/np.sqrt(2)
            v[b,a,0,b,a,1] += 1/np.sqrt(2)
    v = v.ravel()
    return np.outer(v, v.conj())


def test_quantum_switch_with_future_coherence():
    w = switch_process()
    assert np.linalg.matrix_rank(w) == 1
    assert np.linalg.norm(w[::2,1::2]) > 0
    result = QCQCFeasibility([2,2],[2,2], future=4).solve(w)
    assert result['accepted'], result


def test_erratum_negative_future_with_valid_marginal():
    model = QCQCFeasibility([2],[2], future=2)
    marginal = np.eye(4)/2
    good = np.kron(marginal, np.eye(2)/2)
    bad = np.kron(marginal, np.diag([2.,-1.]))
    np.testing.assert_allclose(model.trace(good, model.full, 'F'),
                               model.trace(bad, model.full, 'F'))
    assert model.solve(good)['accepted']
    result = model.solve(bad)
    assert not result['accepted']
    assert result['status'] == 'invalid_nonpositive'


def test_positive_normalized_but_signalling_one_slot_rejected():
    # Correct total trace is not enough: output dependence violates the terminal equation.
    model = QCQCFeasibility([2],[2])
    w = np.kron(np.eye(2)/2, np.diag([2.,0.]))
    result = model.solve(w)
    assert not result['accepted']
    assert result['status'] in ('infeasible', 'infeasible_inaccurate')


def test_ocb_process_not_qcqc():
    eye = np.eye(2)
    x = np.array([[0.,1.],[1.,0.]])
    z = np.diag([1.,-1.])
    w = (np.eye(16) + (np.kron(np.kron(np.kron(eye,z),z),eye)
         + np.kron(np.kron(np.kron(z,eye),x),z))/np.sqrt(2))/4
    assert np.linalg.eigvalsh(w)[0] > -1e-12
    result = QCQCFeasibility([2,2],[2,2]).solve(w)
    assert not result['accepted']
    assert result['status'] in ('infeasible', 'infeasible_inaccurate')


def test_permutation_and_partial_trace_against_index_loops():
    model = QCQCFeasibility([2,3],[3,2], past=2)
    axes = model.axes(frozenset({1}), 0)
    rng = np.random.default_rng(123)
    a = rng.normal(size=(model.size(axes),)*2)
    target = model.axes(frozenset({0,1}))
    result = model.insert_identity(a, axes, 'O0', target)
    dims = tuple(model.dim[t] for t in target)
    for _ in range(100):
        i = tuple(rng.integers(d) for d in dims)
        j = tuple(rng.integers(d) for d in dims)
        si = tuple(i[target.index(t)] for t in axes)
        sj = tuple(j[target.index(t)] for t in axes)
        expected = a[np.ravel_multi_index(si, tuple(model.dim[t] for t in axes)),
                     np.ravel_multi_index(sj, tuple(model.dim[t] for t in axes))]
        if i[target.index('O0')] != j[target.index('O0')]:
            expected = 0
        assert result[np.ravel_multi_index(i,dims), np.ravel_multi_index(j,dims)] == expected
    np.testing.assert_allclose(model.trace(result,target,'O0'),3*a)


def test_complex_switch_and_corrupted_certificate():
    phase = np.tile([1.,1j],32)
    w = phase[:,None]*switch_process()*phase.conj()[None,:]
    model = QCQCFeasibility([2,2],[2,2], future=4)
    result = model.solve(w)
    assert result['accepted'], result
    nodes = {k:v.copy() for k,v in result['nodes'].items()}
    nodes[frozenset(),0] *= 2
    assert model.audit(w,nodes)['max_equation_residual'] > 0.1


def test_identity_wires_nontrivial_past_and_future():
    phi = np.eye(2).ravel()
    v = np.kron(phi,phi)
    model = QCQCFeasibility([2],[2],past=2,future=2)
    assert model.solve(np.outer(v,v))['accepted']


@pytest.mark.parametrize('order', [(0,1,2),(2,1,0),(1,2,0)])
def test_three_slot_nonwhite_fixed_orders(order):
    # In route order: initial |0>, two identity wires, discarded last output.
    phi = np.eye(2).ravel()
    v = np.kron(np.kron([1.,0.],phi),phi)
    route_w = np.kron(np.outer(v,v),np.eye(2))
    source = tuple(a for k in order for a in (f'I{k}',f'O{k}'))
    target = ('I0','O0','I1','O1','I2','O2')
    axes = tuple(source.index(a) for a in target)
    w = route_w.reshape((2,)*12).transpose(axes+tuple(a+6 for a in axes)).reshape(64,64)
    model = QCQCFeasibility([2,2,2],[2,2,2])
    assert model.solve(w)['accepted']


def test_bad_inputs():
    with pytest.raises(ValueError):
        QCQCFeasibility([0],[2])
    model = QCQCFeasibility([2],[2])
    with pytest.raises(ValueError):
        model.solve(np.eye(3))
    w = np.eye(4, dtype=complex)/2
    w[0,1] = 1j
    assert model.solve(w)['status'] == 'invalid_nonhermitian'
