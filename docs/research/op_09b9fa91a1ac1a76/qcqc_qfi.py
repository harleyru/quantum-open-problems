"""Full-gauge hypograph SDP for reduced QC-QC processes at regular points.

Numerical research implementation, not an exact optimality certificate.
Axes are I0,O0,...; Choi vectors use K.T.ravel() (unnormalized).
"""
from itertools import product
import cvxpy as cp
import numpy as np
from qcqc_feasibility import QCQCFeasibility, subsets


def hermitian_basis(r):
    basis = []
    for i in range(r):
        a = np.zeros((r, r), complex)
        a[i, i] = 1
        basis.append(a)
    for i in range(r):
        for j in range(i):
            a = np.zeros((r, r), complex)
            a[i, j] = a[j, i] = 1/np.sqrt(2)
            basis.append(a)
            a = np.zeros((r, r), complex)
            a[i, j], a[j, i] = 1j/np.sqrt(2), -1j/np.sqrt(2)
            basis.append(a)
    return basis


def ad_ensemble(n, u=0.5):
    if not isinstance(n, int) or not 1 <= n <= 3 or not 0 < u < 1:
        raise ValueError('Require 1 <= n <= 3 and 0 < u < 1')
    ks = [np.diag([1, np.sqrt(1-u)]).astype(complex),
          np.array([[0, np.sqrt(u)], [0, 0]], complex)]
    ds = [k @ (-0.5j*np.diag([1, -1])) for k in ks]
    def tensor(vectors):
        result = np.ones(1, complex)
        for v in vectors:
            result = np.kron(result, v.T.ravel())
        return result
    columns, derivatives = [], []
    for indices in product(range(2), repeat=n):
        columns.append(tensor([ks[i] for i in indices]))
        derivatives.append(sum(tensor([ds[i] if k == j else ks[i]
                                      for k, i in enumerate(indices)])
                               for j in range(n)))
    return np.column_stack(columns), np.column_stack(derivatives)


def gram_coefficients(c, dc):
    terms = np.array([dc, *[-1j*c @ h for h in hermitian_basis(c.shape[1])]])
    # G_ab = 4 Re Tr(A_a A_b^dagger S^T); the transpose is essential.
    coefficients = 4*np.einsum('air,bjr->abij', terms, terms.conj())
    return coefficients.reshape(len(terms)**2, -1)


def output_sld(s, c, dc):
    """SLD of a purified strategy output; clip solver-scale negative eigenvalues."""
    vals, vecs = np.linalg.eigh(s.T)
    root = (vecs*np.sqrt(np.maximum(vals, 0)))@vecs.conj().T
    v, dv = root@c, root@dc
    rho = v@v.conj().T
    drho = dv@v.conj().T+v@dv.conj().T
    p, u = np.linalg.eigh(rho)
    derivative = u.conj().T@drho@u
    denom = p[:, None]+p[None, :]
    mask = denom > 1e-9
    return float(np.sum(2*np.abs(derivative[mask])**2/denom[mask])), float(np.trace(rho).real)


def solve(n=1, u=0.5, eps=1e-7, tolerance=2e-5):
    c, dc = ad_ensemble(n, u)
    model = QCQCFeasibility([2]*n, [2]*n)
    size = len(c)
    s = cp.Variable((size, size), hermitian=True)
    nodes = {(used, k): cp.Variable((model.size(model.axes(used, k)),)*2,
                                   hermitian=True)
             for used in subsets(n) for k in range(n) if k not in used}
    coeff = gram_coefficients(c, dc)
    m = 1+c.shape[1]**2
    g = cp.reshape(coeff.real @ cp.vec(cp.real(s), order='C')
                   - coeff.imag @ cp.vec(cp.imag(s), order='C'), (m, m), order='C')
    t = cp.Variable()
    e = np.zeros((m, m)); e[0, 0] = 1
    constraints = [s >> 0, *[v >> 0 for v in nodes.values()]]
    constraints += [a == b for a, b in model.equations(s, nodes, True)]
    constraints += [(g+g.T)/2 - t*e >> 0]
    problem = cp.Problem(cp.Maximize(t), constraints)
    problem.solve(solver='SCS', eps=eps, max_iters=100000)
    result = {'n': n, 'u': u, 'status': problem.status, 'accepted': False}
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE) or s.value is None:
        return result
    values = {key: v.value for key, v in nodes.items()}
    audit = model.audit(s.value, values)
    numeric_g = (coeff @ s.value.ravel()).real.reshape(m, m)
    numeric_g = (numeric_g+numeric_g.T)/2
    q, b = numeric_g[1:, 1:], numeric_g[1:, 0]
    x = -np.linalg.pinv(q, rcond=1e-9) @ b
    h = sum(v*a for v, a in zip(x, hermitian_basis(c.shape[1])))
    d = dc-1j*c@h
    recomputed = float(4*np.trace((d@d.conj().T).T @ s.value).real)
    block_min = float(np.linalg.eigvalsh(numeric_g-float(t.value)*e)[0])
    stationarity = float(np.max(np.abs(q@x+b)))
    sld, normalization = output_sld(s.value, c, dc)
    result.update(audit, objective=float(t.value), gauge_recomputed=recomputed,
                  output_sld_qfi=sld, output_trace=normalization,
                  gauge_stationarity=stationarity, hypograph_min_eigenvalue=block_min,
                  process=s.value, nodes=values, gauge=h)
    result['accepted'] = bool(audit['max_equation_residual'] <= tolerance
                              and audit['min_eigenvalue'] >= -tolerance
                              and block_min >= -tolerance
                              and stationarity <= tolerance
                              and abs(normalization-1) <= tolerance
                              and abs(sld-recomputed) <= tolerance
                              and abs(recomputed-float(t.value)) <= tolerance)
    return result


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, default=1)
    parser.add_argument('--output-dir', help='Write numerical witnesses and summary locally')
    args = parser.parse_args()
    result = solve(args.n)
    summary = {k: v for k, v in result.items()
               if k not in ('process', 'nodes', 'gauge')}
    if args.output_dir:
        from pathlib import Path
        folder = Path(args.output_dir)
        folder.mkdir(parents=True, exist_ok=True)
        (folder/f'qfi_n{args.n}.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='ascii')
        if 'process' in result:
            arrays = {'process': result['process'], 'gauge': result['gauge']}
            arrays.update({f'node_{sum(1 << k for k in used)}_{next_slot}': value
                           for (used, next_slot), value in result['nodes'].items()})
            np.savez_compressed(folder/f'qfi_n{args.n}.npz', **arrays)
    print(json.dumps(summary, indent=2))
