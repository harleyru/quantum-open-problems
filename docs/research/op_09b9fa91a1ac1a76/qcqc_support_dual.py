"""Numerical fixed-gauge support dual. Outputs are candidates, not certificates."""
import cvxpy as cp
import numpy as np
from qcqc_feasibility import QCQCFeasibility, subsets
from qcqc_qfi import ad_ensemble


def slacks(model, y, symbolic=False):
    result = {}
    for used in subsets(model.n):
        for k in range(model.n):
            if k in used:
                continue
            result[used, k] = model.insert_identity(
                y[used], model.axes(used), f'I{k}', model.axes(used, k), symbolic
            ) - model.trace(y[used | {k}], model.axes(used | {k}), f'O{k}', symbolic)
    return result


def solve_support(n, omega):
    model = QCQCFeasibility([2]*n, [2]*n)
    size = model.size(model.full)
    omega = np.asarray(omega, complex)
    if omega.shape != (size, size) or not np.isfinite(omega).all():
        raise ValueError('Invalid performance operator')
    if not np.allclose(omega, omega.conj().T, atol=1e-10, rtol=0):
        raise ValueError('Performance operator must be Hermitian')
    full = frozenset(range(n))
    y = {used: cp.Variable((model.size(model.axes(used)),)*2, hermitian=True)
         for used in subsets(n) if used != full}
    y[full] = cp.Constant(omega)
    constraints = [a >> 0 for a in slacks(model, y, True).values()]
    problem = cp.Problem(cp.Minimize(cp.real(cp.trace(y[frozenset()]))), constraints)
    problem.solve(solver='SCS', eps=1e-8, max_iters=100000)
    result = {'status': problem.status, 'rigorous': False}
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        return result
    values = {key: a.value for key, a in y.items()}
    residuals = slacks(model, values)
    result.update(objective=float(problem.value), matrices=values,
                  min_slack_eigenvalue=min(float(np.linalg.eigvalsh(a)[0])
                                          for a in residuals.values()))
    return result


def solve_joint(n, u=0.5):
    c, dc = ad_ensemble(n, u)
    model = QCQCFeasibility([2]*n, [2]*n)
    y = {used: cp.Variable((model.size(model.axes(used)),)*2, hermitian=True)
         for used in subsets(n)}
    h = cp.Variable((c.shape[1],)*2, hermitian=True)
    d = dc-1j*c@h
    full = frozenset(range(n))
    # Schur complement: Y_all >= 4 conjugate(D) D.T = Omega(h).
    block = cp.bmat([[y[full], 2*cp.conj(d)],
                     [2*d.T, np.eye(c.shape[1])]])
    constraints = [a >> 0 for a in slacks(model, y, True).values()]+[block >> 0]
    problem = cp.Problem(cp.Minimize(cp.real(cp.trace(y[frozenset()]))), constraints)
    problem.solve(solver='SCS', eps=1e-8, max_iters=100000)
    result = {'status': problem.status, 'rigorous': False}
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        return result
    values = {key: a.value for key, a in y.items()}
    dd = dc-1j*c@h.value
    terminal = values[full]-4*(dd@dd.conj().T).T
    result.update(objective=float(problem.value), matrices=values, gauge=h.value,
                  min_slack_eigenvalue=min(float(np.linalg.eigvalsh(a)[0])
                                          for a in slacks(model, values).values()),
                  terminal_min_eigenvalue=float(np.linalg.eigvalsh(terminal)[0]))
    return result


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, required=True, choices=[2, 3])
    parser.add_argument('--folder', default=str(Path(__file__).parent/'qfi_results'))
    parser.add_argument('--joint', action='store_true')
    args = parser.parse_args()
    folder = Path(args.folder)
    if args.joint:
        result = solve_joint(args.n)
        h = result.get('gauge')
    else:
        with np.load(folder/f'qfi_n{args.n}.npz') as data:
            h = data['gauge']
        c, dc = ad_ensemble(args.n)
        d = dc-1j*c@h
        result = solve_support(args.n, 4*(d@d.conj().T).T)
    prefix = 'joint_dual' if args.joint else 'dual'
    summary = {key: v for key, v in result.items() if key not in ('matrices', 'gauge')}
    (folder/f'{prefix}_n{args.n}.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='ascii')
    if 'matrices' in result:
        np.savez_compressed(folder/f'{prefix}_n{args.n}.npz', gauge=h,
                            **{f'y_{sum(1 << k for k in used)}': a
                               for used, a in result['matrices'].items()})
    print(json.dumps(summary, indent=2))
