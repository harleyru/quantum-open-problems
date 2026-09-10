"""Wechs et al., corrected Proposition 7 / Eq. (63), fixed-W feasibility.

Unnormalized Choi convention; axes P,I0,O0,...,I(N-1),O(N-1),F.
Numerical feasibility is not an exact or interval certificate.
"""
from itertools import combinations
from math import prod
import numpy as np
import cvxpy as cp
from scipy import sparse


def subsets(n):
    return [frozenset(c) for size in range(n + 1)
            for c in combinations(range(n), size)]


class QCQCFeasibility:
    def __init__(self, inputs, outputs, past=1, future=1):
        inputs, outputs = tuple(inputs), tuple(outputs)
        dims = (past, *inputs, *outputs, future)
        if not inputs or len(inputs) != len(outputs) or any(
                not isinstance(d, (int, np.integer)) or d < 1 for d in dims):
            raise ValueError('Positive integer dimensions and at least one slot required')
        self.n = len(inputs)
        self.dim = {'P': past, 'F': future}
        for k, (di, do) in enumerate(zip(inputs, outputs)):
            self.dim[f'I{k}'], self.dim[f'O{k}'] = di, do
        self.full = self.axes(frozenset(range(self.n))) + ('F',)

    def axes(self, used, next_slot=None):
        result = ['P']
        for k in range(self.n):
            if k in used:
                result.extend((f'I{k}', f'O{k}'))
            elif k == next_slot:
                result.append(f'I{k}')
        return tuple(result)

    def size(self, axes):
        return prod(self.dim[a] for a in axes)

    def trace(self, matrix, axes, label, symbolic=False):
        dims = [self.dim[a] for a in axes]
        axis = axes.index(label)
        if symbolic:
            return (cp.partial_trace(cp.real(matrix), tuple(dims), axis=axis)
                    + 1j*cp.partial_trace(cp.imag(matrix), tuple(dims), axis=axis))
        tensor = np.asarray(matrix).reshape(dims + dims)
        return np.trace(tensor, axis1=axis, axis2=axis + len(dims)).reshape(
            self.size(axes) // self.dim[label], -1)

    def insert_identity(self, matrix, axes, label, target, symbolic=False):
        source = axes + (label,)
        shape = tuple(self.dim[a] for a in source)
        order = tuple(source.index(a) for a in target)
        indices = np.arange(prod(shape)).reshape(shape).transpose(order).ravel()
        if symbolic:
            perm = sparse.csr_matrix((np.ones(len(indices)),
                                      (np.arange(len(indices)), indices)),
                                     shape=(len(indices), len(indices)))
            return perm @ cp.kron(matrix, np.eye(self.dim[label])) @ perm.T
        extended = np.kron(matrix, np.eye(self.dim[label]))
        return extended[np.ix_(indices, indices)]

    def equations(self, w, nodes, symbolic=False):
        equations = []
        all_slots = frozenset(range(self.n))
        for used in subsets(self.n):
            target = self.axes(used)
            if used == all_slots:
                left = self.trace(w, self.full, 'F', symbolic)
            else:
                left = sum(self.trace(nodes[used, k], self.axes(used, k),
                                      f'I{k}', symbolic)
                           for k in range(self.n) if k not in used)
            right = np.eye(self.dim['P']) if not used else sum(
                self.insert_identity(nodes[used - {k}, k],
                                     self.axes(used - {k}, k), f'O{k}',
                                     target, symbolic) for k in sorted(used))
            equations.append((left, right))
        return equations

    def audit(self, w, nodes):
        matrices = [w, *nodes.values()]
        return {
            'max_equation_residual': max(float(np.max(np.abs(a - b)))
                                         for a, b in self.equations(w, nodes)),
            'min_eigenvalue': min(float(np.linalg.eigvalsh((a + a.conj().T)/2)[0])
                                 for a in matrices),
            'max_hermiticity_residual': max(float(np.max(np.abs(a-a.conj().T)))
                                           for a in matrices),
        }

    def solve(self, w, tolerance=2e-6, solver='SCS', **options):
        if not np.isfinite(tolerance) or tolerance <= 0:
            raise ValueError('Positive finite tolerance required')
        w = np.asarray(w, dtype=complex)
        size = self.size(self.full)
        if w.shape != (size, size) or not np.isfinite(w).all():
            raise ValueError('W has wrong shape or nonfinite entries')
        if np.max(np.abs(w-w.conj().T)) > tolerance:
            return {'status': 'invalid_nonhermitian', 'accepted': False}
        if np.linalg.eigvalsh((w+w.conj().T)/2)[0] < -tolerance:
            return {'status': 'invalid_nonpositive', 'accepted': False}
        nodes = {(used, k): cp.Variable((self.size(self.axes(used, k)),)*2,
                                       hermitian=True)
                 for used in subsets(self.n) for k in range(self.n) if k not in used}
        constraints = [cp.Constant((w+w.conj().T)/2) >> 0]
        constraints += [v >> 0 for v in nodes.values()]
        constraints += [a == b for a, b in self.equations(w, nodes, True)]
        problem = cp.Problem(cp.Minimize(0), constraints)
        if solver == 'SCS':
            options = {'eps': 1e-8, 'max_iters': 30000, **options}
        problem.solve(solver=solver, **options)
        result = {'status': problem.status, 'accepted': False, 'nodes': {}}
        if problem.status in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE) and all(
                v.value is not None for v in nodes.values()):
            values = {key: v.value for key, v in nodes.items()}
            audit = self.audit(w, values)
            result.update(audit, nodes=values)
            result['accepted'] = (audit['max_equation_residual'] <= tolerance
                                  and audit['min_eigenvalue'] >= -tolerance
                                  and audit['max_hermiticity_residual'] <= tolerance)
        return result
