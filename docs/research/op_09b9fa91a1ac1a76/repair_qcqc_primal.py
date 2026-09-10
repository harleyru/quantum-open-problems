"""Exact rational flow repair and interval positivity of QC-QC witnesses."""
from fractions import Fraction as F
from math import factorial
import numpy as np
from mpmath import iv
from qcqc_feasibility import QCQCFeasibility, subsets
from certify_qcqc_dual import positive_ldl


class RationalFlow(QCQCFeasibility):
    def insert_identity(self, matrix, axes, label, target, symbolic=False):
        if symbolic:
            raise ValueError('Exact rational implementation only')
        source = axes+(label,)
        dims = tuple(self.dim[a] for a in source)
        order = tuple(source.index(a) for a in target)
        indices = np.arange(np.prod(dims)).reshape(dims).transpose(order).ravel()
        extended = np.kron(matrix, np.eye(self.dim[label], dtype=int))
        return extended[np.ix_(indices, indices)]


def rational_parts(a):
    re, im = np.empty(a.shape, object), np.empty(a.shape, object)
    for i in range(len(a)):
        re[i, i], im[i, i] = F(str(float(a[i, i].real))), F(0)
        for j in range(i):
            re[i, j] = re[j, i] = F(str(float(a[i, j].real)))
            im[i, j] = F(str(float(a[i, j].imag)))
            im[j, i] = -im[i, j]
    return re, im


def as_interval(re, im):
    def scalar(x):
        x = F(x)
        return iv.mpf(int(x.numerator))/int(x.denominator)
    return np.array([[iv.mpc(scalar(a), scalar(b)) for a, b in zip(rr, ii)]
                     for rr, ii in zip(re, im)], dtype=object)


def repair(path, n, mixing=F(1, 10000)):
    if not 0 < mixing < 1:
        raise ValueError('mixing must lie strictly between zero and one')
    model = RationalFlow([2]*n, [2]*n)
    all_slots = frozenset(range(n))
    with np.load(path) as data:
        pairs = {(used, k): rational_parts(data[f'node_{sum(1 << j for j in used)}_{k}'])
                 for used in subsets(n) for k in range(n) if k not in used}
    parts, process = [], []
    for component in range(2):
        nodes = {key: pair[component].copy() for key, pair in pairs.items()}
        for used in subsets(n):
            target = model.axes(used)
            rhs = (np.array([[F(1-component)]], object) if not used else sum(
                model.insert_identity(nodes[used-{k}, k], model.axes(used-{k}, k),
                                      f'O{k}', target) for k in sorted(used)))
            if used == all_slots:
                process.append(rhs)
                break
            remaining = [k for k in range(n) if k not in used]
            lhs = sum(model.trace(nodes[used, k], model.axes(used, k), f'I{k}')
                      for k in remaining)
            error = rhs-lhs
            for k in remaining:
                nodes[used, k] += model.insert_identity(
                    error*F(1, 2*len(remaining)), target, f'I{k}', model.axes(used, k))
        parts.append(nodes)
    for component in range(2):
        for (used, k), a in parts[component].items():
            m = len(used)
            white = F(factorial(m)*factorial(n-m-1), factorial(n)*2**(m+1))
            parts[component][used, k] = (1-mixing)*a
            if component == 0:
                parts[component][used, k] += mixing*white*np.eye(len(a), dtype=int)
        process[component] = (1-mixing)*process[component]
        if component == 0:
            process[component] += mixing*F(1, 2**n)*np.eye(4**n, dtype=int)
    # Verify all equations over exact rationals; no tolerance enters this gate.
    for component in range(2):
        for index, (a, b) in enumerate(model.equations(process[component], parts[component])):
            if index == 0:
                b = np.array([[F(1-component)]], object)
            if any(value != 0 for value in (a-b).ravel()):
                raise ArithmeticError('Exact flow repair failed')
    return process, parts


def verify(process, parts):
    checks = {str(key): positive_ldl(as_interval(parts[0][key], parts[1][key]))
              for key in parts[0]}
    checks['process'] = positive_ldl(as_interval(*process))
    return {'verified': all(ok for ok, _ in checks.values()), 'checks': checks}


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, choices=[2, 3], required=True)
    args = parser.parse_args()
    iv.dps = 50
    folder = Path(__file__).parent/'qfi_results'
    process, parts = repair(folder/f'qfi_n{args.n}.npz', args.n)
    result = verify(process, parts)
    result.update(n=args.n, mixing='1/10000', exact_flow=True, precision=50)
    arrays = {'process_real': process[0], 'process_imag': process[1]}
    for component, suffix in enumerate(('real', 'imag')):
        arrays.update({f'node_{sum(1 << k for k in used)}_{next_slot}_{suffix}': a
                       for (used, next_slot), a in parts[component].items()})
    # Fraction strings avoid pickle and preserve exact values for independent replay.
    np.savez_compressed(folder/f'repaired_primal_n{args.n}.npz',
                        **{key: np.array([[str(x) for x in row] for row in a])
                           for key, a in arrays.items()})
    (folder/f'repaired_primal_n{args.n}.json').write_text(json.dumps(result, indent=2)+'\n', encoding='ascii')
    print(json.dumps(result, indent=2))
