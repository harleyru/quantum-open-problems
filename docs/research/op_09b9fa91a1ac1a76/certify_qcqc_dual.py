"""Outward interval LDL verification of decimal-rational dual witnesses.

Uses mpmath interval arithmetic; no floating eigenvalues enter acceptance.
The channel is exactly half damping. Failure means unverified, not invalid.
"""
from itertools import product
from fractions import Fraction
import numpy as np
from mpmath import iv
from qcqc_feasibility import QCQCFeasibility
from qcqc_support_dual import slacks


def conj(z):
    return iv.mpc(z.real, -z.imag)


def hermitian_decimal(a):
    out = np.empty(a.shape, object)
    for i in range(len(a)):
        out[i, i] = iv.mpc(str(float(a[i, i].real)), '0')
        for j in range(i):
            out[i, j] = iv.mpc(str(float(a[i, j].real)), str(float(a[i, j].imag)))
            out[j, i] = conj(out[i, j])
    return out


def positive_ldl(a):
    n = len(a)
    l = [[iv.mpc(0) for _ in range(n)] for _ in range(n)]
    pivots = []
    for j in range(n):
        pivot = a[j, j]-sum(l[j][k]*conj(l[j][k])*pivots[k] for k in range(j))
        if not (pivot.real.a > 0 and pivot.imag.a <= 0 <= pivot.imag.b):
            return False, j
        pivots.append(pivot.real)
        for i in range(j+1, n):
            l[i][j] = (a[i, j]-sum(l[i][k]*conj(l[j][k])*pivots[k]
                                   for k in range(j)))/pivot.real
    return True, n


def exact_ensemble(n):
    s = iv.sqrt(iv.mpf('0.5'))
    zero, one = iv.mpc(0), iv.mpc(1)
    local = [[one, zero, zero, iv.mpc(s)], [zero, zero, iv.mpc(s), zero]]
    c = np.empty((4**n, 2**n), object)
    dc = np.empty_like(c)
    for j, indices in enumerate(product(range(2), repeat=n)):
        for row, axes in enumerate(product(range(4), repeat=n)):
            v = one
            for k, axis in zip(indices, axes):
                v = v*local[k][axis]
            c[row, j] = v
            dc[row, j] = v*iv.mpc(0, sum(1 if a >= 2 else -1 for a in axes))/2
    return c, dc


def certify(path, n, margin='0.000001', precision=50):
    iv.dps = precision
    delta = Fraction(margin)
    if delta <= 0:
        raise ValueError('Positive repair margin required')
    model = QCQCFeasibility([2]*n, [2]*n)
    with np.load(path) as data:
        h = hermitian_decimal(data['gauge'])
        y = {}
        for mask in range(2**n):
            used = frozenset(k for k in range(n) if mask & (1 << k))
            a = hermitian_decimal(data[f'y_{mask}'])
            shift = (2**(n-len(used)+1)-1)*iv.mpf(margin)
            for i in range(len(a)):
                a[i, i] += shift
            y[used] = a
        objective = Fraction(str(float(data['y_0'][0, 0].real)))+(2**(n+1)-1)*delta
    c, dc = exact_ensemble(n)
    d = dc-(c@h)*iv.mpc(0, 1)
    omega = np.empty((len(c), len(c)), object)
    for i in range(len(c)):
        for j in range(len(c)):
            omega[i, j] = 4*sum(conj(d[i, k])*d[j, k] for k in range(d.shape[1]))
    matrices = slacks(model, y)
    matrices['terminal'] = y[frozenset(range(n))]-omega
    checks = {str(key): positive_ldl(a) for key, a in matrices.items()}
    return {'verified': all(a for a, _ in checks.values()), 'n': n,
            'precision': precision, 'repair_margin': margin,
            'upper_bound_rational': str(objective), 'checks': checks}


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, choices=[2, 3], required=True)
    args = parser.parse_args()
    folder = Path(__file__).parent/'qfi_results'
    result = certify(folder/f'joint_dual_n{args.n}.npz', args.n)
    (folder/f'certified_dual_n{args.n}.json').write_text(json.dumps(result, indent=2)+'\n', encoding='ascii')
    print(json.dumps(result, indent=2))
