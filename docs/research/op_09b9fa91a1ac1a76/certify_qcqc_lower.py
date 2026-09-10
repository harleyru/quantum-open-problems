"""Interval full-gauge hypograph certificate for an exactly repaired process."""
from fractions import Fraction
import numpy as np
from mpmath import iv
from certify_qcqc_dual import conj, exact_ensemble, positive_ldl
from repair_qcqc_primal import repair, verify, as_interval


def sparse_terms(c, dc):
    def nonzero(z):
        return not (z.real == 0 and z.imag == 0)
    terms = [[(i, k, dc[i, k]) for i in range(len(c)) for k in range(c.shape[1])
              if nonzero(dc[i, k])]]
    r = c.shape[1]
    basis = []
    for a in range(r):
        basis.append([(a, a, iv.mpc(1))])
    for a in range(r):
        for b in range(a):
            basis.append([(a, b, iv.mpc(1)), (b, a, iv.mpc(1))])
            basis.append([(a, b, iv.mpc(0, 1)), (b, a, iv.mpc(0, -1))])
    for h in basis:
        terms.append([(i, k, -iv.mpc(0, 1)*c[i, j]*v)
                      for j, k, v in h for i in range(len(c)) if nonzero(c[i, j])])
    return terms


def hypograph(s, n, lower):
    c, dc = exact_ensemble(n)
    terms = sparse_terms(c, dc)
    g = np.empty((len(terms), len(terms)), object)
    for a, aa in enumerate(terms):
        for b in range(a+1):
            value = 4*sum(v*conj(w)*s[i, j]
                          for i, k, v in aa for j, ell, w in terms[b] if k == ell)
            g[a, b] = g[b, a] = iv.mpc(value.real)
    g[0, 0] -= iv.mpf(lower.numerator)/lower.denominator
    return g


def certify(path, n, lower='4.736'):
    iv.dps = 50
    lower = Fraction(lower)
    process, parts = repair(path, n)
    feasible = verify(process, parts)
    hypograph_check = positive_ldl(hypograph(as_interval(*process), n, lower))
    return {'verified': feasible['verified'] and hypograph_check[0],
            'n': n, 'lower_bound_rational': str(lower),
            'mixing': '1/10000', 'precision': 50,
            'feasibility': feasible, 'hypograph_check': hypograph_check}


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, choices=[2, 3], required=True)
    parser.add_argument('--lower', required=True)
    args = parser.parse_args()
    folder = Path(__file__).parent/'qfi_results'
    result = certify(folder/f'qfi_n{args.n}.npz', args.n, args.lower)
    (folder/f'certified_lower_n{args.n}.json').write_text(json.dumps(result, indent=2)+'\n', encoding='ascii')
    print(json.dumps(result, indent=2))
