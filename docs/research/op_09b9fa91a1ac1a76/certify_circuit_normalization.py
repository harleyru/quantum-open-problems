"""Integer-arithmetic Gram and cumulative physical-normalization bounds.

Bounds compare exactly defined rational circuits and their exact polar
normalizations, not floating replay output and not QFI values.
"""
from fractions import Fraction
import numpy as np


def rational_grid(a, denominator=10**15):
    def component(values):
        return np.array([[int(round(Fraction.from_float(float(x))*denominator))
                          for x in row] for row in values], dtype=object)
    return component(a.real), component(a.imag)


def gram_bound(real, imag, denominator):
    if real.shape != imag.shape:
        raise ValueError('Matching real and imaginary arrays required')
    size = real.shape[1]
    gr = real.T@real+imag.T@imag-denominator**2*np.eye(size, dtype=object)
    gi = real.T@imag-imag.T@real
    rows = [sum(abs(int(x))+abs(int(y)) for x, y in zip(rr, ii))
            for rr, ii in zip(gr, gi)]
    return Fraction(max(rows), denominator**2)


def certify(path, n, denominator=10**15):
    epsilons, arrays = {}, {}
    with np.load(path) as data:
        for mask in range(2**n):
            stacked = np.vstack([data[f'map_{mask}'], data[f'null_{mask}']])
            re, im = rational_grid(stacked, denominator)
            epsilon = gram_bound(re, im, denominator)
            if epsilon >= 1:
                raise ArithmeticError('Gram error too large for the normalization bound')
            epsilons[mask] = epsilon
            arrays[f'real_{mask}'], arrays[f'imag_{mask}'] = re, im
    layers = [max(e for mask, e in epsilons.items() if mask.bit_count() == k)
              for k in range(n+1)]
    growth = Fraction(1)
    for epsilon in layers:
        growth *= 1+epsilon
    cumulative = 3*sum(layers)*growth
    result = {'n': n, 'denominator': denominator, 'integer_arithmetic': True,
              'gram_bounds': {str(mask): str(e) for mask, e in epsilons.items()},
              'layer_bounds': [str(e) for e in layers],
              'output_trace_norm_bound': str(cumulative),
              'output_bound_below_1e_minus_10': cumulative < Fraction(1, 10**10),
              'scope': 'exact rational raw circuit versus exact polar-normalized circuit; not QFI or replay-roundoff bound'}
    return result, arrays


if __name__ == '__main__':
    import json
    from pathlib import Path
    folder = Path(__file__).parent/'qfi_results'
    result, arrays = certify(folder/'completed_circuit_n3.npz', 3)
    np.savez_compressed(folder/'rational_circuit_n3.npz',
                        **{key: np.asarray(a, dtype=str) for key, a in arrays.items()})
    (folder/'circuit_normalization_certificate_n3.json').write_text(
        json.dumps(result, indent=2)+'\n', encoding='ascii')
    print(json.dumps(result, indent=2))
