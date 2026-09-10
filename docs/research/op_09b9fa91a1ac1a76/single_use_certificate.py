"""Exact one-query amplitude-damping metrology certificate.

Symbolic identities use s=sqrt(1-u), 0<s<1. Rational enclosures use only
integer arithmetic and certify F(1/2)=12-8*sqrt(2), not solver output.
"""
from fractions import Fraction
from math import isqrt
import sympy as sp


def symbolic_certificate():
    s = sp.Symbol('s', positive=True)
    h0 = (s-1)/(2*(s+1))
    h1 = sp.Rational(1, 2)
    alpha0 = (sp.Rational(1, 2)+h0)**2
    alpha1 = s**2*(-sp.Rational(1, 2)+h0)**2 + (1-s**2)*(-sp.Rational(1, 2)+h1)**2
    optimum = 4*s**2/(1+s)**2
    p = s/(1+s)
    weight = p+s**2*(1-p)
    classical_fi = 4*p*s**2*(1-p)/weight
    return {
        'alpha0_minus_bound': sp.factor(4*alpha0-optimum),
        'alpha1_minus_bound': sp.factor(4*alpha1-optimum),
        'measurement_fi_minus_bound': sp.factor(classical_fi-optimum),
        'normalization_minus_one': sp.factor(weight+(1-s**2)*(1-p)-1),
    }


def half_noise_interval(digits=30):
    if not isinstance(digits, int) or not 1 <= digits <= 1000:
        raise ValueError('digits must be an integer from 1 to 1000')
    scale = 10**digits
    k = isqrt(2*scale*scale)
    assert k*k < 2*scale*scale < (k+1)*(k+1)
    return Fraction(12)-8*Fraction(k+1, scale), Fraction(12)-8*Fraction(k, scale)


if __name__ == '__main__':
    checks = symbolic_certificate()
    if any(value != 0 for value in checks.values()):
        raise RuntimeError(checks)
    lo, hi = half_noise_interval()
    print('All four symbolic certificate identities vanish exactly.')
    print(f'Exact strict lower bound: {lo}')
    print(f'Exact strict upper bound: {hi}')
    print(f'Interval width: {hi-lo}')
