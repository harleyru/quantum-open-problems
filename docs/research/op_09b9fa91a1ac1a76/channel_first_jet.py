"""Construct a constant-Choi-rank CPTP surrogate with the same first jet.

Floating-point implementation of the accompanying exact matrix proof.
Interior, two-sided differentiability is essential to the kernel condition.
"""
import numpy as np


def first_jet_factor(j, dj, tolerance=1e-10):
    j, dj = np.asarray(j, complex), np.asarray(dj, complex)
    if j.shape != dj.shape or j.ndim != 2 or j.shape[0] != j.shape[1]:
        raise ValueError('Matching square matrices required')
    if not np.isfinite(j).all() or not np.isfinite(dj).all():
        raise ValueError('Finite matrices required')
    if not np.allclose(j, j.conj().T, atol=tolerance, rtol=0) or not np.allclose(
            dj, dj.conj().T, atol=tolerance, rtol=0):
        raise ValueError('Hermitian first jet required')
    p, u = np.linalg.eigh(j)
    if p[0] < -tolerance:
        raise ValueError('Positive Choi matrix required')
    support = p > tolerance
    c = u[:, support]*np.sqrt(p[support])
    if not c.shape[1]:
        raise ValueError('Nonzero Choi matrix required')
    projector = u[:, support]@u[:, support].conj().T
    q = np.eye(len(j))-projector
    if np.max(np.abs(q@dj@q)) > tolerance:
        raise ValueError('Kernel tangent is not compatible with a two-sided PSD curve')
    derivative = (np.eye(len(j))-projector/2)@dj@c@np.diag(1/p[support])
    return c, derivative


def surrogate_kraus(c, dc, input_dim, output_dim, theta):
    operators = [(c[:, k]+theta*dc[:, k]).reshape(input_dim, output_dim).T
                 for k in range(c.shape[1])]
    m = sum(k.conj().T@k for k in operators)
    p, u = np.linalg.eigh(m)
    if p[0] <= 0:
        raise ValueError('Surrogate normalization singular')
    inverse = (u/np.sqrt(p))@u.conj().T
    return [k@inverse for k in operators]


def choi(operators):
    c = np.column_stack([k.T.ravel() for k in operators])
    return c@c.conj().T
