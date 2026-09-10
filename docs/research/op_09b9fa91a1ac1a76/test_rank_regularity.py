"""Pointwise SLD versus smooth-lift cost at an explicit rank change.

These are preparation channels, not QC-QC advantage witnesses.
Analytic derivatives avoid hiding the singularity in finite differences.
"""
import numpy as np
import pytest


def preparation(t, rotate=False):
    columns = np.diag([np.sqrt(1 - t*t), t]).astype(complex)
    tangent = np.diag([-t / np.sqrt(1 - t*t), 1]).astype(complex)
    if rotate:
        psi = np.array([[np.cos(t / 2)], [np.sin(t / 2)]])
        dot = np.array([[-np.sin(t / 2) / 2], [np.cos(t / 2) / 2]])
        tangent, columns = (np.kron(dot, columns) + np.kron(psi, tangent),
                            np.kron(psi, columns))
    return columns, tangent


def state_tangent(columns, tangent):
    return (columns @ columns.conj().T,
            tangent @ columns.conj().T + columns @ tangent.conj().T)


def spectral_sld_qfi(rho, dot):
    values, vectors = np.linalg.eigh(rho)
    matrix = vectors.conj().T @ dot @ vectors
    denominator = values[:, None] + values[None, :]
    keep = denominator > 1e-12
    return float(np.sum(2 * np.abs(matrix[keep])**2 / denominator[keep]))


@pytest.mark.parametrize('rotate,pointwise,limit', [(False, 0., 4.), (True, 1., 5.)])
def test_pointwise_sld_differs_from_smooth_lift(rotate, pointwise, limit):
    columns, tangent = preparation(0., rotate)
    rho, dot = state_tangent(columns, tangent)
    np.testing.assert_allclose(np.trace(rho), 1)
    np.testing.assert_allclose(spectral_sld_qfi(rho, dot), pointwise)
    np.testing.assert_allclose(4 * np.linalg.norm(tangent, 'fro')**2, limit)
    np.testing.assert_allclose(np.vdot(columns, tangent), 0, atol=1e-14)
    for t in (-0.2, -0.01, 0.01, 0.2):
        c, d = preparation(t, rotate)
        expected = 4 / (1 - t*t) + int(rotate)
        np.testing.assert_allclose(spectral_sld_qfi(*state_tangent(c, d)), expected)
        np.testing.assert_allclose(4 * np.linalg.norm(d, 'fro')**2, expected)


@pytest.mark.parametrize('rotate', [False, True])
def test_every_finite_hermitian_gauge_has_unremovable_kernel_tangent(rotate):
    columns, tangent = preparation(0., rotate)
    rng = np.random.default_rng(299)
    for _ in range(20):
        matrix = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        h = (matrix + matrix.conj().T) / 2
        # Kraus rows transform as dK_i -> dK_i - i sum_j h_ij K_j.
        gauged = tangent - 1j * columns @ h.T
        alpha = np.linalg.norm(gauged, 'fro')**2
        overlap = np.vdot(columns, gauged)
        base = 1 + int(rotate) / 4
        np.testing.assert_allclose(alpha, base + h[0, 0].real**2 + abs(h[1, 0])**2)
        np.testing.assert_allclose(alpha - abs(overlap)**2, base + abs(h[1, 0])**2)


@pytest.mark.parametrize('n', [1, 2, 3])
def test_positive_parallel_qfi_does_not_remove_rank_change(n):
    rho, dot = state_tangent(*preparation(0., rotate=True))
    total, derivative = np.ones((1, 1)), np.zeros((1, 1))
    for _ in range(n):
        derivative = np.kron(derivative, rho) + np.kron(total, dot)
        total = np.kron(total, rho)
    np.testing.assert_allclose(spectral_sld_qfi(total, derivative), n, atol=1e-12)
