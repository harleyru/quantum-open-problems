import numpy as np
import pytest
from itertools import product
from ico_projector import depolarize, q_gen


def random_matrix(rng, dimension):
    return rng.normal(size=(dimension, dimension)) + 1j * rng.normal(size=(dimension, dimension))


def test_two_qubit_slots_on_complete_pauli_basis():
    paulis = [np.eye(2), np.array([[0, 1], [1, 0]]),
              np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
    image_dimension = 0
    for indices in product(range(4), repeat=4):
        basis = np.array([[1]])
        for index in indices:
            basis = np.kron(basis, paulis[index])
        # From normalization: each slot must be identity or output-traceless.
        selected = any(indices) and all(
            indices[k + 1] != 0 or indices[k] == 0 for k in (0, 2))
        image_dimension += int(selected)
        np.testing.assert_allclose(q_gen(basis, (2, 2, 2, 2)),
                                   int(selected) * basis, atol=1e-12)
    assert image_dimension == (16 - 4 + 1)**2 - 1


@pytest.mark.parametrize('dims', [(2, 2), (2, 2, 2, 2), (2, 3, 3, 2), (2, 2, 2, 2, 2, 2)])
def test_self_adjoint_projector(dims):
    rng = np.random.default_rng(20260909)
    size = int(np.prod(dims))
    x, y = random_matrix(rng, size), random_matrix(rng, size)
    qx, qy = q_gen(x, dims), q_gen(y, dims)
    np.testing.assert_allclose(q_gen(qx, dims), qx, atol=1e-12)
    np.testing.assert_allclose(np.vdot(x, qy), np.vdot(qx, y), atol=1e-10)
    np.testing.assert_allclose(q_gen(np.eye(size), dims), 0, atol=1e-12)


def test_random_image_kernel_orthogonality():
    rng = np.random.default_rng(91)
    dims = (2, 3, 2, 2)
    size = int(np.prod(dims))
    for _ in range(8):
        x = random_matrix(rng, size)
        y = random_matrix(rng, size)
        image, kernel = q_gen(x, dims), y - q_gen(y, dims)
        np.testing.assert_allclose(np.vdot(image, kernel), 0, atol=1e-10)


@pytest.mark.parametrize('dims', [(2, 2, 2, 2), (2, 3, 3, 2), (1, 2, 3, 2)])
def test_two_party_perpendicular_cancellation(dims):
    rng = np.random.default_rng(42)
    factors = []
    for pair in (dims[:2], dims[2:]):
        raw = random_matrix(rng, int(np.prod(pair)))
        factors.append(raw - depolarize(raw, pair, 1))
    x = np.kron(factors[0], factors[1].conj().T)
    np.testing.assert_allclose(q_gen(x, dims), x, atol=1e-12)
    size = int(np.prod(dims))
    raw = random_matrix(rng, size)
    hermitian = (raw + raw.conj().T) / 2
    kernel = hermitian - q_gen(hermitian, dims)
    kernel -= np.trace(kernel) * np.eye(size) / size
    # Kernel perturbations of white noise give positive normalized processes.
    process = np.eye(size) + kernel / (2 * max(1, np.linalg.norm(kernel, 2)))
    process *= (dims[1] * dims[3]) / np.trace(process)
    assert np.linalg.eigvalsh(process).min() > 0
    np.testing.assert_allclose(q_gen(process, dims), 0, atol=1e-12)
    np.testing.assert_allclose(np.trace(x @ process), 0, atol=1e-10)
