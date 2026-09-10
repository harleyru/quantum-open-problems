"""Nontrivial complex Stinespring tangents exercise the WLP beta convention."""
import numpy as np
import pytest


@pytest.mark.parametrize('input_dim,output_dim,rank', [(2, 2, 2), (2, 3, 2), (3, 2, 3)])
def test_residual_gram_and_beta_convention(input_dim, output_dim, rank):
    rng = np.random.default_rng(1701)
    total = output_dim * rank
    raw = rng.normal(size=(total, input_dim)) + 1j * rng.normal(size=(total, input_dim))
    v, _ = np.linalg.qr(raw)
    raw_h = rng.normal(size=(total, total)) + 1j * rng.normal(size=(total, total))
    h = (raw_h + raw_h.conj().T) / 2
    dv = -1j * h @ v
    beta = dv.conj().T @ v
    alpha = dv.conj().T @ dv
    residual = dv - v @ beta.conj().T
    np.testing.assert_allclose(beta + beta.conj().T, 0, atol=1e-12)
    np.testing.assert_allclose(v.conj().T @ residual, 0, atol=1e-12)
    np.testing.assert_allclose(residual.conj().T @ residual,
                               alpha - beta @ beta.conj().T, atol=1e-12)
    assert np.linalg.eigvalsh(alpha - beta @ beta.conj().T).min() >= -1e-12
    wrong_beta = v.conj().T @ dv
    wrong_residual = dv - v @ wrong_beta.conj().T
    assert np.linalg.norm(v.conj().T @ wrong_residual) > 1e-3
