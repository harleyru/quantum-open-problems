import numpy as np
import pytest
from qcqc_operator_recursion import completeness, derivative_alpha, route_kraus


def amplitude_damping(theta, damping=0.4):
    c, s = np.cos(theta), np.sin(theta)
    # Rotation followed by amplitude damping, with two nonzero Kraus operators.
    u = np.array([[c, -s], [s, c]], dtype=complex)
    du = np.array([[-s, -c], [c, -s]], dtype=complex)
    noise = [np.diag([1, np.sqrt(1 - damping)]),
             np.array([[0, np.sqrt(damping)], [0, 0]])]
    return [(a @ u, a @ du) for a in noise]


@pytest.mark.parametrize('n', [1, 2, 3])
def test_route_labelled_operator_recursion_is_trace_preserving(n):
    operators = route_kraus(amplitude_damping(0.37), n)
    np.testing.assert_allclose(completeness(operators), np.eye(2), atol=1e-12)
    alpha = derivative_alpha(operators)
    assert np.all(np.linalg.eigvalsh(alpha) >= -1e-12)


def test_route_count_and_operator_shapes():
    operators = route_kraus(amplitude_damping(0.37), 3)
    assert len(operators) == 2**3
    for operator, derivative in operators:
        assert operator.shape == (6 * 2, 2)
        assert derivative.shape == operator.shape


def test_product_derivative_matches_finite_difference():
    theta, eps = 0.37, 1e-6
    center = route_kraus(amplitude_damping(theta), 3)
    plus = route_kraus(amplitude_damping(theta + eps), 3)
    minus = route_kraus(amplitude_damping(theta - eps), 3)
    for (_, derivative), (p, _), (m, _) in zip(center, plus, minus):
        np.testing.assert_allclose(derivative, (p - m) / (2 * eps),
                                   atol=2e-9, rtol=2e-9)


def test_retained_control_has_measurable_coherence():
    operators = route_kraus(amplitude_damping(0.37), 2)
    rho = np.ones((2, 2)) / 2
    output = sum(k @ rho @ k.conj().T for k, _ in operators)
    dephased = output.copy()
    dephased[:2, 2:] = 0
    dephased[2:, :2] = 0
    assert np.linalg.norm(output[:2, 2:]) > 1e-3
    plus = np.ones((2, 2)) / 2
    plus_effect = np.kron(plus, np.eye(2))
    coherent_probability = np.trace(plus_effect @ output).real
    classical_probability = np.trace(plus_effect @ dephased).real
    np.testing.assert_allclose(classical_probability, 0.5, atol=1e-12)
    assert abs(coherent_probability - classical_probability) > 1e-3
    # Discarding the route register produces the same target marginal.
    target = output[:2, :2] + output[2:, 2:]
    target_dephased = dephased[:2, :2] + dephased[2:, 2:]
    np.testing.assert_allclose(target, target_dephased, atol=1e-12)
    np.testing.assert_allclose(np.trace(output), 1, atol=1e-12)


def test_noisy_output_derivative_matches_finite_difference():
    rho = np.ones((2, 2)) / 2
    theta, eps = 0.37, 1e-6
    operators = route_kraus(amplitude_damping(theta), 3)
    derivative = sum(d @ rho @ k.conj().T + k @ rho @ d.conj().T
                     for k, d in operators)
    def state(t):
        return sum(k @ rho @ k.conj().T
                   for k, _ in route_kraus(amplitude_damping(t), 3))
    np.testing.assert_allclose(derivative, (state(theta + eps) - state(theta - eps)) / (2 * eps),
                               atol=2e-9, rtol=2e-9)
    np.testing.assert_allclose(np.trace(derivative), 0, atol=1e-12)


def test_history_dependent_unitary_controls_preserve_normalization():
    theta = 0.37
    z = np.diag([1, -1]).astype(complex)
    controls = {
        (frozenset(), 0): z,
        (frozenset({0}), 1): z,
        (frozenset(), 1): z,
        (frozenset({1}), 0): z,
    }
    operators = route_kraus(amplitude_damping(theta), 2, controls)
    np.testing.assert_allclose(completeness(operators), np.eye(2), atol=1e-12)
    rho = np.eye(2) / 2
    state = sum(k @ rho @ k.conj().T for k, _ in operators)
    np.testing.assert_allclose(np.trace(state), 1, atol=1e-12)
