from pathlib import Path
import numpy as np
from complete_qcqc_circuit import completion, apply_channel, load_operations


def test_cptp_on_arbitrary_inputs_and_coherent_support():
    folder = Path(__file__).parent/'qfi_results'
    ops, _ = load_operations(folder, 3)
    rng = np.random.default_rng(811)
    for op in ops.values():
        v, null, diagnostic = completion(op['map'])
        assert diagnostic['tp_residual'] < 1e-12
        assert diagnostic['polar_change'] < 1e-10
        dim = v.shape[1]
        a = rng.normal(size=(dim, 3))+1j*rng.normal(size=(dim, 3))
        rho = a@a.conj().T
        rho /= np.trace(rho)
        output = apply_channel(v, null, rho)
        assert abs(np.trace(output)-1) < 1e-12
        assert np.linalg.eigvalsh(output)[0] > -1e-12
        support = v.conj().T@v
        supported = support@rho@support
        assert np.allclose(apply_channel(v, null, supported), v@supported@v.conj().T, atol=1e-12)


def test_complement_is_required():
    u = np.array([[1, 0]], complex)
    v, null, _ = completion(u)
    rho = np.diag([0, 1])
    assert np.trace(v@rho@v.conj().T) == 0
    assert np.trace(apply_channel(v, null, rho)) == 1
