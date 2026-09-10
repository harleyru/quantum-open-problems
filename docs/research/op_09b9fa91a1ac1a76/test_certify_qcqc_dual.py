from pathlib import Path
import numpy as np
from mpmath import iv
from certify_qcqc_dual import certify, hermitian_decimal, positive_ldl, exact_ensemble
from qcqc_qfi import ad_ensemble
from qcqc_support_dual import solve_joint


def test_interval_ldl_accepts_and_rejects():
    iv.dps = 50
    assert positive_ldl(hermitian_decimal(np.array([[2, 1j], [-1j, 2]])))[0]
    assert not positive_ldl(hermitian_decimal(np.array([[1, 2j], [-2j, 1]])))[0]


def test_exact_channel_ordering():
    c, dc = exact_ensemble(2)
    expected = ad_ensemble(2)
    for a, b in zip((c, dc), expected):
        numeric = np.array([[complex(float(z.real.mid), float(z.imag.mid))
                             for z in row] for row in a])
        assert np.allclose(numeric, b)


def test_joint_one_query():
    result = solve_joint(1)
    assert abs(result['objective']-(12-8*np.sqrt(2))) < 1e-6


def test_saved_certificate_and_corruption(tmp_path):
    source = Path(__file__).parent/'qfi_results/joint_dual_n3.npz'
    assert certify(source, 3)['verified']
    with np.load(source) as data:
        arrays = {key: data[key] for key in data.files}
    arrays['y_0'] = arrays['y_0']-1
    bad = tmp_path/'bad.npz'
    np.savez(bad, **arrays)
    assert not certify(bad, 3)['verified']
