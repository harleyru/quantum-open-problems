"""Full-input CPTP completion and exported-circuit replay (floating point)."""
import json
from pathlib import Path
import numpy as np
from extract_qcqc_circuit import simulate


def completion(u, cutoff=0.5):
    left, singular, right = np.linalg.svd(u, full_matrices=True)
    rank = int(np.sum(singular > cutoff))
    principal = left[:, :rank]@right[:rank, :]
    null = right[rank:, :]
    # K_0=principal, K_j=|0_out><null_j|. The latter reset unreachable inputs.
    residual = principal.conj().T@principal+null.conj().T@null-np.eye(u.shape[1])
    return principal, null, {'rank': rank,
        'tp_residual': float(np.max(np.abs(residual))),
        'polar_change': float(np.max(np.abs(principal-u))),
        'singular_distance_to_partial_isometry': float(max(
            np.max(np.abs(singular[:rank]-1), initial=0),
            np.max(singular[rank:], initial=0)))}


def apply_channel(principal, null, rho):
    output = principal@rho@principal.conj().T
    output[0, 0] += np.trace(null@rho@null.conj().T)
    return output


def load_operations(folder, n, completed=False):
    folder = Path(folder)
    layout = json.loads((folder/f'circuit_n{n}.json').read_text())['layout']
    prefix = 'completed_circuit' if completed else 'circuit'
    operations, complements = {}, {}
    with np.load(folder/f'{prefix}_n{n}.npz') as data:
        for mask, blocks in layout.items():
            used = frozenset(k for k in range(n) if int(mask) & (1 << k))
            operations[used] = {'map': data[f'map_{mask}'], **blocks}
            if completed:
                complements[used] = data[f'null_{mask}']
    return operations, complements


def export_and_replay(folder, n):
    folder = Path(folder)
    original, _ = load_operations(folder, n)
    arrays, diagnostics = {}, {}
    for used, op in original.items():
        mask = sum(1 << k for k in used)
        principal, null, diagnostic = completion(op['map'])
        arrays[f'map_{mask}'], arrays[f'null_{mask}'] = principal, null
        diagnostics[str(mask)] = diagnostic
    np.savez_compressed(folder/f'completed_circuit_n{n}.npz', **arrays)
    loaded, complements = load_operations(folder, n, completed=True)
    rho, drho = simulate(loaded, n)
    expected, dexpected = simulate(original, n)
    p, q = np.linalg.eigh(rho)
    denominator = p[:, None]+p[None, :]
    mask = denominator > 1e-10
    derivative = q.conj().T@drho@q
    result = {'n': n, 'maps': diagnostics,
              'replayed_principal_output_difference': float(np.max(np.abs(rho-expected))),
              'replayed_principal_derivative_difference': float(np.max(np.abs(drho-dexpected))),
              'principal_output_trace': float(np.trace(rho).real),
              'principal_output_qfi': float(np.sum(2*np.abs(derivative[mask])**2/denominator[mask])),
              'verification_level': 'numerical CPTP completion; principal-path replay'}
    (folder/f'completed_circuit_n{n}.json').write_text(json.dumps(result, indent=2)+'\n', encoding='ascii')
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, choices=[2, 3], required=True)
    args = parser.parse_args()
    print(json.dumps(export_and_replay(Path(__file__).parent/'qfi_results', args.n), indent=2))
