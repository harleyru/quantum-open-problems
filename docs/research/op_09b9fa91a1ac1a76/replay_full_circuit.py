"""Density-matrix replay of all completion Kraus branches, without pruning."""
from itertools import product
import numpy as np
from qcqc_feasibility import subsets
from complete_qcqc_circuit import load_operations


def replay(ops, nulls, n, theta=0.0):
    phase = np.diag([np.exp(-0.5j*theta), np.exp(0.5j*theta)])
    ks = [np.diag([1, np.sqrt(0.5)])@phase,
          np.array([[0, np.sqrt(0.5)], [0, 0]])@phase]
    ds = [k@(-0.5j*np.diag([1, -1])) for k in ks]
    rho_total = np.zeros((ops[frozenset(range(n))]['map'].shape[0],)*2, complex)
    drho_total = np.zeros_like(rho_total)
    for labels in product(range(2), repeat=n):
        initial = ops[frozenset()]['map']
        rho = initial@initial.conj().T
        drho = np.zeros_like(rho)
        layout, offset = {}, 0
        for k, width in ops[frozenset()]['outgoing']:
            layout[frozenset(), k] = (offset, width)
            offset += width
        for step in range(1, n+1):
            used_sets = [k for k in subsets(n) if len(k) == step]
            ins, outs, ni, no = {}, {}, 0, 0
            for used in used_sets:
                u = ops[used]['map']
                ins[used], outs[used] = ni, no
                ni += u.shape[1]; no += u.shape[0]
            external = np.zeros((ni, len(rho)), complex)
            dexternal = np.zeros_like(external)
            control = np.zeros((no, ni), complex)
            next_layout = {}
            for used in used_sets:
                op = ops[used]
                start, out = ins[used], outs[used]
                control[out:out+op['map'].shape[0], start:start+op['map'].shape[1]] = op['map']
                current = start
                for k, width in op['incoming']:
                    old, old_width = layout[used-{k}, k]
                    assert old_width == width
                    external[current:current+width, old:old+width] = np.kron(ks[labels[k]], np.eye(width//2))
                    dexternal[current:current+width, old:old+width] = np.kron(ds[labels[k]], np.eye(width//2))
                    current += width
                current = out
                for k, width in op['outgoing']:
                    next_layout[used, k] = (current, width)
                    current += width
            incoming = external@rho@external.conj().T
            dincoming = (dexternal@rho@external.conj().T+external@drho@external.conj().T
                         +external@rho@dexternal.conj().T)
            rho = control@incoming@control.conj().T
            drho = control@dincoming@control.conj().T
            for used in used_sets:
                start, out = ins[used], outs[used]
                null = nulls[used]
                width = null.shape[1]
                block = slice(start, start+width)
                rho[out, out] += np.trace(null@incoming[block, block]@null.conj().T)
                drho[out, out] += np.trace(null@dincoming[block, block]@null.conj().T)
            layout = next_layout
        rho_total += rho
        drho_total += drho
    return rho_total, drho_total


if __name__ == '__main__':
    import json
    from pathlib import Path
    from extract_qcqc_circuit import simulate
    folder = Path(__file__).parent/'qfi_results'
    ops, nulls = load_operations(folder, 3, completed=True)
    rho, drho = replay(ops, nulls, 3)
    principal, derivative = simulate(ops, 3)
    p, q = np.linalg.eigh((rho+rho.conj().T)/2)
    denom = p[:, None]+p[None, :]
    mask = denom > 1e-10
    d = q.conj().T@drho@q
    result = {'n': 3, 'all_completion_branches': True,
              'output_trace': float(np.trace(rho).real),
              'output_min_eigenvalue': float(p[0]),
              'principal_difference': float(np.max(np.abs(rho-principal))),
              'derivative_difference': float(np.max(np.abs(drho-derivative))),
              'qfi': float(np.sum(2*np.abs(d[mask])**2/denom[mask])),
              'verification_level': 'floating-point full-branch replay, not interval error bound'}
    (folder/'full_replay_n3.json').write_text(json.dumps(result, indent=2)+'\n', encoding='ascii')
    print(json.dumps(result, indent=2))
