"""Numerical effective-support circuit extraction from repaired QC-QC nodes."""
import numpy as np
from qcqc_feasibility import QCQCFeasibility, subsets
from repair_qcqc_primal import repair


def root(a):
    vals, vecs = np.linalg.eigh((a+a.conj().T)/2)
    if vals[0] <= 0:
        raise ValueError('This extractor requires positive definite repaired witnesses')
    return (vecs*np.sqrt(vals))@vecs.conj().T


def incoming_factor(model, used, k, factor):
    axes = model.axes(used-{k}, k)
    rank = factor.shape[1]
    tensor = np.einsum('am,ob->aomb', factor, np.eye(2))
    tensor = tensor.reshape(*[model.dim[a] for a in axes], 2, rank, 2)
    source = axes+(f'O{k}', 'memory', 'output')
    target = model.axes(used)+('output', 'memory')
    return tensor.transpose([source.index(a) for a in target]).reshape(model.size(model.axes(used)), -1)


def outgoing_factor(model, used, k, factor):
    axes = model.axes(used, k)
    source = axes+('memory',)
    target = model.axes(used)+(f'I{k}', 'memory')
    return factor.reshape(*[model.dim[a] for a in axes], factor.shape[1]).transpose(
        [source.index(a) for a in target]).reshape(model.size(model.axes(used)), -1)


def extract(process, nodes, n):
    model = QCQCFeasibility([2]*n, [2]*n)
    factors = {key: root(a) for key, a in nodes.items()}
    operations, diagnostics = {}, {}
    full = frozenset(range(n))
    for used in subsets(n):
        incoming = [(k, incoming_factor(model, used, k, factors[used-{k}, k]))
                    for k in sorted(used)]
        outgoing = [(k, outgoing_factor(model, used, k, factors[used, k]))
                    for k in range(n) if k not in used]
        a = np.concatenate([v for _, v in incoming], axis=1) if incoming else np.ones((1, 1))
        b = np.concatenate([v for _, v in outgoing], axis=1) if used != full else root(process)
        inverse = np.linalg.pinv(a.T, rcond=1e-12)
        u = b.T@inverse
        support = a.T@inverse
        operations[used] = {'map': u, 'incoming': [(k, v.shape[1]) for k, v in incoming],
                            'outgoing': [(k, v.shape[1]) for k, v in outgoing]}
        diagnostics[str(sorted(used))] = {
            'input_dimension': a.shape[1], 'output_dimension': b.shape[1],
            'effective_rank': int(np.linalg.matrix_rank(a)),
            'link_residual': float(np.max(np.abs(u@a.T-b.T))),
            'isometry_residual': float(np.max(np.abs(u.conj().T@u-support))),
        }
    return operations, diagnostics


def simulate(operations, n, theta=0.0, channels=None):
    """Insert external Kraus maps, coherently recombine paths, retain slot environments."""
    from itertools import product
    phase = np.diag([np.exp(-0.5j*theta), np.exp(0.5j*theta)])
    ks = [np.diag([1, np.sqrt(0.5)])@phase,
          np.array([[0, np.sqrt(0.5)], [0, 0]])@phase]
    ds = [k@(-0.5j*np.diag([1, -1])) for k in ks]
    if channels is None:
        channels = [(ks, ds) for _ in range(n)]
    columns, derivatives = [], []
    full = frozenset(range(n))
    for labels in product(*[range(len(ks)) for ks, _ in channels]):
        states, tangents = {}, {}
        for used in subsets(n):
            op = operations[used]
            if not used:
                v, dv = np.ones(1), np.zeros(1)
            else:
                chunks, dchunks = [], []
                for k, width in op['incoming']:
                    ks, ds = channels[k]
                    before = states[used-{k}, k].reshape(2, -1)
                    tangent = tangents[used-{k}, k].reshape(2, -1)
                    chunks.append((ks[labels[k]]@before).ravel())
                    dchunks.append((ds[labels[k]]@before+ks[labels[k]]@tangent).ravel())
                v, dv = np.concatenate(chunks), np.concatenate(dchunks)
            out, dout = op['map']@v, op['map']@dv
            if used == full:
                columns.append(out); derivatives.append(dout)
            else:
                start = 0
                for k, width in op['outgoing']:
                    states[used, k] = out[start:start+width]
                    tangents[used, k] = dout[start:start+width]
                    start += width
    v, dv = np.column_stack(columns), np.column_stack(derivatives)
    return v@v.conj().T, dv@v.conj().T+v@dv.conj().T


def load_repaired(path, n):
    process, parts = repair(path, n)
    s = np.asarray(process[0], float)+1j*np.asarray(process[1], float)
    nodes = {key: np.asarray(parts[0][key], float)+1j*np.asarray(parts[1][key], float)
             for key in parts[0]}
    return s, nodes


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path
    from qcqc_qfi import ad_ensemble
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, choices=[2, 3], required=True)
    args = parser.parse_args()
    n = args.n
    folder = Path(__file__).parent/'qfi_results'
    s, nodes = load_repaired(folder/f'qfi_n{n}.npz', n)
    ops, diagnostics = extract(s, nodes, n)
    rho, drho = simulate(ops, n)
    c, dc = ad_ensemble(n)
    v, dv = root(s).T@c, root(s).T@dc
    expected, dexpected = v@v.conj().T, dv@v.conj().T+v@dv.conj().T
    p, q = np.linalg.eigh(rho)
    d = q.conj().T@drho@q
    denom = p[:, None]+p[None, :]
    mask = denom > 1e-10
    result = {'n': n, 'operations': diagnostics,
              'output_residual': float(np.max(np.abs(rho-expected))),
              'derivative_residual': float(np.max(np.abs(drho-dexpected))),
              'output_trace': float(np.trace(rho).real),
              'qfi': float(np.sum(2*np.abs(d[mask])**2/denom[mask])),
              'verification_level': 'numerical effective-support realization'}
    np.savez_compressed(folder/f'circuit_n{n}.npz',
                        **{f'map_{sum(1 << k for k in used)}': op['map'] for used, op in ops.items()})
    result['layout'] = {str(sum(1 << k for k in used)): {'incoming': op['incoming'], 'outgoing': op['outgoing']}
                        for used, op in ops.items()}
    (folder/f'circuit_n{n}.json').write_text(json.dumps(result, indent=2)+'\n', encoding='ascii')
    print(json.dumps(result, indent=2))
