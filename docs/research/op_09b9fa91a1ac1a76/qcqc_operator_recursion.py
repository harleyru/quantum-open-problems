"""Finite operator-valued QC-QC path model.

The route label is retained in the terminal output. This is a valid, explicit
coherent superposition of fixed orders (a QC-QC subclass). Stacking the route
blocks within each Kraus operator retains off-diagonal control coherences in
the output state. Dephasing or discarding that control removes them. This is
not a dynamically routed QC-QC implementation or a process-matrix SDP solver.
"""
from itertools import permutations, product
import numpy as np


def route_kraus(kraus, n, controls=None):
    """Return route-labelled global Kraus operators and their derivatives."""
    routes = list(permutations(range(n)))
    controls = controls or {}
    rank, dim = len(kraus), kraus[0][0].shape[0]
    out = []
    for indices in product(range(rank), repeat=n):
        blocks, dblocks = [], []
        for route in routes:
            value = np.eye(dim, dtype=complex)
            derivative = np.zeros_like(value)
            used = frozenset()
            for slot in route:
                internal = controls.get((used, slot), np.eye(dim, dtype=complex))
                value = internal @ value
                derivative = internal @ derivative
                operator, dot_operator = kraus[indices[slot]]
                derivative = dot_operator @ value + operator @ derivative
                value = operator @ value
                used = used | {slot}
            blocks.append(value / np.sqrt(len(routes)))
            dblocks.append(derivative / np.sqrt(len(routes)))
        out.append((np.vstack(blocks), np.vstack(dblocks)))
    return out


def completeness(operators):
    return sum(op.conj().T @ op for op, _ in operators)


def derivative_alpha(operators):
    return sum(derivative.conj().T @ derivative
               for _, derivative in operators)
