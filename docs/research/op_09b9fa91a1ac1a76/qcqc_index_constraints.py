"""Combinatorial audit of the subset-indexed QC-QC recursion.

This checks the routing-index structure in WLP26 Theorem A. It deliberately
does not claim to solve the operator-valued process-matrix SDP.
"""
from itertools import combinations, permutations


def subsets(n, size):
    return {frozenset(x) for x in combinations(range(n), size)}


def audit(n):
    all_items = frozenset(range(n))
    for size in range(n):
        for used in subsets(n, size):
            # Every unused call can be selected next.
            next_states = {used | {k} for k in all_items - used}
            assert len(next_states) == n - size
            # The outgoing equality's right side has one term per possible
            # last call, and these are exactly the immediate predecessors.
            predecessors = {used - {k} for k in used}
            assert len(predecessors) == size
            if size:
                assert all(pred < used for pred in predecessors)
    # Complete routes are ordered call sequences, not unordered subsets.
    routes = list(permutations(range(n)))
    assert len(routes) == _factorial(n)
    for route in routes:
        prefixes = [frozenset(route[:size]) for size in range(n + 1)]
        assert len(prefixes) == n + 1
        assert len(set(prefixes)) == n + 1
        assert prefixes[0] == frozenset()
        assert prefixes[-1] == all_items
        assert all(prefixes[i] < prefixes[i + 1] for i in range(n))
    return {
        'N': n,
        'subset_nodes': 2**n,
        'complete_routes': len(routes),
        'complete_subset': all_items,
    }


def _factorial(n):
    result = 1
    for value in range(2, n + 1):
        result *= value
    return result


if __name__ == '__main__':
    for n in range(2, 8):
        print(audit(n))
