import math

import pytest

from qcqc_index_constraints import audit


@pytest.mark.parametrize('n', range(2, 8))
def test_subset_indexed_routes(n):
    result = audit(n)
    assert result['subset_nodes'] == 2**n
    assert result['complete_routes'] == math.factorial(n)

