"""Dense trace-and-replace implementation with interleaved input/output axes."""
import numpy as np


def depolarize(matrix, dims, axis):
    """Trace one tensor factor and replace it by its normalized identity."""
    count = len(dims)
    tensor = matrix.reshape(tuple(dims) * 2)
    reduced = np.trace(tensor, axis1=axis, axis2=axis + count)
    remaining = [k for k in range(2 * count) if k not in (axis, axis + count)]
    expanded = np.tensordot(reduced, np.eye(dims[axis]) / dims[axis], axes=0)
    order = remaining + [axis, axis + count]
    return expanded.transpose(np.argsort(order)).reshape(matrix.shape)


def q_gen(matrix, dims):
    """Apply WLP eq:qico_appendix without constructing a superoperator matrix."""
    if not dims or len(dims) % 2 or any(d < 1 for d in dims):
        raise ValueError('dims must contain positive input/output dimension pairs')
    total = int(np.prod(dims))
    matrix = np.asarray(matrix, dtype=complex)
    if matrix.shape != (total, total):
        raise ValueError('matrix shape does not match dims')
    result = matrix.copy()
    for input_axis in range(0, len(dims), 2):
        output_replaced = depolarize(result, dims, input_axis + 1)
        result = result - output_replaced + depolarize(output_replaced, dims, input_axis)
    return result - np.trace(matrix) * np.eye(total) / total
