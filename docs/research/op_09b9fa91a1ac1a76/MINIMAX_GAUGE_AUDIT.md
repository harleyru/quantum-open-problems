# Minimax and Gauge Audit

## Minimax direction

Let `f(S,h)` be the QFI performance functional, with `S` ranging over valid
processes and `h` over global environment gauges. The elementary minimax
inequality gives

```text
max_S min_h f(S,h) <= min_h max_S f(S,h).
```

Restricting the right-hand minimization to product gauges induced by one
single-channel Kraus representation can only increase its value. Therefore
the product-gauge expression is a valid upper bound, not an equality. This
confirms the direction used in WLP `eq:final_derivation_1` and
`eq:final_derivation_2`.

## Positive-semidefinite residual

For a Kraus column `V` and its derivative `dot(V)`, define

```text
beta = dot(V)^dagger V,
R = dot(V) - V beta^dagger.
```

Trace preservation gives `V^dagger V = I` and
`beta + beta^dagger = 0`. Direct multiplication yields

```text
V^dagger R = 0,
R^dagger R = dot(V)^dagger dot(V) - beta beta^dagger.
```

The last operator is a Gram matrix, hence positive semidefinite. In Kraus
components this is exactly
`alpha - beta beta^dagger >= 0`. The derivation is representation-wise and
does not compare alpha and beta from different gauges.

Explicitly, `V^dagger dot(V) = beta^dagger`, so
`V^dagger R = beta^dagger - beta^dagger = 0`. Expanding all four terms
of `R^dagger R` subtracts two copies of `beta beta^dagger` and adds one.
Equivalently, `R = (I - V V^dagger) dot(V)`. An earlier version used
`beta = V^dagger dot(V)` while retaining the opposite-convention residual;
that sign/conjugation error is corrected here. `test_stinespring_residual.py`
contains a negative-control test that rejects the old convention.

## Remaining qualification

This audit validates the minimax direction and residual positivity. It does
not by itself prove the performance-operator expansion, the reduced-process
normalization identities, or parallel attainability. Those are tracked in
`THEOREM_2_3_AUDIT.md` and `PARALLEL_ATTAINABILITY.md`.
