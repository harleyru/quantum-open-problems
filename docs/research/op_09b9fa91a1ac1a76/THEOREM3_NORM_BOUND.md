# Theorem 3 Norm-Bound Check

Let `V` be the Kraus column and decompose its tangent as
`dot V = V beta^dagger + R`, with `V^dagger R = 0`. For a normalized positive
process `S`, define the induced seminorm
`||X||_S^2 = tr(X^dagger S^T X)`. The Gram identities give

```text
||P_j||_S^2 <= ||beta||^2,
||R_j||_S^2 <= ||alpha - beta beta^dagger||,
```

where `P_j` is the parallel component in slot j. The process-kernel
orthogonality removes cross terms between different residual slots, while
the triangle inequality gives

```text
||sum_j P_j||_S <= N ||beta||,
||sum_j R_j||_S <= sqrt(N ||alpha-beta beta^dagger||).
```

Consequently

```text
F_Gen(N) <= 4 ||sum_j(P_j+R_j)||_S^2
          <= 4 [N ||beta|| + sqrt(N ||alpha-beta beta^dagger||)]^2.
```

This independently checks the scaling and coefficient arithmetic in WLP
`eq:bound3`/`eq:ultimate_bound_supp`. The remaining qualifications are that
the process-kernel orthogonality must use the deterministic-process lemma in
`ICO_CANCELLATION_PROOF.md`, and that optimization over gauges is a minimax
upper-bound step rather than an equality without a saddle-point theorem.
