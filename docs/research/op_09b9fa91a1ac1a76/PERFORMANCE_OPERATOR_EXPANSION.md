# Product-Kraus Performance-Operator Expansion

For one channel write the Kraus column as `K` and its derivative as `dot K`.
For N independent uses, the product Kraus operator is
`K_i = K_{i_1} tensor ... tensor K_{i_N}`. Leibniz gives

```text
dot(K_i) = sum_j K_{i_1} tensor ... tensor dot(K_{i_j}) tensor ... tensor K_{i_N}.
```

Expanding `4 sum_i |dot(K_i)><dot(K_i)|` produces N diagonal terms and
`N(N-1)/2` unordered pairs, each appearing with its Hermitian conjugate.
The diagonal term at slot j is

```text
(E^T)^(tensor(j-1)) tensor Omega^(1) tensor (E^T)^(tensor(N-j)),
```

while a pair j<k is

```text
1/4 [E^T ... tensor Lambda_j tensor ... tensor Lambda_k^dagger ... E^T
     + h.c.].
```

The factor `1/4` follows from WLP's definition
`Lambda = 4 sum_i |dot(K_i^dagger)><K_i^dagger|`; the two cross terms in the
Leibniz expansion are precisely the displayed term and its adjoint. This
checks the structure and pair counting in WLP `eq:Omega_expansion_nonlocal`.

The derivation assumes a fixed product Kraus gauge for all N slots. It does
not assert that minimizing over product gauges equals minimizing over all
global gauges; the minimax audit records only the valid upper-bound direction.
Transpose placement depends on the Choi/vectorization convention and must be
kept consistent between `E`, `Lambda`, and the process link product.
