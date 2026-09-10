# Fixed-gauge support dual

Let Omega be a fixed Hermitian performance operator and maximize Tr(Omega S)
over reduced QC-QC processes. For every subset K introduce Hermitian Y_K
on the slot input/output spaces of K; set Y_all=Omega. For every k outside K
require

    Z_K,k = Y_K tensor I_I_k - Tr_O_k Y_(K union {k}) >= 0,

with tensor axes permuted to those of X[K,k]. Minimize Tr(Y_empty).
The empty space is scalar when P=1. No positivity constraint on Y_K itself
is needed. Summing the primal flow equations against Y_K telescopes to

    Tr(Y_empty) - Tr(Omega S) = sum_(K,k) Tr(Z_K,k X[K,k]) >= 0.

This proves weak duality directly, without strong duality assumptions.
For any fixed Hermitian Kraus gauge h, the optimized QFI is at most the
support function for Omega(h), so a feasible Y supplies an upper bound.

`qcqc_support_dual.py` constructs numerical candidates only. Tests check the
adjoint telescoping identity on complex matrices with unequal dimensions,
and verify the identity-operator support against the known process trace.

## Observed gap

At N=3, u=0.5, using the gauge reconstructed from the primal process gives
a support-dual objective about 5.21542160. The primal QFI objective was
4.73731164. The smallest dual slack eigenvalue is about -4.99e-8, so even
the larger number is not yet a rigorous upper bound as stored.

A minimizer of f(S_star,h) need not minimize max_S f(S,h). Degenerate
gauge directions at S_star can change performance on other processes.
This experiment therefore does not refute the primal optimization, but
does rule out treating its reconstructed gauge as an established saddle
point. A joint gauge/dual optimization is the next step toward tight bounds.

For exact certification, round gauge and dual matrices to explicitly
represented numbers, enclose the resulting Omega, and verify every slack
with directed error control. A constructive repair is available: replace
Y_K by Y_K + delta_K I, starting from delta_all=0, and choose

    delta_K >= max_(k outside K) [d_O_k delta_(K union {k}) + epsilon_K,k],

where an independently certified bound establishes Z_K,k >= -epsilon_K,k I.
The repaired slacks are PSD and the objective increases by delta_empty
for P=1. Floating-point eigenvalues alone cannot supply certified epsilon.

Run with `conda run -n quant_dev python research_qcqc/verification/qcqc_support_dual.py --n 3`.
JSON diagnostics and NPZ candidates remain in the local `qfi_results` folder.

## Joint gauge optimization and interval verification

The `--joint` option makes h and Y_all variables and imposes

    [[Y_all, 2 conjugate(D)], [2 D.T, I]] >= 0,
    D = dot C - i C h.

Its Schur complement gives Y_all >= Omega(h). Weak duality still applies:
the extra terminal term Tr((Y_all-Omega) S) is nonnegative. At N=3 the
joint numerical objective is 4.737311635619153, unlike the loose fixed-gauge
candidate. No minimax equality is needed to use a feasible witness as an
upper bound.

`certify_qcqc_dual.py` defines the stored lower-triangular matrix entries
as exact decimal rationals using Python's round-trip float strings. It
reconstructs their conjugate upper triangles and real diagonals. For the
half-damping channel it builds C and dot C from an outward interval for
sqrt(1/2), not from their floating-point approximations. The proof uses
mpmath interval arithmetic at 50 decimal digits and complex Hermitian LDL
elimination; every pivot must have a strictly positive lower endpoint.

For each subset K it adds

    delta_K I = (2^(N-|K|+1)-1) * 10^-6 I.

Each flow slack gains 10^-6 I, and the terminal difference gains 10^-6 I.
Interval verification checks the resulting matrices themselves, so it
does not assume any earlier floating-point eigenvalue estimate is reliable.
At N=3 all twelve flow slacks and the terminal matrix pass, giving

    F_3 <= 4737326635619153 / 1000000000000000
         = 4.737326635619153.

This is an interval-arithmetic upper certificate under the stated regular
channel/QFI conventions, trusted Python integer operations and mpmath's
outward rounding. It is not a formally verified arithmetic kernel. The
numerical primal value 4.7373116403 remains uncertified; do not present
the two values as a rigorous narrow interval yet. Tests reject an indefinite
matrix and a deliberately corrupted root dual variable.

```powershell
conda run -n quant_dev python research_qcqc/verification/qcqc_support_dual.py --n 3 --joint
conda run -n quant_dev python research_qcqc/verification/certify_qcqc_dual.py --n 3
```
