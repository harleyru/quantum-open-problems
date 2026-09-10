# Full-gauge QC-QC QFI optimization

## Scope

This implementation optimizes finite-query reduced QC-QC processes for phase
rotation followed by amplitude damping, with 0 < u < 1 and N <= 3.
It assumes the regular-point purification formula of Wei, Li and Pang,
[Eq. (2)](https://arxiv.org/html/2609.05355v1), and uses the corrected
Wechs constraints documented in `QCQC_SDP_DESIGN.md`.
It does not prove the asymptotic theorem or cover rank-changing points.

## Hypograph derivation

Write the unnormalized Choi ensemble as C and its derivative as D. Each
column uses input-output ordering, K.T.ravel(), with interleaved slot axes.
Expand h = sum_a x_a H_a in an orthonormal Hermitian basis with real x_a.
Include all r^2 basis elements, including imaginary off-diagonal elements.
Set A_0 = D and A_a = -i C H_a. For a reduced process S define

    G_ab(S) = 4 Re Tr[(A_a A_b^dagger)^T S].

Then f(S,x) = [1,x]^T G(S) [1,x]. All entries of G are real affine
functions of S. Write G = [[a,b^T],[b,Q]]. The condition

    G(S) - t e_0 e_0^T >= 0

is equivalent to t <= inf_x f(S,x). The generalized Schur complement
requires Q >= 0, b in range(Q), and a-t >= b^T Q^+ b. Equivalently,
nonnegativity for all [1,x] implies nonnegativity of the homogeneous
quadratic form by scaling and continuity, including vectors with first
coordinate zero. This formulation does not exchange max and min.

Maximize t subject to this LMI and the corrected reduced QC-QC constraints
with P=F=1. Here F=1 labels the reduced operator, not a restriction to
discarding the physical future output. A purification of S supplies an
unrestricted retained future; its partial trace obeys the same constraints,
and the corrected characterization supplies realizability. Extracting a
gate-level realization from that characterization remains separate work.

## Numerical checks

The code independently contracts process constraints with NumPy, reconstructs
the minimizing gauge using Q^+, and evaluates the original performance
operator. It also forms V = sqrt(S^T) C, its derivative, and computes the
SLD QFI of rho = V V^dagger by spectral decomposition. The latter calculation
clips small negative process eigenvalues introduced by the solver and checks
normalization. It is a numerical cross-check, not an exact feasible repair.

The acceptance gate requires process and hypograph positivity within 2e-5,
affine residuals within 2e-5, gauge stationarity, normalization, and agreement
between objective, gauge recomputation and output SLD QFI within 2e-5.
SCS uses eps=1e-7 and at most 100000 iterations. An accepted result is not
a rigorous upper/lower certificate or a bound on distance to the exact optimum.

At u=0.5, initial runs give N=2: 2.2994501356 and N=3: 4.7373116403.
The N=3 value agrees with the previously recorded literature value 4.7373.
Tests compare N=1 with 4(1-u)/(1+sqrt(1-u))^2 at u=0.1,0.5,0.9 and check
the complex-gauge expansion against direct matrix contraction on random PSD S.

## Reproduction

From the workspace root:

```powershell
conda run -n quant_dev python -m pytest -q research_qcqc/verification/test_qcqc_qfi.py
conda run -n quant_dev python research_qcqc/verification/qcqc_qfi.py --n 3 --output-dir research_qcqc/verification/qfi_results
```

The JSON contains scalar diagnostics. The NPZ contains S, the minimizing
gauge, and all subset-indexed auxiliary matrices. A key `node_mask_k` uses
an integer bit mask for the used subset and k for the next slot.

## Remaining work

- Obtain independently audited dual bounds or interval-certified bounds.
- Extract explicit sequential isometries and compare an operational circuit.
- Extend beyond this channel family and test regularity hypotheses.
- Optimize storage before increasing N: the dense Gram coefficient tensor
  scales as (r^2+1)^2 times the square of the process dimension.
- Keep finite-query advantages separate from persistent asymptotic advantage.

All files and witnesses remain local; this work does not modify Draft PR #40.
