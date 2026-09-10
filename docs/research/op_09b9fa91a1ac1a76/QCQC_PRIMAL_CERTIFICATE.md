# Exact-flow primal repair and certified QFI lower bound

For N=3 and half damping the current interval-arithmetic certificates give

    5921/1250 <= F_3(QC-QC) <= 4737326635619153/10^15,
    4.7368 <= F_3(QC-QC) <= 4.737326635619153.

The width is 0.000526635619153. This statement uses the regular-point QFI
formula, corrected QC-QC characterization, and the interval arithmetic
trust assumptions stated in `QCQC_SUPPORT_DUAL.md`. It does not establish
an asymptotic result or provide explicit multi-query circuit isometries.

## Exact affine repair

Read the lower triangles of numerical X[K,k] as decimal rationals and
reconstruct exact Hermitian matrices. Treat real and imaginary components
as separate arrays of Python Fraction values. Visit subsets by increasing
cardinality. For the current subset K let E_K be the required incoming
operator minus the sum of outgoing input traces. For every k outside K add

    E_K tensor I_I_k / (2 (N-|K|))

to X[K,k], with canonical tensor ordering. The sum of its input traces is
exactly E_K. Later repairs never alter already repaired outgoing nodes.
Finally define S by the terminal sum of X tensor output identities.

This gives exact affine feasibility but does not guarantee positivity.
Mix with the white feasible family at weight eta=1/10000:

    X_white[K,k] = |K|! (N-|K|-1)! / (N! 2^(|K|+1)) I,
    S_white = I/2^N.

The implementation then checks every affine equation exactly over rational
numbers, and verifies strict positivity of all X and S by outward interval
LDL at 50 decimal digits. The exported NPZ uses fraction strings rather
than pickled objects or floating approximations.

## QFI hypograph certificate

Construct the exact half-damping Choi ensemble and derivative using an
interval enclosure of sqrt(1/2). Use all r^2 Hermitian gauge basis elements;
the off-diagonal basis is unnormalized to avoid unnecessary square roots.
For A_0=dot C and A_a=-i C H_a form

    G_ab = 4 Re sum_(i,j,k) A_a[i,k] conjugate(A_b[j,k]) S[i,j].

The contraction includes the performance-operator transpose. A positive
matrix G-t e_0 e_0^T implies f(S,h)>=t for every Hermitian gauge h.
At t=5921/1250 the entire 65 by 65 matrix passes interval LDL, so the
repaired process supplies the claimed lower bound. Its positive S and the
locally constant full-column-rank channel ensemble avoid a rank-changing
output purification in this instance. The characterization guarantees a
QC-QC realization of the purified reduced process; extracting its explicit
isometries remains a separate, unfinished obligation.

Tests compare the interval Gram construction with a direct contraction,
check exact Fraction retention, and reject an excessive lower bound 4.74.

```powershell
conda run -n quant_dev python research_qcqc/verification/repair_qcqc_primal.py --n 3
conda run -n quant_dev python research_qcqc/verification/certify_qcqc_lower.py --n 3 --lower 4.7368
```
