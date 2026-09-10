# Exact cumulative physical-normalization bound

## Rational gate definition

For each exported completed map stack its principal rows above its null
rows to form A. This is the compressed Stinespring matrix: principal output
rows share one environment label; each null row occupies output zero with
its own environment label. This embedding preserves inner products.
Round real and imaginary parts to denominator D=10^15, using exact rational
rounding of the binary float value. The resulting integer arrays define a
new exact rational circuit, not a claim that the original floats are exact
isometries. `rational_circuit_n3.npz` contains those integers as strings.

Compute G=A^dagger A-I using integer matrix multiplication, with denominator
D^2. For every row sum |Re G_ij|+|Im G_ij|, and let epsilon be the maximum.
Hermiticity and the induced matrix norm bound imply ||G||_2<=epsilon.
No floating eigenvalues enter this computation.

## Exact normalized realization

When epsilon<1, A^dagger A is positive definite. Define the exact physical
isometry B=A(A^dagger A)^(-1/2). Its existence, exact matrix formula, and the
rational entries of A specify an unambiguous realization, even though its
entries need not be rational. Its associated CPTP channel includes the
same environment labels used in full-branch replay.

Singular-value calculus gives

    ||B-A|| <= epsilon/(1+sqrt(1-epsilon)) <= epsilon.

For any input state with any reference, expand the difference between
A rho A^dagger and B rho B^dagger. The trace norm is at most
(||A||+||B||)||A-B|| <= 3 epsilon. Tracing the environment cannot enlarge
this bound. The raw CP map has induced completely bounded trace norm at
most ||A||^2<=1+epsilon.

## Layer composition

Use the largest block epsilon for each chronological direct-sum control
layer, including preparation and the terminal map. Complete the external
channel calls with their slot-indexed Stinespring environments; they are
isometries on the legal unused-slot control space. Their norm is one.
A telescoping replacement of raw control layers by normalized layers gives
the conservative uniform output-state bound

    ||rho_raw-rho_physical||_1
       <= 3 (sum_l epsilon_l) product_l(1+epsilon_l).

The verifier evaluates the right-hand side as a Fraction and proves it is
less than 10^-10 for the three-query saved witness. The bound includes all
completion branches and arbitrary entangled inputs in the allowed routing
space. It is not an estimate obtained by comparing two floating simulations.

```powershell
conda run -n quant_dev python research_qcqc/verification/certify_circuit_normalization.py
```

## Scope and remaining audit

This certifies proximity of an exact rational raw circuit and an exact
normalized physical circuit. It does not bound the difference from the
unrounded optimized process, or floating BLAS replay roundoff. It also does
not turn a state-distance bound into a QFI bound: QFI near small eigenvalues
requires derivative and spectral-support control. The rigorous QFI interval
continues to refer to the exactly repaired process witness, whose realization
exists by the algebraic extraction proof. A numerical QFI for the exported
rounded gates remains numerical evidence. Final reporting must preserve
these distinctions.
