# QC-QC Feasibility SDP

## Sources and correction

We implement Wechs, Dourdent, Abbott and Branciard, PRX Quantum 2,
030335 (2021), corrected Proposition 7, Eq. (63).
[Corrected arXiv v3](https://arxiv.org/html/2101.08796v3#S5.SS3) and
[2023 erratum](https://doi.org/10.1103/PRXQuantum.4.030901) were checked.
The publisher's one-page erratum requires W itself to be PSD; auxiliary PSD
matrices and marginal equations do not imply this for nontrivial F.
Eq. (63) is unchanged. The sufficiency proof is Appendix B.3.c of v3.

The previous design incorrectly conflated circuit isometries with SDP
variables. V-dagger V=I is quadratic and nonconvex. It is NOT imposed on
optimization variables here. Ancillas and isometries exist by the cited
characterization theorem; extracting them is a separate task.

## Convention and equations

Input W has axes P,I0,O0,...,I(N-1),O(N-1),F, with positive integer
dimensions; P and F may be trivial. Choi vectors are unnormalized:
Tr W = dP product_k dO_k. For every strict subset K and k outside K,
introduce Hermitian PSD X[K,k] on P, input/output pairs of K, and input I_k.
The slot axes are in increasing slot order, with P first. There are
N*2^(N-1) auxiliary matrices; their sizes are
dP*dI_k*product_(j in K)(dI_j*dO_j).

Minimize zero subject to the following affine and PSD constraints:

```text
W >= 0; X[K,k] >= 0;
sum_k Tr_(I_k) X[empty,k] = I_P;
for every nonempty proper K:
  sum_(l outside K) Tr_(I_l) X[K,l]
    = sum_(k in K) X[K\{k},k] tensor I_(O_k);
Tr_F W = sum_k X[Nset\{k},k] tensor I_(O_k).
```

Here Nset={0,...,N-1}. Each identity is unnormalized. Named-axis permutation
aligns summands before addition. W is fixed. This tests full deterministic
QC-QC membership at fixed finite external dimensions, not merely membership
in general ICO or a fixed-order subclass. No link-product transpose is needed
for these equations; the general-ICO projector is not a replacement for them.

## Implementation

`qcqc_feasibility.py` constructs the SDP with CVXPY. Its NumPy audit
recomputes partial traces from tensors, and checks equation residuals,
Hermiticity, and eigenvalues of W and all auxiliary matrices. Identity-axis
permutations are separately compared with explicit index loops in tests.
Real and imaginary partial traces are separated to work around CVXPY's
complex canonicalization limitation.

`solve` returns status, accepted, and available auxiliary witnesses and audit
metrics. Acceptance requires a candidate and all numerical gates to pass.
Default entrywise tolerance is 2e-6; SCS eps is 1e-8. Neither an accepted
floating-point witness nor solver-reported infeasibility is an exact or
interval certificate. Invalid shapes and nonfinite inputs raise ValueError.
Dimension support is general, but computational cost grows exponentially.

## Verification

From the QIOP workspace root:

```text
conda run -n quant_dev python -m pytest -q research_qcqc/verification/test_qcqc_feasibility.py
```

Tests cover unequal dimensions, nontrivial P and F, three-slot white noise,
quantum switch coherence and complex phase, identity wires, deliberately
corrupted witnesses, positive invalid processes and OCB infeasibility.
The OCB fixture is `(I+(IZZI+ZIXZ)/sqrt(2))/4` in AI,AO,BI,BO order:
Oreshkov, Costa and Brukner,
[Nature Communications 3, 1092 (2012)](https://doi.org/10.1038/ncomms2076).

The erratum regression compares white noise with
`Wbad=(I_(I O)/2) tensor diag(2,-1)` for dI=dO=dF=2 and dP=1.
Both have the same future marginal. Wbad is nonpositive and must be rejected,
even though the marginal equations admit the white-noise auxiliaries.

This implementation does not optimize metrological QFI, reconstruct an
explicit dynamic circuit, or independently prove Proposition 7. The existing
route simulator remains a separate subclass experiment.
