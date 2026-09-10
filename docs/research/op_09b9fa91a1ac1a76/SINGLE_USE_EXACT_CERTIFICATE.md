# Exact one-query certificate and circuit

## Statement

For the channel K_0(theta)=diag(1,s) exp(-i theta Z/2) and
K_1(theta)=sqrt(1-s^2)|0><1| exp(-i theta Z/2), with 0<s<1,
the ancilla-assisted one-query optimum is

    F_1 = 4 s^2/(1+s)^2.

One-query QC-QC has no order choice. Parameter-independent preprocessing
and postprocessing cannot exceed the ancilla-assisted channel optimum.
The proof below supplies matching explicit upper and lower witnesses.

## Upper witness

At theta=0 choose the diagonal Kraus gauge

    h = diag((s-1)/(2(1+s)), 1/2).

With D_i = dot K_i - i h_ii K_i, direct multiplication gives

    sum_i D_i^dagger D_i = s^2/(1+s)^2 I.

For any normalized input purification, apply the gauged Stinespring
isometry and retain its environment. Its pure-state QFI is
4(<dot Psi|dot Psi>-|<Psi|dot Psi>|^2), at most
4<dot Psi|dot Psi> = 4 s^2/(1+s)^2. Tracing the environment cannot increase
SLD QFI. Mixed inputs admit a purification, so the bound covers them too.
This is a direct channel upper witness; it does not require SDP strong
duality or an exchange of optimization order.

## Lower witness and physical circuit

Initialize target T and memory R in |00>. Let p=s/(1+s). Apply
R_y(2 arccos sqrt(p)) to T and CNOT from T to R. The resulting probe is

    sqrt(p)|00> + sqrt(1-p)|11>.

Apply the channel once to T; retain R. Its output is

    rho(theta) = |v(theta)><v(theta)| + (1-s^2)(1-p)|01><01|,
    v(theta) = sqrt(p) exp(-i theta/2)|00>
               + s sqrt(1-p) exp(i theta/2)|11>.

At the local operating point theta=0 measure in the orthonormal basis
(|00>+i|11>)/sqrt(2), (|00>-i|11>)/sqrt(2), |01>, |10>.
Writing w=p+s^2(1-p), the first two probabilities are w/2 and their
derivatives have opposite signs with magnitude s sqrt(p(1-p)). The
remaining probabilities have zero derivative. Their classical Fisher
information equals

    4 s^2 p(1-p)/w = 4 s^2/(1+s)^2.

The |10> probability vanishes identically in a neighborhood and contributes
zero. Thus the measurement attains the upper bound. For another known
local operating point, rotate the measurement phase accordingly.
This is a complete one-query preparation and measurement, not an
extraction of the existing multi-query optimized witnesses.

## Exact arithmetic

`single_use_certificate.py` verifies the upper and lower identities with
symbolic rational functions. For u=1-s^2=1/2, F_1=12-8 sqrt(2). Set M=10^30
and k=floor(sqrt(2 M^2)) using integer square root. Integer squaring proves

    12-8(k+1)/M < F_1 < 12-8k/M.

This strict rational interval has width 8/10^30. No floating-point
eigenvalue computation or optimizer status enters this enclosure.

```powershell
conda run -n quant_dev python research_qcqc/verification/single_use_certificate.py
conda run -n quant_dev python -m pytest -q research_qcqc/verification/test_single_use_certificate.py
```

Numerical circuit tests supplement the exact derivation; they do not supply
its rigor. Multi-query dual certificates, generic QC-QC circuit extraction,
and the asymptotic proof remain open work in this project.
