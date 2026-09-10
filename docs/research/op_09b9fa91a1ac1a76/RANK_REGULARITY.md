# Rank-Regularity Scope

The channel family is smooth as a Choi-matrix-valued map, but a globally
differentiable minimal Kraus representation need not exist at a point where
the Choi rank changes. The WLP formulas use a differentiable Kraus lift and
an environment gauge. A constant-rank neighborhood removes this particular
rank-change obstruction, subject to the other theorem hypotheses. Merely
supplying a differentiable nonminimal lift does NOT establish equality of
the lift-gauge cost with pointwise output-state SLD QFI.

At a rank-changing point, the correct invariant object is the derivative of
the Choi operator and the output-state SLD QFI. A fixed nonminimal Kraus lift
can often provide the needed derivative, but this requires an explicit local
construction; minimal-rank gauge formulas cannot simply be assumed across
the singularity.

The problem's phrase "regular parameter value" needs an intrinsic definition,
not an exclusion based only on a chosen Kraus representation. The examples
below distinguish pointwise SLD QFI from a smooth-lift metric. This is a scope
issue, not evidence of a physical QC-QC advantage or a refutation of WLP on
regular points.

## Explicit example

Consider the preparation channel from a one-dimensional input space to a
two-dimensional output space

```text
E_t(rho) = (1-t^2) |0><0| + t^2 |1><1|,
```

for `|t| < 1`. Its Choi matrix is smooth in `t`, while its Choi rank is one
at `t=0` and two away from zero. A fixed two-Kraus representation uses
`K_0(t)=sqrt(1-t^2)|0><*|` and `K_1(t)=t|1><*|`, where `<*|` is the unique
input bra. They satisfy `K_0(t)^dagger K_0(t)+K_1(t)^dagger K_1(t)=1` and
are smooth across the rank change even though the representation is
nonminimal at `t=0`. This demonstrates why a
Choi-smooth family need not have a smooth *minimal-rank* Kraus chart, while a
nonminimal lift can still make the tangent formulas meaningful.

Any application of the WLP/Zhou--Jiang gauge minimization at such a point
must specify whether it uses this fixed lift or a rank-adapted representation;
the two descriptions are related by an isometry only after the environment
dimension is chosen consistently.

## Independent pointwise calculation

Write `rho(t)=diag(1-t^2,t^2)`. The SLD equation is
`dot(rho)=(rho L+L rho)/2`, with `F_SLD=Tr(rho L^2)`. At `t=0`,
`dot(rho)=0`; choosing `L=0` gives `F_SLD(0)=0`. For `t != 0`,
the eigenbasis is fixed and both probabilities are positive, so

```text
F_SLD(t) = (-2t)^2/(1-t^2) + (2t)^2/t^2 = 4/(1-t^2).
```

Thus the point value is zero but its punctured limit is four. The smooth
lift has `K_0(0)=|0>`, `K_1(0)=0`, `dot K_0(0)=0`, and
`dot K_1(0)=|1>`. For any finite Hermitian environment generator `h`,

```text
dot K_0^h = -i h_00 |0>,
dot K_1^h = |1> - i h_10 |0>,
alpha_h = 1 + h_00^2 + |h_10|^2,
beta_h = i h_00,
min_h 4 alpha_h = min_h 4(alpha_h-|beta_h|^2) = 4.
```

Both minima are attained, for example by `h=0`. The new-output component
`|1>` cannot be removed by mixing the Kraus values at zero: their span is
only `span{|0>}`. This is an algebraic proof, not evidence from an optimizer.
It proves that a smooth lift alone is insufficient to identify the gauge
minimum with the pointwise SLD QFI. It does not invalidate its upper-bound
use at this point.

## Positive-QFI variant

Prepare `sigma(t)=|psi(t)><psi(t)| tensor rho(t)`, with
`psi(t)=cos(t/2)|0>+sin(t/2)|1>`. This remains a smooth finite-dimensional
preparation channel and has a smooth two-Kraus lift
`L_i(t)=psi(t) tensor K_i(t)`.

At zero, `dot psi=|1>/2` and `dot rho=0`. The SLD of sigma can be chosen
as the pure-state SLD on the first factor tensored with identity. It gives
`F_SLD,sigma(0)=1`. Away from zero, tensor-product QFI additivity gives
`F_SLD,sigma(t)=1+4/(1-t^2)`, whose limit is five. In the lift calculation,
`dot L_0=|1,0>/2` and `dot L_1=|0,1>` are orthogonal to the span of the
nonzero Kraus values. Therefore

```text
alpha_h = 5/4 + h_00^2 + |h_10|^2,
min_h 4 alpha_h = min_h 4(alpha_h-|beta_h|^2) = 5.
```

For N parallel uses the input space is one-dimensional. Any clean ancilla
is independent of the parameter; a parameter-independent final channel
cannot improve QFI. Keeping all outputs attains
`F_PAR^(N)(0)=N`. One can verify this without a limit argument by using
the sum of the N single-copy SLDs; cross terms have zero expectation.
Hence eventual positive parallel QFI does not by itself exclude this
rank-change mismatch. The example has discontinuous QFI, so it does not
establish a counterexample at a regular parameter value, nor any QC-QC
advantage. A general Choi-level formulation remains to be established.

## Evidence and source

`test_rank_regularity.py` checks the analytic derivatives, the spectral SLD
formula on both sides of zero, the exact Hermitian-gauge identities, and
N=1,2,3 tensor-product pointwise QFI. These are regression tests for the
explicit proof above, not a general rank-regularity theorem.

The rank-change distinction is consistent with D. Safranek,
"Discontinuities of the quantum Fisher information and the Bures metric,"
Physical Review A 95, 052320 (2017),
[doi:10.1103/PhysRevA.95.052320](https://doi.org/10.1103/PhysRevA.95.052320).
