# Constant-rank reduction for pointwise SLD metrology

## Lemma

Let J(t) be a two-sided differentiable finite-dimensional CPTP Choi curve
on an open neighborhood of t=0. There exists a smooth CPTP curve J_tilde(t)
of locally constant Choi rank such that J_tilde(0)=J(0) and
dot J_tilde(0)=dot J(0). The construction does not preserve second derivatives
or continuous-limit QFI; it preserves pointwise SLD QFI.

## Proof

Write J=CC^dagger with C full column rank, let P project onto range(J),
Q=I-P, and D=dot J. For every x in ker(J), the real scalar function
x^dagger J(t)x is nonnegative on both sides of zero and vanishes at zero.
Its derivative is zero. Polarization yields Q D Q=0.

Define

    L=(I-P/2) D C (C^dagger C)^(-1).

Since C(C^dagger C)^(-1)C^dagger=P,

    LC^dagger+CL^dagger = DP+PD-PDP = D,

where the last equality uses QDQ=0. Set C_raw(t)=C+tL; full column rank
persists near zero. Unvectorize its columns as Kraus matrices A_i(t).
Their normalization matrix M(t)=sum_i A_i(t)^dagger A_i(t) obeys M(0)=I
and dot M(0)=0 because Tr_out J=I and Tr_out D=0. Therefore M(t)>0 near
zero and M(t)^(-1/2)=I+O(t^2). Define K_i(t)=A_i(t)M(t)^(-1/2).
They give a smooth CPTP channel, retain the value and derivative of the raw
factor, and preserve its number of linearly independent Kraus operators.
This proves the lemma. No square root of a vanishing Choi eigenvalue is used.

## Equality of metrological objectives at every finite N

For any fixed parameter-independent strategy, output density matrices are
linear in the tensor product of slot Choi matrices. Equal J and dot J
therefore give equal output rho and dot rho for every N. Pointwise SLD QFI
is a function of this pair, so it agrees exactly for the two families,
strategy by strategy. The admissible parallel, QC-QC, and deterministic
process classes do not depend on the family. Taking their suprema preserves
the equality separately for every N, and therefore preserves any defined
asymptotic ratio. N-dependent optimal strategies create no exception.

Thus the constant-rank theorem in `ASYMPTOTIC_SQUEEZE_PROOF.md` transfers
to smooth interior points even when the original Choi rank changes,
provided the cited parallel attainability theorem applies to the surrogate.
The rank-birth example in `RANK_REGULARITY.md` does not contradict this:
its smooth nonminimal lift can overestimate pointwise QFI, while the surrogate
discards rank-birth terms of second order and preserves the true first jet.

## Remaining boundary

The catalog statement uses 'regular' without defining it. The combined
argument covers two-sided smooth interior parameter points with ordinary
pointwise SLD QFI and the stated ancilla-assisted resource model. It does not
cover a one-sided boundary tangent with QDQ nonzero, continuous-limit/Bures
QFI conventions that retain second-order rank-birth information, or changing
channel dimensions. This is a precise scope boundary, not a counterexample.
The parallel attainability theorem remains an explicitly cited dependency;
this document does not reprove its QEC/probe construction.

## Executable evidence

`channel_first_jet.py` implements the factor and CPTP normalization.
Tests check a complex two-Kraus channel, finite-difference derivatives,
constant surrogate rank, an exact symbolic rank-birth example with nonzero
derivative, and rejection of a one-sided rank-birth tangent. The mathematical
proof above, not those finite tests, supplies the all-dimensional claim.
