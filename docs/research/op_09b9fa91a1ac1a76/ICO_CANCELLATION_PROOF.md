# Independent Two-Slot Cancellation Argument

## Assumptions and Convention

Let S be a deterministic process with no retained global future. Use a
Choi convention in which normalization reads Tr[S (J_1 tensor ... tensor
J_N)] = 1 for arbitrary independently chosen CPTP Choi matrices J_k.
Consistent full transposes can be absorbed into the definition of J_k.
Each slot has arbitrary fixed finite input/output dimensions. This is the
deterministic process class used for the reduced strategies, not a
postselected process.

## Proof Without a Projector

First take two slots. Let X and Y be Hermitian with zero output partial
trace in their respective slots. The completely depolarizing Choi matrix
J_0 = I_in tensor I_out / d_out is positive definite and has output
partial trace I_in. Thus J_0 + t X and J_0 + u Y are CPTP for all real
t and u in sufficiently small open intervals about zero.

Normalization holds for every such pair. Expanding it as a polynomial in
t and u, the coefficient of t u must vanish:

```text
Tr[S (X tensor Y)] = 0.                         (1)
```

For non-Hermitian X and Y, their Hermitian real and imaginary parts still
have zero output partial trace. Apply (1) to all four Hermitian pairs and
use complex bilinearity. In particular Y may be an adjoint.

For N slots, insert arbitrary fixed CPTP matrices in the other N-2 slots
before making the same two independent perturbations. Normalization still
holds for all small t,u, so the same coefficient vanishes with those fixed
matrices included. Equivalently, the resulting two-slot functional is
normalized on every CPTP pair. This justifies reduced-process cancellation
without assuming it as a separate lemma.

The argument holds for arbitrary finite input and output dimensions. It
does not require the perturbations to be physical channel derivatives;
zero output partial trace suffices because J_0 is an interior CPTP point.

## Application to WLP

For Lambda_perp = Lambda - D_out(Lambda), the output partial trace is zero.
Equation (1) therefore removes the perpendicular-perpendicular term in
Theorem 2. Theorem 3 residual interference operators satisfy the same
partial-trace condition when sum_i R_i^dagger K_i = 0, so the argument
also applies there, subject to checking vectorization consistently.

This proves the cancellation lemma, not the whole QFI theorem. It does not
verify minimax steps, the Kraus derivative convention at rank changes,
the parallel attainability theorem, or the QC-QC operator recursion.

## Executable Cross-Check

`ico_projector.py` separately implements the trace-and-replace formula
`eq:qico_appendix` in the WLP Supplemental. `test_ico_projector.py` checks
self-adjointness, idempotence, and two-party cancellation for complex
matrices, including unequal input/output dimensions. Random kernel
perturbations of white noise yield positive normalized processes, avoiding
the invalid use of arbitrary density matrices as process matrices.

Run from this directory: `python -m pytest -q test_ico_projector.py`.
