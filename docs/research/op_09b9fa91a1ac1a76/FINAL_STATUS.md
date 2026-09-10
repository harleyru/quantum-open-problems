# Independent Verification Status

Historical five-workstream snapshot. For the current three-obligation goal,
including certificates, circuit extraction and asymptotic proof, use
`GOAL_COMPLETION_AUDIT.md`. The limitations and test counts below describe
the earlier membership-SDP stage, not the current implementation.

Target: Wei, Li, and Pang, arXiv:2609.05355, and problem
`op_09b9fa91a1ac1a76`.

## Results

## Source Anchors

The audit uses the local WLP source anchors below: Theorem 2 and its bound
(`main.tex`, `eq:bound2`), Theorem 3 (`eq:bound3`), the asymptotic upper bound
(`eq:ultimate_asymptotic_upper`), the matching SQL statement (`main.tex`,
lines 1719-1722), the QC-QC recurrence (`eq:iterative_bound`), and the
explicit ICO projector (`eq:qico_appendix`) with its image characterization
(`eq:charaoficq`). These labels are preferred over page numbers because the
local source is the reproducible artifact.

1. **Theorems 2 and 3.** The Kraus gauge decomposition and norm bounds are
   algebraically consistent; the product-Kraus expansion is recorded in
   `PERFORMANCE_OPERATOR_EXPANSION.md`, while the minimax direction and
   residual positivity are independently recorded in
   `MINIMAX_GAUGE_AUDIT.md` and `THEOREM3_NORM_BOUND.md`. Theorem 2's
   interference cancellation follows
   from the deterministic-process normalization argument in
   `ICO_CANCELLATION_PROOF.md`. Theorem 3's residual bound and its SQL/HL
   consequences are source-level audited, not independently rederived in a
   full SDP proof.
2. **ICO orthogonality.** Independently proved for finite-dimensional
   deterministic reduced processes by two independent CPTP perturbations.
   `ico_projector.py` verifies the stated trace-and-replace projector; its
   tests include unequal dimensions and a complete two-qubit-slot Pauli basis.
3. **QC-QC recursion.** The subset-indexed route structure is executable and
   passes `N=2..7`: all `N!` routes, unique prefixes, and correct predecessor
   and successor counts. `qcqc_operator_recursion.py` additionally implements
   a finite route-labelled operator-valued submodel and tests trace
   preservation, product derivatives, and positive-semidefinite alpha. Its
   retained route register preserves measurable off-diagonal coherences; it
   is a coherent superposition of fixed orders, not a dynamically controlled
   QC-QC optimizer or full process SDP implementation. Tests now use two
   nonzero amplitude-damping Kraus operators and compare state derivatives
   to finite differences. The earlier claim that retained labels alone
   remove interference is withdrawn.
4. **Parallel attainability.** Local SDP/QFI calculations reproduce reported
   finite values and extend them. Zhou and Jiang's published result supports
   asymptotic parallel attainability in both SQL and HL regimes, subject to
   matching regularity, QFI, and resource conventions. See
   `PARALLEL_ATTAINABILITY.md`.
5. **Assumptions.** Smooth finite-dimensional repeated channel families,
   output-state QFI, and the parallel/QC-QC class definitions match. The
   eventual-positivity condition is compatible with the SQL squeeze at regular
   SLD points. An explicit rank-changing preparation family has pointwise
   parallel QFI N, while its single-copy smooth-lift minimum is 5 rather
   than the pointwise SLD value 1. Thus smooth lift existence and eventual
   positivity alone do not discharge regularity. This is not a QC-QC
   advantage or a counterexample at a regular point. Exact theorem
   hypotheses remain residual risks. See `ASSUMPTION_AUDIT.md` and
   `RANK_REGULARITY.md`.

The persistent-advantage search log is in `ADVANTAGE_SEARCH.md`; it records
absence of a found example, not an exhaustive disproof independent of WLP.

## Reproduction

```text
python -m pytest -q research_qcqc/verification
```

The parallel SDP tests require the local `quant_dev` conda environment
(`cvxpy 1.9.2`, `pytest`). In that environment the combined verification,
gauge, operator-recursion, rank-regularity, and symmetric-parallel suite
passes 54 tests.

## Fixed-process QC-QC feasibility

`qcqc_feasibility.py` implements the corrected Wechs Proposition 7 / Eq. (63)
with PSD subset auxiliaries and explicit positivity of W. The 2023 erratum
was checked against the publisher PDF. Fourteen dedicated tests cover
quantum switch coherence, complex matrices, unequal dimensions, nontrivial
past/future, three nonwhite three-slot orders, OCB rejection, invalid
normalization structure, and the erratum's missing-positivity failure mode.
The returned witnesses are independently checked through NumPy tensor traces
and eigenvalues. This closes the local fixed-W membership implementation
task, not QFI optimization or the independent asymptotic proof.

The current repository is intentionally uncommitted and unuploaded. Numerical
checks are consistency evidence, not proofs of global optimality.
