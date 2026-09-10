# Independent verification framework

**Research draft, pending independent scientific review.** Read
`REVIEW_STATUS.md` before interpreting any proof or numerical result.

Current results and requirement-level evidence: `GOAL_COMPLETION_AUDIT.md`.
The five-workstream outline below records the original audit organization;
its release-boundary list predates the newer certificates and proof documents.

This local directory audits Wei, Li and Pang, arXiv:2609.05355v1, against
problem `op_09b9fa91a1ac1a76`. It is GitHub-ready English material, but it is
not yet part of the catalog and is intentionally not uploaded.

## Evidence levels

- **Algebraic check:** an identity is evaluated symbolically or to numerical
  precision for explicitly stated finite dimensions.
- **Finite-instance check:** a process constraint or bound is solved for a
  specified N and dimension.
- **Proof audit:** the source's quantifiers, hypotheses, and inference steps
  are checked against the definitions and cited results.
- **Independent proof:** a complete argument is reproduced without relying on
  the source's omitted steps. The project does not claim this level yet.

## Five workstreams

1. Theorems 2 and 3: audit decomposition, gauge changes, norm bounds, and
   asymptotic limit.
2. ICO orthogonality: test partial-trace cancellation for valid finite process
   matrices and audit the extension to arbitrary dimensions.
3. QC-QC recurrence: implement subset-indexed causal constraints and test the
   scalar recurrence on finite instances.
4. Parallel attainability: verify the cited channel-metrology theorem and
   reproduce representative SQL/HL examples.
5. Scope: compare every quantifier and regularity assumption with the problem
   statement and record any mismatch.

Every result must state its evidence level and residual risk.

## Release boundary

The current package is an audit package, not a proof of the full WLP claim.
Before publication, a reviewer must either supply or explicitly accept the
following unresolved boundaries:

- a global minimax/QFI derivation, including the exact SDP domain and gauge
  regularity;
- a general dynamic QC-QC process-matrix construction, beyond the tested
  route-labelled operator submodel;
- a line-by-line identification of WLP's and Zhou--Jiang's QFI and resource
  conventions;
- a Choi-level treatment of parameter points where the minimal Kraus rank
  changes; and
- an exhaustive argument about persistent asymptotic QC-QC advantage.

Passing the executable tests, or finding no counterexample in the search log,
does not discharge these boundaries. They must remain visible in any catalog
contribution and must not be converted into a `Solved` status by implication.

## Evidence index

- `CIRCUIT_NORMALIZATION_BOUND.md`: exact integer Gram bounds and a cumulative
  output-state bound for rational gates versus their exact CPTP normalization;
  not a QFI or floating replay error certificate.

- `ASYMPTOTIC_SQUEEZE_PROOF.md`: CP domination, all-N upper bound and SQL/HL
  squeeze under explicit regularity and published parallel-attainability assumptions.

- `QCQC_CIRCUIT_EXTRACTION.md`: explicit effective-support maps, a full-input
  CPTP completion, and numerical replay with stated completion-branch limits.

- `QCQC_PRIMAL_CERTIFICATE.md`: exact-rational flow repair, interval PSD
  and QFI hypograph checks, and a certified three-query interval at half damping.

- `QCQC_SUPPORT_DUAL.md` and `qcqc_support_dual.py`: fixed-gauge support
  dual, telescoping weak-duality proof, and numerical candidates. These
  are not yet rigorous multi-query upper certificates.

- `SINGLE_USE_EXACT_CERTIFICATE.md` and `single_use_certificate.py`: matching
  exact one-query upper/lower witnesses, explicit probe and measurement,
  and an integer-verified rational enclosure at half damping.

- `QCQC_QFI_SDP.md`, `qcqc_qfi.py`, and `test_qcqc_qfi.py`: finite-query
  full-gauge QC-QC QFI hypograph SDP, independent output-SLD check, and
  analytic single-use regressions. These extend membership testing to
  numerical performance optimization, not an asymptotic proof.

- `FINAL_STATUS.md`: consolidated result and source anchors.
- `THEOREM_2_3_AUDIT.md`, `MINIMAX_GAUGE_AUDIT.md`, and
  `PERFORMANCE_OPERATOR_EXPANSION.md`: theorem-chain checks.
- `ICO_CANCELLATION_PROOF.md`, `ico_projector.py`, and
  `test_ico_projector.py`: cancellation proof and executable projector tests.
- `qcqc_index_constraints.py` and `qcqc_operator_recursion.py`: finite route
  and operator-valued submodel checks; see `QCQC_SCOPE.md` for the boundary
  between this model and general dynamic QC-QC.
- `PARALLEL_ATTAINABILITY.md`, `ASSUMPTION_COVERAGE.md`, and
  `RANK_REGULARITY.md`: external theorem and hypothesis coverage.
- `test_rank_regularity.py`: analytic pointwise-SLD and smooth-lift regression
  checks at a rank-changing preparation channel.
- `ADVANTAGE_SEARCH.md`: finite-use and persistent-advantage search log.
- `QCQC_SCOPE.md`: exact scope of the operator-valued QC-QC submodel.
- `QCQC_SDP_DESIGN.md`, `qcqc_feasibility.py`, and
  `test_qcqc_feasibility.py`: implemented fixed-process QC-QC feasibility
  SDP using corrected Proposition 7, including total-process positivity.
  This is not metrological QFI optimization or circuit extraction.

## Reproducible environment

The core checks use NumPy and pytest. The parallel SDP checks additionally
require CVXPY. QC-QC feasibility also uses CVXPY and SciPy. A tested local environment is the `quant_dev` conda
environment with CVXPY 1.9.2 and pytest. From the repository root:

```text
conda run -n quant_dev python -m pytest -q research_qcqc/verification
conda run -n quant_dev python -m pytest -q research_qcqc/mothe_paper/test_symmetric_parallel.py
conda run -n quant_dev python -m pytest -q research_qcqc/verification/test_rank_regularity.py
```
