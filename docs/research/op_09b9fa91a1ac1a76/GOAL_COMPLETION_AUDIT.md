# Goal evidence audit: 2026-09-10

This is the current status of the three-obligation research goal. Earlier
five-workstream summaries predate the implementation and are historical.
All new research remains local; the catalog checkout and PR are unchanged.

## Finite-query certificates

| N | Half-damping QFI result | Evidence |
| --- | --- | --- |
| 1 | F=12-8 sqrt(2) | `SINGLE_USE_EXACT_CERTIFICATE.md` |
| 2 | 2.2991 <= F <= 2.29945713771398 | `qfi_results/certified_lower_n2.json`, `qfi_results/certified_dual_n2.json` |
| 3 | 4.7368 <= F <= 4.737326635619153 | `QCQC_PRIMAL_CERTIFICATE.md`, `QCQC_SUPPORT_DUAL.md` |

The primal uses exact rational flow repair and a full-gauge hypograph. The
dual uses a joint gauge LMI and a telescoping weak-duality proof. Interval
LDL verifies PSD with outward arithmetic. Tests reject corrupted witnesses
and excessive lower bounds. The trust base includes Python integers and
mpmath intervals, not a formally verified arithmetic kernel.

## Circuit realization

`QCQC_CIRCUIT_EXTRACTION.md` and `extract_qcqc_circuit.py` supply a constructive
support-isometry proof and eight exported maps for the repaired N=3 process.
`complete_qcqc_circuit.py` supplies full-input CPTP completions.
`replay_full_circuit.py` independently loads the exported matrices and
propagates density matrices and derivatives through all completion branches.
The numerical QFI is 4.737052445071333. Tests include distinct complex
unitary calls, finite-difference derivatives, and a nonnegligible completion
branch. `CIRCUIT_NORMALIZATION_BOUND.md` supplies an integer-certified
cumulative state-distance bound below 10^-10 for exact rational gates
versus their exact polar-normalized realization.

This fulfills constructive extraction and independent numerical action/QFI
verification. It does not certify the rounded-gate QFI by interval arithmetic
or bound BLAS replay roundoff. The strict QFI interval refers to the repaired
process; the normalization certificate compares different explicitly defined
circuit objects. These distinctions are necessary in any final report.

## Asymptotic proof

`ASYMPTOTIC_SQUEEZE_PROOF.md` supplies the CP-domination lemma, residual
cancellation, all-N upper bound, and SQL/HL squeeze. It explicitly invokes
Zhou--Jiang's published parallel attainability theorem instead of independently
reconstructing their QEC/probe construction. `FIRST_JET_REDUCTION.md` extends
the argument to rank-changing original families at two-sided smooth interior
points for pointwise SLD QFI. The proof gives ratio one in this stated model.

Excluded conventions are one-sided boundary tangents, continuous-limit/Bures
QFI, postselection, parameter-dependent controls and changing dimensions.
The catalog leaves 'regular' undefined, so publication must state the
covered assumptions. No claim covers every possible interpretation of that
word. No source theorem is replaced by finite numerical evidence.

## Requirement mapping

- Rigorous finite-query upper/lower bounds: exact/interval witnesses for
  N=1,2,3 and constructive primal/dual derivations.
- Circuit/isometry realization: explicit matrices, support proof, full-space
  completion, independent all-branch numerical QFI/action replay.
- Asymptotic question under explicit assumptions: full local upper/squeeze
  and first-jet argument, with published parallel attainability as a named
  mathematical dependency.
- Necessary English files and reproduction: proof notes, Python implementations,
  tests, JSON diagnostics and NPZ witnesses in this local verification folder.
- No upload: clean catalog checkout; no commit, push or PR change this goal.

The stronger optional tasks of a proof-assistant formalization, independent
reproof of the external attainability theorem, and an interval-QFI certificate
for rounded exported gates are not achieved. They must not be implied by
completion of the stated-assumption research package. Publication and external
scientific review remain separate from local implementation and verification.
