# Verification Status

Historical five-workstream snapshot. Current proof and implementation status:
`GOAL_COMPLETION_AUDIT.md`.

This matrix records the current evidence for the five requested checks. It
does not upgrade a source-level audit into an independent proof.

| Workstream | Current result | Evidence level | Remaining limitation |
|---|---|---|---|
| Theorems 2 and 3 | The gauge decomposition, beta cancellation, norm recurrence, and asymptotic implication are internally consistent under the stated hypotheses. | Proof audit | Several general-ICO steps still rely on the arXiv manuscript's lemmas and cited results. |
| ICO orthogonality | `ICO_CANCELLATION_PROOF.md` derives cancellation directly from deterministic CPTP normalization in arbitrary finite dimensions. Dense projector tests check complex matrices, unequal dimensions, and a complete two-qubit-slot Pauli basis. | Independent lemma proof and finite matrix checks | Applies to deterministic reduced processes, not arbitrary positive matrices or postselection. Full QFI theorem and vectorization remain separate audit obligations. |
| QC-QC recursive constraints | Subset nodes, immediate predecessors/successors, ordered routes, and unique prefixes pass for N=2..7. | Finite combinatorial check | This audits the index structure only, not the full operator-valued QC-QC SDP. |
| Parallel attainability | Local parallel-QFI computations reproduce reported N<=3 values and extend the calculation to N=4; symmetric probes improve QFI per use over the tested finite range. | Finite-instance check | These data do not establish super-SQL asymptotic scaling. The external theorem requires full-text hypothesis checking. |
| Assumption coverage | WLP explicitly treats SQL linear and HL quadratic coefficients. See `FINAL_STATUS.md` and `RANK_REGULARITY.md`. | Conditional proof audit | The claimed SQL omission is withdrawn. Rank regularity, SLD conventions, and external attainability hypotheses remain under review. |

## Reproduction

From this directory, run:

```text
python qcqc_index_constraints.py
pytest -q test_qcqc_index_constraints.py
```

All formal claims remain English and must be synchronized with the catalog
JSON/TeX files before any eventual upload.
