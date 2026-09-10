# Research draft: independent review required

This is an AI-assisted proof draft and reproducible evidence package, not an
independently validated resolution of problem op_09b9fa91a1ac1a76. Completion
of the local work goal means delivery of the research artifacts, not external
confirmation that the mathematical argument is correct.

Before any claim that the open problem is solved, an independent reviewer
must audit the all-N upper bound, first-jet reduction, applicability of the
published parallel-attainability theorem, and correspondence with the
problem's resource and pointwise SLD regularity conventions. Passing tests
cannot rule out a gap in a general mathematical proof.

Finite-query interval certificates and exact algebraic checks have a
different evidence level from floating-point circuit replay. The strict QFI
interval applies to the repaired process, not to the rounded exported gates.
The gate-normalization state-distance bound is not a QFI error bound.

The parallel-attainability construction is cited, not independently
reproved. One-sided parameter boundaries and continuous-limit/Bures QFI are
outside the stated argument. See GOAL_COMPLETION_AUDIT.md for details.

Publication to a personal Fork is for review and reproducibility only. It
does not change the authoritative catalog status or imply upstream approval.
Statements in older notes that artifacts remain local describe their
pre-publication history, not the status of this Fork snapshot.

## Run the self-contained snapshot

From this directory, with the requirements installed:

```sh
python -m pytest -q .
python certify_qcqc_dual.py --n 3
python certify_qcqc_lower.py --n 3 --lower 4.7368
python replay_full_circuit.py
```

Source papers and unrelated workspace files are not bundled. Older commands
with research_qcqc/verification prefixes refer to the original local layout;
use the commands above for this snapshot.
