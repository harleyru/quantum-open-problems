# QC-QC Implementation Scope

The executable `qcqc_operator_recursion.py` instantiates the path-Kraus
formula from WLP's `eq:qcqc_kraus_iterative` with identity or
history-dependent unitary internal maps and an explicit route register. The
register stores ordered permutations and the global Kraus operator is a
vertical block stack. Finite-difference tests verify the Leibniz derivative
and completeness tests verify normalization.

This is an operator-valued coherent superposition of fixed orders: off-
diagonal route blocks are retained before the route register is discarded. It
is not the full QC-QC class, because WLP's construction permits arbitrary
history-indexed internal isometries `V^{->k}_{K,k}` with ancillary output
spaces, whereas this implementation only permits history-dependent unitary
controls on the fixed system space. Ancilla extensions, arbitrary internal
isometries, and process-matrix optimization are not implemented. Passing these tests proves neither
the universal Theorem A recurrence nor QC-QC optimality.

The subset-index audit independently verifies allowed unused successors and
unique prefixes for all `N!` routes. The full class-level conclusion is
therefore a source-proof audit plus a tested explicit submodel, not an
independent QC-QC SDP solution for that route simulator.

## Separate process-membership implementation

`qcqc_feasibility.py` now implements fixed-W feasibility for the full
deterministic QC-QC characterization at finite external dimensions, using
Wechs et al.'s corrected Proposition 7. Its PSD auxiliaries range over all
subset nodes; they do not restrict the process to the route simulator above.
`QCQC_SDP_DESIGN.md` records the equations and the 2023 positivity correction.
This numerical membership solver does not implement QFI optimization,
explicit ancilla/isometry extraction, or prove WLP's recurrence.
