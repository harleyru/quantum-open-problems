# Assumption Coverage Audit

WLP treats SQL linear and HL quadratic leading coefficients separately. The
SQL squeeze and positivity argument is summarized in `FINAL_STATUS.md`.

| Requirement | WLP formulation | Assessment |
|---|---|---|
| Smooth one-parameter family | Differentiable channel and Kraus tangent | Covered locally; rank-changing points need special care. |
| Finite-dimensional input/output | Explicit finite-dimensional channel model | Covered. |
| Repeated black-box uses | `E_g` tensor power and identical uses | Covered. |
| Output-state QFI | Final-state QFI optimized through a performance operator | Covered, subject to normalization conventions. |
| Parallel versus QC-QC | `Para subset QC-QC subset Gen` and QC-QC recursion | Covered at the class-definition level. |
| Eventual positive parallel QFI | SQL and HL cases treated separately | Covered conditionally at regular SLD points; see `RANK_REGULARITY.md`. |
| `limsup` ratio | Matched positive SQL or HL leading coefficients give an ordinary limit | Conditional on external attainability and regularity hypotheses. |

The earlier claim that SQL was omitted is withdrawn. Remaining coverage risks
are rank regularity, SLD conventions, and the exact hypotheses of the external
parallel attainability theorem.
