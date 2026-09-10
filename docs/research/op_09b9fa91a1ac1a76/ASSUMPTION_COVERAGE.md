# Assumption Coverage Matrix

| Requirement | WLP formulation | Coverage | Residual risk |
|---|---|---|---|
| One unknown parameter | Single parameter `g` and Kraus tangent `dot K_i` | Covered | None beyond differentiability. |
| Identical black-box uses | `E_g` tensor power with N repeated queries | Covered | Resource normalization must match the problem record. |
| Finite-dimensional channel | Finite input/output spaces and finite Kraus rank | Covered | None in the stated model. |
| Smooth/regular family | Differentiable Kraus/Stinespring tangent | Conditional | A smooth nonminimal lift alone does not guarantee pointwise SLD/gauge equality; an explicit counterexample is proved in `RANK_REGULARITY.md`. |
| Output-state QFI | Performance-operator minimax formulation | Covered conditionally | Choi/vectorization and SLD conventions must agree. |
| Physical QC-QC strategies | `Para subset CS subset QC-QC subset Gen` and Theorem A | Class-level covered | Full coherent process-matrix optimization is not locally implemented. |
| Eventual positive parallel QFI | SQL and HL leading coefficients treated separately | Covered at regular SLD points | Parallel attainability theorem hypotheses must match exactly. |
| `limsup` ratio | Matching positive SQL or HL coefficient gives an ordinary limit | Conditional | This is a theorem-level conclusion, not a finite numerical extrapolation. |
| All channel families | Universal source theorem | Source-level covered | No exhaustive physical family search can replace the universal proof. |

For details on rank-changing points see `RANK_REGULARITY.md`; for the external
attainability dependency see `PARALLEL_ATTAINABILITY.md`.
