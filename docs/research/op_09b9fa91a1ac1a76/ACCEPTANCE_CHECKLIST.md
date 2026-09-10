# Acceptance Checklist

Historical checklist retained for provenance. The current requirement-level
evidence and limitations are in `GOAL_COMPLETION_AUDIT.md`.

| Requirement | Evidence | Status |
|---|---|---|
| Theorem 2/3 minimax and gauge direction | `MINIMAX_GAUGE_AUDIT.md` | Checked; full SDP minimax equality not claimed. |
| Product performance-operator expansion | `PERFORMANCE_OPERATOR_EXPANSION.md` | Algebraically checked. |
| Residual positivity and norm bound | `THEOREM3_NORM_BOUND.md`, `test_stinespring_residual.py` | Checked algebraically and numerically. |
| ICO orthogonality | `ICO_CANCELLATION_PROOF.md`, `ico_projector.py` | Finite-dimensional deterministic-process proof and tests. |
| QC-QC subset recursion | `qcqc_index_constraints.py` | Checked for N=2..7. |
| Operator-valued QC-QC recursion | `qcqc_operator_recursion.py`, `QCQC_SCOPE.md` | Tested explicit submodel; not the full SDP class. |
| Parallel attainability | `PARALLEL_ATTAINABILITY.md` | External theorem mapped conditionally; finite SDP checks reproduced. |
| Rank-changing points | `RANK_REGULARITY.md`, `test_rank_regularity.py` | Explicit proof and seven regression cases distinguish smooth-lift cost from pointwise SLD, including positive parallel QFI. General regularity theorem remains incomplete. |
| Persistent advantage search | `ADVANTAGE_SEARCH.md` | No example found; not an exhaustive disproof. |
| Catalog synchronization | `npm run sync-tex -- --check`, `node site/build.mjs` | Passed; full metadata migration remains pending. |
| Reproducibility | `requirements.txt`, `README.md` | Tested in `quant_dev`; full related suite: 54 passed. |
| Corrected QC-QC feasibility | `qcqc_feasibility.py`, `QCQC_SDP_DESIGN.md` | Implemented Eq. (63), including W PSD; 14 dedicated tests. Numerical membership only, not QFI optimization. |

No item in this checklist authorizes commit, push, upload, or a status change
from `Unsolved` to `Solved` in the catalog.
