# Parallel Attainability Audit

WLP cites Zhou and Jiang, *Asymptotic Theory of Quantum Channel Estimation*,
PRX Quantum 2, 010343 (2021), DOI: 10.1103/PRXQuantum.2.010343. The
published paper states this more precisely. Its Theorem 1 classifies the
entanglement-assisted repeated-channel QFI as `Theta(N^2)` iff HNKS holds and
`Theta(N)` otherwise. Its Theorem 2 gives `lim F_N/N = 4 min_{h: beta=0}
||alpha||` in the SQL case and constructs arbitrarily accurate SDP-computable
parallel probes; its Theorem 3 gives `lim F_N/N^2 = 4 min_h ||beta||^2` in the
HL case and constructs probes attaining that coefficient. The SQL theorem
also states asymptotic equality of sequential and parallel QFI.

## Scope comparison

This supports WLP's use of parallel attainability for the ordinary repeated
channel model, provided the channel family is differentiable and finite
dimensional and the same resource/QFI convention is used. Zhou--Jiang use a
minimal Kraus rank at the parameter value and an arbitrarily large clean
ancilla, with a local Kraus/gauge tangent. At rank-changing points WLP needs
an additional justification identifying the gauge cost with pointwise SLD
QFI; a differentiable nonminimal lift alone is not sufficient. The explicit
counterexample to that identification is in `RANK_REGULARITY.md`. This paper
does not by itself prove WLP's general-ICO theorem.

The citation therefore supports the lower-bound/attainability direction,
not the ICO upper-bound direction. WLP's `min ||beta||` coefficient also
requires the HNKS/SQL case split and the Kraus-gauge identification used in
its own proof.

Reference: https://doi.org/10.1103/PRXQuantum.2.010343
