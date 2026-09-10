# Effective-support circuit extraction

The construction follows the link-inversion approach in Wechs et al.,
[Appendix B.3.3 and Lemma 17](https://arxiv.org/html/2101.08796v3#A2.SS3.SSS3),
with positivity supplied by the corrected characterization and the repaired
primal witness. The implementation currently handles qubit slots and
positive definite repaired nodes, not arbitrary ranks or dimensions.

## Matrix construction and proof

Choose square-root factors R[K,k] of each auxiliary node X[K,k]. For each
used subset K, build a matrix A_K with rows indexed by all open physical
legs of K. Its column blocks, one for each last-called k in K, consist of
the output leg O_k and the memory index of R[K\{k},k]. Insert an identity
wire between the row O_k and the column O_k. Build B_K from R[K,l] for
each next l outside K, moving I_l from the row legs into the column block.

The repaired flow equation gives A_K A_K^dagger = B_K B_K^dagger. At the
empty subset choose A_empty=[1]; at the full subset choose B_all=sqrt(S).
The physical map on the direct-sum column spaces is

    U_K = B_K.T (A_K.T)^+.

The equal Gram matrices imply U_K A_K.T = B_K.T and U_K^dagger U_K=P_K,
where P_K projects onto range(A_K.T). To see this, take thin singular-value
decompositions with the same left singular vectors and singular values,
A=L Sigma R_A^dagger and B=L Sigma R_B^dagger. Substitution gives
U=conjugate(R_B) R_A.T, and hence the stated projector identity.
Thus U_K is an isometry on the reachable input support. This remains
true for arbitrary slot operations, since contraction of open past legs
produces input vectors in that support.

The control labels form direct sums indexed by (used subset, next slot).
The simulator applies a slot Kraus operator to its input leg, then combines
incoming branches coherently before applying U_K. Kraus environment labels
belong to physical slots, not to chronological positions. No branch is
measured during this recombination.

## Three-query instance

For the certified repaired half-damping witness, dimensions are:

| Used slots | Number of maps | Input | Output | Reachable input rank |
| --- | --- | --- | --- | --- |
| 0 | 1 | 1 | 12 | 1 |
| 1 | 3 | 4 | 32 | 4 |
| 2 | 3 | 32 | 64 | 16 |
| 3 | 1 | 192 | 64 | 64 |

The terminal 192-to-64 map is not an isometry on its whole declared input
space. It is an isometry on the 64-dimensional reachable support. A physical
extension to unused input states requires additional output space or a
CPTP completion; the saved matrices do not yet contain that completion.

The circuit simulation produces QFI 4.737052445071347, inside the certified
interval. Output and derivative matrices agree with direct Choi contraction
within 1.5e-15, and the largest effective-isometry residual is below 1.6e-12.
These are floating-point checks, not interval certificates for the extracted
maps. Tests also insert distinct complex unitary channels, compare derivatives
with finite differences, and reject a deliberately misnormalized preparation.

```powershell
conda run -n quant_dev python research_qcqc/verification/extract_qcqc_circuit.py --n 3
```

`circuit_n3.npz` stores the maps; `circuit_n3.json` specifies each direct-sum
block width and diagnostics. The suffix is the integer mask of the used
subset. Within each block the slot leg precedes memory. No file is uploaded.

## Full-input CPTP completion

`complete_qcqc_circuit.py` now supplies a full-input completion. If U has
initial support projector P, choose an orthonormal basis {z_j} for its
orthogonal complement and define Kraus operators

    K_0=U, K_j=|0_out><z_j|.

Their completeness relation is U^dagger U + sum_j |z_j><z_j| = I.
Thus the map is CP and TP on all inputs. For any state supported on P,
including one entangled with a reference, all added Kraus branches vanish.
The construction preserves coherences within the reachable support.
It also gives a Stinespring isometry sum_j K_j tensor |j_env>, if a
coherent dilation rather than a CPTP representation is desired.

For floating matrices the code first uses a polar/SVD projection to set
the nonzero singular values to one. The largest entry change is below
8.6e-13. The exported NPZ stores K_0 and the rows <z_j|; each additional
Kraus operator has its sole nonzero row in output row zero. Output row zero
lies in the first declared outgoing control block. The construction need
not assign a meaningful history to unreachable inputs.

After loading the exported maps again, all local CPTP completeness checks
have residual below 1.4e-15. Random full-input density matrices preserve
trace and positivity, and tests show that omitting the complement loses
trace on an unreachable input. The replayed principal branch output differs
from the original by less than 2.7e-15 and its QFI is 4.737052445071329.

```powershell
conda run -n quant_dev python research_qcqc/verification/complete_qcqc_circuit.py --n 3
```

## All completion branches

`replay_full_circuit.py` now loads the exported completed maps and propagates
density matrices and derivatives through every completion Kraus branch.
At each chronological layer it assembles a coherent direct-sum principal
map. Complement Kraus operators carry distinct environment labels, so their
reset contributions add as density matrices. Summing their contributions as
Tr(N rho N^dagger)|0><0| is exact linear algebra, not branch pruning.
External Kraus labels remain attached to slots throughout the sum.

The full three-query output trace is 1.0000000000000007 and QFI is
4.737052445071333. The maximum difference from principal replay is below
8.4e-17. A separate test deliberately routes half the probability through
a complement and verifies that full replay has trace one while principal
replay has trace one half. Finite differences check the propagated derivative.

```powershell
conda run -n quant_dev python research_qcqc/verification/replay_full_circuit.py
```

Floating-point propagation can yield solver-scale negative eigenvalues
(here about -1.3e-16). A certified accumulated arithmetic/realization error
bound remains unfinished; the full replay is not an interval certificate
for the saved circuit gates. This does not alter the separate exact-flow
and interval-QFI process certificates. The asymptotic proof with its source
dependencies and scope is in `ASYMPTOTIC_SQUEEZE_PROOF.md` and
`FIRST_JET_REDUCTION.md`.
