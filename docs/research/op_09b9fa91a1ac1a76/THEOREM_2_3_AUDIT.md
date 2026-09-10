# Audit of Theorems 2 and 3

Source: W. Wei, Y. Li and S. Pang, arXiv:2609.05355v1, main text and
Supplemental Material, downloaded locally under `wei_li_pang_2026/`.

## Definitions checked

For a differentiable Kraus representation `E_g = {K_i(g)}`, the paper uses

`alpha = sum_i dot(K_i)^dagger dot(K_i)` and
`beta = sum_i dot(K_i)^dagger K_i`.

The completeness derivative implies `beta + beta^dagger = 0`. A Kraus
environment gauge changes the representation and hence changes alpha and
beta; every optimization over these quantities must use one common gauge.
The local regression tests enforce this point.

## Theorem 2 chain

1. The N-use performance operator is expanded into local, parallel-parallel,
   parallel-orthogonal, and orthogonal-orthogonal terms.
2. The paper defines `Lambda_parallel = 4 beta tensor I_out/d` and
   `Lambda_perp = Lambda - Lambda_parallel`.
3. The contraction of the orthogonal-orthogonal term with a valid two-party
   ICO process is asserted to vanish. The paper derives this from the process
   normalization/causality constraints in the supplemental proof.
4. Bounding the remaining local and interference terms yields

   `F_Gen(N) <= min_K [4 N ||alpha|| + N(N-1)
      (4 ||beta||^2 + 2 d ||beta|| ||Lambda_perp||)]`.

5. If a gauge has beta=0, the N^2 term vanishes. This proves that an SQL
   channel under the Hamiltonian-in-Kraus-span condition remains SQL-limited
   for general ICO.

The algebraic definitions and the beta=0 implication are locally checked for
the amplitude-damping rotation family. The universal orthogonal contraction is
independently derived from CPTP normalization in
`ICO_CANCELLATION_PROOF.md` and tested by `test_ico_projector.py`.

## Theorem 3 chain

1. The Stinespring tangent is decomposed as `dot(V)=V beta^dagger + R`, with
   `R^dagger V=0` and residual norm controlled by
   `||alpha-beta beta^dagger||`.
2. The same process-matrix orthogonality removes the cross terms carrying one
   residual factor.
3. The resulting general-ICO bound is

   `F_Gen(N) <= 4 [N ||beta|| +
      sqrt(N ||alpha-beta beta^dagger||)]^2`.

4. The paper invokes known asymptotic parallel attainability to obtain the
   matching coefficient `4 min_K ||beta||^2`.
5. Since `Para subset QC-QC subset Gen`, the squeeze argument gives the same
   leading coefficient for QC-QC whenever the common coefficient is nonzero.

Steps 1 and the positive-semidefinite residual identity are direct algebraic
checks. Step 2 is covered by the independent finite-dimensional normalization
argument, subject to the deterministic-process convention stated there. Step
4 remains the load-bearing external attainability claim; it is supported by
Zhou and Jiang (PRX Quantum 2, 010343 (2021)) but is not reproved here.

## Local evidence and residual risk

- `test_stinespring_residual.py` and the gauge tests under
  `research_qcqc/wei_li_pang_2026/`: residual and same-gauge alpha/beta
  checks for amplitude damping.
- `research_qcqc/mothe_paper/test_symmetric_parallel.py`: independent
  parallel-QFI checks for the tested finite instances.
- `qcqc_index_constraints.py` and `qcqc_operator_recursion.py`: finite
  combinatorial and operator-valued QC-QC recurrence checks.
- The local framework now implements the full reduced QC-QC constraints,
  certifies a finite three-query QFI interval, and extracts circuit maps;
  see `QCQC_PRIMAL_CERTIFICATE.md` and `QCQC_CIRCUIT_EXTRACTION.md`.
  These finite-instance results do not discharge the general asymptotic
  proof obligations. ICO cancellation is independently proved for the stated
  convention; the full asymptotic theorem remains under audit.
