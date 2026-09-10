# Asymptotic squeeze proof under explicit regularity assumptions

Extension: `FIRST_JET_REDUCTION.md` constructs a locally constant-rank CPTP
surrogate at every two-sided smooth interior point. It preserves each finite-N
pointwise SLD objective, extending the theorem below beyond original families
of locally constant Choi rank. One-sided and continuous-limit-QFI conventions
remain outside this extension.

## Statement and scope

Fix a finite-dimensional one-parameter CPTP channel family at theta_0. Assume
a differentiable minimal Kraus family of locally constant rank, a nonzero
channel derivative, and the ancilla-assisted parallel attainability theorem
of Zhou and Jiang in this channel model. Strategies use N identical slots,
parameter-independent controls at the local operating point, no postselection,
and unrestricted finite ancillas at each N. Write F_Para(N), F_QCQC(N) and
F_Gen(N) for the optimized pointwise SLD QFI.

Then F_QCQC(N)/F_Para(N) tends to one. The proof below also bounds the
normalized deterministic-process superset Gen. It uses the published
parallel attainability theorem, not a new independent proof of that theorem.
It does not cover changing dimensions, nonidentical channel uses, unknown
parameter-dependent controls, or rank-changing operating points.

## A CP domination lemma

Absorb Choi transposes into T=S.T, so T>=0 and
Tr[T (J_1 tensor ... tensor J_N)]=1 for independently chosen CPTP Choi
operators. Hold all but slot j fixed at CPTP channels. Let L have Kraus
operators L_a and A=sum_a L_a^dagger L_a<=c I. If c>0, the map L/c is
trace nonincreasing. Choose a fixed output unit vector |0> and Kraus
operators M_b=|0><b| sqrt(I-A/c). Then
sum_b M_b^dagger M_b=I-A/c, so L/c plus M is CPTP.
Positivity of T and the other Choi operators implies

    0 <= Tr[T (... tensor J_L tensor ...)] <= c.

If c=0 then all L_a vanish. This proves the lemma without a matrix
inequality between J_L and a multiple of the original channel Choi matrix.

## Tangent decomposition

Choose any differentiable Kraus gauge and stack its Kraus operators into V.
Set alpha=dot V^dagger dot V and beta=dot V^dagger V. Differentiated trace
preservation gives beta^dagger=-beta. Define

    P=V beta^dagger, R=dot V-P.

Direct multiplication gives V^dagger R=0 and
R^dagger R=alpha-beta beta^dagger>=0. Let P_j and R_j be N-slot Choi
ensemble derivatives with P or R replacing K in slot j. Use the seminorm
||A||_T^2=Tr(A^dagger T A). The CP domination lemma gives

    ||P_j||_T^2 <= ||beta||^2,
    ||R_j||_T^2 <= ||alpha-beta beta^dagger||.

For j!=k, the cross contraction <R_j,R_k>_T contains one local operator
sum_a |R_a>><<K_a| and its opposite-slot adjoint. Their output partial
traces vanish because sum_a K_a^dagger R_a=0 (with the input transpose
fixed by vectorization). All remaining slots contain CPTP Choi operators.
The normalization polynomial argument in `ICO_CANCELLATION_PROOF.md`
therefore gives <R_j,R_k>_T=0, including its complex part.

Consequently, triangle inequality and the orthogonality identity imply

    ||sum_j (P_j+R_j)||_T
       <= N ||beta|| + sqrt(N ||alpha-beta beta^dagger||).

Contracting a purification of the fixed process with the gauged channel
Stinespring maps yields an output purification. The pure-state QFI is at
most four times its squared derivative norm; discarding its environment
cannot increase SLD QFI. Thus, for every fixed gauge and every strategy,

    F_Gen(N) <= 4 [N ||beta||
                         + sqrt(N ||alpha-beta beta^dagger||)]^2.       (A)

Only an upper bound on QFI is used here, not equality with a gauge minimum.
The gauge may be selected before maximizing over strategies.

## SQL squeeze

If some Hermitian gauge has beta=0, define
a_SQL=inf_(h:beta_h=0) ||alpha_h||. Equation (A) gives
F_Gen(N)<=4 N a_SQL. For an infimum, apply the bound to every epsilon-optimal
fixed gauge and then send epsilon to zero; attainment is unnecessary.
Zhou--Jiang's Theorem 2 identifies the parallel limit F_Para(N)/N=4 a_SQL.
Since Para is contained in QC-QC and QC-QC in Gen, all three normalized
limits agree.

The coefficient is positive under the nonzero channel-derivative assumption:
a maximally entangled probe produces the normalized Choi state, whose SLD
QFI is positive when its derivative is nonzero. N independent copies give
F_Para(N)>=N F_Choi(1)>0. Thus division by F_Para(N) is legitimate and the
ratio tends to one. This argument avoids claiming that a zero HL coefficient
leaves the positive SQL case unresolved.

## HL squeeze

If no gauge has beta=0, set b=inf_h ||beta_h||. The image of the finite
dimensional real-linear gauge map is a closed subspace, so its affine
translate has positive distance from zero: b>0.
For each epsilon>0 choose one fixed gauge with ||beta||<=b+epsilon.
Divide (A) by N^2 and take limsup in N before sending epsilon to zero.
The residual norm of that fixed gauge is finite, giving

    limsup_N F_Gen(N)/N^2 <= 4 b^2.

Zhou--Jiang's Theorem 3 supplies the parallel limit 4 b^2. Inclusion again
gives matching limits and ratio one. No N-dependent gauge convergence,
uniform residual bound over gauges, or interchange of limits is required.

## Degenerate and excluded cases

If the channel derivative vanishes, the derivative of its N-fold tensor
power and hence every parameter-independent process output vanishes.
Pointwise SLD QFI is zero for all N; a positive-denominator advantage ratio
is not defined. Rank-changing examples in `RANK_REGULARITY.md` require
separate treatment: smooth lifts alone do not justify the attainability
formula used here. The current argument does not remove that restriction.

## External dependency and evidence

Zhou and Jiang, PRX Quantum 2, 010343 (2021), Theorems 2 and 3, give the
parallel SQL and HL coefficients in the ancilla-assisted repeated-channel
model. Their beta convention differs by an adjoint and/or factor i from
the anti-Hermitian beta here; the zero constraint and operator norm agree.
Primary source: https://journals.aps.org/prxquantum/pdf/10.1103/PRXQuantum.2.010343

The CP domination and fixed-gauge squeeze above complete the argument
conditional on that published attainability input and the stated regular
channel assumptions. This is not an independent reproof of all source
theorems, nor an unconditional resolution of the catalog's intended meaning
of 'regular'. Numerical tests check the contraction identities on complex
two-slot examples; they do not replace this all-N reasoning. No catalog
status or GitHub contribution is changed by this local proof document.
