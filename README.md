# Code and data for "The second row of a random Latin square" working notes

> NOTE (session 8): this README indexes the early-session code only.
> The authoritative, up-to-date index is HANDOFF.md (= HANDOFF_session_8.md)
> and the manifest appendix of second_row_notes.pdf.

Session artifacts, August 2026.

## Exact data
- `cgw_data.py` — exact completion counts C_n(lambda), 4<=n<=11, transcribed from
  Cavenagh–Greenhill–Wanless (RSA 33 (2008), Tables 1–2). Running it verifies the
  data against (a) the published P_n(lambda) fractions and (b) the known totals
  L_n of Latin squares (McKay–Wanless) — both checks pass exactly for all n<=11.

## Analysis of the exact tables
- `fit_additive.py` — additive fits of log C_n(lambda) over cycle counts,
  per-2-cycle chain increments, dTV(P_n,Q_n), max/min spreads.

## Third-row layer (exact theory)
NOTE (attribution): the fusion identity and the reduction of two-permutation
discordant counts to menage numbers turn out to be classical -- Touchard
(Scripta Math. 19, 1953) and Riordan, Introduction to Combinatorial Analysis,
Ch. 8 Sec. 3 (eqs. (24)-(26), Theorem 2). What is new here is the quantitative
per-l-cycle asymptotics and the additivity reading. See Remark 4.4 of the notes.
- `third_row.py` — N^(3) via cycle matching polynomials; brute-force checks
  (n=6,7); fusion identity R3((l,n-l)) = M_n + M_{n-2l} verified for all
  2<=l<=n/2, n<=25; general fusion recursion verified; additive fit at n=30.
- `third_row_asymp.py` — exact ratios to n=400, Richardson extrapolation of
  n^4 rho_2 = 1 + 6/n + 21/n^2 + ... and n^6 rho_3 = 1 + 15/n + 134/n^2 + ...

## Fourth-row layer (exact computation)
- `layer4.c` — N^(4)(sigma) for sigma an n-cycle, (2,n-2)- and (3,n-3)-type,
  n<=11, by summing Ryser–Gray-code permanents over all third rows.
  Results embedded in the notes (Section 5 table).

## Switching identity verification
- `verify_identity.py` — complete enumeration of Latin squares of orders 5, 6
  (first row normalised to id; symbol-relabelling invariance argued in the
  notes), implementing the CGW switch/flip and the joining bookkeeping.
  Exact rational results:
    n=5, (2,3)->(5):    qbar = 2       = 2 C(5)/C(2,3) - 1   (match)
    n=6, (2,4)->(6):    qbar = 10/11   = 2 C(6)/C(2,4) - 1   (match)
    n=6, (2,2,2)->(2,4):qbar = 4/7     = 2 C(2,4)/C(2,2,2)-1 (match, alpha=beta=2)
    n=6, (3,3)->(6):    naive reading gives 5/6 vs required 3/4 (alpha=beta>=3
                        bookkeeping open; see notes, Remark on alpha=beta)
  Also: P(B1) = P(B2) exactly in every case; the row-swap involution exchanges
  the two events instance-by-instance (checked on 7000+ instances).

## Monte Carlo at larger n
- `jm.c` — Jacobson–Matthews sampler, O(1) per move via line-position lists.
  IMPORTANT: sampling must be done on a fixed grid of times, discarding
  improper checkpoints. Emitting "the first proper state after T moves" is
  biased by the improper-excursion flux (correlates with intercalate counts):
  at n=4 it distorts P(type (2,2)) from 1/2 to 0.389.
- `mirror_jm.py` — audited Python mirror of the chain used to localise that bias.
- `jm_stats.py` — measurement pipeline (cycle types vs exact tables; generic-pair
  B probability; post-flip B1/B2 events and qbar for mu=(2,n-2)).
- `run8.txt`, `run11.txt`, `run12.txt`, `run16.txt` — production outputs.

## Notes
The write-up is `../notes/second_row_notes.tex` / `.pdf`.

## Pair calculus (session 2 additions)
- `puncture.py` -- punctured-board Touchard-Riordan calculus (paths+cycles),
  verified against exhaustive enumeration at n=6,7.
- `pair_ie.py` -- the pair inclusion-exclusion identity
  N4(sigma) = sum_j (-1)^j sum_M R(M;sigma)^2, verified exactly at n=6,7.
- `layer4_production.py`, `layer4_big.txt` -- exact T_j (j<=3), truncated
  layer-4 effects to n=18 (validated against exact n<=11 values).
- `three_row_moments.py` -- exact joint intercalate moments of 3xn rectangles
  (certified against brute force at n=7): Poisson(1/2)^3 limit data, the tilt
  identity E[X12]-1/2 = M_{n-4}/(2M_n)(1+O(n^-2)), covariance ~ 1/(2n^2).

## Fifth layer / triple calculus (session 4 additions)
- `triple_ie.py` -- the TRIPLE inclusion-exclusion identity for N5 (two
  equivalent forms: raw M-form and the resummed U-form with per-cell label
  weights -1/-1/-1/+2), verified against brute-force triple enumeration at
  n=5 (both forms) and n=6,7 (both cycle types).
- `n5_direct.c`, `n5_direct2.c` -- exact N5 by direct enumeration (64/128-bit
  cell masks, OpenMP), n<=9, both types; reproduces the known exact N3, N4
  as certification. Key values:
    n=8: N5 = 1539331584 (type (n)),   1541544192 (type (2,n-2))
    n=9: N5 = 1417404643008 (type (n)),1418731458816 (type (2,n-2))
  => r5*(n)_4 = 0.8617 (n=8), 0.9759 (n=9); r5/r4 = 1.023, 1.030.
- `switch5.py` -- the fifth-layer switch decomposition
  DeltaN5 = 3 BrG5 + 3 dG2 + 6 dH (verified by complete enumeration n=6,7),
  the BrG5 leading-term structure, and check3: the U-form evaluated on
  board0 equals P0^(5), and w(empty) = P0^(4) (identities (i),(ii) of the
  notes' stratified bracket).
- `anatomy5.c` -- exact BrG5, N0^(5), G2^(5)-, H^(5)-pieces on board0
  (n<=9; n=9 takes ~1-2 h on 2 threads). At n=8 the pieces sum to 2212608 =
  the direct DeltaN5 (independent cross-check). G5(d1) = G5(d2) exactly
  (Psi-symmetry acts diagonally on triples).
- `layer5_production.py` -- truncated triple-IE series S(Q), Q<=3; converges
  slowly (total coincidence intensity ~6): effect underestimated by 24% at
  n=8, 6% at n=9; use exact anchors instead.
- `anatomy5b.c` -- fast bit-matrix version of anatomy5 (global disjointness
  bitset; n=9 in ~8 min, 378 MB). Verified identical to anatomy5.c at n=8.
  n=9 pieces: 3BrG5+3dG25+6dH5 = 1326815808 = exact DeltaN5 (independent
  cross-check at n=9); G5(d1)=G5(d2)=357845442448 exactly (Psi).
- `n5_types.c` -- exact N3/N4/N5 for arbitrary cycle type, n<=8 (per-type
  layer tables; additivity checks at layers 4,5).

## Marked-pair identity (session 5 additions)
- `pf_orbit.py` -- the CGW cross-switch (Def 3.9 of CGW08) implemented and
  verified exhaustively at n=5 (Latin/type/involution/marked-status checks,
  shifted-pair toggle profile: both shifts k=+-1 toggle deterministically),
  plus splitting-side statistics P(B_i | marked A/B): P(B_i|marked B) = 1/2
  exactly (forced by the cross-switch involution).
- `pf_bridge.py` -- the joining<->splitting bridge, verified
  instance-by-instance at n=5,6: every join lands in X_A (Lemma S), each
  marked-A instance is hit exactly twice (A-join + switched partner), X_B is
  orphaned; J_A = J_B = |X_A| exactly, and qbar = |X_B|/|X_A| for
  (2,3)->(5), (2,4)->(6), (2,2,2)->(2,4), and (3,3)->(6) [= 3/4: resolves
  the alpha=beta>=3 anomaly of rem:albe / TODO 6c].
- `xside.c` -- X-side counts X_A, X_B for ALL types and merges by complete
  enumeration (n<=7): exact match with 2C(lambda)/C(mu)-1 in every case.
- `pf_toggle.py` -- grounding statistics for the re-randomisation programme:
  useful intercalates (completion rows, exactly one marked column) and
  per-switch bit-flip counts, exact at n=5,6.
- `marked_pair_theorem.md` -- working notes for the theorem and programme.
- `section_marked.tex` -- the new notes section (input by the main tex).
