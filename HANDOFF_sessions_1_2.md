# HANDOFF: Cameron's second-row conjecture — state, tools, and prioritized to-do

*Prepared 2026-08-28 at the end of two working sessions. Written to be fully
self-contained for a fresh session with no prior context. The companion
tarball `second_row_session.tar.gz` contains all code, data, and the working
notes (`second_row_notes.tex` / `.pdf`, 20 pp). Read the notes' abstract and
table of contents first; this file is the operational companion.*

---

## 0. The problem in one paragraph

Let L be a uniformly random Latin square of order n with first row = identity,
and let σ (a derangement) be its second row. Cameron conjectured
(blog, "A niggling problem", 24 Jan 2015) that the law of σ tends *very
rapidly* to uniform on derangements. Write N(σ) = #completions of the 2×n
rectangle (id, σ); N depends only on the cycle type λ of σ; set
C_n(λ) = N(σ)/(n−2)!. P_n = law of the cycle type of σ; Q_n = cycle type law
of a uniform derangement. Weak form (Cavenagh–Greenhill–Wanless [CGW08],
Conj. 6.1): d_TV(P_n, Q_n) = o(1). **Still open** (checked to Sept 2025;
Kwan–Petrova–Sawhney arXiv:2509.13125 resolves the *parity* conjecture from
the same blog post but not this one).

## 1. What the two sessions established (status labels matter)

Notation: M_n = ménage numbers (0, 1, 2, 13, 80, 579, 4738, 43387, 439792,
4890741, 59216642, 775596313 for n = 2..13). (n)_4 = n(n−1)(n−2)(n−3).
B₂ = I ∪ P_σ (forbidden board). N^(k) = # of k×n rectangles with rows 1,2 =
(id, σ). "Layer-k effect" r_k = the relative increase of N^(k)/N^(k−1) when σ
has type (2, n−2) instead of (n).

**PROVED (in the notes, complete or near-complete proofs):**
- Fusion identity M_{C_{2a}}·M_{C_{2b}} = M_{C_{2(a+b)}} + x^{2b}M_{C_{2(a−b)}}
  for cycle matching polynomials — classical (Touchard 1953; Riordan,
  *Introduction to Combinatorial Analysis*, Ch. 8 §3, eqs. (24)–(26),
  Theorem 2). Hence N^(3)((ℓ, n−ℓ)) = M_n + M_{n−2ℓ} exactly, and
  r₃ = M_{n−4}/M_n = (1/(n)_4)(1 − 4/n² + O(n⁻³)). Full asymptotic series
  known: n⁴·r₃ = 1 + 6/n + 21/n² + 42/n³ + ...
- Pair-IE identity (notes Prop. 6.1): N^(4)(σ) = Σ_j (−1)^j Σ_{M∈𝓜_j} R(M;σ)²,
  R(M;σ) = # of B₂-avoiding permutations containing the partial matching M,
  computable by the punctured-board calculus (paths + cycles). Verified
  exactly at n = 6, 7.
- The exact switching identity (notes Thm. 6.1 old numbering / §7 in v3):
  for a split of a part α+β into α ≠ β,
  |S(λ)|/|S(μ)| = (γ(λ)/γ(μ))·(1+q̄)/2, i.e. q̄ = 2C(λ)/C(μ) − 1,
  where q̄ = P(B₁) + P(B₂) for two explicit "post-flip B-pair" events.
  Verified by COMPLETE enumeration at n=5 (q̄=2) and n=6 (q̄=10/11 for
  (2,4)→(6); q̄=4/7 for (2,2,2)→(2,4), the α=β=2 case with the coincident
  event counted twice). P(B₁) = P(B₂) exactly (row-swap involution; verified
  instance-by-instance on 7000+ instances).
- Reduction theorem (notes §8): the "B-pair equidistribution" hypothesis
  EQ(δ): |q̄ − 1| ≤ δ(n) for all α≠β adjacent pairs, with δ = o(1/log n),
  implies d_TV(P_n,Q_n) = O(δ log n + n^{−2/3}), hence the weak conjecture,
  hence (with Łuczak–Pyber + Häggkvist–Janssen) almost all loops have
  multiplication group S_n. The routing avoids α=β splits except for
  all-equal-part types, absorbed in tails via CGW's unconditional bounds.
- Three-row joint Poisson theorem (notes Thm. §7 v3): in a uniform 3×n Latin
  rectangle the three pairwise intercalate counts (X₁₂, X₁₃, X₂₃) are
  exchangeable and jointly → independent Poisson(1/2), with O(1/n) rates on
  factorial moments. *Proof currently a careful sketch — see TODO 2.*

**VERIFIED-EXACT (numerics with exact rational arithmetic, no proofs):**
- log C_n(λ) is additive over cycles: at n=11, max residual 1.8e−6 against a
  leading coefficient h(2) = 9.6e−4. No sign(σ) term (< residual).
- h_n(2)·n³ ≈ 1.28–1.33 at n = 10, 11; d_TV(P_n,Q_n)·n³ = 0.402, 0.388 at
  n = 10, 11, matching the prediction (c₂/2)e^{−1/2}·n^{−3} with c₂ ≈ 1.30.
  CAVEAT (notes Remark 5.2): the limit c₂ is NOT pinned by n ≤ 11 — the
  early-layer weights track 1/(n)_4 = n⁻⁴(1+6/n+…), inflating finite-n
  values; plausible range for the limit: [0.8, 1.5].
- Layer-4 effect r₄: exact for n ≤ 11 (from summing ~5×10⁶ Ryser permanents),
  series-estimated to n = 18: r₄·(n)_4 = 0.843 → 0.984 monotone; r₄/r₃ − 1
  decays 18.8% (n=8) → 0.2% (n=18). Strong evidence r₄ = r₃(1 + O(1/n)).
- t₁ = E[#common cells of two independent B₂-avoiding rows] satisfies
  t₁ − n/(n−2) = n(n−2)Var_c(p_c) ≈ 10⁻⁶..10⁻⁷ (n=8..12): the third row's
  cell probabilities are uniform to relative O(n⁻³); exact n=20 profile flat
  to 1e−8 except a single-cell deficit ≈ −1.4 n⁻³ at ω-distance 2
  (τ(i) = σ²(i) disfavored).
- Tilt identity: E[X₁₂] − 1/2 = (M_{n−4}/(2M_n))(1 + O(n⁻²)) — agrees to
  1e−4 relative by n = 13. Cov(X₁₂, X₁₃)·n² = 0.754 → 0.633 (n=8..13),
  extrapolating to ≈ 0.48 + 2.0/n, i.e. Cov ≈ 1/(2n²).
- Exact q̄ tables (from CGW's exact C_n(λ), n ≤ 11): deficits ½(1−q̄)n³
  stratify by the smallest part destroyed: ≈ 1.28 (2-cycles), ≈ 0.02
  (3-cycles), ≈ ±0.003 (≥4). JM Monte Carlo reproduces q̄ at n = 8 to 0.2σ
  (8σ detection of q̄ < 1) and is consistent at n = 11, 12, 16.
- Parity: P(σ even) − 1/2 matches (−1)^{n−1}(n−1)/(2D_n) to ~5% by n = 11.

**KNOWN ISSUES / UNRESOLVED BOOKKEEPING:**
- α=β≥3 joining bookkeeping: our reading of CGW's joining steps gives 5/6 at
  n=6, μ=(3,3), where the identity requires 3/4. Not used anywhere (routing
  avoids it), but unresolved. See notes Remark on α=β.
- The involution lemma (P(B₁)=P(B₂)) lacks a written algebraic proof
  (empirically instance-exact).

## 2. File manifest (tarball) and re-verification commands

Everything is Python 3 (needs `pip install mpmath sympy --break-system-packages`)
plus gcc. All "verify" commands should print PASS/MATCH lines.

| file | what it is | re-verify with |
|---|---|---|
| `second_row_notes.tex/.pdf` | the working notes (v3, 20 pp) | `pdflatex` ×2 (needs xcolor; use `microtype` with `expansion=false` if fonts complain) |
| `cgw_data.py` | exact C_n(λ), 4≤n≤11, from CGW Tables 1–2 | `python3 cgw_data.py` (checks vs published fractions AND L_n) |
| `fit_additive.py` | additive fits, chain increments, dTV | `python3 fit_additive.py` |
| `third_row.py` | N^(3) via matching polys; fusion checks | `python3 third_row.py` (brute force n=6,7; fusion to n=25) |
| `third_row_asymp.py` | exact ratios to n=400 + Richardson | `python3 third_row_asymp.py` (~5 min) |
| `layer4.c` | exact N^(4) for n≤11 (Ryser–Gray sums) | `gcc -O3 -o layer4 layer4.c && ./layer4 11 0` (~2 min) |
| `verify_identity.py` | complete-enumeration check of the switching identity, n=5,6 | `python3 verify_identity.py 5` (instant), `6a`/`6b`/`6c` (~10 min each) |
| `jm.c`, `jm_stats.py` | Jacobson–Matthews sampler + measurements | `./jm 8 2000000 512 42 \| python3 jm_stats.py 8` |
| `puncture.py` | punctured-board calculus | `python3 puncture.py` (brute force n=6,7) |
| `pair_ie.py` | pair-IE identity + T_j sums | `python3 pair_ie.py verify` (~8 min) |
| `layer4_production.py` | fast T_j (j≤3), truncated layer-4 effects | `python3 layer4_production.py 8 9 10 11` (validates vs EXACT dict) |
| `layer4_big.txt` | output for n=14..18 | — |
| `three_row_moments.py` | exact 3-row joint moments | `python3 three_row_moments.py brute7` (certification), then `8 9 10 ...` |

Key data constants for spot-checks: N^(3)((n)) = M_n; N^(4) exact:
n=8: (M=4738, N4=5955474) type (n); (4740, 5960976) type (2,n−2);
n=11: (4890741, 7058091360392) and (4891320, 7059777493824).

## 3. PRIORITIZED TO-DO

### TODO 1 (first priority): prove the layer-4 theorem r₄ = r₃·(1 + O(1/n))

Goal: a rigorous proof of `r₄(n) = (M_{n−4}/M_n)(1 + O(1/n))`, i.e.
`Δ_λ log P(X = 0) = O(n⁻⁵)` where X = #common cells of two independent
B₂-avoiding rows and Δ_λ is the difference between σ-types (2,n−2) and (n).

Route (all pieces exist in embryo):
1. Prove the **uniform punctured estimate**: for any σ ∈ D_n and partial
   matching M with |M| = j ≤ √n avoiding B₂,
   R(M;σ) = (n−j)!·e^{−E/(n−j)}·(1 + O((j+1)²/n)) with E = #cells of the
   punctured board (E = 2n − O(j)). Method: inclusion–exclusion
   R = Σ(−1)^k m_k (N−k)! with Bonferroni truncation; m_k bounds for
   max-degree-2 boards. This is TODO 2's lemma as well — do it once, well.
2. Prove **cell-probability uniformity**: max_c |p_c(n−2) − 1| = O(n⁻³)
   and the profile structure (deviation at ω-distance d is O(n^{−(2d−1)})).
   Method: path fusion — the punctured board for cell (i,c) at distance d
   splits the big cycle into paths of lengths (2d−3, 2n−2d−1); compare
   different splits via the path analogue of the fusion identity
   (p_a·p_b − p_{a'}·p_{b'} identities; p = Chebyshev-U). This also gives the
   exact constant of the d=2 deficit (≈ −1.4 n⁻³ empirically) — a nice
   standalone lemma.
3. Bound Δ_λ of the coincidence-cumulants: Δt₁ = O(n⁻⁶) follows from step 2;
   for t_j, j ≥ 2, show Δt_j = O(t_j · n⁻³ᐨᵉᵖˢ) and sum. The empirical
   decomposition (in `layer4_production.py` output: "Dlog Phi_J" lines) shows
   where the mass sits; the tail (j ≥ 4) matters, so a uniform-in-j argument
   is needed — suggestion: couple the two σ-types' coincidence processes via
   the near-identical cell-probability profiles, or use a log-derivative
   (interpolation) representation in the tilt parameter.
Acceptance: a theorem + proof in the notes replacing the "overwhelming
evidence" paragraph of §6; numerics already in place to confirm constants.
Difficulty: medium. This is the most tractable genuinely-new theorem.

### TODO 2: finish the three-row Poisson theorem write-up

Goal: turn the proof sketch of the joint Poisson(1/2)³ theorem (notes §7)
into a complete proof. Needs: (a) the uniform punctured estimate of TODO 1
step 1; (b) N^(3)(σ) = n!e⁻²(1+O(1/n)) uniformly over types (same lemma with
M = ∅); (c) careful bookkeeping of degenerate plantings (shared indices,
τ-plantings colliding with σ-plantings or with B₂) — each is a relative
O(1/n), count them explicitly; (d) the exchangeability lemma (done).
Also upgrade Propositions "tilt" and "covariance" from sketches to proofs:
- tilt: uses the exact fusion expansion of N^(3)(λ) (the second-order terms
  are exactly computable: N^(3)((2,2,n−4)) = M_n + 2M_{n−4} + M_{n−8}, etc.)
  plus exact-Poisson factorial moments of derangement cycle counts
  (errors exponentially small).
- covariance: the exact count of admissible τ-swap pairs is
  C(n,2) − n + λ₂(σ); make the "mean swap probability" step rigorous with
  the punctured estimate.
Acceptance: proofs with explicit O(·) constants or clean o(1)'s; this is a
self-contained publishable note (check novelty first: search for "intercalates
in Latin rectangles Poisson" — we found nothing, and Riordan's K(3,n) theory
is aggregate-only).
Difficulty: low-medium. Good warm-up task.

### TODO 3: the triple calculus (layer 5)

Goal: generalize the pair-IE identity to three added rows and compute the
layer-5 effect. The identity to derive and verify:
N^(5)(σ) = # of triples (τ, η, θ) of B₂-avoiding rows, pairwise cell-disjoint
= Σ over triples of pairwise-disjoint partial matchings (M_{τη}, M_{τθ}, M_{ηθ})
of (−1)^{|M_{τη}|+|M_{τθ}|+|M_{ηθ}|} · R(M_{τη} ∪ M_{τθ}) · R(M_{τη} ∪ M_{ηθ})
· R(M_{τθ} ∪ M_{ηθ}) — CHECK THIS CAREFULLY: the three coincidence patterns
must be compatible (the union of the two matchings meeting a given row must
itself be a matching), and the IE is over three independent "bad event"
families. Verify exactly at n = 6 against brute force (enumerate triples;
~80³ ≈ 5×10⁵ triples at n=6 — trivial). Then truncate at total order ≤ 3 and
compute the layer-5 effect for types (n) vs (2,n−2) at n = 8..14; compare
with the "aggregate remaining layers" numbers (notes §5: at n=11 the layers
5..11 average 1.51·n⁻⁴ against 1.73–1.76 for layers 3–4). The layer-5 point
is the first direct probe of the profile γ beyond γ(0).
Acceptance: identity verified exactly at n=6; layer-5 effect at n=8
cross-checked against a direct computation (extend `layer4.c` to N^(5) at
n=8: sum per(J−B₂−τ−η) over pairs — ~4738·5955474/4738 ≈ 6×10⁶ permanents of
order 8 — a few minutes in C; or do n=7 exactly).
Difficulty: medium-high (bookkeeping-heavy; high value for the profile).

### TODO 4: the bulk layer profile γ and the constant c₂

Goal: determine γ(t) = lim n⁴·(layer-⌊tn⌋ effect), whence c₂ = ∫₀¹γ.
This is the hardest computational-analytic target. Three attacks:
(a) Dense-regime cluster expansion: for the board B_k = I∪σ∪τ₃∪…∪τ_k
   (k-regular), per(J−B_k) = ∫₀^∞ e^{−t}t^n M_{B_k}(−1/t)dt exactly; the
   λ-sensitive term is the 4-cycle count #C₄(B_k) ⊇ λ₂(σ). Heilmann–Lieb
   zero-freeness gives a convergent cluster expansion of log M_B in the
   relevant region for k ≤ (1−ε)n. Work out the C₄-cluster weight at
   saddle t*(k, n) — this predicts γ(t). Sanity anchors: γ(0) = 1 (proved);
   the empirical layer-4 ≈ layer-3.
(b) Transfer-matrix experiments: circulant boards {0,1,2,…,k−1} with a
   planted defect creating one extra C₄ — permanents of banded circulant
   complements are computable by transfer matrix for n in the hundreds
   (cf. Metropolis–Stein–Stein, "Permanents of cyclic (0,1) matrices");
   measure the per-C₄ boost as a function of band width k. This directly
   measures the cluster weight in a clean ensemble.
(c) Exact layer-5 (TODO 3) + aggregate constraints from exact totals n ≤ 11.
Acceptance: any principled determination of γ beyond t=0, even heuristic
with numerical confirmation; update Remark 5.2 and Conjecture 3.2's constant.
Difficulty: high. Don't start here; do after TODO 1+3.

### TODO 5: attack EQ (B-pair equidistribution) with the KPS machinery

The strategic prize: EQ(o(1/log n)) proves the weak conjecture (notes §8).
This needs its own session. Concretely:
1. Fetch arXiv:2509.13125 (Kwan–Petrova–Sawhney, "Parities in random Latin
   squares") IN FULL (ask the project owner for the PDF if fetching is
   blocked) and absorb: their stable intercalate switchings (§ on
   re-randomisation) and the TRP comparison theorem (their Lemma 4.7).
2. Key observation to exploit (notes §9): an intercalate switch involving
   row 2 and another row composes σ₁,₂ with a transposition = a degenerate
   CGW flip. Their parity theorem shows enough well-separated stable
   intercalates exist to re-randomise a ±1 functional; the target here is
   the *same-cycle indicator* of the two marked rows in a column pair —
   one bit, but of a conditioned structure.
3. Precise target statement: for uniform L ∈ S(λ) and the post-flip column
   pair of the identity, P(B-pair) = 1/2 + o(1/log n). Try: exhibit a
   switching within S(λ) that flips the B-status of the marked pair while
   preserving the measure up to controllable factors (the cross-switch of
   CGW does exactly this — Lemma 3.10 toggles A/B of specific pairs! — so
   the question becomes equidistribution under cross-switch orbits).
   That reframing (cross-switch orbit mixing) may be the cleanest attack
   and doesn't need KPS at all — CHECK IT FIRST: the cross-switch is an
   involution on S(λ) that toggles the B-status of the relevant pair;
   if one can show it is measure-preserving up to 1+o(1/log n) ON AVERAGE
   over instances, EQ follows. (It is exactly measure-preserving as a map;
   the issue is it moves the marked pair's environment — quantify.)
Acceptance: either a proof of EQ(δ) for any δ = o(1), or a precise
identification of the obstruction.
Difficulty: high/research-grade. Highest payoff.

### TODO 6 (cleanups, any session, low effort)

a. Write the algebraic proof of the involution lemma P(B₁)=P(B₂)
   (chase ι = row-swap through switch and flip; empirics guarantee it).
b. Exact constant of the d=2 cell deficit (path-fusion computation; predicts
   the t₁ deviation constants). Falls out of TODO 1 step 2.
c. Resolve the α=β≥3 exact bookkeeping (n=6 μ=(3,3) discrepancy 5/6 vs 3/4).
   Determine the correct joining multiplicities; `verify_identity.py 6b` is
   the test harness — try variants of the "extra edge" rules until the exact
   rational 3/4 is reproduced, then re-derive the rule from the multigraph.
d. Third-order intercalate expectation: prove
   E[#intercalates] = n²/4 − n/4 + Θ(1/n) (ties c₂ to McKay–Wanless-style
   switchings; possibly a standalone paper with TODO 2).
e. Sharp parity conjecture: P(σ even) − 1/2 ~ (−1)^{n−1}(n−1)/(2D_n)
   (KPS-style techniques; harder than it looks).
f. Consider splitting the notes into two papers: (i) "The second row of a
   random Latin square" (conjectures + identity + reduction + numerics);
   (ii) "Intercalates in three-line Latin rectangles" (TODO 2 material).

## 4. Pitfalls (hard-won; do not rediscover)

1. **JM sampler bias**: never emit "the first proper state after T moves" —
   it biases against intercalate-rich squares (at n=4: P((2,2)) = 0.389
   instead of 1/2!). Sample on a fixed time grid and DISCARD improper
   checkpoints (current `jm.c` does this correctly). Validate any sampler
   change against exact cycle-type tables at n=8 and n=11 (χ² should be
   ~dof; biased versions give χ² ≈ 500–700 on 6 dof).
2. **Truncation systematics**: the pair-IE series estimates (J=2, J=3
   average) carry a validated systematic of ≈ −0.4% on the layer-4 effect.
   Signals smaller than that (e.g. r₄/r₃ − 1 for n ≥ 14) are NOT resolved.
3. **Finite-size inflation**: n³h_n(2) ≈ 1.30 at n ≤ 11 does NOT determine
   lim c₂ (Remark 5.2). Do not "confirm c₂ = 1.30" from small n.
4. **α=β≥3**: do not use the exact identity there (see TODO 6c); the
   reduction theorem's routing avoids it — keep it that way unless 6c is
   resolved.
5. Long background runs in the sandbox: launch with
   `setsid nohup bash -c '...' >/dev/null 2>&1 < /dev/null &` and poll files;
   plain `nohup ... &` from the harness shell gets killed. `pkill -f` with a
   pattern matching your own launcher kills the launcher — match exactly
   (e.g. `pkill -f "^\./jm "`).
6. LaTeX in the sandbox: no `lmodern`; use
   `\usepackage[protrusion=true,expansion=false]{microtype}` and `xcolor`.

## 5. Papers on hand / to request

In the original literature tarball (session 1): CGW 2008 (RSA 33, 286–309)
— the central prior work, exact tables source; Häggkvist–Janssen 1996;
Cameron 1992. Uploaded later: Riordan's *Introduction to Combinatorial
Analysis* (full book PDF; Ch. 8 §3–4 are the relevant pages 198–207).
To request from the project owner if needed: Touchard, Scripta Math. 19
(1953) 109–119; Riordan, Scripta Math. 20 (1954) 14–23 (both hard to find;
the book substitutes for most purposes); Kwan–Petrova–Sawhney
arXiv:2509.13125 full PDF (for TODO 5); Moser, Canad. J. Math. 19 (1967)
1011–1017 and Whitehead, J. Austral. Math. Soc. A 28 (1979) 369–377
(for TODO 3/4 context); McKay–Wanless, JCTA 86 (1999) (for TODO 6d).

## 6. Quick-reference numbers

- d_TV(P_n,Q_n): 2.01e−3 (n=7), 1.21e−3, 4.67e−4, 4.02e−4, 2.91e−4 (n=11).
- max/min completion ratio − 1: 1.41e−2 (n=7) … 3.87e−3 (n=11).
- q̄ for (2,n−2)→(n): 0.988372, 0.993345, 0.996822, 0.997363, 0.998080
  (n=7..11); ½(1−q̄)n³: 1.99, 1.70, 1.16, 1.32, 1.28.
- r₃ = M_{n−4}/M_n; r₄ (exact): 5.0153e−4, 3.1340e−4, 1.8715e−4, 1.2049e−4
  (n=8..11); estimates to n=18 in `layer4_big.txt`.
- 3-row moments: E[X_ab]−1/2 = M_{n−4}/(2M_n) to 1e−4 rel. (n≥13);
  n²Cov(X_ab,X_cd) → ≈ 0.48 + 2.0/n; E[C(X,2)] → 0.125.
- Predicted d_TV constant: (c₂/2)e^{−1/2} ≈ 0.394 with c₂ = 1.30 (finite-n).

*End of handoff. The notes (`second_row_notes.tex`) are the living document:
update them, keep every numerical claim backed by a script in the tarball,
and keep the PROVED / VERIFIED / COMPUTED distinction honest.*
