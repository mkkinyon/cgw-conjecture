# HANDOFF: Cameron's second-row conjecture — state after session 4
## >>> TODO 3 IS NOW CLOSED: the layer-5 theorem is proved (same rigor as layer 4), with exact confirmation to n=9. <<<

*Prepared 2026-08-28 (end of the fourth working session). Fully self-contained
for a fresh session. The companion tarball contains all code, data, and the
working notes (`second_row_notes.tex` / `.pdf`, 35 pp). Read the notes'
abstract and ToC first; this file is the operational companion. The previous
handoffs (`HANDOFF_sessions_1_2.md`, `HANDOFF_session_3.md`, kept in the
tarball) remain valid except where superseded below.*

---

## 0. The problem in one paragraph

Let L be a uniformly random Latin square of order n with first row = identity,
and let σ (a derangement) be its second row. Cameron conjectured (blog, "A
niggling problem", 24 Jan 2015) that the law of σ tends *very rapidly* to
uniform on derangements. N(σ) = #completions of the 2×n rectangle (id, σ);
C_n(λ) = N(σ)/(n−2)!; P_n = law of the cycle type of σ; Q_n = that of a
uniform derangement. Weak form (CGW08 Conj 6.1): d_TV(P_n,Q_n) = o(1).
**Still open.** Refined picture (notes §2–3): log C_n additive over cycles,
h_n(2) ≈ c₂n⁻³, d_TV ≍ n⁻³.

## 1. What session 4 added (on top of sessions 1–3)

Session 4 attacked TODO 3 (triple calculus / layer 5) and closed it.
All in the new notes §8 ("The fifth layer: the switch calculus transfers",
label sec:fifth), between the switch-calculus §7 and the three-row section.

- **Triple-IE identity (Prop, proved; verified n=5,6,7 both types):**
  N⁽⁵⁾ = Σ over compatible (M₁,M₂,M₃) of (−1)^{|M₁|+|M₂|+|M₃|}
  R(M₁∪M₂)R(M₁∪M₃)R(M₂∪M₃). KEY RESUMMATION: pairwise compatibility ⇔
  U = M₁∪M₂∪M₃ is a matching; per-cell label weights collapse to
  −1 (each of the three "exactly two constraint sets" patterns) and +2
  (all three); every coincidence cell lies in ≥2 of the constraint sets
  A,B,C. `triple_ie.py`.
- **Exact N⁽⁵⁾ by direct enumeration** (`n5_direct.c` 64-bit, n=8;
  `n5_direct2.c` 128-bit two-word masks + OpenMP, n=9 in ~40 min on
  2 cores; N⁽³⁾,N⁽⁴⁾ reproduced exactly as certification):
  n=8: 1539331584 (type (n)), 1541544192 ((2,n−2));
  n=9: 1417404643008, 1418731458816.
  ⇒ r₅·(n)₄ = 0.8617 (n=8), 0.9759 (n=9); r₅/r₄ = 1.023, 1.030;
  ΔlogN⁽⁵⁾/(3M_{n−4}/M_n) = 1.134, 1.041.
- **Switch decomposition, layer 5 (Lemma, exact; verified n=6,7 brute,
  n=8 vs direct):** ΔN⁽⁵⁾ = 3Br_G⁽⁵⁾ + 3[G₂⁽⁵⁾(e)−G₂⁽⁵⁾(d)] +
  6[H⁽⁵⁾(e)−H⁽⁵⁾(d)] (coefficients 3,3,6 replace layer 4's 2,2,2; same
  pins d₁,d₂,e₁,e₂ on the same board0). At n=8 the pieces (anatomy5.c:
  2649312 + 159840 − 596544) sum to 2212608 = exact ΔN⁽⁵⁾; at n=9
  (anatomy5b.c) 1632687552 + 69975336 − 375847080 = 1326815808 = exact
  ΔN⁽⁵⁾. Ψ-symmetry acts diagonally on triples ⇒ G⁽⁵⁾(d₁)=G⁽⁵⁾(d₂)
  exactly (n=7,8,9: e.g. 357845442448 both at n=9).
- **Stratified bracket (eq:BrG5strata):** Br_G⁽⁵⁾ = Σ_A B_A·w(A) with
  B_A THE SAME four-term layer-4 bracket in the A-slot and w(A) the
  signed spectator weight. Exact identities (verified on board0 at
  n=6,7): (i) Σ_A R₀(A)w(A) = N₀⁽⁵⁾; (ii) w(∅) = N₀⁽⁴⁾,
  so the A=∅ stratum is M_{n−4}N₀⁽⁴⁾ exactly. `switch5.py check3`.
- **Transfer of the strata analysis:** two new observations do all the
  work. (1) Every cell of A lies in B∪C, so line-swap pairings act
  diagonally; pinned-factor equalities (deleted line SETS) are untouched;
  free-factor differences follow by the product rule, one Wronskian per
  factor. (2) Spectator stability: re-indexing the M₃-sum under a pairing
  is the identity off an O(1/n)-measure set (line-measure lemma inside w).
- **Theorem (fifth layer; notes thm:layer5):** ΔlogN⁽⁵⁾ = 3M_{n−4}/M_n +
  O(n⁻⁵), hence r₅ = r₄(1+O(1/n)) = r₃(1+O(1/n)). Same working-note
  rigor as thm:layer4 (referee pass = TODO 1″ below).
- **Remark (rem:layerk):** by induction, Δlog N⁽ᵏ⁾ = (k−2)M_{n−4}/M_n +
  O_k(n⁻⁵) for every fixed k ⇒ the layer profile γ is FLAT AT THE ORIGIN.
- **All types at n=7,8 (rem:alltypes; `n5_types.c`):** exact N⁽³⁾/N⁽⁴⁾/N⁽⁵⁾
  for every derangement type at n≤8. Uniformity of the (k−2)-pattern over
  types; 3-cycles show the predicted M₂=0 near-null (small negative)
  layer-4,5 effects; (4,4) tracks (2,6) at all three layers (three-layer
  confirmation of the M₀=2 convention); (2,3,3) = (2,6)+2·(3,5) to 0.3%
  (layer 4) / 1.3% (layer 5) — additivity over cycles at layers 4,5.
- **Truncated triple-IE series** (`layer5_production.py`, S(Q), Q≤3):
  converges slowly (total coincidence intensity ≈ 6 vs 2 for the pair
  series): effect 24% low at n=8, 6% at n=9 (exact anchors), Q=3 values
  r₅·(n)₄ = 0.922, 0.942, 0.964, 0.970, 0.974 at n=10..14 (each a few %
  low, approaching 1 monotonically). Use only as consistency check.

**Everything from sessions 1–3 still stands.** The layer-4 theorem, the
three-row Poisson package, and all older PROVED/VERIFIED labels are
unchanged; `HANDOFF_session_3.md` §1 remains the reference for those.

## 2. File manifest (additions; everything else as before)

| file | what it is | re-verify with |
|---|---|---|
| `triple_ie.py` | triple-IE identity, both forms | `python3 triple_ie.py verify` (~1 min); `verify7` (~15 min) |
| `n5_direct.c` | exact N⁽⁵⁾, n=8, 64-bit | `gcc -O3 -o n5_direct n5_direct.c && ./n5_direct 8` (1 s) |
| `n5_direct2.c` | exact N⁽⁵⁾, n≤9, 128-bit+OpenMP | `gcc -O3 -fopenmp ...; ./n5_direct2 8` (1 s); `./n5_direct2 9` (~40 min, 2 threads) |
| `n5_types.c` | exact N⁽³⁻⁵⁾ for any type, n≤8 | `./n5_types 8 2 3 3` etc. (instant) |
| `switch5.py` | layer-5 switch decomposition + identities (i),(ii) | `python3 switch5.py 1 2` (~2 min); `python3 -c "import switch5; switch5.check3((6,))"` |
| `anatomy5.c` | exact Br_G⁽⁵⁾, N₀⁽⁵⁾, G₂/H pieces on board0 (reference version) | `gcc -O3 -fopenmp -o anatomy5 anatomy5.c && ./anatomy5 8` (~1 min) |
| `anatomy5b.c` | same, via global disjointness bit-matrix (fast; 378 MB at n=9) | `./anatomy5b 8` (0.5 s, must match anatomy5); `./anatomy5b 9` (~8 min) |
| `layer5_production.py` | truncated triple-IE S(Q), Q≤3 | `python3 layer5_production.py validate` then `8 9 10 11` |
| `HANDOFF_session_3.md` | the previous handoff | — |
| `KPS_2025.pdf` | Kwan–Petrova–Sawhney arXiv:2509.13125v2 (for TODO 5) | — |

Key new exact numbers: see §1 above; plus board0 quantities
N₀ = 791/6203/55000 (n=7/8/9), N₀⁽⁵⁾ = 6632664 / 3934123200 /
3190836569904 (n=7/8/9), Br_G⁽⁵⁾ = 5128 / 883104 / 544229184;
Br_G⁽⁵⁾·N₀/(M_{n−4}N₀⁽⁵⁾) = 0.612, 0.696, 0.722 (n=7,8,9) — approaching
1−c/n with c ≈ 2.5 ≈ 2× the layer-4 constant 1.17 (two spectators).
(`anatomy5b.c` is the fast bit-matrix version: n=9 in ~8 min vs ~10 h
for `anatomy5.c`; identical outputs at n=8.)

## 3. PRIORITIZED TO-DO (updated)

### TODO 1″ (light): referee pass, now covering §7 AND §8

(a) re-derive the case inventories independently (layer 4 as before; for
layer 5 additionally the spectator-stability step — a j-resolved strata
scan for eq:BrG5strata, the analogue of `switch_j1.py`, would be the
safety net and does not yet exist); (b) tighten orphan/crude bounds into
displayed inequalities; (c) the paper split (old TODO 6f) is riper now:
paper (i) could carry layers 3,4,5 as a unified "switch calculus" story.

### TODO 2′ (cleanup): three-row note — unchanged from session 3.

### TODO 4 (profile γ / c₂) — NOW THE MAIN OPEN COMPUTATIONAL TARGET

γ is flat at the origin (rem:layerk). The remaining unknown is the bulk
k ≍ tn. Attacks unchanged from `HANDOFF_sessions_1_2.md` TODO 4 (cluster
expansion at the saddle; banded-circulant transfer-matrix experiments;
aggregate constraints n≤11, now sharpened by exact layer-5 at n≤9). New
leverage: exact all-type data at n≤8 (`n5_types.c`) pins the
ℓ-dependence h⁽ᵏ⁾(ℓ) at k=4,5.

### TODO 5 (EQ / KPS machinery): unchanged — the strategic prize.

**NEW: the KPS paper is now IN THE TARBALL** (`KPS_2025.pdf`,
Kwan–Petrova–Sawhney, "Parities in random Latin squares",
arXiv:2509.13125v2, 17 Sep 2025, 38 pp) — no need to request it.
Orientation from a first skim (verified against v2; the old handoff's
"Lemma 4.7" citation is CORRECT in this version):

- §2.2 "A new approximation lemma": **Lemma 4.7** is the main
  approximation lemma (triangle-removal-process property with prob 1−ε
  ⇒ holds in a random subset of a random Latin square with prob
  1−ε·exp(O(n log²n))); **Corollary 5.5** is the easy-to-apply version
  comparing random Latin squares to Erdős–Rényi random hypergraphs.
- §2.3: cycle switchings (their notation for CGW-type moves; they cite
  the "cross-switch" operations explicitly — our TODO 5 route).
- §2.4: individual intercalate switchings — why single-switch ratio
  estimates (their eq. (2.1), the analogue of our q̄!) reach tails but
  not the bulk; the "infamous upper tail" clustering obstruction.
- §2.5: **multiple intercalate switchings — the re-randomisation**. Key
  device: a stable/canonical collection 𝓘(L) of disjoint intercalates
  with 𝓘(L)=𝓘(L′) under any subset-switch (built via a random sparse
  "template" T of switchable row/column pairs), giving
  P(L→L′)=P(L′→L)=2^{−|𝓘(L)|} exactly; then F₂ linear algebra. The
  collection "robustly spans almost all rows": every set R₀ of ≈log²n
  rows has an intercalate meeting R₀ and its complement. Codimension
  ≈log²n ≪ √n fluctuations (their Remark 2.1 explains why the
  approximation quality matters).
- Main technical theorem: **Theorem 6.4**. Footnote 4: the ε in the
  entry-by-entry enumeration cannot exceed 1/2 — a general barrier for
  TRP-based methods (relevant to how strong an EQ(δ) one can hope for
  from this machinery).

Strategy pointers unchanged (`HANDOFF_sessions_1_2.md` TODO 5): check
the cross-switch orbit reframing FIRST (needs no KPS); if blocked, the
target for KPS-style re-randomisation is the same-cycle indicator of
the marked column pair — one conditioned bit, vs their row-parity bits.
Their stable-collection trick is the plausible transfer: switch a
stable collection of intercalates NOT meeting the marked pair's
environment to re-randomise it within S(λ).

### TODO 6 (cleanups): unchanged from previous handoffs.
(G⁽⁵⁾(d₁)=G⁽⁵⁾(d₂) needs no new work: the diagonal Ψ-action remark in
§8 already is the one-line proof, and it is verified at n=7,8.)

## 4. Pitfalls (all previous ones still apply, plus)

12. **The truncated triple-IE series is NOT the tool for the layer-5
    effect** — its coincidence intensity is 3μ ≈ 6, so Q≤3 partial sums
    are far from converged (effect 24% low at n=8). Trust only exact
    anchors (n≤9) and the switch calculus. If larger-n numerics are
    needed, extend Q or improve the resummation before believing them.
13. **The (k−2) normalisation**: compare ΔlogN⁽ᵏ⁾ against
    (k−2)M_{n−4}/M_n. Sanity: exact ratios 1.000/1.094/1.134 (n=8,
    k=3/4/5) and 1.000/1.023/1.041 (n=9).
14. **n5_direct2 counts ORDERED triples** (×3! vs unordered); N⁽⁵⁾ = # of
    5×n rectangles with rows 1,2 fixed = the ordered count. Keep
    consistent when cross-checking.
15. In `anatomy5.c` board0's rows 2 and n (1-indexed; 1 and n−1 in code)
    keep ONLY the diagonal cell (the σ-cells are the removed switch
    cells); getting this wrong silently shifts N₀ (correct values
    791/6203/55000 at n=7/8/9).

## 5. Papers on hand / to request

As previous handoffs, PLUS: `KPS_2025.pdf` (arXiv:2509.13125v2) is now
in the tarball — the main item previously listed as "to request". Still
outstanding if ever needed: Touchard 1953 / Riordan 1954 Scripta Math.
(the Riordan book substitutes), Moser 1967, Whitehead 1979 (TODO 3/4
context — less urgent now that TODO 3 is closed), McKay–Wanless JCTA
1999 (TODO 6d).

## 6. Quick-reference numbers — previous handoffs plus §§1–2 above.

*The notes are the living document: keep every numerical claim backed by a
script in the tarball, and keep the PROVED / VERIFIED / COMPUTED labels
honest. Session-4 headline: the fifth-layer theorem r₅ = r₃(1+O(1/n))
(by transfer of the switch calculus through a triple-IE), exact N⁽⁵⁾ to
n=9 confirming it to 3%, and all-type tables showing per-cycle
additivity at layers 4 and 5.*
