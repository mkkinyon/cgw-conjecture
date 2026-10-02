# STATUS — CGW Conjecture project (rolling; replaces per-session handoffs)

*Updated 2026-10-01 (session 13, continuous).  The notes
`second_row_notes.tex` (105 pp, 0 errors) are the living document; this file
is the two-page map.  History: `git log`.  Old handoffs `HANDOFF_session_*.md`
are kept for the record; their TODO lists are superseded by this file.*

## ALERT (2026-10-02): gap in CGW08 Lemma 3.12 / Theorem 3.13 (rem:cgwgap)
Splitting case 3 of CGW's proof (cross-switch at the neighbour {ωj,ωj′}, then backflip
at {j,j′}) lands in the WRONG type when min(α,β) ≥ 3: the cross-switch reverses the
ω-arc from ωj to ωj′, which contains j′, so {j,j′} ends at distance 2 (forced by the
proof of their Lemma 3.10; verified independently, s13_cgwcase3.py: 2676/2676 wrong at
n=9,12).  Consequence: CGW's 3/2 direction (⟺ P_X[B] ≤ 2/3 via thm:marked — a weak
form of OUR target) is unproved as published; the 1/2 direction is trivial given
thm:marked.  Conditional on "CGW(3/2)" from now on: eq:cgwA's P_X[A] ≥ 1/3,
prop:shortcycles lower bounds, the unconditional witness bound, unconditional (CP)
c = 2/5, eq:CPcols "≥ 1/3", CGW's P[(n)] ≤ 2n^{−2/3} and Cor 4.5, hence the TAIL BOUND in
thm:reduction — DONE: prop:tailsurvive (refereed) gives E[N_ℓ] ≤ 4/ℓ (ℓ ≤ n/2),
E[(κ)_k] ≤ 2^{k−1}(4+k²/log n)(log n)^k and P[κ ≥ 16 log n] = O(n^{−2} log n) from |J| = 2|X_A|
alone, so thm:reduction is now independent of CGW's unproved direction.  Also: CGW(3/2) ⟸
P_X[Flip=∅] ≤ 1/3 (thm:ladder; sharp at n=5), i.e. our (L) in constant form would repair their theorem.
Everything exact of ours (thm:marked, thm:ladder, offsets, lem:rowfair, lem:legalturn,
prop:orbit, eq:Jp) is unaffected.  No published erratum found (quick search only).

## Goal and chain — REWRITTEN after the ladder identity (§sec:ladder)

CGW weak conjecture d_TV(P_n,Q_n) → 0  ⇐  EQ(o(1/log n)) [thm:reduction]
⇐ |P_X(B) − ½| = o(1/log n) [thm:marked]
⇐ **(L): P_X[no ladder pair flippable] = o(1/log n)** [thm:ladder, cor:ladder].

**thm:ladder (PROVED, exact, verified to all digits at n=7):** with
K₀ = min(|C₁|,|C₂|) (apart) / min(d₁₂,d₂₁) (together) and the ladder pairs
(x_k,y_k) = (π⁻ᵏ(1), π⁻ᵏ(2)), 1 ≤ k < K₀, the j′-trades T_k at the ladder
pairs are commuting involutions of X preserving the ladder; each toggles the
bit iff its pair is flippable (j, j′ in different ρ_{x_k,y_k}-cycles).  Hence
on every orbit of ⟨T_k⟩ the bit is EXACTLY fair unless no ladder pair is
flippable, and
    P[B] − P[A] = P[B, Flip=∅] − P[A, Flip=∅],   |P[B] − ½| ≤ ½ P[Flip=∅],
for X and for every fixed R (quenched).  Data: P[Flip=∅ | K₀] = 2^{−(K₀−1)}
to the resolution of 9000/7200 instances at n=30/50; P[Flip=∅] ≈ 4/n (0.125, 0.079, 0.041 at n=30,50,100).

**(L) is the one open item, and it is now reduced to the witness lemma
(hyp:witness, prop:Lwitness).**  Offset ladders D_c = {(π^{−(c+i)}(1), π^{−i}(2))}
(thm:ladderoffset, exact at n=7 for c=0..3) price the short-arc B-states:
averaging the offset identities over c gives P[B, d₁₂=ℓ, |C₁| ≥ n/2+ℓ, θ ≥ θ₀]
≤ 36/((n−2)θ₀) where θ = fraction of offsets whose diagonal has a flippable
pair (eq:offsetavg); the splice-in switching does |C₁| < n/2+ℓ.  So:
- (L) ⇐ witness lemma: for the ladder (≥ n/log²n disjoint pairs) and for the
  diagonal blocks of a short arc, the number of blocks containing a pair whose
  ρ-cycle through j has length ≤ ℓ′ = log⁵n and avoids j′ is ≥ κ·Σ min(½, ℓ′|B|/n)
  except w.p. o(1/log n) (annealed).  Budget: ℓ₀ = n/log²n; total o(1/log² n).
- Data (witness50.log): witness freq ≈ (ℓ′−1)/(n−1) − ℓ′/n; disjoint pairs
  E[WW′]/p² = 1.00±0.03; Var(G) and Var(star count) = binomial within 5 %.
- **CGW08 Thm 3.13 IS P_X[A] ≥ 1/3** (via thm:marked), and summed over types it
  gives the UNCONDITIONAL witness lower bound with no loss in ℓ:
  (2/3)E[(ℓ+β)N_{ℓ+β}] ≤ E[ℓN_ℓ·βN_β] ≤ 2E[(ℓ+β)N_{ℓ+β}], hence E[N_ℓ] ∈ [2/3, 4]/ℓ
  and P[|Q_q| = ℓ] ∈ [0.66, 4]/n for ℓ ≤ √n/40 (prop:shortcycles, refereed;
  exact at n ≤ 11 from CGW's tables; sharp at n=5). Remaining: the CONDITIONAL
  form.  **The 'path version' of CGW Lemma 3.12 is NOT routine (§sec:crossspace):**
  with columns p,p′ prescribed, rows x,y of a ladder pair have two PATHS in the
  free part, and the cross-path pairs (the only ones whose flip exchanges
  crossed/parallel) cannot be repaired when they are B-pairs (switch would trade
  through p; cross-switch arc passes through p′).  CGW's condition (3) is exactly
  what excludes paths — essential, not cosmetic.  What survives in the cross
  space X (rows 1,2 + columns p,p′): (F) flip identity at legal A cross-path
  pairs; (I) INTERCALATE identity E[I_sep; crossed] = E[I_merge; parallel]
  (intercalates on rows x,z, columns q∈P₁, c∈P₂; always legal, frame-preserving);
  data: I = 0.98·|A||B|/n in every bin (icsw50.log).  (I) is FAIR and gives no
  inequality by itself: the deterministic bound is only I_merge ≤ (ℓ_p−1)(ℓ_p′−1)
  (referee: Z₈×Z₂ counterexample to the earlier min(ℓ_p,ℓ_p′)−1 claim); a constant
  needs two-sided intercalate control in the cross space (i)+(ii)+(iii).
  **CORE PROBLEM (CP), eq:CP:** given rows 1,2 (type+mark) and columns p,p′, for
  rows x,y with four distinct symbols on {p,p′}: P[p,p′ in different
  ρ_{x,y}-cycles] ≥ c.  Unconditional (CP): c = 2/5 (1/3 ≤ P[W] ≤ 3/5+4/(5(n−1)),
  refereed; exact at n=5).  **Columns-only (CP) is CGW transposed**
  (eq:CPcols: = 1/2 apart, ≥ 1/3 together) — the obstacle is the conditioning
  on rows 1,2 (type+mark = X itself), NOT the frame.  thm:ladder = the legal
  half of transposed CGW; the illegal half (switch/cross-switch at {x_k,y_k} in
  Lᵀ) = flip / cross-switch of rows 1,2 at (p,p′): toggles ALL side-1×side-2
  statuses but reverses one side (ladder ↦ anti-diagonal).
- **Orbit method (§sec:orbit, new; PROVED parts: lem:rowfair, lem:legalturn,
  prop:orbit, prop:Lorbit):** in X every row trade of rows ∉{1,2} is legal
  (⇒ exact fairness of A/B status of (q,c) for generic rows x,y when q,c lie in
  different ρ_{x,y}-cycles: P_X = 1/2, lem:rowfair); a column-cycle turn is legal
  iff the cycle meets {1,2} in 0 or 2 rows; turns at a fixed matching M of free
  columns commute.  "Clean coordinates" of a ladder pair (x,y): short cycles
  through x avoiding y, the family, and meeting {1,2} in 0/2.  They generate a
  free (Z/2)-action; conditionally on the orbit the bits of a family of disjoint
  pairs are INDEPENDENT and P[(x,y) crossed | orbit] ≤ 1 − ½ȳ_{x,y}, ȳ = orbit
  average of the fraction of clean coordinates separated by {p,p′}.  Hence
  P[Flip_G = ∅] ≤ E[exp(−½ Σ ȳ)] (prop:orbit), and (L) ⇐ hyp:orbit
  (Σ_{sub-ladder} ȳ ≥ κ|I| w.p. 1−o(1/log n), + block form for offsets; the
  flippability concentration of hyp:witness is no longer needed — Chernoff given
  the orbit).  Data (orbit{50,100}[_fixed].log): column-cycle lengths through x
  uniform; P[separated | clean, ℓ] = 0.19–0.20 = P[separated] (genericity);
  orbit coin 0.478/0.523 at n=100 (greedy matching, |T|≈42); fixed matching:
  |T| ≈ 3, E[ȳ] = 0.19 ⇒ per-pair factor ≈ 0.9, κ ≈ 0.15.
  Residual = genericity (β) + macroscopic arcs (α) + concentration (γ).
  (Π) SOLVED (§sec:orbit (f)): P_η[b=1] ≥ ½·max_Z P[Z crossing]; b ≡ 0 on the cube iff
  P₁,P₂ are in different components of the block graph (vertices = arcs P₁,P₂ and the
  other cycles, edges = chords of T); one coordinate per pair suffices (eq:onecoord):
  P_X[Flip_I=∅] ≤ E[Π(1−½X_k)], X_k = [first clean coordinate of pair k crossing].
  (γ): data (cov{50,100}.log) — orbit weights of different pairs UNCORRELATED
  (Var(Y)/(|I|Var ȳ) = 0.99/1.14, excess explained by Var|I|).
  (α) ⇐ (BL) [interleaving of {1,2},{x,y} in column cycles ≤ 1−c; data 0.19].
  (U) in its simplest instances (§(g), refereed): generic columns ⇒ (U) IS the two-arc
  form of (α) (both arcs/cycles at p,p′ macroscopic); generic rows ⇒ (BL) is TRIVIAL
  (interleaved fraction over row pairs ≤ ½ + O(1/n) deterministically; exchangeability)
  but is lost under the conditioning |Q_p(ρ_{x,y})| = ℓ / frame selection, which is the
  content.  DIRECTION: lower bounds on different-cycles events not conditioned on their
  own toggle ((β), third-row availability, CGW's P_X[A] ≥ 1/3) are gap-direction; upper
  bounds are surviving-direction.  Third-row repairs exist (crossed ladder pair: cross-
  switch of x_k,y_k at (p,p′) always defined, lands in J_B with (d₁₂,d₂₁) → (2k,|C|−2k);
  A-instance: backflip through C₂ → J_A) but have unbounded multiplicity / don't preserve
  the ladder — the precise reason switchings in X cannot give the gap direction.
  ONE FORM: every input is "two pairs of lines interleaved in a third pair" (U):
  c ≤ P[interleaved] ≤ 1−c in X, plus the same-cycle event for generic rows at the
  prescribed columns (availability of the (BL) repairs).  Unconditionally (no type
  conditioning) (U) is EXACT given the cycle type: P = A/3 + 2B ∈ [1/27, 1/3+O(1/n)]
  (conjugation invariance; referee).  All exact identities in X are fair; prop:orbit's
  gain is that fairness compounds over a family.
- Literature (scout, §sec:witnessstatus): the conditional form is not in print;
  Allsop–Morris 2026 give (δ/n)^{|P|} ≤ P[L ⊇ P] ≤ (Δ/n)^{|P|} (exp loss in ℓ;
  δ=Δ=1+o(1) open — that IS the lower-bound clause); KS18's exact ratio method
  (1+O(1/n))/(s+1) for intercalates in two rows is the model to follow.
  Needed: two-row cycle switchings with degrees exact to 1+O(polylog/n), in the
  space of completions of rows 1,2 and columns j,j′ (never touched by the moves).

Superseded by (L) (keep as results, do not work on): (I), (II)/(T), (III),
(F1), (F2), the component question, the pool/master-formula route, the
spectral stationary model (still a standalone theorem: paper/spectral_same_cycle.tex).

## Older map (pre-ladder), kept for orientation

### (T), i.e. (II) the transfer — what it contained

(T) is one statement but it carries three things: (a) the real frame law
π_P(L) under Unif(X) must be close to uniform ON THE TEST FUNCTION
ε·P^Mε (a function of the cycle type and the cycles of u,v) — the quenched
one-point content of the old item (I) is relocated here, not removed;
(b) the pool must be in generic position relative to the frame; (c) R/M
row-sharing rare. Status:
- **Trapping half of (a) PROVED** (prop:trapped, refereed): P_X[row 1 in a
  frame cycle of length ℓ] ≤ 18/(n−ℓ−1), by a Kwan–Sudakov switching whose
  only row-1 move is the CGW switch at the mark (priced by thm:marked).
  Trapped marks cost O(1/k'). The rest of (a) is OPEN.
- (b) OPEN. Exact structure (§s13transfer2): delete columns j,j′; frame =
  leftover graph of the n×(n−2) rectangle; pool = ¼-coin thinning of a
  candidate graph; only local exclusion = π-adjacent rows. Real squares
  (true marks) decorrelate 10–30% FASTER than the model at n=50, half that
  at n=100 (the adjacency exclusion), frame statistics agree.
- **Trade-chain route** (§s13tradechain, §s13comp7; standalone statement
  paper/fairness_question.tex, refereed): symmetric chain on S(R) (trades along
  the ρ_{x,y}-cycle through j′, turns of free frame cycles); EQ(4δ) ⇐ bit
  autocorrelation ≤ δ² ⇐ component fairness Σ|K|(2p(K)−1)²/|S(R)| → 0.
  **Exact enumeration at n=7 (s13_comp7.c, runs/s13/comp7/, 7 instances):
  one giant component (99.1–99.7% of S(R)) with p(K) = P[β=1|R] to 4 decimals;
  the rest is "locked" dust (ℤ₇-like squares, every ρ_{x,y} a 7-cycle, nothing
  flippable; components = 5! row permutations), 0.3–0.8%, carrying the whole
  LHS (0.003–0.008). n=6 is degenerate (row parity invariant).**  So the
  component structure is trivial and the Question for fixed R IS the quenched
  statement P[β=1|R] = ½+o(1); the chain is the tool, not a weakening.
  Rectangle picture (delete columns j,j′): frame cycles = leftover row–symbol
  cycles, turns = orientations, bit = function of the rectangle; flippable ⇔
  the two (x,y)-paths pair "parallel"; for x,y in different cycles a turn
  toggles it, so thinning only bites inside C₁∪C₂ (the bit-changing pairs).
  Within-cycle mark ⇔ (1,2) crossed; cross-cycle marks have an exact
  bit-flipping involution (switch at (1,2) + isotopy) ⇒ p = ½ exactly.
  Needs (F1): P_stationary[no flippable bit-changing pair] = o(1/log² n) —
  a parallel/crossed statement about row-pair paths in a random rectangle,
  target for a KS switching as in prop:trapped — and (F2) mixing of the
  thinned walk vs the "coin walk" on D = {derangements, π(1)≠2, π(2)≠1}
  (uniform on D stationary, P_D[1~2] = ½ exactly; fibres {π⁻¹(1)=a,π⁻¹(2)=b}
  = permutations of n−2 points with no fixed point outside {a,b}, exact ½
  there via π ↦ (a b)∘π).  Right-multiplication trades (through j) are
  equally legal and were included as move set B (no qualitative change).

## Bugs found this session
- Sessions 7–8 real-square pipelines (s7/s8_real_bias.py) used σ's SYMBOLS as column
  indices on unnormalised JM squares: their "marks" were random column pairs.
  Qualitative findings survive; §sec:realbias numbers are for random pairs.
  Fixed in s13_real.py / s13_trade.py (symbol → column via row 1).

## Closed / superseded (do not spend time on these)

- (III) prob:decouple: NOT NEEDED (prop:nodecouple; routing merges m ≥ n/(A log n)).
- (I) in all its session 8–12 forms (prob:mixing one-block start, prob:survival,
  (R1), (C1), (C2″), (S₂), occupation-time lemma, Markov renewal, a_k, t_k(w)):
  closed IN THE MODEL; the quenched frame-law question survives inside (T)(a).  The exact laws (thm:blockcount,
  thm:sizebiased, thm:idleloop, lem:inert, prop:nomino) remain valid results.
- Budget: γ > 1/11 not 4/9 (cor:budget) — now irrelevant since c = 1/2 is proved.
- Dust obstruction: a sup-norm statement; irrelevant to L²/stationary (rem:dustblind).
- Still true and worth keeping: lem:flipsrecur (τ_k < ∞ a.s.); the audit of §s13.

## Do not reopen (proved impossible / false)

Doeblin minorisation (prop:nomino); (V-)uniform ergodicity (cor:noerg);
λ-drift (prop:nodrift); geometric (C2); (R2); transport contraction;
two-block odd perturbation; ρ from burnt-in starts (pitfall 56).

## Possible next steps (owner's choice)

-1. (CP) via hyp:orbit (the current best form of the residual).  Done: (a) the
   one-step switching is fair (useless); the orbit method compounds fairness.
   Open inputs (all of the form (U), see above): (α) arc-trapping in X — DONE
   up to (BL): CGW's lemma for the row pair (x,y) runs inside X with the flips
   blocked exactly when rows 1,2 and rows x,y are INTERLEAVED in the (q,c)-column
   cycles (switch-invariant); eq:Jp gives P_X[|Q_p(ρ_{x,y})| = ℓ] ≤ 3/((1−b_ℓ)(n−ℓ))
   with b_ℓ = interleaved fraction (data 0.19, interleave50.log); (BL): b_ℓ ≤ 1−c is
   again a same-cycle genericity for third rows (the recursion, now explicit);
   for ladder pairs the frame must be fixed too (bookkeeping not done); (β) genericity: P[clean coordinate (q,c) separated | ...] ≥ c —
   note the free row trades (lem:rowfair) re-randomise the clean coordinates
   without touching ρ_{x,y}, but their availability ("q,c in different
   ρ_{z,w}-cycles") is again a same-cycle event; (γ) concentration of Σȳ over
   the sub-ladder (second moment: average covariance o(1/log n) suffices).
   (b) relax rows 1,2 to type + mark: DONE in the sense that the orbit method
   works in X (type+mark), with column turns through both rows 1,2 legal.
   (c) exact identity not pinned by a row pair and a column pair: the orbit
   coordinates are exactly such (pinned by row x and column pair (q,c) only).
0. (superseded) (F1)/(F2).

1. (II): formulate the minimal "coarse product" property of the pool and
   the frame-health property, and test both on real squares at n = 100–400
   against the exact prediction; then the TRP/KPS robustness argument.
2. Publishable checkpoint: thm:marked + thm:master + thm:spectral +
   cor:stationarymodel + the layer theorems is a paper now; thm:spectral
   is a standalone result about the random-transposition walk (check
   literature: two-point/same-cycle statistics of random transpositions).
3. Literature check whether thm:spectral is known.

## Numerics this session

`runs/s13/frame*.log` (frame comparison), `akB_0_p.log` (p = 2.4, moot),
`akB_{6,7}_k40.log` (moot); `comp7/*.log` (exact components n=7, ~1 min each);
`trade{30,50,100}.log` (trade chain on real squares); `f1/ladder{30,50,100}.log` (ladder
statistics), `f1/witness50.log` (witness lemma scoping), `f1/f1_*.log` (flippability of bit-changing pairs: fair, independent), `f1/icsw50.log` (intercalate identity), `orbit/orbit{30,50,100}.log`, `orbit/orbit{50,100}_fixed.log`, `orbit/interleave50.log`, `orbit/cov{50,100}.log` (orbit method / blocked fraction / covariance of orbit weights: `./jm 100 120 300000 9 | python3 s13_orbit.py 100 --inst 3 --lmax 8 --fixed`).  `s13_spectral.py brute N` (N ≤ 7) and
`s13_spectral_big.py N` reproduce the theorem.

## Quickstart

HANDOFF_session_12.md §3 QUICKSTART still green, plus:
`python3 s13_spectral.py brute 6` (exact = brute, all digits);
`python3 s13_spectral_big.py 4000` (sum ν = 1, m·E → 1.9);
`gcc -O3 -o s13_frame s13_frame.c -lm && ./s13_frame 2 32 100000` → 4E[bias²] ≈ 0.115;
`pdflatex second_row_notes.tex` ×2 → 0 errors, 123 pp;
`./jm 9 300 2500 3 | python3 s13_cgwcase3.py 9` → CGW case 3 lands in the wrong type (every instance);
`gcc -O2 -DN=7 -o s13_comp7 s13_comp7.c && ./s13_comp7 0123456 1234560 0 2` → ladder identity exact (~1 min);
`./jm 30 300 30000 1 | python3 s13_ladder.py 30` → P[Flip=∅|K₀] ≈ 2^{−(K₀−1)}.
