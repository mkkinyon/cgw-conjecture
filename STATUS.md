# STATUS — CGW Conjecture project (rolling; replaces per-session handoffs)

*Updated 2026-10-01 (session 13, continuous).  The notes
`second_row_notes.tex` (105 pp, 0 errors) are the living document; this file
is the two-page map.  History: `git log`.  Old handoffs `HANDOFF_session_*.md`
are kept for the record; their TODO lists are superseded by this file.*

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

**(L) is the one open item.**  It needs NO fairness of anything: only that
the K₀−1 disjoint-row ladder pairs are not all crossed.  Plan:
- K₀ ≤ n/log³n costs O(1/log³n) by prop:trapped on A; on B the arc version
  is proved for |C₁| ≤ n/2 (splice-in switching); MISSING: P[B, d₁₂=ℓ,
  |C₁| > n/2, Flip=∅] = O(1/n) (rem:ladderdata(iii)).
- K₀ ≥ n/log³n: it suffices that some ladder pair has a j-cycle of length
  ≤ log⁵n avoiding j′ (a LOCAL witness; P ≈ ℓ/n per pair, expected log² n
  witnesses) — eq:Lloc — a polylog-size substructure statement about
  Θ(n/log³n) disjoint row pairs; target for KS switching / second moment
  (precision o(1/log² n) in covariances).  Data: at n=50, K₀ ≥ 12: a witness
  of length ≤ 10 exists in 93 %, ≤ 20 in 99.6 %.

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

-1. (L): (a) the arc-trapping gap rem:ladderdata(iii); (b) eq:Lloc by a Kwan–Sudakov
   switching in the space of completions of rows 1,2 and columns j,j′ (short
   j-cycles of disjoint ladder pairs); (c) measure covariances of the witness
   events W_k on real squares to see the precision needed.
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
statistics), `f1/f1_*.log` (flippability of bit-changing pairs: fair, independent).  `s13_spectral.py brute N` (N ≤ 7) and
`s13_spectral_big.py N` reproduce the theorem.

## Quickstart

HANDOFF_session_12.md §3 QUICKSTART still green, plus:
`python3 s13_spectral.py brute 6` (exact = brute, all digits);
`python3 s13_spectral_big.py 4000` (sum ν = 1, m·E → 1.9);
`gcc -O3 -o s13_frame s13_frame.c -lm && ./s13_frame 2 32 100000` → 4E[bias²] ≈ 0.115;
`pdflatex second_row_notes.tex` ×2 → 0 errors, 105 pp;
`gcc -O2 -DN=7 -o s13_comp7 s13_comp7.c && ./s13_comp7 0123456 1234560 0 2` → ladder identity exact (~1 min);
`./jm 30 300 30000 1 | python3 s13_ladder.py 30` → P[Flip=∅|K₀] ≈ 2^{−(K₀−1)}.
