# STATUS — CGW Conjecture project (rolling; replaces per-session handoffs)

*Updated 2026-10-01 (session 13, continuous).  The notes
`second_row_notes.tex` (95 pp, 0 errors) are the living document; this file
is the two-page map.  History: `git log`.  Old handoffs `HANDOFF_session_*.md`
are kept for the record; their TODO lists are superseded by this file.*

## Goal and chain

CGW weak conjecture d_TV(P_n,Q_n) → 0  ⇐  EQ(o(1/log n)) [thm:reduction]
⇐ marked bit P(B) = 1/2 + o(1/log n) [thm:marked]
⇐ KPS re-randomisation + rank formula [thm:master] + supply/cond. expander
⇐ prob:annealed: E|bias| = o(1/log n) over orbits, k' ≥ 6 log^11 n chords.

**prob:annealed in the stationary-frame MODEL is PROVED with the sharp
exponent** (§sec:spectralexact): 4E[bias²] ≤ C/k' + Ck'²/N for a uniform
frame and a uniform disjoint pool (so k' ≪ √N; the real pool is a central
element and the restriction can presumably be lifted). It rests on the exact
spectral measure of the same-cycle sign (thm:spectral): ν(λ) = 2c_λ²(1−r_λ)
on shapes (a,b,1^j), c_λ = (a−b+1)/((a+j+1)(b+j)); E_π̄[Φ_m] ≤ 12/m.
The chain is: thm:reduction ⇐ EQ ⇐ thm:marked + [KPS re-randomisation,
thm:master, lem:supply, lem:condexp, prop:nodecouple] ⇐ **(T)**, where
(T): E_{Unif(X)}[ε(π₀R_SM_S) ε(π₀R_TM_T)] = o(1/log² n) over the real orbit
measure. cor:stationarymodel says (T) holds (≈ 2/k') in the model.

## The one open item: (T), i.e. (II) the transfer — what it contains

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
- **Trade-chain route** (§s13tradechain), an alternative that bypasses (b)
  and intercalates entirely: symmetric chain on X preserving S(λ) and the
  mark; EQ(4δ) ⇐ bit autocorrelation along the chain ≤ δ². Merges with a
  free cycle are free; EVERY split and the C₁–C₂ merge are thinned by
  "flippability" (measured ½, uncorrelated with the frame beyond adjacency;
  exact criterion: q* off the ρ_{x,y}-cycle of j). Needs (F1) flippable
  bit-changing pairs in positive proportion for typical squares — merges in
  the apart state AND separating splits in the together state — and (F2)
  mixing of the thinned walk. Numerics: autocorrelation tracks the pure walk
  to the noise floor (n=30). If everything were flippable, EQ(o(1)) would
  follow from Diaconis–Shahshahani on cosets (P[u~v]=½ exactly on each coset).

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

0. (F1) by a switching argument; (F2) by comparison with the full transposition walk.

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
`akB_{6,7}_k40.log` (moot).  `s13_spectral.py brute N` (N ≤ 7) and
`s13_spectral_big.py N` reproduce the theorem.

## Quickstart

HANDOFF_session_12.md §3 QUICKSTART still green, plus:
`python3 s13_spectral.py brute 6` (exact = brute, all digits);
`python3 s13_spectral_big.py 4000` (sum ν = 1, m·E → 1.9);
`gcc -O3 -o s13_frame s13_frame.c -lm && ./s13_frame 2 32 100000` → 4E[bias²] ≈ 0.115;
`pdflatex second_row_notes.tex` ×2 → 0 errors, 95 pp.
