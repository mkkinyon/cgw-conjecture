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

**prob:annealed in the stationary-frame model is PROVED with the sharp
exponent** (§sec:spectralexact, session 13): 4E[bias²] ≤ C/k' + Ck²/N.
It rests on an exact closed form for the spectral measure of the
same-cycle sign under the random-transposition walk (thm:spectral):
ν(λ) = 2c_λ²(1−r_λ) on shapes (a,b,1^j), c_λ = (a−b+1)/((a+j+1)(b+j)).
Consequence: E_π̄[Φ_m] ≤ 12/m uniformly in N (cor:stationarydecay), m·E → 2.

## The one open item: (II) the transfer — three equivalent-ish formulations

Needed: E_{Unif(X)}[ε(π_P(L^S)) ε(π_P(L^T))] = o(1/log² n) (pool subsets S,T).
- **(T1) PROVED** (prop:trapped, refereed): P[row 1 in a frame cycle of length ℓ] ≤ 18/(n−ℓ−1)
  under Unif(X), by a Kwan–Sudakov switching whose only row-1 move is the CGW switch
  at the mark (priced by thm:marked). Trapped marks cost O(1/k′).
- **(T2) OPEN**: the pool (intercalate row pairs) in generic position relative to the
  frame. Exact structure (§s13transfer2): delete columns j,j′; frame = leftover graph
  of the n×(n−2) rectangle; pool = ¼-coin thinning of a candidate graph; only local
  exclusion = π-adjacent rows. Real squares (true marks) decorrelate 10–30% FASTER
  than the model at n=50 (adjacency exclusion), frame stats agree.
- **Trade-chain route (§s13tradechain) OPEN, promising**: symmetric Markov chain on X
  (trades at row pairs ∉{1,2}, turns of cycles avoiding 1,2) preserving S(λ) and the
  mark, acting on the frame by transpositions thinned by "flippability" (prob ½, no
  frame correlation except adjacency). EQ(δ) ⇐ bit autocorrelation along the chain
  ≤ δ². Needs (F1): in the apart state a positive fraction of C_1×C_2 pairs flippable
  (typical squares), (F2): mixing of the thinned walk. Numerics: autocorrelation
  tracks the pure walk to the noise floor (n=30). No intercalates/expander needed.

## Bugs found this session
- Sessions 7–8 real-square pipelines (s7/s8_real_bias.py) used σ's SYMBOLS as column
  indices on unnormalised JM squares: their "marks" were random column pairs.
  Qualitative findings survive; §sec:realbias numbers are for random pairs.
  Fixed in s13_real.py / s13_trade.py (symbol → column via row 1).

## Closed / superseded (do not spend time on these)

- (III) prob:decouple: NOT NEEDED (prop:nodecouple; routing merges m ≥ n/(A log n)).
- (I) in all its session 8–12 forms: prob:mixing (one-block start), prob:survival,
  (R1), (C1), (C2″), (S₂), occupation-time lemma, Markov renewal, a_k, t_k(w):
  all MOOT for the stationary frame.  The exact laws (thm:blockcount,
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
