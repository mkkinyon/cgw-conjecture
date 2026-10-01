# STATUS — CGW Conjecture project (rolling; replaces per-session handoffs)

*Updated 2026-10-01 (session 13, continuous).  The notes
`second_row_notes.tex` are the living document; this file is the two-page
map.  History: `git log`.  Old handoffs `HANDOFF_session_*.md` are kept for
the record; their TODO lists are superseded by this file.*

## Goal and chain

CGW weak conjecture d_TV(P_n,Q_n) → 0  ⇐  EQ(o(1/log n)) [thm:reduction]
⇐ marked bit P(B) = 1/2 + o(1/log n) [thm:marked]
⇐ KPS re-randomisation + rank formula [thm:master] + supply/cond. expander
⇐ prob:annealed  ⇐  (I) prob:mixing  +  (II) transfer iid → orbit  [(III) closed]
(I): prob:mixing (ψ(j,d) ≤ K d^{-2c}, c > 1/11)  ⇐  (S₂) on S'_K with γ > 1/11
      + entrance estimate (lem:srecharge + occupation-time lemma, OPEN)
      + for z=−1: prop:holder with p > 6/5 (crude) — measured p ≈ 2.7.

## Open items (what a proof still needs)

1. **(S₂)** [§sec:s13ltwo]: E_{ν_t}[Φ_m²] ≤ C(K) m^{−2γ}, γ > 1/11, for the
   time-t entrance laws of S'_K.  Dust-blind L² statement; replaces the
   sup-norm (S).  Routes: spectral (H_β) inequality in L²(π̄) [§s13spectral],
   or exact-law route (a_k = O(k^{−6/5}) from thm:blockcount).
2. **Occupation-time lemma** [lem:occupation, OPEN, believed standard].
3. **(R1′)** E_x[τ] ≤ C(K)(1+1/g) on H, and the (1+γ)-moment of τ with
   γ > 1/11 [§s13m1]; TODO 11a's t^{−2} tail would be far more than enough.
4. **(II)(T1)+(T2)** [§s13transfer]: health of the frame and coarse product
   genericity of the intercalate pool under Unif(S(λ)).  HARD, different in
   kind from (I): distributional, cannot be imported by union-of-costs; only
   visible route is a KPS-style "all but ℓ columns" robustness argument in
   the triangle-removal model.  Numerics: true with margin.
5. **(T3)** re-pose (S₂)/entrance from the stationary start (the true frame
   is PD(1)-like, not one-block); budget unaffected, c = min(1/2, γ).

## Closed this session (13)

- Budget corrected: γ > 1/11, not 4/9 (cor:budget).  p > 6/5 not 13/5.
- L² form of the reduction (prop:S2red); dust obstruction does not touch it
  (rem:dustblind).  Plain annealing E[Φ] is NOT sufficient.
- τ_k < ∞ a.s. from any state with ab > 0 (lem:flipsrecur).
- (III) prob:decouple NOT NEEDED: routing chains merge parts m ≥ n/(A log n)
  (prop:nodecouple).  Also: cycle type of σ is a function of the bottom
  rectangle alone (orientations are fair coins).
- Factor-2.2 iid-vs-real-squares discrepancy explained: the real frame is
  stationary (s13_frame.c: one-block/antipodal 0.3/k, random frame 0.9/k,
  real 0.7/k).

## Do not reopen (proved impossible / false)

Doeblin minorisation of the post-flip chain in any form (prop:nomino);
uniform or V-uniform ergodicity of κ₁ (cor:noerg); λ-drift on the full
space (prop:nodrift); geometric (C2) (§s12num: a_k ~ k^{−3}); (R2)
(prop:r2false); transport-metric contraction (§s12left); the two-block
odd-perturbation argument (§s12odd); reading ρ from burnt-in starts
(pitfall 56).

## Running / recent numerics

`runs/s13/akB_0_p.log` (4e8 samples, pin p), `akB_{6,7}_k40.log`
(K-dependence of the tail constant), `frame*.log` (frame comparison).
Handoff numeric task 1 ((C2″) on the unit circle) not yet run.

## Quickstart

See HANDOFF_session_12.md §3 QUICKSTART (all green at start of s13) plus:
`gcc -O3 -o s13_frame s13_frame.c -lm && ./s13_frame 2 32 100000` → E[bias²]k ≈ 0.92.
`pdflatex second_row_notes.tex` ×2 → 0 errors, 91 pp.
