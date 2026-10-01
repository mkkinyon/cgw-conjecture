# HANDOFF: Cameron's second-row conjecture — state after session 6

## >>> SESSION 6 HEADLINE: the MASTER FORMULA (Theorem thm:master, notes §sec:orbit, PROVED + verified on 1.7M real-square subset-checks) — the "non-linear orbit-mixing" bit is EXACTLY an F₂ rank event in a fixed interlacement matrix (Cohn–Lempel + a bordered-rank lemma + passenger/point-splitting constructions). The supply lemma (5′a) and conditional expander (5′b) are PROVED. The orbit-by-orbit mixing lemma of session 5 is REFUTED (too strong: short-arc orbits have bias +1/2 at measure ≍1/k); the open core is RE-SCOPED to an ANNEALED bias bound (Problem prob:annealed), for which generic position gives E|bias| ≈ 0.35·k^{−1/2} — a factor-5 exponent margin over the k^{−1/11} budget. <<<

*Prepared 2026-08-28 (end of session 6). Self-contained; the companion
tarball has all code, data, papers, notes (`second_row_notes.tex`/`.pdf`,
now ~44 pp with §sec:orbit). Previous handoffs remain valid except where
superseded (notably: the "KPS transfer programme" item 4 of
`HANDOFF_session_5.md` — the orbit-wise mixing lemma — is superseded
by the annealed reformulation below; that file is in the tarball).*

---

## 0. The problem in one paragraph

As before (HANDOFF_session_4 §0). Weak conjecture d_TV(P_n,Q_n)=o(1)
STILL OPEN; refined picture h_n(2) ≈ c₂n⁻³ retains merit (owner
directive: CGW Conjecture is the priority). By Theorem M (session 5),
EQ(δ) ⇔ the marked bit [rows 1,2 in one column cycle of
P={j,σ^α(j)}] has law 1/2+O(δ) for uniform L∈S(λ); EQ(o(1/log n)) ⇒
weak conjecture via thm:reduction.

## 1. What session 6 added (all in notes §sec:orbit = `section_orbit.tex`)

### Theorem thm:master (MASTER FORMULA; PROVED + VERIFIED)

Fix an orbit (base square L₀, marked pair P={j,j'}, stable collection
𝓘₀ in rows 3..n, configuration ω uniform on the cube). Then

    bit(ω) = [ x̄_{φ(ω)∪B} ∈ colspace( Ā_{φ(ω)∪B} ) ]   over F₂,

where Ā (alternating) and x̄ are EXPLICIT functions of the positions of
the intercalate rows in the cycles of π₀ = π_P(L₀), B is a frozen
bridge block (#cycles−1), and φ(ω) = lifted chords of active useful
coordinates. Ingredients (each with full proof in §sec:orbit):
- **Cohn–Lempel** (cited, [CL72], re-verified): #cycles(γτ_T) =
  null(A_T)+1 for a single circle and disjoint chords.
- **Bordered-rank lemma** (lem:border): for alternating A, bordering
  by (x,0) changes nullity by +1 iff x ∈ col(A), else −1. (Two-line
  proof; uses zᵀAz=0.)
- **Bit formula** (thm:bitrank): u~v in γτ_S ⇔ x_S ∈ col(A_S), where
  x = separation vector. Proof: adjoin the virtual chord e={u,v} (=
  the CGW switch at P!) and apply CL twice + bordered lemma.
- **Passenger construction** (lem:passenger): multi-cycle π₀ realized
  as single circle γ = (C₁ s₁ C₂ t₁ s₂ C₃ t₂ …) with frozen disjoint
  bridges (s_i,t_i); fresh points ride as passengers; first-return =
  π₀ exactly.
- **Point splitting** (lem:split): two-sided products π₀ρμ (ρ = j-side
  chords copy 1, μ = j'-side chords copy 0, per-side disjoint, sides
  may share rows): double every point; first-return on copies 0 =
  π₀ρμ. Conjugating coordinates (𝓘_{jj'}) contribute both lifts tied
  to one bit. KEY AUTOMATIC FACT: stable intercalates sharing a column
  can't share a row (they'd share a cell ⇒ not isolated) ⇒ per-side
  disjointness is free.
- Left factors are conjugated away using supports ∌ rows 1,2.

Verification: `pf_cl.py` (CL, BIT, MC, SPLIT, SPLIT2 — exhaustive over
all activations, random instances, incl. shared rows + conjugating
coords; conv rule: EARLIER-applied factor → copy 0). `pf_master.py`:
END-TO-END on real squares: n=5 complete (7920 checks), n=6 sampled
(1,706,949 subset-checks, families ≤4 intercalates incl. all four
column-incidence classes): **0 failures**.

### Bias identity + extremal configurations (prop:bias; PROVED)

bias := P_ω(bit)−1/2 = (1/2)E_S[null increment of the virtual chord e].
Also |bias| ≤ (1/2)P(coordinate i doesn't toggle) for ANY single i
(toggle status of i is a function of ω_{−i}).
Exactly solvable (verified `pf_bias.py`):
- **Ladder** (r pairwise-interlacing separating chords, nothing else):
  bit = [|S| even] — F₂-LINEAR parity, bias = 0 exactly ∀r≥1.
- **Nested** (r pairwise-non-interlacing separating): bit=[S=∅],
  bias = 2^{−r}−1/2 (worst case).
- **Cut chord** (one separating chord c crossing NO other chord):
  bit = 0 whenever c is active (e_c ∈ ker with x_c = 1), so EXACTLY
  bias = (1/2)·bias′ − 1/4 ∈ [−1/2, 0], bias′ = conditional bias given
  c inactive. So a cut chord caps bias at 0, and good structure
  elsewhere can't rescue balance: ladder + cut chord = −1/4 exactly
  (verified). (If c is the ONLY separating chord, bias = 0 — the −1/4
  needs a balanced remainder.) Moral: must EXCLUDE bad structure ⇒
  the mixing bound must come from the MEASURE.

### REFUTATION of the session-5 orbit-wise mixing lemma

If rows 1,2 are close in their π₀-cycle and no chord endpoint falls in
the short arc, x=0 and bit ≡ 1: bias = +1/2. Such orbits have measure
≍ 1/k (not o(1/log n)). So "conditional law = 1/2+o(1/log n) for all
but o(1/log n) of orbits" is FALSE in general. The annealed average
survives: short-arc contribution is O(1/k) — self-regularizing, no
join-then-toggle chain needed.

### Problem prob:annealed (the NEW open core)

    E_L | bias(orbit(L), P) | = o(1/log n)
    (L uniform on S(λ) ∩ expander; marks averaged; k ≥ 6log¹¹n supply)

Numerics (`bias_scale.c`, C, MC over subsets): iid-position ensemble:
E|bias| ≈ 0.35 k^{−1/2} (k=8: .124, 16: .099, 32: .075, 64: .044,
128: .032); signed mean = O(1/k)-compatible (≈0 within se ~3e-3);
flat in arc length except arcs with O(1) chord endpoints. ANY bound
k^{−c}, c > 1/11 suffices. Heuristic for signed mean: E null of random
chord systems ~ Harer–Zagier log-growth ⇒ increment ~1/k.

### Lemma lem:supply (5′a; PROVED modulo KPS Lem 8.3 (cited))

H = T∩L an (ℓ,β)-expander (the partial square is called H in the
notes; the letter P is reserved for the marked pair) ⇒ for all but
< tℓ columns c: t stable intercalates through c, pairwise disjoint
elsewhere, avoiding rows {1,2} and column j'. Proof: dead-column rounds + one expander
application per round (five βn-sets built from complements of used
rows/columns/symbols; C* = ℓ dead columns; monotonicity of Def 8.1).
t = 6ℓ ⇒ exceptional ≤ 6log²²n columns.
CAVEAT (Remark rem:bookkeeping(i)): exceptional columns are
square-dependent; for merges with λ_m·m = n^{o(1)} the per-square bad-
mark fraction isn't automatically small. Proposed fix (recorded, not
proved): TEMPLATE AVERAGING — run the whole decomposition with
relabelled template g(T), g uniform column permutation; every step
holds per fixed g; (g,mark)-averaged bad fraction ≤ 6log²²n/(n−1).

### Lemma lem:condexp (5′b; PROVED)

P(expander | S(λ)) ≥ 1 − e^{−ω(n log²n)} for EVERY λ, by direct
division: P(S(λ)) ≥ e^{−O(n log n)} using D_λ ≥ n!e^{−O(n log n)} and
N(ν)/N(λ) ≤ 3ⁿ via CGW Thm 3.13 (ratio ∈ [1/2,3/2] per merge; ≤ n
merges through the one-part type). NO Lemma-4.7 transfer machinery
needed.

## 2. File manifest (additions; all else as in HANDOFF_session_5)

| file | what | re-verify |
|---|---|---|
| `pf_cl.py` | CL/BIT/MC/SPLIT/SPLIT2 formula verification | `python3 pf_cl.py all` (~4 min; expect 0 failures, SPLIT conv0 fails BY DESIGN — records the wrong convention) |
| `pf_master.py` | end-to-end master formula on real squares | `python3 pf_master.py 5 8` (fast); `python3 pf_master.py 6 7 11` (~8 min, 0 failures) |
| `pf_bias.py` | structured bias examples (ladder/nested/cut) | `python3 pf_bias.py struct` (~1 min) |
| `bias_scale.c` | k-scaling of bias, MC over subsets | `gcc -O3 -o bias_scale bias_scale.c -lm && ./bias_scale 32 400 20000` (~1 min) |
| `section_orbit.tex` | notes §sec:orbit (input by main tex) | pdflatex second_row_notes.tex (clean, no undef refs) |
| `pf_cl.log`, `pf_master.log`, `bias_scale.log` | session-6 run logs | — |

Bibliography addition: [CL72] Cohn–Lempel, JCTA 13 (1972) 83–89.

## 3. PRIORITIZED TO-DO

**The strategic target is now Problem prob:annealed (annealed bias
bound). Suggested attack order for session 7:**

### TODO 7a (the mathematical core): position-genericity under the
orbit measure. Two sub-statements:
  (i) P(some available separating chord is a CUT chord, or the
  separating family is nested-degenerate) = o(1/log n) under uniform
  L ∈ S(λ)∩expander — an anti-concentration statement about stable-
  intercalate row positions in the cycles of π_P. Tools: row-3..n
  exchangeability acts by conjugation (NOT template-equivariant — see
  pitfall 22), KSS intercalate density (N ~ n²/4), the expander
  itself (chords with rows in prescribed βn-sets EXIST — use to show
  crossings are plentiful, not to build ladders).
  (ii) Quantitative rank-balance for the resulting random
  interlacement matrices: E|bias| ≤ k^{−c}. Possible routes: second
  moment of the rank event via the kernel-sum identity
  bit = 2^{−null}Σ_{z∈ker}(−1)^{z·x}; Harer–Zagier-type nullity
  computations; interlace-polynomial recursions (pivot/local
  complementation); or a coupling that plants a random ladder.
  START by testing candidate sufficient conditions in pf_bias.py
  (e.g. "every separating chord crossed by ≥ m others" vs bias — the
  cut-chord example shows m=0 fails; is m=1 enough? unknown).

### TODO 7b: make Remark rem:bookkeeping(i) (template averaging) a
proof. Check that KPS Lemma 8.3's template T can be replaced by g(T)
with the SAME failure bound (yes — relabelling symmetry of uniform L),
and write the (g,mark)-averaged assembly of steps 1–5 of §sec:orbit's
summary.

### TODO 7c: real-ensemble grounding at moderate n. Sample S(λ) at
n = 30–100 (MCMC via Jacobson–Matthews — `jm.c` / `jm_stats.py` /
`mirror_jm.py` are IN THE TARBALL from earlier sessions, see
HANDOFF_session_4 manifest — or intercalate switching from
sec:numerics machinery), extract maximal disjoint intercalate families
through a marked column, compute (Ā,x̄) and bias: does E|bias| track
0.35k^{−1/2}? Any cut-chord/nested pathologies at realistic densities?
(This is the empirical version of TODO 7a and will likely reveal the
right genericity condition.)

### TODO 1″ (referee pass): unchanged from session 5, PLUS: referee
§sec:orbit (the proofs are short; check lem:split's convention against
pf_cl SPLIT test, and lem:supply's set-size arithmetic
4tℓ+2t+ℓ+2 ≤ βn/2).

### TODO 2′, TODO 4: unchanged (owner: real merit, lower priority).

## 4. Pitfalls (all previous, plus)

21. In the master formula the EARLIER-applied factor attaches to copy
    0: for π₀∘ρ∘μ (μ applied first), μ-chords → copies 0, ρ-chords →
    copies 1. pf_cl.py's SPLIT test records the failure of the other
    convention (conv0: 122 failures BY DESIGN). π_P orientation:
    π(r) = row of column j' holding symbol L(r,j); j-side switches
    PREcompose (right), j'-side POSTcompose (left).
22. Row-permutation symmetry of S(λ) is NOT template-equivariant:
    permuting rows 3..n preserves uniform-on-S(λ) and conjugates π_P,
    but maps T-stable intercalates to g(T)-stable ones. Any symmetry
    argument must either average over templates (cf. TODO 7b) or
    avoid template-dependence.
23. Orbit-wise mixing (session-5 formulation) is FALSE — do not try
    to prove it; short-arc orbits have bias +1/2 at measure ≍ 1/k.
    The annealed form absorbs them at cost O(1/k).
24. In bias_scale.c the u-slot is 0 and chords never contain u,v
    slots; "arc length" is counted in SLOTS (chord-endpoint
    positions), not rows — when comparing to real data convert via
    endpoint density.
25. A "neither-column" intercalate leaves π_P untouched — but only
    cell-disjointness makes subset-switching well-defined; arbitrary
    families of intercalates (not pairwise cell-disjoint) do NOT give
    a cube. (pf_master uses greedy cell-disjoint families.)

## 5. Papers on hand / to request

**IN THE TARBALL (unchanged from session 5 — CGW_2008.pdf,
KPS_2025.pdf, KSS_large_deviations.pdf, Schramm_2005.pdf,
Kwan_Sudakov_2018.pdf — plus notes PDF). Keep packing all of them.**

**To request if needed (session 7):**
- Cohn–Lempel 1972 (JCTA 13, 83–89) — cited for Theorem thm:CL; the
  statement is re-verified computationally, but the referee pass
  (TODO 1″) may want the original. Also possibly Traldi, "Binary
  nullity, Euler circuits and interlacement graphs" (generalized CL —
  NOT currently needed: the passenger/splitting constructions avoid
  it).
- Arratia–Bollobás–Sorkin interlace polynomial papers + Harer–Zagier
  (1986) IF TODO 7a(ii) goes the nullity-generating-function route.
- KSSS 2202.05088: same conditional status as before (probably NOT
  needed — lem:condexp avoids Lemma 4.7).

## 6. Quick reference

Master formula: bit = [x̄_{φ(ω)∪B} ∈ col Ā_{φ(ω)∪B}]. bias =
(1/2)E[null increment of virtual chord e] = (1/2)E[±1]. Ladder ⇒
parity ⇒ bias 0; nested ⇒ 2^{−r}−1/2; cut chord c ⇒ bias =
(1/2)bias′ − 1/4 (bias′ = given c inactive; ladder remainder ⇒ −1/4).
Generic:
E|bias| ≈ 0.35k^{−1/2}; budget k^{−1/11} at k = 6log¹¹n. q̄ values
unchanged (session-5 list). P(S(λ)) ≥ e^{−O(n log n)} (lem:condexp
proof); CGW Thm 3.13: merge ratio ∈ [1/2, 3/2].

*The notes are the living document; keep PROVED/VERIFIED/COMPUTED
labels honest. Session-6 headline: the orbit bit is a rank condition
(master formula, PROVED+VERIFIED); open core re-scoped to the annealed
bias bound with a factor-5 exponent margin in generic position.*
