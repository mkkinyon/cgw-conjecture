# HANDOFF: Cameron's second-row conjecture / CGW Conjecture — state after session 8

## >>> SESSION 8 HEADLINE: prob:mixing is now a MARKED SPLIT–MERGE
problem with an exact symmetry. PROVED (notes §sec:mixing, new =
`section_mixing.tex`): (1) the fresh-chord operator P IS one step of
uniform split–merge (coagulation–fragmentation) on the measured
partition of the circle into γτ-cycles, with ε = the together-
indicator of u,v (lem:splitmerge; continuum simulator
`s8_splitmerge.c` reproduces the discrete ψ, one-point, and PD(1)
invariants end-to-end); (2) the chain has an exact sign-reversing
involution ι (swap together↔apart keeping the geometry; thm:iota,
verified numerically to MC precision), whence M_d = ε₀·Φ_d(quotient)
with Φ_d = E[(−1)^{#flips}] over an autonomous quotient chain
(cor:quotient) — the sign structure is GONE from the problem; (3) the
stationary gap law is exact: π(g ≤ δ) = 4δ + O(δ²)
(prop:gaplawstat, from uniform-permutation classics; verified); (4)
LOWER BOUND PROVED: ψ_π(d) ≥ c/d (prop:lower) — the measured rate
ψ ≈ 1/d is SHARP, c = 1/2 is the true exponent ceiling of
prob:mixing, margin vs needed c > 1/11 is definitively 5.5×.  The
upper bound is reduced to two scoped sub-problems: prob:gaplaw
(uniform-in-j linear gap law; proved at stationarity, creation-flux
side proved pointwise (lem:recharge), missing only a "no-dust"
escape estimate) and prob:survival (parity decay off small gaps;
renewal-with-infinite-mean endgame mapped out, recharge lemma
proved), plus a spectral route via reversibility (lumped
transposition walk; odd-sector Dirichlet form has a massive term on
flips).  TRANSFER (8b): reduced to a de Finetti mixture-of-iid(θ)
step with O(k²/n) error — the whole split–merge analysis is
θ-covariant, so thm:factor and everything above apply verbatim per
mixture component (§sec:transfer); empirically E[bias²|k]·k = 0.6–0.75
flat over k=8..26, n=30..100 (and n=200: SEE BELOW).  NUMERICS:
iid 4Ebias²·k FLATTENS (k=512: 1.45±0.23, k=1024: 1.17±0.32 vs
k=256: 1.44) — constant ≈1.3–1.45, not log; stationary ψ·d =
0.96/1.06/1.11/0.93/0.87 at d=64..1024, exponent 1.03±0.09, with the
ψ|g₀-scaling collapse confirmed (function of g₀·d alone, transition
at g₀ ≈ 1/d, together/apart bins identical = ι in action); one-point
f(m) → 0.47/m² exactly; n=200: no n-drift. <<<

*Prepared 2026-08-29 (end of session 8). Self-contained; companion
tarball has all code, data, papers, notes (`second_row_notes.tex` now
inputs `section_mixing.tex` = §sec:mixing). Previous handoffs valid
except: session-7's TODO 8a/8b option lists are superseded by the
concrete reductions of §sec:mixing; the session-7 ψ(j,d) table
carries ±0.005 sampling scatter (its inner-MC se were optimistic) —
its exponent conclusion stands. `HANDOFF.md` = this file.*

---

## 0. The problem in one paragraph

As before (HANDOFF_session_4 §0). Weak conjecture d_TV(P_n,Q_n)=o(1)
STILL OPEN; owner directive: CGW Conjecture is the priority. Chain:
thm:marked → programme steps 1–5 → open items now:
  (I)   prob:mixing, REDUCED to prob:gaplaw + prob:survival
        (§sec:mixing; rate 1/d proved sharp from below)
  (II)  transfer iid → orbit, REDUCED to mixture/de Finetti step +
        θ-uniformity (§sec:transfer)
  (III) prob:decouple (analytic framing unchanged; session-8
        empirics now supportive: supply ⊥ class size, see 9d)

## 1. What session 8 added (all in notes §sec:mixing = `section_mixing.tex`)

### Lemma lem:splitmerge (PROVED + VERIFIED end-to-end)
Fresh chord = one uniform split–merge move on the curve partition
(cycles of γτ_W as measured partition of the circle); marked state =
(curve masses; u,v curves; sides (x,y) if together else blocks (p,q));
all transition laws need only uniform variates (full case table in the
lemma).  `s8_splitmerge.c` (continuum, no surgery per step, ~1000×
faster at large d) matches the discrete `s7_twopoint` at fs=1:
ψ(16,4): .1341(22) vs .1305(22); ψ(64,16): .0391(26) vs .0395(26);
one-point f(4): +.0232(4) vs +.0237(5).  PD(1) checks: EΣp²=.4979
(pred .5), E[marked mass|tog]=.6685 (pred 2/3).

### Theorem thm:iota (PROVED + VERIFIED): sign-reversing involution
ι: (together,(x,y),env) ↔ (apart,(p,q)=(x,y),env) commutes with the
kernel; ε∘ι = −ε; hence M_d∘ι = −M_d.  Verified: mirrored starts give
P_A(tog)+P_B(tog)=1 at all t, identical u-marginals (`iota` mode).

### Corollary cor:quotient (PROVED): M_d = ε₀·Φ_d(Γ), ψ = E[Φ²]
Quotient chain Γ=(a,b;env) autonomous; flips (prob 2ab per step,
(a,b)→(α+β, a+b−α−β)) are quotient moves; label = ε₀·(−1)^{#flips}.
Any stationary law is ι-invariant: label ⊥ geometry, P(tog)=1/2 ✓
(measured .5065(29)).

### Prop prop:gaplawstat (PROVED*): π(g≤δ) = 4δ+O(δ²)
π = uniform-permutation limit: given tog, mass s ~ 2s ds, x|s ~ U(0,s);
given apart, (p,q) uniform on the triangle p+q<1.  (*finite-n proof +
routine limit; verified in dyadic bins.)  NOTE marked-total s has
QUADRATIC small mass in both sectors → s-traps contribute only O(1/d²).

### Prop prop:lower (PROVED): ψ_π(d) ≥ c₀/d
Freeze argument: on {g₀ ≤ 1/(Kd)} the g-material is hit by no cut in d
steps w.p. ≥ 1−4/K (two cuts/step, unhit ⇒ frozen), no hit ⇒ no flip ⇒
|M_d| ≥ 1−8/K; mass of the event ≥ c/(Kd) by the gap law. K=16.
CONSEQUENCE: no c > 1/2 attainable; measured rate is exact; any proof
of prob:mixing must tolerate the trapped trajectories.

### Lemma lem:recharge (PROVED): the upper-bound toolkit
(i) flips preserve s = a+b; (ii) post-flip gap quadratically safe:
P(g'≤η|flip from (a,b)) ≤ η²/(ab); (iii) per-step creation flux into
{g≤η} ≤ 4η² from ANY state (channel-by-channel densities ≤ 2η).
Endgame mapped: flips ≈ renewal process with infinite-mean tail
~c/t; alternating renewal cancellation gives Φ_d = Θ(tail(d)) = Θ(1/d);
the two genuine gaps: state-dependence of inter-flip times, and the
escape-rate ("no-dust") estimate — at stationarity dust has Dickman
mass ~1e−8, so the assembly must let dust delay, not veto (dust
coarsens under merges).

### Problems restated (OPEN, sharply scoped)
prob:gaplaw: P_{μ_j}(g₀≤ε) ≤ Cε^κ uniform in j (κ=1 at stationarity
PROVED; creation flux PROVED; only escape/no-dust missing).
prob:survival: E[Φ_d²; g₀≥ε] ≤ C(εd)^{−2α}.  ANY (κ,α) with
ακ/(κ+2α) > 1/11 closes the iid model ((1,1/2) ⇒ c' = 1/4; κ=α>3/10
suffices).  Spectral formulation (§sec:spectral): marked chain =
strong lumping of the reversible transposition walk ⇒ reversible;
ν_ε(top-δ) ≤ Cδ ⟺ upper bound; odd-sector Dirichlet form has massive
flip term E[2ab(F+F')²] forcing low-energy odd functions onto
{g≲δ} (mass ~δ): weighted-Hardy formulation.  CAVEAT: chain is
2-periodic ((−1)^{#blocks} is an exact λ=−1, ι-odd eigenfunction);
work on P² / per parity sector.

### §sec:transfer (8b REDUCED): mixture route
Orbit chords exchangeable given frame ⇒ mixture of iid(θ) with
O(k²/n) TV error (Diaconis–Freedman) + hard-core correction of same
order (KSS-controllable); thm:factor holds per component; the ENTIRE
split–merge analysis is θ-covariant (rates = θ-masses; needs θ
atomless — orbit atoms O(1/n) harmless).  Two open steps: (a)
quantitative exchangeable→mixture for the greedy family, (b)
θ-uniform gaplaw/survival.  8b(γ) DONE from s7 CSVs: E[b²|k]·k flat
0.6–0.75 across k=8–26, n=30–100.

### Numerics (8d)
- iid k=512: 4Eb²·k = 1.448±0.229; k=1024: 1.171±0.324 — the
  0.95→1.44 drift has FLATTENED: constant ≈1.3–1.45, NOT log k.
  8d(b) settled: no log factor.
- Stationary ψ·d: 0.96 / 1.06 / 1.11 / 0.93 / 0.87 at
  d=64/128/256/512/1024 (errors ±0.05→±0.25) — constant ≈1.0, fitted
  exponent 1.03±0.09: pure 1/d.  (Reproducibility: the d=64 point was
  an interactive run with burn=4096, not in any log; from the logged
  points alone `python3 s8_harvest.py` reports slope 1.093 over
  d=128–1024.)  ψ|g₀ is a scaling function φ(g₀·d) alone
  (dyadic profile .06/.23/.49/.69/.82/→1 at u≈.7/.35/.18/.09/.04;
  φ≤.02 for u≥1.4; transition at g₀≈1/d; tog/apart bins identical =
  ι in action; small-u: 1−φ ≈ u^{.5-.6}, unexplained, low stakes).
  CAUTION: a d=2048 point at useful precision is inherently ~9 CPU-h
  (se floor 1/√(cs·fs), work = 2d·cs·fs) — budget it as an early
  launch or skip.
- n=200 real ensemble (s8_inst200.csv, 144+ inst, k=26 cap): E|b| =
  0.118, Eb²·k = 0.87±0.15, signed −0.006±0.014, P(r=0)=0.056±0.019 —
  fixed-k agreement with n=50/100 (0.107–0.116) confirmed at n=200;
  no n-drift in the 1/k law.  NOTE jm.c is ~4× slower at n=200 than
  n=100 (cache): ~8 min/square at thinning 2n³ — plan accordingly
  (a second
  stream yielded only 1 more square and was discarded).
- One-point from clean start: f(m) = .006280(50)/.002951(41)/
  .001713(35)/.000814(29)/.000458(25) at m=8/12/16/24/32; f·m² →
  0.469, local exponent 2.00 on m=24→32: f(m) ≈ 0.47/m² EXACTLY the
  trap-creation prediction (clean start must CREATE its 1/m-trap at
  cost 1/m², lem:recharge(iii); stationary start has the linear atom
  ready-made → 1/d).  8d(c) resolved.

## 2. File manifest (additions)

| file | what | re-verify |
|---|---|---|
| `s8_splitmerge.c` | continuum marked split–merge: psi/psibar/psig/one/traj/iota/gapdist modes | `gcc -O3 -o s8_splitmerge s8_splitmerge.c -lm && ./s8_splitmerge psi 64 16 6000 50 992` (≈.039, matches discrete) |
| `section_mixing.tex` | notes §sec:mixing (input by main tex) | pdflatex second_row_notes.tex (clean) |
| `s8_psibar.log` | stationary ψ̄(d) runs, gap-binned | — |
| `s8_4eb2_big.log` | iid k=512, 1024 | — |
| `s8_real200.log`, `s8_inst200.csv` | n=200 real-ensemble run | — |
| `s8_real_bias.py` | s7_real_bias + mcyc (σ-cycle length) column + k=0 rows | as s7 version |
| `s8_dec{30,50}.csv` | 8c data: corr(k, mcyc) = +0.006±0.018 (n=30, 3200 inst), −0.034±0.022 (n=50, 2000 inst); P(k=0)=0 throughout; E[k|m] flat — supply ⊥ class size empirically (recorded in §sec:vacuity status) | `python3 s8_harvest.py` |
| `s8_harvest.py` | all session-8 summary tables from logs/CSVs | `python3 s8_harvest.py` |
| `s8_onept.log`, `s8_psibar1024.log` | one-point suite; d=1024 ψ̄ | — |

Notes edits: section_orbit.tex (lem:supply β≤1/2 caveat folded in — 
TODO 1‴ DONE), second_row_notes.tex (inputs section_mixing; open-problem
§ updated).

## 3. PRIORITIZED TO-DO (session 9)

**QUICKSTART: `python3 pf_cl.py all` (0 failures; SPLIT conv0 122 BY
DESIGN), `python3 pf_master.py 5 8` (0 failures), then
`gcc -O3 -o s8_splitmerge s8_splitmerge.c -lm && ./s8_splitmerge psi
64 16 6000 50 992` → psi=0.034573 EXACTLY (seeded, deterministic;
the unbiased value of ψ(64,16) is 0.039, see pitfall 33 — the seeded
run's exact output is the reproducibility check, not the best
estimate).  Then launch long numerics first (see 9e).**

### TODO 9a (mathematical core): prove prob:survival.
Route: renewal-with-infinite-mean assembly. Concretely:
(i) first-flip tail from good starts: P_Γ(τ>t) ≤ C/(g₀... t)^α via
occupation bounds (creation flux lem:recharge(iii) + escape);
(ii) recharge (lem:recharge(ii)) ⇒ post-flip states are good w.h.p. ⇒
regeneration; (iii) alternating-renewal cancellation with
state-dependent times: dominate/couple against an iid renewal with
tail c/t; the signed renewal measure Σ(−1)^kG^{*k} has mass 1/2 ⇒
Φ_d = Θ(tail). Watch the 2-periodicity. ANY α > 0 with the κ=1 gap
law gives c' = α/(1+2α) — even α = 1/8 closes the budget.

### TODO 9b: prove the no-dust/escape estimate (finishes prob:gaplaw).
Need: P(escape [0,2ε] per step | g ≤ 2ε, F_t) ≥ cε "usually" —
i.e. a macroscopic merge partner is available except on rare,
short-lived dust episodes. Suggested: Lyapunov on L₁ = largest block
(or on Σp²): merges dominate in dust, coarsening is fast; a
time-fraction bound E[frac of [0,d] with L₁ ≤ 1/8] ≤ small suffices
(dust need only DELAY escapes). Alternatively: prove κ=1/2 version
unconditionally via E[g^{-1/2}] supermartingale-ish bookkeeping
(creation ∫η^{-1/2}·η dη converges) — κ=1/2, α=1/2 gives c'=1/8>1/11 ✓
— possibly the FASTEST route to closing the iid model.

### TODO 9c: transfer (8b cont.): (a) write the exchangeable→mixture
step for the greedy family under the orbit measure (conflict density
via KSS); (b) θ-uniformity of 9a/9b (should be cosmetic: constants
depend only on max atom); (c) numerics: accumulate more n=200
instances if wanted (all 144 current ones sit at the kcap=26 bin,
Eb²·k = 0.87±0.15 already computed — more data or a higher kcap is
the only way to sharpen this).

### TODO 9d: prob:decouple — analytic (L↓ framing). The empirical
part is DONE this session: supply ⊥ class size on real squares
(corr ≈ 0 ± 0.02 at n=30,50; §sec:vacuity status note). Next step is
the orientation-randomness / per-column failure quantification.

### TODO 9e: numerics to launch early (ONLY if wanted — the s8 data
may already suffice; k=1024, the one-point suite, and n=200 are all
DONE, do not rerun them):
- psibar d=2048: work ≈ 2d·cs·fs steps at ~1e7 steps/s/core; a useful
  point (se ≤ ψ/4, i.e. cs·fs ≥ 7e7 at ψ≈5e-4) is ~9 CPU-hours —
  launch at session start or skip (d ≤ 1024 already gives exponent
  1.03±0.09; a d=2048 point mainly guards against a very slow log).
  d ≥ 4096 needs a better estimator or a bigger machine; do NOT
  launch it naively (session 8 burned 5h on a mis-budgeted attempt).
- n=200 --verify run (master formula ground truth at n=200; see §6
  note): ~10 squares with --verify 50 suffices, ~80 min of jm.c.
- Optional: exact one-point f(m) for m ≤ 3 by symbolic integration
  over the split–merge tree (would pin the 0.47 constant's origin;
  HZ-adjacent).

### TODO 1⁗ (referee): §sec:mixing is new — referee it. Also
prop:gaplawstat's finite-n → continuum limit exchange deserves a
written proof (currently "routine"; the lumping/reversibility
statement in §sec:spectral likewise).

## 4. Pitfalls (all previous, plus)

32. The split–merge chain is 2-PERIODIC ((−1)^{#blocks} flips every
    step; it is ι-ODD). Spectral arguments must use P² or parity
    sectors; do not assume aperiodicity. Empirically no λ=−1
    obstruction in ψ (it decays), consistent with phase symmetry.
33. ψ estimation: the s7 discrete table has ±0.005 scatter (inner-MC
    correlations make the printed se optimistic at fs>1). For precise
    ψ use fs=1 (exact se = 1/√cs) or the continuum code with large fs
    AND C-level variance accounted.
34. In s8_splitmerge.c, marked-block splits treat the mark at
    position 0 of the cycle parametrisation (cuts uniform relative to
    the mark): correct because cut positions along a curve are
    uniform ⟂ the mark's position. Don't "optimize" this.
35. The flip move preserves a+b (marked total): a good invariant
    check for any modification of the simulator.
36. Under μ_j (j iid chords), u and v are ALWAYS on curves — there is
    no "u,v on the same point" atom; but the SHORT-ARC atom (session
    7) corresponds to together-states with one tiny side, i.e. the
    g-trap: same object, two languages.

## 5. Papers on hand / to request

**IN THE TARBALL (unchanged): CGW_2008.pdf, KPS_2025.pdf,
KSS_large_deviations.pdf, Schramm_2005.pdf, Kwan_Sudakov_2018.pdf +
notes PDF. Keep packing all of them. (kps.txt kept.)**

**To request if needed (session 9):**
- Diaconis–Freedman "Finite exchangeable sequences" (1980) — for 9c(a).
- Diaconis–Mayer-Wolf–Zeitouni–Zerner (2004, Aldous conjecture /
  uniform split–merge) and/or Berestycki EFC papers — for 9a/9b
  context (Schramm_2005.pdf already packed).
- Dynkin–Lamperti / infinite-mean renewal references (e.g. Bingham–
  Goldie–Teugels ch. 8) — for the alternating-renewal endgame.

## 6. Quick reference

Master formula: unchanged, VERIFIED to n=100 (0/3250 ground-truth
switch checks); at n=200 the bias statistics are computed but the
--verify ground-truth check was NOT run there (cheap 9e item). 
Factorization: 4Ebias² ≤ max ψ(j,d) + 2e^{−k/32}; ψ(j,d) =
E_{μ_j}[Φ_d(Γ)²]; Φ_d = E[(−1)^{#flips}]; flips: prob 2ab, (a,b) →
(α+β, a+b−α−β), preserve a+b. π: P(tog)=1/2; tog: s~2s ds, x|s~U;
apart: (p,q) uniform on triangle. π(g≤δ) = 4δ+O(δ²). ψ_π(d) ≥ c/d
(PROVED) and ≈ 1.05/d (measured). Budget: any c > 1/11, i.e.
ακ/(κ+2α) > 1/11. iid 4Eb²·k → ≈1.45 flat. Real ensemble: E|b| ≈
1.5–1.6× iid at same k; Eb²·k ≈ 0.7; r=0 atom ≈ 0.85/k̄. Ladder ⇒ 0;
nested ⇒ 2^{−r}−1/2; cut chord ⇒ (1/2)bias′−1/4; pendants ⇒
(1−2^{−(m+1)})^r − 1/2. P(S(λ)) ≥ e^{−O(n log n)}; CGW Thm 3.13 merge
ratio ∈ [1/2, 3/2].

*The notes are the living document; keep PROVED/VERIFIED/COMPUTED
labels honest. Session-8 headline: the mixing problem is a marked
split–merge chain with an exact sign symmetry; the 1/d rate is proved
sharp from below; the upper bound is two scoped estimates away, and
the transfer is a de Finetti step away.*
