# HANDOFF: Cameron's second-row conjecture / CGW Conjecture — state after session 9

## >>> SESSION 9 HEADLINE: prob:gaplaw is UNNECESSARY; the label is a
parity functional of (t, #env blocks); an EXACT one-point formula.
PROVED (notes §sec:s9 = `section_mixing_s9.tex`, input by
`section_mixing.tex`): (1) prop:shorthorizon — since j ≤ k ≤ 8d in
thm:factor, the state-independent creation flux (lem:recharge iii)
gives P_{μ_j}(g ≤ ε) ≤ 1[g₀^ini ≤ ε] + 4ε²j, VERIFIED sharp (ratio
1.00±0.05 at j=8, ≤0.93 at j=64,512); cor:kappafree — with
prob:survival at exponent α, ψ(j,d) ≤ 32ε²d + C(εd)^{−2α} ≤ C d^{−α/(1+α)},
so ANY α > 2/9 closes the iid model (α=1/2 ⇒ c'=1/6). The no-dust /
escape estimate is demoted to "sharp-exponent" material. (2) lem:parity —
every step is a split or merge, so B_t ≡ B₀+t and
ε_t = −(−1)^{B₀+t}(−1)^{E_t}, (−1)^{N_t} = (−1)^t(−1)^{E_t−E₀}: the
label is a deterministic function of time and the env block count.
Consequences: (a) f(m) ALTERNATES in sign, f(m) ≈ (−1)^m c₁/m² with
m²|f| → ≈0.5 — session 8 (even m only) misread it as the frozen-label
trap term; the non-alternating part is ≲2e−4 at m=8; (b) coupling
obstruction: any coupling proof of prob:survival must produce envs of
different block-count parity; same-time ι-mirroring is impossible; the
tiny-block swap coupling caps at |Φ_d| ≤ 1−c; (c) the stationary
spectral measure has ≲3% mass near λ=−1 (odd/even ratio 0.97 at d=16,64)
— the alternation is a clean-start phenomenon. (3) thm:exactonepoint —
E[Σc_i²(σ_m)] = n(n+1)/2 + Σ_{k=1}^{n−1}(n−k)(1−2k/(n−1))^m for
σ_m = (uniform n-cycle)·(m transpositions), via hook characters
(⟨χ_{(n−k,1^k)}, Σc²⟩ = (−1)^k(n−k), proved by generating function);
continuum: E[Σp_i(m)²] = 1/2 + 1/(2(m+1)) (m even), 1/2 + 1/(2(m+2))
(m odd); the (u,v)-averaged one-point function is 1/(m+1) resp. 1/(m+2),
= (2m+3)/(2(m+1)(m+2)) + (−1)^m/(2(m+1)(m+2)). VERIFIED vs continuum
simulator to 3 digits and by hand at m=1,2. NUMERICS: first-flip tail
from good starts (stationary env, g₀∈[1/16,1/8)) is t²P(τ>t) → 11.1:
a 1/t² law (subcritical renewal, kernel 2/(s²m²)) — the true α in
prob:survival is plausibly 2; E[Φ²|g₀=u/d] is at the noise floor
(≲1e−3) for all u ≥ 3 in the s8 psibar logs. Spectral route made
concrete on the quotient (twisted kernel P_π = P_nf − P_f, self-adjoint,
Φ_d = P_π^d 1); a Cauchy–Schwarz/segment-smoothing sketch would give
ν([1−δ,1]) ≤ C√δ ⇒ ψ_π(d) ≤ Cd^{−1/2} (c=1/4) AT STATIONARITY; the
passage π → μ_j is the real obstacle for the iid model. <<<

*Prepared 2026-09-03 (end of session 9). Self-contained; companion
tarball has all code, data, papers, notes. Previous handoffs valid
except: HANDOFF_session_8 §1 "(d) one-point exponent ... confirms the
trap-creation mechanism" is INCOMPLETE (see (2a) above); its TODO 9b
(no-dust / κ=1 gap law) is no longer on the critical path (cor:kappafree).
`HANDOFF.md` = this file.*

---

## 0. The problem in one paragraph

As before (HANDOFF_session_4 §0). Weak conjecture d_TV(P_n,Q_n)=o(1)
STILL OPEN; owner directive: CGW Conjecture is the priority. Chain:
thm:marked → programme steps 1–5 → open items now:
  (I)   prob:mixing, REDUCED to prob:survival ALONE (any α > 2/9)
        over the horizon j ≤ 8d; rate 1/d proved sharp from below.
  (II)  transfer iid → orbit, REDUCED to mixture/de Finetti step +
        θ-uniform prob:survival (prop:shorthorizon is θ-covariant).
  (III) prob:decouple (unchanged; s8 empirics supportive).

## 1. What session 9 added (all in `section_mixing_s9.tex`)

### prop:shorthorizon + cor:kappafree (PROVED + VERIFIED)
P_{μ_j}(g_j ≤ ε) ≤ 1[g₀^ini ≤ ε] + 4ε²j (each of the 4 channels enters
(0,ε] w.p. ≤ ε² per step from ANY state). With j ≤ 8d:
ψ ≤ 32ε²d + C(εd)^{−2α}, optimum ε = d^{−(1+2α)/(2+2α)},
ψ ≤ C d^{−α/(1+α)}, c' = α/(2(1+α)) > 1/11 ⟺ α > 2/9.
`s9_survival flux j cs` prints P(g≤2^{−b}) vs 4ε²j (logs s9_flux{64,512}.log).

### lem:parity (PROVED)
ε_t = −(−1)^{B₀+t}(−1)^{E_t}; Φ_d = (−1)^d E[(−1)^{E_d−E₀}]. Smooth
object = (−1)^dΦ_d. One-point at odd m (s9_onept_odd.log): f(7)=−.00753,
f(8)=+.00599, f(9)=−.00526, …, f(17)=−.00149 (se 2.2e−4).

### thm:exactonepoint (PROVED + VERIFIED)
`python3 s9_hooks.py n M` → exact E[Σp²] at finite n (rationals);
n=40: .6007/.5567/.5314 at m=4/8/16 vs continuum .6005/.5561/.5297.
Continuum closed form above. Spectral density of the averaged start is
(1+y)/4 on y=1−2x∈[−1,1] — linear at BOTH ends (trap end y=+1,
alternating end y=−1). Same technique does NOT reach ψ (two-time
quantity is not a class function of the increment) — checked.

### First-flip tail (COMPUTED)
`s9_survival phic lo hi d cs`: stationary start, accept g₀∈[lo,hi),
reports Φ_d, P(τ>d), tail at powers of 2. g₀∈[1/16,1/8): t²P(τ>t) =
29.0/16.1/11.3/11.1 at t=32/64/128/256. `s9_survival tail g0 tmax cs`
(clean start, env = single block): P(τ>t) ~ 1/t². PITFALL: the `phi`
mode (force_gap) is CONFOUNDED when s < g₀ (gap left unforced); use
`phic`.

### Odd/even stationary correlations (COMPUTED)
`s9_survival oddeven d cs fs`: ⟨P^dε,P^{d+1}ε⟩/⟨P^dε,P^dε⟩ = 0.971 (d=16),
0.974 (d=64): no λ≈−1 mass at stationarity beyond ~3%.

### Spectral route, twisted kernel (SKETCH, §sec:twisted)
Gaps: (i) weighted Poincaré for Q² on a segment with weight a'(s−a');
(ii) the −1 end (numerically absent); (iii) π → μ_j, j ≤ 8d — the real
obstacle. Idea for (iii): dust-insensitive comparison of μ_j (j ≥ d/2,
WLOG by lem:freeze) with π: a block of mass w influences d steps w.p.
≤ 2wd.

## 2. File manifest (additions)

| file | what | re-verify |
|---|---|---|
| `section_mixing_s9.tex` | notes §sec:s9 (input by section_mixing.tex) | pdflatex second_row_notes.tex (clean, 0 errors) |
| `s9_survival.c` | quotient-chain diagnostics: flux / tail / tailenv / phi (confounded) / phic / oddeven modes; flip counter NFLIP on s8 kernel | `gcc -O3 -o s9_survival s9_survival.c -lm && ./s9_survival flux 8 200000` (ratios ≈1 at b≥8) |
| `s9_hooks.py` | exact hook-character formula for E[Σc²], finite n | `python3 s9_hooks.py 12 6` → m=0: 1.000000, m=2: 0.671717 (=1/2+1/6+O(1/n)) |
| `s9_flux64.log`, `s9_flux512.log` | flux-bound checks | — |
| `s9_phic.log`, `s9_tail.log`, `s9_oddeven.log`, `s9_onept_odd.log`, `s9_hooks40.log` | session-9 numerics | — |

Notes edits: second_row_notes.tex §sec:open (gaplaw no longer needed);
section_mixing.tex intro, §sec:transfer (θ-uniform survival only),
§sec:mixnum(d) (alternation correction).

## 3. PRIORITIZED TO-DO (session 10)

**QUICKSTART (≈3 min total): `python3 pf_cl.py all` (0 failures; SPLIT
conv0 122 BY DESIGN), `python3 pf_master.py 5 8` (0 failures),
`gcc -O3 -o s8_splitmerge s8_splitmerge.c -lm && ./s8_splitmerge psi 64
16 6000 50 992` → psi=0.034573 EXACTLY, `gcc -O3 -o s9_survival
s9_survival.c -lm && ./s9_survival flux 8 200000` (ratio column ≈1.0±0.1
for b ≥ 8), `python3 s9_hooks.py 12 6` → 1.000000 / 0.638889 / 0.671717
/ 0.585170 / … , `pdflatex second_row_notes.tex` twice (0 errors).**

**START HERE (first hour).** Read §sec:s9 of the notes
(`section_mixing_s9.tex`, ~5 pages): prop:shorthorizon, cor:kappafree,
lem:parity, thm:exactonepoint, §sec:twisted. Then do TODO 10a step 0
below before anything else — it fixes the exact statement to prove.
Do NOT launch long numerics (no d=2048, no n=200 reruns); the machine
had 2 cores in session 9 and everything below is math-first.

### TODO 10a (mathematical core, the ONLY thing standing between the
iid model and a proof): prove prob:survival with ANY α > 2/9.

**Step 0 — write the exact target.** After cor:kappafree the statement
needed is: there are C, α > 2/9 such that for all d, all j ≤ 8d, all
ε ∈ (0,1/4):
    E_{μ_j}[ Φ_d(Γ_0)² ; g_0 ≥ ε ] ≤ C (εd)^{−2α}.
Since |Φ_d| ≤ 1 it is ENOUGH to bound E_{μ_j}[|Φ_d|; g_0 ≥ ε], and one
may additionally discard the s-trap {s_0 ≤ σ} at cost 2σ²j ≤ 16σ²d
(entry into {s ≤ σ} needs a shrink landing a′ ≤ σ − b with b ≤ σ:
probability ≤ σ² per channel, two channels — lem:strap, PROVED in
§sec:s9).
So WLOG a_0, b_0 ≥ ε and s_0 ≥ σ with σ ≍ ε; the budget then reads
ψ ≤ 32ε²d + 16σ²d + C(εd)^{−2α}.

**Step 1 — choose a route.** Recommendation: (A) first, because its
first two gaps are concrete and finite; fall back to (B) if (iii)
looks hopeless after a day.

(A) Spectral, §sec:twisted. Objects: quotient chain, twisted kernel
P_π = P_nf − P_f (self-adjoint on L²(π̃)), Φ_d = P_π^d 1, Dirichlet
form E(F) = ½E[(F−F′)²; nf] + ½E[(F+F′)²; f], flip rate k = 2ab,
normalised flip kernel Q with trapezoid density on the segment
{a′ + b′ = s, env fixed}.
  (i) Prove the 1-D smoothing inequality: for F on (0,s),
      ‖F − Q²F‖²_{w} ≥ c·inf_κ ‖F − κ‖²_{w′} with weights w, w′ ∝
      a′(s−a′) (explicit kernel, half-day; check the constant
      numerically by discretising Q on a grid).
  (ii) Show ν (spectral measure of 1) has no mass in [−1, −1+δ] beyond
      Cδ^{1/2}: either prove E[(−1)^{E}] = 0-type orthogonality in the
      continuum, or run the same argument for ⟨F,(I+P_π)F⟩. Numerics
      say the mass is ≲3% (s9_oddeven.log); confirm with more
      statistics if needed (`./s9_survival oddeven 64 20000 400` ≈ 4 min).
  (iii) π → μ_j (THE CRUX). Everything above gives ψ_π(d) ≤ Cd^{−1/2}
      (c = 1/4 at stationarity). For μ_j: WLOG j ≥ d/2 by lem:freeze
      (ψ(j,d) ≤ ψ(j+d/2, d/2)). Proposed tool: a coupling of μ_j with π
      in which all blocks of mass ≥ w agree except with probability
      p(j,w), plus the stability bound |Φ_d(x) − Φ_d(x′)| ≤ 2d·(mass of
      disagreeing material) (a block of mass w is hit in d steps w.p.
      ≤ 2wd). Needs a quantitative convergence of uniform split–merge
      from one block to PD(1) in the large-block topology —
      Schramm_2005.pdf (packed) is the reference; the exact formula
      thm:exactonepoint (E Σp² → 1/2 at rate 1/(2(m+1))) and its
      higher-moment analogues (TODO 10b(b)) are the first quantitative
      handles. If (iii) fails, the stationary result alone is still
      worth writing up as a theorem (prob:mixing at π with c = 1/4).

(B) Renewal, §sec:parity(c). Objects: first-flip time τ, flip times
τ_1 < τ_2 < …, Φ_d = Σ_k (−1)^k P(N_d = k).
  (i) Prove the first-flip tail from good starts:
      P_Γ(τ > t) ≤ C( e^{−c εσ t} + 1/(σ² t²) ) for a_0, b_0 ≥ ε,
      s_0 ≥ σ. NOTE the 1/t² term is ε-INDEPENDENT (traps at depth
      ~1/t are created from macroscopic states at flux 8η dη and
      survive m steps w.p. ≈ e^{−2ηsm}: kernel K(m) ≈ 2/(s²m²), total
      mass ≈ 4ε/s < 1 — subcritical renewal). Measured: t²P(τ>t) →
      11.1 for g₀ ∈ [1/16,1/8), s ≈ 0.6 (s9_phic.log). Ingredients
      already proved: flux bound lem:recharge(iii), freeze argument of
      prop:lower. Missing: the escape hazard ≥ 2g(s−g) via flips (no
      environment needed! a flip needs only a cut in each marked
      coordinate) — this makes the "no-dust" issue disappear for the
      g-trap; the s-trap needs a merge partner but costs only σ²j.
  (ii) Parity cancellation with state-dependent inter-flip times. The
      obstruction lem:parity(b): NO same-time coupling can work; the
      cancellation must be extracted from the flip geometry. Concrete
      idea: post-flip states (α+β, s−α−β) have a density bounded below
      on the middle of the segment; two consecutive flips from any
      state give a state whose law dominates c·(uniform on the middle
      half of the segment) — a Doeblin minorisation for the post-flip
      chain, which makes the inter-flip times "nearly iid" after two
      flips; then apply the alternating-renewal bound for a regenerative
      process. Aim for α = 1/2 (only α > 2/9 is needed).

**Step 2 — write it up in the notes with honest labels**, and rerun
the quickstart before packing.

### TODO 10b: exact-formula extensions (cheap, HZ-adjacent, ≤ 2 h each).
(a) Antipodal f(m): conjecture m²f(m) → (−1)^m/2. Compute f(m) exactly
for m ≤ 6 by symbolic integration over the split–merge tree (f(1) = 0,
f(2) = 1/12 by hand; f(2) = 1/6 − 1/12 from the shrink and flip
branches) and fit the tail form (−1)^m/(2(m+a)(m+b)).
(b) E[Σp³], E[Σp⁴] along the chain via hook characters of Σc³, Σc⁴
(same s9_hooks.py machinery: replace F(mu)=sum mu_i^2; the generating
function with one marked part is (1+t)^{−1}[Σ_j j^{r−1}(1−(−t)^j)x^j]
·(1+tx)/(1−x)). Gives the exact law of the Simpson index and closes
the drift identity E[ΔL] = 2L² − (7/3)Σp⁴ (L = Σp²) — a first
quantitative handle for 10a(A)(iii).

### TODO 10c: transfer (unchanged from 9c): (a) exchangeable→mixture
step for the greedy family (Diaconis–Freedman, request the paper);
(b) θ-uniform prob:survival only (prop:shorthorizon is θ-covariant).

### TODO 10d: prob:decouple (unchanged from 9d).

### TODO 1⁗⁗ (referee): §sec:s9 is new — referee it, especially the
generating-function step in thm:exactonepoint (one marked part; the
[x^n] extraction) and the ±1 accounting in lem:parity.

### Numerics to have at hand (all already in logs; rerun only if needed):
- E[Φ_d² | g₀ bin] comes from `s8_splitmerge psibar` output (column
  "psi|bin"); `s9_survival phic` gives E[Φ_d | bin] and the τ-tail, NOT
  Φ². Runtimes (2 cores): `phic 0.0625 0.125 256 100000` ≈ 3 min;
  `tail 0.25 4096 200000` ≈ 1 min; `flux 512 400000` ≈ 6 min;
  `oddeven 64 20000 400` ≈ 4 min; `s9_hooks.py 40 16` ≈ 3 s.

## 4. Pitfalls (all previous, plus)

37. The `phi` mode of s9_survival.c (force_gap) leaves the state
    unchanged when s < g₀ — its output at macroscopic g₀ is
    contaminated by unforced small-gap states. Use `phic` (rejection).
38. f(m) alternates in sign. Never sample only even m; the smooth
    object is (−1)^m f(m) (clean start) and (−1)^d Φ_d in general. At
    STATIONARITY there is no alternation (ν has no −1 mass).
39. Hook-character convolution: weight per hook is (−1)^k, NOT
    C(n−1,k)(−1)^k (the d_λ cancels against χ_λ(C_n)/d_λ). m=0 must
    give exactly 1.
40. In thm:factor, j ≤ k and d ≥ k/8 ⇒ j ≤ 8d: the short horizon is
    what makes the flux bound sufficient; a "uniform in j" gap law
    was never needed for the iid model.

## 5. Papers on hand / to request

**IN THE TARBALL (unchanged): CGW_2008.pdf, KPS_2025.pdf,
KSS_large_deviations.pdf, Schramm_2005.pdf, Kwan_Sudakov_2018.pdf +
notes PDF. Keep packing all of them. (kps.txt kept.)**

**To request if needed (session 10):** Diaconis–Freedman 1980 (9c);
Diaconis–Mayer-Wolf–Zeitouni–Zerner 2004 / Berestycki EFC (10a-A(iii));
Diaconis–Shahshahani 1981 (random transposition eigenvalues — used in
thm:exactonepoint, standard).

## 6. Quick reference

Master formula: unchanged, VERIFIED to n=100. Factorization:
4Ebias² ≤ max_{j≤k, d≥k/8} ψ(j,d) + 2e^{−k/32}; ψ(j,d) = E_{μ_j}[Φ_d²];
Φ_d = E[(−1)^{N_d}] = (−1)^d E[(−1)^{E_d−E₀}]; flips prob 2ab,
(a,b)→(α+β, a+b−α−β). Short horizon: P_{μ_j}(g≤ε) ≤ 4ε²j (j≤8d).
Budget: prob:survival with α > 2/9 suffices (c' = α/(2(1+α))).
π: P(tog)=1/2; π(g≤δ) = 4δ+O(δ²); ψ_π(d) ≥ c/d (PROVED), ≈1.0/d
(measured); first-flip tail from good starts ~11/t². Exact: E[Σp²(m)]
= 1/2 + 1/(2(m+1)) (even m), 1/2 + 1/(2(m+2)) (odd m). iid 4Eb²·k → ≈1.45
flat; real ensemble Eb²·k ≈ 0.7. P(S(λ)) ≥ e^{−O(n log n)}; CGW Thm 3.13
merge ratio ∈ [1/2, 3/2].

*The notes are the living document; keep PROVED/VERIFIED/COMPUTED
labels honest. Session-9 headline: the iid model now hangs on ONE
estimate (prob:survival, any α > 2/9), whose mechanism is numerically
a 1/t² first-flip tail; the label is a parity functional of the
environment block count; and the (u,v)-averaged one-point function is
known exactly.*
