# HANDOFF: Cameron's second-row conjecture / CGW Conjecture — state after session 10

## >>> SESSION 10 HEADLINE: the environment CANNOT be avoided in
prob:survival (recharge is essential), but it is now CONTROLLED on the
short horizon: the number of curves B_t after t chords has the EXACT
Harer–Zagier law, E[3^{B_t}] = 3+4t+2t², so B_t ≤ C_γ log d for ALL
t ≤ 9d except w.p. d^{−γ} (cor:horizon) — the "no-dust" obstacle is
closed up to logarithms.  Also PROVED: the exact law of the mass of the
curve through a uniform (or fixed) point at every time m,
(1−(1−2x)^m)dx + [(1+(−1)^m)/(2(m+1))]δ₁ (thm:sizebiased; gives
E[Σp^r(m)] = 1/r + ∫₀¹(1−x^{r−1})(1−2x)^m dx for all r — TODO 10b(b)
done and superseded); a recharge lemma under the block bound
(lem:srecharge) and the log-drift of the marked total (lem:logdrift);
the exact Markov-renewal identity Σ_d Φ_d z^d = (1−z)^{−1}(I+K_z)^{−1}
(I−K_z)1 (prop:genfun) and a conditional decay proposition
(prop:condecay): prob:survival ⇐ (R1) uniform τ^{1+γ}-moments of the
first-flip time on the log-macroscopic set + (R2) uniform invertibility
of I+K_z on the closed unit disc (+ an occupation-time step for s).
NUMERICS: the pure shrink–flip chain (no merges) from (1/4,1/4) decays
like ~1/t² to the noise floor, but from (1/4,0.02) and (0.02,0.02) only
on the ε²t scale — proving that recharge by environment blocks is the
only route (§sec:noenv). <<<

*Prepared 2026-09-03 (end of session 10). Self-contained; companion
tarball has all code, data, papers, notes. Previous handoffs valid
except: HANDOFF_session_9 TODO 10a(A)(iii)'s "quantitative convergence
to PD(1) in the large-block topology" is no longer the tool — the
block-count bound replaces it; TODO 10b(b) is done.  `HANDOFF.md` = this
file.*

---

## 0. The problem in one paragraph

As before (HANDOFF_session_4 §0). Weak conjecture d_TV(P_n,Q_n)=o(1)
STILL OPEN; owner directive: CGW Conjecture is the priority. Chain:
thm:marked → programme steps 1–5 → open items now:
  (I)   prob:mixing, REDUCED to prob:survival (any α > 2/9, j ≤ 8d),
        which after §sec:s10 is REDUCED to (R1)+(R2)+occupation step on
        the log-macroscopic set S'_K = {a,b ≥ 1/(16K)}, K = O(log d),
        with the environment handled by cor:horizon.
  (II)  transfer iid → orbit (unchanged: mixture step + θ-uniform
        survival).
  (III) prob:decouple (unchanged).

## 1. What session 10 added (all in `section_mixing_s10.tex`, §sec:s10)

### §sec:noenv — the environment cannot be avoided (COMPUTED + trivial scaling)
`s10_nomerge a0 b0 tmax S`: the (a,b) chain with merges suppressed.
From (1/4,1/4): Φ = .130/.031/.0070/.0024 at t=8/16/32/64, t²Φ ≈ 8–10,
then noise floor 1e−3 (2e6 samples).  From (1/4,0.02): Φ = .53/.25/.088/
.040 at t=64/256/1024/4096.  From (0.02,0.02): .90/.68/.30/.13 at
t=64/256/1024/8192.  Without merges only ε²t matters (all rates
quadratic; the big coordinate decays like t^{−1/2}), and the budget of
cor:kappafree forces ε ≪ d^{−1/2}, so ε²d ≪ 1: any proof MUST use
merges with environment blocks (recharge).  Logs: s10_nomerge_*.log.

### thm:blockcount (PROVED + VERIFIED) — B_t is Harer–Zagier
The 2t chord endpoints are iid uniform ⇒ a uniform chord diagram ⇒ the
number of curves is the HZ vertex count of a random gluing of a 2t-gon:
  Σ_t E[x^{B_t}] z^t = [((1+z)/(1−z))^x − 1]/(2z),
  P(B_t = m) = [z^{t+1}] (2 artanh z)^m/(2 m!),  B_t ≡ t+1 (mod 2),
  E[2^B] = 2+2t, E[3^B] = 3+4t+2t², E[4^B] = (4/3)(t+1)(t²+2t+3),
  P(B_t = 1) = (1+(−1)^t)/(2(t+1)),  E[B_t] = log t + O(1).
Second proof via hooks: x^{c(σ)} = Σ_λ s_λ(1^x)χ_λ; s_{hook}(1^x) by
hook-content; for integer x only k ≤ x−1 hooks contribute
(`s10_blocks.py n tmax xmax`, exact rationals, extrapolated in 1/(n−1)).
MC check `s10_blocksim.py 6 200000`: all P(B_t=m) agree to ±0.001.

### cor:horizon (PROVED) — block bound over the horizon
P(max_{t≤T} B_t ≥ K) ≤ 3^{1−K}(T+1)³.  With K_γ(T) = ⌈(3+γ)log₃(T+1)⌉+1
the event H_T = {B_t < K_γ ∀t ≤ T} has P ≥ 1 − (T+1)^{−γ}; on it
Σp_i² ≥ 1/K, env has < K blocks, largest env block ≥ (1−s_t)/K.  For
prob:survival (μ_j start, j ≤ 8d, then d steps) use T = 9d.  This is
the no-dust estimate in the form the short horizon needs (at
stationarity B = ∞ and the statement is false — the short horizon is
again what makes the iid model tractable).

### thm:sizebiased (PROVED + VERIFIED) — law of the block of a uniform point
For F on {1..n}: (−1)^k⟨χ_{(n−k,1^k)}, Σ_u F(|C(u)|)⟩ = F(n) − F(k)
(1 ≤ k ≤ n−1), Σ_{j≤n}F(j) (k=0).  Hence P_m(|C(1)| = j) = (1−r_j^m)/n
for j < n, P_m(|C(1)| = n) = (1/n)Σ_k r_k^m, r_k = 1−2k/(n−1).
Continuum: density 1 − (1−2x)^m on (0,1), atom (1+(−1)^m)/(2(m+1)) at 1
(= P(B_m = 1), consistent).  Moments E[Σp^r(m)] = 1/r + ∫₀¹(1−x^{r−1})
(1−2x)^m dx, all r (r=2 is thm:exactonepoint).  Hook identity
(−1)^k⟨χ_k, Σc^r⟩ = n^{r−1} − k^{r−1} checked by brute force r ≤ 6,
n ≤ 9 (`s10_hooks34.py`); law + moments vs MC (`s10_sizebiased.py m S`)
agree to 3 digits.  Closes the drift identity E[ΔL] = 2L² − (7/3)Σp⁴
in expectation.  Small marked blocks: P(mass of u's curve ≤ η) ≤
min(η, mη²).

### lem:srecharge, lem:logdrift (PROVED)
On H with K ≥ 3: from s_0 ≥ ε ≤ 1/(2K), P(max_{t≤T} s_t < 1/(2K)) ≤
e^{−εT/K} + 2ε²T + P(H^c); per coordinate the same with e^{−εT/(3K)}
(merge with the largest env block if s ≤ 1/2; flip landing ≥ 1/(2K)
w.p. ≥ 1/2 if s > 1/2).  Log-drift: E[Δ log s] ≥ (s/K)log(1+1/(2Ks))
− s²/2 > 0 for s ≤ 1/(2K), light-tailed down-jumps.  The occupation-
time consequence (s ≥ 1/(4K) for ≥ half the times of a window, w.h.p.)
is standard but NOT written out.

### prop:genfun, prop:condecay (PROVED as reductions; (R1),(R2) OPEN)
K_zF(x) = E_x[z^τ F(Γ_τ)], ‖K_z‖ ≤ |z|.  Σ_d Φ_d z^d = (1−z)^{−1}(I+K_z)^{−1}
(I−K_z)1 (|z|<1; Neumann series).  G(z) := (I+K_z)^{−1}(I−K_z)1 has
G(1) = 0 (τ < ∞ a.s. on {ab>0}, taken for granted).  If (R1) z ↦ K_z is
C^{1,γ} on the closed disc on a weighted sup-norm space B (suffices:
sup_x E_x[τ^{1+γ}V(Γ_τ)]/V(x) < ∞) and (R2) sup_{|z|≤1}‖(I+K_z)^{−1}‖_B
< ∞, then sup_{x∈S'}|Φ_d(x)| ≤ C d^{−γ}.  (R1) = traps (measured tail
1/t² ⇒ any γ < 1); (R2) = cancellation, needed on the WHOLE circle
(z=1: −1 ∉ spec K_1 with uniform inverse; z ≠ 1 on the circle:
non-lattice condition; z = −1: parity of inter-flip times not
asymptotically deterministic — where the (−1)^d alternation lives).
Doeblin on the fibre: two consecutive flips give a'' ≥ (1/4)·Unif[s/4,
3s/4] in law; the base (s, env) is not regenerated, so B must be a
Meyn–Tweedie weighted space in which the post-flip chain is
V-geometrically ergodic; on H the base has < K blocks.

Reduction bookkeeping (§sec:markovrenewal, last paragraph): S'_K =
{a,b ≥ 1/(16K)}; enter it via lem:srecharge (s reaches 1/(2K) in
T_1 = O(ε^{−1}K log d)) + occupation step + a balanced flip (rate
≥ ε/(4K) on times with s ≥ 1/(4K), lands min ≥ s/4 w.p. ≥ 1/2), total
failure O(d^{−γ}) + O(εK log d); then strong Markov.  NOTE (referee):
applying lem:srecharge to a and b separately does NOT give a common
time (a coordinate at 1/(2K) shrinks on scale K² ≪ ε^{−1}K).

## 2. File manifest (additions)

| file | what | re-verify |
|---|---|---|
| `section_mixing_s10.tex` | notes §sec:s10 (input by section_mixing.tex after s9) | pdflatex second_row_notes.tex twice (0 errors) |
| `s10_nomerge.c` | (a,b) chain without merges; prints Φ_t at powers of 2 | `gcc -O3 -o s10_nomerge s10_nomerge.c -lm && ./s10_nomerge 0.25 0.25 64 200000` → Φ_8 = 0.13±0.002, Φ_64 ≤ 0.006 (noise floor 2e−3 at this sample size; the quoted 0.0024 used 2e6 samples) |
| `s10_blocks.py` | exact pgf E[x^{B_t}] via hook Schur values, continuum by extrapolation | `python3 s10_blocks.py 60 4 3` → 2,4,6,8,10 and 3,9,19,33,51 |
| `s10_blocksim.py` | MC of split–merge block count vs HZ law | `python3 s10_blocksim.py 4 100000` (agreement ±0.003) |
| `s10_hooks34.py` | hook inner products ⟨χ_k, Σc^r⟩, r=2,3,4 (imports s9_hooks) | `python3 s10_hooks34.py 8 2` → rows n^{r−1}−k^{r−1} |
| `s10_sizebiased.py` | MC of size-biased block law + moments vs exact | `python3 s10_sizebiased.py 8 100000` |
| `s10_nomerge_*.log` | no-merge numerics quoted in §sec:noenv | — |

Notes edits: section_mixing.tex (\input s10; intro sentence);
second_row_notes.tex §sec:open (one clause).

## 3. PRIORITIZED TO-DO (session 11)

**QUICKSTART (≈3 min, 2 cores): as session 9 (`python3 pf_cl.py all`
→ 0 failures, SPLIT conv0 122 BY DESIGN; `python3 pf_master.py 5 8` → 0
failures; `gcc -O3 -o s8_splitmerge s8_splitmerge.c -lm && ./s8_splitmerge
psi 64 16 6000 50 992` → psi=0.034573; `gcc -O3 -o s9_survival
s9_survival.c -lm && ./s9_survival flux 8 200000` → ratio ≈1.0±0.1 for
b ≥ 8; `python3 s9_hooks.py 12 6` → 1.000000/0.638889/0.671717)
PLUS `python3 s10_blocks.py 60 4 3` (2,4,6,8,10 / 3,9,19,33,51) and
`python3 s10_blocksim.py 4 100000`, then `pdflatex second_row_notes.tex`
twice (0 errors).**

**THE EXACT TARGET (TODO 10a step 0, now fixed).** Prove: there are
C, γ > 4/9 and K_0 such that for all K ≥ K_0, all states x on the event
H (fewer than K blocks, maintained over the horizon) with
a, b ≥ 1/(16K), and all d ≥ 1,
    |Φ_d(x)| = |E_x[(−1)^{N_d}]| ≤ C K^{c} d^{−γ}        (S)
(any fixed c). By the reduction paragraph at the end of §sec:markovrenewal
(plus the occupation-time step, TODO 11c), (S) gives prob:survival with
α = γ/2 − o(1) > 2/9, hence prob:mixing, hence the iid model. Sufficient
for (S): (R1)+(R2) of prop:condecay on a weighted sup-norm space. Do
NOT try to prove prob:survival for arbitrary environments — §sec:noenv
shows it is false without recharge, and cor:horizon is what makes H
available for free.

**OWNER'S CONCERN AND THE STOPPING RULE FOR SESSION 11 (read this
first; it applies to all future sessions).** The project owner asked at
the end of session 10 whether the residue is really shrinking or merely
being renamed each session. Honest ledger: distinct obstacles HAVE been
removed (EQ → one marked bit, s5; annealed bias → one L² statement, s7;
sign structure, s8; gap law, s9; environment/no-dust, s10). But the
parity-cancellation core has now been FORMULATED THREE TIMES WITHOUT A
PROOF STEP ON IT: "state-dependent inter-flip times" (s8), the spectral
Hardy inequality + coupling obstruction (s9), and condition (R2) (s10).
These are the same object under three names. Session 11's job is (R2)
ITSELF, and it is to be judged by exactly one of these outcomes:

  (A) a PROOF of (R2) at least in the frozen-environment warm-up of
      TODO 11b (fixed block list, merges from it, fibre Doeblin only),
      written in the notes with an honest label; or
  (B) a clean NEGATIVE finding: the cheap numerics of TODO 11d show
      E_x[(−1)^{τ_k}] or the inter-flip parity NOT decaying/mixing as
      (R2) predicts, or the frozen-environment case provably fails —
      in which case say so and STOP working this route.

A fourth reformulation of the cancellation core, however elegant, does
NOT count as progress and must not be the session's deliverable. If by
mid-session neither (A) nor (B) is in sight, write down precisely what
blocks the frozen-environment case and hand off; do not open new
sub-projects (no exact-formula excursions, no transfer/decouple work,
no long numerics) until (A) or (B) is settled. Future sessions: keep
this ledger honest — every handoff must state whether the core moved
or was renamed.

**START HERE.** Read §sec:s10 (~6 pages). The iid model now hangs on
(R1)+(R2) of prop:condecay on the log-macroscopic set, plus the
occupation-time step for s (lem:logdrift ⇒ standard). Everything about
the environment is done by cor:horizon. Do NOT re-derive no-dust
estimates; do NOT run long numerics (2 cores).

### TODO 11a — (R1) (only AFTER (A) or (B) of the stopping rule is
settled; (R1) is trap accounting with proved ingredients and is the
lower-risk piece): the first-flip tail on S'_K (route B(i), now with
a healthy environment). Target: sup_{x∈S'_K} P_x(τ > t) ≤ C(K)/t²
(polynomial in K), or any t^{−(1+γ)} with γ > 4/9 in the weighted
sense sup_x E_x[τ^{1+γ}V(Γ_τ)]/V(x) < ∞. Ingredients all PROVED:
trap-creation flux 4η²/step (lem:recharge iii); a trap of depth g is
left by a flip alone at rate 2g(s−g); s ≥ 1/(4K) most of the time
(lem:logdrift + occupation step — WRITE THIS FIRST, it is standard:
Lyapunov s^{−θ}); merges recharge (lem:srecharge). Assemble as a
subcritical renewal (kernel ≈ 2/(s²m²), §sec:parity(c)). Half a day.

### TODO 11b — (R2): invertibility of I + K_z. THIS IS THE SESSION
(see the stopping rule above). Do 11d (30 min) FIRST as a falsifier,
then the frozen-environment warm-up, then the real base. Steps: (i) write the
post-flip chain on H as a Markov chain on (a,b; s, env) with < K
blocks; (ii) Doeblin on the fibre (two flips ⇒ (1/4)Unif[s/4,3s/4]);
(iii) V-geometric ergodicity of the base (s, env) between flips —
the env is split–merge with < K blocks, mixing in O(K²) steps; find V;
(iv) −1 ∉ spec(K_1) (no sign-alternating invariant function: uses
that τ is not a.s. even); (v) the circle: non-lattice (τ has a
geometric-like component from any state, so E_x[z^τ ·] has modulus
< 1 on "most" of the state, but uniformity is the issue); z = −1
separately. If (iii) is too heavy, first prove (R2) with the env
FROZEN (adversarial but fixed env, merges from a fixed block list) as
a warm-up — the fibre Doeblin alone should give it, and it would show
the cancellation mechanism is real.

### TODO 11c — the occupation-time lemma for s (needed by both 11a and
the reduction paragraph): from lem:logdrift, P(#{t ≤ T: s_t < 1/(4K)}
≥ T/2 | s_0 ≥ 1/(2K), H) ≤ C K^c/T^{?}; polynomial is enough.

### TODO 11d — sanity numerics for (R2) (cheap, ≤ 30 min, on the s8
kernel in s9_survival.c; add a mode "kflip"). From a macroscopic start
x (e.g. stationary env, a,b ∈ [1/8,3/8]), record the flip times
τ_1 < τ_2 < … and print, for k = 1..16: (i) E_x[(−1)^{τ_k}] — this is
(K_{−1}^k 1)(x), and (R2) at z = −1 predicts geometric decay in k;
(ii) E_x[(−1)^{τ_k − τ_{k−1}}] (parity of inter-flip times, should be
bounded away from ±1); (iii) the tail P(τ_k − τ_{k−1} > t) at powers
of 2 for k = 1, 2, 4 (should all be ~1/t², i.e. (R1) holds uniformly
along the post-flip chain, not just at the start). ~1 min at
cs = 100000.

### TODO 11e — (u,v)-averaged ensemble: check whether thm:factor needs
antipodal ψ or the (u,v)-averaged ψ̄ (the "true ensemble" has random
marks); if ψ̄ suffices, ψ̄(j,d) = Σ_k (−1)^k r_k^j ⟨χ_k, h_d⟩ with
h_d(σ) = E_{A,A'}[Σ_{C,C'}|C∩C'|(|C∩C'|−1)]/(n(n−1)) a CLASS function
of σ (two d-chord sets sharing the prefix) — the j-dependence is exact
in hooks; only the d-dependence of ⟨χ_k, h_d⟩ (a stationary
expectation) remains. Session 9 said the technique "does not reach ψ";
re-examine with this formulation (1 h).

### TODO 10c/10d (transfer, decouple): unchanged.

### TODO 11f (referee): §sec:s10 was refereed in-session (one slip
in E[4^B] and the simultaneous-recharge gap were caught and fixed); the
occupation-time step and "τ < ∞ a.s." are flagged as taken for granted.

## 4. Pitfalls (all previous, plus)

41. Without merges, the parity decays only on the ε²t scale; never
    bound the environment from above only.
42. lem:srecharge applied to a and b separately gives large a and large
    b at DIFFERENT times; entering S'_K needs the balanced-flip step.
43. E[4^{B_t}] = 4 + (20/3)t + 4t² + (4/3)t³ (an earlier fit had
    (26/3)t + 2t²; check against t=2: 44).
44. cor:horizon is about the chain from ONE block; at stationarity
    B = ∞ and the block bound is false. It is a short-horizon fact.
45. (R2) is needed on the whole unit circle, not only near z = 1.
46. `s10_blocksim.py` and `s10_sizebiased.py` are pure Python (≈1 s per
    10⁵ trials at t ≤ 8); do not push to t ≥ 64.

## 5. Papers on hand / to request

IN THE TARBALL: CGW_2008.pdf, KPS_2025.pdf, KSS_large_deviations.pdf,
Schramm_2005.pdf, Kwan_Sudakov_2018.pdf + notes PDF. Keep packing.
To request if needed: Alsmeyer, "The Markov renewal theorem and related
results" (for (R2)/Markov renewal with polynomial rates); Meyn–Tweedie
ch. 15–16 (V-geometric ergodicity); Harer–Zagier 1986 / Zagier 1995
(block-count law, standard); Diaconis–Freedman 1980 (transfer).

## 6. Quick reference

Master formula: unchanged. Factorization: 4Ebias² ≤ max_{j≤k, d≥k/8}
ψ(j,d) + 2e^{−k/32}; ψ(j,d) = E_{μ_j}[Φ_d²]; Φ_d = E[(−1)^{N_d}]. Short
horizon: P_{μ_j}(g ≤ ε) ≤ 4ε²j. Budget: prob:survival with α > 2/9.
Block count: E[3^{B_t}] = 3+4t+2t², P(max_{t≤T}B_t ≥ K) ≤ 3^{1−K}(T+1)³.
Size-biased block at time m: density 1−(1−2x)^m, atom (1+(−1)^m)/(2(m+1)).
E[Σp^r(m)] = 1/r + ∫(1−x^{r−1})(1−2x)^m dx. Generating function:
Σ_dΦ_dz^d = (1−z)^{−1}(I+K_z)^{−1}(I−K_z)1. No-merge chain: decay on ε²t.
π: P(tog)=1/2; π(g≤δ) = 4δ+O(δ²); ψ_π(d) ≥ c/d (PROVED), ≈1.0/d.

*The notes are the living document; keep PROVED/VERIFIED/COMPUTED labels
honest. Session-10 headline: the environment is controlled on the
horizon by the Harer–Zagier block count; the survival problem is a
Markov-renewal problem with two named conditions on the log-macroscopic
post-flip chain.*
