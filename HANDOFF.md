# HANDOFF: Cameron's second-row conjecture / CGW Conjecture — state after session 12

## >>> SESSION 12 HEADLINE: outcome (B), sharp.  The minorisation route
(TODO 12b, prop:c1, "uniform ergodicity", "V-geometric ergodicity for
ANY V") is CLOSED by a rigorous obstruction: the environment carries
INERT DUST.  A curve of mass ε is hit with probability ≤ 2ε per step,
so it survives the first n flips untouched w.p. → 1; the n-th post-flip
laws from x_ε are supported on the pairwise-disjoint sets {a curve of
mass exactly ε}; a common minorant ν on Ω_K would have ν(S_i) → 1 for
infinitely many sets S_i of which each configuration lies in ≤ K
(prop:nomino).  Only hypothesis: τ_n < ∞ a.s.  TODO 12a's drift
FAILS on the full space (prop:nodrift; fallback = killed chain, as
instructed).  TODO 12c's criterion does not imply (C2) (§sec:s12odd).
AND: (C2) IS NOT GEOMETRIC.  Measured with 10^8 samples from the
one-curve start, a_k = E_x[(−1)^{τ_k}] = (−1)^k·(2.5–3)·k^{−3}, local
exponent 2.7 ± 0.3 over 8 ≤ k ≤ 20; the same universal tail emerges
after ~8 flips from a 9-curve start.  Session 11's ρ ≈ 0.3 was the
transient above the noise floor.  The polynomial decay is a slow mode
at eigenvalue +1 of κ₁ (the block count B grows like log), not a −1
problem; (C1) is untouched.  What the reduction needs at z = −1 is
only C^{0,γ}, γ > 4/9, and the MEASURED modulus of continuity of the
z = −1 sum along the real axis is |1−w|^{0.75} (prop:holder and its
status).  So the programme survives with (C2) replaced by the
polynomial (C2″); but no Doeblin-type proof of it can exist. <<<

*Prepared 2026-09-15 (end of session 12).  Self-contained; companion
tarball has all code, data, papers, notes (83 pp, 0 errors).
Previous handoffs valid except: HANDOFF_session_11's (C2) is geometric
only in the transient; its TODO 12a/12b/12c are all closed (negative);
its "route (a) uniform ergodicity + parity symmetry" is impossible.
`HANDOFF.md` = this file.*

---

## 0. The problem in one paragraph

As before (HANDOFF_session_4 §0).  Weak conjecture STILL OPEN; owner
directive: CGW Conjecture is the priority.  Chain: thm:marked →
programme steps 1–5 → open items:
  (I)   prob:mixing → prob:survival → (S) on S'_K, which after
        §sec:s11 + §sec:s12 is REDUCED to (R1) + (C1) + (C2″) [+ the
        occupation-time step, TODO 11c], with the whole circle except
        z = ±1 PROVED (thm:idleloop) and z = −1 reduced to z = +1 with
        test function (−1)^E (prop:repair).
  (II)  transfer iid → orbit (unchanged, never examined).
  (III) prob:decouple (unchanged, never examined).

## 0.5 Notation (additions to HANDOFF_session_11 §0.5)

| symbol | meaning |
|---|---|
| a_k(x) | E_x[(−1)^{τ_k}] = (κ_{−1}^k 1)(x) = (−1)^{B(x)} Σ_b (−1)^b P_x[B(Γ̄_{τ_k}) = b]: the parity bias of the block count at the k-th flip.  κ₁^k(−1)^E = (−1)^E (−1)^k a_k |
| y_ε, x_ε | a state y with a sub-arc D of mass ε (inside one curve, no mark) declared a separate curve: "y plus a dust curve" |
| S_ε | {configurations with a curve of mass exactly ε} |
| (C2″) | sup_{x∈S'_K} |E_x[w^{τ_k}(−1)^{E_{τ_k}}]| ≤ C(K) k^{−p}, uniformly for |w| ≤ 1 near w = 1 |
| (M_p) | sup E_y[τ_k] ≤ C k^p over the log-macroscopic landing states (implied by (R1) along the chain; measured E τ_k ≈ 6k) |
| t_k(w) | E_x[w^{τ_k}(−1)^{τ_k}] = (κ_{−w}^k 1)(x), the twisted parity; t_k(1) = a_k |

## 1. What session 12 added (all in `section_mixing_s12.tex`, §sec:s12)

### lem:inert — inert dust coupling (PROVED)
Drive the chains from y and from y_ε by the same uniform points.  Until
a point falls in D, the y_ε-chain is the y-chain with D declared a
separate curve; same flip times; the y_ε-chain has a curve of mass
exactly ε.  P(coupling holds through T) ≥ 1 − 2εT.  (Works in the
positional/cyclic-order model; D is a contiguous segment of its curve.)

### prop:nomino — no state-uniform minorisation, for ANY n (PROVED)
For K ≥ 4, no n, δ > 0, ν on Ω_K with P_x[Γ̄_{τ_n} ∈ ·] ≥ δν for all
x ∈ Ω_K ∩ S'_K (killed or unkilled).  Proof: x_i = (1/4,1/4;{1/2−2^{−i},
2^{−i}}); δν(Ω_K \ S_i) ≤ 2·2^{−i}T_i + P_{x_0}[τ_n > T_i] → 0 with
T_i = 2^{i/2}; so ν(S_i) → 1, but Σ_i ν(S_i) ≤ K (≤ K masses per
configuration).  Consequence: the hypothesis of prop:c1 is NEVER
satisfied; prop:c1 is correct and vacuous.  NOT the anticipated
difficulty (cost K^{−2K}); the obstruction is that no finite sequence
of steps removes a curve that is never hit.

### cor:noerg — no uniform / V-uniform ergodicity (PROVED)
(i) No (π̄, C, ρ<1) with sup_{Ω_K}|κ₁^kF − π̄F| ≤ Cρ^k‖F‖_∞ (the
oscillation of κ₁^k 1_{S_i} over Ω_K → 1 for each k).
(ii) Under (M_p): no V ≥ 1 for which κ₁ is V-uniformly ergodic — V
would have to satisfy V(y_ε) ≥ 1/(2C'ρ^k) for ε ≤ 1/(32Ck^p), making
κ₁V(x_0) = ∞ via "carve dust at step 1, flip at step 2".  For the
killed chain: constants exponential in K (fatal).  This closes the
sentence in §sec:markovrenewal asking for a V-geometrically ergodic
base and route (a) of §sec:s11left.

### prop:nodrift — TODO 12a fails on the full space (PROVED)
V = (Σr⁴)^{−θ}, θ ≤ 1/2: no κ₁V ≤ ϱV + b.  Family x_{m,N}: apart,
a = m, b = m/10, environment in N ≥ 16(T+1)³ max(m^{−1}, m^{−4/3})
equal pieces.  Coupling with an "isolated" process (dust–non-dust
chords suppressed) gives P[τ_k ≤ T] ≥ 3/4 uniformly; deterministically
Σr⁴ ≤ 17m⁴ up to time T (non-dust mass ≤ 2m, dust curves ≤ m^{4/3});
so κ₁^kV ≥ (3/4)17^{−θ}m^{−4θ} while V ≤ m^{−4θ}: contradiction as
m → 0.  On Ω_K the drift is vacuous (V ≤ (2K)^{4θ}); the killing stays
and is harmless (cor:horizon).  Mechanism verified numerically
(`s12_driftfail`).

### §sec:s12odd — TODO 12c's criterion is insufficient (PROVED, trivial)
Inter-flip times 1,2,1,2,… have opposite consecutive parities always,
yet (−1)^{τ_k} has period 4.  The two-block argument needs
‖κ_{−1}²‖_∞ < 1, false exactly (eigenvector (−1)^E).  Dust insertions
are the only invisible odd perturbations and give 1 − O(η) per flip.

### §sec:s12num — (C2) is polynomial (MEASURED, decisive)
`s12_akdecay` (2·10^7 samples, se 2.2e−4) and `s12_akB` (10^8 samples,
se 1e−4), fixed deterministic starts, no censoring:
* one curve, arcs (1/2,1/2): a_k = −.4004, +.1530, −.0674, +.0335,
  −.0186, +.0111, −.0073, +.0048, −.0034, +.0022, −.0019, +.0015 (k ≤ 12);
  10^8 run: |a_k| = .00495 .00349 .00271 .00195 .00142 .00115 .00107
  .00095 .00091 .00079 .00055 .00045 .00035 (k = 8..20), sign (−1)^k;
  k³|a_k| ∈ [2.4, 2.7] for k ≤ 13, [2.8, 3.9] for 14 ≤ k ≤ 20 (2–3σ);
  local exponent 2.7 ± 0.3.  Ratios |a_{k+1}/a_k| rise .38 → .70: NOT
  geometric.
* (0.3,0.3;{0.4}) apart: same tail, k³|a_k| → 2.0–2.4, sign (−1)^{k+1}
  = (−1)^{k+B₀+1}.
* dust inserted into either start (mass 1e−4): a_k identical to 2σ
  (lem:inert in action).
* 9-curve start (0.3,0.3; 7×0.057), 10^8 samples: a_k crosses 0 at
  k = 3, |a_k| GROWS to .00145 at k = 7–8, then decays; 17 consecutive
  sign alternations k = 4..20; k³|a_k| rises from 0.5 to 2–3.  The tail
  is UNIVERSAL (C/k³, C ≈ 2.5–3), reached after a start-dependent
  transient.  17-curve start: |a_k| ≤ 7e−4 for k ≤ 16, transient not over.
* stationary (burnt-in) start: a_k(μ) at noise for k ≥ 5 — because
  the mixture averages the sign (−1)^{B(x)} of eq:akB; sign-blind
  E_μ[a_k²] = (4±5)e−4, consistent with the fixed-start values.  THIS
  IS WHY SESSION 11 SAW ρ ≈ 0.3.
* decomposition by B_{τ_k} (one-curve start): P[B_{τ_k} = 1] ≈ 0.5/k
  (and there τ_k is even, exactly); the law of B_{τ_k} spreads like
  log k; the alternating sum of a lattice law that spreads
  logarithmically decays polynomially.  Mechanism: the block-count
  parity is a SLOW MODE OF κ₁ AT EIGENVALUE +1.  Generating function
  Σ a_k w^k is singular at w = −1 (= eigenvalue −1 of κ_{−1} =
  eigenvalue +1 of κ₁).  Nothing here bears on (C1).

### prop:holder — polynomial decay suffices at z = −1 (PROVED as reduction)
Under (C2″) with exponent p and (M_1), G is C^{0,(p−1)/(p+1)} near
z = −1 (crude bound).  Budget γ > 4/9 ⟺ p > 13/5 on the crude bound.
`s12_twist` (3·10^7 samples, one-curve start, w = .99,.95,.90 real):
(a) exact factorisation t_k(w) = E[w^{τ_k}] a_k FAILS (twisted parity
larger by ×1.2 at |1−w|k = .08, ×2–4 at .4–.8: short paths carry more
memory); (b) uniform bound |t_k(w)| ≤ |a_k| holds at every k, w tested,
so (C2″) holds numerically with the same constant; (c) the SUM
Σ(w) = Σ_k (−1)^k t_k(w) has Σ(1) − Σ(w) = .026, .085, .143 at
|1−w| = .01, .05, .10: modulus |1−w|^{0.75} — far above 4/9, above the
crude 0.5, below the factorised heuristic's min(1,p−1).

## 2. File manifest (additions)

| file | what | re-verify |
|---|---|---|
| `section_mixing_s12.tex` | notes §sec:s12 (input after s11) | `pdflatex second_row_notes.tex` ×2 → 0 errors, 83 pp |
| `s12_state.h` | the state/step code of s11_kflip.c, split out for reuse | — |
| `s12_akdecay.c` | a_k and E_μ[a_k²] from fixed or burnt-in starts, two continuations per start | `gcc -O3 -o s12_akdecay s12_akdecay.c -lm && ./s12_akdecay fixed 0 16 100000 3` → a_1 ≈ −.400, a_4 ≈ +.040 (se .003), 1 s |
| `s12_akB.c` | a_k decomposed by B(Γ̄_{τ_k}); starts 0,1,2,6 (9 curves),7 (17 curves) | `./s12_akB 0 12 3000000 5` → a_8 ≈ +.0047, P[B_8=1] ≈ .063, 1 min |
| `s12_twist.c` | twisted parity t_k(w) and the factorisation check | `./s12_twist 0 8 1000000 41` → t_8(.99) ≈ +.0041 vs a_8 ≈ +.0050 |
| `s12_driftfail.c` | E[(Σr⁴(Γ̄_τ)/Σr⁴(x))^{−θ}] from (m, η; N pieces) | `./s12_driftfail 0.3 0.03 30000 300` → 1.19 (θ=.1) |
| `runs/ak_fixed{0,2,3,4}.log`, `runs/ak_stat_mac.log` | the 2e7 / 5e6 runs quoted | — |
| `runs/akB_{6,7}.log`, `runs/akB_{0,6}_big.log` | the B-decomposition and the 10^8 runs | — |
| `runs/twist_0.log` | the factorisation / modulus run | — |

Notes edits: `section_mixing.tex` (\input s12); `section_mixing_s11.tex`
(SESSION-12 CORRECTION paragraphs at prop:c1 and at the (C2)
measurement).

## 3. PRIORITIZED TO-DO (session 13)

**QUICKSTART (≈6 min, 2 cores; all of it was green at the start of
session 12).**  From the tarball's `second_row_session/`:
```
python3 pf_cl.py all                      # 0 failures; SPLIT conv0 = 122 is BY DESIGN
python3 pf_master.py 5 8                  # 0 failures
gcc -O3 -o s8_splitmerge s8_splitmerge.c -lm && ./s8_splitmerge psi 64 16 6000 50 992   # psi = 0.034573
gcc -O3 -o s9_survival s9_survival.c -lm && ./s9_survival flux 8 200000                 # ratio ≈ 1.0–1.3 for b ≥ 10
python3 s10_blocks.py 60 4 3              # 2,4,6,8,10 / 3,9,19,33,51
python3 s10_blocksim.py 4 100000          # HZ vs MC agree to ±0.003
gcc -O3 -o s11_ztau s11_ztau.c -lm && ./s11_ztau 0.5 1.0 200000 400000 3               # lambda = 0.04153, all rows "ok"
gcc -O3 -o s11_kflip s11_kflip.c -lm && ./s11_kflip kflip 0.125 0.375 16 20000 2000000 7  # parity "HOLDS", a_1 ≈ -0.13
gcc -O3 -o s12_akdecay s12_akdecay.c -lm && ./s12_akdecay fixed 0 16 100000 3          # a_1 ≈ -0.400, a_4 ≈ +0.040 (se 0.003)
gcc -O3 -o s12_akB s12_akB.c -lm && ./s12_akB 0 12 3000000 5                           # a_8 ≈ +0.0047, P[B_8=1] ≈ 0.063
gcc -O3 -o s12_twist s12_twist.c -lm && ./s12_twist 0 8 1000000 41                     # t_8(0.99) ≈ +0.0041 vs a_8 ≈ +0.0050
pdflatex second_row_notes.tex; pdflatex second_row_notes.tex                           # 0 errors, 83 pages
```

**READING ORDER (≈1 h).**  §0.5 of THIS file and of HANDOFF_session_11
(notation); then §sec:s12 of the notes (≈8 pages); then, from earlier
sections, only: prob:survival and the paragraph after it (§sec:mixing),
prop:genfun + prop:condecay + the "Reduction to Problem prob:survival"
paragraph (end of §sec:markovrenewal, s10), lem:parity (s9),
thm:blockcount + thm:sizebiased (section_mixing_s10.tex), thm:factor
(section_annealed.tex) and cor:kappafree (section_mixing_s9.tex).  Do not read
§sec:s11 in full unless you take option (γ).

**THE TARGET, restated.**  (S): |Φ_d(x)| ≤ C K^c d^{−γ}, γ > 4/9, for
x ∈ S'_K ∩ H, K = O(log d).  Sufficient (proved reduction): (R1) +
(C1) + (C2″) + the occupation-time step (TODO 11c).  Proved
unconditionally: everything on the closed disc except z = ±1
(thm:idleloop); z = −1 is z = +1 with test function (−1)^E
(prop:repair).  CLOSED ROUTES — do not reopen, see the DO-NOT list and
pitfalls 54–60: Doeblin minorisation of the post-flip chain in any
form; uniform or V-uniform ergodicity of κ₁ in any weighted sup-norm;
the λ-drift on the full space; the two-block/odd-perturbation argument
in operator-norm form; any argument "uniform over environments".

**STEP 0 (ten minutes, before anything else): TRY TO KILL (C2″).**
Pitfall 53's procedure, now with a second data point: (R2) died by an
eigenvector, geometric (C2) died by a measurement at 100× the
statistics.  For (C2″) the eigenvector test is already done ((−1)^E is
an eigenfunction of κ₁ iff E_x[(−1)^τ] is constant, which it is not).
The remaining cheap kill attempts are: (i) is the exponent p really
> 13/9 (factorised route) or > 13/5 (crude route)?  The measurement
says 2.7 ± 0.3 — pin it with numeric task 3 below before relying on
it; (ii) is the constant C(K) polynomial in K?  Numeric task 2.  If
either fails, record it and the session's finding is that (C2″) is
false; do not then spend the session repairing it.

**OWNER'S CONCERN AND THE LEDGER.**  Session 12's honest entry: the
core did not move toward a proof; it moved toward the truth.  Three
things were established that are not renamings: (a) the whole
Meyn–Tweedie family of approaches is IMPOSSIBLE here, by a two-line
mechanism (inert dust) that no choice of V, ν or n circumvents — the
"base is not regenerated" difficulty of sessions 8–11, finally
identified as structural rather than technical; (b) (C2) as carried
since session 11 is FALSE in its geometric form, and the numerics that
"confirmed" it were reading the transient; (c) the target that
survives is weaker (polynomial (C2″)), is measurably true (p ≈ 2.7 ±
0.3, modulus |1−w|^{0.75} at z = −1), and is exactly what the reduction
consumes.  Against: NO POSITIVE PROOF STEP was made on (C1) or (C2″),
and §sec:s12left says a proof must be "dust-blind" and must use the
environment's actual randomness — the natural transport-metric
coupling is shown NOT to contract, so the tool for a dust-blind proof
does not yet exist.

**SESSION 13: THE THREE OPTIONS, THEIR DELIVERABLES, AND THE STOPPING
RULE.**  The choice between (α), (β), (γ) is the owner's; if the owner
has not said, take (α).  Whichever is taken, session 13 is judged by
exactly one written deliverable, defined below, and by nothing else.
Do not open carried tasks (11a, 11c, 11e, 10c, 10d) to have something
to show.

  (α) THE AUDIT SESSION (recommended first; ≈ one full session).
      Purpose: nothing in the chain thm:marked → programme steps 1–5 →
      prob:survival has been re-checked since it was laid, and two
      carried conditions have now been found false in consecutive
      sessions.  Procedure: for each item below, produce ONE of: a
      proof (written into the notes with a PROVED label), a
      counterexample, or an explicit one-paragraph statement of what
      would make it false and how one would test that.  Items, in
      order of value:
        1. ANNEALED vs QUENCHED — the single most valuable question in
           the project.  Read prob:survival, thm:factor, cor:kappafree
           and the reduction paragraph at the end of §sec:markovrenewal,
           and answer: does the bound the programme needs,
           E_{μ_j}[Φ_d(Γ₀)²; g₀ ≥ ε] ≤ C(εd)^{−2α}, require (S) for
           EVERY environment configuration in S'_K ∩ H, or does it
           suffice to bound E[Φ_d(Γ₀)² | marked pair of Γ₀] with the
           initial environment averaged under μ_j (the two replicas
           still share Γ₀)?  Note that Φ_d(Γ₀)² is a two-replica
           quantity, so "annealed" here means annealed over the shared
           start, not over two independent environments.  If annealed
           suffices, say precisely which statement replaces (S), and
           whether the dust obstruction (which is about uniformity over
           environments) touches it.  Deliverable: a labelled paragraph
           in the notes, "prop:annealedsuffices" or
           "rem:quenchedneeded".
        2. τ < ∞ a.s. on {ab > 0} (taken for granted since s10; used in
           prop:r2false, prop:nomino, prop:nodrift).  Expected: easy
           proof from Σ_t 2a_tb_t = ∞ a.s.; write it.
        3. The occupation-time step (TODO 11c): still unwritten; at
           minimum state it exactly and say what fails if it is false.
        4. The budget arithmetic of cor:kappafree and the short-horizon
           reduction in thm:factor: recompute the exponents by hand,
           including the K^c loss.
        5. (M_1) and (M_p) of §sec:s12 (E τ_k ≤ Ck along the chain):
           state what they need from (R1) and whether TODO 11a gives it.
      STOPPING RULE for (α): item 1 must be answered by the half-way
      point; if it is not, write down exactly where the reduction is
      ambiguous and hand off.  Items 2–5 in the remaining time, in
      order; unfinished items are carried, not faked.

  (β) SCOUTING (two half-sessions, or one if only one branch is
      chosen).  Items (II) transfer iid → orbit and (III) prob:decouple
      have never been examined; their difficulty is UNKNOWN.
      Deliverable per branch: a one-page classification in the notes —
      routine / hard / as hard as (I) — with the first concrete
      sub-problem stated and, if possible, a numerical check that the
      statement is even true.  No proof attempts.

  (γ) CONTINUE ON (I) — only if (α) item 1 says annealed suffices.
      The one route left open by §sec:s12left(iv): prove polynomial
      decay of the block-count parity at flip times from the exact
      fixed-time laws, using a_k = (−1)^{B₀} Σ_b (−1)^b P[B_{τ_k} = b]
      (eq:akB) and treating the flip times as a random sampling of the
      time axis.  Start from the one-curve start, where thm:blockcount
      gives E[x^{B_t}] exactly.  SUCCESS = a_k = O(k^{−p}), p > 13/9,
      from the one-curve start, PROVED, with the argument visibly
      extendable to the μ_j starts.  STOPPING RULE: if no exact
      structure for the flip-time sampling has appeared by the
      half-way point, record what was tried and hand off.

**CHEAP NUMERICS (≈40 min total, 2 cores; do them in the background
during whatever option is taken; all three are in the s12 programs
with trivial edits):**
  1. (C2″) on the unit circle: `s12_twist` with w = e^{iθ}, θ = 0.05,
     0.1, 0.2 (add a complex accumulator, ~10 lines).  Check that
     |Σ(1) − Σ(e^{iθ})| ≍ θ^{0.75} as measured on the real axis.  If
     the modulus on the circle is much worse than on the real axis,
     that is a finding — record it.
  2. K-dependence of the tail constant: `s12_akB` starts 6 (9 curves)
     and 7 (17 curves) to k = 40 with 10^8 samples (≈20 min per start
     per core).  Expected: the same C ≈ 2.5–3 after a transient.  If C
     grows faster than polynomially in the number of curves, (C2″)
     fails the budget — record it.
  3. Pin p: `s12_akB 0 32 400000000` (≈40 min on one core; se 5e−5;
     a_32 ≈ 8e−5 if p = 3).  Fit log|a_k| vs log k over 10 ≤ k ≤ 32.
     Distinguishes p = 2.5 from 3; the budget lines are 13/9 = 1.44
     (factorised) and 13/5 = 2.6 (crude).

**DO NOT, in session 13:**
  1. try to prove the minorisation of prop:c1, in any variant (pitfall 54);
  2. try to prove geometric (C2) (pitfall 55);
  3. look for a Lyapunov function that makes κ₁ V-uniformly ergodic
     (cor:noerg(ii): none exists);
  4. look for an odd perturbation that is exact (§sec:s12odd) — the
     only invisible ones are dust and they are provably too weak;
  5. build a coupling argument on the transport metric between mass
     configurations (§sec:s12left(ii): not contracted);
  6. re-run the s11 stationary-start diagnostic and read ρ ≈ 0.3 from
     it (pitfall 56);
  7. start (γ) before (α) item 1 is answered;
  8. run numerics longer than the three listed above before the
     session's deliverable is settled.

### TODO 11c (carried) — the occupation-time lemma for s.  Unchanged.
### TODO 11a (carried) — (R1).  Unchanged; note it now also supplies
    (M_1), which prop:holder needs.
### TODO 11e, 10c, 10d: unchanged, never opened.

### TODO 12d (referee): §sec:s12 was refereed in-session.  Caught and
fixed: (1) the first draft of prop:nodrift chose T after N although the
law of τ_k depends on N — fixed by coupling with an isolated process
whose τ_k-law is N-free; (2) the first draft of the same proof had N
≥ 16(T+1)² where the coupling-failure bound needs (T+1)³; (3) the first
statement of the Hölder exponent claimed that polynomial decay with
p = 3 gives exponent 1/2 "which suffices" — correct, but the margin
(p > 13/5 against a measured 2.7 ± 0.3) was then honestly recorded as
marginal, and the measured modulus 0.75 added; (4) the initial reading
of the 9-curve-start data as "tail ten times smaller" was replaced,
after the 10^8 run, by "same tail after a transient".  FLAGGED AS
TAKEN FOR GRANTED (unchanged): τ < ∞ a.s.; the occupation-time step.
NEW flag: (M_p) in cor:noerg(ii) and (M_1) in prop:holder are
assumed, not proved (they follow from (R1) along the chain).

## 4. Pitfalls (all previous, plus)

54. NO state-uniform minorisation of the post-flip chain exists on
    Ω_K, for any n, ν, δ > 0 — inert dust (prop:nomino).  Any argument
    that needs P_x(Γ̄_{τ_n} ∈ ·) ≥ δν "for every x" is dead on arrival,
    and so is any that needs (V-)uniform ergodicity of κ₁.
55. (C2) is NOT geometric.  a_k(x) ~ (−1)^k C k^{−3} from low-B starts,
    universal C ≈ 2.5–3, reached after a start-dependent transient.
    Geometric decay is available only from the killing on Ω_K, with
    1 − ρ ~ e^{−cK}: fatal in the budget.
56. Averaging a_k over a burnt-in start CANCELS the tail, because the
    tail carries the sign (−1)^{B(x)} of the start (eq:akB).  Use
    fixed starts, or the sign-blind E_μ[a_k(x)²], to test (C2)-type
    statements.  Session 11's "ρ ≈ 0.3" was produced by this.
57. The block count B is a slow variable (E B_t ≈ log t) whose PARITY
    is fast at fixed times ((−1)^{B_t} deterministic) but slow at flip
    times.  Statements about B mod 2 sampled at random times are where
    the polynomial rates live.
58. A "for every environment" (adversarial) formulation of parity
    decorrelation is FALSE: an adversary choosing the merge gains at
    the (random) merge times can bias each inter-flip parity by a fixed
    amount.  Proofs must use the environment's true randomness.
59. The transport (ℓ¹) distance between mass configurations is NOT
    contracted by the natural coupling: when a dust discrepancy is
    finally hit it becomes a macroscopic one (amplification).  A
    dust-blind contraction, if one exists, must compare coarse laws.
60. A polynomial statement with exponent p gives Hölder exponent
    (p−1)/(p+1) by the crude bound and needs p > 13/5; do not quote
    "p = 3 is enough" without the margin.  The measured modulus (0.75)
    is better than the crude bound, but it is a measurement.

## 5. Papers on hand / to request

Unchanged (CGW_2008, KPS_2025, KSS, Schramm_2005, Kwan_Sudakov_2018 in
the tarball).  Meyn–Tweedie ch. 15–16 is NO LONGER the relevant
reference (cor:noerg).  If (γ) is chosen: Hairer–Mattingly–Scheutzow,
"Asymptotic coupling and a general form of Harris' theorem" (PTRF
2011) — the weak-Harris framework is the only known template for
chains with inert slow modes, but note pitfall 59.  If (α) shows
annealed suffices: nothing new needed — the exact laws are in
§sec:s10.

## 6. Quick reference (additions)

a_k = E_x[(−1)^{τ_k}] = (−1)^{B(x)} Σ_b (−1)^b P_x[B_{τ_k} = b];
κ₁^k(−1)^E = (−1)^E(−1)^k a_k;  t_k(w) = E[w^{τ_k}(−1)^{τ_k}] = κ_{−w}^k 1;
G(−w) = 1 + 2(−1)^E Σ_{k≥1}(−1)^k t_k(w) near z = −1.
Measured: a_k(one curve) ≈ (−1)^k 2.5 k^{−3} (k = 7..13), exponent
2.7 ± 0.3 (k = 8..20); a_1 = −.4004; from (0.3,0.3;.4): a_1 = −.004,
a_2 = −.0606; from 9 curves: a_1 = −.0967, |a_k| peaks .00145 at k = 7–8.
Modulus of Σ(w) at w → 1 (real): |1−w|^{0.75}.  P[B_{τ_k}=1] ≈ .5/k
from one curve.  Dust: hit probability ≤ 2ε/step; created at scale η
at rate ≈ 2η dη/step (arc-refined list).
Everything else: as HANDOFF_session_11 §6.

*The notes are the living document; keep PROVED/VERIFIED/COMPUTED/
MEASURED/OPEN labels honest.  Session-12 headline: the minorisation
route is impossible (inert dust), the geometric form of (C2) is false
(k^{−3}), the polynomial form is true numerically and sufficient, and
the one question that now dominates everything is whether the
programme needs uniformity over environments at all.*
