# HANDOFF: Cameron's second-row conjecture / CGW Conjecture — state after session 11

## >>> SESSION 11 HEADLINE: (R2) of session 10 is FALSE — I+K_z is
never invertible at z=±1 (K_1ε = −ε in the label-carrying bookkeeping;
κ_{−1}[(−1)^E] = −(−1)^E in the label-free one), so prop:condecay was
vacuous.  The repair is proved, not conjectured: the label is a
passenger (lem:autonomy, K_z = κ_z ⊕ (−κ_z)), and the IDLE LOOP —
split a curve, merge the two pieces back, exact probability
λ = (1/3)Σ_j r_j⁴ ≥ 1/(3(B+1)³) — gives UNCONDITIONALLY
   ‖κ_z‖_∞ ≤ (1−λ)/|1−λz²| ≤ 1/(1+λ(1−Re z²))   for all |z| ≤ 1,
hence uniform invertibility of I±κ_z and I−κ_z² on the whole closed
disc away from z=±1, with constant 6K³/(1−Re z²) on the block-bounded
space (thm:idleloop).  That is the non-lattice condition of Markov
renewal theory (pitfall 45, TODO 11b(v)) — PROVED.  z=−1 is not a
separate point: by κ_{−z} = −M̄κ_zM̄ it IS the renewal point z=+1 with
the test function (−1)^E.  What is left is two conditions at one
point on one Markov kernel κ_1 (the post-flip quotient chain):
   (C1) −1 ∉ spec κ_1 with a uniform inverse — PROVED to follow from
        one Doeblin minorisation shared by the n-th and (n+1)-st
        post-flip laws (prop:c1, ‖(I+κ_1)^{−1}‖ ≤ (2n+1)/(2δ));
   (C2) ‖κ_1^k (−1)^E‖_∞ ≤ Cρ^k — the environment block-count parity
        decorrelates geometrically along the post-flip chain.
(C2) = E_x[(−1)^{τ_k}] and is MEASURED: ρ ≈ 0.3, from a macroscopic
AND from a deeply trapped start (s11_kflip).  LEDGER: the core MOVED —
a theorem was proved on it and a stated condition was disproved. <<<

*Prepared 2026-09-15 (end of session 11).  Self-contained; companion
tarball has all code, data, papers, notes.  Previous handoffs valid
except: HANDOFF_session_10's (R2) is FALSE (see §1 below) — every
statement depending on it must be re-read with (C1)+(C2) substituted;
TODO 11b(iv) is disproved, 11b(v) is proved, 11d is done (with a
corrected prediction).  `HANDOFF.md` = this file.*

---

## 0. The problem in one paragraph

As before (HANDOFF_session_4 §0).  Weak conjecture d_TV(P_n,Q_n)=o(1)
STILL OPEN; owner directive: CGW Conjecture is the priority.  Chain:
thm:marked → programme steps 1–5 → open items now:
  (I)   prob:mixing, REDUCED to prob:survival (any α > 2/9, j ≤ 8d),
        which after §sec:s10 + §sec:s11 is REDUCED to (R1) + (C1) +
        (C2) on the log-macroscopic set S'_K, K = O(log d), with the
        environment handled by cor:horizon and the whole circle except
        z=±1 handled by thm:idleloop.
  (II)  transfer iid → orbit (unchanged).
  (III) prob:decouple (unchanged).

## 0.5 Notation used in this handoff (all of it, in one place)

One step of the chain = one fresh uniform chord = pick two uniform
points of the circle; same curve → SPLIT that curve at the two cuts,
different curves → MERGE them.  State and derived quantities:

| symbol | meaning |
|---|---|
| Γ | the marked state: curve masses + which curves carry u and v; when u ∼ v also the two arc masses of the marked curve |
| ε(Γ) | +1 if u ∼ v ("together"), −1 if not ("apart") |
| (p₁,p₂) | the marked pair: the arcs (x,y) when together, the two curve masses (a,b) when apart; s = p₁+p₂, gap g = min(p₁,p₂) |
| Γ̄ | the QUOTIENT state (p₁,p₂; env) — Γ with ε forgotten.  Autonomous (lem:autonomy) |
| B(Γ), E(Γ) | number of curves in total / number carrying neither mark.  B = E + 1 + 1{apart}; B changes by ±1 at every step |
| r(Γ) = (r_j) | ARC-REFINED mass list: the unmarked curve masses together with p₁ and p₂.  Σ r_j = 1, at most B+1 entries |
| flip | a step that changes ε.  Probability 2p₁p₂ from either sector |
| τ, τ_k | first flip time; k-th flip time |
| N_d | number of flips in d steps.  Φ_d(Γ) = E_Γ[(−1)^{N_d}] |
| K_z | first-flip kernel WITH the label: K_zF(x) = E_x[z^τ F(Γ_τ)] |
| κ_z | first-flip kernel on the QUOTIENT: κ_z g(x̄) = E_x̄[z^τ g(Γ̄_τ)].  K_z = κ_z ⊕ (−κ_z); κ_1 is the post-flip chain and is THE object of session 12 |
| λ(Γ) | idle-loop probability = (1/3) Σ_j r_j⁴ ≥ 1/(3(B+1)³) |
| H | the block-bound event {B_t < K ∀ t ≤ 9d}, P(H) ≥ 1 − (9d+1)^{−γ} (cor:horizon) |
| Ω_K | {B ≤ K}, the chain killed at the first split that would make B = K+1.  κ_1 there is SUB-Markov — prop:c1's proof is valid for sub-Markov κ_1 (its residual measures just have mass ≤ 1−δ) |
| S'_K | the log-macroscopic set {a,b ≥ 1/(16K)}, K = O(log d) — the starting set in prob:survival |
| π̄ | a stationary law of the post-flip quotient chain κ_1 on Ω_K (exists once the minorisation of prop:c1 does; does NOT exist on the unkilled space, where B → ∞) |
| M̄, N | multiplication by (−1)^E and by ε |

## 1. What session 11 added (all in `section_mixing_s11.tex`, §sec:s11)

### lem:autonomy — the label is a passenger (PROVED)
The quotient chain Γ̄ = (p₁,p₂; env) is autonomous: in BOTH sectors the
moves are flip (prob 2p₁p₂, (p₁,p₂) → (V, s−V), V ~ U(0,p₁)⋆U(0,p₂)),
shrink p_i → p_i√U (prob p_i²), grow p_i → p_i+q_j (prob 2p_iq_j), env
moves untouched.  (p₁,p₂) = arcs (x,y) when together, blocks (a,b) when
apart.  Consequences: (i) the sector-swap involution J commutes with
the kernel and fixes flip times; (ii) in F = g∘π + ε(f∘π),
K_z = κ_z ⊕ (−κ_z), so spec K_z = ±spec κ_z and
Σ_dΦ_dz^d = (1−z)^{−1}(I+κ_z)^{−1}(I−κ_z)1; (iii) P_π(tog) = 1/2 in
one line.  NOTE the structural payoff of (i): τ has the SAME law from
the two sectors, so the alternating renewal behind Φ_d is symmetric —
that is why Φ_d → 0 instead of → (E₊τ−E₋τ)/(E₊τ+E₋τ).

### lem:conj, prop:r2false — (R2) IS FALSE (PROVED)
N = mult by ε, M̄ = mult by (−1)^E.  Then N K_z N = −K_z (exactly one
flip at τ) and M̄ κ_z M̄ = −κ_{−z} (from lem:parity:
(−1)^τ = −(−1)^{E_τ−E_0}).  Hence spec K_z = −spec K_z and
spec κ_{−z} = −spec κ_z, and
  K_1 ε = −ε,      κ_{−1}[(−1)^E] = −(−1)^E.
So I+K_z has a bounded null vector at z=1 (labels kept) and I+κ_z has
one at z=−1 (labels dropped): prop:condecay is vacuous as stated on
EITHER bookkeeping.  The null vectors are bounded by 1, so no weight,
no smaller state space and no smoother test class removes them.
Session 10 missed it because (I+K_z)^{−1} is only ever applied to
(I−K_z)1, which degenerates at the same points — the real question is
whether the numerator's zero beats the resolvent's blow-up, and that is
not what (R2) says.

### prop:repair — the two singular points, exactly (PROVED)
(I+κ_{−w})^{−1}(I−κ_{−w})1 = M̄(I−κ_w)^{−1}(I+κ_w) M̄1.  So near z=1
the singular factor is (I+κ_z)^{−1} on a numerator that is O(|1−z|);
near z=−1 it is the ORDINARY renewal resolvent (I−κ_w)^{−1}, w→1, on
(I+κ_w)(−1)^E → 2(−1)^E.  The residue at z=−1 is ∝ ⟨π̄,(−1)^E⟩ — its
vanishing is exactly the decay of the alternating part of Φ_d found in
lem:parity(a) ((−1)^d c₁/d²).

### lem:idleloop — EXACT formula for the 2-step return (PROVED+VERIFIED)
With r(Γ) = the ARC-REFINED mass list (unmarked curve masses together
with p₁ and p₂; Σr_j = 1, ≤ B+1 entries),
  λ(Γ) := P(split a curve, then merge back its two pieces)
        = (1/3) Σ_j r_j⁴ ≥ 1/(3(B+1)³),
and up to null sets this is the whole event {Γ̄₂ = Γ̄₀}.  Per entry the
contribution is r⁴/3 = r²·2r²E[V(1−V)], V = |U₁−U₂|, E[V(1−V)] = 1/6;
in the together sector the arc-reassignment factor (x−w)/(s−w) cancels
the merge factor 2w(s−w) exactly.  It is the ARC-refined list because
splitting the marked curve ACROSS its arcs is a flip, not an idle loop.
VERIFIED to 3–4 digits in 3 configurations (`s11_ztau`, 2e6 trials).
Note E[Σp⁴(m)] = 1/4 + ∫(1−x³)(1−2x)^m dx (thm:sizebiased) ⇒ Eλ → 1/12;
the block bound is only needed to make the bound uniform.

### thm:idleloop — the contraction (PROVED+VERIFIED)
From the exact identity (1−λ(x)z²)(κ_zg)(x) = E_x[z^τ g(Γ̄_τ)1_{L^c}]:
  ‖κ_z‖_∞ ≤ sup_x (1−λ(x))/|1−λ(x)z²| ≤ 1/(1+λ_*(1−Re z²)), |z| ≤ 1.
On Ω_K = {B ≤ K} (chain killed at the first split that would make
B = K+1; the killing costs ≤ (9d+1)^{−γ} in Φ_d by cor:horizon):
  ‖κ_z‖_∞ ≤ 1 − min(t₀/(6K³), 1/(2K)),  t₀ = 1−Re z²,
so ‖(I±κ_z)^{−1}‖ ≤ max(6K³/t₀, 2K).  [Boundary layer B = K: EVERY
split kills and a split has probability Σ(mass)² ≥ 1/K, giving the
second branch — this is why no boundary bookkeeping is needed.]
VERIFIED: |E_x[z^τ]| vs the bound at θ/π = 0.1..0.9 from 3 starts
(`s11_ztau`): holds everywhere with margin (θ=π/2, one-curve start:
0.446 vs 0.920).

### prop:c1 — a minorisation suffices for (C1) (PROVED as a reduction)
If for some n, δ>0 and one ν, both P_x(Γ̄_{τ_n} ∈ ·) ≥ δν and
P_x(Γ̄_{τ_{n+1}} ∈ ·) ≥ δν for EVERY x, then I+κ_1 is invertible with
‖(I+κ_1)^{−1}‖ ≤ (2n+1)/(2δ).  (Proof: F = H_j + (−1)^jκ_1^jF with
‖H_j‖ ≤ j‖G‖; take j = n, n+1, multiply by (−1)^n, (−1)^{n+1} and
subtract: the δν terms cancel, 2δm ≤ (2n+1)‖G‖.)  The minorisation is
OPEN and is the classical difficulty: the fibre regenerates (two flips
dominate (1/4)Unif[s/4,3s/4]) but the base (s, env) does not.

### §sec:s11num — TODO 11d, run with a CORRECTED prediction
The handoff-10 prediction "E_x[(−1)^{τ_k}] decays geometrically under
(R2)" is meaningless ((R2) is false), but the quantity IS (C2):
κ_1^k(−1)^E = (−1)^E(−1)^kκ_{−1}^k1 and κ_{−1}^k1(x) = E_x[(−1)^{τ_k}].
Measured (`s11_kflip`, burn 4096, gap-accepted starts, no censoring):
* macroscopic (gap ∈ [1/8,3/8), 1.2e5 samples, se 0.0029):
  E[(−1)^{τ_k}] = −0.1306, +0.0263, −0.0100, +0.0011, then noise floor
  for all k ≥ 4.  Geometric, ratio ≈ 0.2–0.4, NO (−1)^k component.
* trap (gap ∈ [1e−3,1e−2), 4e4 samples, se 0.005): at the noise floor
  from k = 1 on.
* E[(−1)^{τ_k−τ_{k−1}}] → −0.170 ± 0.003 for k ≥ 3 from both starts
  (bounded away from ±1, as needed).
* post-flip chain equilibrates in ≈ 4 flips: from the trap start
  E[gap] at τ_k = .1995,.2375,.2468,.2484,.2491 → macroscopic .2498;
  E[s] = .715,.774,.790,.797,.798 → .799; sign correlations of the gap
  decay .317,.090,.022,.012 (ratio ≈ 0.3).  Evidence for prop:c1 with
  n ≈ 4.
* t²P(τ_k−τ_{k−1} > t) ∈ [15,20] on 16 ≤ t ≤ 128 for k = 2,4,8 from
  BOTH starts: (R1) holds uniformly ALONG the chain — TODO 11d(iii)
  answered yes.
* the assertion B(Γ_t) − B(Γ_0) ≡ t (mod 2) holds on all 1.6e5
  trajectories (hypothesis of lem:conj(ii)).

## 2. File manifest (additions)

| file | what | re-verify |
|---|---|---|
| `section_mixing_s11.tex` | notes §sec:s11 (input by section_mixing.tex after s10) | `pdflatex second_row_notes.tex` ×2 → 0 errors, 75 pp |
| `s11_kflip.c` | post-flip-chain diagnostics: E[(−1)^{τ_k}] (= (C2)), inter-flip parity, gap/s equilibration, dyadic tails; asserts the block-count parity identity | `gcc -O3 -o s11_kflip s11_kflip.c -lm && ./s11_kflip kflip 0.125 0.375 16 20000 2000000 7` → a_1 ≈ −0.13, a_k at noise for k ≥ 4, r_k ≈ −0.17, parity "HOLDS" |
| `s11_ztau.c` | measures λ(x) by exact 2-step state return and checks thm:idleloop's bound on \|E_x[z^τ]\| on a θ-grid | `gcc -O3 -o s11_ztau s11_ztau.c -lm && ./s11_ztau 0.5 1.0 200000 400000 3` → λ = 0.04153 vs exact 0.041667, all rows "ok" |
| `s11_kflip_mac.log`, `s11_kflip_trap.log`, `s11_ztau.log` | the runs quoted in §sec:s11num and the lemma verifications | — |

Notes edits: `section_mixing.tex` (\input s11); `section_mixing_s10.tex`
(a SESSION-11 CORRECTION paragraph before the status note of
prop:condecay); `second_row_notes.tex` (\spec macro).

## 3. PRIORITIZED TO-DO (session 12)

**QUICKSTART (≈4 min, 2 cores): as session 11 (`python3 pf_cl.py all`
→ 0 failures, SPLIT conv0 122 BY DESIGN; `python3 pf_master.py 5 8` →
0 failures; `gcc -O3 -o s8_splitmerge s8_splitmerge.c -lm &&
./s8_splitmerge psi 64 16 6000 50 992` → psi=0.034573; `gcc -O3 -o
s9_survival s9_survival.c -lm && ./s9_survival flux 8 200000` →
ratio ≈ 1.0–1.3 for b ≥ 10;
`python3 s10_blocks.py 60 4 3` → 2,4,6,8,10 / 3,9,19,33,51;
`python3 s10_blocksim.py 4 100000` → HZ vs MC agree to ±0.003) PLUS the two s11 re-verifications
in the manifest above, then `pdflatex second_row_notes.tex` twice
(0 errors, 75 pages).**

**THE EXACT TARGET (unchanged in substance, restated).**  Prove: there
are C, γ > 4/9 and K₀ such that for all K ≥ K₀, all x on the event H
with a,b ≥ 1/(16K), and all d ≥ 1,
    |Φ_d(x)| = |E_x[(−1)^{N_d}]| ≤ C K^c d^{−γ}        (S)
By the reduction paragraph at the end of §sec:markovrenewal plus the
occupation-time step (TODO 11c, still unwritten), (S) gives
prob:survival with α = γ/2 − o(1) > 2/9.  Sufficient for (S):
(R1) + (C1) + (C2).  Everything else on the closed disc is PROVED.

**OWNER'S CONCERN AND THE STOPPING RULE FOR SESSION 12 (read first).**
The session-10 ledger asked whether the residue is shrinking or being
renamed.  Session 11's honest entry: **the core moved.**  Three
distinct things happened that are not reformulations: (a) a stated
condition was DISPROVED by an explicit eigenvector; (b) a THEOREM was
proved about the cancellation object itself (thm:idleloop), with
explicit constants and no hypotheses beyond the block bound, settling
the whole closed disc except two points — this is the first proof step
on the object that sessions 8, 9 and 10 only renamed; (c) the two
points were shown to be one point (κ_{−z} = −M̄κ_zM̄).  What remains,
(C1)+(C2), is genuinely smaller and — this is the test — (C2) is a
SCALAR QUANTITY THAT CAN BE MEASURED, and it was measured (ρ ≈ 0.3).

Session 12 is (C1)+(C2) and is to be judged by exactly one of:

  (A) a PROOF of one of:
      (A1) the minorisation hypothesis of prop:c1 on Ω_K — some
           n = O(1), some ν, some δ = δ(K) with 1/δ(K) POLYNOMIAL in
           K.  This gives (C1) outright, with ‖(I+κ_1)^{−1}‖ ≤
           (2n+1)/(2δ), and (via uniform ergodicity) reduces (C2) to
           ⟨π̄,(−1)^E⟩ = 0;
      (A2) ⟨π̄,(−1)^E⟩ = 0 with a geometric rate, i.e. (C2) directly;
      (A3) a state-uniform lower bound on P(two consecutive inter-flip
           times have opposite parity) — this gives (C2) by TODO 12c.
      Any of the three, written in the notes with an honest label, is
      the session's deliverable; or

  (B) a clean NEGATIVE — and note carefully what does and does not
      qualify, because the budget is generous in K:

      HOW MUCH K-DEPENDENCE IS ALLOWED.  K = O(log d) and the target
      exponent is γ > 4/9 against a bound C K^c d^{−γ} with ANY fixed
      c.  So δ(K) = K^{−100} is FINE, ρ(K) = 1 − K^{−100} is FINE
      (what matters is that (1−ρ(K))^{−1} and 1/δ(K) are polynomial in
      K), and a constant that is exponential in K, e.g. δ = 3^{−K} or
      ρ = 1 − e^{−K}, is FATAL.  Do not report a large power of K as a
      failure.

      (B) therefore means one of: a state-uniform OBSTRUCTION to the
      minorisation (a family of states in Ω_K from which the τ_n and
      τ_{n+1} post-flip laws are mutually singular, or overlap only by
      e^{−cK}); or numerics in which the measured decay rate of
      E_x[(−1)^{τ_k}] degrades faster than polynomially in K over a
      family of starts in S'_K.  The cheap experiment for the latter:
      run `s11_kflip` with the burn-in length and the gap window
      varied so that the environment is progressively more fragmented
      (larger B, hence larger effective K), fit ρ(K), and check
      whether (1−ρ)^{−1} grows polynomially.  ~20 min.  If it does not,
      say so and STOP working this route.

A FIFTH name for the cancellation core does not count.  In particular:
do NOT restate (C2) as "aperiodicity", "non-lattice", "spectral gap"
or "V-geometric ergodicity" and call that progress — those are the
names; the content is the minorisation, the parity average, or the
odd-perturbation bound.

**TIME BUDGET.**  12a first (≈ 2 h).  Then 12b.  At the half-way point
of the session, stop and ask which of (A1)/(A2)/(A3)/(B) is actually
within reach; if none is, write down precisely which uniform lower
bound fails, record it as the session's finding, and hand off — do NOT
start 11a/11c/11e/10c/10d to have something to show.  Those are
carried tasks with known routes and they will still be there.

**START HERE.**

**Step 0 (ten minutes, before anything else): TRY TO DISPROVE (C1) AND
(C2).**  This is the lesson of session 11 turned into a procedure (see
pitfall 53).  (R2) was carried as the centrepiece of a session's brief
and is false by a one-line argument; the project's labels record
whether a claim is PROVED, not whether an OPEN target is ACHIEVABLE.
The disproof that worked was: find a BOUNDED eigenfunction built from a
state function that the dynamics toggles DETERMINISTICALLY.  The three
such functions here are ε (toggles at every flip), (−1)^B (toggles at
every step), and (−1)^E (toggles at every NON-flip step; flips leave E
fixed, since together→apart splits the marked curve into two marked
curves and apart→together merges them).

*Session 11 already ran this check on (C1)/(C2) and it comes out
clean — do not redo it, but do read it, because the method is what you
should apply to any new condition you adopt.*  On the quotient the only
candidate is (−1)^E, and
    κ_1[(−1)^E](x̄) = (−1)^{E(x̄)} · E_x̄[(−1)^{E_τ−E_0}]
                    = −(−1)^{E(x̄)} · E_x̄[(−1)^τ],
using (−1)^{E_τ−E_0} = −(−1)^τ (lem:conj(ii)).  So (−1)^E is an
eigenfunction of κ_1 **iff E_x[(−1)^τ] is constant in x** — and it is
measurably NOT: −0.1306 ± 0.003 from the macroscopic start against
−0.0009 ± 0.005 from the trap start (s11_kflip logs), about 45 standard
errors apart.  So the (R2)-style obstruction does not apply to (C1) or
(C2), and neither is false for the cheap reason.  That is a real, if
small, positive check: it is the first time one of these named
conditions has survived a deliberate attempt to kill it.

Then read §0.5 above (notation), then §sec:s11 of the notes (~7 pages),
then only the following from earlier sections: lem:parity (s9),
cor:horizon and §sec:markovrenewal (s10).

**DO NOT, in session 12:**
  1. try to prove (R2) — it is false (pitfall 47);
  2. re-derive no-dust / block-count estimates — cor:horizon is a black
     box and is used exactly in the form it was proved;
  3. re-prove anything on the circle away from z = ±1 — thm:idleloop
     is unconditional there;
  4. treat z = −1 as a separate point — prop:repair makes it z = +1
     with test function (−1)^E;
  5. look for another EVEN uniform event (pitfall 49) — the idle loop
     is already the best of those, and what is missing is odd;
  6. run long numerics (2 cores) or open 11a/11c/11e/10c/10d before
     (A) or (B) is settled;
  7. report a large POWER of K as a failure (see the budget above).

### TODO 12a — the Lyapunov drift for λ (do this FIRST, ≈ 2 h)
SUCCESS = an inequality κ₁V ≤ ϱV + b with ϱ < 1, b < ∞, for
V = (Σ_j r_j⁴)^{−θ} = (3λ)^{−θ} and SOME θ > 0 (θ may be as small as
you like; it does not enter the final exponent).  Pay-off is double:
it replaces the killing at B = K by a weighted space — removing the
only place where thm:idleloop is not stated on the full state space,
and the only place where κ₁ is sub-Markov — and it supplies the "base"
half of the minorisation of prop:c1.  A finite computation of the same
kind as lem:logdrift.  Ingredients in hand: the exact law of Σp^r at
every time (thm:sizebiased, r = 4 included); Σ r_j⁴ ≥ (Σ r_j²)²/(B+1);
the flip map's trapezoid; the move list in §0.5.
FALLBACK if the drift will not close for any θ: state the obstruction
in one paragraph and continue with the killed chain on Ω_K, which is
enough for everything except the cosmetic point above.  Do not spend
more than 2 h here — 12a is an enabler, not the deliverable.

### TODO 12b — the minorisation of prop:c1 (THIS IS THE SESSION)
SUCCESS = integers n = O(1), a probability measure ν on Ω_K, and
δ = δ(K) with 1/δ polynomial in K, such that FOR EVERY x ∈ Ω_K
    P_x(Γ̄_{τ_n} ∈ ·) ≥ δν(·)   AND   P_x(Γ̄_{τ_{n+1}} ∈ ·) ≥ δν(·).
Three points that are easy to get wrong:
  (i) the two bounds must share ONE ν and ONE δ.  That is the whole
      trick of prop:c1 — subtracting the two identities cancels the
      δν terms and leaves a factor 2δ.
  (ii) "for every x" means every x ∈ Ω_K, traps included.  A bound
      that degrades as g = min(p₁,p₂) → 0 is not enough; but note that
      a trap only makes τ₁ LONG, and the numerics show a long first
      inter-flip time helps (it randomises everything, see the trap
      row in §sec:s11num).  Prove that, do not assume it.
  (iii) the idle loop is NOT the way to convert an n-path into an
      (n+1)-path: it returns Γ̄ exactly but adds 0 flips.  An honest
      extra flip is needed.  (Pitfall 49.)
Structure of the proof, with what is already in hand:
  * fibre: DONE — two consecutive flips give a p₁'' dominating
    (1/4)Unif[s/4,3s/4] (§sec:markovrenewal).
  * s-coordinate: force s into a fixed window using lem:srecharge and
    lem:logdrift (+ 12a).
  * environment: the remaining work.  Drive the env into a
    neighbourhood of a reference configuration (e.g. one curve of mass
    ≈ 1−s, or two of mass ≈ (1−s)/2) by a fixed finite sequence of
    merges, each of probability ≥ (mass product) ≥ K^{−2}; with < K
    curves on Ω_K, ≤ K merges suffice, giving δ ≥ K^{−2K} — which is
    NOT polynomial.  Getting from K^{−2K} to K^{−O(1)} is the real
    content: do it by NOT insisting on a fixed reference configuration
    but on a set of configurations of ν-measure bounded below, or by
    using the O(6) steps between flips (≈ 6n steps for n flips) more
    cleverly than one merge at a time.

### TODO 12c — (C2) directly, if 12b stalls
Show sup_x |E_x[(−1)^{τ_k}]| ≤ Cρ^k.  Needs a state-uniform lower
bound on the probability that two consecutive inter-flip times have
OPPOSITE parity.  The idle loop gives a uniform lag-2 (even)
perturbation and is therefore useless here: this is the one place in
the whole problem where an odd/even asymmetry must be produced rather
than conjugated away.  Candidate odd perturbation: a single
marked-curve shrink inserted before the flip (changes τ by 1, occurs
with probability p_i² ≥ 1/(2(B+1)²)) — but it also changes the state,
so it must be composed with a return; find a 3-step return, or couple
the two resulting states.  SUCCESS = a constant c(K) with 1/c
polynomial such that, from every x ∈ Ω_K,
    P_x[(−1)^{τ_{k+1}−τ_k} = −(−1)^{τ_k−τ_{k−1}}] ≥ c(K)
uniformly in k; this gives (C2) with 1−ρ ≍ c by a standard
two-block argument.  Measured value of the relevant quantity:
E[(−1)^{τ_k−τ_{k−1}}] ≈ −0.17, i.e. c is of order 0.4 in practice, and
it does not drift with k.

SCOPE A NEGATIVE HERE CORRECTLY.  It is quite possible that no
state-uniform odd short return exists at all — every free uniform event
found in eleven sessions is EVEN, and that may reflect something
structural about a chain whose block count has deterministic parity
rather than a failure of ingenuity.  A PROOF of that would be a
genuinely valuable finding and is worth writing up.  But it kills
ROUTE 12c ONLY, not (C2): (C2) can still come from 12b/(A1) via uniform
ergodicity plus ⟨π̄,(−1)^E⟩ = 0, which needs no odd event anywhere.  Do
NOT report an odd-event obstruction as outcome (B) for the session
unless 12b is independently blocked as well; record it as a structural
finding and state explicitly which of the three routes to (C2) it
removes and which it leaves open.  **If neither 12b nor 12c is in sight by
mid-session, that is outcome (B): write down precisely which uniform
lower bound fails and hand off.**

### TODO 11c (carried) — the occupation-time lemma for s
From lem:logdrift: P(#{t ≤ T: s_t < 1/(4K)} ≥ T/2 | s_0 ≥ 1/(2K), H)
≤ CK^c/T^? ; polynomial is enough.  Standard (Lyapunov s^{−θ}); still
NOT written.  Needed by the reduction paragraph and by (R1).

### TODO 11a (carried) — (R1)
sup_{x∈S'_K} P_x(τ > t) ≤ C(K)/t².  Ingredients all proved; assemble
as a subcritical renewal with kernel ≈ 2/(s²m²) (§sec:parity(c)).
Session 11 added: the 1/t² tail now VERIFIED to hold uniformly along
the post-flip chain (k = 2,4,8, both starts), so the statement being
assembled is the right one.  Half a day, lower risk than 12b.

### TODO 11e (carried) — (u,v)-averaged ensemble / hooks for ψ̄
Unchanged; do not open before (A) or (B).

### TODO 10c/10d (transfer, decouple): unchanged.

### TODO 12d (referee): §sec:s11 was refereed in-session.  Caught and
fixed: the first version of lem:idleloop claimed λ ≥ q⁴/6 from the
largest curve and was contradicted by the measurement (0.0415 vs the
claimed 0.1667) — the error was forgetting that splitting the marked
curve across its two arcs is a FLIP; the corrected statement is the
exact arc-refined formula, which matches to 4 digits.  Also caught: a
first attempt to state the repair as a condition on (I−K_z²)^{−1} with
the label kept, which is right but clumsy — the label-free
κ_z is the correct object (lem:autonomy).  FLAGGED AS TAKEN FOR
GRANTED (unchanged from s10): τ < ∞ a.s. on {ab>0}; the
occupation-time step.  NEW flag: thm:idleloop is stated on the KILLED
chain on Ω_K; the killed Φ_d differs from the true one by cor:horizon,
but the resulting sub-Markov κ_1 is what (C1)/(C2) must be proved for.
TODO 12a removes this.

## 4. Pitfalls (all previous, plus)

47. (R2) of prop:condecay is FALSE.  Do not try to prove it, and do
    not quote prop:condecay without substituting (C1)+(C2).
48. Pitfall 45 ("(R2) is needed on the whole unit circle") is now
    SETTLED in the affirmative by thm:idleloop — but only on a state
    space where λ is bounded below, i.e. only with the block bound.
    On the full state space inf λ = 0 and the theorem says nothing.
49. The idle loop changes the number of steps by 2 and the flip index
    by 0.  It is therefore useless for any statement that needs an ODD
    perturbation (TODO 12c).  Every "free" uniform event found so far
    is even; that is the precise residue of the parity problem.
50. λ uses the ARC-REFINED list: for the marked curve in the together
    sector it is (x⁴+y⁴)/3, not s⁴/3.  The error is a factor of up to
    8 and was caught only by the exact numerical check.
51. E_x[(−1)^{τ_k}] is NOT predicted to decay by (R2) (which is
    false); it is predicted to decay by (C2).  Keep the two
    bookkeepings straight — they give DIFFERENT eigenvalues:
      label-carrying:  K_{−1}[(−1)^B] = +(−1)^B   (eigenvalue +1),
      label-free:      κ_{−1}[(−1)^E] = −(−1)^E   (eigenvalue −1).
    Only the second is an obstruction to (R2); the first is harmless,
    and it is why E_x[(−1)^{τ_k}] has no a-priori reason to decay in
    the label-carrying picture.  Read handoff-10's TODO 11d(i) with
    this correction.  (An earlier draft of this handoff conflated the
    two; caught in the session-11 referee pass.)
52. spec K_z = ±spec κ_z: any spectral statement about the
    label-carrying kernel is automatically symmetric under λ ↦ −λ, so
    "−1 ∈ spec" is never informative there.
53. THE LABELS TRACK PROOF, NOT ACHIEVABILITY.  PROVED / VERIFIED /
    COMPUTED / OPEN say whether a claim has been established; nothing
    in the system says whether an OPEN target is TRUE.  (R2) was
    OPEN — correctly labelled — and false by one line, and it survived
    a full session plus an in-session referee pass.  Before adopting
    any new named condition as a session's target, spend ten minutes
    trying to disprove it; the productive move is to look for bounded
    eigenfunctions built from state functions the dynamics toggles
    deterministically (here ε, (−1)^B, (−1)^E).  See START HERE step 0
    for the worked check on (C1)/(C2).

## 5. Papers on hand / to request

IN THE TARBALL: CGW_2008.pdf, KPS_2025.pdf, KSS_large_deviations.pdf,
Schramm_2005.pdf, Kwan_Sudakov_2018.pdf + notes PDF.  Keep packing.
To request if needed: Meyn–Tweedie ch. 15–16 (drift/minorisation — now
the single most relevant reference, for TODO 12a/12b); Alsmeyer, "The
Markov renewal theorem and related results"; Harer–Zagier 1986 /
Zagier 1995; Diaconis–Freedman 1980 (transfer).

## 6. Quick reference

Master formula: unchanged.  Factorization: 4Ebias² ≤ max_{j≤k, d≥k/8}
ψ(j,d) + 2e^{−k/32}; ψ(j,d) = E_{μ_j}[Φ_d²]; Φ_d = E[(−1)^{N_d}].
Short horizon: P_{μ_j}(g ≤ ε) ≤ 4ε²j.  Budget: α > 2/9.
Block count: E[3^{B_t}] = 3+4t+2t², P(max_{t≤T}B_t ≥ K) ≤ 3^{1−K}(T+1)³,
B_t ≡ t+1 (mod 2).  Size-biased block at time m: density 1−(1−2x)^m,
atom (1+(−1)^m)/(2(m+1)); E[Σp^r(m)] = 1/r + ∫(1−x^{r−1})(1−2x)^m dx.
Quotient chain moves (both sectors): flip 2p₁p₂ → (V,s−V),
V ~ U(0,p₁)⋆U(0,p₂); shrink p_i → p_i√U at p_i²; grow p_i → p_i+q_j at
2p_iq_j.  Generating function: Σ_dΦ_dz^d = (1−z)^{−1}(I+κ_z)^{−1}
(I−κ_z)1.  Conjugations: N K_z N = −K_z, M̄ κ_z M̄ = −κ_{−z}.
Idle loop: λ = (1/3)Σr_j⁴ ≥ 1/(3(B+1)³); ‖κ_z‖ ≤ 1/(1+λ(1−Re z²)).
Label-carrying form of the conjugation: M K_z M = K_{−z} with
M = (−1)^B, so K_{−1}[(−1)^B] = +(−1)^B (harmless); (−1)^B = −ε(−1)^E.
π: P(tog)=1/2 (now one line from lem:autonomy); π(g≤δ) = 4δ+O(δ²);
ψ_π(d) ≥ c/d (PROVED), ≈1.0/d.  Annealed one-point check: from the
one-curve start E[Φ_m] = 1/(m+1) (m even), 1/(m+2) (m odd).

*The notes are the living document; keep PROVED/VERIFIED/COMPUTED
labels honest.  Session-11 headline: the stated cancellation condition
was false, the idle loop proves its repair everywhere except at
z = ±1, those two points are one point, and the residue is now a
measurable scalar — measured, and decaying at ρ ≈ 0.3.*

## 7. PROJECT-LEVEL NOTE — FOR THE OWNER, **NOT** SESSION-12 WORK

Session 12 should read this and then ignore it; it is about how
sessions are allocated, not about (C1)/(C2).  Four observations from
stepping back at the end of session 11.

**(a) The foundation has not been audited since it was laid.**  The
(R2) episode is evidence about more than (R2): a target that is false
by a one-line argument survived a session and a referee pass.  The
in-session referee checks the session's OWN new material; nothing
checks the accumulated chain.  Recommend one session, at a time of the
owner's choosing, spent entirely on an adversarial audit of
thm:marked → programme steps 1–5 → prob:survival, with priority on the
items explicitly TAKEN FOR GRANTED rather than proved: τ < ∞ a.s. on
{ab>0}; the occupation-time step; the budget arithmetic of
cor:kappafree; the short-horizon reduction in thm:factor.  Deliverable
per item: a proof, a counterexample, or an explicit statement of what
would make it false.  This is cheap insurance and it gets more
valuable, not less, the closer (I) comes to closing.

**(b) Two of the three open branches have never been examined.**
Seven sessions (5–11) have gone into item (I).  Item (II) (transfer
iid → orbit) and item (III) (prob:decouple) have been carried forward
untouched since at least session 10.  Their difficulty is therefore
UNKNOWN — not low, unknown.  Recommend a half-session of scouting each,
purely to classify them (routine / hard / as hard as (I)), and to do it
BEFORE (I) closes: discovering a second hard problem late is the
expensive order to discover it in.

**(c) The carried "standard" items are accumulating.**  TODO 11a (the
(R1) assembly) and TODO 11c (the occupation-time lemma) are each
labelled "half a day, standard" and each has now been carried two
sessions, because the headline task always eats the session.  Unwritten
standard steps are exactly where errors hide — see (a).  Consider
ring-fencing one session for the carried items alone, with the headline
task explicitly forbidden.

**(d) Honest trend.**  The residue is now FALSIFIABLE for the first
time: (R2) could not be simulated, (C2) is a scalar and was simulated
in twenty minutes.  That genuinely bounds the downside — if this route
is wrong, session 12 or 13 finds out cheaply instead of renaming the
obstacle a fifth time.  Against that: the remaining piece is not
obviously easier than what preceded it (K^{−2K} → K^{−O(1)} is real
content, and it is the same "the base is not regenerated" difficulty
that sessions 8–10 deferred), and the even/odd gap may be structural.
The fair summary is that the project has converted a vague obstacle
into a sharp one on a single branch, while two required branches remain
unexamined and the foundation has not been re-checked.  **The dominant
uncertainty in this project is no longer the piece being actively
worked on.**
