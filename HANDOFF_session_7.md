# HANDOFF: Cameron's second-row conjecture / CGW Conjecture — state after session 7

## >>> SESSION 7 HEADLINE: the TWO-POINT FACTORIZATION (Theorem thm:factor,
notes §sec:anneal; PROVED, three-line conditioning argument) reduces the
iid model of the annealed bias bound to ONE mixing statement:
ψ̄(d) = sup_j E_C[M_d(C)²] ≤ d^{−2c} for some c > 1/11
(Problem prob:mixing), where M_d = P^d ε is d steps of a fresh-chord
averaging Markov operator. Measured: ψ ≈ 0.8/d (c ≈ 0.52 — factor-5
margin persists). TWO REFUTATIONS: (a) pendant stars kill EVERY local
crossing-multiplicity genericity condition (prop:pendant, PROVED); (b)
template averaging is EXACTLY NEUTRAL (prop:vacuity, PROVED) — the
session-6 proposed fix of rem:bookkeeping(i) is withdrawn; the true
remaining bookkeeping item is the supply–mark decoupling
(Problem prob:decouple, NEW, open for classes < log²⁴n columns). REAL
ENSEMBLE (7c): master formula verified against ground truth at
n=30,50,100 (0 mismatches in 3250 subset checks on JM-sampled squares);
real-orbit E|bias| is ~1.5–1.6× the iid value at equal k, with a
short-arc/r=0 atom of measure ≈ 0.85/k (6.5%/3.9%/2.9% at n=30/50/100,
bias exactly +1/2) which accounts for the ENTIRE signed mean;
E[bias²]·k ≈ 0.7 across n — the annealed 1/k second-moment law
survives on real squares. <<<

*Prepared 2026-08-28 (end of session 7). Self-contained; the companion
tarball has all code, data, papers, notes (`second_row_notes.tex`/`.pdf`,
now with §sec:anneal = `section_annealed.tex`). Previous handoffs remain
valid except where superseded (notably: `HANDOFF_session_6.md` §1's
rem:bookkeeping(i) "TEMPLATE AVERAGING" proposal is REFUTED, and its
TODO 7b is resolved-by-refutation into Problem prob:decouple). In this
tarball `HANDOFF.md` = this file; the session-6 handoff is
`HANDOFF_session_6.md`.*

---

## 0. The problem in one paragraph

As before (HANDOFF_session_4 §0). Weak conjecture d_TV(P_n,Q_n)=o(1)
STILL OPEN; owner directive: CGW Conjecture is the priority. Chain
unchanged: thm:marked (EQ(δ) ⇔ marked bit fair to δ/4 on X) →
programme steps 1–5 of §sec:orbit summary → open items now:
  (I)   prob:mixing (iid model of the annealed bound; via thm:factor)
  (II)  transfer iid → true orbit measure (7a(i), untouched this session)
  (III) prob:decouple (supply–mark decoupling for small classes)

## 1. What session 7 added (all in notes §sec:anneal = `section_annealed.tex`)

### Theorem thm:factor + Corollary cor:psienough (PROVED)

iid ensemble (u,v fixed, chords iid): ε_S = ±1 the u~v sign. For S,T
independent uniform subsets of [k], conditioning on the shared chords C:
    E[ε_S ε_T] = E_C[ M_{d1}(C) M_{d2}(C) ],  M_d(C) = E[ε(C ∪ F_d)|C],
so 4E[bias²] = E_multinomial E_C[M M] ≤ max_{j≤k, d≥k/8} ψ(j,d) + 2e^{−k/32},
ψ(j,d) := E_{C_j}[M_d(C_j)²]. Any ψ̄(d) ≤ Kd^{−2c}, c > 1/11 closes the
iid model (E|bias| ≤ √(Ebias²)). KEY: M_d = P^d ε (P = add one fresh
active chord, average) — an L² mixing rate for a coagulation–
fragmentation-type chain (Schramm [Sch05]; one-point fn is Harer–Zagier
[HZ86] material — bibitem added).

### Lemma lem:freeze (PROVED): ψ(j,d+1) ≤ ψ(j+1,d) (conditional Jensen);
hence ψ̄(d) = sup_j ψ(j,d) is non-increasing. Open problem stated as
ψ̄(d) ≤ Kd^{−2c} (prob:mixing).

### Numerics (`s7_twopoint.c`, cycle route — INDEPENDENT of the rank route)
- 4E[bias²]·k = 0.95, 1.08, 1.22, 1.28, 1.27, 1.44 at k=8..256
  (Ebias² ≈ 0.33/k; possible slow log growth — watch at larger k).
- ψ(j,d): j=16: .124/.019/.0017/.0004 at d=4/16/64/256; j=64:
  .174/.041/.0084/.0010; j=256: .197/.050/.0114/.0027; j=1024:
  .064 (d=16), .0128 (d=64). Fit: ψ(256,·) ~ d^{−1.05}; j-dependence
  saturating. c ≈ 0.52.
- one-point E[ε(F_m)]: +.024/+.005/+.001 at m=4/8/16; ≤1e−3 beyond
  (statistically 0 at 4M–8M samples).

### Proposition prop:pendant (PROVED + verified `s7_bias_lab.py pendant`)

r nested separating chords, each crossed by m private nested pendants:
bias = (1−2^{−(m+1)})^r − 1/2 → −1/2 for every fixed m. So "every
separating chord crossed by ≥ m others" NEVER suffices (need m ≳ log r).
Exhaustive scan k≤6 (all matchings × all v-positions; CSVs
`s7_exhaust_k{4,5,6}.csv`): |bias| ≡ 0 only on ladders and on
"r=1 + separating chord crossing NOTHING" (cut-chord identity gives
exact 0 there); all other r≥1 classes reach ≥ 1/16 (all-separating
r=k dense classes) resp. → 1/4 (everything else). ALSO: per-coordinate
toggle route dead: E[min_i P(no toggle)] ≈ .35–.40, NON-decreasing in
k≤12 (`s7_bias_lab.py toggle`).

### Proposition prop:vacuity (PROVED): template averaging is neutral

For every column permutation h: E_{(L,P)~X}[bad(h(T);L,P)] =
E_{(L,P)~X}[bad(T;L,P)]. Proof: column perms fix σ as symbol permutation
(pairs (L(1,c),L(2,c)) move with columns) ⇒ S(λ) and X closed; marks and
the intrinsic exceptional set are covariant; substitute (L,P) =
(h L₀, h P₀). Same for random templates with column-exchangeable law
(covers the KPS iid-Bernoulli template of their Lem 9.2 → 8.3).
CONSEQUENCE: supply bookkeeping is PROVED only for classes with
λ_m·m ≥ log²⁴n (deterministic |E| ≤ 6log²²n); smaller classes = NEW
Problem prob:decouple. Suggested framing (recorded in notes): condition
on bottom (n−2)×n rectangle; class columns = leftover-cycle orientations
(±1 per leftover cycle, tilted by S(λ)); exceptional set is
L↓-measurable UP TO stability's weak dependence on intercalates of T∩L
meeting rows 1,2. NOTE ALSO (recorded): per-column supply CANNOT come
from the expander property (bottoms out at ℓ) nor from lem:condexp-style
division (KPS per-tuple rates too weak vs e^{O(n log n)}).

### TODO 7c DONE: real-ensemble grounding (`s7_real_bias.py` + `jm.c`)

JM-uniform squares, uniform marked instances {j,σ^α(j)}, greedy maximal
cell-disjoint intercalate families through the marked columns (rows≥3),
orbit data exactly as pf_master (same conventions). GROUND TRUTH: master
formula re-verified by actually switching the sampled squares: 0/3250
mismatches total (0/750 n=30 pilot, 0/1250 n=30, 0/750 n=50, 0/500
n=100). Bias data (kcap 26; k>14 by MC over 1500–2000 subsets):

| n | inst | k̄ | E\|b\| | E\|b\| (r>0) | P(r=0) | signed | Eb²·k̄ | iid E\|b\| |
|---|---|---|---|---|---|---|---|---|
| 30 | 4200 | 12.9 | .168 | .145 | .065 | +.038 | 0.67 | .098 |
| 50 | 3300 | 22.6 | .123 | .108 | .039 | +.023 | 0.75 | .074 |
| 100 | 1800 | 26 (cap) | .114 | .101 | .029 | +.013 | 0.73 | .069 |

KEY: (a) Eb²·k̄ ≈ 0.7 across n — the 1/k law holds on real squares
(constant ≈ 2.2× iid; E|bias| ≈ 1.5–1.6× iid at same k — position law
of intercalate rows in π_P cycles is not iid, as expected); (b) at
FIXED k=26 the n=50 and n=100 values agree (.107/.112) — elevation is
a k-level property, not an n-effect; (c) the r=0 atom is ≈ 0.85/k̄ and
its contribution (atom × 1/2) equals the signed mean within error —
the short-arc self-regularization of §sec:annealed measured exactly;
(d) cut-chord frequency 2.1%/1.6%/1.1%, no new pathology class.
Per-instance CSVs: `s7_inst{30,50,100}.csv` with
(n, k, bias, nsep, ncut, ncyc) for regression/pathology hunting.

### Referee pass (TODO 1″ items for §sec:orbit): PASSED with notes
- lem:split convention matches pf_cl conv1 (earlier-applied → copy 0) ✓
- lem:supply arithmetic OK for β ≤ 1/2 (application has β=1/10); for
  general β the hypothesis should read 2tℓ+2 ≤ (1−β)n etc. — cosmetic,
  not yet edited into notes.
- lem:condexp: D_λ, 3ⁿ-path, e^{O(√n)} type-count all check out.
- lem:border/thm:bitrank/lem:passenger: proofs correct.

## 2. File manifest (additions; all else as before)

| file | what | re-verify |
|---|---|---|
| `s7_bias_lab.py` | pendant stars, exhaustive k≤6 scan (CSVs), toggle stats | `python3 s7_bias_lab.py pendant` (~1 min; expect all OK); `exhaust 5` seconds, `exhaust1 6` ~5 min; `toggle` ~10 min |
| `s7_twopoint.c` | cycle-route: 4Ebias², one-point, ψ(j,d) | `gcc -O3 -o s7_twopoint s7_twopoint.c -lm && ./s7_twopoint bias 64 400000` |
| `s7_real_bias.py` | real-ensemble bias via JM + master formula | `gcc -O3 -o jm jm.c && ./jm 30 40 30000 42 \| python3 s7_real_bias.py 30 --verify 30` (expect 0 mismatches) |
| `section_annealed.tex` | notes §sec:anneal (input by main tex) | pdflatex second_row_notes.tex (clean) |
| `s7_*.log`, `s7_exhaust_k*.csv`, `s7_inst*.csv` | session-7 runs/data | — |

Notes edits: section_orbit.tex (rem:bookkeeping(i) withdrawn→pointer,
genericity paragraph updated, summary item 2 split PROVED/OPEN),
second_row_notes.tex (input section_annealed, HZ86 bibitem, manifest
appendix extended, first open problem updated).

## 3. PRIORITIZED TO-DO (session 8)

**QUICKSTART (before building anything): re-verify the stack —
`python3 pf_cl.py all` (0 failures; SPLIT conv0's 122 are BY DESIGN),
`python3 pf_master.py 5 8` (0 failures), and one s7_twopoint sanity
(`./s7_twopoint bias 32 200000` → 4Ebias2 ≈ 0.038). Suggested order of
attack: 8b(γ) + 8d first (hours of background compute — launch early),
then 8a as the session's mathematical core, 8c and 1‴ as secondary.**

### TODO 8a (mathematical core): attack prob:mixing — ψ̄(d) ≤ d^{−2c}.
Concrete routes, roughly in order of promise:
  (i) Direct probabilistic estimate of M_d(C): after d fresh chords,
  P(u~v) − 1/2 given C. The fresh chords perform Schramm-type
  coagulation–fragmentation on cycle space; u,v connectivity after d
  steps should decorrelate at a polynomial rate. Try small-C cases
  first: C = single frozen chord; compute M_d exactly (transfer-matrix
  over the finite "arc-occupancy" state?) — the state space for ONE
  frozen chord + virtual chord is finite-dimensional-ish.
  (ii) Character/map-enumeration: E[ε_S ε_T] is a two-matrix-model
  cycle statistic; Harer–Zagier gives one-point exactly — look for the
  two-point analogue (request Arratia–Bollobás–Sorkin interlace
  polynomial papers + HZ 1986 if going this route).
  (iii) Weaker target that still closes: ANY polynomial decay
  (c > 1/11 = margin 5.5×). E.g. couple d fresh chords so that with
  probability ≥ 1 − d^{−c'} a "randomizing gadget" (e.g. a fresh ladder
  pair straddling u,v, crossed by nothing frozen... careful: pendant
  moral says gadgets must be measure-based) — actually prefer (i).

### TODO 8b: transfer 7a(i) — iid → orbit positions. The factorization
proof only needs conditional iid of fresh chords given shared ones.
Under the orbit measure the chord positions are exchangeable (rows
3..n) but not independent. Options: (α) exchangeable-pair/de Finetti
approximation at k ≪ n; (β) prove the factorization inequality with a
correlation-error term and bound it by intercalate-density
concentration (KSS); (γ) empirical: test the CONCLUSION under the real
ensemble — E[bias²]·k vs k from the s7_inst*.csv data (and an n=200
run). CAUTION (learned this session): the factorization PROOF does not
transfer per-instance — for a fixed family, S∖T and T∖S are disjoint
subsets of the same coordinate pool, not independent; independence came
from position-randomness of distinct chords. So under the orbit
measure the factorization must be run with the SQUARE's randomness
(conditional on the shared intercalates' rows, resample the rest —
an MCMC-with-constraints question), or replaced by a correlation
bound. NOTE: (γ) is cheap — do it first.

### TODO 8c: prob:decouple. Start with the L↓-conditional framing
(recorded in notes §sec:vacuity): quantify (a) how far the exceptional
set is from L↓-measurable (intercalates of T∩L meeting rows 1,2 are
sparse: density ε⁴-ish), (b) orientation-randomness of class columns
given L↓ under S(λ)-tilt. Also: empirical measurement — in the 7c
pipeline, tag each instance with whether its marked columns are
supply-poor (k small) and with its class size λ_m·m: measure the
correlation on real squares.

### TODO 8d: numerics refinements. (a) n=200 run (kcap 26) to see the
n-trend continue; (b) 4Ebias²·k drifted 0.95→1.44 over k=8..256 in the
iid ensemble — run k=512, 1024 (s7_twopoint bias) to decide constant vs
log k (a log factor is harmless for the budget but should be known);
(c) pin the one-point exponent f(d) (currently: decays at least ~d^{-2},
below MC resolution by d=32 — could be computable EXACTLY via
Harer–Zagier two-point/face-connectivity: good warm-up for 8a(ii)).

### TODO 1‴ (referee): §sec:anneal is new — referee it (proofs are
short); also fold the β ≤ 1/2 caveat into lem:supply's statement.

### TODO 2′, 4: unchanged (real merit, lower priority).

## 4. Pitfalls (all previous, plus)

26. Template averaging is a NO-OP for the X-averaged bad-mark fraction
    (prop:vacuity): everything is jointly covariant under column
    relabelling; do not resurrect it. The remaining freedom that is NOT
    killed by covariance: statements conditional on L↓ (bottom
    rectangle), since column relabelling does not fix L↓-fibers.
27. Local genericity conditions on the interlacement graph (crossing
    multiplicities, connectivity, sep-rank) do NOT cap |bias| (pendant
    stars + exhaustive scan). Only ladders and the "r=1, sep chord
    crossing nothing" class are exactly balanced.
28. In s7_twopoint.c the convention is ρ = γ∘τ (τ first); this differs
    from pf_cl's right_mult by relabelling and is fine for ENSEMBLE
    statistics but do not mix the two codes' per-instance bits.
29. ψ-estimation must use PRODUCTS OF TWO INDEPENDENT fresh samples
    per shared C (unbiased for M²); a single-sample square is biased
    upward by the inner-MC variance.
30. In s7_real_bias.py, k is CAPPED (kcap, default 26) and k>kexh
    instances use MC over subsets (se ~ 0.011 at 2000 samples): fine
    for E|bias| but do not read per-instance biases below ~0.03 as
    exact. r=0 orbits have bias +1/2 EXACTLY (no MC noise: still
    computed by MC — noise ~0.011 around 0.5).
31. JM thinning used: 30000/150000 moves at n=30/50; at n=100 two
    independent runs with thinning 10⁶ and 2×10⁵ agreed within error
    (E|b| .116 vs .112, atom 2.95% vs 2.8%) — thinning ≥ 2n³ is safe;
    the 10⁶ run took ~1h wall (jm.c is the bottleneck at n=100).
    jm.c emits raw bytes (n² per square) on stdout.

## 5. Papers on hand / to request

**IN THE TARBALL (unchanged): CGW_2008.pdf, KPS_2025.pdf,
KSS_large_deviations.pdf, Schramm_2005.pdf, Kwan_Sudakov_2018.pdf +
notes PDF. Keep packing all of them. (kps.txt = pdftotext dump, kept.)**

**To request if needed (session 8):**
- Harer–Zagier 1986 and/or Arratia–Bollobás–Sorkin interlace polynomial
  papers — IF TODO 8a goes route (ii).
- Berestycki or Schramm-style coagulation–fragmentation references for
  route (i) (Schramm_2005.pdf is already in the tarball).
- Cohn–Lempel 1972 original (referee luxury; statement re-verified).

## 6. Quick reference

Master formula: bit = [x̄_{φ(ω)∪B} ∈ col Ā_{φ(ω)∪B}]; now verified on
real squares to n=100. bias = (1/2)E[±1]. Factorization:
4Ebias² = E[M_{d1}M_{d2}] over shared config; ψ(j,d) = E[M_d²] ≈ 0.8/d;
ψ(j,d+1) ≤ ψ(j+1,d). Budget: any c > 1/11 with k = 6log¹¹n. Ladder ⇒ 0;
nested ⇒ 2^{−r}−1/2; cut chord ⇒ (1/2)bias′−1/4; pendant stars ⇒
(1−2^{−(m+1)})^r − 1/2. Real ensemble: E|bias| ≈ 1.5–1.6× iid at same k;
r=0 atom ≈ 0.85/k̄ (6.5%/3.9%/2.9% at n=30/50/100). P(S(λ)) ≥
e^{−O(n log n)}; CGW Thm 3.13 merge ratio ∈ [1/2, 3/2]. q̄ values
unchanged.

*The notes are the living document; keep PROVED/VERIFIED/COMPUTED labels
honest. Session-7 headline: annealed second moment factorizes; iid model
reduced to one L² mixing rate with measured exponent 1.05 vs needed
0.18; template averaging refuted (honest bookkeeping restored); master
formula holds on real squares at n=100.*
