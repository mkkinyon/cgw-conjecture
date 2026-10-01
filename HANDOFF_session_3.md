# HANDOFF: Cameron's second-row conjecture — state after session 3
## >>> TODO 1 IS NOW CLOSED: the layer-4 theorem is UNCONDITIONAL. <<<

*Prepared 2026-08-28 (end of the third working session). Fully self-contained
for a fresh session. The companion tarball contains all code, data, and the
working notes (`second_row_notes.tex` / `.pdf`, 31 pp). Read the notes'
abstract and ToC first; this file is the operational companion. The previous
handoff (`HANDOFF_sessions_1_2.md`, kept in the tarball) remains valid except
where superseded below.*

---

## 0. The problem in one paragraph

Let L be a uniformly random Latin square of order n with first row = identity,
and let σ (a derangement) be its second row. Cameron conjectured (blog, "A
niggling problem", 24 Jan 2015) that the law of σ tends *very rapidly* to
uniform on derangements. N(σ) = #completions of the 2×n rectangle (id, σ);
C_n(λ) = N(σ)/(n−2)!; P_n = law of the cycle type of σ; Q_n = that of a
uniform derangement. Weak form (CGW08 Conj 6.1): d_TV(P_n,Q_n) = o(1).
**Still open.** Refined picture (notes §2–3): log C_n additive over cycles,
h_n(2) ≈ c₂n⁻³, d_TV ≍ n⁻³.

## 1. What session 3 added (on top of sessions 1–2)

Session 3 attacked TODO 1 (layer-4 theorem) and TODO 2 (three-row Poisson).

**NEW — the switch calculus (notes §7, new section "An exact switch calculus
for the fourth layer"):**

- **Adjacent representatives.** σ = (12…n) and σ′ = (12)(34…n) have boards
  differing in exactly 4 cells (a 2×2 switch): common board B₀ = two paths
  P¹ (3 edges: c1–r1–c2–r2), P² (2n−5 edges: c3–r3–…–cn–rn); pin cells
  d₁=(2,3), d₂=(n,1) (σ's), e₁=(2,1), e₂=(n,3) (σ′'s). (1-indexed. In the
  code, 0-indexed: d1=(1,2), d2=(n−1,0), e1=(1,0), e2=(n−1,2).)
- **Switch decomposition (Lemma, exact, verified n=6,7):**
  ΔN⁽⁴⁾ = 2Br_G + 2[G₂(e)−G₂(d)] + 2[H(e)−H(d)],
  Br_G = G(d₁)+G(d₂)−G(e₁)−G(e₂), G(z) = #{disjoint good pairs, z ∈ τ}.
- **Path Wronskian (Lemma, proved):** p_a p_b − p_{a−1}p_{b+1} =
  −(−x)^{b+2} p_{a−b−3} (p_e = path matching polys, Chebyshev continuation
  for negative index). Consequences W_C = −x³p_{C−4}, V_C = −x⁴p_{C−5};
  cycle-to-path: M_{C_{2m}} = p_{2m−1} + x p_{2m−3}.
- **Empty bracket (Prop, proved):** Br(∅) = M_{n−4} exactly — a switch-local
  proof of the ℓ=2 fusion identity. Verified n=6..14.
- **Ψ-symmetry (Remark):** Ψ(i,c) = (3−c, 3−i) mod n preserves B₀, swaps
  d₁↔d₂, fixes e₁,e₂; explains the exact G(d₁)=G(d₂) seen at n=8..11.
- **Valuation lemma (proved):** the 4-term bracket polynomial B_M has
  [x⁰]=[x¹]=0 always and [x²] = (1_{c2∉M}−1_{cn∉M})(1_{r3∉M}−1_{r1∉M});
  verified exhaustively (|M|≤2, n=9,10).
- **Master formula (proved):** each pin-difference is a single x^s-shifted
  path product via component Wronskians; generic form
  B_M = −x³p_{C−4}p_D∏mid − x⁴p_{C−5}p_{D−1}∏mid. Verified exactly on
  random M at n=9,11,13.
- **Uniform punctured estimate (Lemma, proved, Bonferroni):** for any
  max-degree-2 board (paths ∪ cycles), ⟨M_F⟩_m = m! e^{−E/m}(1+O(1/m))
  uniformly for E ≤ 2m. Plus weight/moment corollaries. (This was TODO 1
  step 1 AND the missing ingredient (a) of TODO 2 — done once, well.)
- **Generic stratum (Prop, proved):** window-avoiding M give
  ⟨B_M⟩ = (M_{n−4}/M_n)R₀(M)(1+O((j+1)/n)).
- **Pairing lemmas (proved):** classes {one pinned row + one pinned col}
  and {M ∋ pin cell}, and the G₂- and H-differences, are all
  O(n⁻¹(M_{n−4}/M_n)N₀²) — via line-swap pairings with EXACT pinned-factor
  equality (same deleted row/col sets!) and x³-Wronskian free-factor
  differences. Numerically: the two (c-2) classes total +4.277 vs −4.276
  (units M_{n−4}N₀) at n=8.
- **Column-exchange lemma (proved, exact bijection; verified):** for M
  meeting rows 1,2 (col 1 constraint-free): R₀(M∪d₂) − R₀(M∪e₂) =
  −R₀(M∪{(3,1),(n,3)}) EXACTLY; dual for rows 2,3 (sign +, pins
  (1,3),(n,1)); for rows 1,2,3 the bracket vanishes identically.
  This resolves the dangerous x¹-classes: pairing (1,γ)↔(3,γ) has exact
  pinned-factor equality and O(1/n) free-factor Wronskian differences
  (measured ≈ −2/n²). `switch.py check5` verifies all of it.
- **Theorem (layer 4; UNCONDITIONAL, proved in notes §7):**
  Δlog N⁽⁴⁾ = 2M_{n−4}/M_n + O(n⁻⁵), hence r₄ = r₃(1+O(1/n)).
  Every step is an exact identity, an exhaustively verified finite case
  analysis, or a uniform Bonferroni estimate. The crude-bound steps
  (orphan sub-cases, j > n/2 tails) are written at working-note rigor —
  a referee pass is TODO but no gap is known.

**NEW — loops remark (notes Remark 1.x, §1 'Known results'):** a loop's
Cayley table is a REDUCED square, so its second row satisfies σ(1)=2; but
conditioning on σ(1)=2 leaves the cycle-type law exactly P_n (within-class
uniformity + σ(1) uniform on {2..n} on every class), and the TV distance of
the conditioned laws equals d_TV(P_n,Q_n) exactly. So all results transfer
verbatim to loops. Verified by complete enumeration at n=5,6
(`loop_reduction_check.py`).

**NEW — three-row section upgraded to full proofs (TODO 2 essentially DONE):**

- thm:poisson (joint Poisson(1/2)³) now has a complete proof: generic
  configurations via the uniform estimate; collision patterns classified
  with a net-exponent argument (every pattern loses ≥ 1 power of n; a full
  B-pair/C-pair coincidence is impossible since it would force σ(i)=i).
- **Fusion expansion with remainder (Lemma, proved by induction):**
  N⁽³⁾(λ) = M_n + Σ_ℓ λ_ℓ M_{n−2ℓ} + E(λ), |E| ≤ Cκ²(n−8)!.
  Conventions M₀=2, M₁=−1. Numerically max|E|/(κ²(n−8)!) ≈ 0.06 at
  n=12,16,20; E ≥ 0 except type (n/2,n/2) where E=−2. (Also proves
  second-order additivity of log N⁽³⁾.)
- prop:tilt now proved with O(n⁻²): expansion lemma + derangement cycle
  moments E_Q[∏(λ_ℓ)_{s_ℓ}] = ∏ℓ^{−s_ℓ}(1+O(1/(n−m)!)).
- prop:cov now proved: Cov(X₁₂,X₁₃) = (1/2+O(n⁻¹log n))/n², via the exact
  admissible count C(n,2)−n+λ₂ and a JOIN-STABILITY bound
  |Δ(λ)−Δ(μ)| ≤ Cn⁻³ for adjacent types, proved with the switch calculus
  (works for ANY join, not just (2,n−2)↔(n)).

**Everything from sessions 1–2 still stands** (exact switching identity +
reduction theorem; pair-IE; JM numerics; q̄ tables; see
`HANDOFF_sessions_1_2.md` §1, whose PROVED/VERIFIED labels remain accurate).

## 2. File manifest (additions; everything else as before)

All previous files unchanged and still passing. New:

| file | what it is | re-verify with |
|---|---|---|
| `switch.py` | switch calculus: B₀, general punctured counts on any path/cycle board, empty bracket, brute-force checks, ratio scans, column-exchange suite | `python3 switch.py 1`; `python3 switch.py 236 237` (~2 min); `python3 -c "import switch; switch.check5()"` |
| `fusion_lemmas.py` | symbolic Wronskian proof-checks + exact master-formula verification | `python3 fusion_lemmas.py` |
| `valuation_scan.py` | exhaustive valuation-lemma check (all M, j≤2, n=9,10) | `python3 valuation_scan.py` |
| `switch_j1.py` | strata sums S_j vs enumeration; j=1 cell scans; python anatomy(n) | `python3 switch_j1.py signed 8 3` (S₀ must be 12406 = M₄N₀) |
| `anatomy.c` | exact Br_G, N₀⁽⁴⁾, m±, E±[A]: full enumeration + per-τ Ryser | `gcc -O3 -o anatomy anatomy.c && ./anatomy 9` (0.4 s; n=11 ~3 min) |
| `anatomy11.txt` | n=11 anatomy output | — |
| `HANDOFF_sessions_1_2.md` | the previous handoff (still-valid details) | — |

Key new exact numbers: Br_G = 3054 (n=8), 177688 (n=9), 11284750 (n=10),
921606048 (n=11); Br_G·N₀/(M_{n−4}N₀⁽⁴⁾) = 0.8919, 0.8685, 0.8818, 0.8929;
m⁺(E⁺−E⁻)/(M_{n−4}E[A]) ≈ −1.17/n (n=9,10,11); the (c-1,A″=0) "worry
class" by_j at n=8: {2: −4110, 3: −1113, 4: +4465}, total −758.

## 3. PRIORITIZED TO-DO (updated)

### TODO 1′ (light): referee pass on the layer-4 proof

The theorem is unconditional. Remaining polish: (a) re-derive the
case inventory of Lemmas "pairing" and "boundary" independently and
check nothing is missed (the exact class-resolved sums in `switch_j1.py`
are the safety net — any missed class shows up numerically); (b) tighten
the orphan depth-≥2 crude bounds into displayed inequalities; (c)
consider extracting §7 + §4 as the standalone paper "The second row of a
random Latin square" (old TODO 6f(i)).

### TODO 2′ (cleanup): polish the three-row material into a standalone note

Proofs are in. Remaining: (a) make prop:cov's join-stability step a
standalone lemma (it currently leans on §7 across a section boundary);
(b) novelty check ("intercalates in Latin rectangles", "Poisson" — nothing
found so far; Riordan's K(3,n) theory is aggregate-only); (c) consider the
paper split of old TODO 6f.

### TODO 3 (triple calculus / layer 5), TODO 4 (profile γ), TODO 5 (EQ /
KPS machinery), TODO 6 (cleanups): unchanged from `HANDOFF_sessions_1_2.md`.

Notes: for TODO 3 the switch calculus applies verbatim to N⁽⁵⁾ (pair counts
→ triple counts; same pins, same Wronskians — only Lemma "uniform" reuse
plus a triple-IE). For TODO 6a, try deriving the row-swap involution lemma
from the Ψ-involution of notes §7 Remark (they look closely related).

## 4. Pitfalls (all previous ones still apply, plus)

7. **Indexing**: the notes are 1-indexed, all code 0-indexed. Pinned lines:
   rows {2,n}, cols {1,3} (1-indexed) = rows {1,n−1}, cols {0,2} (0-idx).
8. **Bracket conventions in strata sums**: R₀†(M∪z) = R₀(M) if z ∈ M, 0 if
   z's row/col meets M; forgetting this makes pin-cell classes explode
   (the ±854 artifact in the j=1 scan was exactly this).
9. **Do not "improve" the pairing maps**: the pinned-factor equality is
   exact only because the deleted row/col SETS coincide; perturbing the
   map breaks the exactness the argument needs.
10. **Expected sizes**: (E⁺−E⁻)/E[A] is ~ −0.85/n⁴ (constant −1.17 in the
   m⁺-normalisation), not n⁻⁵; the n⁻⁵ arises only after multiplying by
   m⁺ and subtracting inside Br_G. Small-n (n=8) values are off-trend.
11. `switch_j1.py signed` is exponential in j (combinations of all cells);
   keep j ≤ 3–4 and n ≤ 10, or restrict to a class as in the worry-scan.

## 5. Papers on hand / to request — unchanged from previous handoff.

## 6. Quick-reference numbers — previous handoff §6 plus §2 above.

*The notes are the living document: keep every numerical claim backed by a
script in the tarball, and keep the PROVED / VERIFIED / COMPUTED labels
honest. Session-3 headline results: the layer-4 theorem r₄ = r₃(1+O(1/n))
(unconditional, via the new switch calculus) and full proofs for the
three-row Poisson/tilt/covariance package.*
