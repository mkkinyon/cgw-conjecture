# HANDOFF: Cameron's second-row conjecture — state after session 5

## >>> SESSION 5 HEADLINE: the MARKED-PAIR IDENTITY (Theorem M, notes §11 sec:marked, PROVED, self-contained) — the post-flip events B₁,B₂, the joining measure, and the switch tie-break are ELIMINATED from the EQ picture. EQ is now: one bit of one uniform square. TODO 6c (α=β≥3 anomaly) RESOLVED. A precise KPS-transfer programme for that bit is laid out in §sec:kps with quantifiers that fit the o(1/log n) budget. <<<

*Prepared 2026-08-28 (end of session 5). Self-contained; the companion
tarball has all code, data, notes (`second_row_notes.tex`/`.pdf`, ~38 pp).
Previous handoffs (`HANDOFF_sessions_1_2.md`, `HANDOFF_session_3.md`,
`HANDOFF_session_4.md`) remain valid except where superseded.*

---

## 0. The problem in one paragraph

As before (see HANDOFF_session_4.md §0). N(σ) = #completions of the 2×n
rectangle (id,σ); C_n(λ) = N/(n−2)!; weak conjecture d_TV(P_n,Q_n)=o(1)
STILL OPEN; refined picture h_n(2) ≈ c₂n⁻³, d_TV ≍ n⁻³.

## 1. What session 5 added

### Theorem M (marked-pair identity; notes sec:marked; PROVED + verified n≤7)

Setting: first row id ⇒ ω = σ. For a merge μ→λ (split part m = α+β,
α,β ≥ 2, **including α=β**), define the marked instance space
X = {(L',P)}: L' ∈ S(λ), P = {j, σ^α(j)}, j in an m-cycle of σ(L')
(unordered pairs once). X = X_A ⊔ X_B by A/B status of P in L'
(B ⇔ rows 1,2 in the same column cycle of P). Then

    J_A = J_B = |X_A|   and   2 C_n(λ)/C_n(μ) − 1 = q̄ = |X_B| / |X_A|.

Proof ingredients (all short, all verified instance-by-instance n=5,6):
- **Arithmetic:** 2|X|/J = 2C(λ)/C(μ) identically (D_λ/D_μ cancellation;
  both α≠β and α=β conventions).
- **Lemma S:** trading a FULL column cycle of pair P inverts π_P on it;
  the row partition into column cycles is unchanged ⇒ the CGW flip at P
  preserves P's own status ⇒ every join lands in X_A; unflips of X_A are
  A-joins; flip at a fixed pair is an involution.
- **Switch toggle lemma:** the CGW switch at a spanning pair precomposes
  π_P with τ₁₂ ⇒ toggles status; involution on join instances ⇒
  J_A = J_B exactly. Composing bijections: X_A ↔ A-joins ↔ B-joins.
- Corollaries: C(λ)/C(μ) > 1/2 unconditionally; **EQ(δ) ⇔
  P_X(marked pair B) = 1/2 + O(δ)** — a single bit of a single square;
  α=β≥3 needs no special bookkeeping (the old rem:albe anomaly was an
  artifact of the joining-side events): verified (3,3)→(6) gives exactly
  3/4 = 2C₆((6))/C₆((3,3))−1. **TODO 6c CLOSED.**

Verification: `pf_bridge.py` (n=5,6: every structural claim
instance-by-instance), `xside.c` (complete enumeration n≤6) and `xside2.c`
(n=7, one σ-representative per type by conjugation symmetry, ~6.6M
completions per type, 2 s each): X_B/X_A = 2C(λ)/C(μ)−1 EXACTLY in all
cases — n=7: 85/86 [(2,5)→(7)], 143/142 [(3,4)→(7), a B-EXCESS],
107/108 [(2,2,3)→(2,5)], 35/36 [(2,2,3)→(3,4)].
(NOTE: full first-row-normalized enumeration at n=7 is 12.2 BILLION
squares — do NOT rerun `./xside 7`; use xside2.)

### Cross-switch calculus (CGW Def 3.9 implemented; `pf_orbit.py`)

The CGW cross-switch is now implemented and verified at n=5 (Latin,
type-preserving, involution, marked-B preserved). Toggle profile: BOTH
shifted pairs (k=+1 and k=−1 translates) toggle deterministically ⇒
**P(B₁ | marked B) = P(B₂ | marked B) = 1/2 exactly** — explains the
extraordinary post-flip balance of §sec:identity. (Historical note now:
B₁,B₂ eliminated by Theorem M.) POST-SESSION UPGRADE: the CGW paper is
now IN THE TARBALL (`CGW_2008.pdf`; Def 3.9 / Lemma 3.10 on pp.
292–294). Lemma 3.10's verbatim toggle family is {ω^k(j), ω^k'(j')}
with INDEPENDENT indices 1≤k<α, 1≤k'<β (the earlier WebFetch
transcription collapsed k=k'). Both shifted pairs are in the family —
P₊₁ at (1,1), P₋₁ at (α−1,β−1) — so their deterministic toggling is a
THEOREM for all n, not just an n=5 observation. (Marked-status
preservation under the cross-switch remains empirical: 2880/2880 at
n=5.) Equal-shift pairs with k ≥ min(α,β) are not covered — matching
the partial toggle rates seen. CAREFUL: the case-1
trade pairs the arc through j (interior 1≤k<α, with j' = ω^α(j)) with
the interior of the δ-path from row 2 to row 1; pairing it with the
1→2 path breaks Latin-ness (bug found and fixed; see pf_orbit.py).

### KPS transfer programme (notes end of sec:kps; NOT yet a proof)

Target bit: rows 1,2 same column cycle of P = {j, σ^α(j)}, L uniform in
S(λ). Programme:
1. **Conditioning is FREE:** intercalates avoiding rows 1,2 switch
   within S(λ) (σ untouched). Stable collections (KPS §7) restricted to
   rows 3..n re-randomise exactly: orbit-uniform measure 2^{−|𝓘₀|}.
2. **Switch action computed:** an intercalate meeting EXACTLY ONE marked
   column with rows r,r' ≥ 3 acts on π_P by precomposition with τ_{rr'}
   ⇒ flips the bit iff r,r' interleave rows 1,2. Meeting both/neither
   columns ⇒ conjugation ⇒ bit fixed. (`pf_toggle.py` has exact n=5,6
   flip statistics; at n=6 strata are far from asymptotic.)
3. **Quantifiers fit:** KPS Lemma 8.3/Def 8.1 (ℓ=log¹¹n): all but ≤ ℓ
   columns have ≥ 6ℓ disjoint stable intercalates (their R-set
   construction, transposed), with rows/partner-columns prescribable in
   βn-sets (so: avoiding rows 1,2 and column j'). Exceptional mark
   fraction O(log¹¹n/n) = o(1/log n) — INSIDE the EQ budget.
   P(S(λ)) ≥ e^{−O(n log n)} for every λ (worst: all-2-cycles) beats
   their e^{−ω(n log²n)} failure prob ⇒ expander holds conditioned on
   S(λ) for every λ.
4. **THE OPEN CORE — "non-linear orbit-mixing lemma":** bit(ω) =
   [0,1 same cycle of π₀·∏_{ω_i=1}τ_i], ω uniform on the cube, τ_i
   disjoint transpositions from stable intercalates through column j.
   Show conditional bit-law = 1/2 + o(1/log n) for all but o(1/log n)
   of orbits. NOT F₂-linear (KPS's kernel computation doesn't apply).
   Note the adversarial-π₀ issue: fixed τ positions can all be
   non-interleaving; but the expander lets the intercalate be CHOSEN
   with rows in prescribed βn-sets — e.g. R = an arc of the π₀-cycle
   between rows 1,2, R* = its complement ⇒ interleaving guaranteed when
   the relevant arcs have length ≥ βn. Short-cycle/short-arc cases need
   a two-step (join-then-toggle) chain argument. This lemma is the
   session-6+ target.

## 2. File manifest (additions; all else as in HANDOFF_session_4.md)

| file | what | re-verify |
|---|---|---|
| `pf_orbit.py` | cross-switch + X-side event stats | `python3 pf_orbit.py 5` (~1 min) |
| `pf_bridge.py` | bridge bijections, q̄ = X_B/X_A | `python3 pf_bridge.py 5` (fast); `6` (~4 min); `624` (~4 min); 633 via `python3 -c "import pf_bridge; pf_bridge.run(6,(3,3),3,3)"` |
| `xside.c` | X_A,X_B all types/merges, full enumeration n≤6 | `gcc -O3 -o xside xside.c && ./xside 6` (2 s); do NOT run at 7 (12.2e9 squares) |
| `xside2.c` | same per fixed σ-representative, n=7 in seconds | `gcc -O3 -o xside2 xside2.c && ./xside2 7 7` etc. (2 s); certify vs `./xside2 6 6` = 10/11, 3/4 |
| `pf_toggle.py` | useful-intercalate / flip stats | `python3 pf_toggle.py 5`; `6` (~8 min) |
| `marked_pair_theorem.md` | working notes for Theorem M | — |
| `section_marked.tex` | notes §sec:marked (input by main tex) | pdflatex second_row_notes.tex |

Key exact numbers: (2,3)→(5): J_A=J_B=X_A=1440, X_B=2880 (q̄=2);
(2,4)→(6): 1520640 thrice, X_B=1382400 (10/11); (2,2,2)→(2,4):
483840, X_B=276480 (4/7); (3,3)→(6): 829440, X_B=622080 (3/4).

## 3. PRIORITIZED TO-DO (updated)

**PROJECT OWNER DIRECTIVE (end of session 5): session 6 should
prioritize TODO 5′ — the marked-pair/KPS route to CGW's conjecture —
over all other TODOs. Clarification from the owner: this is NOT an
abandonment of the refined n⁻³ picture (TODO 4 etc.), which retains
real merit; it is a judgment that the CGW Conjecture is the more
immediately accessible target right now.**

### TODO 5′ (THE strategic target, now sharply scoped): the non-linear
orbit-mixing lemma (§1 item 4 above). Sub-steps: (a) formalize the
column-transposed version of KPS's R-set construction (their Lemma 8.5
proof) to get "all but log¹¹n columns carry ≥ 6log¹¹n disjoint stable
intercalates avoiding two given rows and one given column"; (b) the
conditional-expander statement given S(λ) (union-of-costs; check their
Lemma 4.7 machinery is not needed — direct conditioning may suffice
since P(S(λ)) is only e^{−O(n log n)} small); (c) the mixing lemma for
the cycle bit under random disjoint-transposition subsets with chosen
anchors; (d) assemble EQ(o(1/log n)) ⇒ weak conjecture via
thm:reduction. Each of (a)–(c) is a self-contained lemma; (c) is the
mathematical core.

### TODO 1″ (referee pass): unchanged (§7 AND §8; the layer-5
j-resolved strata scan safety net still does not exist). ADD: referee
the new §sec:marked (the proofs are short; check the α=β arithmetic
D_λ/D_μ = α²μ_α(μ_α−1)/(2mλ_m) and the unordered-pair conventions
once more against pf_bridge counts).

### TODO 2′ (three-row note), TODO 4 (bulk profile γ / c₂): unchanged.
Note TODO 4 gets new leverage: P_X(B) = q̄/(1+q̄) is exactly computable
from cgw_data to n=11, and the X-side is directly enumerable/samplable
— a cleaner observable for transfer-matrix / cluster expansions than
q̄ itself.

### TODO 6 (cleanups): 6c is CLOSED (Theorem M). 6a (algebraic
involution-lemma proof) is now LOW value (B₁,B₂ eliminated; the
cross-switch section supersedes). Others unchanged.

## 4. Pitfalls (all previous, plus)

16. The cross-switch case-1/case-2 arc↔δ-path pairing (see §1 above);
    re-verify with `pf_orbit.py 5` (toggle profile must be 1.0000 at
    k=1 and k=α+β−1) after any edit.
17. `do_switch` breaks the first-row-id normalization — ALWAYS renorm
    (relabel symbols so row 0 = id) before comparing squares as keys
    (bug found in pf_bridge; fixed with `renorm`).
18. In X-side counts with α=β, ordered marks double-count unordered
    pairs — halve (xside.c does; pf_bridge dedups by frozenset).
19. The F-strata of pf_toggle (flip-count strata) are NOT preserved by
    switching and do NOT cleanly carry the imbalance at small n; don't
    read n=6 strata as asymptotics.
20. Direct downloads (curl/pip) are blocked by egress policy for
    users.monash.edu.au AND arxiv.org; WebFetch summarizes through a
    small model and CANNOT retrieve verbatim math reliably (it
    collapsed Lemma 3.10's independent indices k,k' into k=k').
    All needed papers are now IN THE TARBALL — read them with
    pdftotext or the Read tool; ask the project owner for any others.

## 5. Papers on hand / to request

**IN THE TARBALL (post-session uploads by the project owner — the
session-1 literature had silently dropped out of the tarball chain;
do not let that happen again: pack every provided paper into each
session's outgoing tarball):**
- `CGW_2008.pdf` — Cavenagh–Greenhill–Wanless, RSA 33 (2008). THE
  central prior work: exact tables (§4), switch/flip (Defs 3.2/3.6,
  Lemma 3.7), cross-switch (Def 3.9, Lemma 3.10, pp. 292–294),
  multigraph bookkeeping (§3).
- `KPS_2025.pdf` — Kwan–Petrova–Sawhney, arXiv:2509.13125v2. Stable
  intercalates §7, expander Lemma 8.3/Def 8.1 (ℓ=log¹¹n), master
  theorem §6, approximation Lemma 4.7 / Cor 5.5.
- `KSS_large_deviations.pdf` — Kwan–Sah–Sawhney, Bull. LMS 54 (2022),
  arXiv:2106.11932v3. Intercalate concentration N ~ n²/4 and the
  large-deviation rates (upper tail exp(−Θ(n^{4/3}(log n)^{2/3}))).
  Baseline density facts for stable-intercalate counting.
- `Schramm_2005.pdf` — Schramm, Israel J. Math. 147 (2005),
  compositions of random transpositions: cycle structure after cn
  transpositions → PD(1); coagulation–fragmentation invariance. The
  classical toolbox closest to the orbit-mixing lemma (TODO 5′(c)) —
  our orbit dynamics on π_P is coagulation–fragmentation driven by a
  uniform SUBSET of DISJOINT transpositions with chosen anchors.
- `Kwan_Sudakov_2018.pdf` — Kwan–Sudakov, RSA (2018): intercalates
  and discrepancy; the n²/4 lower bound machinery.

**Still to request if needed:**
- Kwan–Sah–Sawhney–Simkin, "Substructures in Latin squares",
  arXiv:2202.05088 (Israel J. Math. 2023) — the primary
  triangle-removal-process comparison paper (KPS's Lemma 4.7 is their
  new easier variant). Request ONLY if TODO 5′(b) (conditional
  expander given S(λ)) turns out to need the original machinery
  beyond KPS §4–5.
- Riordan book, Häggkvist–Janssen 1996, Cameron 1992 (session-1
  uploads, still absent from the current chain); Touchard 1953 /
  Riordan 1954 / Moser 1967 / Whitehead 1979 / McKay–Wanless JCTA
  1999 — same low-priority status as before.

## 6. Quick reference

q̄ = X_B/X_A values (exact, from cgw_data via Theorem M): n=7:
85/86 (2,5→7), 143/142 (3,4→7) [note >1], 107/108, 35/36; n=11
(2,9→11): 0.998080. P_X(B) = q̄/(1+q̄).

*The notes are the living document; keep PROVED/VERIFIED/COMPUTED labels
honest. Session-5 headline: Theorem M (marked-pair identity) — EQ is now
one bit of one square; the KPS-transfer programme for that bit is scoped
with its single open core lemma.*
