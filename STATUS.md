# STATUS — CGW Conjecture project (rolling; replaces per-session handoffs)

*Updated 2026-10-03 (session 13, continuous; see SEGMENT LOG below for today).  Two refereed documents: `paper/cgw_gap_note.tex`
(8 pp, to send to CGW) and `paper/weak_cgw_conditional.tex` (29 pp, REVISED 2026-10-03 twice: hyp:witness is the primary hypothesis again (§5); hyp:orbitgen / hyp:adaptive / hyp:firstpair are ALTERNATIVE sufficient conditions (§6, incl. new §6(e) first ladder pair); title 'a reduction'; assembled by `cprime_blocks.py` + `witness_subs.py` +
`paper/assemble_conditional.py` + `postprocess_conditional.py`; thm:marked, prop:tailsurvive,
thm:reduction [error O(δ log n + n^{−1+o(1)})], thm:ladder, offsets, prop:trapped, prop:splice (NEW,
refereed: P_X[B,d₁₂=ℓ,|C₁|<n/2+ℓ] ≤ 4/(n−2ℓ+1)), orbit method, hyp:orbit ⇒ (L), What remains).
Referee corrections applied to both (see git log a9d7507 and after).  The notes
`second_row_notes.tex` (139 pp, 0 errors) are the living document; this file
is the two-page map.  History: `git log`.  Old handoffs `HANDOFF_session_*.md`
are kept for the record; their TODO lists are superseded by this file.*

## ASSESSMENT (2026-10-02 evening; own + independent subagent) — read this before the ALERT
**Verdict.** The PROVED base grew and is solid (thm:marked, thm:ladder, offsets, prop:trapped, prop:splice,
prop:tailsurvive, thm:reduction independent of CGW, gap note + constant-2 repair).  The UNPROVED residual has
been reformulated, not reduced: hyp:witness (one rare local event at a frame-selected pair) → hyp:orbit (constant-
strength genericity (β), introduced by the orbit method itself, found gap-direction/unprovable in X) → (U) →
hyp:adaptive (rare local events again, C log log n of them).  hyp:adaptive is of the SAME logical type as
hyp:witness and logically STRONGER (needs many events, not one); the loop (β introduced → removed) was presented
as progress twice.  The weakest sufficient condition we have is still hyp:witness.
**Core obstacle, precisely:** a lower bound (any constant, any rate) on the probability of ANY local event at a
row pair selected by the frame, inside X.  Minimal instance: P_X[(x₁,y₁) flippable] ≥ c, x₁ = π⁻¹(1), y₁ = π⁻¹(2).
Every exact identity in X is fair; gap-direction bounds need a rigid inverse (exists only for min(α,β)=2);
type-conditioning transfers at cost e^{O(n log² n)} (vdW/Brégman) while the heuristic failure of "long ladder,
nothing flippable" is 2^{−n/log² n} — too large by log⁴n in the exponent, and tightening the transfer is
conjecture-strength (the subagent's "cheap conditioning at n^{O(log n)}" used the wrong direction of thm:marked:
the surviving half bounds fine types from ABOVE only).  Data speak to typical λ only; the hypotheses could fail
on rare λ unseen.
**Over-claims to fix:** paper title "conditional proof" → "a reduction"; abstract "test its ingredients up to
n=100" while the regime is out of reach; "(Π) SOLVED", "no genericity constant enters" sold as gains;
"refereed" = by subagents.  **Publishable now:** gap note (+ Prop 3); the identities-plus-reduction paper.
**Next (in order):** (1) DONE 2026-10-03 (hyp:witness primary); (2) DONE 2026-10-03 (fixed-rectangle sampler; no
falsification at n ≤ 100); (3) target P_X[(x₁,y₁) flippable] ≥ c — reduced one-way to P_X[G] ≥ 2c (prop:firstpair), same
obstacle; and P_X[B] itself is controlled by the same G (prop:firstpairbias); (4) see NEXT_STEPS_2026-10-03.md.
**Progress =** a proven lower bound at a frame-selected pair, or a falsification.  Another "lower tail of rare
local events at frame-selected rows" = reformulation, to be labelled as such.

## LITERATURE CHECK 2026-10-03 (MathSciNet + Google Scholar citers of CGW08; KSSS, Allsop–Wanless, DKKS, GMW, CW16, Allsop–Morris, KPS read)
Notes §(m).  (a) Fixed-cell bounds in the literature are rectangle statements (KSSS Lemma 3.3, k ≤ n/4, within-row swaps; DKKS
Thm 1.6, k < n/2, two-sided (1±δ)/n for εn-sparse P, Granet–Joos long switches in the free-pair expander) or open conjectures for
sparse P in full squares (DKKS Conj 1.5; Kelly's (e²+o(1))/n spread).  Our (FC) conditions on two FULL rows (dense): in spread
language P_X[Q ⊆ L] = P[R∪Q]/P[R], and a one-sided numerator bound is useless without P[R ⊆ L] = C(λ)/|L_n| to n^{O(1)}, a weak
CGW.  Allsop–Wanless (full squares, cycle switching): Lemmas 3.8–3.10 are fixed-cell bounds for boundary cells; interior cells get
the trivial bound; the forward-degree count for column cycles avoiding prescribed rows "seems like a difficult task in general"
(their §3.1) — the obstacle in their words.  (b) Nobody states or reproves a row-pair cycle bound independently of CGW §3.
(c) **GMW25 Theorem 9 (hence their average-case Theorem 1 via Theorem 12) uses CGW Cor 4.5 (upper half) and Thm 4.9 — both
affected.  REPAIRED (refereed): prop:tailsurvive(ii) with k = ⌈√n⌉ gives P[κ ≥ 9√n] ≤ exp(−(½−o(1))√n log n) (so Thm 4.9's
STATEMENT survives — the gap note previously said it did not; corrected), and C(finer) ≤ 2C(coarser) iterated gives P_n(λ) ≤
(n/2)Π(2/c_i) (Cor 4.5 upper half with n/2 for n^{1/3}); GMW's proof goes through as stated.**  Gap note §4 and paper rem:gap
updated (commit 828c235).  **Cavenagh–Wanless 2016 (parity equidistribution) read: its Theorem 1.1 uses CGW Cor 4.5 (Lemma 2.1: needs P[(n)] = o(1),
NOT repaired), Thm 4.9 (Lemma 2.2: repaired, modulo writing prop:tailsurvive in Ω(F)), and the UPPER bound of Thm 3.13 at odd
splits min ≥ 3 (Lemma 2.3: the gap direction, outside Prop 3; what they need is only an averaged/any-constant version).  So CW16
currently depends on the gap at two points.  NEW NEGATIVE KNOWLEDGE (refereed): from the surviving material + our identities we
know of no proof that two rows of a random Latin square form a single n-cycle w.p. o(1) (Lemma 4.4's 2n^{−2/3}); second moment
(KS18 + prop:tailsurvive(i)) gives only P[(n)] ≤ 5/6 + o(1); Prop 3 iterated loses 4^k k!.  Gap note §4 updated (8 pp).**
**Google Scholar citers (user's list) + Allsop–Morris 2026 and Kwan–Petrova–Sawhney 2025 read (notes §(n), §(o); refereed).**
*Allsop–Morris (arXiv 2606.18174), Theorem 1.1: two-sided (δ/n)^{|P|} ≤ P[L ⊇ P] ≤ (Δ/n)^{|P|} (δ→1/23, Δ→23) for partial Latin squares P
with |R_P| ≤ αn, |C_P| ≤ βn, 2α+β < 1 — full squares; our R (two full rows) is outside the STATEMENT.  Their tool, the η-switching
η_L(r,r′,c) (displaces the target cell within its row; changes only rows r,r′,r″; reversibility ≤ n preimages if Type One, ≤ 1
otherwise), is LEGAL in X whenever r,r′,r″ ∉ {1,2} ∪ R_Q, with a DETERMINISTIC forward degree n − O(|Q|) — the first switching
available to us in X with that property (all our turns have random forward degrees).  Their Theorem 3.1 (P[Type One | P] ≤ D/n),
which converts reversibility into the per-cell bound, transfers to X in Claim 1 (recursion; bookkeeping α=(2+|R_Q|)/n, β=(|C_Q|+2)/n,
excluded columns C_Q ∪ {p,p′}) and Claims 2–3 (cross-switch segments restricted to contain 0 or 2 of rows 1,2; double count must be
rewritten; constants only), and FAILS at exactly one step: Claim 4(ii) (the column-cycle switch γ_L(c₂,c₃,x) when the row switch
ρ(r₃,x,c₃) hits c₂) in the configuration "exactly one of rows 1,2 on γ_L(c₂,c₃,r₁) and x on the (c₂,c₃)-cycle of the other"; the
other cases are repaired by turning unions of cycles.  Claim 1 alone gives only P_X[Type One | Q] = O(n^{−1/2}), i.e. C/√n per
cell, not (FC).  If (FC) in X held with constant Δ: P_X[ν ≤ M] ≤ 8Δ^{2M+2}/n (notes §(l)(i)), the first term of hyp:firstpair.
Classification: NEGATIVE KNOWLEDGE on one step + a sharp question for Allsop/Morris/Wanless: does Theorem 3.1 hold when P contains
two complete rows the switchings never touch?  AM cite CGW08 only for the §3 switching lemmas (reversibility/effect of cycle and
cross switches) — unaffected by the gap.*
*KPS (arXiv 2509.13125) prove Cameron's parity conjecture (Thm 1.3) via an approximation lemma (triangle-removal process → random
subset of a random Latin square at cost e^{O(n log² n)} on the failure probability) + a canonical switch-invariant family of
"stable intercalates" that meets every large set of rows simultaneously (Lemma 8.3).  Consequences: (i) Thm 1.3(1) ⇒ P[rows 1,2
same parity] = ½ + o(1) ⇒ **P_n[σ even] = ½ + o(1): the parity marginal of the weak conjecture is proved (by them)**; supersedes
CGW Thm 4.3 (gap-affected); (ii) CW16's theorem is a special case of Thm 1.3(4), whose written proof uses CW16 as a black box
(dependence CGW→CW16→KPS 1.3(3)–(4) formally present), Remark 6.6 sketches an independent re-proof — gap note §4 updated
(commit pending); (iii) program sketch (not a result): approximation lemma for completions of a fixed R + stable toggles acting on
flippability bits + re-randomisation; the separation input (β) reappears as "some stable intercalate at row x_k has a separated
column pair".  Classification: literature / WRITING (two gap-note corrections).*

## COLD READ 2026-10-03 (independent subagent, after 6 segments) — verdict and corrections to adopt
**Verdict.** Core obstacle unchanged; the log says so honestly.  New and unconditional today: prop:firstpair, prop:firstpairbias
(found by a referee), lem:nuswitch, the validated sampler.  hyp:firstpair is logically incomparable with hyp:witness; the
residual is now "hyp:witness or hyp:firstpair or ...", weakly smaller as a statement, not smaller in difficulty.  Seg 3 is a
relabeling (eq:minimal is sufficient for nothing); Seg 6 turned a switching-friendly upper bound into a lower tail (neither
easier nor harder).  Recommendation: stop in-house work on the residual; publish (gap note first; human read of §§2–4;
one consultation on the fixed-position question).
**On hyp:firstpair vs hyp:witness.** Better in two defensible ways, undersold in the text: (1) its failure is dominated by SHORT
ARCS, not by a large count: split at M = C log log n, P[Gᶜ] ≲ P[ν ≤ M] (≈ 4M/n; an upper bound on a short structure at the
selected pair) + 0.6^M (a conditional "hazard ≥ c" statement over O(log log n) candidates) — the "Θ(n) candidates" framing in the
paper describes the wrong quantity and should go; (2) its data measure the hypothesised quantity itself (n·P[Gᶜ] flat at 5.5 for
n = 30–150), whereas hyp:witness's ℓ′ = log⁵n exceeds n at every sampled n.  Hidden costs: the ν-tail is the dominant term (R4:
no easier); the candidates are structurally dependent (generated by ρ_{x₁,y₁}; measured dependence positive — the dangerous
direction); admissibility couples G to rows 1,2; the sup over all λ forces uniformity over every rectangle — CHECK whether
thm:reduction needs all λ or only types with O(log n) parts (cheap weakening if so).  Presentation: lead with hyp:firstpair, keep
"neither implies the other; hyp:witness is the formal reference".  State the trivial common weakening
sup min(P_X[Flip=∅], P_X[Gᶜ]) = o(1/log n).
**On the obstacle statement.** "A lower bound on any event at a frame-selected pair" is imprecise (trivial lower bounds exist;
P_X[A] ≥ ¼ is one; R4 shows the UPPER direction stuck too).  Sharper: constant-factor estimates, in either direction, for events on
O(polylog) cells in two fixed rows of a uniform completion of a fixed 2×n rectangle plus O(1) prescribed cells (row exchangeability
makes "frame-selected" exactly two prescribed cells at cost (n−2)(n−3)).  Concretely a forward-degree lemma: for fixed x, c and a
set S of k columns, at least θ_k n rows u have their ρ_{x,u}-cycle through c avoiding S, with failure o(1/log n).  "Frame-selected"
is a red herring once the fixed-cell form is used.  The generic weakness of Latin-square switchings: the only always-legal moves are
trades along whole row cycles (lengths ≈ uniform), so the number of valid switchings is itself a random count needing a lower tail;
one cell per row gets a deterministic 1/n, two cells in a row do not.  Question for an expert: spread-type bounds P[Q ⊆ L] ≤
(C/n)^{|Q|} for small Q in completions of a fixed rectangle (KSS(S)?).
**Process lapses.** Rule 1: "REDUCTION" is not one of the five classes (Seg 3); Segs 4, 6 stacked labels letting PROVED lead.
Rule 3: Seg 4 bolded the implication before "Cost first".  Rule 4: Segs 2–4 had no per-segment "what did not change".  Rule 5:
NEXT_STEPS calling itself the cold read was self-certification.  Paper sentences to fix: §7 "the difficulty is ENTIRELY that the
positions are selected" (overstatement, contradicts the next sentences); §7 "smallest sufficient condition we have" (vs
incomparability) → "simplest"; abstract "reduced to a witness lemma" (four incomparable hypotheses); abstract "lower bounds of the
kind whose proof in CGW has a gap" reads as insinuation; abstract must flag the sampler's unproved connectivity/mixing; item (6)
"Θ(n) candidates"; rem:fixrect "[0.995,1.005]" is a range of point estimates, not a confidence interval.
**Rule amendments adopted:** the classes are PROVED / ONE-WAY REDUCTION (to a statement not known to be easier) / REFORMULATION /
NEGATIVE KNOWLEDGE / DATA / WRITING, exactly one per segment (a segment with several outcomes is split); the cold read is done by
a subagent that has not seen the session, never by a document the lead writes.

## WORKING RULES (added 2026-10-02; every segment reads and obeys these)
1. Every segment report classifies its outcome as exactly one of: PROVED (unconditional), REFORMULATION,
   NEGATIVE KNOWLEDGE, DATA, WRITING.  No superlatives in STATUS or summaries ("headline", "kills", "solved"
   are banned unless the statement is unconditional and refereed).
2. A new hypothesis may REPLACE the reference hypothesis (currently hyp:witness) only if it is proved to be
   implied by it or strictly weaker; otherwise it is listed as an ALTERNATIVE sufficient condition.
3. The cost of any new formulation (rarity, constants, regime, what got worse) is stated in the first paragraph,
   before the benefit.
4. Each report ends with "What did NOT change": the core obstacle restated, and whether this segment touched it.
5. Every ~5 segments, an independent cold-read assessment (subagent given STATUS + paper, asked "shrinking or
   reformulating?"), recorded here.

## SEGMENT LOG 2026-10-03 (classified per the working rules)
**Seg 1 — WRITING.** Paper restructured: hyp:witness primary (§5, with (i) G=0 form, uniform in λ,α,β, row-exchange
transport for d₂₁), hyp:orbitgen / hyp:adaptive as alternatives with corrected offset parts (scaled thresholds, truncated
diagonals); §7 states the minimal instance.  Commit 72c6689.  Did not touch the core obstacle.
**Seg 2 — DATA (no falsification).** `jmfix.c` + `s13_fixrect.py` + `s13_fixstats.py`: JM chain restricted to completions of
a FIXED 2×n rectangle (simple random walk on the restricted move graph; uniform on the component).  Validated: n=5,6 all
completions hit from two starts, χ² within 1σ (36 / 4032 / 5376 states); n=7 P[B|R]=0.4965 vs exact 0.497076.  Types (n),
(n/2,n/2), (n−2,2), (10,10,10), (4,2^13), (6,3^8) at n=30 (60k squares), (50), (25,25), (48,2), (4,2^23), (6,3^14,2) at n=50
(40k), (100), (4,2^48) at n=100 (12k); marks α=2, ≈m/3, m/2: P_X[B] ∈ [0.4974,0.5024] (n≤50; sampling error ≈0.002–0.0025),
hence C(λ)/C(µ) ∈ [0.995,1.005] for every tested pair incl. (4,2^13)→(2^15), (6,3^8)→(3^10); P_X[(x₁,y₁) flippable] ∈
[0.4988,0.5017] (0.4896 at n=7: ½ is not an identity); P[Flip=∅|K₀=k] = 2^{−(k−1)}; witness frequencies at pair 1 identical
across types to ±0.002.  Says: first-moment quantities are type-independent at n ≤ 100.  Says nothing about o(1/log n) tails.
Logs runs/s13/fixrect/; notes §(k).
**Seg 3 — one-way REDUCTION of the minimal instance to a transposed analogue (prop:firstpair, refereed; NOT a lower bound).**
For the first ladder pair alone, μ ≥ 2 always (t=1 is always a candidate) and the adaptive turn needs no short-cycle or
collision condition: with G = {some t<μ: x₁,y₁ in different cycles of the column permutation of {ρ^t p, ρ^t p′}, and the
cycle through x₁ or y₁ meets {1,2} in 0 or 2 rows}, the turn at the least good t is an involution of X∩G toggling flip₁:
|P_X[flip₁] − ½| ≤ ½P_X[Gᶜ], and P_X[A], P_X[B] ≥ ¼P_X[G] for EVERY split.  Data: P[G₁] ≈ 0.41, P[G, t≤8] = 0.81/0.88/0.94
at n=30/50/100, identity P[flip₁∩G]=P[¬flip₁∩G] holds to ±0.002.  Cost/what it is: G is flippability with rows and
columns exchanged at square-selected positions (admissibility clause added; conditioning still on rows 1,2); the implication
is one-way; a lower bound on P_X[G₁] (constant-probability local event, not rare) is the same obstacle.  Referee: fixed
σ′=(q c)σ(q c) (not σ′=σ), the (q c)∘ρ case proved directly (lem:adaptinv does not cover it), regress remark downgraded.
**Seg 4 — DATA on hyp:witness (no falsification) + one PROVED alternative sufficient condition (prop:firstpairbias).**
`s13_wtail.py` on fixed-rectangle samples (n=30: 7.2e5 instances; n=50: 3.6e5 each of (50),(25,25), 2.4e5 of (4,2^23); n=100:
8e4 each of (100),(4,2^48)): P[no witness in first K ladder pairs | K₀−1 ≥ K] = (1−w₁)^K within noise to K = 12/20/40; hazard
flat to 3 digits; P[Flip=∅|K] = 2^{−K} to K=12–15; across diagonals D_c (shared y-rows) P[G=0] and P[G<EG/2] binomial for all
ℓ ≤ 12, all thresholds; identical for extreme types.  Regime of hyp:witness not reachable.  **prop:firstpairbias** (found by
the referee of the recommendations draft; refereed; exact at n=5,6 all classes, n=7 sampled): P_X[B]−P_X[A] =
E[(1_B−1_A)(1_{Fᶜ}−1_F); Gᶜ] ⇒ |P_X[B]−½| ≤ ½P_X[Gᶜ] ⇒ **weak CGW ⇐ hyp:firstpair: sup_{λ,α,β} P_X[Gᶜ] = o(1/log n)**, an
ALTERNATIVE sufficient condition (rule 2: not proved weaker than hyp:witness; neither implies the other).  Cost first: still a
lower tail of a count of square-selected events; decomposes as P[μ ≤ M] (upper bound on a short arc at the frame pair, unproved)
+ P[M disjoint candidates all bad] (the obstacle).  Heuristic/data: per-candidate 0.41, Eμ = 6.3/9.6/18 at n=30/50/100, P_X[Gᶜ]
≈ 5.5/n (0.61, 0.18, 0.11 at n=7,30,50), P[Gᶜ|μ=k] ≈ 0.59^{k−1} (slightly above).  Notes §(k); `s13_firstpair_id.py`;
`NEXT_STEPS_2026-10-03.md` (refereed recommendations: R1 send gap note; R2 put hyp:firstpair + data in the paper, literature
search, post; R3 consultation with the unconditioned question first; R4 one experiment: P_X[Gᶜ] vs n and a switching bound for
P[μ ≤ M], with stopping rule; R5 not: more reformulations, more first-moment sampling, transfer tightening, completion-count side).
**Seg 5 — WRITING (R2a,b).** Paper now 25 pp: §6(e) "The first ladder pair alone" (lem:mu2, prop:firstpair, prop:firstpairbias,
hyp:firstpair, cor:firstpair, cost paragraph, rem:fixrect with the sampler, type-independence, witness-independence and P_X[Gᶜ]
data incl. n=100: 0.055); §7 rewritten (four hypotheses; eq:minimal reduced one-way to P_X[G] ≥ 2c; fixed-cell restatement;
literature sentence hedged, k=2 via KSS21 + exchangeability); Results item (6); item (2) states EQ(δ) correctly; abstract ("three
exact identities", data scope, 5.5/n labelled first-moment); status paragraph; date.  Fresh referee of the revised paper: no
mathematical error in §6(e); fixed: μ clash (arc minimum renamed ν in the paper), cor:firstpair input (uses |q̄−1| ≤ 2s/(1−s)),
rem:fixrect symmetry argument (isotopies fixing rows 1,2, not row exchangeability; "assuming connectivity and mixing"), stale
counts/figures, terminology (ρ = permutation of columns induced by rows; τ_t = δ_{q,c}), t₀−1 = 0 case, "cannot give (L)" scoped
to the first pair, §7 closing includes the ν-tail.  Paper sources: `paper/firstpair_block.py` (new), `cprime_blocks.py` (remains),
`assemble_conditional.py`, `postprocess_conditional.py`.  Did not touch the core obstacle.  Remaining before posting (R2c–e):
literature search (user), human read of §§2–4, the gap-note courtesy window.
**Seg 6 — R4 (one segment, stopping rule applied): DATA + one small PROVED lemma that is a one-way REDUCTION of the ν-term +
NEGATIVE KNOWLEDGE on the switching route.**  Data: n·P_X[Gᶜ] = 5.5, 5.5, 5.5, 5.4 at n = 30, 50, 100, 150 (type (n)); 5.5 for
(4,2^73) at n=150; P[ν=k] ≈ 4/n; P[Gᶜ|ν=k] ≈ 0.6^{k−1} at every n and type.  lem:nuswitch (notes §(l); refereed; verified on
2.5·10⁵ moves): for the staircase event E_k(a,b) (ρ-path from a ∈ {p,p′} first meets {p,p′} at step k at b), trading rows x₁,u
along the ρ_{x₁,u}-cycle W_u through w_{k−1} is always legal, keeps (x₁,y₁), and is rigidly invertible when W_u avoids
{p, w₀..w_{k−2}}; hence P_X[ν ≤ M, f ≥ θn] ≤ 4M/(θn), f = number of admissible rows u.  Where it stalls: the residual
P_X[ν ≤ M, f < θn] is a lower tail of a count of constant-probability events (≈1/(k+1) each) at frame-selected positions —
the same TYPE as the second term of hyp:firstpair (a different event).  Fixed-position route (exchangeability + union
bound over staircases) needs "P[Q ⊆ L | R] ≤ (C/n)^{|Q|} for Q in two rows", which always-legal row trades give only to one
factor 1/n per row; column turns are legal on a set whose size is again a lower tail.  Stopped per the rule.  Scripts
`s13_nutail.py`; data logs runs/s13/fixrect/fp*.log.
**What did NOT change (today):** no proven lower bound at a frame-selected position; both terms of hyp:firstpair are lower
tails of counts of constant-probability events at such positions; the only lower bound remains P_X[A] ≥ ¼ at distance 2.
hyp:witness remains the reference hypothesis; hyp:firstpair is an alternative.  Cold read done (section above).
**Seg 7 — literature (NEGATIVE KNOWLEDGE on one step of a candidate tool + WRITING).**  See LITERATURE CHECK: Allsop–Morris
η-switching legal and rigid in X with deterministic forward degree; their Thm 3.1 transfers except Claim 4(ii) in one configuration;
KPS: parity marginal of the weak conjecture proved by them; gap note §4 (KPS, AM citations) updated.  Did not touch the core
obstacle; located it inside an otherwise working argument (Claim 4(ii)), which is the sharpest consultation question we have.
**Seg 8 — WRITING (R2, cold-read fixes; refereed by a fresh subagent, 15 findings applied).**  Paper 29 pp.  (a) NEW small PROVED
item: rem:fewparts — thm:reduction uses EQ(δ) only along chains from (n) to plain λ with ≤ K = ⌈16 log n⌉ parts, so every
hypothesis is needed only for plain λ with < K parts and plain µ (referee verified both routing cases; K excludes Θ(n)-part types
asymptotically but NOT the sampled types at n ≤ 150, which all have < K parts; (4,2^13)→(2^15), (6,3^8)→(3^10) are excluded only
because α=β); A = 16 made explicit.  (b) Cold-read fixes: abstract (four conditions, none implying another; hazard caveat;
sampler caveat; "whose proof in CGW has a gap" insinuation removed; fragment fixed), Results (2),(6) (M = C log log n split with
the hazard condition stated with ν > t; "Θ(n) candidates" gone; "not shown to be easier"), §6(e) cost paragraph (hazard
reading; data-measurability; hidden costs: ν-tail, positive dependence, admissibility; "simplest" not "smallest"), lem:nuswitch
+ (FC) eq:FC moved into the paper with proof (referee: correct; C ≥ 2 added), rem:fixrect "[0.995,1.005]" as point estimates,
§7 "entirely" removed, min-hypothesis sup min(P_X[Flip=∅], P_X[Gᶜ]) recorded as the weakest condition short of EQ but not
promoted, exchangeability formula P_X[E] = (n−2)(n−3)P_X[E′].  (c) §7 literature paragraph (KSSS Lemma 3.3 incl. full rows
in Q, C = 1+O(k/n), rectangle→square transfer only at e^{O(n log² n)} (McKay–Wanless Prop 4, replacing the KS18 Prop 5
citation); DKKS Thm 1.6 with correct quantifiers (ε,η; γn-sparse); AM Thm 1.1 + η-switching legal in X with deterministic
forward degree; Thm 3.1: Claim 1 transfers, Claims 2–3 "appear to transfer, not written out", Claim 4(ii) configuration; AW
§3.1 quote; KPS parity marginal with rows/columns/symbols).  (d) KPS parity marginal in the intro.  (e) Bibliography sorted;
AM26, AW25, DKKS26, KSSS23, MW99 added.  Did not touch the core obstacle.  Cold read next (user's instruction).

## ALERT (2026-10-02): gap in CGW08 Lemma 3.12 / Theorem 3.13 (rem:cgwgap; paper/cgw_gap_note.tex, refereed)
Splitting case 3 of CGW's proof (cross-switch at the neighbour {ωj,ωj′}, then backflip
at {j,j′}) lands in the WRONG type when min(α,β) ≥ 3: the cross-switch reverses the
ω-arc from ωj to ωj′ (or the complementary arc, by the smallest-symbol rule), which contains j′ (resp. j), so {j,j′} ends at distance 2
(s13_cgwcase3.py: 1006/1006 at n=9, 3098/3098 at n=12 wrong type).  For min(α,β)=2<max the type is
right but the form reversing the LONG arc makes joining rule 4 cross-switch at a different pair, so
G_S ≠ G_J: exact at n=5, (5)→(2,3), all 161280 squares (s13_cgw_multigraph.py, runs/s13/cgw/):
172800 of 1036800 G_S edges absent from G_J (half of case-3 + their switch edges), µ-degrees in G_S
12/18/24 vs 18 asserted.  Only (2,2) is unaffected.  REPAIR for min(α,β)=2 (gap note Prop. 3, refereed;
s13_cgw_min2.py): attributing edges to joining pairs (direct at A-pairs, switch at B-pairs, case-2 switch edge to
the backflipped pair) gives ≤ 4 per pair, hence C(λ)/C(µ) ≤ 2 and **P_X[A] ≥ 1/4 unconditionally for marks at
distance 2 or m−2** (prop:alpha2) — the first gap-direction bound we have; fails for min ≥ 3 by type, not
multiplicity.  §(h),(i): two disjoint row pairs have numerically independent cycle types (n=30,50,100);
prop:genericlong: a generic row pair has a cycle > n/5 w.p. ≥ 1/2 unconditionally; KS/KSS inventory: transfer
lemma e^{O(n log² n)} ⇒ KSS intercalate concentration holds in every class of X; nothing in KS/KSS addresses
frame-selected families.  **§(j) ADAPTIVE COORDINATES (prop:adaptive, refereed, verified at n=30–100):**
take the coordinates of ladder pair k to be (ρ_k^t p, ρ_k^t p′), t ≤ min(T, μ_k−1): they are separated by {p,p′}
BY CONSTRUCTION and orbit-invariant (lem:adaptinv), so (β) disappears: P_X[Flip_I=∅] ≤ E[2^{−|U|}], U = pairs with a
short clean candidate cycle (and non-colliding).  New hypothesis eq:hypadaptive: E[2^{−|U|}; K₀ ≥ ℓ₀] = o(1/log n)
(+ offset version).  Cost: per-pair membership is rare (≈ Tℓ₁/n) — a witness-lemma-type lower tail — or, with
s ≈ T ≈ ℓ₁ ≈ √n, constant per pair but needs arc trapping at scale √n.  Consequence: CGW's 3/2 direction (⟺ P_X[B] ≤ 2/3 via thm:marked — a weak
form of OUR target) is unproved as published; the 1/2 direction is trivial given
thm:marked.  Conditional on "CGW(3/2)" from now on: eq:cgwA's P_X[A] ≥ 1/3,
prop:shortcycles lower bounds, the unconditional witness bound, unconditional (CP)
c = 2/5, eq:CPcols "≥ 1/3", CGW's P[(n)] ≤ 2n^{−2/3} and Cor 4.5, hence the TAIL BOUND in
thm:reduction — DONE: prop:tailsurvive (refereed) gives E[N_ℓ] ≤ 4/ℓ (ℓ ≤ n/2),
E[(κ)_k] ≤ 2^{k−1}(4+k²/log n)(log n)^k and P[κ ≥ 16 log n] = O(n^{−2} log n) from |J| = 2|X_A|
alone, so thm:reduction is now independent of CGW's unproved direction.  Also: CGW(3/2) ⟸
P_X[Flip=∅] ≤ 1/3 (thm:ladder; sharp at n=5), i.e. our (L) in constant form would repair their theorem.
Everything exact of ours (thm:marked, thm:ladder, offsets, lem:rowfair, lem:legalturn,
prop:orbit, eq:Jp) is unaffected.  No published erratum found (quick search only).

## Goal and chain — REWRITTEN after the ladder identity (§sec:ladder)

CGW weak conjecture d_TV(P_n,Q_n) → 0  ⇐  EQ(o(1/log n)) [thm:reduction]
⇐ |P_X(B) − ½| = o(1/log n) [thm:marked]
⇐ **(L): P_X[no ladder pair flippable] = o(1/log n)** [thm:ladder, cor:ladder].

**thm:ladder (PROVED, exact, verified to all digits at n=7):** with
K₀ = min(|C₁|,|C₂|) (apart) / min(d₁₂,d₂₁) (together) and the ladder pairs
(x_k,y_k) = (π⁻ᵏ(1), π⁻ᵏ(2)), 1 ≤ k < K₀, the j′-trades T_k at the ladder
pairs are commuting involutions of X preserving the ladder; each toggles the
bit iff its pair is flippable (j, j′ in different ρ_{x_k,y_k}-cycles).  Hence
on every orbit of ⟨T_k⟩ the bit is EXACTLY fair unless no ladder pair is
flippable, and
    P[B] − P[A] = P[B, Flip=∅] − P[A, Flip=∅],   |P[B] − ½| ≤ ½ P[Flip=∅],
for X and for every fixed R (quenched).  Data: P[Flip=∅ | K₀] = 2^{−(K₀−1)}
to the resolution of 9000/7200 instances at n=30/50; P[Flip=∅] ≈ 4/n (0.125, 0.079, 0.041 at n=30,50,100).

**(L) is the one open item, and it is now reduced to the witness lemma
(hyp:witness, prop:Lwitness).**  Offset ladders D_c = {(π^{−(c+i)}(1), π^{−i}(2))}
(thm:ladderoffset, exact at n=7 for c=0..3) price the short-arc B-states:
averaging the offset identities over c gives P[B, d₁₂=ℓ, |C₁| ≥ n/2+ℓ, θ ≥ θ₀]
≤ 36/((n−2)θ₀) where θ = fraction of offsets whose diagonal has a flippable
pair (eq:offsetavg); the splice-in switching does |C₁| < n/2+ℓ.  So:
- (L) ⇐ witness lemma: for the ladder (≥ n/log²n disjoint pairs) and for the
  diagonal blocks of a short arc, the number of blocks containing a pair whose
  ρ-cycle through j has length ≤ ℓ′ = log⁵n and avoids j′ is ≥ κ·Σ min(½, ℓ′|B|/n)
  except w.p. o(1/log n) (annealed).  Budget: ℓ₀ = n/log²n; total o(1/log² n).
- Data (witness50.log): witness freq ≈ (ℓ′−1)/(n−1) − ℓ′/n; disjoint pairs
  E[WW′]/p² = 1.00±0.03; Var(G) and Var(star count) = binomial within 5 %.
- **CGW08 Thm 3.13 IS P_X[A] ≥ 1/3** (via thm:marked), and summed over types it
  gives the UNCONDITIONAL witness lower bound with no loss in ℓ:
  (2/3)E[(ℓ+β)N_{ℓ+β}] ≤ E[ℓN_ℓ·βN_β] ≤ 2E[(ℓ+β)N_{ℓ+β}], hence E[N_ℓ] ∈ [2/3, 4]/ℓ
  and P[|Q_q| = ℓ] ∈ [0.66, 4]/n for ℓ ≤ √n/40 (prop:shortcycles, refereed;
  exact at n ≤ 11 from CGW's tables; sharp at n=5). Remaining: the CONDITIONAL
  form.  **The 'path version' of CGW Lemma 3.12 is NOT routine (§sec:crossspace):**
  with columns p,p′ prescribed, rows x,y of a ladder pair have two PATHS in the
  free part, and the cross-path pairs (the only ones whose flip exchanges
  crossed/parallel) cannot be repaired when they are B-pairs (switch would trade
  through p; cross-switch arc passes through p′).  CGW's condition (3) is exactly
  what excludes paths — essential, not cosmetic.  What survives in the cross
  space X (rows 1,2 + columns p,p′): (F) flip identity at legal A cross-path
  pairs; (I) INTERCALATE identity E[I_sep; crossed] = E[I_merge; parallel]
  (intercalates on rows x,z, columns q∈P₁, c∈P₂; always legal, frame-preserving);
  data: I = 0.98·|A||B|/n in every bin (icsw50.log).  (I) is FAIR and gives no
  inequality by itself: the deterministic bound is only I_merge ≤ (ℓ_p−1)(ℓ_p′−1)
  (referee: Z₈×Z₂ counterexample to the earlier min(ℓ_p,ℓ_p′)−1 claim); a constant
  needs two-sided intercalate control in the cross space (i)+(ii)+(iii).
  **CORE PROBLEM (CP), eq:CP:** given rows 1,2 (type+mark) and columns p,p′, for
  rows x,y with four distinct symbols on {p,p′}: P[p,p′ in different
  ρ_{x,y}-cycles] ≥ c.  Unconditional (CP): c = 2/5 (1/3 ≤ P[W] ≤ 3/5+4/(5(n−1)),
  refereed; exact at n=5).  **Columns-only (CP) is CGW transposed**
  (eq:CPcols: = 1/2 apart, ≥ 1/3 together) — the obstacle is the conditioning
  on rows 1,2 (type+mark = X itself), NOT the frame.  thm:ladder = the legal
  half of transposed CGW; the illegal half (switch/cross-switch at {x_k,y_k} in
  Lᵀ) = flip / cross-switch of rows 1,2 at (p,p′): toggles ALL side-1×side-2
  statuses but reverses one side (ladder ↦ anti-diagonal).
- **Orbit method (§sec:orbit, new; PROVED parts: lem:rowfair, lem:legalturn,
  prop:orbit, prop:Lorbit):** in X every row trade of rows ∉{1,2} is legal
  (⇒ exact fairness of A/B status of (q,c) for generic rows x,y when q,c lie in
  different ρ_{x,y}-cycles: P_X = 1/2, lem:rowfair); a column-cycle turn is legal
  iff the cycle meets {1,2} in 0 or 2 rows; turns at a fixed matching M of free
  columns commute.  "Clean coordinates" of a ladder pair (x,y): short cycles
  through x avoiding y, the family, and meeting {1,2} in 0/2.  They generate a
  free (Z/2)-action; conditionally on the orbit the bits of a family of disjoint
  pairs are INDEPENDENT and P[(x,y) crossed | orbit] ≤ 1 − ½ȳ_{x,y}, ȳ = orbit
  average of the fraction of clean coordinates separated by {p,p′}.  Hence
  P[Flip_G = ∅] ≤ E[exp(−½ Σ ȳ)] (prop:orbit), and (L) ⇐ hyp:orbit
  (Σ_{sub-ladder} ȳ ≥ κ|I| w.p. 1−o(1/log n), + block form for offsets; the
  flippability concentration of hyp:witness is no longer needed — Chernoff given
  the orbit).  Data (orbit{50,100}[_fixed].log): column-cycle lengths through x
  uniform; P[separated | clean, ℓ] = 0.19–0.20 = P[separated] (genericity);
  orbit coin 0.478/0.523 at n=100 (greedy matching, |T|≈42); fixed matching:
  |T| ≈ 3, E[ȳ] = 0.19 ⇒ per-pair factor ≈ 0.9, κ ≈ 0.15.
  Residual = genericity (β) + macroscopic arcs (α) + concentration (γ).
  (Π) SOLVED (§sec:orbit (f)): P_η[b=1] ≥ ½·max_Z P[Z crossing]; b ≡ 0 on the cube iff
  P₁,P₂ are in different components of the block graph (vertices = arcs P₁,P₂ and the
  other cycles, edges = chords of T); one coordinate per pair suffices (eq:onecoord):
  P_X[Flip_I=∅] ≤ E[Π(1−½X_k)], X_k = [first clean coordinate of pair k crossing].
  (γ): data (cov{50,100}.log) — orbit weights of different pairs UNCORRELATED
  (Var(Y)/(|I|Var ȳ) = 0.99/1.14, excess explained by Var|I|).
  (α) ⇐ (BL) [interleaving of {1,2},{x,y} in column cycles ≤ 1−c; data 0.19].
  (U) in its simplest instances (§(g), refereed): generic columns ⇒ (U) IS the two-arc
  form of (α) (both arcs/cycles at p,p′ macroscopic); generic rows ⇒ (BL) is TRIVIAL
  (interleaved fraction over row pairs ≤ ½ + O(1/n) deterministically; exchangeability)
  but is lost under the conditioning |Q_p(ρ_{x,y})| = ℓ / frame selection, which is the
  content.  DIRECTION: lower bounds on different-cycles events not conditioned on their
  own toggle ((β), third-row availability, CGW's P_X[A] ≥ 1/3) are gap-direction; upper
  bounds are surviving-direction.  Third-row repairs exist (crossed ladder pair: cross-
  switch of x_k,y_k at (p,p′) always defined, lands in J_B with (d₁₂,d₂₁) → (2k,|C|−2k);
  A-instance: backflip through C₂ → J_A) but have unbounded multiplicity / don't preserve
  the ladder — the precise reason switchings in X cannot give the gap direction.
  ONE FORM: every input is "two pairs of lines interleaved in a third pair" (U):
  c ≤ P[interleaved] ≤ 1−c in X, plus the same-cycle event for generic rows at the
  prescribed columns (availability of the (BL) repairs).  Unconditionally (no type
  conditioning) (U) is EXACT given the cycle type: P = A/3 + 2B ∈ [1/27, 1/3+O(1/n)]
  (conjugation invariance; referee).  All exact identities in X are fair; prop:orbit's
  gain is that fairness compounds over a family.
- Literature (scout, §sec:witnessstatus): the conditional form is not in print;
  Allsop–Morris 2026 give (δ/n)^{|P|} ≤ P[L ⊇ P] ≤ (Δ/n)^{|P|} (exp loss in ℓ;
  δ=Δ=1+o(1) open — that IS the lower-bound clause); KS18's exact ratio method
  (1+O(1/n))/(s+1) for intercalates in two rows is the model to follow.
  Needed: two-row cycle switchings with degrees exact to 1+O(polylog/n), in the
  space of completions of rows 1,2 and columns j,j′ (never touched by the moves).

Superseded by (L) (keep as results, do not work on): (I), (II)/(T), (III),
(F1), (F2), the component question, the pool/master-formula route, the
spectral stationary model (still a standalone theorem: paper/spectral_same_cycle.tex).

## Older map (pre-ladder), kept for orientation

### (T), i.e. (II) the transfer — what it contained

(T) is one statement but it carries three things: (a) the real frame law
π_P(L) under Unif(X) must be close to uniform ON THE TEST FUNCTION
ε·P^Mε (a function of the cycle type and the cycles of u,v) — the quenched
one-point content of the old item (I) is relocated here, not removed;
(b) the pool must be in generic position relative to the frame; (c) R/M
row-sharing rare. Status:
- **Trapping half of (a) PROVED** (prop:trapped, refereed): P_X[row 1 in a
  frame cycle of length ℓ] ≤ 18/(n−ℓ−1), by a Kwan–Sudakov switching whose
  only row-1 move is the CGW switch at the mark (priced by thm:marked).
  Trapped marks cost O(1/k'). The rest of (a) is OPEN.
- (b) OPEN. Exact structure (§s13transfer2): delete columns j,j′; frame =
  leftover graph of the n×(n−2) rectangle; pool = ¼-coin thinning of a
  candidate graph; only local exclusion = π-adjacent rows. Real squares
  (true marks) decorrelate 10–30% FASTER than the model at n=50, half that
  at n=100 (the adjacency exclusion), frame statistics agree.
- **Trade-chain route** (§s13tradechain, §s13comp7; standalone statement
  paper/fairness_question.tex, refereed): symmetric chain on S(R) (trades along
  the ρ_{x,y}-cycle through j′, turns of free frame cycles); EQ(4δ) ⇐ bit
  autocorrelation ≤ δ² ⇐ component fairness Σ|K|(2p(K)−1)²/|S(R)| → 0.
  **Exact enumeration at n=7 (s13_comp7.c, runs/s13/comp7/, 7 instances):
  one giant component (99.1–99.7% of S(R)) with p(K) = P[β=1|R] to 4 decimals;
  the rest is "locked" dust (ℤ₇-like squares, every ρ_{x,y} a 7-cycle, nothing
  flippable; components = 5! row permutations), 0.3–0.8%, carrying the whole
  LHS (0.003–0.008). n=6 is degenerate (row parity invariant).**  So the
  component structure is trivial and the Question for fixed R IS the quenched
  statement P[β=1|R] = ½+o(1); the chain is the tool, not a weakening.
  Rectangle picture (delete columns j,j′): frame cycles = leftover row–symbol
  cycles, turns = orientations, bit = function of the rectangle; flippable ⇔
  the two (x,y)-paths pair "parallel"; for x,y in different cycles a turn
  toggles it, so thinning only bites inside C₁∪C₂ (the bit-changing pairs).
  Within-cycle mark ⇔ (1,2) crossed; cross-cycle marks have an exact
  bit-flipping involution (switch at (1,2) + isotopy) ⇒ p = ½ exactly.
  Needs (F1): P_stationary[no flippable bit-changing pair] = o(1/log² n) —
  a parallel/crossed statement about row-pair paths in a random rectangle,
  target for a KS switching as in prop:trapped — and (F2) mixing of the
  thinned walk vs the "coin walk" on D = {derangements, π(1)≠2, π(2)≠1}
  (uniform on D stationary, P_D[1~2] = ½ exactly; fibres {π⁻¹(1)=a,π⁻¹(2)=b}
  = permutations of n−2 points with no fixed point outside {a,b}, exact ½
  there via π ↦ (a b)∘π).  Right-multiplication trades (through j) are
  equally legal and were included as move set B (no qualitative change).

## Bugs found this session
- Sessions 7–8 real-square pipelines (s7/s8_real_bias.py) used σ's SYMBOLS as column
  indices on unnormalised JM squares: their "marks" were random column pairs.
  Qualitative findings survive; §sec:realbias numbers are for random pairs.
  Fixed in s13_real.py / s13_trade.py (symbol → column via row 1).

## Closed / superseded (do not spend time on these)

- (III) prob:decouple: NOT NEEDED (prop:nodecouple; routing merges m ≥ n/(A log n)).
- (I) in all its session 8–12 forms (prob:mixing one-block start, prob:survival,
  (R1), (C1), (C2″), (S₂), occupation-time lemma, Markov renewal, a_k, t_k(w)):
  closed IN THE MODEL; the quenched frame-law question survives inside (T)(a).  The exact laws (thm:blockcount,
  thm:sizebiased, thm:idleloop, lem:inert, prop:nomino) remain valid results.
- Budget: γ > 1/11 not 4/9 (cor:budget) — now irrelevant since c = 1/2 is proved.
- Dust obstruction: a sup-norm statement; irrelevant to L²/stationary (rem:dustblind).
- Still true and worth keeping: lem:flipsrecur (τ_k < ∞ a.s.); the audit of §s13.

## Do not reopen (proved impossible / false)

Doeblin minorisation (prop:nomino); (V-)uniform ergodicity (cor:noerg);
λ-drift (prop:nodrift); geometric (C2); (R2); transport contraction;
two-block odd perturbation; ρ from burnt-in starts (pitfall 56).

## Possible next steps (owner's choice)

-1. (CP) via hyp:orbit (the current best form of the residual).  Done: (a) the
   one-step switching is fair (useless); the orbit method compounds fairness.
   Open inputs (all of the form (U), see above): (α) arc-trapping in X — DONE
   up to (BL): CGW's lemma for the row pair (x,y) runs inside X with the flips
   blocked exactly when rows 1,2 and rows x,y are INTERLEAVED in the (q,c)-column
   cycles (switch-invariant); eq:Jp gives P_X[|Q_p(ρ_{x,y})| = ℓ] ≤ 3/((1−b_ℓ)(n−ℓ))
   with b_ℓ = interleaved fraction (data 0.19, interleave50.log); (BL): b_ℓ ≤ 1−c is
   again a same-cycle genericity for third rows (the recursion, now explicit);
   for ladder pairs the frame must be fixed too (bookkeeping not done); (β) genericity: P[clean coordinate (q,c) separated | ...] ≥ c —
   note the free row trades (lem:rowfair) re-randomise the clean coordinates
   without touching ρ_{x,y}, but their availability ("q,c in different
   ρ_{z,w}-cycles") is again a same-cycle event; (γ) concentration of Σȳ over
   the sub-ladder (second moment: average covariance o(1/log n) suffices).
   (b) relax rows 1,2 to type + mark: DONE in the sense that the orbit method
   works in X (type+mark), with column turns through both rows 1,2 legal.
   (c) exact identity not pinned by a row pair and a column pair: the orbit
   coordinates are exactly such (pinned by row x and column pair (q,c) only).
0. (superseded) (F1)/(F2).

1. (II): formulate the minimal "coarse product" property of the pool and
   the frame-health property, and test both on real squares at n = 100–400
   against the exact prediction; then the TRP/KPS robustness argument.
2. Publishable checkpoint: thm:marked + thm:master + thm:spectral +
   cor:stationarymodel + the layer theorems is a paper now; thm:spectral
   is a standalone result about the random-transposition walk (check
   literature: two-point/same-cycle statistics of random transpositions).
3. Literature check whether thm:spectral is known.

## Numerics this session

`runs/s13/frame*.log` (frame comparison), `akB_0_p.log` (p = 2.4, moot),
`akB_{6,7}_k40.log` (moot); `comp7/*.log` (exact components n=7, ~1 min each);
`trade{30,50,100}.log` (trade chain on real squares); `f1/ladder{30,50,100}.log` (ladder
statistics), `f1/witness50.log` (witness lemma scoping), `f1/f1_*.log` (flippability of bit-changing pairs: fair, independent), `f1/icsw50.log` (intercalate identity), `orbit/orbit{30,50,100}.log`, `orbit/orbit{50,100}_fixed.log`, `orbit/interleave50.log`, `orbit/cov{50,100}.log` (orbit method / blocked fraction / covariance of orbit weights: `./jm 100 120 300000 9 | python3 s13_orbit.py 100 --inst 3 --lmax 8 --fixed`).  `s13_spectral.py brute N` (N ≤ 7) and
`s13_spectral_big.py N` reproduce the theorem.

## Quickstart

HANDOFF_session_12.md §3 QUICKSTART still green, plus:
`python3 s13_spectral.py brute 6` (exact = brute, all digits);
`python3 s13_spectral_big.py 4000` (sum ν = 1, m·E → 1.9);
`gcc -O3 -o s13_frame s13_frame.c -lm && ./s13_frame 2 32 100000` → 4E[bias²] ≈ 0.115;
`pdflatex second_row_notes.tex` ×2 → 0 errors, 123 pp;
`./jm 9 300 2500 3 | python3 s13_cgwcase3.py 9` → CGW case 3 lands in the wrong type (every instance);
`gcc -O2 -DN=7 -o s13_comp7 s13_comp7.c && ./s13_comp7 0123456 1234560 0 2` → ladder identity exact (~1 min);
`./jm 30 300 30000 1 | python3 s13_ladder.py 30` → P[Flip=∅|K₀] ≈ 2^{−(K₀−1)}.
