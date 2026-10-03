# Recommendations on next steps (2026-10-03; draft refereed by an independent subagent, revised)

## Where the project stands (facts)

**Proved and refereed (by subagents, not humans).**  thm:marked; prop:tailsurvive; thm:reduction (weak CGW ⇐
EQ(o(1/log n)), independent of CGW's unproved direction); thm:ladder and offset ladders; prop:trapped; prop:splice;
orbit method (lem:rowfair, lem:legalturn, prop:orbit, lem:adaptinv, prop:adaptive); prop:genericlong; prop:alpha2
(P_X[A] ≥ 1/4 at distance 2); prop:firstpair (single-pair involution); **prop:firstpairbias (new today, found by the
referee of the first draft of this document, verified exactly at n = 5, 6):**

    P_X[B] − P_X[A] = E_X[(1_B − 1_A)(1_{Fᶜ} − 1_F); Gᶜ],   hence   |P_X[B] − ½| ≤ ½ P_X[Gᶜ],

F = {(x₁,y₁) flippable}, G = {some candidate t < μ at the first ladder pair is good}.  So weak CGW follows from
**hyp:firstpair: sup_{λ,α,β} P_X[Gᶜ] = o(1/log n)** — an alternative sufficient condition (not known to be weaker
than hyp:witness; neither implies the other), concerning one row pair, with constant per-candidate probability
(≈0.41 in the data) and Θ(n) candidates (heuristically; Eμ = 6.3, 9.6, 18 at n = 30, 50, 100).  Data: P_X[Gᶜ] ≈
5.5/n (0.61, 0.18, 0.11 at n = 7, 30, 50).  Gap note (CGW08 Lemma 3.12 / Thm 3.13) with constant-2 repair.

**Unproved.**  Any lower bound on the probability of any event at a row or column pair selected by the frame,
inside X.  Both sufficient conditions (hyp:witness, hyp:firstpair) are lower tails of such counts.

**Data (fixed-rectangle sampler, n ≤ 100, including types of P_n-probability e^{−Θ(n)}).**  P_X[B] = ½ ± 0.002;
P_X[flip₁] = ½ ± 0.002; ladder bits are independent fair coins to K = 15–40; witness events are independent
Bernoulli(w₁) to three digits along the ladder (flat hazard) and binomial across the diagonals; identical across
types.  Not reached: the regime of either hypothesis.  Unproved: connectivity of the restricted chain (tested
exactly at n = 5, 6, and against the exact P[B|R] at n = 7).

## Reading of the situation (revised after the referee)

The obstacle is the absence of a technique for lower bounds on square-selected events under two-row
conditioning; the literature's tools (CGW switchings: gap direction only with a rigid inverse, which exists at
min(α,β) = 2 only; Kwan–Sudakov switchings and the KSS process approximation: fixed positions or global counts,
transfer to X at e^{O(n log² n)}) stop there.  The referee's correction to the first draft: the return on the
reformulations is not zero — prop:alpha2 is a real unconditional bound, and prop:firstpairbias replaces a
hypothesis about ≥ n/log²n rare events (plus offset families and a splice bound) by one about Θ(n)
constant-probability events at one pair.  That is a simplification of the *sufficient condition*, not a reduction
of the obstacle, which is unchanged in kind.  The referee also notes that even WITHOUT the two-row conditioning,
a lower bound P[the ρ_{3,4}-cycle through p has length k] ≳ 1/n for 3 ≤ k ≤ log⁵ n is not in the literature as far
as either of us knows (only intercalates, k = 2, are settled); so "the only obstacle is the conditioning" was an
overstatement and has been removed.

## Recommendations, in order

**R1. Send the gap note now** (after the author's own read; `paper/cgw_gap_note.tex`, [author] placeholders).
Self-contained, computer-verified, states what survives, includes the constant-2 repair.  Give the authors ~3
weeks before any preprint that cites the gap.  Cost: an afternoon.  Risk: none.

**R2. Put hyp:firstpair into the paper and the data into the paper, then post.**  Concretely, in this order:
(a) one writing segment: prop:firstpairbias in §4 or §5 (its proof is five lines given thm:ladder's T₁ and
prop:firstpair), eq:hypfirstpair as an alternative sufficient condition next to hyp:witness (rule 2: neither is
proved weaker, so neither replaces the other; state the comparison as heuristics), the §(k) data in one paragraph
(sampler, type-independence, witness independence, P_X[Gᶜ] ≈ 5.5/n), connectivity of the sampler flagged as
unproved; (b) correct stale statements in the abstract/introduction ("typical types only", "P_X[B] = 0.49–0.50
(±0.01)"); (c) a real literature search (arXiv since 2018: two-row cycle structure, errata to CGW08,
Kwan–Sah–Sawhney–Simkin "Substructures in Latin squares"), since STATUS records only a quick one; (d) one human read
of §§2–4; (e) post.  Cost: 2–3 days.  Reason: the exact identities, the CGW-independent reduction, and two cleanly
stated sufficient conditions are the durable output; a public statement of the obstacle is the most likely route
to the missing idea.

**R3. One consultation, combined with R1 where natural** (Wanless and coauthors are natural consultees; otherwise
someone working with the KSS machinery).  Ask the UNCONDITIONED question first, as the referee advises: "for a
uniform Latin square and fixed columns p, p′, is there a lower bound P[the ρ_{3,4}-cycle through p has length k and
avoids p′] ≳ c/n for 3 ≤ k ≤ log⁵ n?" — then the conditioned form (completions of a fixed 2×n rectangle, two cells
(3,p) = L[1][p′], (4,p) = L[2][p′] prescribed; by row-exchangeability this is exactly the frame-selected event).
Cost: one e-mail.  The answer "open even unconditioned" would itself be worth knowing.

**R4. One more targeted in-house experiment, with a stopping rule (replaces the first draft's R4).**
Measure P_X[Gᶜ] with all candidates as a function of n (30, 50, 100, 150) and type, and P_X[Gᶜ | μ = k] against
0.59^{k−1}; and attempt the one piece that is an UPPER bound, P_X[μ ≤ M] ≤ CM/n at the frame-selected pair, by a
splice-type switching (prop:splice's method; note it bounds a different frame event, so this is not a corollary).
Stopping rule: one segment; if P_X[Gᶜ] does not decay like 1/n, or the μ-tail does not yield to a switching, stop
and record.  Expected value: the μ-tail is the only part of either hypothesis in the direction switchings can
handle; proving it would isolate the obstacle to "among M disjoint candidates, not all bad", its smallest
sufficient form.  (Exact enumeration at n = 8 is out of reach: ~1.8·10¹¹ completions per rectangle.)

**R5. Not recommended.**  (a) Further reformulations of the hypothesis by new coordinates or families: each has
been exact and each has landed on the same obstacle; (b) more sampling of first- or second-moment quantities:
done; (c) tightening the transfer lemma (e^{O(n log² n)} → 2^{o(n/log² n)}): conjecture-strength; (d) the
completion-count side as a route to 1 + o(1/log n) (van der Waerden/Brégman give e^{O(n log² n)}; Godsil–McKay-type
asymptotics cover k = o(n^{6/7}) rows, not completions; log L(n) is known to a factor 1 + o(1) only) — usable as a
heuristic for the O(1/n) term at most.

**R6. Process.**  This document counts as the cold read due after five segments (four today: writing, data,
reduction, data+bias).  Next cold read after R2(a).  The working rules held for the day; the one lapse was in
my first draft of this document ("zero return", "the only obstacle is the conditioning"), caught by the referee.
