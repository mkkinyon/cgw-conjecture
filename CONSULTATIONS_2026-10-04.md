# Suggested consultations (as of 2026-10-04)

Five contacts, in a recommended order, each with: why them, what to ask, what to attach, what not to claim, and a draft in your voice. Everything in square brackets is a placeholder or a choice for you. The drafts are deliberately short; the attachments carry the detail. Dates and page counts refer to the snapshots in `paper/versions/`.

Attachments referred to below:

| Document | File | Pages |
|---|---|---|
| Gap note | `paper/versions/cgw_gap_note_2026-10-04a.pdf` | 10 |
| Conditional paper | `paper/versions/weak_cgw_conditional_2026-10-04b.pdf` | 38 |
| Loops note | `paper/versions/loops_note_2026-10-03a.pdf` | 6 |

Recommended order: (1) Wanless/Cavenagh/Greenhill, then (2) Allsop/Morris — the second can go out shortly after the first, since the paper cites the gap note and it is better that the gap note has reached its authors before third parties read about it; (3) Eberhard and then (4) Cameron are independent of the first two and can go any time; (5) is optional.

---

## 1. Wanless, Cavenagh, Greenhill — the gap in Lemma 3.12

**Why them.** Authors of CGW08. The gap is in their proof; they are the people who can confirm it, repair it, or point to a repair, and the ones who should decide about an erratum. Wanless is also the natural reader for the downstream effects on CW16 and GMW25 (he is an author of both).

**What to ask.** (i) Whether they agree that the proof of the upper bound in Lemma 3.12 fails for every split other than (2,2), as the note describes (§3, "case 3 lands in the wrong type"). (ii) Whether they know a repair, or whether the statement should be regarded as open. (iii) Whether they would like to coordinate on an erratum, and how they want the downstream dependence (CW16 Theorem 1.1 at two points; KPS Theorem 1.3(3)–(4) through it) described.

**Attach.** The gap note only. Do not attach the conditional paper in this first message; it can follow if they ask what the note was for.

**Do not claim.** That the statement of Lemma 3.12 is false (the note says it is "expected to be true"); that CW16's theorem is false (it is true, by KPS 1.3(4)); that anything about the CGW conjecture itself has changed. Do not mention the loops note or the conditional paper's hypotheses.

**Draft.**

> Subject: A possible gap in the proof of Lemma 3.12 of "The cycle structure of two rows in a random Latin square"
>
> Dear Ian, Nick and Catherine,
>
> While working on the second-row conjecture from your 2008 paper I had to go through the proof of Lemma 3.12 in detail, and I believe the argument for the upper bound (the direction used in Theorem 3.13) does not go through for splits other than (2,2): the cross-switch of case 3 reverses an arc of the cycle rather than exchanging the two columns, and the resulting square lies in a different type class from the one the inverse count assumes. I have written this up in the attached short note, together with what I could verify survives (the lower bound of Lemma 3.12 and everything that follows from it alone, both assertions of Theorem 4.9, and the (2,2) case, where the upper bound holds with constant 2), and a section on the papers that cite the affected results — in particular Theorem 1.1 of Cavenagh–Wanless 2016, which as far as I can tell depends on the gap at two points, though its statement is now known to be true by other means (Kwan–Petrova–Sawhney, Theorem 1.3(4)); the use in Gill–Mammoliti–Wanless 2025 is repaired in the note.
>
> I would be grateful if you could tell me whether I have misread the argument, and if not, whether you know a repair. I have not been able to find one beyond the (2,2) case, and the note says so. If you think an erratum is appropriate I would of course prefer to coordinate with you on how the dependence of the later papers is described.
>
> I should say plainly that I have no doubt about the statement itself, which the KPS results strongly support; the question is only about the proof.
>
> With best regards,
> [name]

---

## 2. Allsop and Morris (cc Wanless) — the two-complete-rows extension and the lower-bound question

**Why them.** §7 of the conditional paper runs their switching proof of Theorem 3.1 / Lemma 4.1 with a partial Latin square containing two complete rows; they are the only people who can check that at a glance, and §8 asks a question that is phrased entirely in the vocabulary of their Claims 2–4. Wanless in cc because he is Allsop's coauthor on the related rectangle paper and because the question is about the second-row conjecture.

**What to ask.** Three things, separately. (i) Whether §7 is right: in particular that their Claim 2 is a statement in aggregate over the class A_k and not per arc (a referee of ours found a per-arc reading false, and the paper's Claim 2′ is stated in aggregate for that reason), and that keeping every column-cycle switch and cross-switch away from two complete rows costs only constants. (ii) The question: in a uniform completion of a fixed 2×n Latin rectangle, with two further rows x, y pinned by one anchor cell each in a common column p (L(x,p) = L(1,p′), L(y,p) = L(2,p′)) and a column pair (q,c) read off row y (the columns of y holding L(1,p′) and L(x,p′)), is there a switching argument giving P[x, y on different (q,c)-cycles] ≥ c > 0? State what you have: with (q,c) and y fixed, their Claims 2–3 plus row exchangeability give it with 1/30 (Proposition 8.5); with the pair read off y, the length of the cycle through x and y is Θ(n) with the right probability (Proposition 8.3, by an insert/cut switching pinned by x and y), but a positive density of good pairs across the cut between x and y is what their cross-switch count does not deliver, because the cross-switch toggles the interleaving status of the straddling pairs together with their goodness (§8(d)). (iii) Whether they see a move that breaks what §8(d) calls polarisation.

**Attach.** The conditional paper. Point them to §7 and §8; say the rest is context.

**Do not claim.** That the conjecture follows from an answer to (ii): say precisely what would follow (the first step of the hazard hypothesis; the full hypothesis needs all O(log log n) candidates, with the earlier ones preserved, which is a further problem). Do not describe §7 as a result "about" their paper in a way that could read as a correction; it is an extension that uses their proof unchanged except for the restriction of the moves.

**Draft.**

> Subject: Your Theorem 3.1 with two complete rows, and a question in the language of your Claims 2–4
>
> Dear Jack and Patrick,
>
> I have been using your universal probability bounds in work on the Cavenagh–Greenhill–Wanless second-row conjecture, and I would like to ask you two things about the attached draft, in which §7 and §8 are the parts that concern you; the rest is context.
>
> §7 runs your proof of Theorem 3.1 and Lemma 4.1 for a partial Latin square that contains two complete rows, with every column-cycle switch and cross-switch kept away from those rows (the "arcs" of Claims 2 and 3 are cut by the special rows, Claim 2 is applied in aggregate over arcs of equal size — a referee of ours pointed out that a per-arc reading of your Claim 2 is false, and I believe your statement is indeed the aggregate one — and the admissible rows of Claim 4 are those on cycles containing no special row). The outcome is your per-cell bound, with a worse constant, for completions of a fixed 2×n rectangle. I would be very grateful if you could tell me whether you see a problem with this; you are the natural judges and I would rather hear it from you than from a referee later.
>
> The question is this. In a uniform completion of a fixed 2×n rectangle, take two further rows x, y pinned by one cell each in a common column p (L(x,p) = L(1,p′), L(y,p) = L(2,p′)), and let (q,c) be the two columns of row y that hold L(1,p′) and L(x,p′). Is there a switching argument that gives P[x and y lie on different (q,c)-column cycles] ≥ c for an absolute c > 0? What I can do is in §8: with (q,c) and y fixed in advance, your Claims 2–3 plus the exchangeability of y give it with c = 1/30 (Proposition 8.5); with (q,c) read off row y, the cycle through x and y has length Θ(n) with probability 1 − O(δ) at scale δn (Proposition 8.3, by inserting inactive cycles into it and cutting them out, with the inverse pinned by x and y), but the step that fails is a positive density of "good" pairs — row cycle through c avoiding q — across the cut between x and y on that cycle. Your cross-switch count does not give it, because the cross-switch at an interleaving pair toggles the interleaving status of the straddling pairs along with their goodness, and §8(d) describes a configuration that all the moves I know preserve. If you see a move that breaks it, or a reason the question is hopeless by these means, I would like to know either way.
>
> I should be clear about what would follow from a positive answer: only the first step of the hypothesis in §7 (a constant conditional hazard at the first of O(log log n) candidates); the later steps need the earlier candidates preserved, which is a separate problem. So this is a question about your method rather than a request to prove a conjecture.
>
> With best regards,
> [name]
>
> [cc: Ian Wanless]

---

## 3. Eberhard — the loops note

**Why him.** The 2016 sketch (maximal-subgroup remark) is his; the note proves the statement he sketched and credits him. He should decide between coauthorship and attribution before anyone else sees it.

**What to ask.** Whether he would like to be a coauthor, or prefers an acknowledgement and a citation of his sketch; whether the account of his argument in the note is accurate; whether he knows of any prior proof.

**Attach.** The loops note.

**Do not claim.** That the result is "folklore" without qualification (you have decided to leave that sentence for later; it should be settled before sending to Cameron, see 4).

**Draft.**

> Subject: Almost all loops have multiplication group S_n — a short note, and a question about authorship
>
> Dear Sean,
>
> In loop theory it has long been said in passing, and never (as far as I know) proved, that the multiplication group of a random loop is almost always the full symmetric group. Cavenagh, Greenhill and Wanless noted in 2008 that it would follow from their conjecture on the cycle type of two rows of a random Latin square; the attached six-page note proves it unconditionally — the rows of a uniform reduced Latin square generate a group containing A_n with probability 1 − e^{−(log 2 − o(1))n²}, and S_n with probability 1 − 3cⁿ — by the route you sketched in 2016: the maximal-subgroup reduction, Praeger–Saxl and Pyber for the primitive case, a block count for the imprimitive one, and a Latin-square count for the sharp exponent.
>
> Since the argument is yours in outline, I would like you to decide how you prefer to appear: as a coauthor, or with the sketch credited and an acknowledgement. Either is fine with me. I would also be grateful if you could check that the note's account of your argument is accurate, and tell me if you know of an earlier proof.
>
> Once this is settled I intend to let Peter Cameron know, since the question goes back to his remarks on loops and Latin squares.
>
> With best regards,
> [name]

---

## 4. Cameron — the loops note (after 3)

**Why him.** The multiplication-group question and the "second row" formulation of the conjecture both trace to his remarks; he will want to know, and he is a likely reader of both notes.

**What to ask.** Nothing specific; inform, and ask whether he knows an earlier proof or an earlier explicit statement of the belief.

**Attach.** The loops note (in its post-Eberhard form), and optionally the gap note, since he has an interest in the CGW conjecture.

**Do not claim.** Anything about the conjecture itself beyond what the gap note says.

**Draft.**

> Subject: The multiplication group of a random loop
>
> Dear Peter,
>
> You may like to know that the oft-repeated belief that almost every loop has multiplication group S_n now has a proof. Cavenagh, Greenhill and Wanless remarked in 2008 that it would follow from their second-row conjecture; the attached note proves it directly, with the sharp exponent, along lines Sean Eberhard sketched in 2016 [and with him as coauthor / with his sketch credited]. I would be glad to hear of any earlier statement or proof of the belief that you know of.
>
> [Optional: Separately, in the course of work on the second-row conjecture itself I found what appears to be a gap in the proof of the upper bound in Lemma 3.12 of the 2008 paper; I have written to the authors, and attach the note in case it is of interest. The statement is not in doubt.]
>
> With best regards,
> [name]

---

## 5. (Optional, later) Kwan, Petrova, Sawhney — Remark 6.6 and the dependence of CW16 on the gap

**Why them.** The gap note says KPS's written proof of Theorem 1.3(3)–(4) uses CW16 as a black box, and that their Remark 6.6 sketches an independent re-proof of the Cavenagh–Wanless parity theorem, which the note has not checked. If the sketch is complete, CW16's theorem — and with it KPS 1.3(3)–(4) — stops depending on the gap, which simplifies the "Other affected papers" section. This is a concrete, bounded question they can answer in a line.

**What to ask.** Whether Remark 6.6 is a complete argument for the mod-2 equidistribution, so that Theorem 1.3(3)–(4) can be read as independent of CW16. (Do not ask them about the hazard question; the previous cold read's assessment is that their approximation lemma, with its e^{O(n log² n)} cost, is not the tool for a constant-probability event at O(1) selected positions, and asking would be asking them to do the work of §8.)

**Attach.** The gap note, after the CGW authors have had it.

**Draft.**

> Subject: A question about Remark 6.6 of "Parities in random Latin squares"
>
> Dear Matthew, Kalina and Mehtaab,
>
> I have found what I believe is a gap in the proof of Lemma 3.12 of Cavenagh–Greenhill–Wanless (2008), described in the attached note, which I have sent to the authors. Among the results that depend on it is Theorem 1.1 of Cavenagh–Wanless (2016), which your Theorem 1.3(4) contains, so the statement is safe; but your written proof of parts (3)–(4) uses that theorem as a black box, and your Remark 6.6 sketches an independent proof of it by your own method. May I ask whether you regard that sketch as complete, so that I can describe parts (3)–(4) as independent of the 2008 lemma? The note currently says the dependence is "presumably removable" and that I have not checked the sketch; I would rather state it correctly.
>
> With best regards,
> [name]

---

## Notes on timing and tone

- Send 1 before anything else goes public (the conditional paper cites the gap note).
- 2 can go within days of 1. If Allsop or Morris answer (i) negatively, §7 needs revisiting before the paper is posted, so it is worth asking before posting.
- 3 and 4 are independent of the conjecture and can go now; the "folklore" sentence in the loops note should be settled before 4.
- 5 is optional and best sent after 1 has been acknowledged.
- None of the drafts states or implies that the conjecture is proved, that any hypothesis is close to proved, or that the gap makes any published statement false. Keep it that way in follow-ups.
