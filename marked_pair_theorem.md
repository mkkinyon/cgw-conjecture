# The marked-pair reformulation (session 5 working notes)

Setting: first row id, so the CGW column map satisfies omega = sigma
(1∘omega(j) = 2∘j with row 1 the identity gives omega(j) = sigma(j)).

Fix mu with parts alpha != beta, lambda = mu with alpha,beta merged.

**Splitting instance space** X = {(L', P)}: L' in S(lambda), P = {j, j'}
an unordered column pair with j' = sigma^alpha(j), j in a merged
(alpha+beta)-cycle of sigma(L').  |X| = |S(lambda)| * lambda_{a+b} * (a+b).
X_A / X_B = subsets where P is an A-/B-pair of L' (rows 1,2 in
different/same column cycle of P).

**Joining instance space** J = {(L, P)}: L in S(mu), P = {j,j'} spanning a
chosen alpha-cycle and beta-cycle.  J_A / J_B by the status of P in L.

## Lemma S (full-cycle trades preserve pair status)
Trading (swapping the two columns' entries along) a FULL column cycle C of
the pair P inverts pi_P on C and fixes it elsewhere; the partition of rows
into column cycles of P is unchanged.  Hence the A/B status of P is
invariant under any full-column-cycle trade at P — in particular under the
CGW flip at P (which trades the full column cycle through row 2).
Proof: for r in C, new pi(r) = old pi^{-1}(r) (direct computation); cycles
are inverted in place, partition unchanged.

## Corollary 1 (joins land in X_A; unflips are A-joins)
The flip at an A-pair keeps it an A-pair (Lemma S).  So:
- join(L, P) (= flip after switch-if-B, the switch making P an A-pair)
  always produces a marked-A instance:  image(join) ⊆ X_A.
- unflip(x) for x in X_A (flip at the marked A-pair; splits the merged
  cycle into an alpha- and beta-cycle) produces (L, P) with P an A-pair
  of L: every unflip is an A-join.

## Lemma T (switch toggle on J)
The CGW switch at P (trade the row cycle through j or j', canonical
min-symbol choice) is an involution on J (Lemma 3.8 of CGW; the traded
cycle's symbol set is preserved, so the min-symbol rule re-selects it)
that TOGGLES the status of P and preserves the instance data (it inverts
one sigma-cycle in place; type and cycle column-sets preserved).
Hence J_A = J_B exactly.

## Theorem M (marked-pair identity)
    J_A = J_B = |X_A|,   hence   J = 2|X_A|,
and with the CGW multigraph balance  2|X| = J(1 + qbar):

    qbar(mu -> lambda) = |X_B| / |X_A|            (exact)

i.e.  C_n(lambda)/C_n(mu) = (1+qbar)/2 = |X|/(2|X_A|), and

    EQ(delta)  <=>  | |X_B|/|X_A| - 1 | <= delta
               <=>  P_X(marked pair is B) = 1/2 + O(delta).

Proof of J_A = |X_A|: x -> unflip(x) maps X_A into A-joins (Cor 1);
flip∘flip = id at a fixed pair (the flipped cycle re-contains row 2 and is
re-selected), so join(unflip(x)) = x and unflip(join(L,P)) = (L,P) for
A-joins: bijection.  Composing with the switch bijection (Lemma T) maps
X_A bijectively onto B-joins: x -> switch(unflip(x)); its join-image is
flip(switch(switch(unflip x))) = flip(unflip x) = x.  QED

## Consequence: the conjecture in one-square form
With first row id and uniform L in S(lambda):

    P( rows 1,2 lie in the same column cycle of {j, sigma^alpha(j)} )
      = qbar/(1+qbar) = 1/2 - Theta(n^-3)   (conjecturally),

j marked uniformly in a merged (alpha+beta)-cycle.  The post-flip events
B1, B2 and the joining measure have been eliminated: EQ is now a
statement about ONE conditioned bit of ONE uniform square.

## Status of verification
- n=5 (2,3)->(5): J_A=J_B=|X_A|=1440, |X_B|=2880, qbar=2 ✓ exact.
  unflip status all-A ✓, roundtrip ✓, join images all marked-A ✓,
  each marked-A instance hit exactly twice (one A-join, one B-join) ✓,
  X_B entirely orphaned ✓.
- n=6 (2,4)->(6) [expect qbar=10/11] and (2,2,2)->(2,4) [expect 4/7,
  alpha=beta=2 with doubled conventions]: RUNNING (pf_bridge.py).
- Cross-switch (CGW Def 3.9) implemented in pf_orbit.py; at n=5 it is an
  involution on X_B preserving the marked-B status and toggling the
  status of BOTH shifted pairs P_{+1}, P_{-1} deterministically
  (k=1 and k=alpha+beta-1 toggle 2880/2880; middle k's only sometimes).
  This gives P(B1 | marked B) = P(B2 | marked B) = 1/2 exactly — now
  mostly of historical interest since B1,B2 are eliminated, but the
  cross-switch remains the only involution acting on X_B.

## The remaining core problem
|X_B|/|X_A| -> 1.  The status-toggling switch exists on J (two row
cycles) but NOT on X (one row cycle); the obstruction to toggling the
marked status inside S(lambda) IS the conjecture.  Natural attacks:
(1) local analysis of pi_P around rows 1,2 given the sigma-marking
    (pair calculus; completion-uniform measure);
(2) KPS re-randomisation of the single bit 1_B(P) via stable
    intercalates avoiding the P-environment (now with no conditioning
    on a join);
(3) quantify the defect of candidate near-toggles on X, e.g.
    partial-cycle trades at P patched elsewhere (the failure
    configurations should have measure O(n^-2) x imbalance O(n^-1)).
