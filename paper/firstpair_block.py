"""firstpair_block.py -- paper paragraph (e) of the orbit section: the first ladder pair alone
(prop:firstpair, prop:firstpairbias, hyp:firstpair, data).  Source: notes section (k)."""

firstpair_block = r"""

\paragraph{(e) The first ladder pair alone.}  The adaptive coordinates of (c$'$)
specialise to the single pair $(x_1,y_1)=(\pi^{-1}(1),\pi^{-1}(2))$ with none of the side
conditions of Proposition~\ref{prop:adaptive} (there are no other pairs to protect), and
the resulting involution combines with the first ladder trade $T_1$ into a bound on
$\Pr_X[B]-\Pr_X[A]$ by an event selected by that one pair.  Since $\alpha,\beta\ge2$ the
pair always consists of rows outside $\{1,2\}$: $x_1=2$ would mean $p'=\sigma(p)$ and
$y_1=1$ would mean $p=\sigma(p')$ (so $K_0\ge2$ always).  Write $\rho=\rho_{x_1,y_1}$, $F$
for the event that $(x_1,y_1)$ is flippable (parallel), and $\nu$ for the minimum of the
two arc lengths from $p$ to $p'$ and from $p'$ to $p$ (crossed case) or of $|Q_p|,|Q_{p'}|$
(parallel case) --- the quantity written $\mu_{x_1,y_1}$ in (c$'$); we use $\nu$ here
because $\mu$ is the split type elsewhere in the paper.

\begin{lemma}\label{lem:mu2}
$\nu\ge2$ always.
\end{lemma}
\begin{proof}
$\rho(p)=p'$ would mean $L(y_1,p')=L(x_1,p)=L(1,p')$, impossible in a column; $\rho(p')=p$
would mean $L(x_1,p')=L(y_1,p)=L(2,p')$, again impossible in a column; and $\rho$ has no
fixed points since $x_1\ne y_1$.
\end{proof}

For $1\le t<\nu$ the column pair $\{q_t,c_t\}=\{\rho^t(p),\rho^t(p')\}$ (written
$q^t,c^t$ in (c$'$)) is free, consists of two distinct columns, and is separated by the
mark, as observed in (c$'$).  Let $\tau_t=\delta_{q_t,c_t}$ be the permutation of rows
induced by the column pair ($\tau_t(r)$ is the row of $L(r,q_t)$ in column $c_t$) and
$Z_t(r)$ its cycle through row $r$.  Call $t$ \emph{good} if $x_1$ and
$y_1$ lie in different cycles of $\tau_t$ and at least one of $Z_t(x_1)$, $Z_t(y_1)$ meets
$\{1,2\}$ in $0$ or $2$ rows (so that its turn is legal); let
\[
G=\{\text{some }t<\nu\text{ is good}\},
\]
and on $G$ let $t_0$ be the least good $t$ and $Z$ the admissible cycle at $t_0$ (the one
through $x_1$ if it is admissible, otherwise the one through $y_1$).

\begin{proposition}[single-pair involution]\label{prop:firstpair}
Turning $Z$ (exchanging the entries of columns $q_{t_0},c_{t_0}$ in the rows of $Z$) is an
involution of $X\cap G$ that preserves the frame $\pi$ (hence $x_1,y_1$, $A$ and $B$),
toggles $F$, and fixes $t_0$ and $Z$.  Consequently
$\Pr_X[F\cap G]=\Pr_X[F^c\cap G]$ and
$\bigl|\Pr_X[F]-\tfrac12\bigr|\le\tfrac12\Pr_X[G^c]$.
\end{proposition}

\begin{proof}
Columns $p,p'$ are untouched, so $\pi$, $x_1,y_1$, the mark and the rows outside $Z$ are
unchanged.  Rows $1,2$ are both in $Z$ or both outside; if both are in, they become
$L_1\circ(q\,c)$ and $L_2\circ(q\,c)$ (row $1$ is then no longer the identity, which the
unnormalised definition of $X$ allows), so the new row permutation is
$(q\,c)\,\sigma\,(q\,c)$, of the type of $\sigma$, and as $q,c\notin\{p,p'\}$ it still has
$p'$ at distance $\alpha$ from $p$: the turned square is in $X$.  Exactly one of $x_1,y_1$
is in $Z$ (they lie in different $\tau_{t_0}$-cycles), so $\rho'=\rho\circ(q\,c)$ if
$x_1\in Z$ and $\rho'=(q\,c)\circ\rho$ if $y_1\in Z$.  In the first case
Lemma~\ref{lem:adaptinv} applies with $t=t_0$.  In the second,
$(q\,c)\circ\rho=\rho\circ(u\,u')$ with $u=\rho^{t_0-1}(p)$, $u'=\rho^{t_0-1}(p')$
(possibly $t_0=1$ and $\{u,u'\}=\{p,p'\}$), and the computation in the proof of
Lemma~\ref{lem:adaptinv}, read with $t_0-1$ in place of $t$ (for $t_0=1$ the paths are
empty and the computation is immediate), shows that the $\rho$-paths from $p$ to $u$ and
from $p'$ to $u'$ are unchanged, that $\rho'(u)=c$ and $\rho'(u')=q$ (the
ordered pair at step $t_0$ is reversed, the unordered pair is $\{q,c\}$), and that crossed
with arcs $(a,b)$ becomes parallel with cycle lengths $(b,a)$ and conversely.  In both
cases $\nu'=\nu$, the candidate pairs are the same for $t\le t_0$, and $F$ toggles.  The
column pairs with $t<t_0$ are disjoint from $\{q,c\}$, so their $\tau_t$, hence their
badness, are unchanged.  At $t_0$ the turn replaces $\tau_{t_0}$ by $\tau_{t_0}$ with the
cycle $Z$ reversed (or by the inverse of that, if the ordered pair is reversed), so $x_1,y_1$
are still in different cycles, $Z$ is still a cycle with the same intersection with
$\{1,2\}$, and the other cycle through $\{x_1,y_1\}$ is unchanged: $t_0$ and $Z$ are the
same in the turned square, and the turn is an involution.  The identity follows, and
$\Pr[F]\ge\tfrac12\Pr[G]$, $\Pr[F^c]\ge\tfrac12\Pr[G]$ give the bound.
\end{proof}

\begin{proposition}[the first ladder pair controls the bias]\label{prop:firstpairbias}
For every $n,\lambda,\alpha,\beta$,
\[
\Pr_X[B]-\Pr_X[A]=\E_X\bigl[(\mathbf 1_B-\mathbf 1_A)(\mathbf 1_{F^c}-\mathbf 1_F);\,G^c\bigr];
\]
consequently $\bigl|\Pr_X[B]-\tfrac12\bigr|\le\tfrac12\Pr_X[G^c]$,
$\Pr_X[A],\Pr_X[B]\ge\tfrac12\Pr_X[G]$, and
$|\bar q(\mu\to\lambda)-1|\le2\Pr_X[G^c]/(1-\Pr_X[G^c])$.
\end{proposition}

\begin{proof}
Write $B_{FG}=\Pr_X[B\cap F\cap G]$ and so on.  The trade $T_1$ of
Theorem~\ref{thm:ladder} ($k=1$) is an involution of $X\cap F$ whose frame is
$(x_1\,y_1)\circ\pi$, so it exchanges $A$ and $B$:
$B_{FG}+B_{FG^c}=A_{FG}+A_{FG^c}$.  The turn of Proposition~\ref{prop:firstpair} is an
involution of $X\cap G$ exchanging $F$ and $F^c$ and preserving $\pi$, hence $A$ and $B$,
which are functions of $\pi$: $B_{FG}=B_{F^cG}$ and $A_{FG}=A_{F^cG}$.  Hence
\begin{align*}
\Pr[B]-\Pr[A]&=(B_{F^cG}+B_{F^cG^c})-(A_{F^cG}+A_{F^cG^c})
=(B_{FG}-A_{FG})+(B_{F^cG^c}-A_{F^cG^c})\\
&=(A_{FG^c}-B_{FG^c})+(B_{F^cG^c}-A_{F^cG^c}),
\end{align*}
which is the displayed expectation; its absolute value is at most $\Pr[G^c]$.  ($T_1$ need
not preserve $G$.)  The last inequality is $|\bar q-1|=|\Pr[B]-\Pr[A]|/\Pr[A]$ with
$\Pr[A]\ge\tfrac12(1-\Pr[G^c])$.
\end{proof}

\begin{hypothesis}[first-pair hypothesis]\label{hyp:firstpair}
$\displaystyle\sup_{\lambda,\alpha,\beta}\Pr_{X(n,\lambda;\alpha,\beta)}[G^c]=o(1/\log n)$.
\end{hypothesis}

\begin{corollary}\label{cor:firstpair}
Hypothesis~\ref{hyp:firstpair} implies Conjecture~\ref{conj:cgw}, with
$\dtv(\PP_n,\QQ_n)=O(\log n\cdot\sup\Pr_X[G^c]+n^{-1+o(1)})$.
\end{corollary}
\begin{proof}
With $s=\sup\Pr_X[G^c]$, Proposition~\ref{prop:firstpairbias} gives
$|\bar q(\mu\to\lambda)-1|\le2s/(1-s)$ for every adjacent pair, i.e.\ $\mathrm{EQ}(\delta)$
with $\delta=2s/(1-s)=o(1/\log n)$, and Theorem~\ref{thm:reduction} applies.
\end{proof}

\paragraph{What the first-pair hypothesis asks, and what it costs.}  The event $G^c$ is
selected by the first ladder pair but is not local to it: it involves the column cycles of
the pairs $\{q_t,c_t\}$ through all rows, and rows $1,2$ through admissibility.  It is
again a lower tail of a count of events at square-selected positions --- the obstacle is
unchanged in kind --- and it decomposes as
\[
\Pr_X[G^c]\le\Pr_X[\nu\le M]+\Pr_X[\nu>M,\ t=1,\dots,M\text{ all bad}].
\]
The first term is an upper bound on a short arc or cycle of $\rho_{x_1,y_1}$ at $p,p'$;
nothing of the kind is proved at a frame-selected pair (Proposition~\ref{prop:splice}
bounds a different frame event), and heuristically it is $\approx4M/n$.  The second term
--- among $M$ candidates on disjoint free column pairs not all are bad --- is the
obstacle itself.  Compared with Hypothesis~\ref{hyp:witness}, at the level of heuristics
and data only: one row pair instead of $\ge n/\log^2n$; per-event probability $\approx0.4$
(the random-permutation value for ``$x_1\not\sim y_1$ under $\tau_t$ and one of the two
cycles has even intersection with $\{1,2\}$'') instead of $\ell'/n$; expected count
$\approx0.4\,\E\nu$ with $\E\nu=6.3,9.6,18$ at $n=30,50,100$ in the data
(Remark~\ref{rem:fixrect}), instead of $\ge\log^2n$ ($\log^3n$ at the parameters of
Proposition~\ref{prop:Lwitness}); no offset families, no splice bound, no conditioning on
$K_0$.  Neither hypothesis is known to imply the other: $\Pr_X[F^c]\ge\tfrac12\Pr_X[G]$, so
the first pair alone cannot give (L), and Hypothesis~\ref{hyp:witness} says nothing about
$G$.  So Hypothesis~\ref{hyp:firstpair} is a third alternative sufficient condition, the one
with the fewest moving parts; the reference hypothesis of the paper remains
Hypothesis~\ref{hyp:witness}.  The minimal instance \eqref{eq:minimal} of
\S\ref{sec:remains} follows from Hypothesis~\ref{hyp:firstpair} (with $c=\tfrac12-o(1)$)
but not evidently from Hypothesis~\ref{hyp:witness}, and is sufficient for neither.

\begin{remark}[data: completions of a fixed rectangle]\label{rem:fixrect}
To test the hypotheses on types that uniform sampling never produces we ran the
Jacobson--Matthews chain restricted to the completions of a fixed $2\times n$ rectangle
(\texttt{jmfix.c}: moves touching rows $1,2$ are never proposed from a proper state and are
rejected from an improper one; the walk is a simple random walk on the restricted move
graph, every proper state having restricted degree $(n-2)n(n-3)$, so its stationary law is
uniform on the proper states of the component of the start).  Connectivity of the restricted
graph is not proved; at $n=5,6$ (types $(5)$, $(6)$, $(2,2,2)$; $36$, $4032$, $5376$
completions, enumerated exactly) every completion is reached from two different starts with
$\chi^2$ statistics within one standard deviation of uniformity, and at $n=7$ the sampled
$\Pr[B\mid R]=0.4965$ agrees with the exact $0.497076$ over all $6566400$ completions of
$R=(\mathrm{id},(1\,2\cdots7))$.  Mixing is not analysed; runs with different thinning
and starts agree within noise.  Any two rectangles of type $\lambda$ are isotopic by a
column permutation and a symbol relabelling fixing rows $1,2$, which preserve completion
counts and carry a mark of class $(m,\alpha)$ to one of the same class; every event
considered here is defined through the cycle structures of row and column permutations and
is invariant under such isotopies.  Hence $\Pr_X$ of such an event equals its probability
under a uniform completion of any one rectangle $R$ of type $\lambda$ with any one mark of
the class, and averaging over the marks of the class in $R$ (which are equivalent under the
isotopies preserving $R$) only reduces the variance.  Findings (assuming connectivity and
mixing of the sampler) ($n=30$: $6\cdot10^4$ squares per type; $n=50$: $4\cdot10^4$;
$n=100$: $1.2\cdot10^4$; types $(n)$, $(n/2,n/2)$, $(n-2,2)$, $(10,10,10)$, $(4,2^{13})$,
$(6,3^8)$, $(4,2^{23})$, $(6,3^{14},2)$, $(4,2^{48})$; marks at $\alpha=2$,
$\approx m/3$, $\lfloor m/2\rfloor$):
(i) $\Pr_X[B]\in[0.4974,0.5024]$ in every class at $n\le50$ (sampling error
$\approx0.0025$), i.e.\ by Theorem~\ref{thm:marked}
$\CC_n(\lambda)/\CC_n(\mu)\in[0.995,1.005]$ for every tested adjacent pair, including
$(4,2^{13})\to(2^{15})$ and $(6,3^8)\to(3^{10})$, types of $\PP_n$-probability
$e^{-\Theta(n)}$; $[0.495,0.504]$ at $n=100$ (error $\approx0.005$).
(ii) $\Pr_X[F]\in[0.4988,0.5017]$ in every class at $n\le50$ ($0.4896$ at $n=7$, so
$\tfrac12$ is not an identity).
(iii) Along the ladder, $\Pr[\text{no witness among the first }K\text{ pairs}\mid K_0-1\ge K]$
equals $(1-w_1)^K$ within sampling error for $K\le12,20,40$ at $n=30,50,100$ and
thresholds $\ell'=4,8,16,24$, where $w_1$ is the first-pair witness rate ($0.100$, $0.060$,
$0.030$ at $\ell'=4$); the hazard $\Pr[\text{witness at }k\mid\text{none before }k]$ is
flat in $k$ to three digits; $\Pr[\text{none of the first }K\text{ ladder pairs flippable}\mid K_0-1\ge K]=2^{-K}$
within noise to $K=12$--$15$ ($7\cdot10^5$ instances at $n=30$); across the diagonals $D_c$ of
$B$-instances with $d_{12}=\ell\le n/4$, which share their $y$-rows, the number of
diagonals carrying a witness vanishes with the binomial probability and has binomial
lower tails, for all $\ell\le12$.
All of this is identical for the extreme types.
(iv) For Hypothesis~\ref{hyp:firstpair}, with all candidates $t<\nu$: $\Pr_X[G^c]=0.61,
0.18, 0.11, 0.055$ at $n=7,30,50,100$ (type $(n)$, $\alpha=2$; $0.056$ for $\alpha=50$ at
$n=100$), i.e.\ $\approx5.5/n$ for $n\ge30$; $\E\nu=6.3,9.6,18$ and $\Pr[\nu=k]\approx4/n$
for small $k$; $\Pr[G^c\mid\nu=k]=0.59,0.37,0.22,0.13,0.08$ for $k=2,\dots,6$ at $n=30$
against $0.59^{k-1}$ (a slight positive dependence among the candidates, which share the
rows $x_1,y_1$); $t=1$ is good with probability $0.41$ in every class at $n\le100$.  Other
types were measured only with $t\le8$: $\Pr[\text{some }t\le8\text{ good}]=0.81,0.88,0.94$
at $n=30,50,100$ in every class.  The identity of
Proposition~\ref{prop:firstpairbias} and its two ingredients hold to all digits on the
exact enumerations at $n=5,6$ (every class) and within sampling error at $n=7$ ($-0.0052$
and $-0.0058$, each $\pm0.001$, against the exact $-0.005848$).
None of these data reach the regime of any hypothesis ($\ell'\ge\log^5n$, failure
$o(1/\log n)$); they exclude a visible dependence between the events, or a visible effect
of the type, at the first- and second-moment level for $n\le100$.
\end{remark}
"""
