"""cprime_blocks.py -- fresh text for the conditional paper: adaptive coordinates (section 5(c')), the hypothesis,
prop:Lorbit, data, and the new section 6.  Imported by assemble_conditional.py."""

cprime_intro = r"""
\paragraph{(c$'$) Adaptive coordinates: the toggling is automatic.}
Proposition~\ref{prop:orbit} needs the clean coordinates of a ladder pair to be
\emph{separated} by $\{p,p'\}$ with positive probability, a genericity input we cannot prove
(\S\ref{sec:remains}).  It disappears if the coordinates are chosen along $\rho_{x,y}$ from
the mark.  For a pair $(x,y)$ of the family and $t\ge1$ let
\[
(q^t,c^t):=\bigl(\rho_{x,y}^{\,t}(p),\ \rho_{x,y}^{\,t}(p')\bigr),
\]
and let $\mu_{x,y}:=\min(d(p\to p'),d(p'\to p))$ if $(x,y)$ is crossed and
$\min(|Q_p|,|Q_{p'}|)$ if parallel.  For $t<\mu_{x,y}$ the pair $(q^t,c^t)$ is free
($q^t\ne c^t$, both $\notin\{p,p'\}$) and separated by $\{p,p'\}$ in $\rho_{x,y}$ by
construction: on different arcs when crossed, one in $Q_p$ and one in $Q_{p'}$ when parallel.
So the turn of its column cycle $Z$ through $x$ (with $y\notin Z$) toggles the status of
$(x,y)$ whenever it is legal (Lemma~\ref{lem:legalturn}).  What makes these coordinates
usable in the orbit method is that they are invariant under their own turn.  (We keep (c)
for comparison: its orbit weight is what the data of (d) measure, and it is where the
genericity constant of \S\ref{sec:remains} lives.)
"""

cprime_rest = r"""
\begin{hypothesis}[short clean cycles at the ladder rows]\label{hyp:adaptive}
There are $T,\ell_1=O(\log^2n)$ such that, uniformly in $\lambda,\alpha,\beta$, with
$\ell_0=n/\log^2n$ and $U$ as in Proposition~\ref{prop:adaptive}:
\begin{itemize}[leftmargin=2em]
\item[(i)] for the sub-ladder $I$ of the first $\min(K_0-1,\lceil n/\log^2n\rceil)$ ladder
pairs, $\E_X\bigl[2^{-|U|};\ K_0\ge\ell_0\bigr]=o(1/\log n)$;
\item[(ii)] for a $B$-instance with $d_{12}=\ell<\ell_0$ and $|C_1|\ge n/2+\ell$, with
$S=\{2\ell t:\ 0\le t\le(n/2-\ell)/(2\ell)\}$ (offsets at gaps $2\ell$, so that the
$x$-rows $\pi^{-(c+i)}(1)$, $c\in S$, $1\le i<\ell$, are pairwise distinct and
$|S|\ge n/(8\ell)$) and $U$ taken for the family $\bigcup_{c\in S}D_c$: the event that fewer
than half of the diagonals $D_c$, $c\in S$, meet $U$ has probability $o(1/\log n)$ (as a
single event over all $\ell<\ell_0$); and symmetrically for $d_{21}$.
\end{itemize}
\end{hypothesis}

\begin{proposition}\label{prop:Lorbit}
Hypothesis~\ref{hyp:adaptive} implies \textup{(L)}, hence the weak CGW conjecture.
\end{proposition}

\begin{proof}
By Theorem~\ref{thm:ladder},
$|\Pr[B]-\Pr[A]|\le\Pr[\mathrm{Flip}=\emptyset]\le\Pr[K_0\ge\ell_0,\mathrm{Flip}=\emptyset]+\Pr[A,K_0<\ell_0]+\Pr[B,K_0<\ell_0]$.
\emph{Long ladders.}  $K_0$ is a function of the frame, hence orbit-invariant, and by
Proposition~\ref{prop:adaptive} $\Pr[\mathrm{Flip}_I=\emptyset\mid\text{orbit}]\le2^{-|U|}$,
so $\Pr[K_0\ge\ell_0,\mathrm{Flip}=\emptyset]\le\E[2^{-|U|};K_0\ge\ell_0]=o(1/\log n)$ by (i).
\emph{Short ladders on $A$.}  $\Pr[A,K_0<\ell_0]\le72\ell_0/n$ by Proposition~\ref{prop:trapped}.
\emph{Short arcs on $B$.}  Fix $\ell<\ell_0$ and consider $B$-instances with $d_{12}=\ell$.
The part $|C_1|<n/2+\ell$ costs at most $8/n$ by Proposition~\ref{prop:splice}.  On
$|C_1|\ge n/2+\ell$ apply Proposition~\ref{prop:adaptive} to the family
$\bigcup_{c\in S}D_c$: its $x$-rows lie in the disjoint depth intervals $[c+1,c+\ell-1]$,
$c\in S$, on the frame arc from $2$ to $1$, and its $y$-rows $\pi^{-i}(2)$ on the arc from
$1$ to $2$, so the family is as required.  Given the orbit the statuses of all pairs are
independent and $T_k$ toggles pair $k\in U$, so a diagonal $D_c$ with $U\cap D_c\ne\emptyset$
is flippable (contains a parallel pair) with conditional probability $\ge\tfrac12$, and the
diagonals' flippabilities are independent (they are disjoint sets of pairs).  Off the
exceptional event of (ii), at least $|S|/2$ diagonals meet $U$, so the number of flippable
offsets $c\in S$ has conditional mean $\mu\ge|S|/4$, and by Chernoff it is $\ge|S|/8$ except
with conditional probability $\le e^{-\mu/8}\le e^{-|S|/32}\le e^{-\log^2n/256}$ (as
$|S|\ge n/(8\ell_0)=\log^2n/8$), which is summable over $\ell<\ell_0$ to $o(1/n)$.  The
identity behind \eqref{eq:offsetavg} holds for the average over $c\in S$ as well as over all
$c$, so with $\theta_0=\tfrac18$,
$\sum_{\ell<\ell_0}\Pr[B,d_{12}=\ell,|C_1|\ge n/2+\ell,\theta_S\ge\theta_0]\le36\ell_0/((n-2)\theta_0)=O(1/\log^2n)$.
The case $d_{21}=\ell$ is symmetric.  Altogether $\Pr[\mathrm{Flip}=\emptyset]=o(1/\log n)$.
\end{proof}

\paragraph{What the hypothesis says, and what it costs.}  A ladder pair enters $U$ when one
of its $T$ candidate column pairs has a short column cycle through $x_k$ that is clean ---
a configuration of $O(\log^2n)$ cells ($4t$ cells of rows $x_k,y_k$ to locate the
candidate, and the cycle itself).  Heuristically a candidate's column cycle has length
$\le\ell_1$ with probability $\approx\ell_1/n$ (the data of \S(d) show the lengths spread
like a uniform variable on $\{2,\dots,n\}$), is clean with probability
$\approx e^{-2s\ell_1/n}$ for a family of $s$ pairs, and the chosen pair avoids the
$\approx2Ts$ candidate columns of the other pairs with probability $\approx(1-2Ts/n)^2$, so
\[
\E|U|\;\approx\;s\cdot\frac{T\ell_1}{n}\,e^{-2s\ell_1/n}\Bigl(1-\frac{2Ts}n\Bigr)^2\;\asymp\;\log^2n
\qquad\text{for }s=\frac n{\log^2n},\ T=\frac{\log^2n}4,\ \ell_1=\frac{\log^2n}2,
\]
against the $C\log\log n$ ($C>1/\log2$) that (i) needs; the constant is small
($\approx1/80$), so the heuristic excess over the threshold is a large-$n$ statement.  Compared with the fixed matching
of (c), which offers $\approx n/2$ coordinates per pair of which a constant number are short
and clean and a fraction $\approx0.19$ of those separated, the adaptive rule offers only
$T=O(\log^2n)$ coordinates per pair, all separated: the genericity constant is gone, and in
exchange membership in $U$ is a rare event per pair, so (i) is a lower-tail statement about
rare local events at rows selected by the square --- of the same shape as a statement that
some ladder pair has a short $\rho_{x_k,y_k}$-cycle through $p$ avoiding $p'$ (which would
make it flippable outright), now with a column cycle of a free pair through $x_k$ in its
place.  There is a second regime:
with $s\approx T\approx\ell_1\approx\sqrt n$ a pair enters $U$ with probability bounded
below and $\E|U|\asymp\sqrt n$, so (i) becomes a second-moment statement about $\sqrt n$
events of constant probability, at the price of the arcs of $\rho_k$ at $p,p'$ exceeding
$\sqrt n$ for a constant fraction of the pairs.  In both regimes the lengths of the arcs of $\rho_k$ at $p,p'$ --- which the
fixed-matching version needed to be macroscopic for a positive orbit weight --- enter only
through the number of candidates $\min(T,\mu_k-1)$.  The collision rule in $U$ (the chosen pair must avoid the candidate
columns of \emph{all} other pairs) is what makes $U$ orbit-invariant --- a turn can shrink
the candidate set of an excluded pair, so the obvious greedy refinement is not invariant ---
and it is crude: it is why the asymptotic regime $2Ts/n\le\tfrac12$ is out of reach of
sampling at $n\le100$.

\paragraph{(d) Data.}  All scripts and logs are in the repository \cite{repo}, under
\texttt{s13\_*.py} and \texttt{runs/s13/orbit/}.  \emph{Adaptive coordinates.}  On
true-marked instances at $n=30,50,100$ (family the first $\min(K_0-1,n/8)$ ladder pairs,
$T=6$, $\ell_1=6,8,10$; $2400$, $3000$, $900$ instances), for a random $k\in U$ the turn of
$Z_k$ was applied and checked to toggle the status of $(x_k,y_k)$, to leave the status,
$\mu_l$, $t_0(l)$, $C_l$ of every other pair and the set $U$ unchanged, and to preserve the
type of rows $1,2$ and the frame ($682$, $697$, $215$ turns, no failure; $444$ and $160$
more with families of $12$ and $30$ pairs); the bound
$\Pr[\mathrm{Flip}_I=\emptyset]\le\E[2^{-|U|}]$ holds in every $K_0$ class.  At these sizes
$30\%$, $23\%$, $16\%$ of the ladder pairs have a good candidate and $13\%$, $5\%$, $3\%$
are in $U$ ($8\%$ and $1.4\%$ with the family of $30$ pairs at $n=100$): $2Ts\ge n$ there,
so the collision rule removes most pairs.  \emph{Column-cycle lengths
and the fixed-matching coordinates.}  The length of the column cycle of a free pair through
a ladder row is spread like a uniform variable on $\{2,\dots,n\}$ (about $1/n$ per length;
the turn is legal for $98\%$ of the $2$-cycles and $64\%$ of the long cycles, and
avoids $y$ for $97\%$ resp.\ $27\%$ of them, at $n=100$).  For the fixed matching of
(c): $\Pr[(q,c)\text{ separated}\mid\text{clean, length }\ell]=0.193$--$0.201$ for
$\ell=2,\dots,8$ against $0.195$ over all column pairs ($n=100$), i.e.\ no detectable
dependence on the cycle length; with $|\mathcal T|\approx42$ coordinates the orbit average
of the parallel bit is $0.478$ from a crossed and $0.523$ from a parallel start ($n=100$);
with the fixed matching of consecutive free columns $|\mathcal T|=1.7$ resp.\ $3.0$ clean
coordinates per pair ($n=50,100$) and $\E[\bar y]=0.159$ resp.\ $0.189$; and the orbit
weights of different ladder pairs are uncorrelated to within sampling error
($\mathrm{Var}(Y_I)/(|I|\mathrm{Var}\,\bar y)=0.99$ and $1.14$ at $n=50,100$, the excess
coming from the variation of $|I|$).  The ladder data of Remark~\ref{rem:ladderdata}
complete the picture: $\Pr[\mathrm{Flip}=\emptyset]\approx4/n$, as if the ladder coins
were fair and independent.
"""

remains = r"""
\section{What remains}\label{sec:remains}

Hypothesis~\ref{hyp:adaptive} is a statement about short cycles at the frame-selected rows:
that enough ladder pairs $(x_k,y_k)$ have, at one of the prescribed column pairs
$(\rho_k^t(p),\rho_k^t(p'))$, a column cycle through $x_k$ of length $O(\log^2n)$ that
avoids the other ladder rows and meets $\{1,2\}$ in $0$ or $2$ rows; and that the arcs of
$\rho_k$ at $p,p'$ are not shorter than $T$ for most $k$.  Both are local (polylogarithmic)
configurations whose heuristic frequency is, for large $n$, far above what is needed; the difficulty is
entirely that the rows are selected by the frame of the mark and the whole space is
conditioned on the type of rows $1,2$.  We record what the available tools do and do not
give.

\emph{Unconditional facts.}  For two generic rows of a uniform square,
Proposition~\ref{prop:tailsurvive}(i) gives, by Markov's inequality on the number of
columns on short cycles, $\Pr[\text{all cycles of }\rho_{x,y}\text{ have length}\le\varepsilon n]\le2\varepsilon/(1-\varepsilon)$;
so a generic pair has a cycle longer than $n/5$ with probability $\ge\tfrac12$, and the
column-cycle analogue holds by transposition.  Kwan and Sudakov's bound
$\Pr[N_2\ge t]=e^{-\Omega(t\log t)}$ on the number of $2$-cycles of two fixed rows
\cite[Thm.~3]{KS18} is proved by joining switchings of the same kind as those of
Proposition~\ref{prop:tailsurvive} and is likewise an upper bound on short cycles.  The
number of completions of any $2\times n$ rectangle is within $e^{O(n\log^2n)}$ of any other
(van der Waerden and Br\'egman; \cite[Prop.~5]{KS18}), so any property of uniform squares
with failure probability $e^{-\omega(n\log^2n)}$ --- such as the intercalate concentration
of Kwan, Sah and Sawhney \cite{KSS21} --- holds conditionally on any fixed rows $1,2$ and
mark, hence in $X$ for every $\lambda,\alpha,\beta$; but the natural failure probability of a
statement about $m\le n$ prescribed row pairs is $e^{-\Omega(m)}$, which this transfer does
not reach.

\emph{Why the fixed-matching version is harder.}  With a fixed matching the orbit weight
$\bar y_k$ of Proposition~\ref{prop:orbit} needs a lower bound on the probability that a
generic column pair is interleaved with $\{p,p'\}$ in $\rho_{x_k,y_k}$ --- an instance of
the statement ``two given pairs of lines of a random Latin square are interleaved in the
cycle structure of a third pair with probability in $[c,1-c]$'' inside $X$.  Without the
conditioning its upper half is elementary and its lower half reduces to macroscopic cycles
(by conjugation invariance the interleaving probability given the cycle type of $\rho_{x,y}$
is $\tfrac13\sum_i(m_i)_4/(n)_4+2\sum_{i\ne j}(m_i)_2(m_j)_2/(n)_4\le\tfrac13+O(1/n)$, and
$\ge\tfrac13\Pr[W\mid\text{type}]^3-O(1/n)$ by Jensen, $W$ the event that two given columns
share a cycle).  Inside $X$ the exact tools are the free trades of rows other than $1,2$
(which preserve interleaving) and the legal column turns (which are legal exactly off the
interleaved configurations), and every identity they produce is fair: it equates two
expectations with random weights whose lower bounds are again of the same form.  Lower
bounds on ``different cycles'' events are, in the language of \cite{CGW08}, of the type of
the direction of their Lemma~3.12 whose proof has the gap: they require a repair of the
same-cycle configurations with a rigid inverse, and the repairs available in $X$ (trades or
cross-switches of a free pair of third rows) have unbounded multiplicity.  The one case in
which we found a rigid inverse is CGW's own case 3 with $\min(\alpha,\beta)=2$, which gives
$\Pr_X[A]\ge\tfrac14$ for marks at distance $2$ (Remark~\ref{rem:gap}); it does not extend
to other distances for a reason of type rather than multiplicity.  The adaptive
coordinates sidestep the question: the coordinate is placed where it is separated by
construction, and the price is paid in rarity rather than in a constant we cannot prove.
"""
