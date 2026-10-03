"""cprime_blocks.py -- fresh text for the conditional paper: adaptive coordinates (section 5(c')), the hypothesis,
prop:Lorbit, data, and the new section 6.  Imported by assemble_conditional.py."""

cprime_intro = r"""
\begin{hypothesis}[orbit genericity; fixed matching]\label{hyp:orbitgen}
Fix the matching $\mathcal M$ of the free columns in consecutive pairs (in the natural
order of the columns) and a constant $\ell_1$ as in (c),
and for a family $\mathcal G$ and $D\subseteq\mathcal G$ write
$Y_D=\sum_{(x,y)\in D}\bar y_{x,y}$ (orbit weights computed with cleanliness with respect to
all rows of $\mathcal G$).  There is a constant $\kappa_0>0$ such that, with
$\ell_0=n/\log^2n$ and uniformly in $\lambda,\alpha,\beta$: (i) for the sub-ladder $I$ of
the first $\min(K_0-1,\lceil n/8\rceil)$ ladder pairs,
$\E_X\bigl[e^{-Y_I/2};\ K_0\ge\ell_0\bigr]=o(1/\log n)$; (ii) for a $B$-instance, with
$\ell:=d_{12}<\ell_0$ and $|C_1|\ge n/2+\ell$, $S=\{2\ell t:\ 0\le t\le(n/2-\ell)/(2\ell)\}$
and $\mathcal G=\bigcup_{c\in S}D_c$: the event that fewer than $\kappa_0|S|$ of the
diagonals $D_c$, $c\in S$, have $Y_{D_c}\ge\kappa_0(\ell-1)$ has probability $o(1/\log n)$
(one event, not one per $\ell$); and symmetrically for $d_{21}$.  (A fixed small fraction,
not a half: with a fixed matching the orbit weight of a single pair vanishes with
probability $\approx0.7$, \S(d), so for $\ell=2$ only a minority of the one-pair diagonals can
qualify.)
\end{hypothesis}
\noindent Since $\bar y_{x,y}$ is typically of constant size (\S(d): $0.16$--$0.19$ on
average with the fixed matching), (i) asks that the sum of $\approx n/8$ such weights not collapse; its content is a
\emph{genericity constant} --- that a clean coordinate is separated by the mark with
probability bounded below --- which none of the exact tools of this paper provides
(\S\ref{sec:remains}).

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
There are $T,\ell_1=O(\log^2n)$ and a constant $\kappa_2>0$ such that, uniformly in
$\lambda,\alpha,\beta$, with $\ell_0=n/(16\log^2n)$, $M=\lfloor n/\log^4n\rfloor$ and $U$ as
in Proposition~\ref{prop:adaptive}:
\begin{itemize}[leftmargin=2em]
\item[(i)] for the sub-ladder $I$ of the first $\min(K_0-1,\lceil\ell_0\rceil)$ ladder
pairs, $\E_X\bigl[2^{-|U|};\ K_0\ge\ell_0\bigr]=o(1/\log n)$;
\item[(ii)] for a $B$-instance, with $\ell:=d_{12}<\ell_0$ and $|C_1|\ge n/2+\ell$: let
$S=\{2\ell t:\ 0\le t\le(n/2-\ell)/(2\ell)\}$ (offsets at gaps $2\ell$, so that the $x$-rows
of the diagonals $D_c$, $c\in S$, are pairwise distinct; $|S|\ge n/(8\ell)\ge2\log^2n$),
let $S'$ be the first $m_\ell:=\min\bigl(|S|,\lceil n/(\log^2n\cdot\min(\ell-1,M))\rceil\bigr)$
offsets of $S$, let $D_c^{(M)}$ be the first $\min(\ell-1,M)$ pairs of $D_c$, and take $U$
for the family $\mathcal F_\ell=\bigcup_{c\in S'}D_c^{(M)}$ (which has at most
$2n/\log^2n$ pairs).  Then the event that fewer than $\kappa_2\log^2n$ of the diagonals
$D_c^{(M)}$, $c\in S'$, meet $U$ has probability $o(1/\log n)$ (one event, not one per
$\ell$); and symmetrically for $d_{21}$.
\end{itemize}
\end{hypothesis}
\noindent The scaling in (ii) is forced: a diagonal of $\ell-1$ pairs meets $U$ with
probability $\approx\min(\ell-1,M)\,T\ell_1/n$, so for short arcs only a small fraction of
the diagonals can be expected to meet $U$, and the family must be kept to
$O(n/\log^2n)$ pairs for the cleanliness and collision conditions of
Proposition~\ref{prop:adaptive} to have constant cost; $|S'|\min(\ell-1,M)\approx n/\log^2n$
makes the expected number of diagonals meeting $U$ of order $\log^2n$ for every $\ell$.

\begin{proposition}\label{prop:Lorbit}
Each of Hypotheses~\ref{hyp:orbitgen} and~\ref{hyp:adaptive} implies \textup{(L)}, hence
the weak CGW conjecture.
\end{proposition}

\begin{proof}
By Theorem~\ref{thm:ladder},
$|\Pr[B]-\Pr[A]|\le\Pr[\mathrm{Flip}=\emptyset]\le\Pr[K_0\ge\ell_0,\mathrm{Flip}=\emptyset]+\Pr[A,K_0<\ell_0]+\Pr[B,K_0<\ell_0]$,
with $\ell_0$ as in the hypothesis used.  \emph{Long ladders.}  $K_0$ is a function of the
frame, hence orbit-invariant, and by Proposition~\ref{prop:orbit} resp.\
Proposition~\ref{prop:adaptive} $\Pr[\mathrm{Flip}_I=\emptyset\mid\text{orbit}]\le e^{-Y_I/2}$
resp.\ $2^{-|U|}$, so $\Pr[K_0\ge\ell_0,\mathrm{Flip}=\emptyset]=o(1/\log n)$ by (i).
\emph{Short ladders on $A$.}  $\Pr[A,K_0<\ell_0]\le72\ell_0/n$ by Proposition~\ref{prop:trapped}.
\emph{Short arcs on $B$.}  Fix $\ell<\ell_0$ and consider $B$-instances with $d_{12}=\ell$;
the case $d_{21}=\ell$ reduces to it by the involution exchanging rows $1,2$, as in
Proposition~\ref{prop:Lwitness}.  The part $|C_1|<n/2+\ell$ costs at most $8/n$ by
Proposition~\ref{prop:splice}.  On $|C_1|\ge n/2+\ell$ the family of (ii) is as required by
the propositions: its $x$-rows lie in the disjoint depth intervals $[c+1,c+\ell-1]$, $c\in S$,
on the frame arc from $2$ to $1$, and its $y$-rows $\pi^{-i}(2)$ on the arc from $1$ to $2$.
Given the orbit the statuses of all pairs are independent, and the diagonals'
flippabilities are independent (they are disjoint sets of pairs).  Under
Hypothesis~\ref{hyp:orbitgen}, a diagonal with $Y_{D_c}\ge\kappa_0(\ell-1)$ is flippable
with conditional probability $\ge1-e^{-\kappa_0(\ell-1)/2}\ge4\theta_0$,
$\theta_0:=\tfrac14(1-e^{-\kappa_0/2})$; off the exceptional event of (ii) at least
$\kappa_0|S|$ diagonals qualify, so the number of flippable offsets has conditional mean
$\ge4\theta_0\kappa_0|S|$ and the fraction $\theta_S$ is $\ge2\theta_0\kappa_0=:\theta_0'$
except with conditional probability $\le e^{-\theta_0\kappa_0|S|/2}\le e^{-\theta_0\kappa_0\log^2n/16}$
(as $|S|\ge n/(8\ell_0)=\log^2n/8$).
Under Hypothesis~\ref{hyp:adaptive}, a diagonal $D_c^{(M)}$ meeting $U$ is flippable with
conditional probability $\ge\tfrac12$ (the status of a pair $k\in U$ is
$\mathrm{bit}_k\oplus\varepsilon_k$ with $\varepsilon$ uniform on the orbit); off the
exceptional event of (ii) at least $\kappa_2\log^2n$ diagonals of $S'$ meet $U$, so the
number of flippable offsets in $S'$ has conditional mean $\mu\ge\tfrac12\kappa_2\log^2n$ and
is $\ge\mu/2$ except with conditional probability $\le e^{-\mu/8}\le e^{-\kappa_2\log^2n/16}$;
hence the fraction $\theta_{S'}$ is $\ge\theta_0(\ell):=\kappa_2\log^2n/(4m_\ell)$.  In both
cases the conditional failure probabilities are summable over $\ell<\ell_0$ to $o(1/n)$.
The derivation of \eqref{eq:offsetavg} from \eqref{eq:offsetprice} applies to the average
over any set of offsets determined by $\ell=d_{12}$ alone, such as $S$ and $S'$, so
\[
\sum_{\ell<\ell_0}\Pr\bigl[B,d_{12}=\ell,|C_1|\ge n/2+\ell,\ \theta_{S}\ge\theta_0(\ell)\bigr]
\le\sum_{\ell<\ell_0}\frac{36}{(n-2)\,\theta_0(\ell)}
\]
(with $S'$ and $\theta_{S'}$ in the adaptive case), which is
$36\ell_0/((n-2)\theta_0')=O(1/\log^2n)$ for Hypothesis~\ref{hyp:orbitgen}, and
for Hypothesis~\ref{hyp:adaptive}, with $m_\ell\le n/(\log^2n\min(\ell-1,M))+1$,
\begin{multline*}
\sum_{\ell<\ell_0}\frac{144\,m_\ell}{(n-2)\kappa_2\log^2n}
\le\frac{144}{\kappa_2\log^4n}\Bigl(\sum_{\ell-1\le M}\frac1{\ell-1}+\sum_{M<\ell-1<\ell_0}\frac1M\Bigr)+\frac{144\,\ell_0}{(n-2)\kappa_2\log^2n}\\
=O\Bigl(\frac{\log n}{\log^4n}+\frac{\ell_0}{M\log^4n}+\frac1{\log^4n}\Bigr)=O\Bigl(\frac1{\log^2n}\Bigr).
\end{multline*}
The exceptional events of (ii) over all $\ell$ form a single event of probability
$o(1/\log n)$ ($\ell=d_{12}$ is a frame statistic).  Altogether
$\Pr[\mathrm{Flip}=\emptyset]=o(1/\log n)$.
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
\qquad\text{for }s=\frac n{16\log^2n},\ T=\frac{\log^2n}4,\ \ell_1=\frac{\log^2n}2,
\]
against the $C\log\log n$ ($C>1/\log2$) that (i) needs; the constant is small
($\approx1/145$), so the heuristic excess over the threshold is a large-$n$ statement.  Compared with the fixed matching
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

\paragraph{Comparison with the witness lemma.}  Hypothesis~\ref{hyp:witness} asks for
\emph{one} witness on a long ladder (a short $\rho_{x_k,y_k}$-cycle through $p$ avoiding
$p'$, which makes the pair flippable with no further argument), at a per-pair rate
$\approx\ell'/n$.  Hypothesis~\ref{hyp:adaptive} asks for $C\log\log n$ good candidates at
a per-pair rate $\approx T\ell_1/n$, and Hypothesis~\ref{hyp:orbitgen} asks for a
genericity constant.  Neither orbit hypothesis is implied by the witness lemma or implies
it (we see no implication in either direction); the witness lemma is the least demanding
heuristically and is the reference hypothesis of this paper.  What the orbit method contributes is the mechanism ---
exact fairness compounded over commuting involutions into exponential decay --- and, with
the adaptive coordinates, the observation that the mechanism needs only rare local
witnesses of a second kind (column cycles of a free pair through $x_k$) rather than a
constant.  All three hypotheses share the same unproved core, stated in
\S\ref{sec:remains}.

\paragraph{(d) Data.}  All scripts and logs are in the repository \cite{repo}, under
\texttt{s13\_*.py} and \texttt{runs/s13/orbit/}.  \emph{Adaptive coordinates.}  On
instances with a uniformly random mark at $n=30,50,100$ (family the first $\min(K_0-1,n/8)$ ladder pairs,
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

Hypotheses~\ref{hyp:witness}, \ref{hyp:orbitgen}, \ref{hyp:adaptive}
and~\ref{hyp:firstpair} are statements about configurations at positions selected by the
frame: a short $\rho_{x_k,y_k}$-cycle through $p$ avoiding $p'$; a clean column pair
separated by the mark; a short clean column cycle at a prescribed column pair; a column pair
along $\rho_{x_1,y_1}$ whose column permutation separates $x_1$ from $y_1$.  Their heuristic
frequencies are, for large $n$, far above what is needed; the difficulty is entirely that
the positions are selected by the frame of the mark and the whole space is conditioned on
the type of rows $1,2$.  The minimal instance of the obstacle is a lower bound on a single
frame-selected pair:
\begin{equation}\label{eq:minimal}
\Pr_X\bigl[(x_1,y_1)\text{ flippable}\bigr]\ge c,\qquad x_1=\pi^{-1}(1),\ y_1=\pi^{-1}(2),
\end{equation}
which no identity of this paper gives unconditionally: flippability of a ladder pair is a
class invariant of the ladder trades, and Proposition~\ref{prop:firstpair} reduces
\eqref{eq:minimal} one way to $\Pr_X[G]\ge2c$, a lower bound for an event of the same kind
with rows and columns exchanged (the column permutation of a pair of columns selected from
rows $x_1,y_1$ separating the frame-selected rows $x_1,y_1$), which is not known to be
easier.  By row exchangeability \eqref{eq:minimal} is a statement at fixed positions under
conditioning: it is the probability that $p,p'$ lie in different cycles of $\rho_{3,4}$
given rows $1,2$ and the two cells $L(3,p)=L(1,p')$, $L(4,p)=L(2,p')$.  Even without the
conditioning on rows $1,2$ we know of no lower bound in the literature for the related
quantity $\Pr[\text{the }\rho_{3,4}\text{-cycle through }p\text{ has length }k\text{ and avoids }p']$,
$3\le k\le\log^5n$; the one case we know to be settled is $k=2$, where the concentration
of the number of intercalates at $n^2/4$ \cite{KSS21} and the exchangeability of row
pairs and of columns give $\approx1/n$.  We record what
the available tools do and do not give.

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
construction, and the price is paid in rarity rather than in a constant we cannot prove;
at the first ladder pair alone (\S\ref{sec:orbit}(e)) the rarity disappears as well, and
what is left is an upper bound on a short arc at the frame-selected pair together with the
lower tail of the number of good candidates among the $\nu-1$ column pairs along
$\rho_{x_1,y_1}$ --- the smallest sufficient condition we have, and still a lower bound at
square-selected positions.
"""
