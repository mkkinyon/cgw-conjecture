"""Paper section 'The hazard term: an exact identity, the first candidate, and fixed positions' (sec:hazardid).
Source: notes sections (r), (s), (t) of section_ladder_s13.tex, as corrected after the referees and cold read 3 (2026-10-04)."""

hazard_block = r"""
\section{The hazard term: an exact identity, the first candidate, and fixed positions}\label{sec:hazardid}

This section records what we can prove about Hypothesis~\ref{hyp:hazard} with the tools of
\S\ref{sec:fc}: an exact switching identity for the conditional hazard, which locates the
difficulty on both of its sides; a lower bound, at scale $n$, on the length of the column
cycle and of its arcs at the first candidate, conditionally on $x_1,y_1$ sharing a cycle there
(the joined part of ``bad'');
and the statement that the same machinery \emph{does} give at fixed positions, which shows
exactly what the selection of the candidate columns by the square costs.  None of it
proves the hypothesis at any candidate, and \S\ref{sec:hazardid}(d) says why the first
candidate is not reached.  Throughout, $R$ and the mark are fixed and all probabilities are
in $S(R)$; the statements for $\Pr_X$ follow by averaging.  The special rows are the frame rows $1,2$ and
$x_1,y_1$; a row is \emph{free} if it is not special; $\tau_t$, $Z_t(\cdot)$ and ``good'' are as in
\S\ref{sec:orbit}(e), and $D_t=Z_t(x_1)$, its rows listed cyclically reading from column $q_t$ as in
\S\ref{sec:fc}.

\subsection*{(a) The switching identity}

For $j\ge1$ let $E_j=\{\nu>j,\ \text{candidates }1,\dots,j\text{ bad}\}$ and
$F_j=\{\nu>j,\ 1,\dots,j-1\text{ bad},\ j\text{ good}\}$, and let
$h_j=\Pr[F_j]/(\Pr[E_j]+\Pr[F_j])$ be the conditional hazard at $j$.  A bad candidate is
\emph{joined} if $y_1\in D_j$ and \emph{separated-inadmissible} otherwise.  For $L\in E_j$
with $j$ joined, a \emph{toggle} is a row switch $\rho_L(z,z',c_j)$ of two free rows
$z,z'\in D_j$ which interleave $x_1,y_1$ in the cyclic order of $D_j$ (exactly one of
$x_1,y_1$ lies strictly between them) and whose footprint, the row cycle of $(z,z')$
through $c_j$, avoids $q_j$.  By \cite[Lemma~2.4(1)]{AM26} it separates $x_1$ from $y_1$
on the $(q_j,c_j)$-cycles; it changes only rows $z,z'$, so it stays in $S(R)$ and leaves
$\rho_{x_1,y_1}$, $\nu$ and every candidate unchanged.  It is \emph{ok} if candidate $j$ is good afterwards (separated and
admissible) and \emph{preserving} if moreover the result lies in $F_j$ (candidates
$1,\dots,j-1$ still bad).  For $L'\in F_j$ a \emph{merge} is a row switch
$\rho_{L'}(z,z',c_j)$ with $z$ a free row of $Z_j(x_1)$ and $z'$ a free row of
$Z_j(y_1)$ whose footprint avoids $q_j$; it joins the two cycles with $z,z'$ at the
junctions (\cite[Lemma~2.4(2)]{AM26}), and it is \emph{preserving} if the result lies in
$E_j$.  Let $f_j(L)$ and $b_j(L')$ be the numbers of preserving toggles and merges ($f_j=0$ when
$j$ is separated-inadmissible).

\begin{proposition}[the hazard identity]\label{prop:hazardid}
For every $j\ge1$ with $\Pr[E_j]\Pr[F_j]>0$,
\[
\frac{\Pr[E_j]}{\Pr[F_j]}=\frac{\mathbb E[\,b_j\mid F_j]}{\mathbb E[\,f_j\mid E_j]},
\qquad b_j\le|Z_j(x_1)|\,|Z_j(y_1)|\le n^2/4 .
\]
In particular $h_j\ge c$ iff $\mathbb E[b_j\mid F_j]\le\frac{1-c}{c}\,\mathbb E[f_j\mid E_j]$.
\end{proposition}

\begin{proof}
A preserving toggle of $L$ at $(z,z')$ is undone by the row switch $\rho_{L'}(z,z',c_j)$
(\cite[Fact~2.1]{AM26}: the same set of cells), which is a preserving merge of $L'$ with
$z\in Z_j(x_1)$, $z'\in Z_j(y_1)$; conversely a preserving merge of $L'$ at $(z,z')$ is
undone by a preserving toggle of its result at the same pair (the merged cycle has $z,z'$
at the junctions, hence interleaving, and the footprint is the same row cycle, avoiding
$q_j$; a switch whose footprint contained $q_j$ would conjugate the column permutation and
change no candidate's status).  So the pairs $(L,\text{preserving toggle})$ and
$(L',\text{preserving merge})$ are in bijection, and double counting gives
$|E_j|\,\mathbb E[f_j\mid E_j]=|F_j|\,\mathbb E[b_j\mid F_j]$.
\end{proof}

\begin{remark}[what the identity shows]\label{rem:hazardid}
Both sides are counts of local moves inside $S(R)$, so the hazard at $j$ is a comparison
of two cycle statistics in two conditioned ensembles.  With the trivial bound
$b_j\le n^2/4$ one needs $\mathbb E[f_j\mid E_j]\ge c'n^2$: the arcs of $D_j$ between
$x_1$ and $y_1$ must have length $\Theta(n)$, and a positive density of the pairs across
them must be good and preserving.  In the data (Remark~\ref{rem:fixrect}; \cite{repo},
\texttt{s13\_hazard.py}) the preserving fraction decays like $(1-\varepsilon)^{j-1}$ with
$\varepsilon\approx0.05$ at $n\le100$ (an ok toggle flips an earlier candidate $s$ when
its footprint, of length $\approx0.35n$, meets exactly one of $q_s,c_s$ and the pair
interleaves $x_1,y_1$ on $D_s$), so that the bound the chain of conditional hazards gives on
$\Pr[E_M]$ with the trivial backward bound stays above a positive constant for every $M$
(if $\varepsilon(n)$ stays bounded below, which the data do not decide); the identity forces
$\mathbb E[b_j\mid F_j]$ to decay at the same rate, and it does (the measured ratio of
the two sides is $\approx1.5$ at $j\le6$, matching $\Pr[E_j]/\Pr[F_j]$ for $j\le4$).  So what fails is
not the move family but the two obvious bounds on its two sides; a direct comparison of
the two conditional expectations is not excluded, and we do not have one.  The data are
for the type $(n)$ only.
\end{remark}

\subsection*{(b) Lengths at the first candidate}

Let $E_1^{\rm j}=\{x_1\sim y_1\text{ on }D_1\}$ (the joined part of $E_1$; $\nu\ge2$
always by Lemma~\ref{lem:mu2}), $\ell=|D_1|$, and $a,a'$ with $a+a'=\ell$ the numbers of
rows of the two arcs of $D_1$ between $x_1$ and $y_1$ (each arc counted with one of its
endpoints).  A $(q_1,c_1)$-cycle is \emph{inactive} if it contains no special row, and
$g(L)$ is the number of rows on inactive cycles.

\begin{proposition}[cycle and arc lengths at the first candidate]\label{prop:lengths}
Let $n\ge\Delta^{32}$, $8\le k\le n/100$ and $1/n\le\delta\le1/100$.  Then
\[
\Pr\bigl[\ell\in[k,2k]\ \big|\ E_1^{\rm j}\bigr]\le400\,\frac kn,
\qquad
\Pr\bigl[\min(a,a')\in[k,2k],\ \ell\ge\delta n\ \big|\ E_1^{\rm j}\bigr]\le6000\,\frac{k}{\delta^2n},
\]
and consequently, with $k_0=\lfloor\log n/(4\log\Delta)\rfloor$,
\[
\Pr\bigl[\ell\le\delta n\text{ or }\min(a,a')\le\delta^3n\ \big|\ E_1^{\rm j}\bigr]
\le48000\,\delta+\frac{5(\Delta^5+k_0\Delta^2)}{\sqrt n\,\Pr[E_1^{\rm j}]}.
\]
\end{proposition}

The mechanism is that of Claim~$4'$: $\Omega(kn)$ moves out of the class, $O(k^2)$ into
any square, the inverse pinned by the two marked rows $x_1,y_1$, which must both lie in a
short piece; the dyadic sum $\sum_kk/n$ over $k\le\delta n$ is $O(\delta)$, which is the
right order ($\Pr[\ell\le\delta n]\approx\delta$ for a random permutation).  Nothing here
bounds a specific cycle from below; the proposition says that a cycle through two marked
rows is unlikely to be of intermediate length.  The constants are not optimised and the
statement is purely asymptotic ($n\ge\Delta^{32}$ comes from matching the dyadic range with
the fixed-cell range at $k_0\ge8$).  If $\Pr[E_1^{\rm j}]$ is small the first-candidate
hazard is large anyway, so the denominator is harmless.

\begin{proof}
\emph{(i) The cycle length.}  Fix $k$ and let $\mathcal C_k=E_1^{\rm j}\cap\{\ell\in[k,2k]\}$.
\emph{Moves.}  For $L\in\mathcal C_k$, a free row $z\in D_1$ and a row $z'$ on an inactive
cycle $C$: if $\rho_L(z,z',c_1)$ avoids $q_1$, switch it; otherwise first switch the column
cycle $C$ (columns $q_1,c_1$, rows of $C$) to obtain $L''$, in which $\rho_{L''}(z,z',c_1)$
avoids $q_1$ (\cite[Lemma~2.3(1)]{AM26} after re-basing the row cycle at a third of its
columns, which exist since $z,z'$ lie on different $(q_1,c_1)$-cycles; directly as in the
proof of Claim~$4'$, since $C\not\ni z$), and then switch it.  The column switch changes
only rows of $C$, which are not special; the row switch only rows $z,z'$.  So rows
$1,2,x_1,y_1$ are unchanged and the result $L'$ is in $S(R)$ with the same mark, first
pair, $\rho_{x_1,y_1}$, $\nu$ and candidates.  By \cite[Lemma~2.4(2)]{AM26} $D_1$ and $C$
merge into one cycle $D'=D_1\cup C$ and every other cycle is unchanged; directly, column
$q_1$ is unchanged and the $c_1$-entries of $z,z'$ are exchanged, so the column permutation
is composed with the transposition of $L[z,c_1],L[z',c_1]$, which splices $C$ into $D_1$
immediately after $z$ in the list.  Hence $L'\in E_1^{\rm j}$ with $\ell'=\ell+|C|\ge\ell+2$
and the cyclic order of the special rows on the cycle preserved.  There are at least
$(\ell-4)\,g(L)$ moves from $L$.
\emph{Inverse count.}  If $L'$ arises from $(L,z,z')$ then $L$ is recovered from $(L',z,z')$
by switching $\rho_{L'}(z,z',c_1)$ (\cite[Fact~2.1]{AM26}) and then, in the two-step case,
switching the $(q_1,c_1)$-cycle through $z'$; so $L'$ has at most $2N$ preimages, $N$ the
number of unordered pairs $\{z,z'\}$ of free rows on $D'$ such that switching
$\rho_{L'}(z,z',c_1)$ splits $D'$ into a piece $P\ni x_1,y_1$ with $|P|\in[k,2k]$ and a
piece $S$.  By \cite[Lemma~2.4(1)]{AM26} $S$ consists of the rows strictly between $z$ and
$z'$ in the list together with $z'$, so $P=D'\setminus S$ is a contiguous arc of $D'$
ending at $z$, which determines $\{z,z'\}$; and a contiguous arc of length $m$ containing
$x_1$ can be placed in at most $m$ ways.  Hence $N\le\sum_{m=k}^{2k}m\le3k^2$ and $L'$ has
at most $6k^2$ preimages.
\emph{Lower bound on $g$.}  The proof of Claim~$3'$ applies verbatim with special set
$\{1,2,x_1,y_1\}$, $Q=\emptyset$ and the class $\mathcal C_k$: its moves (row switches of
good same-arc pairs on the cycles through rows $1,2$ other than $D_1$, and the
cross-switches of Claim~$2'$ on those cycles) touch neither $D_1$ nor a special row, so
$\ell$, $E_1^{\rm j}$, the mark and the candidates are preserved.  The free rows off $D_1$
number at least $n-2k-4$ and lie on the switched active cycles or on inactive cycles; the
forward count is at least $\frac1{27}((n-2k-g-13)_+)^2$ in aggregate, the backward count at
most $n\,g(L')$, and Jensen's inequality gives $\mathbb E[g\mid\mathcal C_k]\ge(n-2k)/31$
for $k\le n/100$ and $n$ large.
\emph{Conclusion.}  $|\mathcal C_k|\,(k-4)\,\mathbb E[g\mid\mathcal C_k]\le6k^2\,|E_1^{\rm j}|$,
so $\Pr[\ell\in[k,2k]\mid E_1^{\rm j}]\le6k^2\cdot31/((k-4)(n-2k))\le400k/n$ for $8\le k\le n/100$.

\emph{(ii) The arcs.}  Let $\mathcal C'_k=E_1^{\rm j}\cap\{a\in[k,2k],\ \ell\ge\delta n\}$,
$a$ the arc from $x_1$ to $y_1$ in the list order (the other arc by symmetry).  If
$k>\delta n/4$ or $\delta n<240$ the bound is $\ge1$; assume otherwise.  The other arc has
$\ell-a\ge\delta n/2$ rows and at most two frame rows, so it contains a special-row-free
stretch of at least $\delta n/6-1$ free rows; let $B=B(L)$ be the longest (ties by
position), $b=|B|$.  \emph{Transfer moves.}  For a good same-arc pair $(w,w')$ in $B$
(footprint avoiding $q_1$) and a free row $z$ of the $x_1$--$y_1$ arc (at least $k-3$
choices): switch $\rho_L(w,w',c_1)$, which by \cite[Lemma~2.4(1)]{AM26} cuts the rows
strictly between $w$ and $w'$ together with $w'$ off as a new cycle $S$ containing no
special row, hence inactive, leaving $x_1,y_1$ and their arc untouched; then perform the
insert move of (i) at the pair $(z,z')$ with $z':=w'$.  The result is in $E_1^{\rm j}$ with
$\ell'=\ell$, $a'=a+|S|\ge a+2$, the same cyclic order of special rows, and $S$ spliced in
after $z$ and ending at $w'$ (in the two-step case the column switch reverses $S$ as a
cycle and the splice still ends at $z'=w'$).  Claim~$2'$, applied to the class of triples
$(L,D_1,B(L))$ with $L\in\mathcal C'_k$ and $|B(L)|=b$ --- closed under its cross-switches,
which act inside $B$, touch no special row, keep $B$ as a set (\cite[Remark~2.7]{AM26})
and change neither $\ell$ nor $a$ --- gives at least $(b-3)^2/9\ge\delta^2n^2/400$ good
pairs in $B$ in aggregate, hence at least $(k-3)\delta^2n^2/400$ moves per square of
$\mathcal C'_k$ on average.  \emph{Inverse count.}  From $L'$: the segment $S$ to be cut
out of the $x_1$--$y_1$ arc leaves a prefix and a suffix of that arc of total length in
$[k,2k]$, at most $(2k)^2$ choices, and determines $(z,z')$ ($z$ the row before $S$,
$z'=w'$ its last row); undoing the insert is the row switch at $(z,z')$ and, in the
two-step case, the column switch of the cycle through $z'$ (two alternatives); undoing the
cut is the row switch at $(w,w')$ with $w'$ known and $w$ one of at most $n$ rows.  So at
most $8k^2n$ preimages, and
$\Pr[a\in[k,2k],\ \ell\ge\delta n\mid E_1^{\rm j}]\le8k^2n\cdot400/((k-3)\delta^2n^2)\le6000k/(\delta^2n)$
for $k\ge8$.

\emph{(iii) Short lengths.}  $\{\ell=m\}\cap E_1^{\rm j}$ is the union, over the positions
of $x_1,y_1$ ($n^2$), of $q_1,c_1$ ($n^2$), of the $m-2$ other rows of the cycle
($n^{m-2}$), of the position of $y_1$ on it ($m-1$) and of the $m-1$ free symbols
($n^{m-1}$; $L[y_1,q_1]=L[1,p']$ is forced), of fixed-cell events with $2m+2$ cells (the
cycle and the anchors $(x_1,p),(y_1,p)$); by Corollary~\ref{cor:FCtwo} and the
exchangeability of rows $3,\dots,n$, as in Corollary~\ref{cor:nutail}, its probability is
at most $(m-1)\Delta^{2m+2}/n$.  For $\{a=m\}\cap E_1^{\rm j}$ the path from $x_1$ to $y_1$
has $m+1$ rows and $2m+2$ cells, and $c_1$ is pinned by the further cells $(x_1,p')$ and
$(y_1,c_1)$: $2m+5$ cells, $n^{2m+4}$ placements, probability at most $\Delta^{2m+5}/n$.
Summing over $m\le k_0$ (at most $2k_0\Delta^2/\sqrt n$ and $4\Delta^5/\sqrt n$, using
$\Delta^{2k_0}\le\sqrt n$), (i) over the dyadic $k\in[k_0,\delta n]$ ($\le800\delta$) and
(ii) over the dyadic $k\in[k_0,\delta^3n]$ for each arc ($\le2\cdot12000\delta$) gives
the consequence.
\end{proof}

\begin{remark}
The same argument, with the frame rows (when they lie on $D_1$) as the pinning rows,
bounds the lengths of the stretches of $D_1$ between consecutive special rows below by
$\delta^3n$ with
probability $1-O(\delta)$; and it applies to $D_j$ for any $j<\nu$ as a statement about
$\Pr[\,\cdot\mid x_1\sim y_1\text{ at }j]$, but not about $\Pr[\,\cdot\mid E_j]$, since its
row switches have footprints of length $\Theta(n)$ that change the status of earlier
candidates.  On sampled completions at $n=30,50,100$ every sampled insert and transfer
move behaved as stated and was undone by the stated inverse, the preimage counts $N$ were
at most $3k^2$, $\mathbb E[g/(n-\ell)]$ was at least $0.26$ in every bin of $\ell/n$
(bound $1/31$), and $\Pr[\ell\in[k,2k)\mid E_1^{\rm j}]$ tracked the permutation value
$3k^2/n^2$ \cite{repo}.
\end{remark}

\subsection*{(c) Fixed positions}

\begin{proposition}[two fixed rows at a fixed column pair]\label{prop:fixedpos}
Let $R$ be a $2\times n$ Latin rectangle (rows $1,2$), $n\ge504$, $q\ne c$ columns and
$x\ne y$ rows in $\{3,\dots,n\}$.  Then
\[
\Pr_{S(R)}[x\text{ and }y\text{ lie on the same }(q,c)\text{-column cycle}]\ \le\ 1-\tfrac1{30}+\tfrac3n .
\]
For a uniform Latin square of order $n\ge502$ and arbitrary rows $x\ne y$ the bound is
$1-\tfrac1{30}$.
\end{proposition}

\begin{proof}
Let $D_x$ be the $(q,c)$-cycle through $x$ and $\ell_x=|D_x|$.  For a permutation $\tau$ of
$\{3,\dots,n\}$ fixing $x$, $L\mapsto\tau L$ is a bijection of $S(R)$ with $y\in D_x(L)$
iff $\tau y\in D_x(\tau L)$; averaging over $\tau$,
$\Pr[y\in D_x]=\mathbb E\,|D_x\cap\{3,\dots,n\}\setminus\{x\}|/(n-3)\le\mathbb E[\ell_x-1]/(n-3)$.
The proof of Claim~$3'$ applies verbatim with special set $\{1,2,x\}$, $Q=\emptyset$ and
the class $S(R)$ (at most three arcs; none of its moves touches a special row), and gives
$\mathbb E[g]\ge f/30$ for the number $g$ of rows on cycles containing no special row,
$f=n-3\ge501$.  Since $\ell_x\le n-g$, $\mathbb E[\ell_x-1]\le n-1-(n-3)/30$.  For a
uniform square take the special set $\{x\}$, $f=n-1$, and average over all $\tau$ fixing
$x$.
\end{proof}

\begin{corollary}\label{cor:longest}
For a uniform Latin square and a fixed column pair, with $M$ the length of the longest
cycle of the pair, $\mathbb E[n-M]\ge(n-1)/60$ and
$\Pr[n-M\ge n/120]\ge1/120-O(1/n)$; by transposition the same holds for the longest cycle
of a fixed row pair.
\end{corollary}

\begin{proof}
$\sum_C|C|(n-|C|)=\sum_x(n-\ell_x)\le2n(n-M)$, and $\mathbb E\sum_x(n-\ell_x)\ge n(n-1)/30$
by the proposition.
\end{proof}

This is a statement in the ``coarse'' direction that the surviving material of
\cite{CGW08} does not give (\cite{gapnote}, \S4): with probability bounded below, the rows
off the longest cycle of a column pair are a positive fraction of all rows.  It is not an
improvement for the type $(n)$ itself, for which the intercalate argument of
\cite{gapnote} gives $5/6+o(1)$, and it does not approach $\Pr[(n)]=o(1)$.

\subsection*{(d) What the first candidate still needs}

Proposition~\ref{prop:fixedpos} is the first-candidate statement with $(q,c)$ and $y$
\emph{fixed}: it needs no pair across any cut, because $y$ is exchangeable with the other
free rows and ``is $y$ on $D_x$'' becomes ``how long is $D_x$'', which Claim~$3'$
answers.  In the ladder, $y_1$ is selected by the cell $(y_1,p)$ and the candidate
columns $q_1$, $c_1$ --- the columns of row $y_1$ holding $L(1,p')$ and $L(x_1,p')$ --- are
functions of row $y_1$; conditioning on the columns destroys the exchangeability of
$y_1$, and conditioning on less leaves the columns random.  So by
Proposition~\ref{prop:hazardid} with $j=1$ and $b_1\le n^2/4$, the first-candidate hazard
$h_1\ge c$ follows from $\mathbb E[f_1\mid E_1]\ge c'n^2$, and
Proposition~\ref{prop:lengths} with the remark after it supplies the lengths ($\ge\delta^6n^2$ interleaving pairs whose split is admissible, with
probability $1-O(\delta)$ on $E_1^{\rm j}$) but not the goodness: a positive density of the
pairs across the $x_1$--$y_1$ cut must have footprint avoiding $q_1$.  This does not
follow from Claim~$2'$: the cross-switch at a bad interleaving pair reverses the segment
between the two rows containing one of $x_1,y_1$ and no frame row (legal, and it preserves
$\nu$ and every candidate as an unordered pair), and by \cite[Lemma~2.6]{AM26} it toggles the goodness of
every pair with exactly one row in that segment --- but it toggles their interleaving
status as well, since the inner row and the inner special row exchange their order; so the
Kind-2 image of a bad interleaving pair is a good \emph{non-interleaving} pair, outside the
family, and with same-stretch auxiliary pairs the Kind-1 image leaves the family instead.
More than that --- and this whole paragraph is heuristic, not a theorem --- the
configuration in which a pair of free rows of $D_1$ at distance $\ge2$
is good iff it does not interleave $x_1,y_1$ (adjacent pairs $\{z,\tau_1(z)\}$ are always
bad), with inactive cycles carrying an orientation-dependent affinity to one arc, is
preserved, up to the $O(n)$ pairs containing a switched row, by every move used in this
section and in \S\ref{sec:fc} --- cross-switches, inserts, cuts, transfers and column
switches of inactive cycles --- checked move by move, not proved; and on it $f_1=0$.  We
therefore see no counting argument from these moves that gives the goodness across the
cut.
(The truth is the opposite: in the data the good fraction is $0.50$ among interleaving and
among non-adjacent same-arc pairs alike, at every distance $\ge2$ along $D_1$, with
adjacent pairs never good \cite{repo}.)  The separated-inadmissible part of $E_1$ ($15$--$17\%$
of the bad candidates in the data, pooled over $j$), on which $f_1=0$, needs a toggle that moves a
frame row off the cycle of $x_1$, which we have not analysed.  So the obstacle of
\S\ref{sec:remains} has, at the first candidate, the exact form: a positive density of
good pairs across the $x_1$--$y_1$ cut, without the exchangeability of $y_1$; and at later
candidates, in addition, the preservation of the earlier ones.
"""
