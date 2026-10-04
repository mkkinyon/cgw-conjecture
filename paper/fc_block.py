"""fc_block.py -- paper section: the Allsop-Morris per-cell bound with two complete rows, (FC) in S(R), the nu-tail,
and the hazard hypothesis.  Source: notes section (q) (two referee passes).  Imported by assemble_conditional.py."""

fc_block = r"""

\section{The fixed-cell bound in completions of a rectangle}\label{sec:fc}

Throughout this section $R$ is a fixed $2\times n$ Latin rectangle occupying rows $1,2$,
$Q$ is a partial Latin square with cells in rows other than $1,2$ such that $R\cup Q$ is a
partial Latin square, and $L$ is uniform on $\{L\in\mathcal L_n:L\supseteq R\cup Q\}$; we
write $S(R)$ for the completions of $R$.  Nothing in this section moves rows $1,2$, and
the mark plays no role.  (In this section $\alpha,\beta$ denote the density parameters of
\cite{AM26}, not the mark parameters; the mark parameters return in
Corollary~\ref{cor:nutail}.  Their statistic $\nu(L)$ is not used; our $\nu$ is the arc
minimum of \S\ref{sec:orbit}(e).)  We use the switchings of Allsop and Morris \cite{AM26} and follow
their notation (\S\S2--4 there): $\rho_L(r,r',c)$ is the row cycle of rows $r,r'$ through
column $c$, $\gamma_L(c,c',r)$ the column cycle of columns $c,c'$ through row $r$,
$\zeta_L(c,c',r,r')$ the cross-switch, and $\eta_L(r,r',c)$ the $\eta$-switch, which
replaces the entries $(r,c,s),(r',c,s'),(r,c',s')$ --- where $s=L[r,c]$, $s'=L[r',c]$ and
$c'$ is the column with $L[r,c']=s'$ --- by $(r,c,s'),(r',c,s),(r,c',s)$ and otherwise
changes only rows $r'$ and $r''$, the row with $L[r'',c']=s$; it is of \emph{Type One} if
$r''=r'$ (then $\rho_L(r,r',c)$ is an intercalate), and by their Lemma~2.8 a square has
at most $n$ preimages under an $\eta$-switch of Type One and at most one otherwise.  Their
Theorem~1.1 bounds $\Pr[L\supseteq P]$ for partial Latin squares $P$ occupying at most
$\alpha n$ rows and $\beta n$ columns with $2\alpha+\beta<1$, and so does not apply to
$P=R\cup Q$, which occupies every column.  Their proof does, with changes we describe in
full.  Call rows $1,2$ \emph{special} and write $E=\{1,2\}\cup R_Q$.

\begin{theorem}[\cite{AM26}, Theorem~3.1, with two complete rows]\label{thm:AM31two}
Let $R,Q$ be as above with $|R_Q|+6\le\alpha n$, $|C_Q|\le\beta n$ and $2\alpha+\beta<1$.
Let $r,r'\notin\{1,2\}$ be distinct rows and $c$ a column, and suppose that $Q$ has no
cell in row $r'$ outside column $c$.  Then, for $n\ge n_0(\alpha,\beta)$,
\[
\Pr\bigl[\eta_L(r,r',c)\text{ is of Type One}\;\big|\;L\supseteq R\cup Q\bigr]\le\frac{D}{n},
\qquad D=D(\alpha,\beta)\to62\text{ as }\alpha,\beta\to0 .
\]
\end{theorem}

The constant is not optimised.  The proof is that of \cite{AM26} with three changes: the
row-exclusion set $R_P$ becomes $E$ and the column-exclusion set $C_P$ becomes $C_Q$ (the
cells of $R$ need no exclusion because no switch below changes rows $1,2$); the
cross-switches of their Claim~2 and the row switches of their Claim~3 are restricted to
pairs of rows in a common \emph{arc} of a column cycle, so that no column part of a switch
passes through a special row; and the statistic $\nu(L)$ of their Claims~3--4 is replaced by
one that also counts the rows on the $(c_2,c_3)$-cycles through the special rows, so that the
column-cycle switch of their Claim~4 is only ever applied to a cycle avoiding rows $1,2$.

\paragraph{The recursion and Claim~1.}  For $L\supseteq R\cup Q$ put $r_0=r$, $c_0=c$,
$s_0=L[r_0,c_0]$, $r_1=r'$, and let $c_1$ be the column with $L[r_0,c_1]=L[r_1,c_0]$.
Recursively, for $j\ge1$, $s_j=L[r_j,c_j]$; $r_{j+1}$ is the row with $L[r_{j+1},c_j]=s_{j-1}$
and $c_{j+1}$ the column with $L[r_j,c_{j+1}]=s_{j-1}$; $X_j$ is the set of squares with
$s_j=s_{j-1}$, and $Y_j$ the set with $s_j\ne s_{j-1}$, $r_{j+1}\notin E\cup\{r_\ell:\ell\le j\}$
and $c_{j+1}\notin C_Q\cup\{c_\ell:\ell\le j\}$ (with $Y_0$ all squares); $p_{j-1}=\Pr[X_j\mid Y_{j-1}]$,
so that $p_0$ is the probability to be bounded.  Their Claim~1,
\[
p_{j-1}\le\frac{np_j+1}{n(1+p_j-2\alpha-\beta)-3j-2}\qquad(1\le j<n(1-2\alpha-\beta)/3),
\]
holds verbatim.  Its switch $\eta_L(r_j,x,c_j)$, for $L\in X_j$ and a row $x\notin E\cup\{r_\ell\}$
with $y$ (the column with $L[r_j,y]=L[x,c_j]$) outside $C_Q\cup\{c_\ell\}$ and $x'$ (the row
with $L[x',y]=s_{j-1}$) outside $E\cup\{r_\ell\}$, changes only rows $r_j,x,x'$, all outside
$\{1,2\}$, so rows $1,2$ are untouched; the cells of $Q$ it could touch are $(r_j,y)$,
excluded, $(r_j,c_j)$, which is not a cell of $Q$ ($c_j\notin C_Q$ for $j\ge2$, and for $j=1$
the only cell of $Q$ allowed in row $r_1$ is in column $c_0\ne c_1$), and cells in rows
$x,x'$, excluded.  The maps $x\mapsto y\mapsto x'$ are injective, so at most $|C_Q|+j+1$
choices of $x$ are lost to $y$ and at most $|E|+j+1$ to $x'$, giving at least
$n-2|E|-|C_Q|-3j-3\ge n(1-2\alpha-\beta)-3j-3$ switches from each $L\in X_j$ into $Y_j$,
and the reversibility count is theirs.

\paragraph{Setting for Claims~2--4.}  For $L\in X_3$ let $Q(L)$ consist of the cells of $Q$
and the ten structure cells $\{(r_\ell,c_\ell),(r_{\ell+1},c_\ell),(r_\ell,c_{\ell+1}):\ell\le2\}\cup\{(r_3,c_3)\}$;
a row is \emph{free} if it is not in $E\cup\{r_0,r_1,r_2,r_3\}$, and $f\ge n-|E|-4\ge n(1-\alpha)$
is the number of free rows.  (The row $r_0$ need not be free and may lie anywhere; its cells
in columns $c_2,c_3$ are not in $Q(L)$ and are never read by the recursion, which reads
only the ten structure cells.)  Consider the partition of the rows into $(c_2,c_3)$-column
cycles.  As in \cite{AM26}: $\gamma_L=\gamma_L(c_2,c_3,r_1)$ is one of them; $\{r_2,r_3\}$ is
another (an intercalate, since $s_3=s_2$); and $c_2,c_3\notin C_Q$.  Call a cycle
\emph{active} if it is $\gamma_L$ or contains a special row (so there are one to three
active cycles, and $\{r_2,r_3\}$ is never active), let $\mathcal S(L)$ be the set of free
rows on active cycles, and $g(L)=f-|\mathcal S(L)|$ the number of free rows on inactive
cycles.  An inactive cycle other than $\{r_2,r_3\}$ avoids rows $1,2$ and $r_1,r_2,r_3$, so
switching it changes no cell of $Q(L)$ or $R$.

Give each active cycle $D$ (a cycle, not to be confused with the constant $D$) a \emph{reference row}: $r_1$ if $D=\gamma_L$, otherwise row $1$
if $1\in D$, otherwise row $2$.  List the rows of $D$ as $x'_1,\dots,x'_m$ in the cyclic
order from the reference row $x'_1$ read from column $c_2$ (the row after $x$ is the row
whose $c_2$-entry is the $c_3$-entry of $x$, as on p.~16 of \cite{AM26}), and let the free
rows of $D$ in this order be $x_1,\dots,x_d$.  The special rows of $D$ other than the
reference row cut the list into \emph{arcs}; a pair $(i,j)$, $i<j$, of free positions is a
\emph{same-arc pair} if no special row and no reference row lies strictly between $x_i$
and $x_j$.  The free rows of the active cycles are the disjoint union of the free rows of
the arcs, and there are at most three arcs in all (if $\gamma_L$ contains both special rows
it has three arcs and is the only active cycle; if one, it has two and the other special
row's cycle has one; if none, $\gamma_L$ has one and the special rows lie on one cycle with
two arcs or on two cycles with one each).  As in \cite{AM26}, $(L,i,j)$ is \emph{good} if
$\rho_L(x_i,x_j,c_3)$ does not hit column $c_2$, and \emph{bad} otherwise.

\paragraph{Claim~$2'$ (good same-arc pairs, in aggregate).}  For $a\ge4$ let $\mathcal A_a$
be the set of triples $(L,D,A)$ with $L\in X_3$, $D$ an active cycle of $L$ and $A$ an arc
of $D$ with exactly $a$ free rows, and $N_a=|\mathcal A_a|$.  Then
\[
\sum_{(L,D,A)\in\mathcal A_a}\#\{\text{good pairs }(i,j)\text{ with }x_i,x_j\in A\}
\ \ge\ N_a\,\frac{a(a-1)(a-3)}{9(a-2)}\ \ge\ N_a\,\frac{(a-3)^2}{9}.
\]
This is a statement about the class $\mathcal A_a$, as their Claim~2 is about the class
$A_k$; a single arc may have no good pair.

\emph{Proof.}  This is their proof run inside the arcs of size $a$.  For $(i,j)$ in $A$ let
$W_A(i,j)$ be the set of pairs $(u,v)$ in $A$ with exactly one of $u,v$ strictly between $i$
and $j$ and $\{u,v\}\cap\{i,j\}=\emptyset$; then $|W_A(i,j)|=(j-i-1)(a-(j-i)-1)$ in arc
coordinates, $\sum_{(i,j)\subseteq A}|W_A(i,j)|=a(a-1)(a-2)(a-3)/12$ and $|W_A|\le(a-2)^2/4$
(their (5)--(6) with $n-k$ replaced by $a$).  For a bad $(L,i_0,j_0)$ in $A$ and
$(u,v)\in W_A(i_0,j_0)$: if $(L,u,v)$ is good, map to $(L,u,v)$ (Kind~1); if it is bad,
cross-switch on $\zeta_L(c_2,c_3,x_u,x_v)$ (Kind~2).  The cross-switch is well defined
($(L,u,v)$ is bad, and $x_u,x_v$ are not adjacent in the list since a position of
$\{i_0,j_0\}$ separates them).  It changes the rows $x_u,x_v$, which are free, and the
entries in columns $c_2,c_3$ of the rows strictly between $x_u$ and $x_v$ in the list,
which are rows of $D$ that are neither special nor the reference row (the pair is
same-arc) nor $r_2,r_3$ (not on $D$); such rows may include $r_0$ or rows of $R_Q$, whose
cells in columns $c_2,c_3\notin C_Q$ are not in $Q(L)\cup R$.  So it changes no cell of
$R$ or $Q(L)$, and $L'\in X_3$ with the same $r_\ell,c_\ell,s_\ell$.  By their Lemma~2.6 the
row set of $D$ is unchanged, and by their Remark~2.7 the list of $L'$ is that of $L$ with
the segment strictly between $x_u$ and $x_v$ reversed; so the special and reference rows
keep their positions, $A$ is an arc of $L'$ with the same free rows, $(u,v)\in W_A(i',j')$
for the new positions $(i',j')$ of $x_{i_0},x_{j_0}$, and $(L',i',j')$ is good by their
Lemma~2.6(2); hence $(L',D,A)\in\mathcal A_a$.  The inverse counts are theirs: a good
$(L',D',A',i',j')$ has at most $|W_{A'}(i',j')|$ Kind-1 preimages (the interleaving
relation is symmetric) and at most $|W_{A'}(i',j')|$ Kind-2 preimages (one for each
$(u,v)\in W_{A'}(i',j')$, by their Lemma~2.5 and Remark~2.7, the inverse cross-switch
again lying inside $A'$).  Summing over $\mathcal A_a$, $\sum_{\text{bad}}|W|\le2\sum_{\text{good}}|W|$,
whence $N_a\,a(a-1)(a-2)(a-3)/12=\sum_{\text{all}}|W|\le3\sum_{\text{good}}|W|\le3\cdot\#\{\text{good}\}\cdot(a-2)^2/4$.
The last inequality of the claim is $a(a-1)\ge(a-2)(a-3)$.  $\square$

\paragraph{Claim~$3'$ (free rows off the active cycles).}  $\sum_{L\in X_3}g(L)\ge f|X_3|/30$ for $f\ge501$.

\emph{Proof.}  Double count the row switches $\rho_L(x_i,x_j,c_3)$ over $L\in X_3$, active
cycles $D$ of $L$ and good same-arc pairs $(i,j)$ on $D$.  \emph{Forward.}  Such a switch
exchanges the entries of the free rows $x_i,x_j$ on the columns of their row cycle through
$c_3$, which does not contain $c_2$; it changes no cell of $R$ or $Q(L)$, so $L'\in X_3$
with the same $r_\ell,c_\ell,s_\ell$.  The $(c_2,c_3)$-column permutation (the $c_2$-entry
of a row $\mapsto$ its $c_3$-entry) changes only at rows $x_i,x_j$, whose $c_3$-entries are
exchanged: it is composed with the transposition of those two symbols, which lie on one
cycle.  So every cycle other than $D$ is unchanged and $D$ splits into two: one consists
of the rows outside the segment strictly between $x_i$ and $x_j$ together with $x_i$, and
contains the reference row and every special row of $D$; the other, $S$, consists of the
rows strictly between together with $x_j$, and contains no special row because $(i,j)$ is
a same-arc pair (this is their Lemma~2.4(1), with the orientation of our list).  Hence the
active cycles of $L'$ are those of $L$ with $D$ replaced by $D\setminus S$ (if $D=\gamma_L$
then $\gamma_{L'}=D\setminus S\ni r_1$; otherwise $\gamma_{L'}=\gamma_L$), $S$ is inactive in
$L'$, $x_i\in\mathcal S(L')$, $x_j\in S$, and $g(L')=g(L)+|S\cap\text{free}|$.  Let $F(L)$
be the number of such switches from $L$ (distinct pairs of rows give distinct switches,
and a pair of rows lies on one cycle).  Summing Claim~$2'$ over $a$, with
$h(a)=((a-3)_+)^2/9$,
\[
\sum_{L\in X_3}F(L)\ \ge\ \sum_{a\ge4}\ \sum_{(L,D,A)\in\mathcal A_a}\#\{\text{good pairs in }A\}
\ \ge\ \sum_{a\ge4}N_ah(a)\ =\ \sum_{L\in X_3}\ \sum_{A}h(|A|),
\]
and for each $L$, since the arcs partition $\mathcal S(L)$ into at most three parts,
$\sum_Ah(|A|)\ge\frac1{27}(\sum_A(|A|-3)_+)^2\ge\frac1{27}((f-g(L)-9)_+)^2$.
\emph{Backward.}  If $L'$ arises from $(L,D,i,j)$ then $L$ is obtained from $L'$ by
switching $\rho_{L'}(x_i,x_j,c_3)$ (a row cycle switch is an involution), so $L$ is
determined by the unordered pair $\{x_i,x_j\}$, of which one row lies in $\mathcal S(L')$
and the other on an inactive cycle of $L'$; the number of preimages of $L'$ is at most
$|\mathcal S(L')|\,g(L')\le f\,g(L')$.  \emph{Conclusion.}  With
$\bar g=|X_3|^{-1}\sum_Lg(L)$, Jensen's inequality for the convex function
$t\mapsto((f-t-9)_+)^2$ gives
$|X_3|((f-\bar g-9)_+)^2/27\le\sum_L\sum_Ah(|A|)\le\sum_LF(L)\le f\sum_{L'}g(L')=f|X_3|\bar g$.
If $\bar g<f/30$ then $(f-\bar g-9)^2>(29f/30-9)^2\ge0.9f^2>27f\bar g$ for $f\ge501$, a
contradiction.  $\square$

\paragraph{Claim~$4'$.}  $p_2\le60/f\le60/(n(1-\alpha))$ for $f\ge501$.

\emph{Proof.}  Their proof with the set of admissible rows changed.  For $L\in X_3$ let
$x$ be a free row on an inactive cycle ($g(L)$ choices) and $\rho=\rho_L(r_3,x,c_3)$.  If
$\rho$ does not hit column $c_2$, switch on $\rho$ to obtain $L'$; otherwise first switch
on the $(c_2,c_3)$-cycle through $x$ to obtain $L''$, then on $\rho_{L''}(r_3,x,c_3)$ to
obtain $L'$.  The cycle through $x$ is inactive, so it avoids rows $1,2$, $r_1$ (on
$\gamma_L$) and $r_2,r_3$ (the intercalate), and $c_2,c_3\notin C_Q$: switching it changes
no cell of $R$ or $Q(L)$, so $L''\in X_3$ with the same data; and $\rho_{L''}(r_3,x,c_3)$
does not hit $c_2$, because the switch composes the row permutation of $(x,r_3)$ with the
transposition of $L[x,c_2],L[x,c_3]$, two symbols of one cycle, which splits that cycle
with the two symbols, hence the columns $c_2,c_3$, on different cycles (their
Lemma~2.3(1)).  In either case the final switch involves the rows $r_3,x$, both outside
$\{1,2\}\cup R_Q$, along a row cycle avoiding $c_2$; it changes no cell of $R$ or of
$Q(L)\setminus\{(r_3,c_3)\}$ and replaces $s_3=s_2$ by $L[x,c_3]\ne s_2$ (or $L[x,c_2]$),
so $L'\in Y_2\setminus X_3$, as in \cite{AM26}.  This gives at least $\sum_Lg(L)\ge f|X_3|/30$
switches from $X_3$ to $Y_2\setminus X_3$.  Their inverse count applies unchanged: from
$L'$ the row $r''$ with $L'[r'',c_3]=s_2$ equals $x$, and in each of the two cases $L$ is
determined, by switching back $\rho_{L'}(r_3,r'',c_3)$ and, in the second case, then the
$(c_2,c_3)$-cycle of $L''$ through $r''$ (the same set of rows as before); so each $L'$ has
at most two preimages.  Hence $|Y_2\setminus X_3|\ge f|X_3|/60$ and
$p_2=|X_3|/|Y_2|\le1/(1+f/60)\le60/f$.  $\square$

\begin{proof}[Proof of Theorem~\ref{thm:AM31two}]
Claim~1 twice from $p_2\le60/(n(1-\alpha))$, as in (16) of \cite{AM26}: with
$\phi_j(x)=(nx+1)/(n(1+x-2\alpha-\beta)-3j-2)$, which is non-decreasing for $j=o(n)$,
$p_1\le\phi_2(p_2)\le(60/(1-\alpha)+1)/(n(1-2\alpha-\beta)-8)$ and
$p_0\le\phi_1(p_1)\le(np_1+1)/(n(1-2\alpha-\beta)-5)$, so $p_0\le D(\alpha,\beta)/n$ with
$D(\alpha,\beta)=(60/(1-\alpha)+2+o(1))/(1-2\alpha-\beta)^2\to62$.
\end{proof}

\begin{corollary}[per-cell bound and (FC) in $S(R)$]\label{cor:FCtwo}
With $R,Q$ as in Theorem~\ref{thm:AM31two}, for a row $r\notin\{1,2\}$ (in $R_Q$ or not) and
any $(c,s)$ such that $Q$ has no cell $(r,c,\cdot)$,
$\Pr[L[r,c]=s\mid L\supseteq R\cup Q]\le\Delta/n$ with $\Delta=(D+1)/(1-2\alpha-\beta)+o(1)$,
which is $\le66$ for $\alpha,\beta$ small.  Consequently, for every partial Latin square
$Q$ on rows other than $1,2$ with $|R_Q|,|C_Q|\le\varepsilon n$ ($\varepsilon$ a small
absolute constant), \eqref{eq:FC} holds with $C=\Delta$:
$\Pr[Q\subseteq L\mid L\supseteq R]\le(\Delta/n)^{|Q|}$.
\end{corollary}

\begin{proof}
The first statement is Lemma~4.1 of \cite{AM26} with $\mathcal R=\{r\}$: its switch
$\eta_L(r,r_1,c)$, $r_1\notin E\cup\{r\}$, avoids the cells of $R\cup Q$ when $c_1\notin C_Q$
and $r_2\notin E$ (at most $|C_Q|+|E|$ excluded choices of $r_1$) and changes only rows
$r,r_1,r_2\notin\{1,2\}$; the images $L'$ all have $r'(L')=r_1\notin E$ (the row with
$L'[r_1,c]=s$), so in their sum (19) only $r'\notin\mathcal R\cup E$ occur, and for these
the inverse count uses Theorem~\ref{thm:AM31two} for $Q'=Q\cup\{(r',c,s),(r,c,s')\}$ with
$r'$ as its $r'$: the only cell of $Q'$ in row $r'$ is in column $c$, as the theorem
allows, and $\alpha,\beta$ increase by $O(1/n)$.  The second statement follows by
conditioning on the cells of $Q$ one at a time.
\end{proof}

\begin{corollary}[the $\nu$-tail]\label{cor:nutail}
For every $\lambda,\alpha,\beta$ and every $M$ with $2M+2\le\varepsilon n$,
$\Pr_X[\nu\le M]\le8\Delta^{2M+2}/n$.  In particular
$\Pr_X[\nu\le C\log\log n]=O\bigl((\log n)^{2C\log\Delta}/n\bigr)=o(1/\log n)$ for every $C$.
\end{corollary}

\begin{proof}
Corollary~\ref{cor:FCtwo} is \eqref{eq:FC} with $C=\Delta$, in every $S(R)$; the
computation after \eqref{eq:FC} in \S\ref{sec:orbit}(e) gives the bound.
\end{proof}

So the first term of the first-pair condition is settled, and the weak conjecture follows
from the second alone:

\begin{hypothesis}[hazard hypothesis]\label{hyp:hazard}
For some absolute $C$ and $M=M(n)\le C\log\log n$,
$\displaystyle\sup_{\lambda,\alpha,\beta}\Pr_X\bigl[\nu>M,\ \text{the candidates }t=1,\dots,M\text{ are all bad}\bigr]=o(1/\log n)$.
\end{hypothesis}

\begin{corollary}\label{cor:hazard}
Hypothesis~\ref{hyp:hazard} implies Hypothesis~\ref{hyp:firstpair}, hence
Conjecture~\ref{conj:cgw}.  It holds with $M=\lceil C\log\log n\rceil$ if, for $t\le M$,
the conditional probability that candidate $t$ is bad, given $\nu>t$ and that candidates
$1,\dots,t-1$ are bad, is at most $1-c$ for an absolute $c>0$, and $C\log(1/(1-c))>1$.
\end{corollary}

\begin{proof}
$\Pr_X[G^c]\le\Pr_X[\nu\le M]+\Pr_X[\nu>M,\ t\le M\text{ all bad}]$, and the first term is
$o(1/\log n)$ by Corollary~\ref{cor:nutail}; under the conditional bound the second is at
most $(1-c)^M=(\log n)^{-C\log(1/(1-c))}$.
\end{proof}

\begin{remark}[what changed, and what did not]\label{rem:fcchanged}
Hypothesis~\ref{hyp:hazard} is a statement about $O(\log\log n)$ events of constant
probability at column pairs selected by the frame-selected row pair; it is a lower bound
(each candidate is good with conditional probability $\ge c$) and is not a fixed-cell
event, so Corollary~\ref{cor:FCtwo} does not reach it.  The obstacle of
\S\ref{sec:remains} is therefore unchanged in kind but reduced in scope: the upper-bound
half --- short structures at a frame-selected pair --- is now available, through the
$\eta$-switching of \cite{AM26}, whose forward degree is deterministic; only lower bounds
at square-selected positions remain, and \S\ref{sec:hazardid} shows which: at the first
candidate, lengths at scale $n$ are available (Proposition~\ref{prop:lengths}), good
pairs across the cut between $x_1$ and $y_1$ are not.  In the data (Remark~\ref{rem:fixrect}) the
conditional probabilities of Corollary~\ref{cor:hazard} (that a candidate is bad) are
$\approx0.6$ at every $n$ and type, i.e.\ $c\approx0.4$.  Where the obstruction of our earlier analysis went: with the admissible rows of
their Claim~4 taken to be all free rows off $\gamma_L$, the column-cycle switch through
$x$ fails to be legal exactly when the cycle through $x$ contains one special row; the
statistic $g$ restricts to rows whose cycles contain none, and Claim~$3'$ shows that
these are $\Omega(n)$ on average, by cutting all active cycles at once with the cuts
confined to arcs.  Everything stays in $S(R)$; neither the mark nor the isotopy argument
of Remark~\ref{rem:fixrect} is used.  The constants ($D\to62$, $\Delta\le66$, against their
$22$ and $23$) are not optimised and are irrelevant for $M=O(\log\log n)$.  The argument
was checked by two independent automated referees (the first found a claim stated per arc
that is true only in aggregate, corrected above) and numerically on sampled completions
at $n=30,50,100$: every restricted cross-switch, cut and Claim-$4'$ move behaved as
stated, the aggregate of Claim~$2'$ held for every arc size with a factor about $2$ to
spare, and the mean of $g/f$ was $0.22$--$0.25$ against the bound $1/30$ \cite{repo}.
\end{remark}
"""
