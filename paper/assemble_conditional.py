"""Assemble paper/weak_cgw_conditional.tex from the working notes (exact extraction of refereed blocks
plus fresh connective text).  Run from /home/claude/cgw: python3 paper/assemble_conditional.py"""
import re

def lines(path):
    return open(path).read().split('\n')

M = lines('section_marked.tex')
N = lines('second_row_notes.tex')
Lad = lines('section_ladder_s13.tex')
Mix = lines('section_mixing_s13.tex')

def block(src, a, b):            # 1-based inclusive
    return '\n'.join(src[a-1:b])

def find(src, pat, start=1):
    for i in range(start-1, len(src)):
        if pat in src[i]: return i+1
    raise KeyError(pat)

def sub(text, pairs):
    for a, b in pairs:
        assert a in text, a[:60]
        text = text.replace(a, b)
    return text

# ---------- blocks ----------
marked = block(M, 17, 159)
marked = sub(marked, [
 (r"Throughout, squares are normalised to have first row $\mathrm{id}$, so"+"\n"+r"the column map is $\omega=\sigma$ ($1\circ\omega(j)=2\circ j$ with the"+"\n"+r"first row the identity gives $\omega(j)=\sigma(j)$).",
  r"Throughout, squares are normalised to have first row $\mathrm{id}$; CGW's column map"+"\n"+r"$\omega$, defined by $L(1,\omega(j))=L(2,j)$, is then $\omega(j)=L(2,j)$, and we take"+"\n"+r"$\sigma=\omega$ (only the cycle type matters, and it is the same for $\omega$ and $\omega^{-1}$)."),
])

hyp_eq = block(N, find(N, r'\begin{hypothesis}[B-pair equidistribution'), find(N, r'\begin{hypothesis}[B-pair equidistribution')+3)
hyp_eq = sub(hyp_eq, [(r"For every pair $\mu\to\lambda$ as above with $\alpha\neq\beta$,"
                       "\n$\\bigl|\\bar q(\\mu\\to\\lambda)-1\\bigr|\\le\\delta(n)$.",
                       r"For every $\lambda\vdash n$ with a part $m=\alpha+\beta$, $\alpha\ne\beta$, $\alpha,\beta\ge2$, and $\mu$ the split type,"
                       "\n$\\bigl|2\\CC_n(\\lambda)/\\CC_n(\\mu)-1\\bigr|\\le\\delta(n)$.")])

t0 = find(N, r'\begin{theorem}\label{thm:reduction}')
t1 = find(N, r'\end{proof}', t0)
reduction = block(N, t0, t1)
reduction = sub(reduction, [
 (r"\ref{rem:cgwgap}", r"\ref{rem:gap}"),
 (r"Conjecture~\ref{conj:cgw}, hence also that almost all loops have"+"\n"+r"multiplication group $S_n$ and trivial character theory.",
  r"Conjecture~\ref{conj:cgw}."),
 (r"For each step of such a chain, Theorem~\ref{thm:identity} and", r"For each step of such a chain, Theorem~\ref{thm:marked} and"),
 (r"  The loops statement follows"+"\n"+r"from \eqref{eq:setTV}, Remark~\ref{rem:reduced} (the passage to reduced"+"\n"+r"squares is exact), and \S\ref{sec:known} as in \cite[Lemma 6.2]{CGW08}.", ""),
])

p0 = find(Lad, r'\begin{proposition}[the number of cycles, from the surviving direction only]')
p1 = find(Lad, r'\end{proof}', p0)
tail = block(Lad, p0, p1)

ladder = block(Lad, 9, 137)
ladder = sub(ladder, [
 (r"combinatorial statement.  Notation as in \S\ref{sec:s13trapped}:", r"Notation:"),
 (r"through $j'$ (\S\ref{sec:s13trapped}).", r"through $j'$."),
])
# strip the leading 'combinatorial statement.' sentence fragment
ladder = ladder.replace("Notation:\n$(L,P)", "Notation:\n$(L,P)")

off0 = find(Lad, r'\subsection{Offset ladders, and what (L) needs}')
off1 = find(Lad, r'\paragraph{The single remaining input.}') - 1
offsets = block(Lad, off0+1, off1)
offsets = sub(offsets, [
 (r"(the argument of Remark~\ref{rem:ladderdata}(iii) with"+"\n"+r"$n-|C_1|>n/2-\ell$ free rows)", r"(a splice-in switching with $n-|C_1|>n/2-\ell$ free rows, Remark~\ref{rem:splice})"),
])

tr0 = find(Mix, r'\begin{proposition}[short frame cycles through a marked row are rare]')
tr1 = find(Mix, r'\end{proof}', tr0)
trapped = block(Mix, tr0, tr1)

o_a0 = find(Lad, r'\paragraph{(a) Transposition: which half of the conditioning is the')
o_b0 = find(Lad, r'\paragraph{(b) Row trades are free in $X$; column turns are')
o_c0 = find(Lad, r'\paragraph{(c) The product bound.}')
o_w0 = find(Lad, r'\paragraph{What the hypothesis says.}')
o_d0 = find(Lad, r'\paragraph{(d) Data.}')
o_dp = find(Lad, r"\paragraph{(d$'$) The concentration input $(\gamma)$: data.}")
o_e0 = find(Lad, r'\paragraph{(e) The arc-trapping input $(\alpha)$')
orbit_a = block(Lad, o_a0, o_b0-1)
orbit_a = sub(orbit_a, [
 (r"\eqref{eq:cgwA}:", r"CGW's Theorem~3.13 (the inequality $\Pr_X[A]\ge\tfrac13$, whose published proof has a gap for $\min(\alpha,\beta)\ge3$, Remark~\ref{rem:gap}):"),
])
orbit_b = block(Lad, o_b0, o_c0-1)
orbit_b = sub(orbit_b, [
 (r"The effect on"+"\n"+r"$\rho_{x,y}$ is the computation of \S\ref{sec:crossspace}; separation"+"\n"+r"is preserved by the arc computation there.",
  r"The effect on $\rho_{x,y}$: if $x\in Z$ and $y\notin Z$ the entries of row $x$ at $q,c$ are exchanged, so $\rho'_{x,y}=\rho_{x,y}\circ(q\,c)$, and symmetrically; multiplying by a transposition of two points on one cycle splits it at those points, so $p,p'$ are separated iff $q,c$ lie on different arcs between them, and multiplying by a transposition of two points on different cycles merges the cycles, so $p,p'$ are joined iff $q,c$ lie one in each of their cycles; in both cases $q,c$ lie one on each side afterwards."),
])
orbit_c = block(Lad, o_c0, o_w0-1)
orbit_c = sub(orbit_c, [
 (r"Follow the proof of Proposition~\ref{prop:Lwitness}.  Long ladders:", r"Put $\ell_0=n/\log^2n$.  By Theorem~\ref{thm:ladder}, $|\Pr[B]-\Pr[A]|\le\Pr[\mathrm{Flip}=\emptyset]\le\Pr[K_0\ge\ell_0,\mathrm{Flip}=\emptyset]+\Pr[A,K_0<\ell_0]+\Pr[B,K_0<\ell_0]$.  Long ladders:"),
 (r"the parts $|C_1|<n/2+\ell$"+"\n"+r"and $K_0<\ell_0$ on $A$ are as before.",
  r"the part $|C_1|<n/2+\ell$ costs $\le8/n$ per $\ell$ by the splice-in bound, and $\Pr[A,K_0<\ell_0]\le72\ell_0/n$ by Proposition~\ref{prop:trapped}; the case $d_{21}=\ell$ is symmetric.  Altogether $\Pr[\mathrm{Flip}=\emptyset]=o(1/\log n)$."),
])
orbit_w = block(Lad, o_w0, o_d0-1)
orbit_w = sub(orbit_w, [
 (r"(the arc-trapping input (ii) of"+"\n"+r"\S\ref{sec:crossspace})", r"(the arc-trapping input $(\alpha)$ of \S\ref{sec:remains})"),
 (r"This is the recursion noted in \S\ref{sec:crossspace}: every", r"This is the recursion: every"),
 (r"(Lemma~\ref{lem:rowfair}, (F), (I), Proposition~\ref{prop:orbit})", r"(Lemma~\ref{lem:rowfair}, Proposition~\ref{prop:orbit}, and the one-step switchings of the marked space)"),
])
orbit_d = block(Lad, o_d0, o_e0-1)
orbit_d = sub(orbit_d, [
 (r"\paragraph{(d) Data.}  \texttt{s13\_orbit.py}"+"\n"+r"(\texttt{runs/s13/orbit/orbit\{30,50,100\}.log}; $600$, $600$ and",
  r"\paragraph{(d) Data.}  All scripts and logs referred to below are in the repository \cite{repo}, under \texttt{s13\_*.py} and \texttt{runs/s13/orbit/}.  The orbit test ($600$, $600$ and"),
 (r"Proposition~\ref{prop:orbit} (\texttt{--fixed};"+"\n"+r"\texttt{orbit\{50,100\}\_fixed.log}): $|\mathcal T|=1.7$",
  r"Proposition~\ref{prop:orbit}: $|\mathcal T|=1.7$"),
 (r"\texttt{s13\_cov.py} (\texttt{runs/s13/orbit/cov\{50,100\}.log}; $543$ and"+"\n"+r"$286$ instances) computes,",
  r"The covariance test ($543$ and $286$ instances at $n=50,100$) computes,"),
 (r"a"+"\n"+r"referee's re-run with fresh seeds", r"an"+"\n"+r"independent re-run with fresh seeds"),
])

# (f) block: eq:Pi and the block graph characterization
f0 = find(Lad, r'\paragraph{(f) The permutation problem $(\Pi)$')
f1 = find(Lad, r'Every input of Hypothesis~\ref{hyp:orbit} is now of one form.') - 1
pi_block = block(Lad, f0, f1)
pi_block = sub(pi_block, [
 (r"(``separated''"+"\n"+r"in \S(c) and Lemma~\ref{lem:legalturn})", r"(``separated'' in Lemma~\ref{lem:legalturn})"),
 (r"the"+"\n"+r"adversarial examples of \S(c) are instances of the disconnected case"+"\n"+r"(the distance-$2$ chords straddling $p$ or $p'$ join $P_1$ to $P_2$"+"\n"+r"and must be excluded there).", r"a family of chords at $\rho_0$-distance $2$, none straddling $p$ or $p'$, is a disconnected example."),
])

preamble = r"""\documentclass[11pt]{article}
\usepackage[margin=1.1in]{geometry}
\usepackage{amsmath,amssymb,amsthm,mathtools}
\usepackage{booktabs}
\usepackage{xcolor}
\usepackage[colorlinks=true,linkcolor=blue!50!black,citecolor=blue!50!black,urlcolor=blue!50!black]{hyperref}
\usepackage{enumitem}
\newtheorem{theorem}{Theorem}[section]
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{corollary}[theorem]{Corollary}
\newtheorem{conjecture}[theorem]{Conjecture}
\newtheorem{hypothesis}[theorem]{Hypothesis}
\theoremstyle{definition}
\newtheorem{definition}[theorem]{Definition}
\newtheorem{remark}[theorem]{Remark}
\DeclareMathOperator{\sgn}{sgn}
\newcommand{\dtv}{d_{\mathrm{TV}}}
\newcommand{\E}{\mathbb{E}}
\renewcommand{\Pr}{\mathbb{P}}
\newcommand{\DD}{\mathcal{D}}
\newcommand{\PP}{P}
\newcommand{\QQ}{Q}
\newcommand{\CC}{C}

\title{The cycle type of two rows of a random Latin square:\\ an exact ladder identity and a conditional proof of the weak Cavenagh--Greenhill--Wanless conjecture}
\author{[author]}
\date{2 October 2026 --- draft}

\begin{document}
\maketitle

\begin{abstract}
Let $L$ be a uniformly random Latin square of order $n$ and $\sigma$ the permutation of columns
carrying its first row to its second.  Cavenagh, Greenhill and Wanless conjectured that the
cycle type of $\sigma$ is asymptotically that of a uniform derangement, in total variation.
We prove three exact identities --- the marked-pair identity, which expresses the ratio of
completion counts of adjacent cycle types as the odds of one bit (``rows $1,2$ lie in a
common column cycle of a marked column pair'') in a single uniform square; a reduction
theorem, by which the conjecture follows from that bit being fair to $o(1/\log n)$; and
the ladder identity, by which the bit is exactly fair except on the squares in which no
``ladder pair'' of rows is flippable --- and an orbit method by which the probability of
the exceptional set is bounded by a product of conditionally independent factors.  The
conjecture is thereby reduced to a genericity hypothesis about the cycle structure of two
rows at four columns (Hypothesis~\ref{hyp:orbit}), which we state precisely, test
numerically up to $n=100$, and locate within the known switching technology.  Along the
way we observe that the published proof of the upper bound in Lemma~3.12 of
Cavenagh--Greenhill--Wanless has a gap for $\min(\alpha,\beta)\ge3$; our results do not
depend on that bound, and we re-derive from the surviving direction the tail estimate on
the number of cycles that the reduction theorem needs.
\end{abstract}

\section{Introduction}

Let $\PP_n$ be the law of the cycle type of $\sigma=\sigma_{1,2}(L)$ for a uniform Latin
square $L$ of order $n$, and $\QQ_n$ the law of the cycle type of a uniform derangement of
$[n]$.  Cavenagh, Greenhill and Wanless \cite[Conj.~6.1]{CGW08} conjectured:

\begin{conjecture}\label{conj:cgw}
$\dtv(\PP_n,\QQ_n)=o(1)$ as $n\to\infty$.
\end{conjecture}

For $\lambda\vdash n$ with all parts $\ge2$ let $\CC_n(\lambda)$ be the number of
completions of a fixed $2\times n$ Latin rectangle of cycle type $\lambda$ to a Latin
square, divided by $(n-2)!$ (the count does not depend on the rectangle chosen, since
rectangles of the same type are isotopic by a row-preserving isotopy), and let
$\gamma(\lambda)$ be the number of permutations of type $\lambda$.  Then
$\PP_n(\lambda)\propto\gamma(\lambda)\CC_n(\lambda)$ and
$\QQ_n(\lambda)\propto\gamma(\lambda)$, so the conjecture says that $\CC_n$ is nearly
constant on the types that carry most of the mass.  CGW proved
$\CC_n(\lambda)/\CC_n(\mu)\in[\tfrac12,\tfrac32]$ for adjacent types (one part of
$\lambda$ split into two parts of $\mu$, both $\ge2$), and used it to show that the two
laws agree within a factor $n^{3/2}$ on every type; the upper bound $\tfrac32$, and
with it the upper half of the latter statement, is affected by the gap described in
Remark~\ref{rem:gap} and \cite{gapnote}.

\paragraph{Results.}  Throughout, a \emph{mark} is a column pair $P=\{j,j'\}$ with
$j'=\sigma^\alpha(j)$ on an $m$-cycle of $\sigma$, $m=\alpha+\beta$, $\alpha,\beta\ge2$;
$X=X(n,\lambda;\alpha,\beta)$ is the set of marked squares of type $\lambda$; the
\emph{frame} of $P$ is the permutation of rows $\pi(r)=r'$ with $L(r',j')=L(r,j)$; and
$P$ is a \emph{$B$-pair} if rows $1,2$ lie in one cycle of $\pi$, an \emph{$A$-pair}
otherwise.  Our results are:
\begin{enumerate}[leftmargin=2em]
\item[(1)] \emph{Marked-pair identity} (Theorem~\ref{thm:marked}):
$2\CC_n(\lambda)/\CC_n(\mu)-1=\Pr_X[B]/\Pr_X[A]$ exactly, for every $n,\lambda,\alpha,\beta$.
\item[(2)] \emph{Reduction} (Theorem~\ref{thm:reduction}): if
$|\Pr_X[B]-\tfrac12|\le\delta(n)/4$ for all adjacent pairs with $\alpha\ne\beta$, then
$\dtv(\PP_n,\QQ_n)=O(\delta\log n+n^{-1+o(1)})$.  The tail estimate it needs on the number
of cycles of $\sigma$ is proved from the trivial half of (1) alone
(Proposition~\ref{prop:tailsurvive}).
\item[(3)] \emph{Ladder identity} (Theorem~\ref{thm:ladder}): with the ladder pairs
$(x_k,y_k)=(\pi^{-k}(1),\pi^{-k}(2))$, $1\le k<K_0$, and ``flippable'' meaning that $j,j'$
lie in different cycles of the column permutation of rows $x_k,y_k$,
\[
\Pr_X[B]-\Pr_X[A]=\Pr_X[B,\mathrm{Flip}=\emptyset]-\Pr_X[A,\mathrm{Flip}=\emptyset],
\qquad |\Pr_X[B]-\tfrac12|\le\tfrac12\Pr_X[\mathrm{Flip}=\emptyset].
\]
So Conjecture~\ref{conj:cgw} follows from (L): $\Pr_X[\mathrm{Flip}=\emptyset]=o(1/\log n)$.
\item[(4)] \emph{Orbit method} (Proposition~\ref{prop:orbit}): the flippabilities of a
family of ladder pairs are conditionally independent given the orbit of a group of
commuting column-cycle turns, and $\Pr_X[\mathrm{Flip}_I=\emptyset]\le\E_X\exp(-\tfrac12\sum_{k\in I}\bar y_k)$,
where $\bar y_k$ is an explicit ``orbit weight''.  Hence (L), and the conjecture, follow
from Hypothesis~\ref{hyp:orbit} (Proposition~\ref{prop:Lorbit}).
\end{enumerate}
Hypothesis~\ref{hyp:orbit} is a statement of the form ``two pairs of lines of a random
Latin square are interleaved in the cycle structure of a third pair with probability
bounded away from $0$ and $1$'', under conditioning on the cycle type of rows $1,2$.
Without that conditioning it is elementary (\S\ref{sec:remains}); the data give the
interleaving probability $0.19$ in every version we measured, and $\Pr_X[B]=0.50\pm0.01$ at
$n=30,50,100$.  Section~\ref{sec:remains} explains why the switching methods available
inside $X$ --- all of which we have pushed as far as they go --- produce only fair
identities and cannot by themselves give the hypothesis.

\paragraph{A gap in CGW's Lemma 3.12.}  Case 3 of the splitting procedure in
\cite{CGW08} (cross-switch at $\{\omega j,\omega j'\}$, then backflip at $\{j,j'\}$) does
not land in $S(\mu,F)$ when $\min(\alpha,\beta)\ge3$: the cross-switch reverses the
$\omega$-arc from $\omega j$ to $\omega j'$, which contains $j'$, so $\{j,j'\}$ ends at
distance $2$.  We give the details in a separate note \cite{gapnote}; here we only record
(Remark~\ref{rem:gap}) which published statements depend on it.  Nothing in the present
paper does: the lower bound $\tfrac12$ of CGW's lemma is the trivial half of
Theorem~\ref{thm:marked}, and the upper bound $\tfrac32$ is, by the same theorem,
$\Pr_X[B]\le\tfrac23$ --- a weak form of what we are after.  By Theorem~\ref{thm:ladder}
it would follow from $\Pr_X[\mathrm{Flip}\ne\emptyset]\ge\tfrac13$.

\paragraph{Conventions.}  $L(r,c)$ is the entry in row $r$, column $c$.  For rows $x,y$,
$\rho_{x,y}$ is the permutation of columns with $L(y,\rho_{x,y}(q))=L(x,q)$, so that
$\sigma=\rho_{1,2}$ up to inversion; a \emph{row cycle} of $(x,y)$ is a cycle of
$\rho_{x,y}$.  For columns $q,c$, $\delta_{q,c}$ is the permutation of rows with
$L(\delta(r),c)=L(r,q)$, and its cycles are the \emph{column cycles} of $(q,c)$.  A
\emph{trade} of rows $x,y$ along a row cycle swaps their entries on the columns of that
cycle; a \emph{turn} of a column cycle swaps the entries of columns $q,c$ in its rows.  Both
are Latin trades (the result is again a Latin square).  For $\lambda\vdash n$, $S(\lambda)$ is
the set of Latin squares of order $n$ with first row the identity whose rows $1,2$ have
cycle type $\lambda$, so $|S(\lambda)|=\gamma(\lambda)(n-2)!\,\CC_n(\lambda)$.
Jacobson--Matthews sampling
\cite{JM96} is used for all numerical data; the scripts and logs are in the project
repository.

\section{The marked-pair identity}\label{sec:marked}

"""

body = []
body.append(preamble)
body.append(marked)
body.append(r"""

\section{The reduction theorem}\label{sec:reduction}

""")
body.append(hyp_eq)
body.append("\n\nBy Theorem~\\ref{thm:marked}, $\\mathrm{EQ}(\\delta)$ is the statement $|\\Pr_X[B]-\\tfrac12|\\le\\tfrac\\delta4(1+o(1))$ for all adjacent pairs with $\\alpha\\ne\\beta$.  The proof below needs a tail bound on the number $\\kappa$ of cycles of $\\sigma$; we first prove it from the trivial half of Theorem~\\ref{thm:marked} alone, so that nothing depends on the upper bound of \\cite[Lemma~3.12]{CGW08} (Remark~\\ref{rem:gap}).\n\n")
body.append(tail)
body.append("\n\n")
body.append(reduction)
body.append(r"""

\begin{remark}[the gap in CGW's Lemma 3.12]\label{rem:gap}
In \cite{CGW08}, the direction $|S(\lambda,F)|N_\lambda\ge\tfrac12|S(\mu,F)|N_\mu$ of
Lemma~3.12 uses only that each splitting pair carries at most two edges and each joining
pair at least one; by Theorem~\ref{thm:marked} it is the trivial inequality $|X|\ge|X_A|$.
The direction $\le\tfrac32$ requires every splitting pair of every $\lambda$-square to carry
exactly two edges, and case~3 of the splitting procedure fails to supply them when
$\min(\alpha,\beta)\ge3$ \cite{gapnote}.  Consequently the upper bound of CGW's Theorem~3.13,
the lower bound of their Lemma~4.1, the parity bounds of Theorem~4.3, the bound
$\Pr[\sigma\text{ is an $n$-cycle}]\le2n^{-2/3}$ of Lemma~4.4 and Corollary~4.5 are not
established by the published argument.  The present paper uses none of them:
Proposition~\ref{prop:tailsurvive} replaces Corollary~4.5 in the proof of
Theorem~\ref{thm:reduction}.  The exact completion counts tabulated in \cite{CGW08} give
$\CC_n(\lambda)/\CC_n(\mu)$ within $0.4\%$ of $1$ for every adjacent pair at $7\le n\le11$,
so the statement itself is not in doubt.
\end{remark}

\section{The ladder identity}\label{sec:ladder}

""")
body.append(ladder)
body.append(r"""

\begin{remark}[data]\label{rem:ladderdata}
On Jacobson--Matthews samples with a uniformly random mark ($9000$, $7200$, $3600$ ladder-pair
instances at $n=30,50,100$): $\Pr[\mathrm{Flip}=\emptyset\mid K_0]=2^{-(K_0-1)}$ within
statistical error for every $K_0$; $\Pr[\mathrm{Flip}=\emptyset]=0.125,0.079,0.041\approx4/n$;
the ladder pairs are flippable with frequency $0.5002$ ($62082$ pairs at $n=50$); and
$\Pr[B]-\Pr[A]=-0.022\pm0.011$, $+0.008\pm0.012$, $-0.030\pm0.018$.
\end{remark}

\subsection{Offset ladders and short arcs}\label{sec:offsets}

""")
body.append(offsets)
body.append(r"""

\begin{remark}[the splice-in bound]\label{rem:splice}
The bound $\Pr[B,d_{12}=\ell,|C_1|<n/2+\ell]\le4/(n-2\ell+1)$ is the $B$-state version of
the following switching, which we state and prove for the cycle of row $1$.
\end{remark}

""")
body.append(trapped)
body.append(r"""

\section{The orbit method}\label{sec:orbit}

Throughout this section $X=X(n,\lambda;\alpha,\beta)$, and ``legal'' means: preserves the
type of rows $1,2$ and the mark, i.e.\ maps $X$ to itself.

""")
body.append(orbit_a); body.append("\n\n")
body.append(orbit_b); body.append("\n\n")
body.append(orbit_c); body.append("\n\n")
body.append(orbit_w); body.append("\n\n")
body.append(orbit_d)
body.append(r"""

\section{What remains}\label{sec:remains}

Hypothesis~\ref{hyp:orbit} has three inputs: $(\alpha)$ the cycles of $\rho_{x_k,y_k}$
through $j,j'$ are macroscopic for a constant fraction of the ladder pairs (both arcs when
they coincide, both cycles when they do not); $(\beta)$ \emph{genericity}: the clean
coordinates of a ladder pair are separated by $\{j,j'\}$ with probability bounded below;
$(\gamma)$ concentration of $\sum_k\bar y_k$ over the family.  Numerically $(\gamma)$ is an
independence statement (the orbit weights of different ladder pairs are uncorrelated to
within sampling error), $(\beta)$ holds with the same constant $0.19$--$0.20$ at every cycle
length, and $(\alpha)$ holds with the lengths spread uniformly.  All three are instances of
one statement:
\[
\textup{(U)}\qquad c\;\le\;\Pr\bigl[\text{two given pairs of lines are interleaved in the cycle structure of a third pair}\bigr]\;\le\;1-c
\]
in $X$ --- where two column pairs are interleaved in a row pair $(x,y)$ if they alternate
on a common cycle of $\rho_{x,y}$ or lie crosswise in two of its cycles, and dually for
row pairs in a column pair.  The pure permutation part of the problem is settled:

""")
body.append(pi_block)
body.append(r"""

Three facts locate the difficulty.  (i) Without conditioning on the type of rows $1,2$,
(U) for four generic columns in a generic row pair is elementary: by conjugation invariance
the four columns are uniform points given the cycle type of $\rho_{x,y}$, and the
interleaving probability is $\tfrac13\sum_i(m_i)_4/(n)_4+2\sum_{i\ne j}(m_i)_2(m_j)_2/(n)_4$
over the cycle lengths $m_i$, which lies in $[\tfrac1{27}-O(1/n),\tfrac13+O(1/n)]$.
(ii) For generic rows the dual statement is trivial: the number of row pairs interleaved with
$\{1,2\}$ in a column pair is at most $(n-2)^2/4$, half of all pairs, deterministically, and
rows other than $1,2$ are exchangeable in $X$.  What is needed is the same for a row pair
selected by its own row-cycle structure (the ladder pair, at the columns where its cycle
through $j$ meets a joining pair), and exchangeability is lost exactly there.  (iii)
Inside $X$ the exact tools are the free trades of rows other than $1,2$ (which preserve
interleaving) and the column-cycle turns (which are legal exactly off the interleaved
configurations), and every identity they produce is fair: it equates two expectations with
random weights whose lower bounds are again statements of the form (U).  Lower bounds on
``different cycles'' events are, in the language of \cite{CGW08}, of the type of the
direction of their Lemma~3.12 whose proof has the gap: they require a canonical
bounded-multiplicity repair of the same-cycle configurations, and the repairs available in
$X$ (trades or cross-switches of a free pair of third rows) have unbounded multiplicity.
What the orbit method adds to this is that the fairness compounds: a constant per-pair
genericity gives exponential decay in the family size, which is why a hypothesis of
constant strength suffices for a conclusion of strength $o(1/\log n)$.

\begin{thebibliography}{9}
\bibitem{CGW08} N.~J. Cavenagh, C. Greenhill and I.~M. Wanless, The cycle structure of two
rows in a random Latin square, \emph{Random Structures Algorithms} 33 (2008), 286--309.
\bibitem{gapnote} [author], A gap in the proof of Lemma 3.12 of Cavenagh--Greenhill--Wanless, note, October 2026.
\bibitem{JM96} M.~T. Jacobson and P. Matthews, Generating uniformly distributed random Latin
squares, \emph{J. Combin. Des.} 4 (1996), 405--437.
\bibitem{KS18} M. Kwan and B. Sudakov, Intercalates and discrepancy in random Latin squares,
\emph{Random Structures Algorithms} 52 (2018), 181--196.
\bibitem{repo} [author], \texttt{cgw-conjecture}: scripts, logs and notes, \url{https://github.com/mkkinyon/cgw-conjecture}, 2026.
\end{thebibliography}
\end{document}
""")
open('paper/weak_cgw_conditional.tex', 'w').write(''.join(body))
print('written')
