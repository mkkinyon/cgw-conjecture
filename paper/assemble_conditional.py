"""Assemble paper/weak_cgw_conditional.tex from the working notes (exact extraction of refereed blocks
plus fresh connective text).  Run from /home/claude/cgw: python3 paper/assemble_conditional.py"""
import re

def lines(path):
    return open(path).read().split('\n')

import sys as _s; _s.path.insert(0, 'paper')
from cprime_blocks import cprime_intro, cprime_rest, remains
from firstpair_block import firstpair_block
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
 (r"and let $Q$ be the row cycle through one column of"+"\n"+r"$P$ (the two row cycles are distinct).", r"and let $Q$ be the row cycle through one column of"+"\n"+r"$P$ (the two row cycles are distinct), selected by CGW's min-symbol rule."),
 (r"Throughout, squares are normalised to have first row $\mathrm{id}$, so"+"\n"+r"the column map is $\omega=\sigma$ ($1\circ\omega(j)=2\circ j$ with the"+"\n"+r"first row the identity gives $\omega(j)=\sigma(j)$).",
  r"Throughout, squares are normalised to have first row $\mathrm{id}$; CGW's column map"+"\n"+r"$\omega$, defined by $L(1,\omega(j))=L(2,j)$, is then $\omega(j)=L(2,j)$, and we take"+"\n"+r"$\sigma=\omega$ (only the cycle type matters, and it is the same for $\omega$ and $\omega^{-1}$)."),
])

hyp_eq = block(N, find(N, r'\begin{hypothesis}[B-pair equidistribution'), find(N, r'\begin{hypothesis}[B-pair equidistribution')+3)
hyp_eq = sub(hyp_eq, [(r"For every pair $\mu\to\lambda$ as above with $\alpha\neq\beta$,"
                       "\n$\\bigl|\\bar q(\\mu\\to\\lambda)-1\\bigr|\\le\\delta(n)$.",
                       r"For every $\lambda\vdash n$ with a part $m=\alpha+\beta$, $\alpha\ne\beta$, $\alpha,\beta\ge2$, and $\mu$ the split type,"
                       "\n$\\bigl|\\bar q(\\mu\\to\\lambda)-1\\bigr|\\le\\delta(n)$, where $\\bar q(\\mu\\to\\lambda):=2\\CC_n(\\lambda)/\\CC_n(\\mu)-1=\\Pr_X[B]/\\Pr_X[A]$ by Theorem~\\ref{thm:marked}.")])

t0 = find(N, r'\begin{theorem}\label{thm:reduction}')
t1 = find(N, r'\end{proof}', t0)
reduction = block(N, t0, t1)
reduction = sub(reduction, [
 (r"\ref{rem:cgwgap}", r"\ref{rem:gap}"),
 (r"Conjecture~\ref{conj:cgw}, hence also that almost all loops have"+"\n"+r"multiplication group $S_n$ and trivial character theory.",
  r"Conjecture~\ref{conj:cgw}."),
 (r"For each step of such a chain, Theorem~\ref{thm:identity} and", r"For each step of such a chain, Theorem~\ref{thm:marked} and"),
 (r"  The loops statement follows"+"\n"+r"from \eqref{eq:setTV}, Remark~\ref{rem:reduced} (the passage to reduced"+"\n"+r"squares is exact), and \S\ref{sec:known} as in \cite[Lemma 6.2]{CGW08}.", ""),
 (r"\dtv(\PP_n,\QQ_n)\;=\;O\bigl(\delta\,\log n\;+\;n^{-2/3}\bigr).", r"\dtv(\PP_n,\QQ_n)\;=\;O\bigl(\delta\,\log n\;+\;n^{-1+o(1)}\bigr)."),
 (r"=\frac{\gamma(\lambda)}{\gamma((n))}\,(1+\varepsilon)^{\pm(\kappa-1)},"+"\n"+r"\qquad |\varepsilon|\le\delta,",
  r"=\frac{\gamma(\lambda)}{\gamma((n))}\prod_{i=1}^{\kappa-1}(1+\varepsilon_i)^{\pm1},"+"\n"+r"\qquad |\varepsilon_i|\le\delta,"),
 (r"$\PP_n(\lambda)=\QQ_n(\lambda)\,Z^{-1}(1+O(K\delta))$ where $Z>0$ is a", r"$\PP_n(\lambda)=\QQ_n(\lambda)\,Z^{-1}(1+O(K\delta))$ (the bound is vacuous unless $K\delta=O(1)$, which is the case of interest) where $Z>0$ is a"),
])
# remove the "original argument, kept for the record" passage
i0 = reduction.index(r"; the"+"\n"+r"original argument, kept for the record,")
i1 = reduction.index(r"On the complement of $\mathcal B$")
reduction = reduction[:i0] + ".\n\n" + reduction[i1:]
reduction = sub(reduction, [
 (r"shows $Z^{-1}=1+O(K\delta)+O(n^{-2/3})$, whence"+"\n"+r"$\sum_{\lambda\notin\mathcal B}|\PP_n-\QQ_n|=O(K\delta)+O(n^{-2/3})$.",
  r"shows $Z^{-1}=1+O(K\delta)+O(n^{-1+o(1)})$, whence"+"\n"+r"$\sum_{\lambda\notin\mathcal B}|\PP_n-\QQ_n|=O(K\delta)+O(n^{-1+o(1)})$."),
 (r"$\dtv(\PP_n,\QQ_n)=O(\delta\log n+n^{-2/3})$ (with"+"\n"+r"Proposition~\ref{prop:tailsurvive} in place of \cite[Cor.~4.5]{CGW08}"+"\n"+r"the error term becomes $O(\delta\log n+n^{-1+o(1)})$).",
  r"$\dtv(\PP_n,\QQ_n)=O(\delta\log n+n^{-1+o(1)})$."),
])

p0 = find(Lad, r'\begin{proposition}[the number of cycles, from the surviving direction only]')
p1 = find(Lad, r'\end{proof}', p0)
tail = block(Lad, p0, p1)
tail = sub(tail, [(r"$\E[\ell N_\ell]\le2n/(n-\ell)$ for every $\ell$;", r"$\E[\ell N_\ell]\le2n/(n-\ell)$ for every $\ell<n$;")])

ladder = block(Lad, 9, 137)
ladder = sub(ladder, [
 (r"n^{-2/3}", r"n^{-1+o(1)}"),
 (r"The same identity holds with $X$ replaced by $S(R)\times\{P\}$ for any"+"\n"+r"fixed admissible pair of first rows $R$ and mark $P$.", r"The same identity holds with $X$ replaced by the set of squares with a fixed pair of"+"\n"+r"rows $1,2$ (of type $\lambda$) and a fixed mark $P$."),
 (r"combinatorial statement.  Notation as in \S\ref{sec:s13trapped}:", r"Notation:"),
 (r"through $j'$ (\S\ref{sec:s13trapped}).", r"through $j'$."),
])
# strip the leading 'combinatorial statement.' sentence fragment
ladder = ladder.replace("Notation:\n$(L,P)", "Notation:\n$(L,P)")

off0 = find(Lad, r'\subsection{Offset ladders, and what (L) needs}')
off1 = find(Lad, r'\paragraph{The single remaining input.}') - 1
offsets = block(Lad, off0+1, off1)
offsets = sub(offsets, [
 (r"(the argument of Remark~\ref{rem:ladderdata}(iii) with"+"\n"+r"$n-|C_1|>n/2-\ell$ free rows)", r"(Proposition~\ref{prop:splice} below)"),
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
 (r"\eqref{eq:cgwA}:", r"CGW's Theorem~3.13 (the inequality $\Pr_X[A]\ge\tfrac13$, whose published proof has a gap for every split other than $(2,2)$, Remark~\ref{rem:gap}):"),
])
orbit_b = block(Lad, o_b0, o_c0-1)
orbit_b = sub(orbit_b, [
 (r"The effect on"+"\n"+r"$\rho_{x,y}$ is the computation of \S\ref{sec:crossspace}; separation"+"\n"+r"is preserved by the arc computation there.",
  r"The effect on $\rho_{x,y}$: if $x\in Z$ and $y\notin Z$ the entries of row $x$ at $q,c$ are exchanged, so $\rho'_{x,y}=\rho_{x,y}\circ(q\,c)$, and symmetrically; multiplying by a transposition of two points on one cycle splits it at those points, so $p,p'$ are separated iff $q,c$ lie on different arcs between them, and multiplying by a transposition of two points on different cycles merges the cycles, so $p,p'$ are joined iff $q,c$ lie one in each of their cycles; in both cases $q,c$ lie one on each side afterwards."),
])
o_h0 = find(Lad, r'\begin{hypothesis}[orbit genericity]')
orbit_c = block(Lad, o_c0, o_h0-1)
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

# witness block
w0 = find(Lad, r'\paragraph{The single remaining input.}')
w1 = find(Lad, r'\subsection{Status of the witness lemma: literature and data}') - 1
witness = block(Lad, w0, w1)
witness = witness.replace(r'\kappa', r'\kappa_1')
from witness_subs import WITNESS_SUBS
witness = sub(witness, WITNESS_SUBS)

# (j) block: adaptive coordinates
j0 = find(Lad, r'\begin{lemma}[the adaptive coordinate is orbit-invariant]')
j1 = find(Lad, r'\end{proof}', find(Lad, r'\begin{proposition}[orbit bound with adaptive coordinates]'))
adaptive = block(Lad, j0, j1)

# (f) block: eq:Pi and the block graph characterization
f0 = find(Lad, r'\paragraph{(f) The permutation problem $(\Pi)$')
f1 = find(Lad, r'Every input of Hypothesis~\ref{hyp:orbit} is now of one form.') - 1
pi_block = block(Lad, f0, f1).replace(r'\paragraph{(f) The permutation problem', r'\paragraph{(e) The permutation problem')
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

\title{The cycle type of two rows of a random Latin square:\\ exact identities and a reduction of the weak Cavenagh--Greenhill--Wanless conjecture}
\author{[author]}
\date{3 October 2026 --- draft}

\begin{document}
\maketitle

\begin{abstract}
Let $L$ be a uniformly random Latin square of order $n$ and $\sigma$ the permutation of columns
carrying its first row to its second.  Cavenagh, Greenhill and Wanless conjectured that the
cycle type of $\sigma$ is asymptotically that of a uniform derangement, in total variation.
We prove three exact identities and a reduction theorem: the marked-pair identity, which
expresses the ratio of completion counts of adjacent cycle types as the odds of one bit
(``rows $1,2$ lie in a common column cycle of a marked column pair'') in a single uniform
square; the reduction theorem, by which the conjecture follows from that bit being fair to
$o(1/\log n)$; and the ladder identity, by which the bit is exactly fair except on the
squares in which no ``ladder pair'' of rows is flippable.  An orbit method bounds the probability of
the exceptional set by a product of conditionally independent factors.  The
conjecture thereby follows from any one of four sufficient conditions, none of which is
known to imply another, all concerning the probability of a configuration at rows or
columns selected by the square inside the space conditioned on the type of rows $1,2$, and
all needed only for types with $O(\log n)$ parts.  The reference condition is a ``witness
lemma'' (Hypothesis~\ref{hyp:witness}): with probability $1-o(1/\log n)$ some ladder pair
of a long ladder has a short row cycle through one marked column avoiding the other.  The
simplest to state comes from a third exact identity at the first ladder pair alone, which
expresses the bias $\Pr_X[B]-\Pr_X[A]$ as an expectation over the event that none of the
column pairs along the row permutation of that pair admits a legal toggling turn; the
conjecture follows if that event has probability $o(1/\log n)$, and, given a constant
conditional hazard along the row permutation, its failure is dominated by a short arc at
the marked columns rather than by a large count of events.  We
state these precisely, report the numerical evidence (up to $n=150$, including, through a
sampler of completions of a fixed rectangle whose connectivity and mixing are not proved,
types of exponentially small probability), and locate the remaining difficulty: every
exact identity we have found inside the conditioned space is fair, and what is missing is
a constant-factor estimate for events on a few cells in two fixed rows of a uniform
completion of a fixed $2\times n$ rectangle, outside the range of the known fixed-cell
bounds for Latin squares and rectangles.  Along the
way we observe that the published proof of the upper bound in Lemma~3.12 of
Cavenagh--Greenhill--Wanless has a gap for every split other than $(2,2)$; our results do not
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
laws agree within a factor $n^{3/2}$ on every type with at most $\tfrac65\log n$ parts; the upper bound $\tfrac32$, and
with it the upper half of the latter statement, is affected by the gap described in
Remark~\ref{rem:gap} and \cite{gapnote}.  One marginal of the conjecture is known: by
Kwan, Petrova and Sawhney's law of large numbers for the number of odd rows
\cite[Thm.~1.3(1)]{KPS25} and the exchangeability of rows,
$\PP_n(\sigma\text{ even})=\tfrac12+o(1)$, as for uniform derangements
(\S\ref{sec:remains}).

\paragraph{Status of the arguments.}  The proofs below have been checked by automated
referees (independent language-model agents with access to the sources and the data) and,
where feasible, by exhaustive computation at $n\le7$ and by sampling (this includes
\S\ref{sec:orbit}(e) and Remark~\ref{rem:fewparts}, added last); no human referee has read
them.  The literature statements of \S\ref{sec:remains} were checked against the cited
papers by the same means.

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
$|\Pr_X[B]/\Pr_X[A]-1|\le\delta(n)$ for all adjacent pairs with $\alpha\ne\beta$ in which
$\mu$ has at most $\lceil16\log n\rceil$ parts
(equivalently $|\Pr_X[B]-\tfrac12|\le\tfrac14\delta(1+o(1))$), then
$\dtv(\PP_n,\QQ_n)=O(\delta\log n+n^{-1+o(1)})$.  The tail estimate it needs on the number
of cycles of $\sigma$ is proved from the trivial half of (1) alone
(Proposition~\ref{prop:tailsurvive}).
\item[(3)] \emph{Ladder identity} (Theorem~\ref{thm:ladder}): with the ladder pairs
$(x_k,y_k)=(\pi^{-k}(1),\pi^{-k}(2))$, $1\le k<K_0$, and ``flippable'' meaning that $j,j'$
lie in different cycles of the permutation $\rho_{x_k,y_k}$ of columns induced by rows $x_k,y_k$,
\[
\Pr_X[B]-\Pr_X[A]=\Pr_X[B,\mathrm{Flip}=\emptyset]-\Pr_X[A,\mathrm{Flip}=\emptyset],
\qquad |\Pr_X[B]-\tfrac12|\le\tfrac12\Pr_X[\mathrm{Flip}=\emptyset].
\]
So Conjecture~\ref{conj:cgw} follows from (L): $\Pr_X[\mathrm{Flip}=\emptyset]=o(1/\log n)$.
\item[(4)] \emph{Witness lemma} (Proposition~\ref{prop:Lwitness}): (L), hence the
conjecture, follows from Hypothesis~\ref{hyp:witness}: with probability $1-o(1/\log n)$,
a long ladder contains a pair $(x_k,y_k)$ whose $\rho_{x_k,y_k}$-cycle through $j$ has
length $\le\log^5n$ and avoids $j'$ (such a pair is flippable outright), and, for the offset
diagonals of short-arc instances, that a constant fraction of the expected number of
diagonals carry one.  This is the reference hypothesis of the paper.
\item[(5)] \emph{Orbit method} (Propositions~\ref{prop:orbit} and~\ref{prop:adaptive}):
the flippabilities of a family of ladder pairs are conditionally independent given the
orbit of a group of commuting column-cycle turns.  With a fixed matching of columns as
coordinates this gives $\Pr_X[\mathrm{Flip}_I=\emptyset]\le\E_X\exp(-\tfrac12\sum_{k\in I}\bar y_k)$
with an explicit ``orbit weight'' $\bar y_k$; with the \emph{adaptive} coordinates
$(\rho_k^t(j),\rho_k^t(j'))$ --- the $t$-th columns along the permutation $\rho_k=\rho_{x_k,y_k}$
of the $k$-th ladder pair from the marked columns --- every coordinate is separated by the
mark, so its turn toggles whenever it is legal and clean, and
$\Pr_X[\mathrm{Flip}_I=\emptyset]\le\E_X[2^{-|U|}]$, where $U$ is the set of ladder pairs
having a short clean candidate cycle (subject to a collision condition).  These give two
alternative sufficient conditions (Hypotheses~\ref{hyp:orbitgen} and~\ref{hyp:adaptive}),
for which we see no implication to or from Hypothesis~\ref{hyp:witness}: the first needs a
genericity constant that the switchings available inside $X$ cannot provide, the second
needs $C\log\log n$ rare local events on a long ladder where the witness lemma needs one.
\item[(6)] \emph{The first ladder pair alone} (Propositions~\ref{prop:firstpair}
and~\ref{prop:firstpairbias}): at the first pair $(x_1,y_1)$ the adaptive turn needs no
side conditions, and combined with the first ladder trade it gives
\[
\Pr_X[B]-\Pr_X[A]=\E_X\bigl[(\mathbf 1_B-\mathbf 1_A)(\mathbf 1_{F^c}-\mathbf 1_F);\,G^c\bigr],
\qquad\bigl|\Pr_X[B]-\tfrac12\bigr|\le\tfrac12\Pr_X[G^c],
\]
where $F$ is the flippability of $(x_1,y_1)$ and $G$ the event that one of the $\nu-1$
column pairs $\{\rho^t(j),\rho^t(j')\}$ along $\rho=\rho_{x_1,y_1}$ ($\nu$ = the length of
the shorter arc or cycle of $\rho$ at the mark) induces a permutation of rows that separates
$x_1$ from $y_1$ with a legal cycle.  So the conjecture also follows from
Hypothesis~\ref{hyp:firstpair}: $\sup_{\lambda,\alpha,\beta}\Pr_X[G^c]=o(1/\log n)$.  With
$M=C\log\log n$, $\Pr_X[G^c]\le\Pr_X[\nu\le M]+\Pr_X[\nu>M,\text{ the first }M\text{ candidates bad}]$,
and the second term is $o(1/\log n)$ for $C$ large (depending on $c$) as soon as the
conditional probability that candidate $t\le M$ is bad, given $\nu>t$ and that its
predecessors are bad, is at most $1-c$; so the condition is
an upper bound on a short arc at the frame-selected pair together with a constant
conditional hazard over $O(\log\log n)$ events, in place of a long ladder of rare events,
and it is the one hypothesis whose quantity the data measure directly ($n\Pr_X[G^c]$ is
$5.5$ for $30\le n\le150$, a first-moment fact that says nothing about the rate).  Again
no implication to or from Hypothesis~\ref{hyp:witness} is known; the $\nu$-tail is not
shown to be easier than the rest (Lemma~\ref{lem:nuswitch}).
\end{enumerate}
All four hypotheses are of one logical type: a constant-factor estimate for the
probability of a configuration at rows or columns selected by the frame of the mark, inside
the space conditioned on the type of rows $1,2$; by Remark~\ref{rem:fewparts} each is
needed only for types with at most $\lceil16\log n\rceil$ parts.  Their heuristic margins
are large; none of the exact tools of this paper gives an estimate of that kind, and the
fixed-cell bounds in the literature for Latin squares and rectangles stop short of the
conditioning on two full rows (\S\ref{sec:remains}, where we also record how far the most
recent of them, Allsop and Morris's switching argument \cite{AM26}, goes in the conditioned
space).  The data (Jacobson--Matthews samples, and a sampler of completions of a fixed
$2\times n$ rectangle that reaches types of $\PP_n$-probability $e^{-\Theta(n)}$,
Remark~\ref{rem:fixrect}; its connectivity and mixing are not proved) give point estimates
of $\Pr_X[B]$ within $0.003$ of $\tfrac12$ in every tested class at $n\le50$ and witness
events that are independent to three digits along the ladder; they do not reach the regime
of any of the hypotheses, and the empirical sizes quoted for the hypothesised quantities
are first-moment information only.

\paragraph{A gap in CGW's Lemma 3.12.}  Case 3 of the splitting procedure in
\cite{CGW08} (cross-switch at $\{\omega j,\omega j'\}$, then backflip at $\{j,j'\}$) does
not behave as claimed: the cross-switch reverses the $\omega$-arc from $\omega j$ to
$\omega j'$ (or the complementary arc, depending on the smallest corner symbol), which
contains $j'$ (resp.\ $j$), so $\{j,j'\}$ ends at distance $2$.  For $\min(\alpha,\beta)\ge3$
the backflip then lands in the wrong type; for $\min(\alpha,\beta)=2<\max(\alpha,\beta)$ the
type is right but the edge is not an edge of the joining multigraph, so the identity
$G_S=G_J$ on which the double count rests fails (already at $n=5$, where we computed both
multigraphs exactly).  We give the details in a separate note \cite{gapnote}; here we only record
(Remark~\ref{rem:gap}) which published statements depend on it.  Nothing in the present
paper does: the lower bound $\tfrac12$ of CGW's lemma is the trivial half of
Theorem~\ref{thm:marked}, and the upper bound $\tfrac32$ is, by the same theorem,
$\Pr_X[B]\le\tfrac23$ --- a weak form of what we are after.  By Theorem~\ref{thm:ladder}
it would follow from $\Pr_X[\mathrm{Flip}=\emptyset]\le\tfrac13$.

\paragraph{Conventions.}  $L(r,c)$ is the entry in row $r$, column $c$.  For rows $x,y$,
$\rho_{x,y}$ is the permutation of columns with $L(y,\rho_{x,y}(q))=L(x,q)$, so that
$\sigma=\rho_{1,2}$ up to inversion; a \emph{row cycle} of $(x,y)$ is a cycle of
$\rho_{x,y}$.  For columns $q,c$, $\delta_{q,c}$ is the permutation of rows with
$L(\delta(r),c)=L(r,q)$, and its cycles are the \emph{column cycles} of $(q,c)$.  A
\emph{trade} of rows $x,y$ along a row cycle swaps their entries on the columns of that
cycle; a \emph{turn} of a column cycle swaps the entries of columns $q,c$ in its rows.  Both
are Latin trades (the result is again a Latin square).  For $\lambda\vdash n$, $S(\lambda)$ is
the set of Latin squares of order $n$ with first row the identity whose rows $1,2$ have
cycle type $\lambda$, so $|S(\lambda)|=\gamma(\lambda)(n-2)!\,\CC_n(\lambda)$; we also
write $D_\lambda=\gamma(\lambda)$ for the number of derangements of type $\lambda$.  A
column pair $\{j,j'\}$ is \emph{admissible} if $j'\notin\{j,\sigma(j),\sigma^{-1}(j)\}$
(CGW's condition for the pair operations; every marked pair is admissible since
$\alpha,\beta\ge2$).  Following CGW, the \emph{flip} at an $A$-pair whose columns lie in
different $\sigma$-cycles, and the \emph{backflip} (or \emph{unflip}) at an $A$-pair whose
columns lie in one $\sigma$-cycle, is the turn of the column cycle of the pair through row
$2$ (it joins, resp.\ splits, the $\sigma$-cycles involved); the \emph{switch} at a pair
whose columns lie in different $\sigma$-cycles is the trade of rows $1,2$ along the row
cycle through one of the two columns, selected by CGW's rule (the cycle whose columns
carry the smaller minimum symbol).  We call Conjecture~\ref{conj:cgw}, which asks only
for convergence in total variation, the \emph{weak} conjecture; the stronger statements
suggested by Cameron (``tends very rapidly'') and by the data (pointwise ratios
$1+o(1)$) are not addressed here.
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
reduction = sub(reduction, [
 (r"Let $K=\lceil A\log n\rceil$ with $A$ a large absolute constant, and let", r"Let $K=\lceil A\log n\rceil$ with $A=16$, and let"),
 (r"Assume $\mathrm{EQ}(\delta)$ with $\delta=\delta(n)\to0$.",
  r"Assume $\mathrm{EQ}(\delta)$ with $\delta=\delta(n)\to0$, or only its instances in which $\mu$ has at most $\lceil16\log n\rceil$ parts."),
])
body.append(reduction)
body.append(r"""

\begin{remark}[which types are needed]\label{rem:fewparts}
The proof uses $\mathrm{EQ}(\delta)$ only along the chains from $(n)$ to the plain types
with at most $K=\lceil16\log n\rceil$ parts, i.e.\ for adjacent pairs $\mu\to\lambda$ with
$\kappa(\mu)\le K$ (and $\alpha\ne\beta$; every intermediate type of a chain is plain,
since the residual part always exceeds every part already split off, so $\lambda$ and
$\mu$ are both plain).  Consequently each of the sufficient conditions below
(Hypotheses~\ref{hyp:witness}, \ref{hyp:orbitgen}, \ref{hyp:adaptive},
\ref{hyp:firstpair}), which are stated as suprema over all $\lambda,\alpha,\beta$, is
needed only for plain types $\lambda$ with fewer than $K$ parts and plain splits $\mu$.
Asymptotically this excludes the types with $\Theta(n)$ parts, such as the families
$(4,2^{(n-4)/2})$ of Remark~\ref{rem:fixrect} (at the sizes we sample, $n\le150$, every
tested type still has fewer than $\lceil16\log n\rceil$ parts; the pairs
$(4,2^{13})\to(2^{15})$ and $(6,3^8)\to(3^{10})$ there are not needed because $\alpha=\beta$
and $\mu$ is not plain); it does not remove the conditioning itself, which is on two full
rows whatever the type.
\end{remark}

\begin{remark}[the gap in CGW's Lemma 3.12]\label{rem:gap}
In \cite{CGW08}, the direction $|S(\lambda,F)|N_\lambda\ge\tfrac12|S(\mu,F)|N_\mu$ of
Lemma~3.12 uses only that each splitting pair carries at most two edges and each joining
pair at least one; by Theorem~\ref{thm:marked} it is the trivial inequality $|X|\ge|X_A|$.
The direction $\le\tfrac32$ requires every splitting pair of every $\lambda$-square to carry
exactly two edges that are also edges of the joining multigraph, and case~3 of the splitting
procedure fails to supply them for every split other than $(2,2)$ \cite{gapnote}.
Consequently the upper bound of CGW's Theorem~3.13, the lower bound of their Lemma~4.1,
the parity bounds of Theorem~4.3, the bound $\Pr[\sigma\text{ is an $n$-cycle}]\le2n^{-2/3}$
of Lemma~4.4, both halves of Corollary~4.5, and the results of their Section~4.2 that rest
on these (Corollary~4.6, Lemma~4.7, Theorems~4.9 and~4.11, Corollary~4.10 --- the last
superseded by \cite{KS18}) are not established by the published argument.  The
\emph{statements} of Theorem~4.9 (at most $9\sqrt n$ cycles with probability
$1-o(e^{-\sqrt n/4})$) and of the upper half of Corollary~4.5 (with $n/2$ in place of
$n^{1/3}$) do follow from the surviving direction: Proposition~\ref{prop:tailsurvive}(ii)
with $k=\sqrt n$ gives $\Pr[\kappa\ge9\sqrt n]\le\exp(-(\tfrac12-o(1))\sqrt n\log n)$, and
$\CC_n(\lambda)\le2^{\kappa-1}\CC_n((n))$ gives $\PP_n(\lambda)\le\tfrac n2\prod_i(2/c_i)$
for a type with parts $c_i$; this repairs the use of both in \cite[Thm.~9]{GMW25}.  The
surviving direction gives no upper bound on coarse types: in particular we know of no proof
that two rows of a random Latin square form a single $n$-cycle with probability $o(1)$
(Lemma~4.4's $2n^{-2/3}$; the second-moment method with \cite{KS18} and
Proposition~\ref{prop:tailsurvive}(i) gives only $\le\tfrac56+o(1)$), and the parity theorem
of \cite{CW16} uses the unproved direction at odd splits (details in \cite{gapnote}); that
theorem is expected to be true independently of the gap, being a special case of
\cite[Thm.~1.3(4)]{KPS25}, whose written proof uses \cite{CW16} as a black box but whose
Remark~6.6 sketches an independent re-proof, which we have not checked.  For
$\min(\alpha,\beta)=2$ the count can be corrected with the constant $2$ in place of $\tfrac32$
\cite[Prop.~3]{gapnote}: for marks at distance $2$ or $m-2$ on an $m$-cycle,
$\CC_n(\lambda)/\CC_n(\mu)\le2$, i.e.\ $\Pr_X[A]\ge\tfrac14$ unconditionally.  The present
paper uses none of the affected statements of \cite{CGW08}:
Proposition~\ref{prop:tailsurvive} replaces Corollary~4.5 in the proof of
Theorem~\ref{thm:reduction}.  The exact completion counts tabulated in \cite{CGW08} give
$\CC_n(\lambda)/\CC_n(\mu)\in[0.985,1.0035]$ for every adjacent pair at $7\le n\le11$ (and
$=\tfrac32$ at $n=5$), so the statement itself is not in doubt.
\end{remark}

\section{The ladder identity}\label{sec:ladder}

""")
body.append(ladder)
body.append(r"""

\begin{remark}[data]\label{rem:ladderdata}
On Jacobson--Matthews samples with a uniformly random mark ($9000$, $7200$, $3000$ ladder-pair
instances at $n=30,50,100$; \texttt{runs/s13/f1/}): $\Pr[\mathrm{Flip}=\emptyset\mid K_0]=2^{-(K_0-1)}$ within
statistical error for $K_0\le6$ and consistent with it beyond; $\Pr[\mathrm{Flip}=\emptyset]=0.125,0.079,0.041\approx4/n$;
the ladder pairs are flippable with frequency $0.5002$ ($62082$ pairs at $n=50$); and
$\Pr[B]-\Pr[A]=-0.022\pm0.011$, $+0.008\pm0.012$, $-0.030\pm0.018$.  On completions of a
fixed rectangle (Remark~\ref{rem:fixrect}; $7\cdot10^5$ marked instances at $n=30$)
$\Pr[\mathrm{Flip}=\emptyset\mid K_0-1\ge K]=2^{-K}$ within noise to $K=12$--$15$, in
every tested type.
\end{remark}

\subsection{Offset ladders and short arcs}\label{sec:offsets}

""")
body.append(offsets)
body.append(r"""

""")
body.append(trapped)
body.append(r"""

\begin{proposition}[the splice-in bound]\label{prop:splice}
For $2\le\ell\le n/4$,
$\Pr_X[B,\ d_{12}=\ell,\ |C_1|<n/2+\ell]\le4/(n-2\ell+1)\le8/n$, and the same with
$d_{21}$ in place of $d_{12}$.
\end{proposition}

\begin{proof}
Let $E$ be the event.  \emph{Forward moves.}  From $(L,P)\in E$ choose $y=\pi^i(1)$ with
$1\le i\le\ell-1$ (a row strictly inside the arc from $1$ to $2$; $\ell-1$ choices) and a row
$x\notin C_1$ ($n-|C_1|\ge(n-2\ell+1)/2$ choices), and let $C_x$ be the frame cycle of $x$,
which contains neither $1$ nor $2$.  If $(x,y)$ is not flippable, turn $C_x$: this trades
columns $j,j'$ on the rows of $C_x$ only, so it preserves $X$, inverts $\pi$ on $C_x$
(Lemma~\ref{lem:cycletrade}), fixes row $y$, and replaces $\rho_{x,y}$ by
$\rho_{x,y}\circ(j\,j')$, which makes $(x,y)$ flippable.  Then apply the $j'$-trade at
$(x,y)$, which replaces the frame by $(x\,y)\circ\pi^{\pm}$ ($\pi^\pm=\pi$, or $\pi$ inverted
on $C_x$): the arc from $1$ to $2$ becomes
$1\to\dots\to\pi^{-1}(y)\to x\to\dots\to y\to\dots\to2$ with $C_x$ occupying
positions $i,\dots,i+|C_x|-1$, so $d_{12}$ grows by $|C_x|$, rows $1,2$ are untouched, and the result
is again in $X_B$ with the same mark.  \emph{Backward moves.}  Given the target
$(L',P)\in X_B$ and $\ell$, a preimage is determined by the position $i\in\{1,\dots,\ell-1\}$
of $x$ on the arc of $L'$ from $1$ to $2$ (then $|C_x|=d'_{12}-\ell$ and
$y=\pi'^{\,i+|C_x|}(1)$ are determined), together with the binary choice of whether $C_x$ was turned; undoing is the
same $j'$-trade (the pair stays flippable, $\rho'_{x,y}=\rho_{x,y}^{-1}$ on the traded
cycle) followed by the turn.  So each target has at most $2(\ell-1)$ preimages while each
source has at least $(\ell-1)(n-2\ell+1)/2$ forward moves, whence
$|E|\,(\ell-1)(n-2\ell+1)/2\le2(\ell-1)\,|X|$.
\end{proof}
""")
body.append(r"""

\section{The witness hypothesis}\label{sec:witness}

Notation as in \S\ref{sec:ladder}: $(L,P)\in X$, $P=\{j,j'\}$, frame $\pi$, ladder pairs
$(x_k,y_k)$, $K_0$, and for $B$-instances the offset diagonals $D_c$ of
\S\ref{sec:offsets}.  The statement below is the reference hypothesis of this paper --- the
least demanding sufficient condition for (L) that we have; the orbit method of
\S\ref{sec:orbit} gives two others.

""")
body.append(witness)
body.append(r"""

\begin{remark}[data]\label{rem:witnessdata}
On Jacobson--Matthews samples with a uniformly random mark (\texttt{runs/s13/f1/}), the
length of the $j$-cycle of $\rho_{x_k,y_k}$ for a ladder pair is spread like a uniform
variable on $\{2,\dots,n\}$: $\Pr[\le\ell]=0.021,0.081,0.18,0.38$ for $\ell=2,5,10,20$ at
$n=50$ against $(\ell-1)/(n-1)=0.020,0.082,0.18,0.39$, and $0.010,0.041,0.092,0.19$ at
$n=100$; among instances with $K_0\ge n/4$ a flippable ladder pair with $j$-cycle of length
$\le10$ exists in $93\%$ ($n=50$) and $95\%$ ($n=100$) of the instances, and one of length
$\le20$ in $99.6\%$ and $99.8\%$.  These samples see only typical types $\lambda$; the
hypothesis is required uniformly in $\lambda$.  The fixed-rectangle sampler of
Remark~\ref{rem:fixrect} reaches the rare types: there the witness events along the ladder
are independent Bernoulli variables to three digits (flat hazard,
$\Pr[\text{no witness in }K\text{ pairs}]=(1-w_1)^K$ to $K=40$ at $n=100$) and the witness
counts across the offset diagonals have binomial lower tails, identically for $(n)$ and for
$(4,2^{(n-4)/2})$.  Classes with all cycles of $\rho_{1,2}$ of length $\le3$ carry no
admissible mark and do not enter.
\end{remark}

\section{The orbit method: alternative sufficient conditions}\label{sec:orbit}

Throughout this section $X=X(n,\lambda;\alpha,\beta)$, and ``legal'' means: preserves the
type of rows $1,2$ and the mark, i.e.\ maps $X$ to itself.  The method bounds
$\Pr_X[\mathrm{Flip}_I=\emptyset]$ by a product of conditionally independent factors; the
two hypotheses it leads to are compared with Hypothesis~\ref{hyp:witness} at the end of the
section.

""")
body.append(orbit_a); body.append("\n\n")
body.append(orbit_b); body.append("\n\n")
body.append(orbit_c); body.append("\n\n")
body.append(cprime_intro); body.append("\n\n")
body.append(adaptive); body.append("\n\n")
body.append(cprime_rest)
body.append(firstpair_block)
body.append(remains)
body.append(r"""

\begin{thebibliography}{9}
\bibitem{AM26} J. Allsop and P. Morris, Universal probability bounds for partial Latin squares,
arXiv:2606.18174 (2026).
\bibitem{AW25} J. Allsop and I.~M. Wanless, Subsquares in random Latin rectangles,
\emph{Combinatorica} 45 (2025), art.~29.
\bibitem{CGW08} N.~J. Cavenagh, C. Greenhill and I.~M. Wanless, The cycle structure of two
rows in a random Latin square, \emph{Random Structures Algorithms} 33 (2008), 286--309.
\bibitem{CW16} N.~J. Cavenagh and I.~M. Wanless, There are asymptotically the same number of
Latin squares of each parity, \emph{Bull. Aust. Math. Soc.} 94 (2016), 187--194.
\bibitem{DKKS26} A. Divoux, T. Kelly, C. Kennedy and J. Sidhu, Subsquares in random Latin
squares and rectangles, \emph{J. Combin. Des.} 34 (2026), 184--197.
\bibitem{gapnote} [author], A gap in the proof of Lemma 3.12 of Cavenagh--Greenhill--Wanless, note, October 2026.
\bibitem{GMW25} M.~J. Gill, A. Mammoliti and I.~M. Wanless, Canonical labeling of Latin
squares in average-case polynomial time, \emph{Random Structures Algorithms} 66 (2025), e70015.
\bibitem{JM96} M.~T. Jacobson and P. Matthews, Generating uniformly distributed random Latin
squares, \emph{J. Combin. Des.} 4 (1996), 405--437.
\bibitem{MW99} B.~D. McKay and I.~M. Wanless, Most Latin squares have many subsquares,
\emph{J. Combin. Theory Ser.~A} 86 (1999), 322--347.
\bibitem{KPS25} M. Kwan, K. Petrova and M. Sawhney, Parities in random Latin squares,
arXiv:2509.13125 (2025).
\bibitem{KS18} M. Kwan and B. Sudakov, Intercalates and discrepancy in random Latin squares,
\emph{Random Structures Algorithms} 52 (2018), 181--196.
\bibitem{KSS21} M.~Kwan, A.~Sah and M.~Sawhney, Large deviations in random Latin squares,
\emph{Bull. London Math. Soc.} 54 (2022), 1420--1438.
\bibitem{KSSS23} M.~Kwan, A.~Sah, M.~Sawhney and M.~Simkin, Substructures in Latin squares,
\emph{Israel J. Math.} 256 (2023), 363--416.
\bibitem{repo} [author], \texttt{cgw-conjecture}: scripts, logs and notes, \url{https://github.com/mkkinyon/cgw-conjecture}, 2026.
\end{thebibliography}
\end{document}
""")
import sys; sys.path.insert(0, 'paper')
from postprocess_conditional import postprocess
open('paper/weak_cgw_conditional.tex', 'w').write(postprocess(''.join(body)))
print('written')
