"""witness_subs.py -- substitutions applied to the witness block (hyp:witness, prop:Lwitness, remark) extracted
from section_ladder_s13.tex for the conditional paper.  Applied after the rename kappa -> kappa_1."""

WITNESS_SUBS = [
 (r"\paragraph{The single remaining input.}  Say that a pair of rows", r"Say that a pair of rows"),
 (r"(a cycle of length $\le\log^4n$ of a row-pair permutation through a",
  r"(a cycle of length $\le\log^5n$ of a row-pair permutation through a"),
 (r"(Remark~\ref{rem:ladderdata}(iv)), so a witness has probability" + "\n" + r"$\approx\ell'/n$.",
  r"(Remark~\ref{rem:witnessdata} below), so a witness has probability" + "\n" + r"$\approx\ell'/n$."),
 (r"witness probability is $\approx\ell'/n$, Remark~\ref{rem:ladderdata}(iv),",
  r"witness probability is $\approx\ell'/n$, Remark~\ref{rem:witnessdata},"),
 (r"The part with $|C_1|<n/2+\ell$ costs $\le8/n$.",
  r"The part with $|C_1|<n/2+\ell$ costs $\le8/n$ (Proposition~\ref{prop:splice})."),
 (r"""given column), in the space of completions of two rows and two
columns, for families of pairs selected by the frame.  It asks for a
lower bound on their expected number and a lower-tail
concentration (a second moment with relative variance $o(1/\log n)$
suffices) --- the territory of
the switching method \cite{KS18}: a switching that shortens or lengthens
the $j$-cycle of $\rho_{x,y}$ by a bounded amount, applied in the
conditional space, compares the counts of squares with
$j$-cycle lengths $m$ and $m+1$ and gives the uniform spread, and the
same switchings applied to two pairs at once control the covariance.""",
  r"""given column), in the space of completions of two rows and two
columns, for families of pairs selected by the frame.  It asks for a
lower bound on their expected number and a lower-tail
concentration (a second moment with relative variance $o(1/\log n)$
suffices).  For two \emph{fixed} rows, counts of such substructures
are what the switching method of \cite{KS18} controls; what it does
not address is a family of row pairs selected by the square, and that
is where the difficulty sits (\S\ref{sec:remains})."""),
 (r"""Let $(L,P)$ be uniform on $X$.  Let $\mathcal G$ be one of the
following families of \emph{blocks} of pairs of rows, determined by
the frame: (i) the ladder, each pair a block of its own; (ii) for a
$B$-instance with $d_{12}=\ell\le n/4$ and $|C_1|\ge n/2+\ell$, the
diagonals $D_c$, $0\le c\le n/2-\ell$, each a block of $\ell-1$ pairs.
Let $G_{\mathcal G}$ be the number of blocks containing a pair that
carries a witness, and
$\mu_{\mathcal G}=\sum_{\mathcal B\in\mathcal G}\min\bigl(\tfrac12,\ell'|\mathcal B|/n\bigr)$.
Then there is a constant $\kappa_1>0$ such that
\[
\Pr\bigl[G_{\mathcal G}<\kappa_1\,\mu_{\mathcal G}\ \text{ and }\ \mu_{\mathcal G}\ge\log^2n\bigr]=o(1/\log n).
\]""",
  r"""Let $(L,P)$ be uniform on $X$, uniformly in $\lambda,\alpha,\beta$.  Let
$\mathcal G$ be one of the following families of \emph{blocks} of pairs
of rows, determined by the frame: (i) the ladder, each pair a block of
its own; (ii) for a $B$-instance, with $\ell:=d_{12}$, provided
$\ell\le n/4$ and $|C_1|\ge n/2+\ell$, the diagonals $D_c$,
$0\le c\le n/2-\ell$, each a block of $\ell-1$ pairs (one family, not
one per $\ell$).  Let $G_{\mathcal G}$ be the number of blocks
containing a pair that carries a witness, and
$\mu_{\mathcal G}=\sum_{\mathcal B\in\mathcal G}\min\bigl(\tfrac12,\ell'|\mathcal B|/n\bigr)$.
Then, for the ladder (i),
\[
\Pr\bigl[G_{\mathcal G}=0\ \text{ and }\ \mu_{\mathcal G}\ge\log^2n\bigr]=o(1/\log n),
\]
and for the diagonals (ii) there is a constant $\kappa_1>0$ such that
\[
\Pr\bigl[G_{\mathcal G}<\kappa_1\,\mu_{\mathcal G}\ \text{ and }\ \mu_{\mathcal G}\ge\log^2n\bigr]=o(1/\log n).
\]"""),
 (r"""$\mu\ge(\ell_0-1)\ell'/n\ge\log^2n$, so $G\ge\kappa_1\mu\ge1$ except on an""",
  r"""$\mu\ge(\ell_0-1)\ell'/n\ge\log^2n$, so $G\ge1$ except on an"""),
 (r"""we treat $d_{12}=\ell$, and $d_{21}=\ell$ is the same argument with the
roles of rows $1,2$ exchanged (Proposition~\ref{prop:trapped} covers
$C_2$, and the offset families become those with the $y$-row deeper).""",
  r"""we treat $d_{12}=\ell$; the case $d_{21}=\ell$ reduces to it by exchanging
rows $1$ and $2$ (with the symbol relabelling that restores the first
row): this replaces $\sigma$ by $\sigma^{-1}$, keeps the mark $\{j,j'\}$
with $j'=\sigma^{-\beta}(j)$, so maps $X(n,\lambda;\alpha,\beta)$ onto
$X(n,\lambda;\beta,\alpha)$ measure-preservingly, conjugates the frame
by $(1\,2)$, hence exchanges $d_{12}$ with $d_{21}$, preserves the
ladder as a set of unordered pairs and carries the diagonals with the
$x$-row deeper to those with the $y$-row deeper, and leaves every
$\rho_{x,y}$ with $x,y\notin\{1,2\}$ unchanged; so the $d_{21}=\ell$ case
for $(\lambda,\alpha,\beta)$ is the $d_{12}=\ell$ case of the hypothesis
for $(\lambda,\beta,\alpha)$, which ``uniformly in $\lambda,\alpha,\beta$''
covers (for $\alpha=\beta$ the hypothesis is read for the ordered mark)."""),
 (r"""\[
\sum_{\ell<\ell_0}\Pr[B,d_{12}=\ell,|C_1|\ge n/2+\ell,\theta\ge\theta_0(\ell)]
\le\sum_{\ell-1\le n/(2\ell')}\frac{36\,n}{(n-2)\,\kappa_1\,(\ell-1)\ell'}
+\sum_{n/(2\ell')<\ell-1<\ell_0}\frac{72}{(n-2)\kappa_1}
=O\Bigl(\frac{\log n}{\ell'}+\frac{\ell_0}{n}\Bigr)=O\Bigl(\frac1{\log^2n}\Bigr),
\]""",
  r"""\begin{multline*}
\sum_{\ell<\ell_0}\Pr[B,d_{12}=\ell,|C_1|\ge n/2+\ell,\theta\ge\theta_0(\ell)]
\le\sum_{\ell-1\le n/(2\ell')}\frac{36\,n}{(n-2)\,\kappa_1\,(\ell-1)\ell'}
+\sum_{n/(2\ell')<\ell-1<\ell_0}\frac{72}{(n-2)\kappa_1}\\
=O\Bigl(\frac{\log n}{\ell'}+\frac{\ell_0}{n}\Bigr)=O\Bigl(\frac1{\log^2n}\Bigr),
\end{multline*}"""),
 (r"""are one admissible choice; the proof needs
$\ell_0\ell'/n\ge\log^2n$, $\ell_0/n=o(1/\log n)$ and
$\log n/\ell'=o(1/\log n)$.""",
  r"""are one admissible choice; the proof needs
$\ell_0\ell'/n\ge\log^2n$, $\ell'\ge4\log^2n$, $\ell_0/n=o(1/\log n)$ and
$\log n/\ell'=o(1/\log n)$."""),
]
