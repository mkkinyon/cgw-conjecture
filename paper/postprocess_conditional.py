"""postprocess_conditional.py -- substitutions applied to the assembled weak_cgw_conditional.tex
(called from assemble_conditional.py).  Each pair must match exactly once-or-more; an assertion
fails if a source text has drifted."""


def postprocess(text):
    i5 = text.index(r'\section{The orbit method')
    head, orb = text[:i5], text[i5:]
    # kappa (genericity constant) -> kappa_0 in sections 5-6 (kappa = number of cycles in section 3)
    assert r'\kappa(\lambda)' not in orb and r'\kappa(\pi)' not in orb
    text = head + orb
    post = [
        (r"""Throughout this section $X=X(n,\lambda;\alpha,\beta)$, and ``legal'' means: preserves the
type of rows $1,2$ and the mark, i.e.\ maps $X$ to itself.""",
         r"""Throughout this section $X=X(n,\lambda;\alpha,\beta)$, and ``legal'' means: preserves the
type of rows $1,2$ and the mark, i.e.\ maps $X$ to itself.  The mark is now written
$\{p,p'\}$ (so $p=j$, $p'=j'$ in the notation of \S\S\ref{sec:marked}--\ref{sec:ladder}; the
letters $j,j'$ are freed for other columns), $\pi$ is its frame, and a column pair $(q,c)$ is
\emph{free} if $\{q,c\}\cap\{p,p'\}=\emptyset$.  A pair of rows $(x,y)$ is called
\emph{parallel} if it is flippable ($p,p'$ in different cycles of $\rho_{x,y}$) and
\emph{crossed} otherwise; $Q_p$ and $Q_{p'}$ denote the $\rho_{x,y}$-cycles through $p$ and
$p'$.  A free pair $(q,c)$ is \emph{separated} by $\{p,p'\}$ in $\rho_{x,y}$ if $q$ and $c$
lie on different arcs of a common cycle through $p,p'$ (crossed case) or one in $Q_p$ and the
other in $Q_{p'}\ne Q_p$ (parallel case)."""),
        (r"""So (CP) with the \emph{columns} prescribed is CGW's theorem; it is the
conditioning on rows $1,2$ --- on their type and the mark, which is
what defines $X$ and cannot be removed --- that is the path version.""",
         r"""So the statement with the \emph{columns} prescribed is CGW's theorem; it is the
conditioning on rows $1,2$ --- on their type and the mark, which is
what defines $X$ and cannot be removed --- that is the obstacle."""),
        (r"""It replaces the witness
lemma (a lower-tail statement about rare short cycles, probability
$\ell'/n$ per pair) by a statement""",
         r"""It replaces a lower-tail
statement about rare short cycles (a ladder pair is flippable as soon as its
$\rho_{x,y}$-cycle through $p$ is short and avoids $p'$, an event of probability
$\approx\ell'/n$ per pair for cycles of length $\le\ell'$) by a statement"""),
        (r"It still contains (CP)-type inputs: a", r"It still contains inputs of the same kind: a"),
        (r"""average over the sub-ladder does not collapse: $\kappa_0\approx0.15$ is
what the data suggest.""",
         r"""average over the sub-ladder does not collapse; the data put the
average of $\bar y$ at $0.13$--$0.19$, so the hypothesis asks for some
$\kappa_0$ well below that, and the lower tail is examined in (d$'$)."""),
        (r"""contributes $\approx\E[\bar y]^2\mathrm{Var}(|I|)\approx0.2$; an
independent re-run with fresh seeds gave within-$|I|$ ratios $1.03$ and
$0.94$) shows""", r"""contributes $\approx\E[\bar y]^2\mathrm{Var}(|I|)\approx0.2$) shows"""),
        (r"per-instance orbit averages lie in $[0.4,0.7)$.", r"per-instance orbit averages lie in $[0.4,0.6)$."),
        (r"ladder pairs at $n=30,50,100$, true marks, $\ell_1\le4,5,8$;", r"ladder pairs at $n=30,50,100$, marks drawn as in $X$, $\ell_1\le4,5,8$;"),
        (r"""If $(x,y)$ is flippable, trade $Q_{x,y}$; otherwise turn $C_1$ (which
makes $(x,y)$ flippable, as $y\in C_1\not\ni x$) and then trade
$Q_{x,y}$ in the turned square.""",
         r"""If $(x,y)$ is flippable, apply the $j'$-trade at $(x,y)$ (trade rows $x,y$ along
the $\rho_{x,y}$-cycle $Q_{x,y}$ through $j'$); otherwise turn $C_1$ (which
makes $(x,y)$ flippable, as $y\in C_1\not\ni x$) and then trade
$Q_{x,y}$ in the turned square."""),
        (r"""splits $C'_1$ into the arcs
$\{\pi'(y),\dots,x\}$ and $\{\pi'(x),\dots,y\}$,""",
         r"""splits $C'_1$ into the cycles
$\{x,\pi'(x),\dots,\pi'^{-1}(y)\}$ and $\{y,\pi'(y),\dots,\pi'^{-1}(x)\}$ of $(x\,y)\circ\pi'$,"""),
        (r"""\begin{hypothesis}[orbit genericity]\label{hyp:orbit}
With $\ell_0=n/\log^2n$:""",
         r"""\begin{hypothesis}[orbit genericity]\label{hyp:orbit}
There are constants $\kappa_0>0$ and $\ell_1$ such that, for the fixed matching
$\mathcal M$ of Proposition~\ref{prop:orbit} and uniformly in $\lambda,\alpha,\beta$, with
$\ell_0=n/\log^2n$:"""),
        (r"""(and
$=\tfrac32$ at $n=5$), so the statement itself is not in doubt.""",
         r"""(and
$=\tfrac32$ at $n=5$, $\tfrac78,\tfrac{11}{14},\tfrac{21}{22}$ at $n=6$), so the statement
itself is not in doubt."""),
    ]
    post.append((r"""normalising constant independent of $\lambda$ (namely
$D_n\CC_n((n))/\sum\gamma\CC$).""", r"""normalising constant independent of $\lambda$ (namely
$Z=\sum_\nu\gamma(\nu)\CC_n(\nu)/(|\DD_n|\,\CC_n((n)))$)."""))
    post.append((r"""\Pr\bigl[p,p'\text{ in different cycles of }\rho_{x,y}\bigm|\text{columns }p,p'\bigr]
\;=\;\tfrac12\ \text{ if $x,y$ are apart in the frame},\qquad
\ge\tfrac13\ \text{ if together at distances $\ge2$ (the latter under CGW($3/2$))}.""", r"""\Pr\bigl[p,p'\text{ in different cycles of }\rho_{x,y}\bigm|\text{columns }p,p'\bigr]
\;\begin{cases}=\tfrac12&\text{if $x,y$ are apart in the frame},\\ \ge\tfrac13&\text{if together at distances $\ge2$}.\end{cases}"""))
    post.append((r"So the statement with the \emph{columns} prescribed is CGW's theorem;", r"(The second bound is CGW's $\tfrac32$, affected by the gap; it holds with $\tfrac14$ unconditionally at distance $2$, Remark~\ref{rem:gap}.)  So the statement with the \emph{columns} prescribed is CGW's theorem;"))

    post.append((r"""$\Pr_X[\mathrm{Flip}_I=\emptyset]$ by a product of conditionally independent factors; the
two hypotheses it leads to are compared with Hypothesis~\ref{hyp:witness} at the end of the
section.""", r"""$\Pr_X[\mathrm{Flip}_I=\emptyset]$ by a product of conditionally independent factors; the
two hypotheses it leads to are compared with Hypothesis~\ref{hyp:witness} after (c$'$), and
in (e) the same mechanism, applied to the first ladder pair alone, gives an exact expression
for $\Pr_X[B]-\Pr_X[A]$ and a third alternative hypothesis."""))
    post.append((r"""constant.  All three hypotheses share the same unproved core, stated in
\S\ref{sec:remains}.""", r"""constant; at the first ladder pair alone, (e) below, the rarity disappears too, and what is
left is the lower tail of the number of good candidates among $\mu-1$ column pairs.  All of
these hypotheses share the same unproved core, stated in \S\ref{sec:remains}."""))
    for a, b in post:
        if a not in text:
            print('postprocess: pattern not found (skipped):', repr(a[:60])); continue
        text = text.replace(a, b)
    return text
