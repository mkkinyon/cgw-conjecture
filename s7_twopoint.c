/* s7_twopoint.c -- two-point structure of the orbit bias (session 7).
 *
 * Continuum ensemble: u at 0.0, v at 0.5 on the unit circle; chords
 * have iid Uniform[0,1) endpoint pairs (a.s. disjoint).  For an ACTIVE
 * chord set T, define (cycle route, NO rank computation)
 *
 *     eps(T) = +1 if u ~ v in gamma tau_T, else -1
 *
 * (gamma = circle successor on the active endpoints + u,v).  By the
 * verified bit formula this equals the rank-event sign, so the cycle
 * route independently cross-checks the rank route (mode "bias").
 *
 * Key identity (session-7 factorization lemma): for S,T uniform
 * independent subsets of k iid chords, with j = |S cap T|,
 * d1 = |S\T|, d2 = |T\S|:
 *
 *   E_pos[eps_S eps_T] = E_C [ M_{d1}(C) M_{d2}(C) ],
 *   M_d(C) := E[ eps(C u F) | C ],  F = d fresh iid chords,
 *
 * because eps_S and eps_T are conditionally independent given the
 * positions C of the shared chords.  Hence
 *
 *   4 E[bias^2] = E_{(j,d1,d2)~multinomial} E_C[M_{d1} M_{d2}]
 *
 * and everything reduces to  psi(j,d) := E_C[ M_d(C)^2 ],
 * estimated unbiasedly by products of two independent fresh samples.
 *
 * Modes:
 *   bias k trials          4*E[bias^2] direct (two indep subset samples)
 *                          + E[eps_S] (signed one-point, subset-avgd)
 *   psi  j d cs fs         psi(j,d) via cs C-samples x fs fresh-pairs
 *   one  m trials          one-point f(m) = E[eps_m], all m active
 *
 * Usage: ./s7_twopoint mode args... [seed]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXP 4100  /* max points on circle = 2*chords + 2 */

static unsigned long long rng_s[2];
static inline unsigned long long rng_next(void) {
    unsigned long long s1 = rng_s[0], s0 = rng_s[1];
    rng_s[0] = s0;
    s1 ^= s1 << 23;
    rng_s[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
    return rng_s[1] + s0;
}
static inline double rng_u(void) {
    return (rng_next() >> 11) * (1.0 / 9007199254740992.0);
}

/* points: coords[i], partner[i] (index of chord partner, -1 for u,v).
 * eps for the active set = all listed chords.  n = number of points. */
static double coords[MAXP];
static int order[MAXP], rankp[MAXP], nxt[MAXP], part[MAXP];

static int cmp_pts(const void *x, const void *y) {
    double a = coords[*(const int *)x], b = coords[*(const int *)y];
    return (a > b) - (a < b);
}

/* build from chord coordinate array ch[2*nc] (pairs), u=0.0, v=0.5;
 * returns eps = +-1 */
static int eps_of(const double *ch, int nc) {
    int n = 2 * nc + 2;
    coords[0] = 0.0;  part[0] = -1;       /* u = point 0 */
    coords[1] = 0.5;  part[1] = -1;       /* v = point 1 */
    for (int c = 0; c < nc; c++) {
        coords[2 + 2 * c] = ch[2 * c];
        coords[3 + 2 * c] = ch[2 * c + 1];
        part[2 + 2 * c] = 3 + 2 * c;
        part[3 + 2 * c] = 2 + 2 * c;
    }
    for (int i = 0; i < n; i++) order[i] = i;
    qsort(order, n, sizeof(int), cmp_pts);
    for (int i = 0; i < n; i++) rankp[order[i]] = i;
    /* gamma: successor on circle; rho = gamma o tau: rho(p) =
       gamma(partner(p)) for chord points, gamma(p) for u,v */
    for (int i = 0; i < n; i++) {
        int p = order[i];
        int q = (part[p] >= 0) ? part[p] : p;   /* tau */
        nxt[p] = order[(rankp[q] + 1) % n];     /* gamma after tau */
    }
    /* wait: rho(p) = gamma(tau(p)); computed above with tau applied
       first at p itself.  walk the rho-cycle of u (=0), look for v(=1) */
    int p = 0, found = 0, steps = 0;
    do {
        p = nxt[p];
        if (p == 1) found = 1;
        if (++steps > n + 2) { fprintf(stderr, "cycle overflow\n"); exit(1); }
    } while (p != 0);
    return found ? 1 : -1;
}

static void fresh(double *ch, int nc) {
    for (int i = 0; i < 2 * nc; i++) ch[i] = rng_u();
}

int main(int argc, char **argv) {
    const char *mode = argc > 1 ? argv[1] : "bias";
    unsigned seed = 12345;
    if (!strcmp(mode, "bias") || !strcmp(mode, "one")) {
        if (argc > 4) seed = (unsigned)atoi(argv[4]);
    } else if (argc > 6) seed = (unsigned)atoi(argv[6]);
    rng_s[0] = seed * 2654435761ULL + 1;
    rng_s[1] = seed ^ 0x9e3779b97f4a7c15ULL;
    for (int i = 0; i < 64; i++) rng_next();

    static double ch[2 * (MAXP / 2)], sub1[2 * (MAXP / 2)], sub2[2 * (MAXP / 2)];

    if (!strcmp(mode, "bias")) {
        /* 4*E[bias^2] via two independent subset samples per config;
           also E[eps] (signed). */
        int k = argc > 2 ? atoi(argv[2]) : 64;
        long trials = argc > 3 ? atol(argv[3]) : 200000;
        double s2 = 0, s1 = 0;
        for (long t = 0; t < trials; t++) {
            fresh(ch, k);
            int n1 = 0, n2 = 0;
            for (int i = 0; i < k; i++) {
                if (rng_next() & 1) { sub1[2 * n1] = ch[2 * i]; sub1[2 * n1 + 1] = ch[2 * i + 1]; n1++; }
                if (rng_next() & 1) { sub2[2 * n2] = ch[2 * i]; sub2[2 * n2 + 1] = ch[2 * i + 1]; n2++; }
            }
            int e1 = eps_of(sub1, n1), e2 = eps_of(sub2, n2);
            s2 += e1 * e2;  s1 += e1;
        }
        printf("bias k=%d trials=%ld: 4Ebias2=%.6f (se %.6f)  Eeps=%+.6f\n",
               k, trials, s2 / trials, 2.0 / sqrt((double)trials),
               s1 / trials);
    } else if (!strcmp(mode, "psi")) {
        /* psi(j,d) = E_C[M_d(C)^2], unbiased: per C, average over fs
           pairs of independent fresh d-sets of eps1*eps2. */
        int j = argc > 2 ? atoi(argv[2]) : 0;
        int d = argc > 3 ? atoi(argv[3]) : 16;
        long cs = argc > 4 ? atol(argv[4]) : 2000;
        long fs = argc > 5 ? atol(argv[5]) : 50;
        double acc = 0;
        double accM = 0;  /* also E_C[M_d] = one-point f(j+d)-ish check */
        for (long t = 0; t < cs; t++) {
            fresh(ch, j);  /* shared block C */
            double inner = 0, innerM = 0;
            for (long r = 0; r < fs; r++) {
                memcpy(sub1, ch, sizeof(double) * 2 * j);
                memcpy(sub2, ch, sizeof(double) * 2 * j);
                fresh(sub1 + 2 * j, d);
                fresh(sub2 + 2 * j, d);
                int e1 = eps_of(sub1, j + d), e2 = eps_of(sub2, j + d);
                inner += e1 * e2;  innerM += 0.5 * (e1 + e2);
            }
            acc += inner / fs;  accM += innerM / fs;
        }
        printf("psi j=%3d d=%4d cs=%ld fs=%ld: psi=%.6f (se~%.6f)  Ef=%+.6f\n",
               j, d, cs, fs, acc / cs, 1.0 / sqrt((double)cs * fs),
               accM / cs);
    } else if (!strcmp(mode, "one")) {
        /* one-point f(m) = E[eps], all m chords active */
        int m = argc > 2 ? atoi(argv[2]) : 16;
        long trials = argc > 3 ? atol(argv[3]) : 1000000;
        double s = 0;
        for (long t = 0; t < trials; t++) {
            fresh(ch, m);
            s += eps_of(ch, m);
        }
        printf("one m=%4d trials=%ld: Eeps=%+.6f (se %.6f)\n",
               m, trials, s / trials, 1.0 / sqrt((double)trials));
    } else {
        fprintf(stderr, "unknown mode %s\n", mode);
        return 1;
    }
    return 0;
}
