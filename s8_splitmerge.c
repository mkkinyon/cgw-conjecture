/* s8_splitmerge.c -- continuum split-merge (coagulation-fragmentation)
 * representation of the fresh-chord averaging operator P (session 8).
 *
 * IDENTIFICATION (Lemma lem:splitmerge, session 8): on the circle of
 * unit length with marked points u=0.0, v=0.5, the cycles of
 * gamma tau_W partition the circle (Lebesgue) into "curves"; adding
 * one fresh chord with iid Uniform endpoints (c1,c2) performs ONE
 * uniform split-merge move on this measured partition:
 *   - c1,c2 in the same curve (prob |B|^2): split it along the two
 *     cut points, positions uniform along the curve;
 *   - c1,c2 in different curves: merge them (cut-and-join), insertion
 *     positions uniform along each.
 * eps(W) = +1 iff u,v lie on the same curve.  The relevant marked
 * state is: (together: sides x,y = curve-lengths of the two u->v
 * arcs) or (apart: sizes p,q of the curves through u,v); all updates
 * need only uniform variates (see session-8 notes).
 *
 * Modes:
 *   psi  j d cs fs [seed]   psi(j,d)=E_C[M_d(C)^2]: C = j iid chords
 *                           (exact discrete surgery), then two indep
 *                           split-merge continuations of d steps;
 *                           matches s7_twopoint psi interface.
 *   psig j d cs fs [seed]   same, but binned by the initial "gap"
 *                           g0 = min(x,y) (together) / min(p,q)
 *                           (apart), log2 bins; also by parity.
 *   one  m trials [seed]    one-point f(m)=E[eps(F_m)] from the empty
 *                           start (single curve, x=y=1/2).
 *   traj d trials [seed]    from empty start: E[eps_d], E[M?]  --
 *                           quick sanity of the chain itself.
 *
 * gcc -O3 -o s8_splitmerge s8_splitmerge.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

/* ----------------- marked-partition state ----------------- */
/* blocks 0..nb-1 with sizes sz[]; ub, vb = indices of u's / v's
 * block; if ub==vb, (sx,sy) = sides (u->v, v->u) along that curve.
 * Free-list-less: blocks are swapped-from-end on deletion. */
typedef struct {
    int nb, ub, vb, cap;
    double *sz;
    double sx, sy;      /* valid iff ub==vb */
} St;

static void st_init(St *S, int cap) {
    S->cap = cap;
    S->sz = malloc(sizeof(double) * cap);
}
static void st_copy(St *D, const St *S) {
    D->nb = S->nb; D->ub = S->ub; D->vb = S->vb;
    D->sx = S->sx; D->sy = S->sy;
    memcpy(D->sz, S->sz, sizeof(double) * S->nb);
}

/* sample a block index with prob proportional to size (linear scan;
 * nb stays modest, and this keeps the code transparent) */
static inline int pick_block(const St *S) {
    double r = rng_u();
    double acc = 0;
    for (int i = 0; i < S->nb; i++) {
        acc += S->sz[i];
        if (r < acc) return i;
    }
    return S->nb - 1;
}

static inline void del_block(St *S, int i) {
    int last = S->nb - 1;
    S->sz[i] = S->sz[last];
    if (S->ub == last) S->ub = i;
    if (S->vb == last) S->vb = i;
    S->nb--;
}

/* one split-merge step (one fresh chord) */
static void step(St *S) {
    int i1 = pick_block(S), i2 = pick_block(S);
    if (i1 == i2) {                      /* split */
        double s = S->sz[i1];
        if (i1 == S->ub && i1 == S->vb) {       /* together block */
            double x = S->sx, y = S->sy;
            double c1 = rng_u() * s, c2 = rng_u() * s;
            /* positions along cycle: u at 0, v at x */
            int c1x = c1 < x, c2x = c2 < x;
            if (c1x && c2x) {                   /* both cuts in x-side */
                double piece = fabs(c1 - c2);
                S->sz[i1] = s - piece;  S->sx = x - piece;
                S->sz[S->nb++] = piece;
            } else if (!c1x && !c2x) {          /* both in y-side */
                double piece = fabs(c1 - c2);
                S->sz[i1] = s - piece;  S->sy = y - piece;
                S->sz[S->nb++] = piece;
            } else {                            /* FLIP: separate u,v */
                double t1 = c1x ? c1 : c2;      /* in (0,x)  */
                double t2 = c1x ? c2 : c1;      /* in (x,s)  */
                double q = t2 - t1;             /* v's piece */
                S->sz[i1] = s - q;              /* u keeps i1 */
                S->sz[S->nb] = q;
                S->vb = S->nb++;
            }
        } else {                                /* <=1 mark in block */
            double t1 = rng_u() * s, t2 = rng_u() * s;
            double piece = fabs(t1 - t2);
            /* if a mark is in the block it sits at position 0, i.e.
             * in the complement piece s-piece (cuts uniform rel. mark);
             * for unmarked blocks the labelling is irrelevant */
            S->sz[i1] = s - piece;
            S->sz[S->nb++] = piece;
        }
    } else {                              /* merge i1,i2 */
        int hasu1 = (i1 == S->ub), hasv1 = (i1 == S->vb);
        int hasu2 = (i2 == S->ub), hasv2 = (i2 == S->vb);
        double p = S->sz[i1], q = S->sz[i2];
        if ((hasu1 && hasv2) || (hasv1 && hasu2)) {   /* FLIP: join u,v */
            double pu = hasu1 ? p : q;   /* u's old block size */
            double pv = hasu1 ? q : p;
            double a = rng_u() * pu, b = rng_u() * pv;
            S->sz[i1] = p + q;
            S->ub = S->vb = i1;
            S->sx = a + b;  S->sy = (p + q) - (a + b);
            del_block(S, i2);
        } else if ((hasu1 && hasv1)) {      /* together-block + env */
            double s = p;
            double c = rng_u() * s;
            if (c < S->sx) S->sx += q; else S->sy += q;
            S->sz[i1] = p + q;
            del_block(S, i2);
        } else if ((hasu2 && hasv2)) {
            double s = q;
            double c = rng_u() * s;
            if (c < S->sx) S->sx += p; else S->sy += p;
            S->sz[i2] = p + q;
            del_block(S, i1);
        } else {                            /* at most one mark each,
                                               not u-and-v across */
            S->sz[i1] = p + q;
            if (S->ub == i2) S->ub = i1;
            if (S->vb == i2) S->vb = i1;
            del_block(S, i2);
        }
    }
}

static inline int eps_st(const St *S) { return S->ub == S->vb ? 1 : -1; }

/* ------------- exact initial state from j iid chords ------------- */
/* discrete surgery on 2j+2 points; u=0.0 (pt 0), v=0.5 (pt 1);
 * computes the curve partition with arc measures, fills St */
#define MAXP 8300
static double coords[MAXP];
static int order[MAXP], rankp[MAXP], nxt[MAXP], part[MAXP], cyc[MAXP];

static int cmp_pts(const void *x, const void *y) {
    double a = coords[*(const int *)x], b = coords[*(const int *)y];
    return (a > b) - (a < b);
}

static void init_from_chords(St *S, int j) {
    int n = 2 * j + 2;
    coords[0] = 0.0;  part[0] = -1;
    coords[1] = 0.5;  part[1] = -1;
    for (int c = 0; c < j; c++) {
        coords[2 + 2 * c] = rng_u();
        coords[3 + 2 * c] = rng_u();
        part[2 + 2 * c] = 3 + 2 * c;
        part[3 + 2 * c] = 2 + 2 * c;
    }
    for (int i = 0; i < n; i++) order[i] = i;
    qsort(order, n, sizeof(int), cmp_pts);
    for (int i = 0; i < n; i++) rankp[order[i]] = i;
    for (int i = 0; i < n; i++) {
        int p = order[i];
        int q = (part[p] >= 0) ? part[p] : p;
        nxt[p] = order[(rankp[q] + 1) % n];   /* rho = gamma o tau */
    }
    /* cycles of rho; measure of a cycle = sum over its points p of
     * the arc [p, circle-successor(p)) */
    for (int i = 0; i < n; i++) cyc[i] = -1;
    S->nb = 0;
    for (int i = 0; i < n; i++) {
        if (cyc[i] >= 0) continue;
        int id = S->nb++;
        double m = 0;
        int p = i;
        do {
            cyc[p] = id;
            double a0 = coords[p];
            double a1 = coords[order[(rankp[p] + 1) % n]];
            m += (p == order[n - 1]) ? (1.0 - a0 + coords[order[0]])
                                     : (a1 - a0);
            p = nxt[p];
        } while (p != i);
        S->sz[id] = m;
    }
    S->ub = cyc[0];
    S->vb = cyc[1];
    if (S->ub == S->vb) {
        /* side x = measure from u to v along the curve */
        double x = 0;
        int p = 0;
        while (p != 1) {
            double a0 = coords[p];
            double a1 = coords[order[(rankp[p] + 1) % n]];
            x += (p == order[n - 1]) ? (1.0 - a0 + coords[order[0]])
                                     : (a1 - a0);
            p = nxt[p];
        }
        S->sx = x;
        S->sy = S->sz[S->ub] - x;
    }
}

int main(int argc, char **argv) {
    const char *mode = argc > 1 ? argv[1] : "psi";
    unsigned seed = 12345;
    if (!strcmp(mode, "one") || !strcmp(mode, "traj")) {
        if (argc > 4) seed = (unsigned)atoi(argv[4]);
    } else if (argc > 6) seed = (unsigned)atoi(argv[6]);
    rng_s[0] = seed * 2654435761ULL + 1;
    rng_s[1] = seed ^ 0x9e3779b97f4a7c15ULL;
    for (int i = 0; i < 64; i++) rng_next();

    if (!strcmp(mode, "psi") || !strcmp(mode, "psig")) {
        int j = argc > 2 ? atoi(argv[2]) : 16;
        int d = argc > 3 ? atoi(argv[3]) : 16;
        long cs = argc > 4 ? atol(argv[4]) : 2000;
        long fs = argc > 5 ? atol(argv[5]) : 50;
        int binned = !strcmp(mode, "psig");
        St C, A, B;
        int cap = 2 * j + 2 * d + 64;
        st_init(&C, cap); st_init(&A, cap); st_init(&B, cap);
        double acc = 0, accM = 0;
        /* gap bins: bin b = floor(-log2(g0)), 0..29; plus together/apart */
        double bacc[2][30] = {{0}}; long bcnt[2][30] = {{0}};
        for (long t = 0; t < cs; t++) {
            init_from_chords(&C, j);
            int tog = (C.ub == C.vb);
            double g0 = tog ? (C.sx < C.sy ? C.sx : C.sy)
                            : (C.sz[C.ub] < C.sz[C.vb] ? C.sz[C.ub]
                                                       : C.sz[C.vb]);
            int b = (int)(-log2(g0 > 1e-9 ? g0 : 1e-9));
            if (b < 0) b = 0;
            if (b > 29) b = 29;
            double inner = 0, innerM = 0;
            for (long r = 0; r < fs; r++) {
                st_copy(&A, &C); st_copy(&B, &C);
                for (int s = 0; s < d; s++) step(&A);
                for (int s = 0; s < d; s++) step(&B);
                int e1 = eps_st(&A), e2 = eps_st(&B);
                inner += e1 * e2;  innerM += 0.5 * (e1 + e2);
            }
            acc += inner / fs;  accM += innerM / fs;
            bacc[tog][b] += inner / fs;  bcnt[tog][b]++;
        }
        printf("SM psi j=%3d d=%4d cs=%ld fs=%ld: psi=%.6f (se~%.6f)"
               "  Ef=%+.6f\n", j, d, cs, fs, acc / cs,
               1.0 / sqrt((double)cs * fs), accM / cs);
        if (binned) {
            printf("  gap bins (g0 in 2^-(b+1)..2^-b): "
                   "[state] b mass psi|bin contrib\n");
            for (int tg = 0; tg < 2; tg++)
                for (int b = 0; b < 30; b++)
                    if (bcnt[tg][b])
                        printf("  [%s] %2d  %.4f  %+.4f  %+.5f\n",
                               tg ? "tog " : "apart", b,
                               (double)bcnt[tg][b] / cs,
                               bacc[tg][b] / bcnt[tg][b],
                               bacc[tg][b] / cs);
        }
    } else if (!strcmp(mode, "psibar")) {
        /* psibar d cs fs burn [seed]: stationary-start psi (j -> inf):
         * burn split-merge steps from the empty start, then two indep
         * continuations of d steps; gap-binned like psig. */
        int d = argc > 2 ? atoi(argv[2]) : 256;
        long cs = argc > 3 ? atol(argv[3]) : 20000;
        long fs = argc > 4 ? atol(argv[4]) : 100;
        int burn = argc > 5 ? atoi(argv[5]) : 16384;
        if (argc > 6) { /* reseed with trailing seed */
            seed = (unsigned)atoi(argv[6]);
            rng_s[0] = seed * 2654435761ULL + 1;
            rng_s[1] = seed ^ 0x9e3779b97f4a7c15ULL;
            for (int i = 0; i < 64; i++) rng_next();
        }
        St C, A, B;
        int cap = 2 * burn + 2 * d + 64;
        st_init(&C, cap); st_init(&A, cap); st_init(&B, cap);
        double acc = 0;
        double bacc[2][40] = {{0}}; long bcnt[2][40] = {{0}};
        for (long t = 0; t < cs; t++) {
            C.nb = 1; C.sz[0] = 1.0; C.ub = C.vb = 0;
            C.sx = 0.5; C.sy = 0.5;
            for (int s = 0; s < burn; s++) step(&C);
            int tog = (C.ub == C.vb);
            double g0 = tog ? (C.sx < C.sy ? C.sx : C.sy)
                            : (C.sz[C.ub] < C.sz[C.vb] ? C.sz[C.ub]
                                                       : C.sz[C.vb]);
            int b = (int)(-log2(g0 > 1e-11 ? g0 : 1e-11));
            if (b < 0) b = 0;
            if (b > 39) b = 39;
            double inner = 0;
            for (long r = 0; r < fs; r++) {
                st_copy(&A, &C); st_copy(&B, &C);
                for (int s = 0; s < d; s++) step(&A);
                for (int s = 0; s < d; s++) step(&B);
                inner += eps_st(&A) * eps_st(&B);
            }
            acc += inner / fs;
            bacc[tog][b] += inner / fs;  bcnt[tog][b]++;
        }
        printf("SM psibar d=%5d cs=%ld fs=%ld burn=%d: psi=%.3e "
               "psi*d=%.4f (se_floor~%.1e)\n", d, cs, fs, burn,
               acc / cs, acc / cs * d, 1.0 / sqrt((double)cs * fs));
        printf("  gap bins: [state] b mass psi|bin contrib\n");
        for (int tg = 0; tg < 2; tg++)
            for (int b = 0; b < 40; b++)
                if (bcnt[tg][b])
                    printf("  [%s] %2d  %.5f  %+.4f  %+.3e\n",
                           tg ? "tog " : "apart", b,
                           (double)bcnt[tg][b] / cs,
                           bacc[tg][b] / bcnt[tg][b],
                           bacc[tg][b] / cs);
        fflush(stdout);
    } else if (!strcmp(mode, "iota")) {
        /* iota x y z t trials [seed]: test the label-flip symmetry.
         * Start A: together, sides (x,y), one env block z (x+y+z=1).
         * Start B: apart, blocks x (u), y (v), env z.
         * Claim: P_A(tog at t) + P_B(tog at t) = 1 for every t,
         * and the geometry marginals agree.  Prints both + sum. */
        double x = argc > 2 ? atof(argv[2]) : 0.2;
        double y = argc > 3 ? atof(argv[3]) : 0.5;
        double z = 1.0 - x - y;
        int tmax = argc > 4 ? atoi(argv[4]) : 8;
        long trials = argc > 5 ? atol(argv[5]) : 1000000;
        if (argc > 6) {
            seed = (unsigned)atoi(argv[6]);
            rng_s[0] = seed * 2654435761ULL + 1;
            rng_s[1] = seed ^ 0x9e3779b97f4a7c15ULL;
            for (int i = 0; i < 64; i++) rng_next();
        }
        St S;
        st_init(&S, 2 * tmax + 64);
        for (int t = 1; t <= tmax; t <<= 1) {
            double pa = 0, pb = 0, ga = 0, gb = 0;
            for (long r = 0; r < trials; r++) {
                /* A: together */
                S.nb = z > 1e-12 ? 2 : 1;
                S.sz[0] = x + y; S.ub = S.vb = 0; S.sx = x; S.sy = y;
                if (z > 1e-12) S.sz[1] = z;
                for (int s = 0; s < t; s++) step(&S);
                pa += (eps_st(&S) + 1) / 2;
                ga += (S.ub == S.vb) ? S.sx : S.sz[S.ub];
                /* B: apart */
                S.nb = z > 1e-12 ? 3 : 2;
                S.sz[0] = x; S.sz[1] = y; S.ub = 0; S.vb = 1;
                if (z > 1e-12) S.sz[2] = z;
                for (int s = 0; s < t; s++) step(&S);
                pb += (eps_st(&S) + 1) / 2;
                gb += (S.ub == S.vb) ? S.sx : S.sz[S.ub];
            }
            printf("iota x=%.2f y=%.2f z=%.2f t=%4d: P_A(tog)=%.5f "
                   "P_B(tog)=%.5f  sum=%.5f (want 1; se %.5f)  "
                   "Eu_A=%.5f Eu_B=%.5f\n", x, y, z, t,
                   pa / trials, pb / trials, (pa + pb) / trials,
                   sqrt(0.5 / trials), ga / trials, gb / trials);
        }
    } else if (!strcmp(mode, "gapdist")) {
        /* gapdist burn cs [seed]: stationary law of the marked state:
         * P(together), and histogram of gap g (log2 bins), plus the
         * together-block total size and x,y marginals */
        int burn = argc > 2 ? atoi(argv[2]) : 8192;
        long cs = argc > 3 ? atol(argv[3]) : 100000;
        if (argc > 4) {
            seed = (unsigned)atoi(argv[4]);
            rng_s[0] = seed * 2654435761ULL + 1;
            rng_s[1] = seed ^ 0x9e3779b97f4a7c15ULL;
            for (int i = 0; i < 64; i++) rng_next();
        }
        St C;
        st_init(&C, 2 * burn + 64);
        long ntog = 0;
        long hist[2][40] = {{0}};
        double esz = 0;
        for (long t = 0; t < cs; t++) {
            C.nb = 1; C.sz[0] = 1.0; C.ub = C.vb = 0;
            C.sx = 0.5; C.sy = 0.5;
            for (int s = 0; s < burn; s++) step(&C);
            int tog = (C.ub == C.vb);
            ntog += tog;
            double g = tog ? (C.sx < C.sy ? C.sx : C.sy)
                           : (C.sz[C.ub] < C.sz[C.vb] ? C.sz[C.ub]
                                                      : C.sz[C.vb]);
            esz += tog ? C.sz[C.ub] : 0;
            int b = (int)(-log2(g > 1e-11 ? g : 1e-11));
            if (b < 0) b = 0;
            if (b > 39) b = 39;
            hist[tog][b]++;
        }
        printf("SM gapdist burn=%d cs=%ld: P(tog)=%.4f  "
               "E[sz|tog]=%.4f\n", burn, cs, (double)ntog / cs,
               esz / (ntog ? ntog : 1));
        printf("  b  P(apart,b)  P(tog,b)   (g in 2^-(b+1)..2^-b)\n");
        for (int b = 0; b < 40; b++)
            if (hist[0][b] || hist[1][b])
                printf("  %2d  %.5f    %.5f\n", b,
                       (double)hist[0][b] / cs, (double)hist[1][b] / cs);
    } else if (!strcmp(mode, "one")) {
        int m = argc > 2 ? atoi(argv[2]) : 16;
        long trials = argc > 3 ? atol(argv[3]) : 1000000;
        St S;
        st_init(&S, 2 * m + 64);
        double s = 0;
        for (long t = 0; t < trials; t++) {
            S.nb = 1; S.sz[0] = 1.0; S.ub = S.vb = 0;
            S.sx = 0.5; S.sy = 0.5;
            for (int i = 0; i < m; i++) step(&S);
            s += eps_st(&S);
        }
        printf("SM one m=%4d trials=%ld: Eeps=%+.6f (se %.6f)\n",
               m, trials, s / trials, 1.0 / sqrt((double)trials));
    } else if (!strcmp(mode, "traj")) {
        /* block-count / Sum p_i^2 diagnostics */
        int d = argc > 2 ? atoi(argv[2]) : 1000;
        long trials = argc > 3 ? atol(argv[3]) : 100;
        St S;
        st_init(&S, 2 * d + 64);
        double nb = 0, p2 = 0;
        for (long t = 0; t < trials; t++) {
            S.nb = 1; S.sz[0] = 1.0; S.ub = S.vb = 0;
            S.sx = 0.5; S.sy = 0.5;
            for (int i = 0; i < d; i++) step(&S);
            nb += S.nb;
            double q = 0;
            for (int i = 0; i < S.nb; i++) q += S.sz[i] * S.sz[i];
            p2 += q;
        }
        printf("SM traj d=%d: E#blocks=%.2f  E sum p^2=%.4f\n",
               d, nb / trials, p2 / trials);
    } else {
        fprintf(stderr, "unknown mode %s\n", mode);
        return 1;
    }
    return 0;
}
