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
static long NFLIP=0;
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
                NFLIP++;
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
            NFLIP++;
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
