/* bias_scale.c -- k-scaling of the orbit bias for random chord systems.
 *
 * Ensemble: circle with 2k+2 slots; u at slot 0, v at slot vslot
 * (default antipodal k+1); the other 2k slots carry a uniform random
 * perfect matching (k chords).  For each instance, estimate
 *
 *   bias = P_{S uniform}[ x_S in colspace(A_S) ] - 1/2
 *
 * by Monte Carlo over S (or exhaustively if k <= 20), using the
 * verified bordered-rank formula:
 *   bit(S) = 1  iff  rank([A_S | x_S]) == rank(A_S).
 *
 * Interlacement of chords {a<b},{c<d}: a<c<b<d or c<a<d<b.
 * Separation: exactly one of u,v strictly inside (a,b).
 *
 * Usage: ./bias_scale k trials ssamples [vslot]
 * Output: per-instance biases -> summary (mean, E|.|, quantiles).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXK 192
#define W ((MAXK + 2 + 63) / 64)

typedef unsigned long long u64;

static u64 A[MAXK][W];     /* adjacency rows (bit j of row i) */
static u64 xvec[W];        /* separation vector */

static unsigned long long rng_s[2];
static inline u64 rng_next(void) {
    u64 s1 = rng_s[0], s0 = rng_s[1];
    rng_s[0] = s0;
    s1 ^= s1 << 23;
    rng_s[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
    return rng_s[1] + s0;
}
static inline double rng_u(void) { return (rng_next() >> 11) * (1.0 / 9007199254740992.0); }

/* rank over F2 of rows[0..m-1] each W words, ncols columns */
static int f2rank(u64 rows[][W], int m, int ncols) {
    int rank = 0;
    for (int col = 0; col < ncols && rank < m; col++) {
        int w = col >> 6, b = col & 63;
        int piv = -1;
        for (int r = rank; r < m; r++)
            if ((rows[r][w] >> b) & 1ULL) { piv = r; break; }
        if (piv < 0) continue;
        for (int t = 0; t < W; t++) {
            u64 tmp = rows[rank][t]; rows[rank][t] = rows[piv][t]; rows[piv][t] = tmp;
        }
        for (int r = 0; r < m; r++) {
            if (r != rank && ((rows[r][w] >> b) & 1ULL))
                for (int t = 0; t < W; t++) rows[r][t] ^= rows[rank][t];
        }
        rank++;
    }
    return rank;
}

static u64 work[MAXK][W];

/* bit for activated set given as list idx[0..s-1] */
static int bit_of(int *idx, int s, int k) {
    /* build A_S rows restricted to S columns, plus x_S in column s */
    for (int a = 0; a < s; a++) {
        memset(work[a], 0, sizeof(work[a]));
        int i = idx[a];
        for (int b = 0; b < s; b++) {
            int j = idx[b];
            if ((A[i][j >> 6] >> (j & 63)) & 1ULL)
                work[a][b >> 6] |= 1ULL << (b & 63);
        }
    }
    /* rank(A_S) */
    static u64 tmp[MAXK][W];
    memcpy(tmp, work, sizeof(u64) * s * W);
    int rA = f2rank(tmp, s, s);
    /* augmented */
    for (int a = 0; a < s; a++) {
        int i = idx[a];
        memcpy(tmp[a], work[a], sizeof(u64) * W);
        if ((xvec[i >> 6] >> (i & 63)) & 1ULL)
            tmp[a][s >> 6] |= 1ULL << (s & 63);
    }
    int rAx = f2rank(tmp, s, s + 1);
    return rA == rAx;
}

int main(int argc, char **argv) {
    int k = argc > 1 ? atoi(argv[1]) : 24;
    int trials = argc > 2 ? atoi(argv[2]) : 200;
    int ssamp = argc > 3 ? atoi(argv[3]) : 20000;
    int N = 2 * k + 2;
    int vslot = argc > 4 ? atoi(argv[4]) : (N / 2);
    unsigned seed = argc > 5 ? (unsigned)atoi(argv[5]) : 12345u;
    rng_s[0] = seed * 2654435761ULL + 1; rng_s[1] = seed ^ 0x9e3779b97f4a7c15ULL;
    for (int i = 0; i < 64; i++) rng_next();

    double *biases = malloc(sizeof(double) * trials);
    int slots[2 * MAXK + 2];
    int chA[MAXK], chB[MAXK];
    int idx[MAXK];

    for (int t = 0; t < trials; t++) {
        /* slots other than 0 (=u) and vslot (=v) */
        int m = 0;
        for (int spos = 0; spos < N; spos++)
            if (spos != 0 && spos != vslot) slots[m++] = spos;
        /* random perfect matching: shuffle then pair consecutive */
        for (int i2 = m - 1; i2 > 0; i2--) {
            int j2 = (int)(rng_u() * (i2 + 1));
            int tt = slots[i2]; slots[i2] = slots[j2]; slots[j2] = tt;
        }
        for (int c = 0; c < k; c++) {
            int a = slots[2 * c], b = slots[2 * c + 1];
            if (a > b) { int tt = a; a = b; b = tt; }
            chA[c] = a; chB[c] = b;
        }
        /* interlacement + separation */
        memset(A, 0, sizeof(A));
        memset(xvec, 0, sizeof(xvec));
        for (int i2 = 0; i2 < k; i2++) {
            for (int j2 = i2 + 1; j2 < k; j2++) {
                int a = chA[i2], b = chB[i2], c = chA[j2], d = chB[j2];
                int inter = (a < c && c < b && b < d) || (c < a && a < d && d < b);
                if (inter) {
                    A[i2][j2 >> 6] |= 1ULL << (j2 & 63);
                    A[j2][i2 >> 6] |= 1ULL << (i2 & 63);
                }
            }
            int a = chA[i2], b = chB[i2];
            int inu = (a < 0 && 0 < b);           /* u at slot 0: never inside (a,b) since a>=1 */
            int inv = (a < vslot && vslot < b);
            if (inu != inv) xvec[i2 >> 6] |= 1ULL << (i2 & 63);
        }
        /* bias */
        long hits = 0, tot = 0;
        if (k <= 20 && (1L << k) <= ssamp * 2L) {
            for (long Sb = 0; Sb < (1L << k); Sb++) {
                int s = 0;
                for (int i2 = 0; i2 < k; i2++) if ((Sb >> i2) & 1) idx[s++] = i2;
                hits += bit_of(idx, s, k);
                tot++;
            }
        } else {
            for (int rep = 0; rep < ssamp; rep++) {
                int s = 0;
                u64 r = 0;
                for (int i2 = 0; i2 < k; i2++) {
                    if ((i2 & 63) == 0) r = rng_next();
                    if ((r >> (i2 & 63)) & 1ULL) idx[s++] = i2;
                }
                hits += bit_of(idx, s, k);
                tot++;
            }
        }
        biases[t] = (double)hits / tot - 0.5;
    }
    /* summary */
    double mean = 0, meanabs = 0;
    for (int t = 0; t < trials; t++) { mean += biases[t]; meanabs += fabs(biases[t]); }
    mean /= trials; meanabs /= trials;
    /* sort |bias| */
    for (int i2 = 1; i2 < trials; i2++) {
        double key = fabs(biases[i2]); int j2 = i2 - 1;
        double keys = biases[i2];
        while (j2 >= 0 && fabs(biases[j2]) > key) { biases[j2 + 1] = biases[j2]; j2--; }
        biases[j2 + 1] = keys;
    }
    double q50 = fabs(biases[trials / 2]), q90 = fabs(biases[(int)(0.9 * trials)]);
    double mx = fabs(biases[trials - 1]);
    double se = meanabs / sqrt((double)trials);
    printf("k=%3d vslot=%3d trials=%d ssamp=%d: mean %+.5f (se~%.5f) "
           "E|b| %.5f  med %.5f  q90 %.5f  max %.5f\n",
           k, vslot, trials, ssamp, mean, se, meanabs, q50, q90, mx);
    free(biases);
    return 0;
}
