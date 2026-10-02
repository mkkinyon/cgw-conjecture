/* s13_comp7.c -- exact component enumeration of the trade/turn chain on S(R), n = 7
 * (session 13, paper/fairness_question.tex).
 *
 * usage: ./s13_comp7 row1 row2 j jp
 *   row1,row2: the two fixed rows as digit strings over 0..6 (e.g. 0123456 1234560)
 *   j,jp:      the two distinguished columns (0-based)
 *
 * Enumerates every completion L of R (rows 2..6 free), computes the frame pi
 * (L[pi(r)][jp] = L[r][j]), the bit beta = [0 ~ 1 in pi], and runs union-find over
 *   A: trades at (x,y), x<y in 2..6, along the rho_{x,y}-cycle through jp (always),
 *      and turns of pi-cycles avoiding rows 0,1;
 *   B: A plus trades along the rho_{x,y}-cycle through j when (x,y) is flippable
 *      (j not on the cycle through jp; otherwise the two cycles coincide).
 * Reports the components of A and of B with size, p(K) = fraction with beta=1,
 * and the distribution of (row parity, column parity, symbol parity, sign pi).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifndef N
#define N 7
#endif
typedef struct { uint64_t hi, lo; } key_t_;

static int R[2][N], J, JP;
static key_t_ *keys; static uint32_t nk = 0, cap = 0;
static uint32_t *parA, *parB;

/* ---------- enumeration of completions ---------- */
static int L[N][N], colused[N];

static void pack(int M[N][N], uint64_t *hi, uint64_t *lo) {
    unsigned __int128 k = 0;
    for (int r = 2; r < N - 1; r++) for (int q = 0; q < N; q++) k = k * 8 + M[r][q];
    *hi = (uint64_t)(k >> 64); *lo = (uint64_t)k;
}
static void store(void) {
    uint64_t hi, lo; pack(L, &hi, &lo);
    if (nk == cap) { cap = cap ? 2 * cap : (1u << 20); keys = realloc(keys, cap * sizeof(key_t_)); }
    keys[nk].hi = hi; keys[nk].lo = lo; nk++;
}
static void rec(int r, int q, int rowused) {
    if (r == N) { store(); return; }
    if (q == N) { rec(r + 1, 0, 0); return; }
    int avail = ((1 << N) - 1) & ~rowused & ~colused[q];
    while (avail) {
        int s = __builtin_ctz(avail); avail &= avail - 1;
        L[r][q] = s; colused[q] |= 1 << s;
        rec(r, q + 1, rowused | (1 << s));
        colused[q] &= ~(1 << s);
    }
}
static int cmpkey(const void *a, const void *b) {
    const key_t_ *x = a, *y = b;
    if (x->hi != y->hi) return x->hi < y->hi ? -1 : 1;
    if (x->lo != y->lo) return x->lo < y->lo ? -1 : 1;
    return 0;
}
static uint32_t lookup(uint64_t hi, uint64_t lo) {
    uint32_t a = 0, b = nk;
    while (a < b) {
        uint32_t m = (a + b) / 2;
        if (keys[m].hi < hi || (keys[m].hi == hi && keys[m].lo < lo)) a = m + 1; else b = m;
    }
    if (a >= nk || keys[a].hi != hi || keys[a].lo != lo) { fprintf(stderr, "lookup failed\n"); exit(1); }
    return a;
}
static void decode(uint32_t i, int M[N][N]) {
    unsigned __int128 k = ((unsigned __int128)keys[i].hi << 64) | keys[i].lo;
    for (int r = N - 2; r >= 2; r--) for (int q = N - 1; q >= 0; q--) { M[r][q] = (int)(k & 7); k >>= 3; }
    for (int q = 0; q < N; q++) {
        M[0][q] = R[0][q]; M[1][q] = R[1][q];
        int used = 0; for (int r = 0; r < N - 1; r++) used |= 1 << M[r][q];
        M[N - 1][q] = __builtin_ctz(((1 << N) - 1) & ~used);
    }
}
static uint32_t encode(int M[N][N]) {
    uint64_t hi, lo; pack(M, &hi, &lo);
    return lookup(hi, lo);
}

/* ---------- union-find ---------- */
static uint32_t find(uint32_t *par, uint32_t a) {
    while (par[a] != a) { par[a] = par[par[a]]; a = par[a]; }
    return a;
}
static void unite(uint32_t *par, uint32_t a, uint32_t b) {
    a = find(par, a); b = find(par, b); if (a != b) par[a] = b;
}

/* ---------- square utilities ---------- */
static void frame(int M[N][N], int pi[N]) {
    int rowof[N]; for (int r = 0; r < N; r++) rowof[M[r][JP]] = r;
    for (int r = 0; r < N; r++) pi[r] = rowof[M[r][J]];
}
static int sgnperm(const int *p, int n) {
    int seen = 0, s = 1;
    for (int i = 0; i < n; i++) if (!(seen >> i & 1)) {
        int l = 0, k = i; while (!(seen >> k & 1)) { seen |= 1 << k; k = p[k]; l++; }
        if (l % 2 == 0) s = -s;
    }
    return s;
}
static int same_cycle(const int *pi, int u, int v) {
    int x = u; for (;;) { x = pi[x]; if (x == v) return 1; if (x == u) return 0; }
}
/* columns of the rho_{x,y}-cycle through q0; returns length, fills cyc[] */
static int rho_cycle(int M[N][N], int x, int y, int q0, int *cyc) {
    int pos[N]; for (int q = 0; q < N; q++) pos[M[y][q]] = q;
    int len = 0, q = q0;
    do { cyc[len++] = q; q = pos[M[x][q]]; } while (q != q0);
    return len;
}

int main(int argc, char **argv) {
    if (argc != 5) { fprintf(stderr, "usage: %s row1 row2 j jp\n", argv[0]); return 1; }
    for (int q = 0; q < N; q++) { R[0][q] = argv[1][q] - '0'; R[1][q] = argv[2][q] - '0'; }
    J = atoi(argv[3]); JP = atoi(argv[4]);
    for (int q = 0; q < N; q++) { L[0][q] = R[0][q]; L[1][q] = R[1][q]; colused[q] = (1 << R[0][q]) | (1 << R[1][q]); }
    /* sigma (column permutation of R: R[1][sigma(q)] = R[0][q]) and hypothesis (1) */
    int pos1[N]; for (int q = 0; q < N; q++) pos1[R[1][q]] = q;
    int sigma[N]; for (int q = 0; q < N; q++) sigma[q] = pos1[R[0][q]];
    int ok = (R[0][J] != R[1][JP]) && (R[0][JP] != R[1][J]);
    printf("R = %s / %s, j=%d jp=%d, sigma(j)=%d sigma(jp)=%d, hypothesis (1): %s\n",
           argv[1], argv[2], J, JP, sigma[J], sigma[JP], ok ? "holds" : "FAILS");
    { int d = 0, q = J; while (q != JP && d <= N) { q = sigma[q]; d++; }
      if (d <= N) printf("  jp = sigma^%d(j) (same sigma-cycle)\n", d); else printf("  j, jp in different sigma-cycles\n"); }

    rec(2, 0, 0);
    qsort(keys, nk, sizeof(key_t_), cmpkey);
    printf("|S(R)| = %u completions\n", nk); fflush(stdout);

    parA = malloc(nk * sizeof(uint32_t)); parB = malloc(nk * sizeof(uint32_t));
    uint8_t *beta = malloc(nk);
    for (uint32_t i = 0; i < nk; i++) parA[i] = parB[i] = i;

    long nflip = 0, nnon = 0, nturn = 0, nturn_odd = 0;
#define NOFF 4
    long lad_in[NOFF][2] = {{0}}, lad_out[NOFF][2] = {{0}}, lad_none[NOFF][2] = {{0}}, lad_K0hist[N + 1] = {0}, lad_none_K0[N + 1] = {0};
    for (uint32_t i = 0; i < nk; i++) {
        int M[N][N], pi[N], cyc[N];
        decode(i, M); frame(M, pi);
        beta[i] = same_cycle(pi, 0, 1);
        for (int c = 0; c < NOFF; c++) {   /* diagonal families D_c: pairs (pi^-(c+i)(1), pi^-i(2)), 1 <= i < K_c */
            int inv[N]; for (int r = 0; r < N; r++) inv[pi[r]] = r;
            int c0 = 0, r = 0; do { c0++; r = pi[r]; } while (r != 0);
            int c1 = 0; r = 1; do { c1++; r = pi[r]; } while (r != 1);
            int Kc;
            if (beta[i]) { int d01 = 0; r = 0; while (r != 1) { r = pi[r]; d01++; } int d10 = c0 - d01; Kc = (d10 - c) < d01 ? (d10 - c) : d01; }
            else Kc = (c0 - c) < c1 ? (c0 - c) : c1;
            if (Kc < 2) { lad_out[c][beta[i]]++; continue; }   /* outside the domain D_c */
            int anyflip = 0, x = 0, y = 1;
            for (int t = 0; t < c; t++) x = inv[x];
            for (int k = 1; k < Kc; k++) {
                x = inv[x]; y = inv[y];
                int len = rho_cycle(M, x, y, JP, cyc), hasj = 0;
                for (int t = 0; t < len; t++) if (cyc[t] == J) hasj = 1;
                if (!hasj) { anyflip = 1; break; }
            }
            lad_in[c][beta[i]]++; if (!anyflip) lad_none[c][beta[i]]++;
            if (c == 0) { lad_K0hist[Kc]++; if (!anyflip) lad_none_K0[Kc]++; }
        }
    }
    printf("moves: %ld flippable trades, %ld non-flippable, %ld turns (%ld of odd cycles)\n", nflip, nnon, nturn, nturn_odd);
    double pb = 0; for (uint32_t i = 0; i < nk; i++) pb += beta[i]; pb /= nk;
    printf("P[beta=1 | R] = %.6f\n", pb);
    for (int c = 0; c < NOFF; c++)
        printf("diagonal c=%d: P[B,D]-P[A,D] = %+.6f ; P[B,D,Flip=0]-P[A,D,Flip=0] = %+.6f ; P[D^c] = %.4f (B:%ld A:%ld) ; P[D,Flip=0] = %.4f\n", c,
               (double)(lad_in[c][1] - lad_in[c][0]) / nk, (double)(lad_none[c][1] - lad_none[c][0]) / nk,
               (double)(lad_out[c][0] + lad_out[c][1]) / nk, lad_out[c][1], lad_out[c][0], (double)(lad_none[c][0] + lad_none[c][1]) / nk);
    printf("  K0 histogram (K0: count, Flip=0 count):"); for (int k = 0; k <= N; k++) if (lad_K0hist[k]) printf(" %d:%ld,%ld", k, lad_K0hist[k], lad_none_K0[k]); printf("\n");

    for (int which = 0; which < 2; which++) {
        uint32_t *par = which ? parB : parA;
        /* roots */
        uint32_t *root = malloc(nk * sizeof(uint32_t));
        for (uint32_t i = 0; i < nk; i++) root[i] = find(par, i);
        /* component ids for all roots */
        int32_t *cid = malloc(nk * sizeof(int32_t)); int nr = 0;
        for (uint32_t i = 0; i < nk; i++) cid[i] = -1;
        for (uint32_t i = 0; i < nk; i++) if (root[i] == i) cid[i] = nr++;
        long *sz = calloc(nr, sizeof(long)), *nb = calloc(nr, sizeof(long)), (*hist)[16] = calloc(nr, sizeof *hist);
        long *rowpar1 = calloc(nr, sizeof(long));
        for (uint32_t i = 0; i < nk; i++) {
            int c = cid[root[i]];
            int M[N][N], pi[N]; decode(i, M); frame(M, pi);
            int rp = 0, cp = 0, sp = 0, col[N], pos[N][N];
            for (int r = 0; r < N; r++) { if (sgnperm(M[r], N) < 0) rp ^= 1; for (int q = 0; q < N; q++) pos[r][M[r][q]] = q; }
            for (int q = 0; q < N; q++) { for (int r = 0; r < N; r++) col[r] = M[r][q]; if (sgnperm(col, N) < 0) cp ^= 1; }
            for (int s = 0; s < N; s++) { for (int r = 0; r < N; r++) col[r] = pos[r][s]; if (sgnperm(col, N) < 0) sp ^= 1; }
            int spi = sgnperm(pi, N) < 0;
            sz[c]++; nb[c] += beta[i]; hist[c][rp * 8 + cp * 4 + sp * 2 + spi]++; rowpar1[c] += rp;
        }
        int *ord = malloc(nr * sizeof(int)); for (int t = 0; t < nr; t++) ord[t] = t;
        /* simple selection of the 12 largest; the rest summarised */
        for (int a = 0; a < nr && a < 12; a++) for (int b = a + 1; b < nr; b++) if (sz[ord[b]] > sz[ord[a]]) { int t = ord[a]; ord[a] = ord[b]; ord[b] = t; }
        double lhs = 0, lhs_small = 0; long n_small = 0, sz_small = 0, maxsmall = 0;
        for (int t = 0; t < nr; t++) { double p = (double)nb[t] / sz[t]; lhs += sz[t] * (2 * p - 1) * (2 * p - 1); }
        for (int t = 12; t < nr; t++) { int c = ord[t]; double p = (double)nb[c] / sz[c]; lhs_small += sz[c] * (2 * p - 1) * (2 * p - 1); n_small++; sz_small += sz[c]; if (sz[c] > maxsmall) maxsmall = sz[c]; }
        printf("\n== move set %s: %d components, LHS of Question = %.6f\n", which ? "B (both trades + turns)" : "A (jp-trades + turns)", nr, lhs / nk);
        for (int t = 0; t < nr && t < 12; t++) {
            int c = ord[t];
            printf("  |K| = %8ld (%.5f of S(R))  p(K) = %.4f  row-parity-1 fraction = %.3f  (row,col,sym,sgnpi) bins:",
                   sz[c], (double)sz[c] / nk, (double)nb[c] / sz[c], (double)rowpar1[c] / sz[c]);
            for (int h = 0; h < 16; h++) if (hist[c][h]) printf(" %d%d%d%d:%ld", h >> 3 & 1, h >> 2 & 1, h >> 1 & 1, h & 1, hist[c][h]);
            printf("\n");
        }
        if (which == 0 && getenv("DUMP")) {
            int *done = calloc(nr, sizeof(int));
            for (uint32_t i = 0; i < nk; i++) {
                int c = cid[root[i]]; if (sz[c] > 2000 || done[c]) continue; done[c] = 1;
                int M[N][N], pi[N], cyc[N]; decode(i, M); frame(M, pi);
                printf("DUMP size=%ld p=%.3f square=", sz[c], (double)nb[c] / sz[c]);
                for (int r = 0; r < N; r++) { for (int q = 0; q < N; q++) printf("%d", M[r][q]); printf(r < N - 1 ? "/" : ""); }
                printf(" pi="); for (int r = 0; r < N; r++) printf("%d", pi[r]);
                int nf = 0; printf(" rho-cycle-lengths-through-jp:");
                for (int x = 0; x < N; x++) for (int y = x + 1; y < N; y++) {
                    int len = rho_cycle(M, x, y, JP, cyc), hasj = 0;
                    for (int t = 0; t < len; t++) if (cyc[t] == J) hasj = 1;
                    printf("%d%s", len, hasj ? "" : "f"); if (!hasj && x >= 2) nf++;
                }
                printf(" flippable(x,y>=2)=%d\n", nf);
            }
            free(done);
        }
        if (n_small) printf("  + %ld smaller components, total size %ld (%.5f of S(R)), largest %ld, contributing %.6f to the LHS\n",
                            n_small, sz_small, (double)sz_small / nk, maxsmall, lhs_small / nk);
        free(cid); free(ord);
        free(root); free(sz); free(nb); free(hist); free(rowpar1);
    }
    return 0;
}
