/* Layer-5 anatomy on board0 via a global disjointness bit-matrix.
   Same outputs as anatomy5.c but O(N0^2/64) word-ops per tau row:
   D[a] = bitset over b of "masks[a] and masks[b] cell-disjoint";
   L_tau = D[tau];  #ordered pairs in L_tau = sum_{i in L} popcount(D[i]&L).
   Memory: N0^2/8 bytes (378 MB at n=9).  Usage: ./anatomy5b 9  (~minutes)
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static int n;
static int forb1[16], forb2[16];
typedef struct { uint64_t lo, hi; } m2;
static m2 *masks;
static int nmask, cap;
static int *pincode;
static int D1r, D1c, D2r, D2c, E1r, E1c, E2r, E2c;

static void dfs(int row, unsigned usedcols, uint64_t lo, uint64_t hi) {
    if (row == n) {
        if (nmask == cap) { cap *= 2; masks = realloc(masks, cap * sizeof(m2)); }
        masks[nmask].lo = lo; masks[nmask].hi = hi; nmask++;
        return;
    }
    for (int c = 0; c < n; c++) {
        if (c == forb1[row] || c == forb2[row]) continue;
        if (usedcols & (1u << c)) continue;
        int b = row * n + c;
        dfs(row + 1, usedcols | (1u << c),
            b < 64 ? lo | ((uint64_t)1 << b) : lo,
            b < 64 ? hi : hi | ((uint64_t)1 << (b - 64)));
    }
}

static inline int bitat(m2 m, int r, int c) {
    int b = r * n + c;
    return b < 64 ? (int)((m.lo >> b) & 1) : (int)((m.hi >> (b - 64)) & 1);
}

int main(int argc, char **argv) {
    n = (argc > 1) ? atoi(argv[1]) : 9;
    for (int i = 0; i < n; i++) { forb1[i] = i; forb2[i] = (i + 1) % n; }
    forb2[1] = 1; forb2[n - 1] = n - 1;
    D1r = 1; D1c = 2; D2r = n - 1; D2c = 0;
    E1r = 1; E1c = 0; E2r = n - 1; E2c = 2;
    cap = 1 << 16; masks = malloc(cap * sizeof(m2)); nmask = 0;
    dfs(0, 0, 0, 0);
    long long N0 = nmask;
    pincode = malloc(nmask * sizeof(int));
    for (int a = 0; a < nmask; a++) {
        int pc = 0;
        if (bitat(masks[a], D1r, D1c)) pc |= 1;
        if (bitat(masks[a], D2r, D2c)) pc |= 2;
        if (bitat(masks[a], E1r, E1c)) pc |= 4;
        if (bitat(masks[a], E2r, E2c)) pc |= 8;
        pincode[a] = pc;
    }
    int W = (nmask + 63) / 64;
    uint64_t *D = calloc((size_t)nmask * W, sizeof(uint64_t));
    if (!D) { fprintf(stderr, "alloc failed\n"); return 1; }
#pragma omp parallel for schedule(dynamic, 64)
    for (int a = 0; a < nmask; a++) {
        uint64_t *row = D + (size_t)a * W;
        uint64_t alo = masks[a].lo, ahi = masks[a].hi;
        for (int b = 0; b < nmask; b++)
            if (!((alo & masks[b].lo) | (ahi & masks[b].hi)))
                row[b >> 6] |= (uint64_t)1 << (b & 63);
    }
    printf("n=%d  N0=%lld  (matrix built)\n", n, N0); fflush(stdout);

    long long P05 = 0, G5d1 = 0, G5d2 = 0, G5e1 = 0, G5e2 = 0;
    long long G25d = 0, G25e = 0, H5d = 0, H5e = 0;
#pragma omp parallel for schedule(dynamic, 16) reduction(+:P05,G5d1,G5d2,G5e1,G5e2,G25d,G25e,H5d,H5e)
    for (int a = 0; a < nmask; a++) {
        uint64_t *L = D + (size_t)a * W;
        int pc = pincode[a];
        long long cnt = 0, cntD2 = 0, cntE2 = 0;
        for (int w = 0; w < W; w++) {
            uint64_t bits = L[w];
            while (bits) {
                int i = w * 64 + __builtin_ctzll(bits);
                bits &= bits - 1;
                uint64_t *Di = D + (size_t)i * W;
                long long di = 0;
                for (int v = 0; v < W; v++)
                    di += __builtin_popcountll(Di[v] & L[v]);
                cnt += di;
                if (pincode[i] & 2) cntD2 += di;
                if (pincode[i] & 8) cntE2 += di;
            }
        }
        P05 += cnt;
        if (pc & 1) { G5d1 += cnt; H5d += cntD2; }
        if (pc & 2) G5d2 += cnt;
        if (pc & 4) { G5e1 += cnt; H5e += cntE2; }
        if (pc & 8) G5e2 += cnt;
        if ((pc & 3) == 3) G25d += cnt;
        if ((pc & 12) == 12) G25e += cnt;
    }
    long long BrG5 = G5d1 + G5d2 - G5e1 - G5e2;
    printf("P05=%lld\n", P05);
    printf("G5(d1)=%lld G5(d2)=%lld G5(e1)=%lld G5(e2)=%lld\n", G5d1, G5d2, G5e1, G5e2);
    printf("BrG5=%lld\n", BrG5);
    printf("G25(d)=%lld G25(e)=%lld  H5(d)=%lld H5(e)=%lld\n", G25d, G25e, H5d, H5e);
    printf("3BrG5+3dG25+6dH5 = %lld   (should equal DeltaN5)\n",
           3 * BrG5 + 3 * (G25e - G25d) + 6 * (H5e - H5d));
    return 0;
}
