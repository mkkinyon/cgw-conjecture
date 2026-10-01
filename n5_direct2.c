/* Exact N5(sigma) for n = 8, 9 (two-word cell masks, OpenMP).
   N5 = # ordered triples of B2-avoiding perms, pairwise cell-disjoint
      = 3 * sum_tau #{ordered pairs (eta,theta) in L_tau, idx > tau, disjoint}
   where L_tau = avoiders with index > tau disjoint from tau.
   Usage: ./n5_direct2 9        (prints both types; n=9 ~ 1 h on 2 cores) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static int n;
static int sigma[16];
typedef struct { uint64_t lo, hi; } m2;
static m2 *masks;
static int nmask, cap;

static void dfs(int row, unsigned usedcols, uint64_t lo, uint64_t hi) {
    if (row == n) {
        if (nmask == cap) { cap *= 2; masks = realloc(masks, cap * sizeof(m2)); }
        masks[nmask].lo = lo; masks[nmask].hi = hi; nmask++;
        return;
    }
    for (int c = 0; c < n; c++) {
        if (c == row || c == sigma[row]) continue;
        if (usedcols & (1u << c)) continue;
        int b = row * n + c;
        dfs(row + 1, usedcols | (1u << c),
            b < 64 ? lo | ((uint64_t)1 << b) : lo,
            b < 64 ? hi : hi | ((uint64_t)1 << (b - 64)));
    }
}

static void make_sigma_type(int type) {
    if (type == 0) { for (int i = 0; i < n; i++) sigma[i] = (i + 1) % n; }
    else {
        sigma[0] = 1; sigma[1] = 0;
        for (int i = 2; i < n; i++) sigma[i] = 2 + (i - 1) % (n - 2);
    }
}

int main(int argc, char **argv) {
    n = (argc > 1) ? atoi(argv[1]) : 9;
    if (n * n > 128) { fprintf(stderr, "n too large\n"); return 1; }
    for (int type = 0; type < 2; type++) {
        make_sigma_type(type);
        cap = 1 << 16; masks = malloc(cap * sizeof(m2)); nmask = 0;
        dfs(0, 0, 0, 0);
        long long N0 = nmask, N4 = 0, N5part = 0;
#pragma omp parallel for schedule(dynamic, 16) reduction(+:N4, N5part)
        for (int a = 0; a < nmask; a++) {
            uint64_t alo = masks[a].lo, ahi = masks[a].hi;
            /* full L for N4; restricted (idx > a) for N5 */
            int nl = 0;
            m2 *L = malloc((nmask - a) * sizeof(m2));
            for (int b = 0; b < nmask; b++) {
                if ((alo & masks[b].lo) | (ahi & masks[b].hi)) continue;
                N4++;
                if (b > a) L[nl++] = masks[b];
            }
            long long pairs = 0;
            for (int i = 0; i < nl; i++) {
                uint64_t ilo = L[i].lo, ihi = L[i].hi;
                for (int j = i + 1; j < nl; j++)
                    if (!((ilo & L[j].lo) | (ihi & L[j].hi))) pairs++;
            }
            N5part += 2 * pairs;
            free(L);
        }
        printf("n=%d type=%s: N0=%lld  N4=%lld  N5=%lld\n",
               n, type ? "(2,n-2)" : "(n)", N0, N4, 3 * N5part);
        fflush(stdout);
        free(masks);
    }
    return 0;
}
