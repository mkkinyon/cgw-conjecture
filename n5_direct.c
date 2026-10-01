/* Exact N5(sigma) at n = 8 (fits 64-bit cell masks): number of ordered
   triples (tau, eta, theta) of B2-avoiding permutations, pairwise
   cell-disjoint.  Method: enumerate avoiders as cell bitmasks; for each
   tau collect the disjoint sublist L_tau; count ordered disjoint pairs
   within L_tau.  Also prints N0, N4 (ordered disjoint pairs) as
   cross-checks against known exact values.

   Usage: ./n5_direct 8      (both types; ~1 min)                     */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int n;
static int sigma[16];
static uint64_t masks[200000];
static int nmask;

static int perm[16];
static void dfs(int row, unsigned usedcols, uint64_t m) {
    if (row == n) { masks[nmask++] = m; return; }
    for (int c = 0; c < n; c++) {
        if (c == row || c == sigma[row]) continue;
        if (usedcols & (1u << c)) continue;
        dfs(row + 1, usedcols | (1u << c), m | ((uint64_t)1 << (row * n + c)));
    }
}

static void make_sigma_type(int type) {
    /* type 0: single n-cycle; type 1: (2, n-2) */
    if (type == 0) {
        for (int i = 0; i < n; i++) sigma[i] = (i + 1) % n;
    } else {
        sigma[0] = 1; sigma[1] = 0;
        for (int i = 2; i < n; i++) sigma[i] = 2 + (i - 2 + 1) % (n - 2);
    }
}

int main(int argc, char **argv) {
    n = (argc > 1) ? atoi(argv[1]) : 8;
    if (n * n > 64) { fprintf(stderr, "n too large for 64-bit masks\n"); return 1; }
    uint64_t *L = malloc(sizeof(uint64_t) * 200000);
    for (int type = 0; type < 2; type++) {
        make_sigma_type(type);
        nmask = 0;
        dfs(0, 0, 0);
        long long N0 = nmask;
        long long N4 = 0, N5 = 0;
        for (int a = 0; a < nmask; a++) {
            uint64_t ma = masks[a];
            int nl = 0;
            for (int b = 0; b < nmask; b++)
                if (!(ma & masks[b])) L[nl++] = masks[b];
            N4 += nl;
            long long pairs = 0;
            for (int i = 0; i < nl; i++) {
                uint64_t mi = L[i];
                for (int j = i + 1; j < nl; j++)
                    if (!(mi & L[j])) pairs++;
            }
            N5 += 2 * pairs; /* ordered (eta, theta) */
        }
        printf("n=%d type=%s: N0=%lld  N4=%lld  N5=%lld\n",
               n, type ? "(2,n-2)" : "(n)", N0, N4, N5);
        fflush(stdout);
    }
    free(L);
    return 0;
}
