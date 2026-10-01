/* Exact N5 for arbitrary derangement cycle type at small n (<=8: fast).
   Usage: ./n5_types n l1 l2 ...   (cycle lengths summing to n, all >=2) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static int n;
static int sigma[16];
static uint64_t masks[200000];
static int nmask;

static void dfs(int row, unsigned usedcols, uint64_t m) {
    if (row == n) { masks[nmask++] = m; return; }
    for (int c = 0; c < n; c++) {
        if (c == row || c == sigma[row]) continue;
        if (usedcols & (1u << c)) continue;
        dfs(row + 1, usedcols | (1u << c), m | ((uint64_t)1 << (row * n + c)));
    }
}

int main(int argc, char **argv) {
    n = atoi(argv[1]);
    if (n * n > 64) { fprintf(stderr, "n too large\n"); return 1; }
    int pos = 0;
    for (int a = 2; a < argc; a++) {
        int l = atoi(argv[a]);
        for (int i = 0; i < l; i++) sigma[pos + i] = pos + (i + 1) % l;
        pos += l;
    }
    if (pos != n) { fprintf(stderr, "type does not sum to n\n"); return 1; }
    nmask = 0;
    dfs(0, 0, 0);
    long long N0 = nmask, N4 = 0, N5 = 0;
    uint64_t *L = malloc(sizeof(uint64_t) * nmask);
    for (int a = 0; a < nmask; a++) {
        uint64_t ma = masks[a];
        int nl = 0;
        for (int b = 0; b < nmask; b++)
            if (!(ma & masks[b])) { L[nl++] = masks[b]; N4++; }
        long long pairs = 0;
        for (int i = 0; i < nl; i++) {
            uint64_t mi = L[i];
            for (int j = i + 1; j < nl; j++)
                if (!(mi & L[j])) pairs++;
        }
        N5 += 2 * pairs;
    }
    printf("n=%d type=(", n);
    for (int a = 2; a < argc; a++) printf("%s%s", argv[a], a + 1 < argc ? "," : "");
    printf("): N3=%lld N4=%lld N5=%lld\n", N0, N4, N5);
    free(L);
    return 0;
}
