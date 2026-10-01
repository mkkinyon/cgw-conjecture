/* xside.c: enumerate all Latin squares of order n (first row = identity),
   and for each derangement type lambda, each cycle length m >= 4 of sigma,
   each alpha in 2..m-2, accumulate the counts of marked pairs
   {j, sigma^alpha(j)} (j in an m-cycle) that are A-pairs / B-pairs.

   Theorem M check:  X_B / X_A  ==  2 C(lambda)/C(mu) - 1   where mu is
   lambda with one m-part split into (alpha, m-alpha).

   Usage: ./xside n     (n <= 7 practical)
   Output lines: type  m  alpha  XA  XB
   (alpha ranges 2..m/2; pairs are unordered: for alpha < m-alpha each
   unordered pair counted once via all j; for m even, alpha=m/2 each pair
   arises twice over j -- we halve those counts.)                       */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n;
static int L[9][9];
static unsigned colmask[9];
static long long acc[2][40]; /* [status][slot] */
/* slot encoding: we build a table of (type-index, m, alpha) at runtime */

#define MAXTYPES 16
static int ntypes = 0;
static int typekey[MAXTYPES]; /* encoded partition */
static long long XA[MAXTYPES][9][9], XB[MAXTYPES][9][9];

static int encode_type(int *cnt) { /* cnt[l] = #cycles of length l */
    int key = 0;
    for (int l = 2; l <= n; l++) key = key * 9 + cnt[l];
    return key;
}
static int type_index(int *cnt) {
    int key = encode_type(cnt);
    for (int i = 0; i < ntypes; i++) if (typekey[i] == key) return i;
    typekey[ntypes] = key;
    return ntypes++;
}
static void decode_print(int key) {
    int cnt[10];
    for (int l = n; l >= 2; l--) { cnt[l] = key % 9; key /= 9; }
    printf("(");
    int first = 1;
    for (int l = 2; l <= n; l++)
        for (int c = 0; c < cnt[l]; c++) {
            if (!first) printf(",");
            printf("%d", l); first = 0;
        }
    printf(")");
}

static void process(void) {
    int sig[9], cyclen[9], cycid[9], cid = 0;
    for (int j = 0; j < n; j++) sig[j] = L[1][j];
    int seen[9] = {0};
    int cyccnt[10] = {0};
    for (int j = 0; j < n; j++) if (!seen[j]) {
        int len = 0, x = j;
        while (!seen[x]) { seen[x] = 1; cycid[x] = cid; x = sig[x]; len++; }
        cyclen[cid++] = len;
        cyccnt[len]++;
    }
    int ti = type_index(cyccnt);
    /* position of symbol s in column jp: pos[s] per column computed lazily */
    for (int j0 = 0; j0 < n; j0++) {
        int m = cyclen[cycid[j0]];
        if (m < 4) continue;
        for (int alpha = 2; alpha <= m - 2; alpha++) {
            if (alpha > m - alpha) break;
            int jp = j0;
            for (int t = 0; t < alpha; t++) jp = sig[jp];
            /* status of pair {j0, jp}: rows 0,1 same column cycle? */
            int pos_jp[9];
            for (int i = 0; i < n; i++) pos_jp[L[i][jp]] = i;
            int i = 0, same = 0;
            for (;;) {
                i = pos_jp[L[i][j0]];
                if (i == 0) break;
                if (i == 1) { same = 1; break; }
            }
            if (same) XB[ti][m][alpha]++;
            else XA[ti][m][alpha]++;
        }
    }
}

static void rec(int r, int c, unsigned used) {
    if (r == n) { process(); return; }
    if (c == n) { rec(r + 1, 0, 0); return; }
    unsigned avail = ((1u << n) - 1) & ~used & ~colmask[c];
    /* row 1 must be a derangement: L[1][c] != c */
    while (avail) {
        unsigned b = avail & (-avail);
        avail ^= b;
        int s = __builtin_ctz(b);
        if (r == 1 && s == c) continue;
        L[r][c] = s;
        colmask[c] |= b;
        rec(r, c + 1, used | b);
        colmask[c] ^= b;
    }
}

int main(int argc, char **argv) {
    n = argc > 1 ? atoi(argv[1]) : 6;
    for (int c = 0; c < n; c++) { L[0][c] = c; colmask[c] = 1u << c; }
    rec(1, 0, 0);
    for (int ti = 0; ti < ntypes; ti++)
        for (int m = 4; m <= n; m++)
            for (int a = 2; a + a <= m; a++)
                if (XA[ti][m][a] + XB[ti][m][a]) {
                    long long xa = XA[ti][m][a], xb = XB[ti][m][a];
                    if (a + a == m) { xa /= 2; xb /= 2; }
                    decode_print(typekey[ti]);
                    printf("  m=%d alpha=%d  XA=%lld XB=%lld  XB/XA=%.9f\n",
                           m, a, xa, xb, (double)xb / (double)xa);
                }
    return 0;
}
