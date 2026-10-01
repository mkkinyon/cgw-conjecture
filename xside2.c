/* xside2.c: X-side counts with sigma FIXED to a canonical representative
   of a given cycle type (justified by symbol/column conjugation symmetry:
   the ratio X_B/X_A per type depends only on the type).

   Enumerate all completions (rows 3..n) of the 2xn rectangle (id, sigma0),
   and for each marked pair {j, sigma0^alpha(j)} (j in an m-cycle, m>=4,
   2 <= alpha <= m/2) accumulate A/B status counts.

   Usage: ./xside2 n p1 p2 ...   (parts of the cycle type, sum = n)
   e.g.   ./xside2 7 7        (type (7))
          ./xside2 7 2 5      (type (2,5))
   Output: per (m,alpha): XA XB XB/XA, plus total completions N(sigma).  */
#include <stdio.h>
#include <stdlib.h>

static int n, np, parts[12];
static int sig[12];          /* sigma0 */
static int L[12][12];
static unsigned colmask[12];
static long long XA[12][12], XB[12][12], NC = 0;

/* marked pair list, precomputed from sigma0 */
static int nmark;
static int mj[144], mjp[144], mm[144], ma[144];

static void process(void) {
    NC++;
    for (int t = 0; t < nmark; t++) {
        int j = mj[t], jp = mjp[t];
        int pos_jp[12];
        for (int i = 0; i < n; i++) pos_jp[L[i][jp]] = i;
        int i = 0, same = 0;
        for (;;) {
            i = pos_jp[L[i][j]];
            if (i == 0) break;
            if (i == 1) { same = 1; break; }
        }
        if (same) XB[mm[t]][ma[t]]++; else XA[mm[t]][ma[t]]++;
    }
}

static void rec(int r, int c, unsigned used) {
    if (r == n) { process(); return; }
    if (c == n) { rec(r + 1, 0, 0); return; }
    unsigned avail = ((1u << n) - 1) & ~used & ~colmask[c];
    while (avail) {
        unsigned b = avail & (-avail);
        avail ^= b;
        L[r][c] = __builtin_ctz(b);
        colmask[c] |= b;
        rec(r, c + 1, used | b);
        colmask[c] ^= b;
    }
}

int main(int argc, char **argv) {
    n = atoi(argv[1]);
    np = argc - 2;
    int s = 0;
    for (int i = 0; i < np; i++) { parts[i] = atoi(argv[i + 2]); s += parts[i]; }
    if (s != n) { fprintf(stderr, "parts must sum to n\n"); return 1; }
    /* canonical sigma0: consecutive blocks, each a cyclic shift */
    int base = 0;
    for (int i = 0; i < np; i++) {
        int p = parts[i];
        for (int k = 0; k < p; k++) sig[base + k] = base + (k + 1) % p;
        base += p;
    }
    /* marked pairs */
    nmark = 0;
    base = 0;
    for (int i = 0; i < np; i++) {
        int m = parts[i];
        if (m >= 4) {
            for (int a = 2; a + a <= m; a++) {
                int cnt = (a + a == m) ? m / 2 : m; /* unordered pairs */
                for (int k = 0; k < cnt; k++) {
                    int j = base + k;
                    int jp = j;
                    for (int t = 0; t < a; t++) jp = sig[jp];
                    mj[nmark] = j; mjp[nmark] = jp;
                    mm[nmark] = m; ma[nmark] = a;
                    nmark++;
                }
            }
        }
        base += parts[i];
    }
    /* rows 0,1 fixed */
    for (int c = 0; c < n; c++) {
        L[0][c] = c; L[1][c] = sig[c];
        colmask[c] = (1u << c) | (1u << sig[c]);
    }
    rec(2, 0, 0);
    printf("n=%d type=(", n);
    for (int i = 0; i < np; i++) printf("%s%d", i ? "," : "", parts[i]);
    printf(")  completions=%lld\n", NC);
    for (int m = 4; m <= n; m++)
        for (int a = 2; a + a <= m; a++)
            if (XA[m][a] + XB[m][a])
                printf("  m=%d alpha=%d  XA=%lld XB=%lld  XB/XA=%.9f\n",
                       m, a, XA[m][a], XB[m][a],
                       (double)XB[m][a] / (double)XA[m][a]);
    return 0;
}
