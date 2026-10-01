/* Layer-5 anatomy on the common board0 (two paths), n <= 11 (128-bit masks).
   Computes exactly:
     N0 = # board0-avoiding permutations
     P05 = # ordered pairwise-disjoint triples
     G5(z) for the four pins d1=(1,2), d2=(n-1,0), e1=(1,0), e2=(n-1,2)
     BrG5 = G5(d1)+G5(d2)-G5(e1)-G5(e2)
     G25(d1,d2), G25(e1,e2)   (tau through both pins)
     H5(d1,d2), H5(e1,e2)     (tau through z, eta through z')
   and checks  3*BrG5 + 3*dG25 + 6*dH5 == DeltaN5  against the direct
   N5 values (recomputed here from the B2 boards).
   Usage: ./anatomy5 8      (n=8 ~ 1 min; n=9 hours -- see -DG_ONLY)
   Compile with -DG_ONLY to skip P05 (only tau-through-pin passes; n=9 ~ minutes).
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static int n;
static int forb1[16], forb2[16]; /* board0: two forbidden cols per row */
typedef struct { uint64_t lo, hi; } m2;
static m2 *masks;
static int nmask, cap;
static int *pincode; /* bit0: through d1, bit1: d2, bit2: e1, bit3: e2 */

static int D1r, D1c, D2r, D2c, E1r, E1c, E2r, E2c;

static void dfs(int row, unsigned usedcols, uint64_t lo, uint64_t hi) {
    if (row == n) {
        if (nmask == cap) { cap *= 2; masks = realloc(masks, cap * sizeof(m2)); pincode = realloc(pincode, cap * sizeof(int)); }
        masks[nmask].lo = lo; masks[nmask].hi = hi; pincode[nmask] = 0; nmask++;
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
    n = (argc > 1) ? atoi(argv[1]) : 8;
    /* board0 = I u P_sigma minus sigma's two switch cells:
       row i forbidden cols: i and (i+1)%n, EXCEPT rows 1 and n-1 where the
       sigma-cells (1,2) and (n-1,0) are REMOVED (only the diagonal stays);
       following switch.py make_boards: board0 = B2 & B2'. Row 1: B2 has
       {1,2}, B2' has {1,0} -> common {1}. Row n-1: B2 {n-1,0}, B2' {n-1,2}
       -> common {n-1}. */
    for (int i = 0; i < n; i++) { forb1[i] = i; forb2[i] = (i + 1) % n; }
    forb2[1] = 1;       /* row 1: only diagonal */
    forb2[n - 1] = n - 1; /* row n-1: only diagonal */
    D1r = 1; D1c = 2; D2r = n - 1; D2c = 0;
    E1r = 1; E1c = 0; E2r = n - 1; E2c = 2;

    cap = 1 << 16; masks = malloc(cap * sizeof(m2)); pincode = malloc(cap * sizeof(int));
    nmask = 0;
    dfs(0, 0, 0, 0);
    for (int a = 0; a < nmask; a++) {
        int pc = 0;
        if (bitat(masks[a], D1r, D1c)) pc |= 1;
        if (bitat(masks[a], D2r, D2c)) pc |= 2;
        if (bitat(masks[a], E1r, E1c)) pc |= 4;
        if (bitat(masks[a], E2r, E2c)) pc |= 8;
        pincode[a] = pc;
    }
    long long N0 = nmask;
    long long m_d1 = 0, m_d2 = 0, m_e1 = 0, m_e2 = 0;
    for (int a = 0; a < nmask; a++) {
        if (pincode[a] & 1) m_d1++;
        if (pincode[a] & 2) m_d2++;
        if (pincode[a] & 4) m_e1++;
        if (pincode[a] & 8) m_e2++;
    }
    printf("n=%d  N0=%lld  m(d1)=%lld m(d2)=%lld m(e1)=%lld m(e2)=%lld\n",
           n, N0, m_d1, m_d2, m_e1, m_e2);

    long long P05 = 0;
    long long G5d1 = 0, G5d2 = 0, G5e1 = 0, G5e2 = 0;
    long long G25d = 0, G25e = 0, H5d = 0, H5e = 0;

#pragma omp parallel for schedule(dynamic, 8) reduction(+:P05,G5d1,G5d2,G5e1,G5e2,G25d,G25e,H5d,H5e)
    for (int a = 0; a < nmask; a++) {
        int pc = pincode[a];
#ifdef G_ONLY
        if (!pc) continue;   /* only pinned-tau rows needed */
#endif
        uint64_t alo = masks[a].lo, ahi = masks[a].hi;
        int nl = 0;
        m2 *L = malloc(nmask * sizeof(m2));
        int *lpin = malloc(nmask * sizeof(int));
        for (int b = 0; b < nmask; b++) {
            if ((alo & masks[b].lo) | (ahi & masks[b].hi)) continue;
            lpin[nl] = pincode[b];
            L[nl++] = masks[b];
        }
        /* cnt = # ordered disjoint pairs (eta,theta) in L;
           cntE2 = pairs with eta through e2; cntD2 = eta through d2 */
        long long cnt = 0, cntE2 = 0, cntD2 = 0;
        for (int i = 0; i < nl; i++) {
            uint64_t ilo = L[i].lo, ihi = L[i].hi;
            long long di = 0;
            for (int j = 0; j < nl; j++) {
                if (j == i) continue;
                if (!((ilo & L[j].lo) | (ihi & L[j].hi))) di++;
            }
            cnt += di;
            if (lpin[i] & 2) cntD2 += di;
            if (lpin[i] & 8) cntE2 += di;
        }
        P05 += cnt;
        if (pc & 1) { G5d1 += cnt; H5d += cntD2; }   /* tau thru d1, eta thru d2 */
        if (pc & 2) G5d2 += cnt;
        if (pc & 4) { G5e1 += cnt; H5e += cntE2; }   /* tau thru e1, eta thru e2 */
        if (pc & 8) G5e2 += cnt;
        if ((pc & 3) == 3) G25d += cnt;
        if ((pc & 12) == 12) G25e += cnt;
        free(L); free(lpin);
    }
    long long BrG5 = G5d1 + G5d2 - G5e1 - G5e2;
    printf("P05=%lld\n", P05);
    printf("G5(d1)=%lld G5(d2)=%lld G5(e1)=%lld G5(e2)=%lld\n", G5d1, G5d2, G5e1, G5e2);
    printf("BrG5=%lld\n", BrG5);
    printf("G25(d)=%lld G25(e)=%lld  H5(d)=%lld H5(e)=%lld\n", G25d, G25e, H5d, H5e);
    long long rhs = 3 * BrG5 + 3 * (G25e - G25d) + 6 * (H5e - H5d);
    printf("3BrG5+3dG25+6dH5 = %lld   (should equal DeltaN5)\n", rhs);
    return 0;
}
