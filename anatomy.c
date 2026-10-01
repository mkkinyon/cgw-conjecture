/* Exact anatomy of BrG at order n:
   board0 = I u P_sigma minus cells (1,2),(n-1,0)  [0-indexed],
   sigma = cycle i -> i+1 mod n.
   Enumerate all board0-avoiding permutations tau; for each compute
   A(tau) = per(J - board0 - P_tau) by Ryser/Gray; accumulate
     N0, sumA, and G(d) for d in {(1,2),(n-1,0),(1,0),(n-1,2)}.
   Usage: ./anatomy n
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef unsigned long long u64;
typedef __int128 i128;

static int n;
static u64 allowed[32];   /* allowed[i] = bitmask of admissible cols for row i (complement of board0 row) */
static u64 rowmask_perm[32]; /* current tau as col per row */
static i128 N0 = 0, sumA = 0;
static i128 Gd1 = 0, Gd2 = 0, Ge1 = 0, Ge2 = 0;
static i128 m_d1 = 0, m_d2 = 0, m_e1 = 0, m_e2 = 0;

/* Ryser with Gray code: per of matrix with row masks rm[i] (bit c = entry 1) */
static i128 per_ryser(u64 *rm) {
    int N = n;
    i128 total = 0;
    long long rowsum[32];
    for (int i = 0; i < N; i++) rowsum[i] = 0;
    u64 gray = 0;
    long long sign = -1; /* will flip; standard: sum over nonempty subsets */
    /* iterate Gray code over 1..2^N-1 */
    u64 lim = 1ULL << N;
    for (u64 k = 1; k < lim; k++) {
        u64 g = k ^ (k >> 1);
        u64 diff = g ^ gray;   /* one bit */
        int b = __builtin_ctzll(diff);
        if (g & diff) { for (int i = 0; i < N; i++) rowsum[i] += (rm[i] >> b) & 1ULL; }
        else          { for (int i = 0; i < N; i++) rowsum[i] -= (rm[i] >> b) & 1ULL; }
        gray = g;
        i128 prod = 1;
        for (int i = 0; i < N; i++) { prod *= rowsum[i]; if (!prod) break; }
        int bits = __builtin_popcountll(g);
        if ((N - bits) & 1) total -= prod; else total += prod;
    }
    return total;
}

static int taucol[32];

static void rec(int row, u64 used) {
    if (row == n) {
        /* compute A(tau) */
        u64 rm[32];
        for (int i = 0; i < n; i++) rm[i] = allowed[i] & ~(1ULL << taucol[i]);
        i128 a = per_ryser(rm);
        N0++; sumA += a;
        if (taucol[1] == 2) { Gd1 += a; m_d1++; }
        if (taucol[n-1] == 0) { Gd2 += a; m_d2++; }
        if (taucol[1] == 0) { Ge1 += a; m_e1++; }
        if (taucol[n-1] == 2) { Ge2 += a; m_e2++; }
        return;
    }
    u64 av = allowed[row] & ~used;
    while (av) {
        int c = __builtin_ctzll(av);
        av &= av - 1;
        taucol[row] = c;
        rec(row + 1, used | (1ULL << c));
    }
}

static void p128(const char *lbl, i128 v) {
    char buf[64]; int p = 63; buf[p] = 0;
    int neg = v < 0; if (neg) v = -v;
    if (!v) buf[--p] = '0';
    while (v) { buf[--p] = '0' + (int)(v % 10); v /= 10; }
    if (neg) buf[--p] = '-';
    printf("%s%s", lbl, buf + p);
}

int main(int argc, char **argv) {
    n = atoi(argv[1]);
    for (int i = 0; i < n; i++) {
        u64 forb = (1ULL << i);
        int s = (i + 1) % n;
        if (!(i == 1 || i == n - 1)) forb |= (1ULL << s); /* sigma-cells except the two removed */
        allowed[i] = ((1ULL << n) - 1) & ~forb;
    }
    rec(0, 0);
    printf("n=%d\n", n);
    p128(" N0=", N0); p128("  sumA=", sumA); printf("\n");
    p128(" m+=", m_d1 + m_d2); p128(" m-=", m_e1 + m_e2); printf("\n");
    p128(" G(d1)=", Gd1); p128(" G(d2)=", Gd2); p128(" G(e1)=", Ge1); p128(" G(e2)=", Ge2); printf("\n");
    p128(" BrG=", Gd1 + Gd2 - Ge1 - Ge2); printf("\n");
    double EA = (double)sumA / (double)N0;
    double EpA = (double)(Gd1 + Gd2) / (double)(m_d1 + m_d2);
    double EmA = (double)(Ge1 + Ge2) / (double)(m_e1 + m_e2);
    printf(" E[A]=%.6f E+[A]=%.6f E-[A]=%.6f\n", EA, EpA, EmA);
    printf(" (E+-E-)/E=%.6e  (E--E)/E=%.6e\n", (EpA - EmA) / EA, (EmA - EA) / EA);
    return 0;
}
