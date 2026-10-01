/* Exact N4(sigma) = number of 4 x n Latin rectangles with first two rows
 * (id, sigma):  N4 = sum over third rows tau discordant with {id,sigma} of
 * per(complement of the 3-row rectangle).
 *
 * Usage: ./layer4 n type   where type 0 = single n-cycle, 1 = (2, n-2),
 *                                2 = (3, n-3)
 * Prints N3 (number of taus) and N4.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef unsigned __int128 u128;

static int n;
static int sigma[16];
static int tau[16];
static uint16_t avail[16];    /* symbols available for tau at column i (bitmask) */
static uint16_t colmask3[16]; /* forbidden symbols per column for 4th row: id, sigma, tau */
static long long N3 = 0;
static u128 N4 = 0;

/* Ryser with Gray code: per of n x n 0-1 matrix given as row masks (allowed columns).
 * Here we need per of matrix M[i][j] = allowed symbol j in column i.
 * We compute over subsets of columns(symbols).  */
static int64_t ryser(const uint16_t *rowmask)
{
    /* per(A) = (-1)^n * sum_{S subset} (-1)^{|S|} prod_i (sum_{j in S} A[i][j]) */
    int64_t total = 0;
    int cnt[16];
    for (int i = 0; i < n; i++) cnt[i] = 0;
    int64_t sign = (n % 2 == 0) ? 1 : -1;
    uint32_t nsub = 1u << n;
    int popc = 0;
    for (uint32_t g = 1; g < nsub; g++) {
        uint32_t gray = g ^ (g >> 1);
        uint32_t prev = (g - 1) ^ ((g - 1) >> 1);
        uint32_t diff = gray ^ prev;
        int j = __builtin_ctz(diff);
        int add = (gray >> j) & 1;
        popc += add ? 1 : -1;
        for (int i = 0; i < n; i++)
            cnt[i] += ((rowmask[i] >> j) & 1) ? (add ? 1 : -1) : 0;
        int64_t prod = 1;
        for (int i = 0; i < n; i++) {
            prod *= cnt[i];
            if (!prod) break;
        }
        total += ((popc & 1) ? -1 : 1) * prod;
    }
    return sign * total;
}

static void rec(int col, uint16_t used)
{
    if (col == n) {
        N3++;
        /* build allowed mask for 4th row: complement of {i, sigma[i], tau[i]} */
        uint16_t rm[16];
        uint16_t full = (uint16_t)((1u << n) - 1);
        for (int i = 0; i < n; i++)
            rm[i] = full & ~((uint16_t)((1u << i) | (1u << sigma[i]) | (1u << tau[i])));
        N4 += (u128)ryser(rm);
        return;
    }
    uint16_t a = avail[col] & ~used;
    while (a) {
        int s = __builtin_ctz(a);
        a &= a - 1;
        tau[col] = s;
        rec(col + 1, used | (uint16_t)(1u << s));
    }
}

int main(int argc, char **argv)
{
    n = atoi(argv[1]);
    int type = atoi(argv[2]);
    /* build sigma */
    if (type == 0) {
        for (int i = 0; i < n; i++) sigma[i] = (i + 1) % n;
    } else if (type == 1) {
        sigma[0] = 1; sigma[1] = 0;
        for (int i = 2; i < n; i++) sigma[i] = (i - 2 + 1) % (n - 2) + 2;
    } else {
        sigma[0] = 1; sigma[1] = 2; sigma[2] = 0;
        for (int i = 3; i < n; i++) sigma[i] = (i - 3 + 1) % (n - 3) + 3;
    }
    uint16_t full = (uint16_t)((1u << n) - 1);
    for (int i = 0; i < n; i++)
        avail[i] = full & ~((uint16_t)((1u << i) | (1u << sigma[i])));
    rec(0, 0);
    unsigned long long hi = (unsigned long long)(N4 >> 64);
    unsigned long long lo = (unsigned long long)N4;
    printf("n=%d type=%d N3=%lld N4_hi=%llu N4_lo=%llu\n", n, type, N3, hi, lo);
    return 0;
}
