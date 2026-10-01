/* Jacobson-Matthews uniform sampler for Latin squares, O(1) per move.
 * Emits sampled squares (proper states only) as raw bytes (n*n symbols per
 * sample) to stdout.
 *
 * Usage: ./jm n nsamples thin seed
 *
 * State: incidence tensor T in {0,1,-1} (one -1 cell iff improper), plus for
 * every line (row,sym), (col,sym), (row,col) the list of its +1 cells
 * (at most 2 entries; exactly 2 on the three lines through the -1 cell when
 * improper, and on no others).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int n;
static int8_t *T;
#define TT(r,c,s) T[((r)*n+(c))*n+(s)]

/* line lists: L_rs[r][s] -> cols, L_cs[c][s] -> rows, L_rc[r][c] -> syms */
static uint8_t *rs_e, *rs_n;   /* rs_e[(r*n+s)*2 + k], rs_n[r*n+s] */
static uint8_t *cs_e, *cs_n;
static uint8_t *rc_e, *rc_n;

static uint64_t rng_state;
static inline uint64_t rng(void){          /* splitmix64 */
    uint64_t z = (rng_state += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}
static inline int rnd(int m){ return (int)(rng() % (uint64_t)m); }

static inline void ladd(uint8_t *e, uint8_t *cnt, int idx, int v){
    e[idx*2 + cnt[idx]] = (uint8_t)v;
    cnt[idx]++;
}
static inline void lrem(uint8_t *e, uint8_t *cnt, int idx, int v){
    if (e[idx*2] == v){ e[idx*2] = e[idx*2+1]; cnt[idx]--; }
    else { cnt[idx]--; }   /* v was in slot 1 */
}

static int neg_r, neg_c, neg_s, improper;

static inline void inc_cell(int r, int c, int s){
    int8_t t = ++TT(r,c,s);
    if (t == 1){            /* 0 -> 1: add to lists */
        ladd(rs_e, rs_n, r*n+s, c);
        ladd(cs_e, cs_n, c*n+s, r);
        ladd(rc_e, rc_n, r*n+c, s);
    }
    /* -1 -> 0: nothing to add (negative cell was not listed) */
}
static inline void dec_cell(int r, int c, int s){
    int8_t t = --TT(r,c,s);
    if (t == 0){            /* 1 -> 0: remove from lists */
        lrem(rs_e, rs_n, r*n+s, c);
        lrem(cs_e, cs_n, c*n+s, r);
        lrem(rc_e, rc_n, r*n+c, s);
    } else {                /* 0 -> -1: new negative cell */
        neg_r = r; neg_c = c; neg_s = s;
    }
}

static void move(void){
    int r, c, s, r0, c0, s0;
    if (!improper){
        r = rnd(n); c = rnd(n);
        s0 = rc_e[(r*n+c)*2];
        do { s = rnd(n); } while (s == s0);
        c0 = rs_e[(r*n+s)*2];
        r0 = cs_e[(c*n+s)*2];
    } else {
        r = neg_r; c = neg_c; s = neg_s;
        r0 = cs_e[(c*n+s)*2 + (cs_n[c*n+s] > 1 ? rnd(2) : 0)];
        c0 = rs_e[(r*n+s)*2 + (rs_n[r*n+s] > 1 ? rnd(2) : 0)];
        s0 = rc_e[(r*n+c)*2 + (rc_n[r*n+c] > 1 ? rnd(2) : 0)];
    }
    improper = 0;
    inc_cell(r,c,s);
    dec_cell(r,c0,s);
    dec_cell(r0,c,s);
    dec_cell(r,c,s0);
    inc_cell(r,c0,s0);
    inc_cell(r0,c,s0);
    inc_cell(r0,c0,s);
    dec_cell(r0,c0,s0);
    if (TT(r0,c0,s0) == -1) improper = 1;
}

#ifndef MAIN_EXCLUDED
int main(int argc, char **argv){
    n = atoi(argv[1]);
    long nsamples = atol(argv[2]);
    long thin = atol(argv[3]);
    rng_state = strtoull(argv[4], NULL, 10) * 2654435761u + 88172645463325252ull;
    T = calloc((size_t)n*n*n, 1);
    rs_e = calloc((size_t)n*n*2,1); rs_n = calloc((size_t)n*n,1);
    cs_e = calloc((size_t)n*n*2,1); cs_n = calloc((size_t)n*n,1);
    rc_e = calloc((size_t)n*n*2,1); rc_n = calloc((size_t)n*n,1);
    for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++){
            int s = (r + c) % n;
            TT(r,c,s) = 1;
            ladd(rs_e, rs_n, r*n+s, c);
            ladd(cs_e, cs_n, c*n+s, r);
            ladd(rc_e, rc_n, r*n+c, s);
        }
    improper = 0;
    for (long i = 0; i < 200*thin; i++) move();
    uint8_t *buf = malloc((size_t)n*n);
    for (long smp = 0; smp < nsamples; smp++){
        /* fixed-grid sampling: advance in windows of exactly `thin` moves and
           discard improper checkpoints (unbiased; emission time never depends
           on the state). */
        do {
            for (long i = 0; i < thin; i++) move();
        } while (improper);
        for (int r = 0; r < n; r++)
            for (int c = 0; c < n; c++)
                buf[r*n+c] = rc_e[(r*n+c)*2];
        fwrite(buf, 1, (size_t)n*n, stdout);
    }
    return 0;
}
#endif
