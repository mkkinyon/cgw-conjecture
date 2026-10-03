/* Jacobson-Matthews sampler restricted to completions of a FIXED 2 x n rectangle (rows 0 and 1 never change).
 * The walk chooses uniformly among the JM moves that touch neither row 0 nor row 1; the move graph is undirected
 * and every proper state has the same restricted degree (n-2) n (n-3), so the stationary law restricted to proper
 * states is uniform on the connected component of the start.  Connectivity of the restricted graph is NOT proved;
 * test against exact enumeration (n = 6, 7) before trusting it.
 * Usage: ./jmfix n nsamples thin seed < initial_square.bin    (n*n bytes, row-major; rows 0,1 are the rectangle)
 * Output: nsamples proper squares, n*n bytes each. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static int n; static int8_t *T;
#define TT(r,c,s) T[((r)*n+(c))*n+(s)]
static uint8_t *rs_e,*rs_n,*cs_e,*cs_n,*rc_e,*rc_n;
static uint64_t rng_state;
static inline uint64_t rng(void){ uint64_t z=(rng_state+=0x9E3779B97F4A7C15ull); z=(z^(z>>30))*0xBF58476D1CE4E5B9ull; z=(z^(z>>27))*0x94D049BB133111EBull; return z^(z>>31);}
static inline int rnd(int m){ return (int)(rng()%(uint64_t)m); }
static inline void ladd(uint8_t*e,uint8_t*cnt,int idx,int v){ e[idx*2+cnt[idx]]=(uint8_t)v; cnt[idx]++; }
static inline void lrem(uint8_t*e,uint8_t*cnt,int idx,int v){ if(e[idx*2]==v){ e[idx*2]=e[idx*2+1]; cnt[idx]--; } else cnt[idx]--; }
static int neg_r,neg_c,neg_s,improper;
static inline void inc_cell(int r,int c,int s){ int8_t t=++TT(r,c,s); if(t==1){ ladd(rs_e,rs_n,r*n+s,c); ladd(cs_e,cs_n,c*n+s,r); ladd(rc_e,rc_n,r*n+c,s);} }
static inline void dec_cell(int r,int c,int s){ int8_t t=--TT(r,c,s); if(t==0){ lrem(rs_e,rs_n,r*n+s,c); lrem(cs_e,cs_n,c*n+s,r); lrem(rc_e,rc_n,r*n+c,s);} else { neg_r=r; neg_c=c; neg_s=s; } }
static long rejected=0;
static void move(void){
    int r,c,s,r0,c0,s0;
    if(!improper){
        r=2+rnd(n-2); c=rnd(n);
        s0=rc_e[(r*n+c)*2];
        int a=rc_e[(0*n+c)*2], b=rc_e[(1*n+c)*2];
        do { s=rnd(n); } while(s==s0||s==a||s==b);         /* s not in rows 0,1 of column c => r0 >= 2 */
        c0=rs_e[(r*n+s)*2]; r0=cs_e[(c*n+s)*2];
    } else {
        r=neg_r; c=neg_c; s=neg_s;                          /* r >= 2 by induction */
        int k=cs_n[c*n+s]; int cand[2]; int m=0;
        for(int i=0;i<k;i++){ int rr=cs_e[(c*n+s)*2+i]; if(rr>=2) cand[m++]=rr; }
        /* uniform over the ORIGINAL choice set, reject if forbidden: keeps the walk a simple random walk on the restricted graph */
        int pick=cs_e[(c*n+s)*2+(k>1?rnd(2):0)];
        if(pick<2){ rejected++; return; }
        r0=pick;
        c0=rs_e[(r*n+s)*2+(rs_n[r*n+s]>1?rnd(2):0)];
        s0=rc_e[(r*n+c)*2+(rc_n[r*n+c]>1?rnd(2):0)];
    }
    improper=0;
    inc_cell(r,c,s); dec_cell(r,c0,s); dec_cell(r0,c,s); dec_cell(r,c,s0);
    inc_cell(r,c0,s0); inc_cell(r0,c,s0); inc_cell(r0,c0,s); dec_cell(r0,c0,s0);
    if(TT(r0,c0,s0)==-1) improper=1;
}
int main(int argc,char**argv){
    n=atoi(argv[1]); long nsamples=atol(argv[2]); long thin=atol(argv[3]);
    rng_state=strtoull(argv[4],NULL,10)*2654435761u+88172645463325252ull;
    T=calloc((size_t)n*n*n,1);
    rs_e=calloc((size_t)n*n*2,1); rs_n=calloc((size_t)n*n,1); cs_e=calloc((size_t)n*n*2,1); cs_n=calloc((size_t)n*n,1); rc_e=calloc((size_t)n*n*2,1); rc_n=calloc((size_t)n*n,1);
    uint8_t *buf=malloc((size_t)n*n);
    if(fread(buf,1,(size_t)n*n,stdin)!=(size_t)n*n){ fprintf(stderr,"need n*n bytes of initial square\n"); return 1; }
    for(int r=0;r<n;r++) for(int c=0;c<n;c++){ int s=buf[r*n+c]; TT(r,c,s)=1; ladd(rs_e,rs_n,r*n+s,c); ladd(cs_e,cs_n,c*n+s,r); ladd(rc_e,rc_n,r*n+c,s); }
    improper=0;
    for(long i=0;i<200*thin;i++) move();
    for(long smp=0;smp<nsamples;smp++){
        do { for(long i=0;i<thin;i++) move(); } while(improper);
        for(int r=0;r<n;r++) for(int c=0;c<n;c++) buf[r*n+c]=rc_e[(r*n+c)*2];
        fwrite(buf,1,(size_t)n*n,stdout);
    }
    fprintf(stderr,"jmfix: rejected improper proposals = %ld\n",rejected);
    return 0;
}
