#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static unsigned long long rng_s[2];
static inline unsigned long long rng_next(void) {
    unsigned long long s1 = rng_s[0], s0 = rng_s[1];
    rng_s[0] = s0;
    s1 ^= s1 << 23;
    rng_s[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
    return rng_s[1] + s0;
}
static inline double rng_u(void) {
    return (rng_next() >> 11) * (1.0 / 9007199254740992.0);
}

/* ----------------- marked-partition state ----------------- */
/* blocks 0..nb-1 with sizes sz[]; ub, vb = indices of u's / v's
 * block; if ub==vb, (sx,sy) = sides (u->v, v->u) along that curve.
 * Free-list-less: blocks are swapped-from-end on deletion. */
typedef struct {
    int nb, ub, vb, cap;
    double *sz;
    double sx, sy;      /* valid iff ub==vb */
} St;

static void st_init(St *S, int cap) {
    S->cap = cap;
    S->sz = malloc(sizeof(double) * cap);
}
static void st_copy(St *D, const St *S) {
    D->nb = S->nb; D->ub = S->ub; D->vb = S->vb;
    D->sx = S->sx; D->sy = S->sy;
    memcpy(D->sz, S->sz, sizeof(double) * S->nb);
}

/* sample a block index with prob proportional to size (linear scan;
 * nb stays modest, and this keeps the code transparent) */
static inline int pick_block(const St *S) {
    double r = rng_u();
    double acc = 0;
    for (int i = 0; i < S->nb; i++) {
        acc += S->sz[i];
        if (r < acc) return i;
    }
    return S->nb - 1;
}

static inline void del_block(St *S, int i) {
    int last = S->nb - 1;
    S->sz[i] = S->sz[last];
    if (S->ub == last) S->ub = i;
    if (S->vb == last) S->vb = i;
    S->nb--;
}

/* one split-merge step (one fresh chord) */
static long NFLIP=0;
static void step(St *S) {
    int i1 = pick_block(S), i2 = pick_block(S);
    if (i1 == i2) {                      /* split */
        double s = S->sz[i1];
        if (i1 == S->ub && i1 == S->vb) {       /* together block */
            double x = S->sx, y = S->sy;
            double c1 = rng_u() * s, c2 = rng_u() * s;
            /* positions along cycle: u at 0, v at x */
            int c1x = c1 < x, c2x = c2 < x;
            if (c1x && c2x) {                   /* both cuts in x-side */
                double piece = fabs(c1 - c2);
                S->sz[i1] = s - piece;  S->sx = x - piece;
                S->sz[S->nb++] = piece;
            } else if (!c1x && !c2x) {          /* both in y-side */
                double piece = fabs(c1 - c2);
                S->sz[i1] = s - piece;  S->sy = y - piece;
                S->sz[S->nb++] = piece;
            } else {                            /* FLIP: separate u,v */
                NFLIP++;
                double t1 = c1x ? c1 : c2;      /* in (0,x)  */
                double t2 = c1x ? c2 : c1;      /* in (x,s)  */
                double q = t2 - t1;             /* v's piece */
                S->sz[i1] = s - q;              /* u keeps i1 */
                S->sz[S->nb] = q;
                S->vb = S->nb++;
            }
        } else {                                /* <=1 mark in block */
            double t1 = rng_u() * s, t2 = rng_u() * s;
            double piece = fabs(t1 - t2);
            /* if a mark is in the block it sits at position 0, i.e.
             * in the complement piece s-piece (cuts uniform rel. mark);
             * for unmarked blocks the labelling is irrelevant */
            S->sz[i1] = s - piece;
            S->sz[S->nb++] = piece;
        }
    } else {                              /* merge i1,i2 */
        int hasu1 = (i1 == S->ub), hasv1 = (i1 == S->vb);
        int hasu2 = (i2 == S->ub), hasv2 = (i2 == S->vb);
        double p = S->sz[i1], q = S->sz[i2];
        if ((hasu1 && hasv2) || (hasv1 && hasu2)) {   /* FLIP: join u,v */
            NFLIP++;
            double pu = hasu1 ? p : q;   /* u's old block size */
            double pv = hasu1 ? q : p;
            double a = rng_u() * pu, b = rng_u() * pv;
            S->sz[i1] = p + q;
            S->ub = S->vb = i1;
            S->sx = a + b;  S->sy = (p + q) - (a + b);
            del_block(S, i2);
        } else if ((hasu1 && hasv1)) {      /* together-block + env */
            double s = p;
            double c = rng_u() * s;
            if (c < S->sx) S->sx += q; else S->sy += q;
            S->sz[i1] = p + q;
            del_block(S, i2);
        } else if ((hasu2 && hasv2)) {
            double s = q;
            double c = rng_u() * s;
            if (c < S->sx) S->sx += p; else S->sy += p;
            S->sz[i2] = p + q;
            del_block(S, i1);
        } else {                            /* at most one mark each,
                                               not u-and-v across */
            S->sz[i1] = p + q;
            if (S->ub == i2) S->ub = i1;
            if (S->vb == i2) S->vb = i1;
            del_block(S, i2);
        }
    }
}

static inline int eps_st(const St *S) { return S->ub == S->vb ? 1 : -1; }

/* ------------- exact initial state from j iid chords ------------- */
/* discrete surgery on 2j+2 points; u=0.0 (pt 0), v=0.5 (pt 1);
 * computes the curve partition with arc measures, fills St */
#define MAXP 8300
static double coords[MAXP];
static int order[MAXP], rankp[MAXP], nxt[MAXP], part[MAXP], cyc[MAXP];

static int cmp_pts(const void *x, const void *y) {
    double a = coords[*(const int *)x], b = coords[*(const int *)y];
    return (a > b) - (a < b);
}

static void init_from_chords(St *S, int j) {
    int n = 2 * j + 2;
    coords[0] = 0.0;  part[0] = -1;
    coords[1] = 0.5;  part[1] = -1;
    for (int c = 0; c < j; c++) {
        coords[2 + 2 * c] = rng_u();
        coords[3 + 2 * c] = rng_u();
        part[2 + 2 * c] = 3 + 2 * c;
        part[3 + 2 * c] = 2 + 2 * c;
    }
    for (int i = 0; i < n; i++) order[i] = i;
    qsort(order, n, sizeof(int), cmp_pts);
    for (int i = 0; i < n; i++) rankp[order[i]] = i;
    for (int i = 0; i < n; i++) {
        int p = order[i];
        int q = (part[p] >= 0) ? part[p] : p;
        nxt[p] = order[(rankp[q] + 1) % n];   /* rho = gamma o tau */
    }
    /* cycles of rho; measure of a cycle = sum over its points p of
     * the arc [p, circle-successor(p)) */
    for (int i = 0; i < n; i++) cyc[i] = -1;
    S->nb = 0;
    for (int i = 0; i < n; i++) {
        if (cyc[i] >= 0) continue;
        int id = S->nb++;
        double m = 0;
        int p = i;
        do {
            cyc[p] = id;
            double a0 = coords[p];
            double a1 = coords[order[(rankp[p] + 1) % n]];
            m += (p == order[n - 1]) ? (1.0 - a0 + coords[order[0]])
                                     : (a1 - a0);
            p = nxt[p];
        } while (p != i);
        S->sz[id] = m;
    }
    S->ub = cyc[0];
    S->vb = cyc[1];
    if (S->ub == S->vb) {
        /* side x = measure from u to v along the curve */
        double x = 0;
        int p = 0;
        while (p != 1) {
            double a0 = coords[p];
            double a1 = coords[order[(rankp[p] + 1) % n]];
            x += (p == order[n - 1]) ? (1.0 - a0 + coords[order[0]])
                                     : (a1 - a0);
            p = nxt[p];
        }
        S->sx = x;
        S->sy = S->sz[S->ub] - x;
    }
}


/* s11_kflip.c -- session 11: diagnostics for condition (R2) along the
 * post-flip chain (TODO 11d, with the corrected prediction).
 *
 * K_z F(x) = E_x[z^tau F(Gamma_tau)].  (R2) asks that -1 is not in the
 * spectrum of K_z, uniformly for |z|<=1.  Session 11 proves
 * K_{-z} = M K_z M with M = multiplication by (-1)^{B(Gamma)} (B = total
 * block count), because B changes by exactly +-1 at every step, so
 * (-1)^tau = (-1)^{B(Gamma_tau)-B(x)}.  Hence spec(K_{-1}) = spec(K_1)
 * and the right diagnostic at z=-1 is NOT "a_k = E_x[(-1)^{tau_k}]
 * decays geometrically" (it need not: K_{-1} has eigenvalue +1 with
 * eigenfunction (-1)^B) but "a_k converges without a (-1)^k
 * oscillation".  An oscillation a_k ~ c(-1)^k with |c| bounded away
 * from 0 would mean -1 in spec(K_1): the NEGATIVE finding.
 *
 * mode: kflip lo hi kmax cs maxsteps [seed]
 *   stationary-ish start (burn 4096 steps from one block), accepted if
 *   gap in [lo,hi).  Then run to kmax flips (censor at maxsteps steps).
 *   prints, for k=1..kmax:
 *     a_k   = E[(-1)^{tau_k}]                       = (K_{-1}^k 1)(x)
 *     s_k   = E[(-1)^{B(Gamma_{tau_k})}]            (= (-1)^{B_0} a_k, exact identity)
 *     r_k   = E[(-1)^{tau_k - tau_{k-1}}]           (inter-flip parity)
 *     c_k   = E[sig_0 sig_k], sig_k = (-1)^{B(Gamma_{tau_k})}
 *     mean and median of tau_k - tau_{k-1}, and E[gap] at tau_k.
 *   and the tail P(tau_k - tau_{k-1} > t) at powers of two for k=1,2,4,8.
 */

static void stat_start(St *C, int burn){
    C->nb=1; C->sz[0]=1.0; C->ub=C->vb=0; C->sx=0.5; C->sy=0.5;
    for(int s=0;s<burn;s++) step(C);
}
static double gap(const St *C){
    return C->ub==C->vb ? (C->sx<C->sy?C->sx:C->sy)
        : (C->sz[C->ub]<C->sz[C->vb]?C->sz[C->ub]:C->sz[C->vb]);
}
static double stot(const St *C){
    return C->ub==C->vb ? C->sz[C->ub] : C->sz[C->ub]+C->sz[C->vb];
}
static int parity_ok = 1;
int main2(int argc, char **argv);
int main(int argc, char **argv) { return main2(argc, argv); }
int main2(int argc,char**argv){
    if(argc<6){ fprintf(stderr,"usage: %s kflip lo hi kmax cs maxsteps [seed]\n",argv[0]); return 1; }
    double lo=atof(argv[2]), hi=atof(argv[3]);
    int kmax=atoi(argv[4]); long cs=atol(argv[5]); long maxsteps=atol(argv[6]);
    unsigned seed = (argc>7)?(unsigned)atoi(argv[7]):12345;
    rng_s[0]=seed*2654435761ULL+1; rng_s[1]=seed^0x9e3779b97f4a7c15ULL;
    for(int i=0;i<64;i++) rng_next();

    St C; st_init(&C, 2*(4096+maxsteps)+64);
    double *a=calloc(kmax+1,sizeof(double)), *sg=calloc(kmax+1,sizeof(double));
    double *r=calloc(kmax+1,sizeof(double)), *cc=calloc(kmax+1,sizeof(double));
    double *mt=calloc(kmax+1,sizeof(double)), *gp=calloc(kmax+1,sizeof(double));
    double *st_=calloc(kmax+1,sizeof(double));
    double *hg=calloc(kmax+1,sizeof(double)), *hs=calloc(kmax+1,sizeof(double));
    double *hgm=calloc(kmax+1,sizeof(double)), *hsm=calloc(kmax+1,sizeof(double));
    long *nk=calloc(kmax+1,sizeof(long));
    /* tails of tau_k - tau_{k-1} for k in {1,2,4,8}: 32 dyadic bins */
    long tail[4][40]; memset(tail,0,sizeof(tail));
    long tailn[4]={0,0,0,0};
    int kidx[4]={1,2,4,8};
    long done=0, censored=0, tried=0;
    long *tauv=calloc(kmax+2,sizeof(long));
    for(long t=0;t<cs;){
        stat_start(&C,4096); tried++;
        double g=gap(&C); if(g<lo||g>=hi) continue;
        t++;
        int B0=C.nb; int sig0 = (B0&1)?-1:1;
        int hg0 = (gap(&C)>0.25)?1:-1; int hs0 = (stot(&C)>0.80)?1:-1;
        long steps=0; int k=0; long prev=0;
        NFLIP=0;
        while(k<kmax && steps<maxsteps){
            int before=NFLIP; step(&C); steps++;
            if(NFLIP!=before){
                k++; tauv[k]=steps;
                long dt=steps-prev; prev=steps;
                int Bk=C.nb;
                /* parity identity check: B_k - B_0 == steps mod 2 */
                if(((Bk-B0)&1)!=(int)(steps&1)) parity_ok=0;
                int sgk=(Bk&1)?-1:1;
                int ptau=(steps&1)?-1:1;
                a[k]+=ptau; sg[k]+=sgk; r[k]+= (dt&1)?-1:1; cc[k]+=sig0*sgk;
                mt[k]+=(double)dt; gp[k]+=gap(&C); st_[k]+=stot(&C); nk[k]++;
                { int hgk=(gap(&C)>0.25)?1:-1, hsk=(stot(&C)>0.80)?1:-1;
                  hg[k]+=hg0*hgk; hs[k]+=hs0*hsk; hgm[k]+=hgk; hsm[k]+=hsk; }
                for(int q=0;q<4;q++) if(kidx[q]==k){
                    tailn[q]++; int b=0; long x=dt; while(x>1){x>>=1;b++;}
                    if(b>39)b=39; tail[q][b]++;
                }
            }
        }
        if(k<kmax) censored++;
        done++;
    }
    printf("# kflip gap in [%g,%g) cs=%ld maxsteps=%ld accept=%.4f censored=%ld parity_identity_B-B0==t_mod2: %s\n",
        lo,hi,cs,maxsteps,(double)cs/tried,censored,parity_ok?"HOLDS":"FAILS");
    printf("#  k     n_k      a_k=E(-1)^tau_k   s_k=E(-1)^{B_k}   r_k=E(-1)^{dtau}   c_k=E[sig0 sig_k]   E[dtau]   E[gap_k]  E[s_k]   Cgap_k   Cs_k\n");
    for(int k=1;k<=kmax;k++){
        long n=nk[k]; if(!n) continue;
        printf("  %3d %8ld   %+.5f          %+.5f          %+.5f          %+.5f       %9.2f  %.4f  %.4f  %+.5f  %+.5f\n",
            k,n,a[k]/n,sg[k]/n,r[k]/n,cc[k]/n,mt[k]/n,gp[k]/n,st_[k]/n,
            hg[k]/n-(hgm[k]/n)*0.0, hs[k]/n-(hsm[k]/n)*0.0);
    }
    printf("# se ~ %.5f\n",1.0/sqrt((double)cs));
    for(int q=0;q<4;q++){
        if(!tailn[q]) continue;
        printf("# tail of tau_k - tau_{k-1}, k=%d (n=%ld):  t   P(dtau>=t)  t^2*P\n",kidx[q],tailn[q]);
        long cum=0;
        for(int b=39;b>=0;b--){ cum+=tail[q][b]; long tt=1L<<b;
            if(cum>0 && b<=20) printf("     %8ld  %.4e  %8.2f\n",tt,(double)cum/tailn[q],(double)cum/tailn[q]*tt*tt); }
    }
    return 0;
}
