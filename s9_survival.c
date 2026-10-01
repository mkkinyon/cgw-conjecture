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


/* s9_survival.c -- session 9: quotient-chain diagnostics (flip counting).
 * Modes:
 *  flux j cs [seed]      P_{mu_j}(g<=eps) for eps=2^-b, vs bound 4 eps^2 j
 *  tail g0 tmax cs [seed]  first-flip tail P(tau>t) from start
 *                        together, sides (g0, 1/2-g0)?? -- we use the
 *                        start: single block, x=g0, y=1-g0 (clean env).
 *  tailenv g0 tmax cs [seed] start from stationary env (burn 4096),
 *                        then force gap: reset marked coords to
 *                        (g0, s-g0) keeping env.
 *  phi g0 d cs [seed]    Phi_d = E[(-1)^N_d] and P(tau>d) from the
 *                        stationary-env start with forced gap g0.
 */
static void stat_start(St *C, int burn){
    C->nb=1; C->sz[0]=1.0; C->ub=C->vb=0; C->sx=0.5; C->sy=0.5;
    for(int s=0;s<burn;s++) step(C);
}
/* force the marked geometry to gap g0 keeping total s and env */
static void force_gap(St *C, double g0){
    if(C->ub==C->vb){ double s=C->sz[C->ub]; if(g0<s){C->sx=g0; C->sy=s-g0;} }
    else { double s=C->sz[C->ub]+C->sz[C->vb]; if(g0<s){C->sz[C->ub]=g0; C->sz[C->vb]=s-g0;} }
}
static double gap(const St *C){
    return C->ub==C->vb ? (C->sx<C->sy?C->sx:C->sy)
        : (C->sz[C->ub]<C->sz[C->vb]?C->sz[C->ub]:C->sz[C->vb]);
}
static double stot(const St *C){
    return C->ub==C->vb ? C->sz[C->ub] : C->sz[C->ub]+C->sz[C->vb];
}
int main(int argc,char**argv){
    const char*mode=argc>1?argv[1]:"flux";
    unsigned seed=777;
    if(argc>5 && !strcmp(mode,"tail")) seed=(unsigned)atoi(argv[5]);
    else if(argc>5 && !strcmp(mode,"tailenv")) seed=(unsigned)atoi(argv[5]);
    else if(argc>6 && !strcmp(mode,"phic")) seed=(unsigned)atoi(argv[6]);
    else if(argc>5 && !strcmp(mode,"phi")) seed=(unsigned)atoi(argv[5]);
    else if(argc>4 && !strcmp(mode,"flux")) seed=(unsigned)atoi(argv[4]);
    rng_s[0]=seed*2654435761ULL+1; rng_s[1]=seed^0x9e3779b97f4a7c15ULL;
    for(int i=0;i<64;i++) rng_next();
    if(!strcmp(mode,"flux")){
        int j=atoi(argv[2]); long cs=atol(argv[3]);
        St C; st_init(&C,2*j+64);
        long hist[40]={0}; double sm[40]={0}; long shist[40]={0};
        for(long t=0;t<cs;t++){
            C.nb=1;C.sz[0]=1;C.ub=C.vb=0;C.sx=C.sy=0.5;
            for(int s=0;s<j;s++) step(&C);
            double g=gap(&C), s=stot(&C);
            int b=(int)(-log2(g>1e-11?g:1e-11)); if(b<0)b=0; if(b>39)b=39;
            hist[b]++;
            int bs=(int)(-log2(s>1e-11?s:1e-11)); if(bs<0)bs=0; if(bs>39)bs=39;
            shist[bs]++;
        }
        printf("flux j=%d cs=%ld:  eps=2^-b   P(g<=eps)   4eps^2 j   ratio | P(s<=eps)  2eps^2 j\n",j,cs);
        long cum=0, scum=0;
        for(int b=39;b>=0;b--){ cum+=hist[b]; scum+=shist[b];
            double eps=pow(2.0,-b); double bd=4*eps*eps*j, sbd=2*eps*eps*j;
            if(cum>0 && b>=1) printf("  %2d  %.3e  %.3e  %.3f | %.3e  %.3e\n",b,(double)cum/cs,bd,(double)cum/cs/bd,(double)scum/cs,sbd);
        }
    } else if(!strcmp(mode,"tail")||!strcmp(mode,"tailenv")){
        double g0=atof(argv[2]); int tmax=atoi(argv[3]); long cs=atol(argv[4]);
        int env=!strcmp(mode,"tailenv");
        St C; st_init(&C,2*tmax+8192+64);
        long *surv=calloc(tmax+1,sizeof(long));
        for(long t=0;t<cs;t++){
            if(env){ stat_start(&C,4096); force_gap(&C,g0); }
            else { C.nb=1;C.sz[0]=1;C.ub=C.vb=0;C.sx=g0;C.sy=1-g0; }
            NFLIP=0; int s=0;
            for(;s<tmax && NFLIP==0;s++) step(&C);
            /* tau = s if flipped at step s (1-based), else > tmax */
            if(NFLIP==0) surv[tmax]++; else surv[s-1]++; /* tau=s */
        }
        /* P(tau>t) */
        long cum=0; printf("%s g0=%g cs=%ld:  t  P(tau>t)  t*P  (g0 t)*P\n",mode,g0,cs);
        for(int t=tmax;t>=0;t--){ cum+=surv[t];
            if(t==tmax||(t&(t-1))==0) if(t>0) printf("  %6d  %.4e  %.4f  %.4f\n",t,(double)cum/cs,(double)cum/cs*t,(double)cum/cs*t*g0);
        }
        /* note: cum at index t counts tau>=t+1 ... fix: surv[s-1] means tau=s => tau>t iff s>t iff index>=t */
    } else if(!strcmp(mode,"phi")){
        double g0=atof(argv[2]); int d=atoi(argv[3]); long cs=atol(argv[4]);
        St C; st_init(&C,2*d+8192+64);
        double acc=0; long notau=0;
        for(long t=0;t<cs;t++){
            stat_start(&C,4096); force_gap(&C,g0);
            NFLIP=0; for(int s=0;s<d;s++) step(&C);
            acc += (NFLIP&1)?-1:1; notau += (NFLIP==0);
        }
        printf("phi g0=%g d=%d cs=%ld: Phi=%.5f (se %.5f)  P(tau>d)=%.5f  ratio=%.3f  g0*d=%.2f\n",
            g0,d,cs,acc/cs,1.0/sqrt((double)cs),(double)notau/cs,acc/cs/((double)notau/cs+1e-30),g0*d);
    } else if(!strcmp(mode,"phic")){
        /* phic lo hi d cs [seed]: stationary start (burn 4096), accept if
         * gap in [lo,hi); report Phi_d, P(tau>d), Phi restricted to tau<=d,
         * tail P(tau>t) at powers of 2, and E[s]. */
        double lo=atof(argv[2]), hi=atof(argv[3]); int d=atoi(argv[4]); long cs=atol(argv[5]);
        if(argc>6){ seed=(unsigned)atoi(argv[6]); rng_s[0]=seed*2654435761ULL+1; rng_s[1]=seed^0x9e3779b97f4a7c15ULL; for(int i=0;i<64;i++) rng_next(); }
        St C; st_init(&C,2*d+8192+64);
        double acc=0, acc_after=0, es=0; long notau=0, tried=0;
        long *surv=calloc(d+2,sizeof(long));
        for(long t=0;t<cs;){
            stat_start(&C,4096); tried++;
            double g=gap(&C); if(g<lo||g>=hi) continue;
            t++; es+=stot(&C);
            NFLIP=0; int tau=d+1;
            for(int s=0;s<d;s++){ step(&C); if(tau>d && NFLIP) tau=s+1; }
            double sg=(NFLIP&1)?-1:1; acc+=sg;
            if(NFLIP==0) notau++; else acc_after+=sg;
            surv[tau<=d?tau:d+1]++;
        }
        printf("phic gap in [%g,%g) d=%d cs=%ld (accept %.4f) E[s]=%.3f: Phi=%.5f (se %.5f) P(tau>d)=%.5f Phi|tau<=d part=%.5f  g_lo*d=%.2f\n",
            lo,hi,d,cs,(double)cs/tried,es/cs,acc/cs,1.0/sqrt((double)cs),(double)notau/cs,acc_after/cs,lo*d);
        long cum=0; for(int t=d+1;t>=1;t--){ cum+=surv[t]; int tt=t-1; if(tt>0 && (tt&(tt-1))==0) printf("   P(tau>%d)=%.4e  t^2P=%.3f  (lo t)^2 P=%.4f\n",tt,(double)cum/cs,(double)cum/cs*tt*tt,(double)cum/cs*tt*tt*lo*lo); }
    } else if(!strcmp(mode,"oddeven")){
        /* oddeven d cs fs [seed]: stationary start; A runs d steps, B runs d
         * steps (even total 2d) or d+1 steps (odd total 2d+1):
         * reports <P^d eps,P^d eps> and <P^d eps,P^{d+1} eps>. */
        int d=atoi(argv[2]); long cs=atol(argv[3]); long fs=atol(argv[4]);
        if(argc>5){ seed=(unsigned)atoi(argv[5]); rng_s[0]=seed*2654435761ULL+1; rng_s[1]=seed^0x9e3779b97f4a7c15ULL; for(int i=0;i<64;i++) rng_next(); }
        St C,A,B; int cap=2*(4096+d)+64; st_init(&C,cap); st_init(&A,cap); st_init(&B,cap);
        double ev=0, od=0;
        for(long t=0;t<cs;t++){
            stat_start(&C,4096);
            double ie=0, io=0;
            for(long r=0;r<fs;r++){
                st_copy(&A,&C); st_copy(&B,&C);
                for(int s=0;s<d;s++){ step(&A); step(&B); }
                int ea=eps_st(&A), eb=eps_st(&B);
                ie+=ea*eb; step(&B); io+=ea*eps_st(&B);
            }
            ev+=ie/fs; od+=io/fs;
        }
        printf("oddeven d=%d cs=%ld fs=%ld: <P^d e,P^d e>=%.4e  <P^d e,P^{d+1} e>=%.4e  (se floor %.1e)  ratio odd/even=%.3f\n",
            d,cs,fs,ev/cs,od/cs,1.0/sqrt((double)cs*fs),od/ev);
    }
    return 0;
}
