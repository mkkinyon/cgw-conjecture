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


/* s11_ztau.c -- session 11: VERIFY the idle-loop contraction theorem.
 *
 * Theorem (s11).  Let L be the event that steps 1,2 are a split of some
 * block followed by the merge of the two pieces it created (so Gamma_2 =
 * Gamma_0 and no flip occurred), lambda(x) = P_x(L).  Then for every
 * |z| <= 1 and every bounded F,
 *      |K_z F(x)| <= (1-lambda(x)) ||F|| / |1 - lambda(x) z^2| ,
 * hence  ||K_z|| <= 1/(1 + lambda_min (1-Re z^2)/(1-lambda_min)).
 * Also proved: lambda(x) >= (1/6) sum_i q_i^4 >= 1/(6 B(x)^3).
 *
 * This program measures, from a DETERMINISTIC start x:
 *   lambda(x)              (by running 2 steps and testing Gamma_2 == x)
 *   |E_x[z^tau]| = |K_z 1(x)|  for z = e^{i theta} on a grid,
 * and prints them against the theoretical bound and the a-priori bound
 * with lambda replaced by 1/(6 B^3).
 *
 * usage: s11_ztau g0 s0 cs maxsteps [seed]
 *   start: marked block (together) of mass s0 with arcs (g0, s0-g0),
 *          plus one environment block of mass 1-s0 (omitted if s0==1).
 */
static int same_state(const St *A, const St *B_){
    if(A->nb!=B_->nb) return 0;
    if((A->ub==A->vb)!=(B_->ub==B_->vb)) return 0;
    if(A->ub==A->vb){
        double x1=A->sx<A->sy?A->sx:A->sy, x2=B_->sx<B_->sy?B_->sx:B_->sy;
        if(fabs(x1-x2)>1e-12) return 0;
        if(fabs(A->sz[A->ub]-B_->sz[B_->ub])>1e-12) return 0;
    } else {
        double a1=A->sz[A->ub], b1=A->sz[A->vb];
        double a2=B_->sz[B_->ub], b2=B_->sz[B_->vb];
        double lo1=a1<b1?a1:b1, hi1=a1<b1?b1:a1;
        double lo2=a2<b2?a2:b2, hi2=a2<b2?b2:a2;
        if(fabs(lo1-lo2)>1e-12||fabs(hi1-hi2)>1e-12) return 0;
    }
    /* compare the multiset of unmarked masses */
    double u1[4096],u2[4096]; int n1=0,n2=0;
    for(int i=0;i<A->nb;i++) if(i!=A->ub&&i!=A->vb) u1[n1++]=A->sz[i];
    for(int i=0;i<B_->nb;i++) if(i!=B_->ub&&i!=B_->vb) u2[n2++]=B_->sz[i];
    if(n1!=n2) return 0;
    for(int i=0;i<n1;i++) for(int j=i+1;j<n1;j++){ if(u1[j]<u1[i]){double t=u1[i];u1[i]=u1[j];u1[j]=t;} }
    for(int i=0;i<n2;i++) for(int j=i+1;j<n2;j++){ if(u2[j]<u2[i]){double t=u2[i];u2[i]=u2[j];u2[j]=t;} }
    for(int i=0;i<n1;i++) if(fabs(u1[i]-u2[i])>1e-12) return 0;
    return 1;
}
static void set_start(St *C,double g0,double s0){
    C->nb=1; C->sz[0]=s0; C->ub=C->vb=0; C->sx=g0; C->sy=s0-g0;
    if(s0<1.0){ C->sz[C->nb++]=1.0-s0; }
}
#define NTH 9
int main(int argc,char**argv){
    if(argc<5){fprintf(stderr,"usage: %s g0 s0 cs maxsteps [seed]\n",argv[0]);return 1;}
    double g0=atof(argv[1]), s0=atof(argv[2]);
    long cs=atol(argv[3]), maxsteps=atol(argv[4]);
    unsigned seed=(argc>5)?(unsigned)atoi(argv[5]):999;
    rng_s[0]=seed*2654435761ULL+1; rng_s[1]=seed^0x9e3779b97f4a7c15ULL;
    for(int i=0;i<64;i++) rng_next();
    St C,X; st_init(&C,2*maxsteps+64); st_init(&X,2*maxsteps+64);
    set_start(&X,g0,s0);
    int B0=X.nb;
    /* (a) lambda(x) = P(Gamma_2 == x) */
    long hit=0; long lam_cs = cs*10;
    for(long t=0;t<lam_cs;t++){ st_copy(&C,&X); C.cap=X.cap; NFLIP=0; step(&C); step(&C);
        if(NFLIP==0 && same_state(&C,&X)) hit++; }
    double lam=(double)hit/lam_cs;
    /* (b) tau distribution -> E[z^tau] */
    double re[NTH]={0},im[NTH]={0}; long cens=0; double etau=0;
    for(long t=0;t<cs;t++){
        st_copy(&C,&X); NFLIP=0; long s=0;
        while(s<maxsteps && NFLIP==0){ step(&C); s++; }
        if(NFLIP==0){cens++; continue;}
        etau+=s;
        for(int q=0;q<NTH;q++){ double th=M_PI*(q+1)/(double)(NTH+1);
            re[q]+=cos(th*s); im[q]+=sin(th*s); }
    }
    long n=cs-cens;
    printf("# start: together, s0=%g arcs (%g,%g), env block %g ; B0=%d\n",s0,g0,s0-g0,1-s0,B0);
    double lam_exact=(1.0/3.0)*(pow(g0,4)+pow(s0-g0,4)+(s0<1?pow(1-s0,4):0.0));
    int nparts=B0+((X.ub==X.vb)?1:0);
    double lam_lb=1.0/(3.0*(double)nparts*nparts*nparts);
    printf("# lambda(x): measured=%.6f (%ld/%ld)  exact formula (1/3)sum r^4 = %.6f  a-priori 1/(3(B+1)^3)=%.6f\n",
        lam,hit,lam_cs,lam_exact,lam_lb);
    printf("# E[tau]=%.3f  censored=%ld/%ld\n",etau/n,cens,cs);
    printf("#  theta/pi    |E[z^tau]|    bound_measured_lambda   bound_apriori_1/(3(B+1)^3)\n");
    for(int q=0;q<NTH;q++){
        double th=M_PI*(q+1)/(double)(NTH+1);
        double R=re[q]/n, I=im[q]/n, mod=sqrt(R*R+I*I);
        double t2=1-cos(2*th);
        double bm=(1-lam)/sqrt((1-lam*cos(2*th))*(1-lam*cos(2*th))+lam*lam*sin(2*th)*sin(2*th));
        double la=lam_lb;
        double ba=(1-la)/sqrt((1-la*cos(2*th))*(1-la*cos(2*th))+la*la*sin(2*th)*sin(2*th));
        printf("   %.4f      %.5f       %.5f   %s          %.5f   %s   (1-Re z^2=%.3f)\n",
            th/M_PI,mod,bm,(mod<=bm+3e-3)?"ok":"VIOLATED",ba,(mod<=ba+3e-3)?"ok":"VIOLATED",t2);
    }
    return 0;
}
