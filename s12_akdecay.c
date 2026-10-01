/* s12_akdecay.c -- session 12: is (C2) geometric or polynomial?
 * a_k(x) = E_x[(-1)^{tau_k}] ; (C2) <=> sup_x |a_k(x)| <= C rho^k.
 * mode stat lo hi kmax nsamp [seed]:  burn 4096 from one block, then walk:
 *   at each accepted start x (gap in [lo,hi)) run TWO independent
 *   continuations to kmax flips; a_k = mean of (-1)^{tau_k} (copy 1),
 *   b_k = mean of product (both copies) = E_mu[a_k(x)^2] (sign-blind).
 *   Next start = end state of copy 1.
 * mode fixed which kmax nsamp [seed]: deterministic start
 *   which=0: one curve arcs (.5,.5); 1: together s=.6 arcs (.1,.5) env .4;
 *   2: apart a=b=.3 env .4; 3: as 0 but arcs (.5,.5-eps) env eps=1e-4 (dust);
 *   4: apart a=b=.3 env {.4-3eps, eps,eps,eps}, eps=1e-4;
 *   5: apart a=b=.3 env 40 curves of .01 (coarse fragmentation)
 */
#include "s12_state.h"
static void stat_start(St *C, int burn){
    C->nb=1; C->sz[0]=1.0; C->ub=C->vb=0; C->sx=0.5; C->sy=0.5;
    for(int s=0;s<burn;s++) step(C);
}
static double gap(const St *C){
    return C->ub==C->vb ? (C->sx<C->sy?C->sx:C->sy)
        : (C->sz[C->ub]<C->sz[C->vb]?C->sz[C->ub]:C->sz[C->vb]);
}
static void fixed_start(St *C,int which){
    double eps=1e-4;
    C->nb=0;
    switch(which){
    case 0: C->nb=1;C->sz[0]=1;C->ub=C->vb=0;C->sx=.5;C->sy=.5;break;
    case 1: C->nb=2;C->sz[0]=.6;C->sz[1]=.4;C->ub=C->vb=0;C->sx=.1;C->sy=.5;break;
    case 2: C->nb=3;C->sz[0]=.3;C->sz[1]=.3;C->sz[2]=.4;C->ub=0;C->vb=1;break;
    case 3: C->nb=2;C->sz[0]=1-eps;C->sz[1]=eps;C->ub=C->vb=0;C->sx=.5;C->sy=.5-eps;break;
    case 4: C->nb=6;C->sz[0]=.3;C->sz[1]=.3;C->sz[2]=.4-3*eps;C->sz[3]=C->sz[4]=C->sz[5]=eps;C->ub=0;C->vb=1;break;
    case 5: C->nb=42;C->sz[0]=.3;C->sz[1]=.3;for(int i=2;i<42;i++)C->sz[i]=.01;C->ub=0;C->vb=1;break;
    }
}
/* run to kmax flips; fill par[k] = (-1)^{tau_k}; return steps */
static long run_flips(St *C,int kmax,int *par,long maxsteps){
    long steps=0;int k=0;
    while(k<kmax && steps<maxsteps){ long before=NFLIP; step(C); steps++;
        if(NFLIP!=before){ k++; par[k]=(steps&1)?-1:1; } }
    for(int j=k+1;j<=kmax;j++) par[j]=0;
    return steps;
}
int main(int argc,char**argv){
    if(argc<5){fprintf(stderr,"usage: stat lo hi kmax nsamp [seed] | fixed which kmax nsamp [seed]\n");return 1;}
    int fixed=!strcmp(argv[1],"fixed");
    double lo=0,hi=0;int which=0,kmax;long ns;unsigned seed;
    if(fixed){which=atoi(argv[2]);kmax=atoi(argv[3]);ns=atol(argv[4]);seed=argc>5?atoi(argv[5]):1;}
    else{lo=atof(argv[2]);hi=atof(argv[3]);kmax=atoi(argv[4]);ns=atol(argv[5]);seed=argc>6?atoi(argv[6]):1;}
    rng_s[0]=seed*2654435761ULL+1; rng_s[1]=seed^0x9e3779b97f4a7c15ULL; for(int i=0;i<64;i++)rng_next();
    long maxsteps=4000000;
    St C,C1,C2; st_init(&C,2*(8192+maxsteps)+64); st_init(&C1,C.cap); st_init(&C2,C.cap);
    double *a=calloc(kmax+1,8),*b=calloc(kmax+1,8); long *n=calloc(kmax+1,sizeof(long));
    int *p1=calloc(kmax+2,sizeof(int)),*p2=calloc(kmax+2,sizeof(int));
    long tried=0,done=0;
    if(!fixed) stat_start(&C,4096);
    while(done<ns){
        if(fixed) fixed_start(&C,which);
        else { double g=gap(&C); tried++; if(g<lo||g>=hi){ step(&C); continue; } }
        st_copy(&C1,&C); st_copy(&C2,&C);
        run_flips(&C1,kmax,p1,maxsteps); run_flips(&C2,kmax,p2,maxsteps);
        for(int k=1;k<=kmax;k++){ if(p1[k]&&p2[k]){ a[k]+=p1[k]; b[k]+=p1[k]*p2[k]; n[k]++; } }
        done++;
        if(!fixed){ st_copy(&C,&C1); for(int s=0;s<32;s++) step(&C); }
    }
    printf("# %s %s nsamp=%ld seed=%u\n",argv[1],argv[2],ns,seed);
    printf("#  k    n        a_k=E(-1)^tau_k     se        b_k=E[a_k(x)^2]   se\n");
    for(int k=1;k<=kmax;k++){ double ak=a[k]/n[k], bk=b[k]/n[k];
        printf("  %3d %9ld   %+.5f   %.5f    %+.5f   %.5f\n",k,n[k],ak,1/sqrt((double)n[k]),bk,sqrt((1-bk*bk)/n[k])); }
    return 0;
}
