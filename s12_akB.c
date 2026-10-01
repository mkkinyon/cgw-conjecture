/* s12_akB.c -- decompose a_k(x)=E_x[(-1)^{tau_k}] by the block count B at tau_k.
 * usage: which kmax nsamp [seed]   (which as in s12_akdecay fixed starts) */
#include "s12_state.h"
static void fixed_start(St *C,int which){
    double eps=1e-4; C->nb=0;
    switch(which){
    case 0: C->nb=1;C->sz[0]=1;C->ub=C->vb=0;C->sx=.5;C->sy=.5;break;
    case 1: C->nb=2;C->sz[0]=.6;C->sz[1]=.4;C->ub=C->vb=0;C->sx=.1;C->sy=.5;break;
    case 2: C->nb=3;C->sz[0]=.3;C->sz[1]=.3;C->sz[2]=.4;C->ub=0;C->vb=1;break;
    case 6: C->nb=9;C->sz[0]=.3;C->sz[1]=.3;for(int i=2;i<9;i++)C->sz[i]=.4/7;C->ub=0;C->vb=1;break;
    case 7: C->nb=17;C->sz[0]=.3;C->sz[1]=.3;for(int i=2;i<17;i++)C->sz[i]=.4/15;C->ub=0;C->vb=1;break;
    }
}
#define NB 8
int main(int argc,char**argv){
    int which=atoi(argv[1]),kmax=atoi(argv[2]); long ns=atol(argv[3]); unsigned seed=argc>4?atoi(argv[4]):1;
    rng_s[0]=seed*2654435761ULL+1; rng_s[1]=seed^0x9e3779b97f4a7c15ULL; for(int i=0;i<64;i++)rng_next();
    St C; st_init(&C,1<<22);
    double *a=calloc((kmax+1)*NB,8); long *n=calloc((kmax+1)*NB,sizeof(long)); double *tot=calloc(kmax+1,8);
    for(long s=0;s<ns;s++){
        fixed_start(&C,which); long steps=0;int k=0;
        while(k<kmax){ long before=NFLIP; step(&C); steps++;
            if(NFLIP!=before){ k++; int p=(steps&1)?-1:1; int b=C.nb; if(b>NB) b=NB; b--;
                a[k*NB+b]+=p; n[k*NB+b]++; tot[k]+=p; } }
    }
    printf("# which=%d nsamp=%ld: per k: a_k total, then for B=1..7,>=8: P(B_k=b)  E[(-1)^tau_k; B_k=b]\n",which,ns);
    for(int k=1;k<=kmax;k++){ printf("k=%2d a_k=%+.5f |",k,tot[k]/ns);
        for(int b=0;b<NB;b++) printf(" B%s%d: %.4f %+.5f |",b==NB-1?">=":"=",b+1,(double)n[k*NB+b]/ns,a[k*NB+b]/ns);
        printf("\n"); }
    return 0;
}
