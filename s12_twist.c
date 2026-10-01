/* s12_twist.c -- factorisation check for (C2''): from a fixed start, compare
 *   t_k(w)=E[w^{tau_k}(-1)^{tau_k}]  with  E[w^{tau_k}] * E[(-1)^{tau_k}],  w in {0.99,0.95,0.9}
 * usage: which kmax nsamp [seed]  (which: 0 one curve, 2 three-curve apart) */
#include "s12_state.h"
static void fixed_start(St *C,int which){
    switch(which){
    case 0: C->nb=1;C->sz[0]=1;C->ub=C->vb=0;C->sx=.5;C->sy=.5;break;
    case 2: C->nb=3;C->sz[0]=.3;C->sz[1]=.3;C->sz[2]=.4;C->ub=0;C->vb=1;break;
    }
}
#define NW 3
int main(int argc,char**argv){
    int which=atoi(argv[1]),kmax=atoi(argv[2]); long ns=atol(argv[3]); unsigned seed=argc>4?atoi(argv[4]):1;
    double W[NW]={0.99,0.95,0.90};
    rng_s[0]=seed*2654435761ULL+1; rng_s[1]=seed^0x9e3779b97f4a7c15ULL; for(int i=0;i<64;i++)rng_next();
    St C; st_init(&C,1<<22);
    double *a=calloc(kmax+1,8), *ew=calloc((kmax+1)*NW,8), *tw=calloc((kmax+1)*NW,8);
    for(long s=0;s<ns;s++){
        fixed_start(&C,which); long steps=0;int k=0;
        while(k<kmax){ long before=NFLIP; step(&C); steps++;
            if(NFLIP!=before){ k++; int p=(steps&1)?-1:1; a[k]+=p;
                for(int j=0;j<NW;j++){ double wt=pow(W[j],(double)steps); ew[k*NW+j]+=wt; tw[k*NW+j]+=wt*p; } } }
    }
    printf("# twist which=%d nsamp=%ld se~%.1e\n#  k   a_k      | per w: E[w^tau (-1)^tau]  E[w^tau]*a_k   E[w^tau]\n",which,ns,1/sqrt((double)ns));
    for(int k=1;k<=kmax;k++){ printf("  %2d %+.5f |",k,a[k]/ns);
        for(int j=0;j<NW;j++) printf("  w=%.2f: %+.5f  %+.5f  %.4f |",W[j],tw[k*NW+j]/ns,(ew[k*NW+j]/ns)*(a[k]/ns),ew[k*NW+j]/ns);
        printf("\n"); }
    return 0;
}
