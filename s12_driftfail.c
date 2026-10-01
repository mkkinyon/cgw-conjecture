/* s12_driftfail.c -- TODO 12a obstruction check: from x = (a=m, b=eta apart; env = N pieces
 * of mass (1-m-eta)/N), measure R = sum r^4(Gamma_tau)/sum r^4(x) and E[R^{-theta}]. */
#include "s12_state.h"
static double S4(const St*C){ double s=0; for(int i=0;i<C->nb;i++){ if(i==C->ub&&i==C->vb){ s+=pow(C->sx,4)+pow(C->sy,4);} else s+=pow(C->sz[i],4);} return s; }
int main(int argc,char**argv){
    double m=atof(argv[1]), eta=atof(argv[2]); int N=atoi(argv[3]); long ns=atol(argv[4]);
    rng_s[0]=77; rng_s[1]=0x9e3779b97f4a7c15ULL; for(int i=0;i<64;i++)rng_next();
    St C; st_init(&C, 4*N+1000000);
    double sumR=0, sumV1=0, sumV3=0, sumtau=0, sums=0; double th1=0.1, th3=0.3;
    for(long s=0;s<ns;s++){
        C.nb=N+2; C.sz[0]=m; C.sz[1]=eta; for(int i=2;i<N+2;i++) C.sz[i]=(1-m-eta)/N; C.ub=0; C.vb=1;
        double S0=S4(&C); long steps=0; long before=NFLIP;
        while(NFLIP==before){ step(&C); steps++; }
        double S1=S4(&C), R=S1/S0; sumR+=R; sumV1+=pow(R,-th1); sumV3+=pow(R,-th3); sumtau+=steps; sums+=C.sz[C.ub];
    }
    printf("m=%g eta=%g N=%d ns=%ld: E[tau]=%.1f E[s_tau]=%.4f E[R]=%.4f E[R^-0.1]=%.3f E[R^-0.3]=%.3f  (drift needs E[R^-theta]<1)\n",m,eta,N,ns,sumtau/ns,sums/ns,sumR/ns,sumV1/ns,sumV3/ns);
    return 0;
}
