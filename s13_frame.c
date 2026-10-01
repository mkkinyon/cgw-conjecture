/* s13_frame.c -- does the FRAME (one block vs random permutation) explain the
 * real-square bias constant?   Session 13 scouting of branch (II).
 *
 * Discrete model on N slots.  Frame pi_0 is either
 *   mode 0: the N-cycle, u = 0, v = N/2 (the antipodal one-block ensemble of
 *           bias_scale.c / s7_twopoint.c),
 *   mode 1: the N-cycle, u = 0, v uniform (one block, stationary gap law),
 *   mode 2: a uniform random permutation of N slots, u = 0, v uniform
 *           (the PD(1)-like frame of a real orbit: the cycle structure of the
 *           pair permutation pi_P with rows 1,2 marked).
 * Pool: k chords = k uniform random pairs of distinct slots, pairwise disjoint
 * and avoiding u,v (resampled otherwise).  For S a subset, pi_S = pi_0 o tau_S
 * and eps_S = +1 iff u ~ v in pi_S.
 *   4 E[bias^2] = E[eps_S eps_T]   (S,T independent uniform subsets)
 *   E|bias|      from ssamp subsets per instance.
 * usage: ./s13_frame mode k ninst [ssamp] [N] [seed]
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef unsigned long long u64;
static u64 rng_s[2];
static inline u64 rng_next(void){u64 s1=rng_s[0],s0=rng_s[1];rng_s[0]=s0;s1^=s1<<23;rng_s[1]=s1^s0^(s1>>18)^(s0>>5);return rng_s[1]+s0;}
static inline double rng_u(void){return (rng_next()>>11)*(1.0/9007199254740992.0);}
static int *pi0,*tau,*used; static int *ca,*cb;
static int N;
static int same_cycle(int u,int v){ /* in pi0 o tau */
    int x=u; do{ x=pi0[tau[x]]; if(x==v) return 1; }while(x!=u); return 0; }
int main(int argc,char**argv){
    int mode=atoi(argv[1]),k=atoi(argv[2]); long ninst=atol(argv[3]);
    int ssamp=argc>4?atoi(argv[4]):0; N=argc>5?atoi(argv[5]):4096; unsigned seed=argc>6?atoi(argv[6]):7;
    rng_s[0]=seed*2654435761ULL+11; rng_s[1]=seed^0x9e3779b97f4a7c15ULL; for(int i=0;i<64;i++)rng_next();
    pi0=malloc(N*sizeof(int)); tau=malloc(N*sizeof(int)); used=malloc(N*sizeof(int));
    ca=malloc(k*sizeof(int)); cb=malloc(k*sizeof(int));
    double sum2=0, sumabs=0, sumsgn=0; long n2=0;
    for(long it=0;it<ninst;it++){
        int u=0,v;
        if(mode==2){ for(int i=0;i<N;i++)pi0[i]=i; for(int i=N-1;i>0;i--){int j=rng_next()%(i+1);int t=pi0[i];pi0[i]=pi0[j];pi0[j]=t;} }
        else { for(int i=0;i<N;i++)pi0[i]=(i+1)%N; }
        v=(mode==0)?N/2:1+(int)(rng_next()%(N-1));
        for(int i=0;i<N;i++){used[i]=0;tau[i]=i;} used[u]=used[v]=1;
        for(int c=0;c<k;c++){ int a,b; do{a=rng_next()%N;}while(used[a]); used[a]=1; do{b=rng_next()%N;}while(used[b]); used[b]=1; ca[c]=a;cb[c]=b; }
        /* two independent subsets */
        int e[2];
        for(int r=0;r<2;r++){ for(int c=0;c<k;c++){ if(rng_next()&1){tau[ca[c]]=cb[c];tau[cb[c]]=ca[c];} else {tau[ca[c]]=ca[c];tau[cb[c]]=cb[c];} }
            e[r]=same_cycle(u,v)?1:-1; }
        sum2+=e[0]*e[1]; n2++;
        if(ssamp>0){ double acc=0; for(int s=0;s<ssamp;s++){ for(int c=0;c<k;c++){ if(rng_next()&1){tau[ca[c]]=cb[c];tau[cb[c]]=ca[c];} else {tau[ca[c]]=ca[c];tau[cb[c]]=cb[c];} } acc+=same_cycle(u,v)?1:-1; }
            double bias=acc/(2.0*ssamp); sumabs+=fabs(bias); sumsgn+=bias; }
    }
    double m2=sum2/n2; double se=sqrt((1-m2*m2)/n2);
    printf("mode=%d k=%d N=%d ninst=%ld: 4E[bias^2]=%.5f(%.5f)  E[bias^2]*k=%.4f(%.4f)",mode,k,N,ninst,m2,se,m2*k/4,se*k/4);
    if(ssamp>0) printf("  E|bias|=%.4f  mean=%+.4f (ssamp=%d, MC floor ~%.3f)",sumabs/ninst,sumsgn/ninst,ssamp,0.5/sqrt(ssamp));
    printf("\n"); return 0;
}
