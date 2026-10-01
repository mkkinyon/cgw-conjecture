// s10_nomerge.c : parity decay E[(-1)^N_t] for the marked pair (a,b) with NO merges
// (adversarial dust environment).  Moves per step: shrink a (prob a^2), shrink b (b^2),
// flip (2ab), nothing (else).  Shrink: a -> a*sqrt(U).  Flip: (a,b)->(al+be, s-al-be).
// usage: ./s10_nomerge a0 b0 tmax samples [seed]
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
static uint64_t st=88172645463325252ULL;
static inline uint64_t xr(void){st^=st<<7;st^=st>>9;return st;}
static inline double ur(void){return (xr()>>11)*(1.0/9007199254740992.0);}
int main(int argc,char**argv){
  double a0=atof(argv[1]),b0=atof(argv[2]); long tmax=atol(argv[3]); long S=atol(argv[4]);
  if(argc>5) st^=atoll(argv[5])*0x9E3779B97F4A7C15ULL;
  int nb=0; long tb[64]; for(long t=1;t<=tmax;t*=2) tb[nb++]=t;
  double *acc=calloc(nb,sizeof(double)); double *acc2=calloc(nb,sizeof(double));
  for(long s=0;s<S;s++){
    double a=a0,b=b0; int N=0; int ib=0;
    for(long t=1;t<=tmax;t++){
      double u=ur();
      if(u<a*a){ a*=sqrt(ur()); }
      else if(u<a*a+b*b){ b*=sqrt(ur()); }
      else if(u<a*a+b*b+2*a*b){ double al=a*ur(), be=b*ur(); double ss=a+b; a=al+be; b=ss-a; N++; }
      if(t==tb[ib]){ acc[ib]+=(N&1)?-1:1; acc2[ib]+=1; ib++; if(ib==nb)break; }
    }
  }
  printf("# a0=%g b0=%g samples=%ld\n",a0,b0,S);
  for(int i=0;i<nb;i++){ double m=acc[i]/S; printf("t=%8ld  Phi=%+.5f  se=%.5f  t*Phi=%+.4f  t^2*Phi=%+.4f\n",tb[i],m,sqrt((1-m*m)/S),tb[i]*m,(double)tb[i]*tb[i]*m); }
  return 0;
}
