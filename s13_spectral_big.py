"""Exact spectral measure of the same-cycle sign (session 13), closed form:
  lambda=(a,b,1^j), a>=b>=1, a+b+j=N:
  c_lambda = (a-b+1)/((a+j+1)(b+j)),  r_lambda = [a(a-1)+b(b-3)-j(j+3)]/(N(N-1)),
  nu_lambda = 2 c_lambda^2 (1-r_lambda);  sum nu = 1.
  E_pi[Phi_m] = <eps,P^m eps> = sum nu r^m ;  psi_pi(d) = ||P^d eps||^2 = E_pi[Phi_{2d}].
usage: python3 s13_spectral_big.py N [m1 m2 ...]"""
import sys, numpy as np
def measure(N):
    b=np.arange(1,N,dtype=np.float64)[:,None]; j=np.arange(0,N,dtype=np.float64)[None,:]
    a=N-b-j; ok=a>=b
    c=np.where(ok,(a-b+1)/((a+j+1)*(b+j)),0.0)
    r=np.where(ok,(a*(a-1)+b*(b-3)-j*(j+3))/(N*(N-1)),0.0)
    nu=np.where(ok,2*c*c*(1-r),0.0)
    return r[ok],nu[ok]
if __name__=="__main__":
    N=int(sys.argv[1]); ms=[int(x) for x in sys.argv[2:]] or [1,2,4,8,16,32,64,128,256,512,1024,2048]
    r,nu=measure(N); print("N=%d shapes=%d sum nu=%.12f"%(N,len(nu),nu.sum()))
    for m in ms:
        if m>N: break
        phi=(nu*r**m).sum(); print("  m=%5d E_pi[Phi_m]=%+.7f  m*E=%.4f"%(m,phi,m*phi))
    for delta in [0.2,0.1,0.05,0.02,0.01,0.005,0.002]:
        m1=nu[r>=1-delta].sum(); m2=nu[r<=-1+delta].sum()
        print("  delta=%.3f nu[1-d,1]=%.5f ratio=%.3f  nu[-1,-1+d]=%.2e"%(delta,m1,m1/delta,m2))
