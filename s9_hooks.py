"""s9_hooks.py -- exact finite-n formula for E[sum c_i^2] of sigma_m = gamma * tau_1...tau_m
(gamma uniform n-cycle, tau_i uniform transpositions), via hook characters:
  E F(sigma_m) = sum_k d_k rho_k r_k^m <chi_k, F>,  hooks (n-k,1^k):
  d_k = C(n-1,k), rho_k = (-1)^k, r_k = 1 - 2k/(n-1),
  <chi_k,F> = sum_mu chi_k(mu) F(mu) / z_mu,  F(mu)=sum mu_i^2,
  sum_k chi_k(mu) t^k = prod_i (1-(-t)^mu_i) / (1+t).
Output: E[sum p_i^2] = E[sum c_i^2]/n^2 for m=0..M, compared with 1/2 + ... """
import sys
from fractions import Fraction
from math import comb, factorial
from collections import Counter

def partitions(n, mx=None):
    if mx is None: mx=n
    if n==0: yield (); return
    for k in range(min(n,mx),0,-1):
        for p in partitions(n-k,k): yield (k,)+p

def polymul(a,b,deg):
    c=[0]*(deg+1)
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if i+j>deg: break
            c[i+j]+=x*y
    return c

def hook_chars(mu,n):
    # numerator prod (1-(-t)^mu_i), then divide by (1+t)
    num=[1]+[0]*n
    for m in mu:
        f=[0]*(n+1); f[0]=1; f[m]=-((-1)**m)
        num=polymul(num,f,n)
    # divide by 1+t: q with (1+t)q = num
    q=[0]*n
    q[0]=num[0]
    for i in range(1,n): q[i]=num[i]-q[i-1]
    assert num[n]-q[n-1]==0
    return q  # q[k] = chi_{(n-k,1^k)}(mu)

def zmu(mu):
    z=1
    for part,c in Counter(mu).items(): z*=part**c*factorial(c)
    return z

def inner(n):
    ip=[Fraction(0)]*n
    for mu in partitions(n):
        ch=hook_chars(mu,n); F=sum(x*x for x in mu); z=zmu(mu)
        for k in range(n): ip[k]+=Fraction(ch[k]*F, z)
    return ip

def E_sum_c2(n,M):
    ip=inner(n)
    out=[]
    for m in range(M+1):
        s=Fraction(0)
        for k in range(n):
            r=Fraction(n-1-2*k, n-1)
            s+=(-1)**k*r**m*ip[k]
        out.append(s/n**2)
    return out

if __name__=='__main__':
    n=int(sys.argv[1]); M=int(sys.argv[2])
    res=E_sum_c2(n,M)
    for m,v in enumerate(res):
        print(f"n={n} m={m:3d}  E[sum p^2]={float(v):.6f}   -1/2 = {float(v)-0.5:+.6f}  (m^2*(..)={(float(v)-0.5)*m*m:+.4f})")
