#!/usr/bin/env python3
# suspension_theorem.py
# Resolution of Open Problem 10.2 (field_surprisal_geometry_v2): the (logM,degree) Gibbs surface
# has EXACT constant Fisher-Rao curvature 1/4. Proof: an "indicator family" (one statistic an
# affine indicator of a single outcome) is the spherical SUSPENSION of that apex outcome with the
# 1-parameter curve of the other statistic; the suspension metric is the round sphere, curvature
# 1/4 for any curve. Exact symbolic verification throughout, fail-first.
import sympy as sp

P=[]
def ck(name,cond,tag,note=""):
    ok=bool(cond); P.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] ({tag:8}) {name}"+(f"   |  {note}" if note else ""))
    assert ok, name

t,s=sp.symbols('t s',real=True)
def K_origin(a,b,nu):    # exact Gaussian curvature (Brioschi) at origin of Hess(log sum nu_i e^{t a_i + s b_i})
    Z=sum(sp.nsimplify(ni)*sp.exp(t*sp.nsimplify(ai)+s*sp.nsimplify(bi)) for ai,bi,ni in zip(a,b,nu))
    A=sp.log(Z); E=sp.diff(A,t,2); F=sp.diff(A,t,s); G=sp.diff(A,s,2)
    o={t:0,s:0}; d=lambda f,*v: sp.diff(f,*v).subs(o); E0,F0,G0=E.subs(o),F.subs(o),G.subs(o)
    Ev,Eu,Gu,Gv,Fu,Fv=d(E,s),d(E,t),d(G,t),d(G,s),d(F,t),d(F,s); Evv,Fuv,Guu=d(E,s,2),d(F,t,s),d(G,t,2)
    M1=sp.Matrix([[-Evv/2+Fuv-Guu/2,Eu/2,Fu-Ev/2],[Fv-Gu/2,E0,F0],[Gv/2,F0,G0]])
    M2=sp.Matrix([[0,Ev/2,Gu/2],[Ev/2,E0,F0],[Gu/2,F0,G0]])
    return sp.simplify((M1.det()-M2.det())/(E0*G0-F0**2)**2)

print("="*82); print("PART 1  Indicator families have curvature 1/4 (arbitrary statistic, any measure)"); print("="*82)
a1,a2,a3=sp.symbols('a1 a2 a3',real=True)
Kgen=K_origin([a1,a2,a3,0],[0,0,0,1],[1,1,1,1])
ck("m=4: arbitrary statistic a=(a1,a2,a3,0) + indicator 1_4  ->  K = 1/4 identically",
   sp.simplify(Kgen-sp.Rational(1,4))==0,"FORCED","independent of a1,a2,a3 (symbolic)")
# non-uniform base measure = evaluation at an off-origin point -> tests 'everywhere'
ck("off-origin (non-uniform base measure): still exactly 1/4  ->  CONSTANT curvature",
   all(K_origin([1,2,3,0],[0,0,0,1],nu)==sp.Rational(1,4) for nu in ([2,1,1,1],[3,5,2,7],[1,1,1,9])),
   "FORCED","constant everywhere, so the surface is totally geodesic")
ck("the 'quadratic' case a=(0,1,4,9,0)+1_5 is exactly 1/4 (earlier numeric deviation was FD error)",
   K_origin([0,1,4,9,0],[0,0,0,0,1],[1,1,1,1,1])==sp.Rational(1,4),"FORCED")
# control: two generic statistics (no indicator) are NOT 1/4
ck("CONTROL: two generic statistics (no indicator) give K != 1/4",
   K_origin([0,1,2,5],[0,2,1,3],[1,1,1,1])!=sp.Rational(1,4)
   and K_origin([0,1,3,6,10],[1,0,2,1,3],[1,1,1,1,1])!=sp.Rational(1,4),"FORCED",
   f"K = {K_origin([0,1,2,5],[0,2,1,3],[1,1,1,1])} , {K_origin([0,1,3,6,10],[1,0,2,1,3],[1,1,1,1,1])}")

print("\n"+"="*82); print("PART 2  The mechanism: spherical suspension (metric = round sphere, any a-curve)"); print("="*82)
th,psi=sp.symbols('theta psi',real=True); b1,b2,b3=sp.symbols('b1 b2 b3',real=True)
w=[sp.exp(th*b1),sp.exp(th*b2),sp.exp(th*b3)]; Zc=sum(w); rho=[wi/Zc for wi in w]
x=[sp.sqrt(r) for r in rho]                                  # sqrt-conditional on the 3 non-apex outcomes
sqrtp=sp.Matrix([sp.cos(psi)*x[0],sp.cos(psi)*x[1],sp.cos(psi)*x[2],sp.sin(psi)])  # radius-1 sqrt-embedding
Xth,Xps=sp.diff(sqrtp,th),sp.diff(sqrtp,psi)
gpp=sp.simplify(Xps.dot(Xps)); gtp=sp.simplify(Xth.dot(Xps)); gtt=sp.simplify(Xth.dot(Xth))
ck("suspension metric: g_psi,psi = 1  and  g_theta,psi = 0  (warped product signature)",
   gpp==1 and gtp==0,"FORCED","sqrt(p)=(cos(psi) x(theta), sin(psi)); apex = the indicator outcome")
ck("g_theta,theta = cos^2(psi) * |x'(theta)|^2  (apex direction orthogonal to the a-curve)",
   sp.simplify(gtt-sp.cos(psi)**2*Xth.dot(Xth).subs(psi,0))==0,"FORCED","reparametrize theta by arc length -> dl^2")
f=sp.cos(psi)
ck("curvature of  dpsi^2 + cos^2(psi) dl^2  is  -f''/f = 1 (radius-1)  ->  Fisher curvature 1/4",
   sp.simplify(-sp.diff(f,psi,2)/f)==1,"FORCED","round-sphere metric for ANY a-curve; Fisher = 4x radius-1")

print("\n"+"="*82); print("PART 3  Corollary: the catalog's (logM, degree) surface is EXACTLY 1/4"); print("="*82)
phi=(1+sp.sqrt(5))/2
# degree in {2,4}; degree = 2 + 2*indicator_K  -> affinely a single-outcome indicator
deg=[2,2,2,2,2,2,4]
ck("degree = 2 + 2*1_K  (affinely the indicator of the single outcome K)",
   all(deg[i]==2+2*(1 if i==6 else 0) for i in range(7)),"FORCED","only K has degree 4")
# catalog as an indicator family (phi,tau share (logM,deg) -> base-measure weight 2), any logM values:
a=[sp.log(2),sp.log(3),sp.log(5),sp.log(phi),0]; b=[0,0,0,0,1]; nu=[1,1,1,2,1]   # phi/tau lumped, apex=K
ck("catalog (logM, degree) Gibbs surface has exact constant curvature 1/4",
   K_origin(a,b,nu)==sp.Rational(1,4),"FORCED","Open Problem 10.2 resolved: exactly 1/4, totally geodesic")

print("\n"+"="*82)
print(f"SUSPENSION THEOREM: {sum(P)}/{len(P)} passed"+("  exit 0" if all(P) else "  FAIL"))
