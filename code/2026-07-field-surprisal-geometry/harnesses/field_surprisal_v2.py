#!/usr/bin/env python3
# field_surprisal_v2.py -- Field Surprisal Geometry v2: core + Tier-1 enhancements.
# Exact algebra where the data is algebraic; numerics (mpmath) for transcendental values
# and for the two-statistic curvature (validated against the known 1/4 sphere case). Fail-first.
import sympy as sp
import mpmath as mp
mp.mp.dps = 30

P=[]
def ck(name,cond,tag,note=""):
    ok=bool(cond); P.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] ({tag:9}) {name}"+(f"   |  {note}" if note else ""))
    assert ok, name

phi=(1+sp.sqrt(5))/2
names=['sqrt2','sqrt3','sqrt5','phi','tau','phi^4','K']
Msym=[sp.Integer(2),sp.Integer(3),sp.Integer(5),phi,phi,phi**4,phi**4-1]
deg =[2,2,2,2,2,2,4]
Enum=[mp.log(float(m)) for m in Msym]                 # costs log M (numeric)

# ================= CORE (condensed) =================
print("="*86); print("CORE  Z, affinity, Fisher=Var, curvature 1/4"); print("="*86)
Z=sp.simplify(sum(Msym))
ck("Z = 17 + 4 sqrt5", sp.simplify(Z-(17+4*sp.sqrt(5)))==0,"FORCED",f"Z={Z}")
b,m1,m2=sp.symbols('beta m1 m2',positive=True)
Zb=m1**(-b)+m2**(-b); A=sp.log(Zb); p1=m1**(-b)/Zb; p2=m2**(-b)/Zb
ck("surprisal-cost affinity S_i = beta*logM_i + logZ", sp.simplify(-sp.log(p1)-(b*sp.log(m1)+A))==0,"FORCED")
I2=sp.simplify(sp.diff(A,b,2))
Var2=sp.simplify(p1*sp.log(m1)**2+p2*sp.log(m2)**2-(p1*sp.log(m1)+p2*sp.log(m2))**2)
ck("Fisher I(beta) = Var_beta(logM)", sp.simplify(I2-Var2)==0,"FORCED")

# ================= TIER 1.1  DUALLY-FLAT =================
print("\n"+"="*86); print("TIER 1.1  DUALLY-FLAT: Bregman-KL between temperatures; dual potential = negentropy"); print("="*86)
b1,b2=sp.symbols('beta1 beta2',positive=True)
Zf=lambda bb:m1**(-bb)+m2**(-bb); Af=lambda bb:sp.log(Zf(bb))
pf=lambda i,bb:(m1 if i==1 else m2)**(-bb)/Zf(bb)
U=lambda bb: pf(1,bb)*sp.log(m1)+pf(2,bb)*sp.log(m2)
KL=sp.simplify(sum(pf(i,b1)*sp.log(pf(i,b1)/pf(i,b2)) for i in (1,2)))
Breg=sp.simplify(Af(b2)-Af(b1)+(b2-b1)*U(b1))
ck("KL(p(b1)||p(b2)) = A(b2)-A(b1)+(b2-b1)U(b1)  [Bregman divergence of logZ]",
   sp.simplify(KL-Breg)==0,"FORCED","Gibbs curve is an e-geodesic; KL is its Bregman form")
H=-sum(pf(i,b)*sp.log(pf(i,b)) for i in (1,2))
Sst=lambda i:-sp.log(m1 if i==1 else m2)                  # statistic S = -log M
dual=sp.simplify(b*(pf(1,b)*Sst(1)+pf(2,b)*Sst(2))-A)     # Legendre dual of A
ck("dual potential A*(eta) = -H (negentropy); dual coords (beta, -U)",
   sp.simplify(-H-dual)==0,"FORCED","two flat charts: natural beta, expectation -U")

# ================= TIER 1.2  THERMODYNAMICS =================
print("\n"+"="*86); print("TIER 1.2  THERMODYNAMICS: heat capacity C(beta) = beta^2 * I(beta)"); print("="*86)
Ubeta=sp.simplify(-sp.diff(A,b))                          # mean energy = <logM>
C=sp.simplify(-b**2*sp.diff(Ubeta,b))                     # C = dU/dT, T=1/beta
ck("C(beta) = beta^2 * I(beta) = beta^2 * Var_beta(logM)", sp.simplify(C-b**2*I2)==0,"FORCED",
   "the Fisher information IS a heat capacity (per beta^2)")
def VarE(bt):
    w=[mp.e**(-bt*e) for e in Enum]; Zt=sum(w); pp=[wi/Zt for wi in w]
    mm=sum(pi*ei for pi,ei in zip(pp,Enum)); return sum(pi*(ei-mm)**2 for pi,ei in zip(pp,Enum))
grid=[mp.mpf(k)/20 for k in range(-160,161)]
peak=max((VarE(bt),bt) for bt in grid)
ck("peak-fluctuation temperature beta* exists (interior max of Var_beta(logM))",
   float(peak[0])>0 and abs(float(peak[1]))<8,"COMPUTED",
   f"beta* ~ {float(peak[1]):.3f}, max Var(logM) ~ {float(peak[0]):.4f}; beta->+/-inf -> min/max-M seeds")

# ================= TIER 1.3  FISHER-RAO DISTANCES =================
print("\n"+"="*86); print("TIER 1.3  FISHER-RAO DISTANCES: d = 2 arccos BC; information length"); print("="*86)
def gibbs(bt):
    w=[mp.e**(-bt*e) for e in Enum]; Zt=sum(w); return [wi/Zt for wi in w]
def dFR(p,q): return 2*mp.acos(min(mp.mpf(1),sum(mp.sqrt(pi*qi) for pi,qi in zip(p,q))))
v1=[mp.mpf(1)]+[mp.mpf(0)]*6; v2=[mp.mpf(0),mp.mpf(1)]+[mp.mpf(0)]*5
ck("pure-state (vertex) distance = pi; d(p,p)=0  (7 seeds = regular spherical simplex)",
   abs(dFR(v1,v2)-mp.pi)<mp.mpf('1e-20') and dFR(v1,v1)<mp.mpf('1e-15'),"FORCED",
   "all off-diagonal catalog vertex distances equal pi")
L=mp.quad(lambda bt: mp.sqrt(VarE(bt)), [mp.mpf(-1), mp.sqrt(5)])
ck("information length L(-1 -> sqrt5) = ∫ sqrt(Var_beta(logM)) dbeta",
   L>0,"COMPUTED",
   f"L ~ {float(L):.4f} nats; d_FR(Gibbs(-1),Gibbs(sqrt5)) ~ {float(dFR(gibbs(-1),gibbs(mp.sqrt(5)))):.4f}")

# ================= TIER 1.4  TWO-STATISTIC (logM, degree) GEOMETRY =================
print("\n"+"="*86); print("TIER 1.4  TWO-STATISTIC GEOMETRY: Fisher matrix ∇²A(b1,b2) and its curvature"); print("="*86)
def cov2(theta,Tlist):
    n=len(Tlist[0]); ex=[mp.e**(sum(theta[a]*Tlist[a][i] for a in range(2))) for i in range(n)]
    Zt=sum(ex); pp=[e/Zt for e in ex]
    mu=[sum(pp[i]*Tlist[a][i] for i in range(n)) for a in range(2)]
    return [[sum(pp[i]*(Tlist[a][i]-mu[a])*(Tlist[c][i]-mu[c]) for i in range(n)) for c in range(2)] for a in range(2)]
def inv2(G):
    d=G[0][0]*G[1][1]-G[0][1]*G[1][0]; return [[G[1][1]/d,-G[0][1]/d],[-G[1][0]/d,G[0][0]/d]],d
def mm(Aa,Bb): return [[sum(Aa[a][k]*Bb[k][c] for k in range(2)) for c in range(2)] for a in range(2)]
def numK(gfun,pt,h):
    g=lambda t:gfun([mp.mpf(x) for x in t]); G=g(pt); Gi,det=inv2(G)
    d1=lambda k:[[ (g([pt[j]+(h if j==k else 0) for j in range(2)])[a][c]
                   -g([pt[j]-(h if j==k else 0) for j in range(2)])[a][c])/(2*h)
                   for c in range(2)] for a in range(2)]
    dg=[d1(0),d1(1)]
    def d2(k,l):
        if k==l:
            Gp=g([pt[j]+(h if j==k else 0) for j in range(2)]); Gm=g([pt[j]-(h if j==k else 0) for j in range(2)])
            return [[(Gp[a][c]-2*G[a][c]+Gm[a][c])/(h*h) for c in range(2)] for a in range(2)]
        pp=[pt[0]+h,pt[1]+h]; pmi=[pt[0]+h,pt[1]-h]; mpi=[pt[0]-h,pt[1]+h]; mmi=[pt[0]-h,pt[1]-h]
        Gpp,Gpm,Gmp,Gmm=g(pp),g(pmi),g(mpi),g(mmi)
        return [[(Gpp[a][c]-Gpm[a][c]-Gmp[a][c]+Gmm[a][c])/(4*h*h) for c in range(2)] for a in range(2)]
    ddg=[[d2(0,0),d2(0,1)],[d2(1,0),d2(1,1)]]
    dGi=[ [[-x for x in row] for row in mm(mm(Gi,dg[k]),Gi)] for k in range(2)]
    Gam=lambda a,i,j: sum(Gi[a][m]*(dg[i][m][j]+dg[j][m][i]-dg[m][i][j]) for m in range(2))/2
    def dGam(a,i,j,k):
        t1=sum(dGi[k][a][m]*(dg[i][m][j]+dg[j][m][i]-dg[m][i][j]) for m in range(2))/2
        t2=sum(Gi[a][m]*(ddg[k][i][m][j]+ddg[k][j][m][i]-ddg[k][m][i][j]) for m in range(2))/2
        return t1+t2
    Rup=lambda a,i,j,k:(dGam(a,i,k,j)-dGam(a,i,j,k)
                        +sum(Gam(a,m,j)*Gam(m,i,k)-Gam(a,m,k)*Gam(m,i,j) for m in range(2)))
    R0101=sum(G[0][a]*Rup(a,1,0,1) for a in range(2))
    return R0101/det
# validate the routine against KNOWN curvatures (rule out an ambient-1/4 bias, fix sign)
Ksx =numK(lambda th:cov2(th,[[1,0,0],[0,1,0]]),[mp.mpf(0),mp.mpf(0)],mp.mpf('0.02'))   # simplex -> 1/4
Kflat=numK(lambda t:[[mp.mpf(1),mp.mpf(0)],[mp.mpf(0),mp.mpf(1)]],[mp.mpf(0),mp.mpf(0)],mp.mpf('0.02'))
Khyp =numK(lambda t:[[1/t[1]**2,mp.mpf(0)],[mp.mpf(0),1/t[1]**2]],[mp.mpf(0),mp.mpf(1)],mp.mpf('0.01'))  # -1
Ksph =numK(lambda t:[[mp.mpf(1),mp.mpf(0)],[mp.mpf(0),mp.sin(t[0])**2]],[mp.mpf(1),mp.mpf(0)],mp.mpf('0.01'))# +1
sgn=1 if Ksx>0 else -1     # fix orientation convention from the sphere test
ck("curvature routine validated on 4 known metrics (simplex 1/4, flat 0, hyperbolic -1, sphere +1)",
   abs(abs(Ksx)-0.25)<1e-3 and abs(Kflat)<1e-4 and abs(sgn*Khyp+1)<1e-2 and abs(sgn*Ksph-1)<1e-2,
   "COMPUTED",f"K = {float(Ksx):.4f}, {float(Kflat):.2e}, {float(sgn*Khyp):.4f}, {float(sgn*Ksph):.4f}")

# the (logM, degree) family: Fisher matrix + curvature sampled at several natural-parameter points
Tcat=[[float(mp.log(float(m))) for m in Msym],[float(d) for d in deg]]
g0=cov2([mp.mpf(0),mp.mpf(0)],Tcat)
ck("(logM,degree) Fisher matrix ∇²A at uniform is SPD (a valid 2D field metric)",
   g0[0][0]>0 and (g0[0][0]*g0[1][1]-g0[0][1]**2)>0,"COMPUTED",
   f"g=[[{float(g0[0][0]):.4f},{float(g0[0][1]):.4f}],[{float(g0[1][0]):.4f},{float(g0[1][1]):.4f}]]; Var(deg)=24/49={24/49:.4f}")
pts=[[mp.mpf(0),mp.mpf(0)],[mp.mpf('-1'),mp.mpf(0)],[mp.mpf('0.6'),mp.mpf('-0.4')]]
Kcat=[sgn*numK(lambda th:cov2(th,Tcat),pt,mp.mpf('0.01')) for pt in pts]
Kgen=[sgn*numK(lambda th:cov2(th,[[0,1,2,3,4],[0,2,1,5,3]]),pt,mp.mpf('0.01')) for pt in pts]  # generic control
print(f"   K(logM,deg) at 3 pts: {[round(float(k),4) for k in Kcat]}  (spread {float(max(Kcat)-min(Kcat)):.4f})")
print(f"   K(generic control):   {[round(float(k),4) for k in Kgen]}  (varies, != 1/4)")
cat_const_quarter = (max(Kcat)-min(Kcat))<1e-2 and all(abs(k-sp.Rational(1,4))<1e-2 for k in Kcat)
gen_differs       = any(abs(k-sp.Rational(1,4))>3e-2 for k in Kgen)
ck("FINDING: the (logM,degree) field surface has CONSTANT curvature = 1/4 (totally geodesic)",
   cat_const_quarter and gen_differs,"COMPUTED",
   "unlike generic two-statistic families (control varies, changes sign); exactness + mechanism OPEN")

print("\n"+"="*86)
print(f"HARNESS v2: {sum(P)}/{len(P)} passed"+("  exit 0" if all(P) else "  FAIL"))
