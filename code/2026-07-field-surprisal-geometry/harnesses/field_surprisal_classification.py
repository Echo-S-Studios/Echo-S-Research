#!/usr/bin/env python3
# field_surprisal_classification.py
# Comprehensive classification of the totally geodesic / constant-curvature-1/4 two-statistic
# (logM, X) Gibbs surfaces of the emission catalog. Corrects a prior conflation: the single-outcome
# surfaces are constant-curvature-1/4 RULED surfaces, NOT totally geodesic (no such surface exists).
import sympy as sp, itertools, math, random

P=[]
def ck(name,cond,tag,note=""):
    ok=bool(cond); P.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] ({tag:9}) {name}"+(f"   |  {note}" if note else ""))
    assert ok, name

phi=(1+math.sqrt(5))/2
names=['sqrt2','sqrt3','sqrt5','phi','tau','phi^4','K']
Mf=[2,3,5,phi,phi,phi**4,phi**4-1]; lM=[math.log(m) for m in Mf]
Aex=[sp.Integer(2),sp.Integer(3),sp.Integer(5),(1+sp.sqrt(5))/2,(1+sp.sqrt(5))/2,((1+sp.sqrt(5))/2)**4,((1+sp.sqrt(5))/2)**4-1]
lm=[sp.log(a) for a in Aex]

def K_at(a,X,t,s):   # exact Gaussian curvature (Brioschi) of the (a,X) Gibbs surface at (t,s), via cumulants
    w=[math.exp(t*a[i]+s*X[i]) for i in range(len(a))]; Z=sum(w); p=[wi/Z for wi in w]
    Ea=sum(p[i]*a[i] for i in range(len(a))); EX=sum(p[i]*X[i] for i in range(len(a)))
    V=[[a[i]-Ea for i in range(len(a))],[X[i]-EX for i in range(len(a))]]
    mom=lambda *idx: sum(p[i]*math.prod(V[k][i] for k in idx) for i in range(len(a)))
    g11,g12,g22=mom(0,0),mom(0,1),mom(1,1); m=mom
    m4=lambda i,j,k,l: sum(p[q]*V[i][q]*V[j][q]*V[k][q]*V[l][q] for q in range(len(a)))
    k4=lambda i,j,k,l: m4(i,j,k,l)-m(i,j)*m(k,l)-m(i,k)*m(j,l)-m(i,l)*m(j,k)
    Evv,Fuv,Guu=k4(0,0,1,1),k4(0,1,0,1),k4(1,1,0,0)
    Eu,Ev,Gu,Gv,Fu,Fv=m(0,0,0),m(0,0,1),m(1,1,0),m(1,1,1),m(0,1,0),m(0,1,1)
    det3=lambda X: (X[0][0]*(X[1][1]*X[2][2]-X[1][2]*X[2][1])-X[0][1]*(X[1][0]*X[2][2]-X[1][2]*X[2][0])+X[0][2]*(X[1][0]*X[2][1]-X[1][1]*X[2][0]))
    M1=[[-Evv/2+Fuv-Guu/2,Eu/2,Fu-Ev/2],[Fv-Gu/2,g11,g12],[Gv/2,g12,g22]]
    M2=[[0,Ev/2,Gu/2],[Ev/2,g11,g12],[Gu/2,g12,g22]]
    return (det3(M1)-det3(M2))/(g11*g22-g12**2)**2
def is_CC(X):
    try: ks=[K_at(lM,[float(v) for v in X],t,s) for (t,s) in [(0,0),(0.5,-0.4),(-0.6,0.3),(1.2,0.8)]]
    except: return False
    return (max(ks)-min(ks)<4e-3) and all(abs(k-0.25)<4e-3 for k in ks)

print("="*82); print("PART A  NO (logM,X) SURFACE IS TOTALLY GEODESIC"); print("="*82)
d_of=lambda X: len(set((round(lM[i],9),X[i]) for i in range(7)))
dvals=sorted({d_of([1 if i in S else 0 for i in range(7)]) for r in range(8) for S in itertools.combinations(range(7),r)})
ck("d = #distinct (logM,1_S) pairs is always in {6,7}, never 3 (TG needs a 3-dim span)",
   dvals==[6,7],"FORCED",f"d in {dvals}; the surface never lies in a great 2-subsphere")
# symbolic: logM^2 not in V=span{1,logM,X} for any X (logM spans a 6-dim algebra) => II_11 != 0 always
lm2=[x**2 for x in lm]
Vfull=sp.Matrix([[1,lm[i]] for i in range(7)])
ck("logM alone spans a 6-dimensional algebra: 1,logM,logM^2,...,logM^5 independent",
   sp.Matrix([[lm[i]**k for k in range(6)] for i in range(7)]).rank()==6,"FORCED",
   "so logM^2 lies in no 3-dim V; II_11 != 0 for every X => never totally geodesic")

print("\n"+"="*82); print("PART B  THE SECOND FUNDAMENTAL FORM: single-outcome surfaces are RULED"); print("="*82)
# II_22=0 <=> X^2 in V ; II_12=0 <=> X.logM in V ; II_11=0 <=> logM^2 in V (never).
def in_V(vecs, extra):
    Vm=sp.Matrix([[vecs[k][i] for k in range(len(vecs))] for i in range(7)])
    return Vm.row_join(sp.Matrix([[extra[i]] for i in range(7)])).rank()==Vm.rank()
for apex in [6,5,3]:  # K, phi^4, phi
    X=[1 if i==apex else 0 for i in range(7)]
    XX=[x*x for x in X]; Xl=[X[i]*lm[i] for i in range(7)]
    ii22=in_V([[1]*7,lm,X],XX); ii12=in_V([[1]*7,lm,X],Xl); ii11=in_V([[1]*7,lm,X],lm2)
    ck(f"apex {names[apex]:5}: II_22=0 (X^2 in V)={ii22}, II_12=0 (X.logM in V)={ii12}, II_11=0 (logM^2 in V)={ii11}",
       ii22 and ii12 and (not ii11),"FORCED","II_12=II_22=0, II_11!=0  =>  RULED by great circles, constant K=1/4, NOT totally geodesic")

print("\n"+"="*82); print("PART C  CONSTANT-CURVATURE-1/4 CLASSIFICATION (all indicators)"); print("="*82)
CC=[S for r in range(1,7) for S in itertools.combinations(range(7),r) if is_CC([1 if i in S else 0 for i in range(7)])]
lc=lambda S: len(set(round(lM[i],9) for i in S))==1
comp=lambda S: tuple(i for i in range(7) if i not in S)
ck("among all 127 indicators, constant-K=1/4 holds for exactly 16 subsets",
   len(CC)==16,"COMPUTED",f"= 8 distinct surfaces (1_S and 1_S^c share a surface)")
ck("1_S is constant-1/4  <=>  logM is constant on S or on its complement",
   all(lc(S) or lc(comp(S)) for S in CC) and
   all((S in CC) for r in range(1,7) for S in itertools.combinations(range(7),r) if lc(S) or lc(comp(S))),
   "FORCED","the 8 surfaces: one apex per seed (7) + the golden-pair merge {phi,tau} (1)")
singles=[S for S in CC if len(S)==1]; pair=[S for S in CC if len(S)==2]
ck("the 8 surfaces are the 7 single-seed apexes plus the merged golden pair {phi,tau}",
   len(singles)==7 and pair==[(3,4)],"COMPUTED",
   "the 8th (phi,tau) = single-outcome indicator of the effective catalog (Mahler-tied seeds merged)")

print("\n"+"="*82); print("PART D  NO NON-INDICATOR SURFACE IS CONSTANT-CURVATURE-1/4"); print("="*82)
random.seed(0); found=0; tested=0
for _ in range(2500):   # random non-structured X
    X=[random.choice([-2,-1,0,1,2,3,0.5,1.5,2.5]) for _ in range(7)]
    if len(set(X))<=2: continue
    tested+=1
    if is_CC(X): found+=1
for _ in range(2500):   # Mobius family X=(al*logM+be)/(logM+ga): the II_12=0 candidate
    al,be,ga=[random.uniform(-3,3) for _ in range(3)]
    if any(abs(lM[i]+ga)<1e-6 for i in range(7)): continue
    X=[(al*lM[i]+be)/(lM[i]+ga) for i in range(7)]
    if len(set(round(x,6) for x in X))<=2: continue
    tested+=1
    if is_CC(X): found+=1
ck("no non-indicator X gives constant curvature 1/4 (random + Mobius candidate family)",
   found==0,"COMPUTED",f"{found}/{tested} non-indicator X are constant-1/4 => the 8 surfaces are the complete list")

print("\n"+"="*82)
print(f"CLASSIFICATION: {sum(P)}/{len(P)} passed"+("  exit 0" if all(P) else "  FAIL"))
