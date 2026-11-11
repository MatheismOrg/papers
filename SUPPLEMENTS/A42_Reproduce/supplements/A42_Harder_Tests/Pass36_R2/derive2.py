import sympy as sp, pickle
t=sp.symbols('t'); k,b=sp.symbols('k b')
a=sp.Function('a')(t); N=sp.Function('N')(t)
ad,add=sp.diff(a,t),sp.diff(a,t,2); Nd=sp.diff(N,t)
R=6*(add/(a*N**2)-ad*Nd/(a*N**3)+ad**2/(a**2*N**2)+k/a**2)
EL=lambda L: sp.diff(L,N)-sp.diff(sp.diff(L,Nd),t)
E_EH=sp.simplify(EL(N*a**3*R/2).subs(N,1).doit())          # kappa=1
E_R2=sp.simplify(EL(N*a**3*b*R**2/2).subs(N,1).doit())
H=ad/a; Phi=H**2+k/a**2
print('E_EH/(3 a^3 Phi) =', sp.simplify(E_EH/(3*a**3*Phi)))
# constraint: E_EH + E_R2 - a^3 rho_m = 0 ; divide by 3 a^3:  Phi + E_R2/(3a^3) = rho_m/3  => Phi = rho_m/3 - E_R2/(3a^3)
Q=sp.simplify(-E_R2/(3*a**3))           # the R^2 term on the matter side, normalized like c/a^3
Qd=sp.diff(Q,t)
P=sp.simplify(sp.cancel(Qd/H))
print('P polynomial in derivatives (no division by adot):', not sp.denom(sp.together(P)).has(ad))
a0,a1,a2,a3,a4=sp.symbols('a0 a1 a2 a3 a4')
sub={sp.diff(a,t,4):a4,sp.diff(a,t,3):a3,sp.diff(a,t,2):a2,sp.diff(a,t):a1}
Qs=sp.simplify(Q.subs(sub).subs(a,a0)); Ps=sp.simplify(P.subs(sub).subs(a,a0))
print('Q =',sp.factor(Qs)); print('P coeff of a4 =',sp.simplify(sp.diff(Ps,a4)))
# GR+R^2 sanity: known Starobinsky flat-space check, de Sitter a=exp(h t) should give Q=0 for k=0
h=sp.symbols('h',positive=True)
print('Q on flat de Sitter (should be 0):', sp.simplify(Qs.subs({k:0,a1:h*a0,a2:h**2*a0,a3:h**3*a0})))
pickle.dump((sp.srepr(Qs),sp.srepr(Ps)),open('QP.pkl','wb'))
