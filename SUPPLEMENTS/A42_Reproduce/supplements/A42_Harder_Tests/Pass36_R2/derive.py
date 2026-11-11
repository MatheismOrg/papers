import sympy as sp
t=sp.symbols('t'); k,beta,kap=sp.symbols('k beta kappa',positive=True)
a=sp.Function('a')(t); N=sp.Function('N')(t)
ad=sp.diff(a,t); add=sp.diff(a,t,2); Nd=sp.diff(N,t)
# Ricci scalar for ds^2=-N^2dt^2+a^2 dSigma_k^2
R=6*(add/(a*N**2)-ad*Nd/(a*N**3)+ad**2/(a**2*N**2)+k/a**2)
L_EH=N*a**3*R/(2*kap); L_R2=N*a**3*beta*R**2/(2*kap)
def EL_N(L):  # Euler-Lagrange wrt N (L depends on N, Nd)
    return sp.diff(L,N)-sp.diff(sp.diff(L,Nd),t)
E_EH=sp.simplify(EL_N(L_EH).subs(N,1).doit())
E_R2=sp.simplify(EL_N(L_R2).subs(N,1).doit())
# matter: L_m=-N a^3 rho  => EL_N = -a^3 rho.  Constraint: E_EH+E_R2 - a^3 rho = 0
H=ad/a; Phi=H**2+k/a**2
Rs=R.subs(N,1).doit(); Rd=sp.diff(Rs,t); Rdd=sp.diff(Rs,t,2)
print('EH check (should be 0):', sp.simplify(E_EH - (-3*a**3*Phi/kap)))
# my effective R^2 density: kappa*rho_R2 = beta*(R^2/2 - 6 H Rdot - 6 R Phi); in constraint form: -(3/kap)a^3 Phi + a^3 rho_R2 ... 
# EH+R2 constraint: E_EH + E_R2 = a^3 rho_m  ->  -3a^3Phi/kap + E_R2 = a^3 rho_m  ->  3Phi = kap rho_m + (kap/a^3) E_R2
rhoR2_from_action=sp.simplify(E_R2/a**3)         # this is rho_R2 (times 1) in units where 3Phi = kap(rho_m + rho_R2)
mine=beta*(Rs**2/2-6*H*Rd-6*Rs*Phi)/kap
print('density check (should be 0):', sp.simplify(rhoR2_from_action-mine))
rho=mine; rhop=-(2*beta*Rdd-2*beta*H*Rd+4*beta*Rs*(sp.diff(H,t)-k/a**2))/kap
print('conservation check (should be 0):', sp.simplify(sp.diff(rho,t)+3*H*rhop))
