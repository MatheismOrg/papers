"""SA73 T49 exploration (all data, no protocol): two-halo term for the KiDS stacks. Matter correlated with each halo adds
dSigma_2h(R) = b(M) * dSigma_matter(R): Tinker et al. (2010) bias at M200c, linear matter correlation (colossus, the
A19planck cosmology), z = Z0. Nothing fitted. Prints chi2 by stack with a multiplier f on the term (f = 1 is the full term)."""
import numpy as np, json
from scipy.linalg import solve_triangular
from scipy.interpolate import interp1d
from colossus.lss import bias as cb
import harness_T45 as M
H=M.H; S=M.S; K=H.K; B=M.B; P=M.P
B.FROZEN['abc']=np.array([11.7481,.8178,-.7119]); th=np.array([11.339])
cosmo=S.COSMO if hasattr(S,'COSMO') else __import__('colossus.cosmology.cosmology',fromlist=['x']).getCurrent()
h=cosmo.h; z=float(K.Z0); print('z0',z,'h',h)
# dSigma_matter(R) in comoving (Msun/h)/(Mpc/h)^2 -> physical Msun/pc^2
rr=np.logspace(-2.5,2.3,400); xi=cosmo.correlationFunction(rr,z)          # r in Mpc/h comoving
lxi=interp1d(np.log(rr),xi,bounds_error=False,fill_value=(xi[0],0.))
rho_m=cosmo.rho_m(0.)*1e9                                                   # Msun h^2 / Mpc^3 comoving
def sigma(Rc):
    l=np.logspace(-3,2.2,600); r=np.sqrt(Rc**2+l**2); return 2*rho_m*np.trapezoid(lxi(np.log(r))*l,np.log(l))
Rg=np.logspace(-3,1,160); Sg=np.array([sigma(R) for R in Rg])
cum=np.concatenate([[0],np.cumsum(.5*(Sg[1:]*Rg[1:]**2+Sg[:-1]*Rg[:-1]**2)*np.diff(np.log(Rg)))])+.5*Sg[0]*Rg[0]**2
dS=2*cum/Rg**2-Sg                                                           # h Msun/Mpc^2 comoving
def dsig_phys(Rphys_mpc):                                                    # R physical Mpc -> Msun/pc^2 physical
    Rc=Rphys_mpc*(1+z)*h; return np.interp(np.log(Rc),np.log(Rg),dS)*h*(1+z)**2/1e12
gs=H.disks()
for drop in (False,True):
    L=H.Lens(gs,drop)
    for f in (0.,.5,1.):
        tot=[]
        for j,st in enumerate(L.stacks):
            base=H.esd(M.MODEL,th,st)
            m=P.G.mass(th,st.ms,st.mp,st.re,st.z,0.); b=cb.haloBias(m*h,model='tinker10',z=z,mdef='200c')
            two=st.w@(b[:,None]*np.sum(st.aw*dsig_phys(st.R),axis=2))
            e=(base+f*two)[L.masks[j]]; zz=solve_triangular(L.L[j],L.y[j]-e,lower=True); tot.append(round(float(zz@zz),2))
            if f==1. and not drop: print('  stack',j+1,'bias',round(float(st.w@b),2),'2h',np.round(two[L.masks[j]],2),'model',np.round(base[L.masks[j]],2),'data',np.round(L.y[j],2))
        print('strict' if drop else 'primary','f',f,tot,round(sum(tot),2))
