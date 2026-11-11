"""SA73 T49 exploration (all data): baryon closure as the cluster halo itself. For each CLASH cluster solve for M200c such
that the NFW total mass inside the outermost radius r_o equals M_bar(<r_o) / (f_b Y), Y = baryon depletion there; then
(a) closure alone, (b) geometric combination with the formula's expected halo M_E (weights 1/sigma^2, sigma_E = 0.25 dex,
sigma_closure = 0.15 dex assumed). Also checks T48's sum form with baryons extrapolated toward r200 (aperture test)."""
import numpy as np
from scipy.optimize import brentq
import harness_T45 as M
H=M.H; S=M.S; K=H.K; B=M.B; B.FROZEN['abc']=np.array([11.7481,.8178,-.7119]); th=np.array([11.339])
a=H.arr([r for r in H.cluster_rows() if r['kind']!='MaNGA_BCG']); names=np.unique(a['name'])
mb=a['gb']/H.ACC*a['r']**2/H.G
def nfw_in(m,z,r):
    c=K.c_dm14(np.array([m]),np.array([z]))[0]; rv=(H.G*m/(100*H.hz(np.array([z]))[0]**2))**(1/3); return m*H.ff(r*c/rv)/H.ff(c), rv
def gp(m,k):
    z=a['z'][k]; c=K.c_dm14(m,z); rv=(H.G*m/(100*H.hz(z)**2))**(1/3)
    return a['gb'][k]+H.G*(1-H.FB)*m*H.ff(a['r'][k]*c/rv)/H.ff(c)/a['r'][k]**2*H.ACC
def score(mh,lab):
    out={}
    for kind in ('CLASH_centre','CLASH_outer'):
        ms=[];mn=[]
        for n in names:
            k=(a['name']==n)&(a['kind']==kind); er=np.log10(a['go'][k]/gp(np.full(k.sum(),mh[n]),k)); ms.append(np.mean(er**2)); mn.append(er.mean())
        out[kind]=(round(float(np.sqrt(np.mean(ms))),4),round(float(np.mean(mn)),3))
    print(lab,out,'median logM',round(float(np.median(np.log10(list(mh.values())))),2))
ME={}; RO={}; MBO={}; SL={}
for n in names:
    k=a['name']==n; i=np.argmax(a['r'][k]); RO[n]=a['r'][k][i]; MBO[n]=mb[k][i]
    ME[n]=float(S.m_exp(a['ms'][k][:1],a['z'][k][:1],np.array([5e14]))[0])
    o=np.argsort(a['r'][k]); SL[n]=np.log(mb[k][o][-1]/mb[k][o][-2])/np.log(a['r'][k][o][-1]/a['r'][k][o][-2])
score(ME,'M_E alone')
print('outer baryon slope dlnM/dlnr median',np.median(list(SL.values())).round(2), 'r_o median',np.median(list(RO.values())))
for Y in (1.0,0.9,0.8,0.7):
    MC={}
    for n in names:
        z=a['z'][a['name']==n][0]; MC[n]=brentq(lambda lm: np.log(nfw_in(10**lm,z,RO[n])[0]*H.FB*Y/MBO[n]),13,16.5); MC[n]=10**MC[n]
    score(MC,f'closure alone Y={Y}')
    score({n:10**((np.log10(ME[n])/.25**2+np.log10(MC[n])/.15**2)/(1/.25**2+1/.15**2)) for n in names},f'  combined with M_E, Y={Y}')
# aperture test of T48's sum form: baryons extrapolated to x * r_o with each cluster's own slope
for x in (1.,1.5,2.,3.):
    score({n:ME[n]+MBO[n]*x**SL[n]/H.FB for n in names},f'T48 sum form, baryons to {x} r_o')
