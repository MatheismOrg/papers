"""Check that the compact rule A42_Step_Rule.py reproduces the supplied Section 8 steps exactly,
and count its length as A38 counts the Standard Model statement (raw LZMA2, preset 9 extreme)."""
import sys, json, lzma, importlib.util, numpy as np
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from A42_Step_Rule import step
import os; SUP=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','A42_Supplied_Physics')+'/'
def mod(path,name):
    s=importlib.util.spec_from_file_location(name,SUP+path); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
mc=mod('massive/massive_curved.py','mc'); hz=mod('horizon/horizon_test.py','hz'); gr=mod('gauge/gauge_ring.py','gr')
out={}
# 1. Massive Dirac on a curved static 1+1 background (Section 8.2), periodic ring
for (n,amp,m,c0) in [(128,.25,2.,.9),(256,.25,.5,.9),(128,0.,2.,.9)]:
    ref,rec=mc.local_circuit(n,amp,m,c0)
    dx,x,N,Nmid,y,Y=mc.geometry(n,amp); steps=rec['steps']; dt=rec['dt']
    psi=mc.initial(n,amp,c0); F=[list(psi[:,0]),list(psi[:,1])]
    V=[list(Nmid),list(-Nmid)]; M=list(.5*dt*m*N)
    for s in range(steps): step(F,V,None,M,dt,dx,s)
    out[f'massive_n{n}_a{amp}_m{m}']=float(np.linalg.norm(np.array(F).T-ref))
# 2. Transport across a Schwarzschild horizon, Painleve-Gullstrand (Section 8.4), open interval
for (n,s_,c0) in [(256,-1,3.),(512,1,1.5)]:
    a=(hz.HIGH-hz.LOW)/n; r=hz.LOW+(np.arange(n)+.5)*a; mid=r[:-1]+a/2
    steps=int(np.ceil(hz.T/(.4*a))); steps+=steps%2; dt=hz.T/steps
    psi=np.sqrt(a)*hz.initial(r,c0); psi/=np.linalg.norm(psi)
    th=np.arcsin(dt*hz.vel(mid,s_)/(2*a)); ref=psi.copy()
    for t in range(steps):
        p=t%2
        for par,frac in [(p,.5),(1-p,1.),(p,.5)]:
            j=np.arange(par,n-1,2); c=np.cos(frac*th[j]); b=np.sin(frac*th[j]); u=ref[j].copy(); v=ref[j+1].copy(); ref[j]=c*u-b*v; ref[j+1]=b*u+c*v
    F=[list(psi)]; V=[list(hz.vel(mid,s_))]
    for t in range(steps): step(F,V,None,None,dt,a,t)
    out[f'horizon_n{n}_s{s_}']=float(np.linalg.norm(np.array(F[0])-ref))
# 3. Prescribed electromagnetic link phases (Section 8.2 gauge test), periodic ring
n=128; dx=2*np.pi/n; x=np.arange(n)*dx; dt=dx*.8
links=np.exp(1j*(1.1/n+.2*(np.sin(x+dx)-np.sin(x)))); theta=np.arcsin(dt*np.ones(n)/(2*dx))
psi=np.exp(1j*3*x)*np.exp(-(x-2)**2); psi=psi/np.linalg.norm(psi); ref=psi.copy()
for t in range(40): ref=gr.step(ref,theta,links,t%2)
F=[list(psi)]
for t in range(40): step(F,[list(np.ones(n))],list(links),None,dt,dx,t)
out['gauge_links_n128']=float(np.linalg.norm(np.array(F[0])-ref))
src=open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'A42_Step_Rule.py'),'rb').read()
raw=lzma.compress(src,format=lzma.FORMAT_RAW,filters=[{'id':lzma.FILTER_LZMA2,'preset':9|lzma.PRESET_EXTREME}])
out['rule_source_bytes']=len(src); out['rule_source_bits']=8*len(src)
out['rule_raw_lzma2_bytes']=len(raw); out['rule_raw_lzma2_bits']=8*len(raw)
import zlib, bz2; out['zlib9_bits']=8*len(zlib.compress(src,9)); out['bz2_bits']=8*len(bz2.compress(src,9))
print(json.dumps(out,indent=1)); json.dump(out,open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'A42_Step_Rule_Check.json'),'w'),indent=1)
