import json, numpy as np
import harness_T48 as M8
import compare46 as C
H=C.H; S=C.S; K=H.K
cl=[r for r in H.cluster_rows() if r['kind']!='MaNGA_BCG']; a=H.arr(cl)
def g_lcdm_gas(a):
    m=S.m_exp(a['ms'],a['z'],S.thresholds(a))+a['mp']/H.FB; c=K.c_dm14(m,a['z']); rv=(H.G*m/(100*H.hz(a['z'])**2))**(1/3)
    return a['gb']+H.G*(1-H.FB)*m*H.ff(a['r']*c/rv)/H.ff(c)/a['r']**2*H.ACC
recs=json.load(open('res48_MOGF0IESB_primary.json'))['records']; rng=np.random.default_rng(20261109)
import harness_T45 as M5
th=np.array([11.339]); M5.B.FROZEN['abc']=np.array([11.7481,.8178,-.7119])
mh=M5.P.G.mass(th,a['ms'],a['mp'],a['re'],a['z'],S.thresholds(a)); print('formula median log halo', np.median(np.log10(mh)).round(3))
import csv
U={r['host']:float(r['M200c_1e14_Msun_h70inv'])*1e14/(0.671/0.7) for r in csv.DictReader(open(H.PKG/'data/Clusters/umetsu2016_table2.csv'))}
d=[np.log10(U[n])-np.log10(mh[a['name']==n][0]) for n in np.unique(a['name']) if n in U]; print('median Umetsu - formula', np.median(d).round(3), len(d))
er0=np.log10(a['go']/H.gpred('MOGF0IESB',th,a)); print('mean outer offset', er0[a['kind']=='CLASH_outer'].mean().round(3))
for p in ('CLASH_centre','CLASH_outer'):
    pooled={}
    for r in recs:
        for h,v in r['losses'][p].items(): pooled.setdefault(h,[]).append(v)
    hosts=sorted(pooled); f=np.array([np.mean(pooled[h]) for h in hosts]); s=a['kind']==p
    for nm,gf in (('MOND',C.g_mond),('LCDM',C.g_lcdm),('LCDM+gas',g_lcdm_gas)):
        er=np.log10(a['go']/gf(a)); o=np.array([np.mean(er[s&(a['name']==h)]**2) for h in hosts])
        idx=rng.integers(0,len(hosts),(4000,len(hosts))); dd=np.sqrt(f[idx].mean(1))-np.sqrt(o[idx].mean(1))
        print(p,nm,round(float(np.sqrt(f.mean())),4),round(float(np.sqrt(o.mean())),4),round(float(np.sqrt(f.mean())-np.sqrt(o.mean())),4),np.round(np.percentile(dd,[2.5,97.5]),4))
