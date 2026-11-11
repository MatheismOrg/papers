"""SA73 T47 exploration (all data, no protocol): baryon closure for clusters. In the formula's first term M_p/f_b, M_p is
'the baryons seen'. For CLASH, take M_p = all baryons inside the outermost measured radius (g_bar r^2/G there, stars + hot
gas) instead of the BCG's stars alone. Nothing fitted."""
import numpy as np, json
import harness_T45 as M
H=M.H; B=M.B; B.FROZEN['abc']=np.array([11.7481,.8178,-.7119]); th=np.array([11.339])
cl=[r for r in H.cluster_rows() if r['kind']!='MaNGA_BCG']; a=H.arr(cl)
def score(a,lab):
    er=np.log10(a['go']/H.gpred(M.MODEL,th,a)); out={}
    for k in ('CLASH_centre','CLASH_outer'):
        s=a['kind']==k; hs=np.unique(a['name'][s])
        out[k]=(round(float(np.sqrt(np.mean([np.mean(er[s&(a['name']==h)]**2) for h in hs]))),4), round(float(np.mean(er[s])),4))
    print(lab,out); return er
e0=score(a,'T45')
for rmaxcap in (None,):
    b={k:np.array(v,copy=True) for k,v in a.items()}; mb=a['gb']/H.ACC*a['r']**2/H.G
    for n in np.unique(a['name']):
        k=a['name']==n; i=np.argmax(a['r'][k]); b['mp'][k]=max(mb[k][i], a['ms'][k][0])
    e1=score(b,'M_p = baryons inside outermost radius')
    print('median log M_p/M*', np.median(np.log10(b['mp']/b['ms'])).round(2), 'median log(M_p/f_b)', np.median(np.log10(b['mp']/H.FB)).round(2))
