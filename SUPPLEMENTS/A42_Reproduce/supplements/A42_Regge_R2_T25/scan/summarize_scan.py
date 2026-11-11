import numpy as np, glob, re
rows=[]
for fn in sorted(glob.glob('L3_k*_*.dat')):
    m=re.match(r'L3_k([0-9.]+)_(flat|rough)\.dat',fn); k=float(m.group(1)); st=m.group(2)
    try: a=np.loadtxt(fn)
    except Exception: continue
    if len(a)<4000: continue
    th=a[a[:,1]==0]; me=a[a[:,1]==1]
    rows.append((k,st,th[:500,4].mean(),me[:1000,4].mean(),me[1000:,4].mean(),me[:,21].mean(),me[:,9].mean(),me[:,5].mean(),(2*me[:,2]-k*me[:,3]).mean()/(15*81), me[:,27].mean()/(50*81)))
print(" k     start  R(th0-500)  R(meas1)  R(meas2)  std_d  maxq    meanq  scaleId  R2/hinge")
for r in rows: print(f"{r[0]:.2f} {r[1]:6s} {r[2]:9.3f} {r[3]:9.4f} {r[4]:9.4f} {r[5]:6.3f} {r[6]:7.2f} {r[7]:7.3f} {r[8]:7.4f} {r[9]:8.3f}")
