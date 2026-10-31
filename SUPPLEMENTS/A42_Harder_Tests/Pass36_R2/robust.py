import sys, math, numpy as np, json
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import r2_bounce as R
from scipy.integrate import solve_ivp
from pass36_model import turning_points, y_from_a, Fp, F
c=R.C
for b in (2e-3,5e-3):
  for meth,rt in (('Radau',1e-12),('LSODA',1e-11)):
    amin,amax=turning_points(c,R.LAM); a=(amin+amax)/2; y=y_from_a(a,c,R.LAM); a1=-math.sqrt(a*a*y-1); a2=a*(y-3*c/a**3/(2*Fp(y,R.LAM)))
    a3=-(R.Q_no3(a,a1,a2,b)-(F(y,R.LAM)-c/a**3))/R.cQ3(a,a1,b)
    def edge(t,z): return 1-1e-9-(1+z[1]**2)/z[0]**2
    edge.terminal=True
    def back(t,z): return z[0]-a
    back.direction=1; back.terminal=True
    try:
        s=solve_ivp(lambda t,z:R.rhs4(t,z,b),(0,20),[a,a1,a2,a3],method=meth,rtol=rt,atol=1e-14,events=[edge,back])
        print(b,meth,'edge' if len(s.t_events[0]) else ('returned' if len(s.t_events[1]) else 'neither'), round(float(s.y[0].min()),5), round(float(s.t[-1]),4))
    except Exception as e: print(b,meth,'exception',e)
