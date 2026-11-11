from math import asin,cos,sin
def step(F,V,L,M,dt,dx,s):
 def coin():
  if M:
   u,w=F
   for i,a in enumerate(M):c,q=cos(a),sin(a);u[i],w[i]=c*u[i]-1j*q*w[i],c*w[i]-1j*q*u[i]
 def bonds(p,f):
  for u,v in zip(F,V):
   n=len(u)
   for j in range(p,len(v),2):
    k=(j+1)%n;t=f*asin(dt*v[j]/(2*dx));c,q=cos(t),sin(t);U=L[j] if L else 1
    u[j],u[k]=c*u[j]-q*U.conjugate()*u[k],q*U*u[j]+c*u[k]
 p=s%2;coin();bonds(p,.5);bonds(1-p,1);bonds(p,.5);coin()
