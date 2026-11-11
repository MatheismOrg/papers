"""Closed FLRW NPQT + free conformal 4D matter; total Box R term set to zero.

All evolution variables use alpha=1 length/time units. lam=eta/alpha.
rad is the independent radiation/state contribution to Friedmann's RHS.
This homogeneous effective model is not a fundamental spacetime rule.
"""
import math
import numpy as np
from scipy.optimize import brentq

def F(y,lam):
    return y/np.sqrt(1-y*y)-lam*y*y

def Fp(y,lam):
    return (1-y*y)**(-1.5)-2*lam*y

def critical_dust(lam):
    # Turning points satisfy c=x^3 F(1/x^2); d/dx=0 determines y.
    y=brentq(lambda y:1-3*y*y+lam*y*(1-y*y)**1.5,.5,1-1e-12)
    x=1/math.sqrt(y)
    return x**3*F(y,lam),x

def turning_points(c,lam,rad=0.):
    from scipy.optimize import minimize_scalar
    def ct(x):return x**3*F(1/x**2,lam)-rad/x
    opt=minimize_scalar(ct,bounds=(1+1e-9,5),method='bounded',options={'xatol':1e-13})
    xm=opt.x
    if c<=opt.fun:raise ValueError('no distinct turning-point pair')
    lo=brentq(lambda x:ct(x)-c,1+1e-12,xm,xtol=1e-13)
    hi=brentq(lambda x:ct(x)-c,xm,max(10,3*c+rad+10),xtol=1e-13)
    return lo,hi

def y_from_a(a,c,lam,rad=0.):
    density=c/a**3+rad/a**4
    return brentq(lambda y:F(y,lam)-density,0,1-1e-13,xtol=1e-14)

def rhs(t,state,c,lam,rad=0.):
    a,v=state;y=(1+v*v)/a**2
    if not 0<=y<1:raise ValueError('outside real action branch')
    fp=Fp(y,lam)
    if fp<=0:raise ValueError('algebraic response branch lost')
    D=y-(3*c/a**3+4*rad/a**4)/(2*fp)
    return np.array([v,a*D])

def quantum_stress(y,D,kappa,state_density=0.):
    rho=kappa*y*y+state_density
    p=kappa*(y*y-4*y*D)/3+state_density/3
    return rho,p

def potential(a,c,lam,rad=0.):
    return a*a*y_from_a(a,c,lam,rad)-1

def discrete_step(a,v,dt,c,lam,rad=0.):
    # Ordinary implicit-midpoint time step, shared between all test cases.
    # No imposed bounce or velocity-sign reversal; nonlinear solve is part
    # of this numerical step, not an asserted microscopic law.
    from scipy.optimize import root
    def eq(z):return z-np.array([a,v])-dt*rhs(0,(z+np.array([a,v]))/2,c,lam,rad)
    start=np.array([a,v])+dt*rhs(0,[a,v],c,lam,rad)
    sol=root(eq,start,tol=1e-11)
    if np.max(abs(eq(sol.x)))>1e-10:raise RuntimeError('implicit step did not solve')
    return sol.x
