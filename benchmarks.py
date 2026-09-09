"""DTLZ benchmarks and Pareto-front reference sets used by the pilot study."""
from __future__ import annotations
import math
import numpy as np


def n_var(problem: str, m: int) -> int:
    p=problem.upper()
    if p=='DTLZ1': return m+4
    if p in ('DTLZ2','DTLZ4'): return m+9
    raise ValueError(problem)


def evaluate(problem: str, x: np.ndarray, m: int) -> np.ndarray:
    p=problem.upper(); x=np.asarray(x,float); n=len(x); k=n-m+1
    if p=='DTLZ1':
        tail=x[m-1:]
        g=100.0*(k + np.sum((tail-0.5)**2 - np.cos(20*math.pi*(tail-0.5))))
        f=np.empty(m)
        for i in range(m):
            v=0.5*(1+g)
            for j in range(m-i-1): v*=x[j]
            if i>0: v*=1-x[m-i-1]
            f[i]=v
        return f
    if p in ('DTLZ2','DTLZ4'):
        alpha=100.0 if p=='DTLZ4' else 1.0
        g=np.sum((x[m-1:]-0.5)**2)
        xx=x[:m-1]**alpha
        f=np.empty(m)
        for i in range(m):
            v=1+g
            for j in range(m-i-1): v*=math.cos(xx[j]*math.pi/2)
            if i>0: v*=math.sin(xx[m-i-1]*math.pi/2)
            f[i]=v
        return f
    raise ValueError(problem)


def ideal(problem: str, m: int):
    return np.zeros(m)


def pareto_reference(problem: str, m: int, n: int=2001) -> np.ndarray:
    p=problem.upper()
    if m==2:
        t=np.linspace(0,1,n)
        if p=='DTLZ1': return np.c_[0.5*t,0.5*(1-t)]
        if p in ('DTLZ2','DTLZ4'): return np.c_[np.cos(t*math.pi/2),np.sin(t*math.pi/2)]
    if m==3:
        # quasi-uniform simplex lattice for DTLZ1; spherical angles for DTLZ2/4
        h=max(10,int(math.sqrt(2*n)))
        pts=[]
        if p=='DTLZ1':
            for i in range(h+1):
                for j in range(h+1-i):
                    k=h-i-j
                    pts.append((0.5*i/h,0.5*j/h,0.5*k/h))
        elif p in ('DTLZ2','DTLZ4'):
            # map a fine simplex lattice through normalized positive vectors
            for i in range(h+1):
                for j in range(h+1-i):
                    k=h-i-j
                    v=np.array([i+0.5,j+0.5,k+0.5],float)
                    v/=np.linalg.norm(v)
                    pts.append(tuple(v))
        else: raise ValueError(problem)
        return np.asarray(pts,float)
    raise ValueError((problem,m))
