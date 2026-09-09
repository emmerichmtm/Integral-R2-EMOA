"""ZDT/DTLZ benchmarks and Pareto-front references used by IR2-EMOA.

The :func:`ideal` function returns the *exact componentwise ideal point* of the
benchmark, not an artificially strictly dominating point.  Consequently some
attained objective values may have zero loss.  The Integral R2 implementation
maps such zero losses to an actual reciprocal ``+inf``.
"""
from __future__ import annotations
import math
import numpy as np


def _bisect_root(fun, lo: float, hi: float, iterations: int = 90) -> float:
    flo = fun(lo); fhi = fun(hi)
    if flo == 0.0: return lo
    if fhi == 0.0: return hi
    if flo * fhi > 0.0:
        raise ValueError("root is not bracketed")
    for _ in range(iterations):
        mid = 0.5*(lo+hi); fm = fun(mid)
        if flo * fm <= 0.0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5*(lo+hi)


def _zdt3_f2_pf(x: float) -> float:
    return 1.0 - math.sqrt(x) - x*math.sin(10.0*math.pi*x)


def _zdt3_ideal_f2() -> float:
    # Global PF minimum is the stationary point in this final nondominated arc.
    def d(x: float) -> float:
        return (-0.5/math.sqrt(x)
                - math.sin(10.0*math.pi*x)
                - 10.0*math.pi*x*math.cos(10.0*math.pi*x))
    xstar = _bisect_root(d, 0.84, 0.87)
    return _zdt3_f2_pf(xstar)


def _dtlz7_q(x: float) -> float:
    return x*(1.0 + math.sin(3.0*math.pi*x))


def _dtlz7_dq(x: float) -> float:
    return 1.0 + math.sin(3.0*math.pi*x) + 3.0*math.pi*x*math.cos(3.0*math.pi*x)


def _dtlz7_pf_intervals():
    """Return the two 1-D record intervals that generate the DTLZ7 PF.

    For g=1 the last objective is additive in q(x)=x(1+sin(3*pi*x)).
    A coordinate value is Pareto-relevant exactly while q is a new prefix
    maximum.  This gives [0,a] U [b,c], where a and c are the first and
    last local maxima and q(b)=q(a) on the rising branch after x=1/2.
    """
    a = _bisect_root(_dtlz7_dq, 0.20, 0.30)
    c = _bisect_root(_dtlz7_dq, 5.0/6.0, 0.90)
    qa = _dtlz7_q(a)
    b = _bisect_root(lambda x: _dtlz7_q(x)-qa, 0.50, c)
    return a,b,c


ZDT3_IDEAL_F2 = _zdt3_ideal_f2()
DTLZ7_A, DTLZ7_B, DTLZ7_C = _dtlz7_pf_intervals()
DTLZ7_QMAX = _dtlz7_q(DTLZ7_C)


def n_var(problem: str, m: int) -> int:
    p=problem.upper()
    if p in ('ZDT1','ZDT2','ZDT3'):
        if m != 2: raise ValueError(f"{p} is bi-objective")
        return 30
    if p=='DTLZ1': return m+4
    if p in ('DTLZ2','DTLZ4'): return m+9
    if p=='DTLZ7': return m+19  # k=20
    raise ValueError(problem)


def evaluate(problem: str, x: np.ndarray, m: int) -> np.ndarray:
    p=problem.upper(); x=np.asarray(x,float); n=len(x)
    if p in ('ZDT1','ZDT2','ZDT3'):
        if m != 2: raise ValueError(f"{p} is bi-objective")
        f1=x[0]
        g=1.0 + 9.0*np.sum(x[1:])/(n-1)
        ratio=f1/g
        if p=='ZDT1': h=1.0-math.sqrt(ratio)
        elif p=='ZDT2': h=1.0-ratio*ratio
        else: h=1.0-math.sqrt(ratio)-ratio*math.sin(10.0*math.pi*f1)
        return np.array([f1,g*h],float)

    k=n-m+1
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
    if p=='DTLZ7':
        tail=x[m-1:]
        g=1.0 + 9.0*np.sum(tail)/k
        f=np.empty(m)
        f[:m-1]=x[:m-1]
        h=m - np.sum((f[:m-1]/(1.0+g))*(1.0+np.sin(3.0*math.pi*f[:m-1])))
        f[m-1]=(1.0+g)*h
        return f
    raise ValueError(problem)


def ideal(problem: str, m: int):
    """Return the exact componentwise ideal point (objective infima)."""
    p=problem.upper()
    if p in ('ZDT1','ZDT2'):
        if m != 2: raise ValueError(f"{p} is bi-objective")
        return np.array([0.0,0.0])
    if p=='ZDT3':
        if m != 2: raise ValueError("ZDT3 is bi-objective")
        return np.array([0.0,ZDT3_IDEAL_F2])
    if p in ('DTLZ1','DTLZ2','DTLZ4'):
        return np.zeros(m)
    if p=='DTLZ7':
        # On the PF g=1.  The last objective is
        # 2m - sum_i q(f_i), q(x)=x(1+sin(3*pi*x)).
        z=np.zeros(m)
        z[-1]=2.0*m-(m-1)*DTLZ7_QMAX
        return z
    raise ValueError(problem)


def _nondominated_rows(A: np.ndarray, atol: float=0.0) -> np.ndarray:
    """Filter a small dense 2D reference sample to its nondominated rows."""
    A=np.asarray(A,float)
    if A.shape[1] != 2: raise ValueError("2D helper")
    # Sort by f1 ascending, retain points whose f2 is a new strict minimum.
    order=np.lexsort((A[:,1],A[:,0]))
    out=[]; best=math.inf
    for i in order:
        y=A[i,1]
        if y < best-atol:
            out.append(A[i]); best=y
    return np.asarray(out,float)


def pareto_reference(problem: str, m: int, n: int=2001) -> np.ndarray:
    p=problem.upper()
    if p in ('ZDT1','ZDT2','ZDT3'):
        if m != 2: raise ValueError(f"{p} is bi-objective")
        t=np.linspace(0.0,1.0,n)
        if p=='ZDT1': return np.c_[t,1.0-np.sqrt(t)]
        if p=='ZDT2': return np.c_[t,1.0-t*t]
        raw=np.c_[t, 1.0-np.sqrt(t)-t*np.sin(10.0*math.pi*t)]
        return _nondominated_rows(raw)

    if m==2:
        t=np.linspace(0,1,n)
        if p=='DTLZ1': return np.c_[0.5*t,0.5*(1-t)]
        if p in ('DTLZ2','DTLZ4'): return np.c_[np.cos(t*math.pi/2),np.sin(t*math.pi/2)]
        if p=='DTLZ7':
            n1=max(2,n//2)
            n2=max(2,n-n1)
            t=np.r_[np.linspace(0.0,DTLZ7_A,n1,endpoint=True),
                    np.linspace(DTLZ7_B,DTLZ7_C,n2,endpoint=True)]
            return np.c_[t, 4.0-t*(1.0+np.sin(3.0*math.pi*t))]
    if m==3:
        # quasi-uniform simplex lattice for DTLZ1; spherical points for DTLZ2/4
        h=max(10,int(math.sqrt(2*n)))
        pts=[]
        if p=='DTLZ1':
            for i in range(h+1):
                for j in range(h+1-i):
                    k=h-i-j
                    pts.append((0.5*i/h,0.5*j/h,0.5*k/h))
        elif p in ('DTLZ2','DTLZ4'):
            for i in range(h+1):
                for j in range(h+1-i):
                    k=h-i-j
                    v=np.array([i+0.5,j+0.5,k+0.5],float)
                    v/=np.linalg.norm(v)
                    pts.append(tuple(v))
        elif p=='DTLZ7':
            # Four disconnected Pareto patches = Cartesian product of the
            # two Pareto-relevant intervals [0,a] and [b,c].
            per=max(8,int(math.ceil(math.sqrt(max(4,n)/4.0))))
            intervals=((0.0,DTLZ7_A),(DTLZ7_B,DTLZ7_C))
            for lo1,hi1 in intervals:
                for lo2,hi2 in intervals:
                    for f1 in np.linspace(lo1,hi1,per,endpoint=True):
                        for f2 in np.linspace(lo2,hi2,per,endpoint=True):
                            f3 = 6.0 - _dtlz7_q(float(f1)) - _dtlz7_q(float(f2))
                            pts.append((float(f1),float(f2),float(f3)))
        else: raise ValueError((problem,m))
        return np.asarray(pts,float)
    raise ValueError((problem,m))
