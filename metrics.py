from __future__ import annotations
import numpy as np


def _nearest(A,B):
    A=np.asarray(A,float); B=np.asarray(B,float)
    d2=((A[:,None,:]-B[None,:,:])**2).sum(axis=2)
    return np.sqrt(d2.min(axis=1))


def convergence(A, reference):
    return float(_nearest(A,reference).mean())


def delta_p(A, reference, p=2):
    d1=_nearest(A,reference); d2=_nearest(reference,A)
    gd=float((np.mean(d1**p))**(1/p)); igd=float((np.mean(d2**p))**(1/p))
    return max(gd,igd)
