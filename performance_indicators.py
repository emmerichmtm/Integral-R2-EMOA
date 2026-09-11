"""Performance indicators used in the IR2-EMOA EMO 2027 study.

Reports exact Integral R2, exact hypervolume (2D/3D), and Delta_p.
Integral R2 is evaluated with the perspective-mapping weighted-complement
formulation. In dimension m the Jacobian weight is (sum x)^(-(m+1)); the
normalized uniform measure contributes the factor (m-1)!.
"""
from __future__ import annotations
import math
import numpy as np
from integral_r2 import (
    reciprocal_points,
    weighted_rect_integral_2d,
    weighted_box_integral_3d,
)


def delta_p(A, R, p=2):
    A=np.asarray(A,float); R=np.asarray(R,float)
    d_ar=np.sqrt(((A[:,None,:]-R[None,:,:])**2).sum(axis=2)).min(axis=1)
    d_ra=np.sqrt(((R[:,None,:]-A[None,:,:])**2).sum(axis=2)).min(axis=1)
    gd=(np.mean(d_ar**p))**(1.0/p)
    igd=(np.mean(d_ra**p))**(1.0/p)
    return float(max(gd,igd))


def hypervolume(A, ref):
    """Exact HV for 2D/3D minimization by coordinate-cell decomposition."""
    A=np.asarray(A,float); ref=np.asarray(ref,float)
    A=A[np.all(A <= ref,axis=1)]
    if len(A)==0:
        return 0.0
    m=A.shape[1]
    if m==2:
        A=A[np.argsort(A[:,0])]
        hv=0.0; y=ref[1]
        for x,yy in A:
            if yy < y:
                hv += max(0.0,ref[0]-x)*(y-yy)
                y=yy
        return float(hv)
    if m==3:
        xs=sorted(set(A[:,0].tolist()+[float(ref[0])]))
        hv=0.0
        for i in range(len(xs)-1):
            x0,x1=xs[i],xs[i+1]
            if x1<=x0: continue
            S=A[A[:,0] <= x0+1e-15][:,1:]
            if len(S): hv += (x1-x0)*hypervolume(S,ref[1:])
        return float(hv)
    raise ValueError('2D/3D only')


def _finite_axis_bounds(values):
    vals=sorted({0.0,*(float(v) for v in values if math.isfinite(float(v)))})
    vals.append(math.inf)
    return vals


def _integral_r2_2d(A, ideal):
    corners=reciprocal_points(A,ideal)
    xb=_finite_axis_bounds(q[0] for q in corners)
    total=0.0
    for i in range(len(xb)-1):
        x0,x1=xb[i],xb[i+1]
        # Coverage is constant in the open x-cell. A reciprocal anchored box
        # spans the whole cell iff its x extent reaches x1.
        ymax=0.0
        for bx,by in corners:
            if bx >= x1 and by > ymax:
                ymax=by
        if math.isinf(ymax):
            continue
        total += weighted_rect_integral_2d((x0,x1,ymax,math.inf))
    return float(total)  # (m-1)! = 1


def _integral_r2_3d(A, ideal):
    corners=reciprocal_points(A,ideal)
    xb=_finite_axis_bounds(q[0] for q in corners)
    yb=_finite_axis_bounds(q[1] for q in corners)
    total=0.0
    # For every x-y cell the union of reciprocal anchored boxes covers z from
    # zero up to one zmax. The complement above zmax is integrated exactly.
    for i in range(len(xb)-1):
        x0,x1=xb[i],xb[i+1]
        for j in range(len(yb)-1):
            y0,y1=yb[j],yb[j+1]
            zmax=0.0
            for bx,by,bz in corners:
                if bx >= x1 and by >= y1 and bz > zmax:
                    zmax=bz
            if math.isinf(zmax):
                continue
            total += weighted_box_integral_3d((x0,x1,y0,y1,zmax,math.inf))
    return float(2.0*total)  # normalized simplex measure: (3-1)! = 2


def ir2_value(A, ideal):
    """Exact Integral R2 via perspective mapping in two or three objectives."""
    A=np.asarray(A,float); ideal=np.asarray(ideal,float)
    if A.ndim != 2 or A.shape[1] != len(ideal):
        raise ValueError('dimension mismatch')
    if A.shape[1]==2:
        return _integral_r2_2d(A,ideal)
    if A.shape[1]==3:
        return _integral_r2_3d(A,ideal)
    raise ValueError('2D/3D only')
