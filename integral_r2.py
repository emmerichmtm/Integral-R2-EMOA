"""Exact Integral R2 deletion losses for 2D and 3D.

Minimization convention. Objective vectors are translated by an ideal point z*.
Zero loss is allowed and is mapped to an actual +infinity reciprocal coordinate.
No epsilon or finite surrogate is used.

The deletion loss of point i is the Jacobian-weighted exclusive volume of the
reciprocal anchored box B_i.  In 3D the Jacobian density is
    (x+y+z)^(-4)
and the mixed antiderivative is -1/(6(x+y+z)).
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math
from typing import List, Sequence, Tuple

Point = Tuple[float, ...]
Rect2 = Tuple[float,float,float,float]
Box3 = Tuple[float,float,float,float,float,float]

@dataclass(frozen=True)
class ExclusiveCell2D:
    owner: int
    rect: Rect2

@dataclass(frozen=True)
class ExclusiveCell3D:
    owner: int
    box: Box3


def _loss_vectors(points: Sequence[Sequence[float]], ideal: Sequence[float], m: int) -> List[Point]:
    if len(points) == 0:
        raise ValueError("at least one point is required")
    if len(ideal) != m:
        raise ValueError("ideal has wrong dimension")
    z = tuple(float(v) for v in ideal)
    if not all(math.isfinite(v) for v in z):
        raise ValueError("ideal must be finite")
    out: List[Point] = []
    for i,a in enumerate(points):
        if len(a) != m:
            raise ValueError(f"point {i} has wrong dimension")
        aa = tuple(float(v) for v in a)
        if not all(math.isfinite(v) for v in aa):
            raise ValueError("points must be finite")
        p = tuple(aa[j]-z[j] for j in range(m))
        # tolerate only tiny negative floating error at the ideal boundary
        p = tuple(0.0 if -1e-14 <= v < 0.0 else v for v in p)
        if any(v < 0.0 for v in p):
            raise ValueError(f"point {i} is better than supplied ideal: loss={p}")
        out.append(p)
    return out


def _reciprocal(losses: Sequence[Point]) -> List[Point]:
    return [tuple(math.inf if v == 0.0 else 1.0/v for v in p) for p in losses]


def reciprocal_points(points: Sequence[Sequence[float]], ideal: Sequence[float]) -> List[Point]:
    m = len(ideal)
    return _reciprocal(_loss_vectors(points, ideal, m))


def _F2(x: float, y: float) -> float:
    if math.isinf(x) or math.isinf(y):
        return 0.0
    s = x+y
    if s == 0.0:
        raise ValueError("nondegenerate weighted rectangle contains singular origin")
    return 1.0/(2.0*s)


def weighted_rect_integral_2d(rect: Rect2) -> float:
    x0,x1,y0,y1 = rect
    if x0 == x1 or y0 == y1:
        return 0.0
    if x1 < x0 or y1 < y0:
        raise ValueError("invalid rectangle")
    xs=(x0,x1); ys=(y0,y1)
    s=0.0
    for ex,ey in product((0,1), repeat=2):
        sign = -1.0 if (2-(ex+ey)) % 2 else 1.0
        s += sign*_F2(xs[ex],ys[ey])
    return 0.0 if -1e-14 < s < 0.0 else s


def _F3(x: float,y: float,z: float) -> float:
    if math.isinf(x) or math.isinf(y) or math.isinf(z):
        return 0.0
    s=x+y+z
    if s == 0.0:
        raise ValueError("nondegenerate weighted box contains singular origin")
    return -1.0/(6.0*s)


def weighted_box_integral_3d(box: Box3) -> float:
    x0,x1,y0,y1,z0,z1 = box
    # Degeneracy must be handled before any corner evaluation; inf==inf is valid.
    if x0 == x1 or y0 == y1 or z0 == z1:
        return 0.0
    if x1 < x0 or y1 < y0 or z1 < z0:
        raise ValueError("invalid box")
    xs=(x0,x1); ys=(y0,y1); zs=(z0,z1)
    terms=[]
    for ex,ey,ez in product((0,1), repeat=3):
        sign = -1.0 if (3-(ex+ey+ez)) % 2 else 1.0
        terms.append(sign*_F3(xs[ex],ys[ey],zs[ez]))
    s=math.fsum(terms)
    if s < 0.0 and abs(s) < 1e-13*max(1.0, math.fsum(abs(t) for t in terms)):
        return 0.0
    if s < 0.0:
        raise ArithmeticError(f"negative weighted integral {s} for {box}")
    return s


class _TopTwo:
    __slots__=("top","owner","count","second")
    def __init__(self):
        self.top=-math.inf; self.owner=-1; self.count=0; self.second=-math.inf
    def add(self,y:float,owner:int):
        if y > self.top:
            self.second=self.top; self.top=y; self.owner=owner; self.count=1
        elif y == self.top:
            self.count += 1; self.owner=-1
        elif y > self.second:
            self.second=y
    def unique_band(self):
        if self.count != 1 or self.owner < 0:
            return None
        ylo = 0.0 if self.second == -math.inf else self.second
        if self.top <= ylo:
            return None
        return self.owner,ylo,self.top


def _groups_desc(corners: Sequence[Point], axis: int=0):
    d={}
    for i,q in enumerate(corners):
        d.setdefault(q[axis],[]).append(i)
    return sorted(d.items(), key=lambda kv: kv[0], reverse=True)


def exclusive_cells_2d(points: Sequence[Sequence[float]], ideal: Sequence[float]) -> List[ExclusiveCell2D]:
    corners=_reciprocal(_loss_vectors(points,ideal,2))
    groups=_groups_desc(corners,0)
    cells=[]; top2=_TopTwo()
    for g,(x_hi,owners) in enumerate(groups):
        for i in owners:
            top2.add(corners[i][1],i)
        x_lo = groups[g+1][0] if g+1 < len(groups) else 0.0
        band=top2.unique_band()
        if band is not None and x_hi != x_lo:
            owner,ylo,yhi=band
            if yhi != ylo:
                cells.append(ExclusiveCell2D(owner,(x_lo,x_hi,ylo,yhi)))
    return cells


def integral_r2_contributions_2d(points: Sequence[Sequence[float]], ideal: Sequence[float]) -> List[float]:
    n=len(points); _loss_vectors(points,ideal,2)
    if n == 1: return [math.inf]
    out=[0.0]*n
    for c in exclusive_cells_2d(points,ideal):
        out[c.owner]+=weighted_rect_integral_2d(c.rect)
    return [0.0 if abs(v)<1e-15 else v for v in out]


def exclusive_cells_3d(points: Sequence[Sequence[float]], ideal: Sequence[float]) -> List[ExclusiveCell3D]:
    """Exact O(n^2) z-slab/x-sweep exclusive-cell decomposition.

    Supports actual +infinity reciprocal coordinates.  The sorted z levels may
    therefore end in +infinity; the last slab is semi-infinite and the weighted
    corner formula handles it by its limiting value.
    """
    corners=_reciprocal(_loss_vectors(points,ideal,3))
    n=len(corners)
    z_levels=sorted({0.0,*(q[2] for q in corners)})
    x_groups=_groups_desc(corners,0)
    cells=[]
    for k in range(len(z_levels)-1):
        z_lo,z_hi=z_levels[k],z_levels[k+1]
        if z_lo == z_hi: continue
        top2=_TopTwo(); active_groups=[]
        for x,owners in x_groups:
            active=[i for i in owners if corners[i][2] >= z_hi]
            if active: active_groups.append((x,active))
        for g,(x_hi,owners) in enumerate(active_groups):
            for i in owners:
                top2.add(corners[i][1],i)
            x_lo = active_groups[g+1][0] if g+1 < len(active_groups) else 0.0
            band=top2.unique_band()
            if band is not None and x_hi != x_lo:
                owner,ylo,yhi=band
                if yhi != ylo:
                    cells.append(ExclusiveCell3D(owner,(x_lo,x_hi,ylo,yhi,z_lo,z_hi)))
    return cells


def integral_r2_contributions_3d(points: Sequence[Sequence[float]], ideal: Sequence[float]) -> List[float]:
    n=len(points); _loss_vectors(points,ideal,3)
    if n == 1: return [math.inf]
    out=[0.0]*n
    for c in exclusive_cells_3d(points,ideal):
        out[c.owner]+=weighted_box_integral_3d(c.box)
    return [0.0 if abs(v)<1e-15 else v for v in out]


def integral_r2_contributions(points: Sequence[Sequence[float]], ideal: Sequence[float]) -> List[float]:
    m=len(ideal)
    if m==2: return integral_r2_contributions_2d(points,ideal)
    if m==3: return integral_r2_contributions_3d(points,ideal)
    raise ValueError("only 2D and 3D are supported")


def finite_r2_contributions(points: Sequence[Sequence[float]], ideal: Sequence[float], weights: Sequence[Sequence[float]]) -> List[float]:
    """Deletion losses for the finite-weight unary R2 indicator in O(n|W|)."""
    m=len(ideal); losses=_loss_vectors(points,ideal,m); n=len(losses)
    if n==1: return [math.inf]
    out=[0.0]*n
    for w in weights:
        if len(w)!=m: raise ValueError("weight dimension mismatch")
        vals=[max(float(w[j])*losses[i][j] for j in range(m)) for i in range(n)]
        order=sorted(range(n), key=vals.__getitem__)
        i0=order[0]; v0=vals[i0]; v1=vals[order[1]]
        # If tied for best, deleting any tied point changes nothing.
        if v1 > v0:
            out[i0] += v1-v0
    inv=1.0/len(weights)
    return [v*inv for v in out]
