"""IR2-EMOA and close steady-state comparators for 2D/3D continuous MOO."""
from __future__ import annotations
import math
import numpy as np
from integral_r2 import integral_r2_contributions, finite_r2_contributions
from benchmarks import evaluate, ideal, n_var


def dominates(a,b):
    return np.all(a <= b) and np.any(a < b)


def nondominated_sort(F):
    """Fast nondominated sorting for the small steady-state population.

    The dominance matrix is vectorized; this avoids millions of tiny NumPy
    calls from pairwise Python loops during long single-shot runs.
    """
    A=np.asarray(F,float); n=len(A)
    le=np.all(A[:,None,:] <= A[None,:,:], axis=2)
    lt=np.any(A[:,None,:] < A[None,:,:], axis=2)
    dom=le & lt
    cnt=dom.sum(axis=0).astype(int)
    fronts=[]; current=np.flatnonzero(cnt==0).tolist()
    assigned=np.zeros(n,dtype=bool)
    while current:
        fronts.append(current); assigned[current]=True
        dec=dom[current].sum(axis=0).astype(int)
        cnt-=dec
        current=np.flatnonzero((cnt==0) & (~assigned)).tolist()
    rank=np.empty(n,dtype=int)
    for r,fr in enumerate(fronts): rank[fr]=r
    return fronts,rank


def _exclusive_cells_anchored_3d(extents):
    # identical ownership geometry to Integral R2, but ordinary volume measure
    n=len(extents); z_levels=sorted({0.0,*(q[2] for q in extents)})
    groups={}
    for i,q in enumerate(extents): groups.setdefault(q[0],[]).append(i)
    x_groups=sorted(groups.items(),key=lambda kv:kv[0],reverse=True)
    out=[]
    for k in range(len(z_levels)-1):
        z0,z1=z_levels[k],z_levels[k+1]
        if z1<=z0: continue
        top=-math.inf; second=-math.inf; owner=-1; count=0
        active=[]
        for x,inds in x_groups:
            ai=[i for i in inds if extents[i][2] >= z1]
            if ai: active.append((x,ai))
        for g,(x1,inds) in enumerate(active):
            for i in inds:
                y=extents[i][1]
                if y>top: second=top; top=y; owner=i; count=1
                elif y==top: count+=1; owner=-1
                elif y>second: second=y
            x0=active[g+1][0] if g+1<len(active) else 0.0
            y0=0.0 if second==-math.inf else second
            if count==1 and owner>=0 and x1>x0 and top>y0:
                out.append((owner,(x0,x1,y0,top,z0,z1)))
    return out


def hypervolume_contributions(front):
    """Exact HVC after affine front normalization, ref=(1.1,...,1.1).

    Normalization is local to the front to avoid problem-specific arbitrary
    fixed scales.  This remains a dystopian reference-point construction.
    """
    A=np.asarray(front,float); n,m=A.shape
    if n==1: return [math.inf]
    lo=A.min(axis=0); hi=A.max(axis=0); span=np.where(hi>lo,hi-lo,1.0)
    B=(A-lo)/span
    ref=np.full(m,1.1)
    ext=np.maximum(ref-B,0.0)
    out=np.zeros(n)
    if m==2:
        # anchored rectangle exclusive bands
        groups={}
        for i,q in enumerate(ext): groups.setdefault(q[0],[]).append(i)
        groups=sorted(groups.items(),key=lambda kv:kv[0],reverse=True)
        top=-math.inf; second=-math.inf; owner=-1; count=0
        for g,(x1,inds) in enumerate(groups):
            for i in inds:
                y=ext[i,1]
                if y>top: second=top; top=y; owner=i; count=1
                elif y==top: count+=1; owner=-1
                elif y>second: second=y
            x0=groups[g+1][0] if g+1<len(groups) else 0.0
            y0=0.0 if second==-math.inf else second
            if count==1 and owner>=0 and x1>x0 and top>y0:
                out[owner]+=(x1-x0)*(top-y0)
    elif m==3:
        for owner,b in _exclusive_cells_anchored_3d([tuple(q) for q in ext]):
            x0,x1,y0,y1,z0,z1=b
            out[owner]+=(x1-x0)*(y1-y0)*(z1-z0)
    else: raise ValueError("2D/3D only")
    return out.tolist()


def lattice_weights(m, target=101):
    if m==2:
        h=max(2,target-1)
        return [(i/h,1-i/h) for i in range(h+1)]
    if m==3:
        # smallest H yielding roughly target vectors
        h=1
        while (h+2)*(h+1)//2 < target: h+=1
        return [(i/h,j/h,(h-i-j)/h) for i in range(h+1) for j in range(h+1-i)]
    raise ValueError


def sbx(p1,p2,rng,eta=15.0,pc=0.9):
    n=len(p1)
    if rng.random() > pc: return p1.copy()
    c=np.empty(n)
    for i in range(n):
        x1,x2=p1[i],p2[i]
        if rng.random() <= 0.5 and abs(x1-x2)>1e-14:
            yl=0.0; yu=1.0
            if x1>x2: x1,x2=x2,x1
            rand=rng.random()
            beta=1+2*(x1-yl)/(x2-x1); alpha=2-beta**(-(eta+1))
            if rand<=1/alpha: betaq=(rand*alpha)**(1/(eta+1))
            else: betaq=(1/(2-rand*alpha))**(1/(eta+1))
            child=0.5*((x1+x2)-betaq*(x2-x1))
            if rng.random() <= 0.5:
                beta=1+2*(yu-x2)/(x2-x1); alpha=2-beta**(-(eta+1))
                if rand<=1/alpha: betaq=(rand*alpha)**(1/(eta+1))
                else: betaq=(1/(2-rand*alpha))**(1/(eta+1))
                child=0.5*((x1+x2)+betaq*(x2-x1))
            c[i]=min(1.0,max(0.0,child))
        else: c[i]=p1[i]
    return c


def polynomial_mutation(x,rng,eta=20.0,pm=None):
    y=x.copy(); n=len(y); pm=(1/n if pm is None else pm)
    for i in range(n):
        if rng.random() <= pm:
            val=y[i]; delta1=val; delta2=1-val; rand=rng.random(); mut_pow=1/(eta+1)
            if rand<=0.5:
                xy=1-delta1; v=2*rand+(1-2*rand)*(xy**(eta+1)); dq=v**mut_pow-1
            else:
                xy=1-delta2; v=2*(1-rand)+2*(rand-0.5)*(xy**(eta+1)); dq=1-v**mut_pow
            y[i]=min(1.0,max(0.0,val+dq))
    return y


def run(problem='DTLZ2', m=2, algorithm='ir2', pop_size=30, evaluations=1500, seed=1, initial_X=None, n_variables=None):
    rng=np.random.default_rng(seed); nv=(n_var(problem,m) if n_variables is None else int(n_variables)); z=ideal(problem,m)
    X=(rng.random((pop_size,nv)) if initial_X is None else np.asarray(initial_X,float).copy())
    if X.shape != (pop_size,nv):
        raise ValueError(f'initial_X has shape {X.shape}, expected {(pop_size,nv)}')
    F=np.asarray([evaluate(problem,x,m) for x in X])
    weights=lattice_weights(m,101 if m==2 else 120)
    fe=pop_size
    fronts,rank=nondominated_sort(F)
    while fe < evaluations:
        def tournament():
            a,b=rng.integers(0,pop_size,size=2)
            if rank[a] < rank[b]: return a
            if rank[b] < rank[a]: return b
            return a if rng.random()<0.5 else b
        p1,p2=tournament(),tournament()
        child=polynomial_mutation(sbx(X[p1],X[p2],rng),rng)
        fc=evaluate(problem,child,m); fe+=1
        X=np.vstack([X,child]); F=np.vstack([F,fc])
        fronts,rank_aug=nondominated_sort(F); worst=fronts[-1]
        if len(worst)==1:
            kill=worst[0]
        else:
            FW=[tuple(F[i]) for i in worst]
            if algorithm=='ir2': c=integral_r2_contributions(FW,tuple(z))
            elif algorithm=='r2': c=finite_r2_contributions(FW,tuple(z),weights)
            elif algorithm=='sms': c=hypervolume_contributions(FW)
            else: raise ValueError(algorithm)
            # least contribution; deterministic index tie-break
            local=min(range(len(worst)), key=lambda k:(c[k],worst[k]))
            kill=worst[local]
        X=np.delete(X,kill,axis=0); F=np.delete(F,kill,axis=0)
        # Only the last front loses a point, so surviving ranks do not change.
        rank=np.delete(rank_aug,kill)
    fronts,_=nondominated_sort(F)
    nd=fronts[0]
    return X[nd],F[nd]
