import sys; sys.path.insert(0, sys.argv[1])
import numpy as np, gapcheck as G
def dense_only(e, n=54):
    br=e.breaks(); return max(e.g(np.linspace(a,b,n)).max() for a,b in zip(br[:-1],br[1:]))
def lm_polish(e, n=17):
    """split, n uniform per piece, polish EVERY local max of the samples."""
    br,best=e.breaks(),0.0
    for a,b in zip(br[:-1],br[1:]):
        t=np.linspace(a,b,n); v=e.g(t); best=max(best,v.max()); h=t[1]-t[0]
        for i in range(n):
            if (i==0 or v[i]>=v[i-1]) and (i==n-1 or v[i]>=v[i+1]): best=max(best,G.polish(e,a,b,t[i],h))
    return best
for name,kw in [("bspline20",dict(kind="bspline",fit_knots=20)),("cyl30",dict(kind="cyl",fit_knots=30))]:
    tv0=G.truth(G.Edge(**kw))
    for wd in (None,5e-4,1e-4,3e-5):
        e=G.Edge(**(kw if wd is None else dict(kw,bump=(3*tv0,0.4137,wd)))); tv=G.truth(e); out=[]
        for mn,m in [("dense54",dense_only),("lmpol17",lm_polish),("lmpol33",lambda e: lm_polish(e,33)),("occt",G.m_occt_like),("cheb1e-5",lambda e: G.m_cheb(e,1e-5))]:
            e.calls=0; out.append(f"{mn}={m(e)/tv:.3f}({e.calls})")
        print(f"{name} bump={wd} true={tv:.3e} | "+"  ".join(out))
