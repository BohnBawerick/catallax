"""Compare edge-gap checker designs on synthetic edges with a known worst gap.

g(t) = C(t) - S(p(t)). C is a least-squares fit of S(p(t)), so the gap is a real fit residual.
Truth = 400k uniform samples + bounded polish around the best one.
"""
import numpy as np
from numpy.polynomial import chebyshev as Ch
from scipy.interpolate import BSpline, make_lsq_spline
from scipy.optimize import brentq, minimize_scalar

rng = np.random.default_rng(7)
K = 3


def clamped(nint, a=0.0, b=1.0):
    inner = np.linspace(a, b, nint + 2)[1:-1]
    return np.r_[[a] * (K + 1), inner, [b] * (K + 1)]


class Edge:
    def __init__(s, kind, fit_knots, bump=None):
        s.kind, s.calls = kind, 0
        if kind == "bspline":
            s.tu, s.tv = clamped(4), clamped(4)
            n = len(s.tu) - K - 1
            gu, gv = np.meshgrid(np.linspace(0, 100, n), np.linspace(0, 100, n), indexing="ij")
            s.P = np.stack([gu, gv, 15 * np.random.default_rng(7).standard_normal((n, n))], -1)
            pc = np.c_[np.linspace(0.05, 0.95, 8), 0.5 + 0.35 * np.sin(np.linspace(0, 5, 8))]
        else:  # cylinder R=10, u angle, v height
            s.tu = s.tv = np.array([0.0, 1.0])
            pc = np.c_[np.linspace(0, 3.0, 8), 10 + 4 * np.sin(np.linspace(0, 6, 8))]
        s.tp = clamped(len(pc) - K - 1)
        s.p = BSpline(s.tp, pc, K)
        tt = np.linspace(0, 1, 4000)
        w = np.ones_like(tt); w[[0, -1]] = 1e5  # fitter pins the vertices
        s.C = make_lsq_spline(tt, s.surf(s.p(tt)), clamped(fit_knots), K, w=w)
        s.bump = bump

    def surf(s, uv):
        u, v = uv[..., 0], uv[..., 1]
        if s.kind == "cyl":
            return np.stack([10 * np.cos(u), 10 * np.sin(u), v], -1)
        Nu = BSpline.design_matrix(np.clip(u, 0, 1), s.tu, K).toarray()
        Nv = BSpline.design_matrix(np.clip(v, 0, 1), s.tv, K).toarray()
        return np.einsum("ni,nj,ijk->nk", Nu, Nv, s.P)

    def g(s, t):  # gap norm, counts true evaluations
        t = np.atleast_1d(np.asarray(t, float))
        s.calls += t.size
        d = s.C(t) - s.surf(s.p(t))
        if s.bump:
            a, c, w = s.bump
            d[:, 2] += a * np.exp(-(((t - c) / w) ** 2))
        return np.linalg.norm(d, axis=1)

    def breaks(s):
        """Split points: knots of C and p, plus where p crosses surface knot lines."""
        b = set(np.unique(s.C.t)) | set(np.unique(s.tp))
        tt = np.linspace(0, 1, 4001)
        uv = s.p(tt)
        for dim, kn in ((0, s.tu), (1, s.tv)):
            for k in np.unique(kn)[1:-1]:
                f = uv[:, dim] - k
                for i in np.nonzero(np.sign(f[:-1]) != np.sign(f[1:]))[0]:
                    b.add(brentq(lambda x: s.p(x)[dim] - k, tt[i], tt[i + 1]))
        return np.array(sorted(x for x in b if 0 <= x <= 1))


def polish(e, a, b, x0, h):
    lo, hi = max(a, x0 - h), min(b, x0 + h)
    r = minimize_scalar(lambda x: -e.g(x)[0], bounds=(lo, hi), method="bounded",
                        options={"xatol": 1e-10})
    return max(-r.fun, e.g(x0)[0])


def m_uniform(e, n=23):
    return e.g(np.linspace(0, 1, n)).max()


def m_golden(e, n=64):  # low-discrepancy, deterministic
    return e.g((np.arange(1, n + 1) * 0.6180339887498949) % 1.0).max()


def m_split_k(e, k=8):
    br = e.breaks()
    t = np.unique(np.concatenate([np.linspace(a, b, k + 1) for a, b in zip(br[:-1], br[1:])]))
    return e.g(t).max()


def m_split_polish(e, k=8):
    """split+8, then polish the best sample of each piece with bounded Brent."""
    br, best = e.breaks(), 0.0
    for a, b in zip(br[:-1], br[1:]):
        t = np.linspace(a, b, k + 1)
        v = e.g(t)
        best = max(best, polish(e, a, b, t[v.argmax()], t[1] - t[0]))
    return best


def m_occt_like(e, deg=18):
    """Split, 3*degree uniform samples per piece, polish every local max with bounded Brent."""
    br, best = e.breaks(), 0.0
    for a, b in zip(br[:-1], br[1:]):
        t = np.linspace(a, b, 3 * deg)
        v = e.g(t)
        best = max(best, v.max())
        h = t[1] - t[0]
        for i in range(len(t)):
            if (i == 0 or v[i] >= v[i - 1]) and (i == len(t) - 1 or v[i] >= v[i + 1]):
                best = max(best, polish(e, a, b, t[i], h))
    return best


def cheb_pts(n, a, b):
    return 0.5 * (a + b) + 0.5 * (b - a) * np.cos(np.pi * np.arange(n + 1) / n)


def vals2coeffs(v):
    n = len(v) - 1
    ext = np.r_[v, v[-2:0:-1]]
    c = np.real(np.fft.fft(ext))[: n + 1] / n
    c[0] /= 2
    c[-1] /= 2
    return c


def m_cheb(e, rel_chop=1e-3, nmax=129):
    """Split; per piece adaptive Chebyshev interpolation of each component of g on nested grids
    17, 33, 65, 129; stop when the coefficient tail is below rel_chop * max|g| on the piece;
    two off-grid sample tests; max of the interpolant, then polish on the true g."""
    br, best = e.breaks(), 0.0
    for a, b in zip(br[:-1], br[1:]):
        n = 16
        while True:
            t = cheb_pts(n, a, b)
            e.calls += t.size  # vector g evaluation, counted once per point
            d = e.C(t) - e.surf(e.p(t))
            if e.bump:
                A, c0, w = e.bump
                d[:, 2] += A * np.exp(-(((t - c0) / w) ** 2))
            co = np.stack([vals2coeffs(d[:, k]) for k in range(3)], -1)
            scale = np.linalg.norm(d, axis=1).max() + 1e-15
            tail = np.abs(co[-4:]).max()
            x = 0.5 * (a + b) + 0.5 * (b - a) * np.array([0.3183098861837907, -0.7071067811865476])
            xs = (2 * x - a - b) / (b - a)
            interp = np.linalg.norm(np.stack([Ch.chebval(xs, co[:, k]) for k in range(3)], -1), axis=1)
            sample_ok = np.all(np.abs(interp - e.g(x)) <= 10 * rel_chop * scale)
            if (tail <= rel_chop * scale and sample_ok) or n >= nmax - 1:
                break
            n *= 2
        xx = np.linspace(-1, 1, 4001)
        mag = np.linalg.norm(np.stack([Ch.chebval(xx, co[:, k]) for k in range(3)], -1), axis=1)
        i = mag.argmax()
        x0 = 0.5 * (a + b) + 0.5 * (b - a) * xx[i]
        best = max(best, scale, polish(e, a, b, x0, (b - a) / n))
    return best


def truth(e):
    t = np.linspace(0, 1, 400001)
    v = e.g(t)
    i = v.argmax()
    return polish(e, 0, 1, t[i], 3e-6)


CASES = [
    ("bspline fit 5 knots", dict(kind="bspline", fit_knots=5)),
    ("bspline fit 20 knots", dict(kind="bspline", fit_knots=20)),
    ("cylinder fit 6 knots", dict(kind="cyl", fit_knots=6)),
    ("cylinder fit 30 knots", dict(kind="cyl", fit_knots=30)),
]
METHODS = [("uniform23", m_uniform), ("golden64", m_golden), ("split+8", m_split_k),
           ("split+3deg+Brent", m_occt_like), ("split+cheb+polish", m_cheb),
           ("split+8+polish", lambda e: m_split_polish(e))]

if __name__ == "__main__":
    rows = []
    for name, kw in CASES:
        base = Edge(**kw)
        tv = truth(base)
        rows.append((name, kw, tv))
        # narrow defect: Gaussian bump of 3x the fit residual, width 0.002, off any knot
        for wd in (2e-3, 5e-4, 1e-4):
            rows.append((name + f" + bump w={wd:g}", dict(kw, bump=(3 * tv, 0.4137, wd)), None))
    for name, kw, tv in rows:
        e = Edge(**kw)
        tv = truth(e)
        out = []
        for mname, m in METHODS:
            e.calls = 0
            val = m(e)
            out.append(f"{mname}={val / tv:.3f}({e.calls})")
        print(f"{name:34s} true={tv:.3e} | " + "  ".join(out))
