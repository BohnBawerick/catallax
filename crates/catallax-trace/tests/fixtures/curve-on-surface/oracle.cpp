// Test-only oracle driver for the catallax-trace curve-on-surface gate.
// Reads one request line on stdin and prints one response line in the token
// format of ../foundation/oracle.cpp. A surface is `du dv nu nv rational
// knots_u... knots_v... (x y z w)...` with every knot listed, Euclidean
// coordinates, weights, and v varying fastest. Catallax never compiles or
// links this file.
#include "StdAfx.h"
#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmCrvOnSurf.h>
#include <SmCurveCache.h>
#include <SmGlobalSolver.h>
#include <SmHermiteCurve.h>
#include <unistd.h>
#include <SmCacheMgr.h>
#include <SmContext.h>
#include <SmNurbsSrf.h>
#include <SmPseudoBox.h>
#include <SmSurfaceCache.h>
#include <SmTol.h>
#include <SmTree.h>
#include <SmTrimSrfCache.h>
#include <nurbs.h>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static std::istringstream in;

static double rd() {
    std::string token;
    if (!(in >> token)) { std::fprintf(stderr, "missing input\n"); std::exit(3); }
    uint64_t bits = std::strtoull(token.c_str(), nullptr, 16);
    double value;
    std::memcpy(&value, &bits, sizeof value);
    return value;
}
static long ri() {
    long value;
    if (!(in >> value)) { std::fprintf(stderr, "missing integer\n"); std::exit(3); }
    return value;
}
static void od(double value) {
    uint64_t bits;
    std::memcpy(&bits, &value, sizeof bits);
    std::printf(" d%016" PRIx64, bits);
}
static void oi(long value) { std::printf(" i%ld", value); }
static void o3(const SmVector3d &v) { od(v.x); od(v.y); od(v.z); }
static bool status(SmStatus s) { oi((long)s); return s == SM_SUCCESS; }

// Unique knots and multiplicities of a full knot vector.
static void unique(const std::vector<double> &knots, SmTArray<double> &values, SmTArray<ULONG> &mults) {
    for (double k : knots) {
        if (values.GetSize() > 0 && values[values.GetSize() - 1] == k) mults[mults.GetSize() - 1]++;
        else { values.Add(k); mults.Add(1); }
    }
}

SmBSplineSurface *rsurface(SmContext &context) {
    long du = ri(), dv = ri(), nu = ri(), nv = ri(), rational = ri();
    std::vector<double> ku(nu + du + 1), kv(nv + dv + 1);
    for (double &k : ku) k = rd();
    for (double &k : kv) k = rd();
    SmTArray<SmPoint3d> points;
    SmTArray<double> weights;
    for (long i = 0; i < nu * nv; ++i) {
        double x = rd(), y = rd(), z = rd(), w = rd();
        points.Add(SmPoint3d(x, y, z));
        weights.Add(w);
    }
    SmTArray<double> vu, vv;
    SmTArray<ULONG> mu, mv;
    unique(ku, vu, mu);
    unique(kv, vv, mv);
    SmBSplineSurface *s = NULL;
    SmStatus st = SmBSplineSurface::CreateCanonical(context, du, dv, points, SM_SF_UNSPECIFIED, mu, mv, vu, vv,
                                                    SM_KT_UNSPECIFIED, rational ? &weights : NULL, NULL, s);
    if (st != SM_SUCCESS || !s) { std::fprintf(stderr, "surface create failed %ld\n", (long)st); std::exit(4); }
    return s;
}

// `dim degree count knots... (x y z w)...`, every value as stored.
static gw_CURVE *read_curve(long &dim) {
    dim = ri();
    long p = ri(), n = ri();
    std::vector<double> u(n + p + 1);
    for (auto &v : u) v = rd();
    gw_CURVE *raw = sm_AllocateNurbCurve(n - 1, p, u.size() - 1);
    for (size_t i = 0; i < u.size(); ++i) raw->knt->U[i] = u[i];
    for (long i = 0; i < n; ++i) {
        raw->pol->Pw[i].x = rd(); raw->pol->Pw[i].y = rd();
        raw->pol->Pw[i].z = rd(); raw->pol->Pw[i].w = rd();
    }
    return raw;
}


// SmGlobalSolver with the branch test exposed.
struct Probe : SmGlobalSolver {
    SmBoolean branch(SmTreeNode **b) { return BranchMayContainAnswers(b); }
};

static void collect(SmTreeNode *node, std::vector<SmTreeNode *> &out) {
    if (!node) return;
    out.push_back(node);
    collect(node->m_pChild1, out);
    collect(node->m_pChild2, out);
}

static SmExtent3d read_box() {
    double a = rd(), b = rd(), c = rd(), d = rd(), e = rd(), f = rd();
    SmExtent3d box;
    box.m_vMin = SmPoint3d(a, b, c);
    box.m_vMax = SmPoint3d(d, e, f);
    return box;
}

int main() {
    SmContext context;
    std::string line, op;
    std::getline(std::cin, line);
    in.str(line);
    in >> op;
    std::printf("%s", op.c_str());
    // surface, then the parameter-space curve: 0 Hermite (p1 d1 p2 d2, planar) or 1 B-spline
    SmBSplineSurface &surface = *rsurface(context);
    long kind = ri();
    SmCurve *uv = NULL;
    if (kind == 0) {
        double a = rd(), b = rd(), c = rd(), d = rd(), e = rd(), f = rd(), g = rd(), h = rd();
        SmHermiteCurve *hermite = new SmHermiteCurve(SmPoint2d(a, b), SmVector2d(c, d), SmPoint2d(e, f), SmVector2d(g, h), 2);
        hermite->SetContext(&context);
        uv = hermite;
    } else {
        long dim;
        gw_CURVE *rawc = read_curve(dim);
        uv = new SmBSplineCurve(rawc, dim, TRUE, &context);
    }
    SmCrvOnSurf cos(*uv, surface, NULL, 0, &context);
    if (op == "cos-eval") {
        double t = rd();
        long left = ri();
        SmPoint3d p;
        if (status(cos.EvaluatePoint(t, p))) { od(p.x); od(p.y); od(p.z); }
        SmPoint2d q;
        if (status(cos.EvaluateUVPoint(t, q))) { od(q.x); od(q.y); }
        SmVector2d d[2];
        SmVector3d n;
        if (status(cos.EvaluateSurfaceNormal(t, left ? TRUE : FALSE, d, n))) {
            od(d[0].x); od(d[0].y); od(d[1].x); od(d[1].y); od(n.x); od(n.y); od(n.z);
        }
    } else if (op == "cos-node" || op == "cos-branch") {
        long index = ri();
        SmCurveCache *cache = (SmCurveCache *)SmCacheMgr::GetOrCreateObjectCache(SM_OC_CURVE, &cos);
        if (!cache) { oi(-1); std::printf("\n"); _exit(0); }
        std::vector<SmTreeNode *> nodes;
        collect(cache->GetTree()->GetTopNode(), nodes);
        if (op == "cos-node") {
            oi((long)nodes.size());
            if (index >= 0 && index < (long)nodes.size()) {
                SmTreeNode *n = nodes[index];
                od(n->m_sBBox.m_vMin.x); od(n->m_sBBox.m_vMin.y); od(n->m_sBBox.m_vMin.z);
                od(n->m_sBBox.m_vMax.x); od(n->m_sBBox.m_vMax.y); od(n->m_sBBox.m_vMax.z);
                SmBezierAux1d *aux = (SmBezierAux1d *)n->m_pData;
                od(aux->GetInterval().GetMin()); od(aux->GetInterval().GetMax());
            }
        } else {
            SmExtent3d box1 = read_box();
            double tolerance = rd();
            SmTreeNode first;
            first.m_sBBox = box1;
            SmTreeNode *branch[SM_GS_MAX_TREES] = {&first, nodes[index], NULL};
            Probe probe;
            probe.SetSolverOperation(SM_SO_SURFACE_PROJECTED_INTERSECT);
            probe.SetNumTrees(2);
            probe.Set3DTolerance(tolerance);
            oi(probe.branch(branch));
        }
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    std::fflush(stdout);
    _exit(0);
}
