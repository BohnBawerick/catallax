// Test-only oracle driver for the catallax-trace surface cache gate.
// Reads one request line on stdin and prints one response line in the token
// format of ../foundation/oracle.cpp. A surface is `du dv nu nv rational
// knots_u... knots_v... (x y z w)...` with every knot listed, Euclidean
// coordinates, weights, and v varying fastest. Catallax never compiles or
// links this file.
#include "StdAfx.h"
#include <SmBSplineSurface.h>
#include <SmCacheMgr.h>
#include <SmContext.h>
#include <SmNurbsSrf.h>
#include <SmPseudoBox.h>
#include <SmPolarBox.h>
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

static void o3e(const SmExtent3d &e) { o3(e.GetMin()); o3(e.GetMax()); }
static void o2e(const SmExtent2d &e) { od(e.GetMin().x); od(e.GetMin().y); od(e.GetMax().x); od(e.GetMax().y); }
static void opseudo(const SmPseudoBox &b) {
    for (ULONG i = 0; i < 3; ++i) o3(b.GetBasis(i));
    SmExtent1d a, c, d;
    b.GetIntervals(a, c, d);
    od(a.GetMin()); od(a.GetMax()); od(c.GetMin()); od(c.GetMax()); od(d.GetMin()); od(d.GetMax());
}

// Preorder: aux type, box, then patch or head data, then children.
static void dump_node(const SmTreeNode *node) {
    oi((long)node->m_eAuxDataType);
    o3e(node->m_sBBox);
    if (node->m_eAuxDataType == SM_AD_BEZIER_SURFACE) {
        SmBezierPatch *patch = (SmBezierPatch *)node->m_pData;
        o2e(patch->GetUVDomain());
        opseudo(patch->GetPseudoBox());
        oi(patch->m_eUCurveConts[0]); oi(patch->m_eUCurveConts[1]);
        oi(patch->m_eVCurveConts[0]); oi(patch->m_eVCurveConts[1]);
        oi(patch->m_eSplitDir);
        oi(patch->m_bChordHeightSatisfied ? 1 : 0); oi(patch->m_bAngleTolSatisfied ? 1 : 0);
    } else {
        SmBezierAux2d *aux = (SmBezierAux2d *)node->m_pData;
        o2e(aux->GetUVDomain());
        oi(aux->m_eSplitDir);
    }
    oi(node->m_pChild1 ? 1 : 0);
    if (node->m_pChild1) { dump_node(node->m_pChild1); dump_node(node->m_pChild2); }
}

static void opolar(const SmPolarBox &b) {
    oi((long)b.GetPolarBoxState());
    if (b.GetPolarBoxState() == SM_PS_UNINITIALIZED || b.GetPolarBoxState() == SM_PS_UNBOUNDED) return;
    for (ULONG i = 0; i < 3; ++i) o3(*b.GetBasis(i));
    if (b.GetPolarBoxState() != SM_PS_SINGLE_POINT && b.GetPolarBoxState() != SM_PS_POINT_WITH_BASIS)
        o2e(b.GetDomain());
}

// The leaves in preorder, each with the polar box GetPolarBox builds on
// demand (which fills in the ancestors it walks through, in this order).
static void polar_leaves(const SmTreeNode *node, const SmSurface *surface) {
    if (node->m_pChild1) {
        polar_leaves(node->m_pChild1, surface);
        polar_leaves(node->m_pChild2, surface);
        return;
    }
    SmBezierAux2d *aux = (SmBezierAux2d *)node->m_pData;
    // Every node's surface is the original one (checked when the fixture was
    // made), so each box is built from it over the node's domain.
    for (const SmTreeNode *up = node; up; up = up->m_pParent)
        if (((const SmBezierAux2d *)up->m_pData)->mBA_pSurface != surface) { std::fprintf(stderr, "own surface\n"); std::exit(5); }
    opolar(aux->GetPolarBox());
    for (const SmVector3d &v : {SmVector3d(1, 0, 0), SmVector3d(0, 1, 0), SmVector3d(0, 0, 1), SmVector3d(1, 1, 1),
                                SmVector3d(0.3, -0.2, 0.05)})
        oi(aux->GetPolarBox().HasPerpendicularToVector(v) ? 1 : 0);
}

static void conts(const SmTArray<SmContinuityType> &c) {
    oi((long)c.GetSize());
    for (ULONG i = 0; i < c.GetSize(); ++i) oi(c[i]);
}

int main() {
    SmContext context;
    std::string line;
    std::getline(std::cin, line);
    in.str(line);
    std::string op;
    in >> op;
    std::printf("%s", op.c_str());
    SmBSplineSurface *surface = rsurface(context);
    if (op == "cache") {
        SmTrimSrfCache *cache = (SmTrimSrfCache *)SmCacheMgr::GetOrCreateObjectCache(SM_OC_TRIMSRF, surface);
        if (!cache) { oi(SM_ERR); std::printf("\n"); return 0; }
        oi(SM_SUCCESS);
        conts(cache->GetUContinuitiesArray());
        conts(cache->GetVContinuitiesArray());
        dump_node(cache->GetTree()->GetTopNode());
    } else if (op == "polar") {
        SmTrimSrfCache *cache = (SmTrimSrfCache *)SmCacheMgr::GetOrCreateObjectCache(SM_OC_TRIMSRF, surface);
        if (!cache) { oi(SM_ERR); std::printf("\n"); return 0; }
        oi(SM_SUCCESS);
        polar_leaves(cache->GetTree()->GetTopNode(), surface);
    } else if (op == "polarbox") {
        // The polar box over a domain, without a guess and with the whole
        // surface's box as the guess.
        double a = rd(), b = rd(), c = rd(), d = rd();
        SmExtent2d domain(SmPoint2d(a, b), SmPoint2d(c, d));
        SmPolarBox whole, plain, guided;
        if (!status(surface->CalculateBoundingBox(surface->GetNaturalUVDomain(), NULL, NULL, &whole))) { std::printf("\n"); return 0; }
        opolar(whole);
        if (status(surface->CalculateBoundingBox(domain, NULL, NULL, &plain))) opolar(plain);
        if (status(surface->CalculateBoundingBox(domain, NULL, NULL, &guided, NULL, &whole))) opolar(guided);
    } else if (op == "refuse") {
        // The port refuses this cache; SMLib's tree stays in the raw evidence.
        SmTrimSrfCache *cache = (SmTrimSrfCache *)SmCacheMgr::GetOrCreateObjectCache(SM_OC_TRIMSRF, surface);
        std::fprintf(stderr, "cache %s\n", cache ? "built" : "failed");
        oi(-1);
    } else if (op == "net") {
        // sm_ComputeNetConstants on the control net, as the cache asks:
        // chords with the diagonal, then angles alone.
        SmZoneTol3d zone = SmTol::GetZoneTol3d(surface);
        double u = 0, v = 0, uv = 0, ua = 0, va = 0;
        gw_SURFACE *net = surface->GetOrCreateGwNurbPointer();
        if (status(sm_ComputeNetConstants(NULL, NULL, net, zone, &u, &v, NULL, NULL, &uv))) { od(u); od(v); od(uv); }
        if (status(sm_ComputeNetConstants(NULL, NULL, net, zone, NULL, NULL, &ua, &va, NULL))) { od(ua); od(va); }
    } else if (op == "cont") {
        for (int axis = 0; axis < 2; ++axis) {
            SmContinuityType minimum = SM_CT_UNDEFINED;
            SmTArray<SmContinuityType> c;
            if (status(surface->CalculateContinuities(axis == 0 ? SM_SP_U : SM_SP_V, minimum, c))) {
                oi(minimum);
                conts(c);
            }
        }
    } else if (op == "planes") {
        SmPatchBoundaryBoundingPlanes planes(surface, surface->GetNaturalUVDomain());
        oi(SM_SUCCESS);
        for (ULONG i = 0; i < 4; ++i) {
            SmBoolean exists = FALSE, degenerate = FALSE;
            SmPoint3d base;
            SmVector3d normal, *samples = NULL;
            double tolerance = 0.0;
            if (!status(planes.GetNthBoundingPlane(i, exists, degenerate, base, normal, samples, tolerance))) continue;
            oi(exists ? 1 : 0); oi(degenerate ? 1 : 0);
            if (exists) {
                o3(base); o3(normal); od(tolerance);
                for (int k = 0; k < 5; ++k) o3(samples[k]);
            }
        }
    } else if (op == "bbox") {
        SmExtent3d box;
        SmPseudoBox pseudo;
        if (status(surface->CalculateBoundingBox(surface->GetNaturalUVDomain(), &box, &pseudo, NULL))) {
            o3e(box);
            opseudo(pseudo);
        }
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
