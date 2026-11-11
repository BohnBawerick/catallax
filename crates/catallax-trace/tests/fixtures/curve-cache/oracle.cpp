// Test-only oracle driver for the catallax-trace curve cache gate.
// Reads one request line on stdin: `cache dim degree n_points weights?
// knots... points...` and prints the cache SMLib builds for that curve
// through its cache manager, in the token format of ../foundation/oracle.cpp.
// Catallax never compiles or links this file.
#include "StdAfx.h"
#include <SmBSplineCurve.h>
#include <nurbs.h>
#include <SmCacheMgr.h>
#include <SmContext.h>
#include <SmCurveCache.h>
#include <SmPseudoBox.h>
#include <SmTol.h>
#include <SmTree.h>
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

static SmCurveCache *g_cache;

static void dump_span(const SmBezierAux1d *span) {
    oi((long)span->m_lLeft); oi((long)span->m_lRight);
    oi((long)span->m_aeConts[0]); oi((long)span->m_aeConts[1]);
    od(span->m_dMaxLength); od(span->m_dMaxChordHeightSquared); od(span->m_dMaxTurningAngleDeg);
    od(((SmBezierAux1d *)span)->GetInterval().GetMin()); od(((SmBezierAux1d *)span)->GetInterval().GetMax());
}

// Preorder walk: box, then the leaf span or parent index the node's data is.
static void dump_node(SmTreeNode *node) {
    o3(node->m_sBBox.GetMin()); o3(node->m_sBBox.GetMax());
    long kind = 0, index = -1;
    for (ULONG i = 0; i < g_cache->GetSpanCount(); ++i)
        if (node->m_pData == g_cache->GetAt(i)) { kind = 1; index = (long)i; }
    if (kind == 0 && node->m_pData)
        for (ULONG i = 0; i + 1 < g_cache->GetSpanCount(); ++i)
            if (node->m_pData == g_cache->GetParentAt(i)) { kind = 2; index = (long)i; }
    oi(kind); oi(index); oi(node->m_pChild1 ? 1 : 0);
    if (node->m_pChild1) { dump_node(node->m_pChild1); dump_node(node->m_pChild2); }
}

int main() {
    SmContext context;
    std::string line;
    std::getline(std::cin, line);
    in.str(line);
    std::string op;
    in >> op;
    std::printf("%s", op.c_str());
    if (op != "cache") { std::fprintf(stderr, "unknown operation %s\n", op.c_str()); return 2; }
    long degree = ri(), count = ri(), rational = ri();
    std::vector<double> knots(count + degree + 1);
    for (double &k : knots) k = rd();
    std::vector<gw_CPOINT> points(count);
    for (gw_CPOINT &p : points) {
        p.x = rd(); p.y = rd(); p.z = rd();
        double w = rd();
        p.w = rational ? w : NL_NOW;
        if (rational) { p.x *= w; p.y *= w; p.z *= w; }
    }
    SmBSplineCurve curve(3, count - 1, count + degree, (short)degree, points.data(), knots.data(), &context);
    od(SmTol::GetApproxTol3d(&curve));
    g_cache = (SmCurveCache *)SmCacheMgr::GetOrCreateObjectCache(SM_OC_CURVE, &curve);
    if (!g_cache) { oi(0); std::printf("\n"); return 0; }
    oi(1);
    const SmTArray<SmContinuityType> &continuities = g_cache->GetContinuitiesArray();
    oi((long)continuities.GetSize());
    for (ULONG i = 0; i < continuities.GetSize(); ++i) oi((long)continuities[i]);
    oi((long)g_cache->GetSpanCount());
    for (ULONG i = 0; i < g_cache->GetSpanCount(); ++i) {
        SmBezierSpan *span = g_cache->GetAt(i);
        dump_span(span);
        SmPseudoBox &box = span->GetPseudoBox();
        for (ULONG b = 0; b < 3; ++b) o3(box.GetBasis(b));
        SmExtent1d a, bb, c;
        box.GetIntervals(a, bb, c);
        od(a.GetMin()); od(a.GetMax()); od(bb.GetMin()); od(bb.GetMax()); od(c.GetMin()); od(c.GetMax());
    }
    for (ULONG i = 0; i + 1 < g_cache->GetSpanCount(); ++i) dump_span(g_cache->GetParentAt(i));
    dump_node(g_cache->GetTree()->GetTopNode());
    std::printf("\n");
    return 0;
}
