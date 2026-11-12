// Test-only oracle driver for the catallax-trace fillet-trace gate.
// Reads one request line on stdin: an operation name and its inputs as
// 16-digit hexadecimal binary64 bit patterns or decimal integers. Prints one
// response line: tokens `d<hex>` for binary64 values and `i<decimal>` for
// integers, statuses, booleans and enumerations. After a failing status
// only the status is printed. Catallax never compiles or links this file.
#include "StdAfx.h"
#include <SmContext.h>
#include <SmExtent1d.h>
#include <SmVector3d.h>
#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmNurbsSrf.h>
#include <SmExtent2d.h>
#include <SmOffsetSurface.h>
#include <SmFilletStandardSolver.h>
#include <SmSurfaceIntersector.h>
#include <SmFilletIntersector.h>
#include <SmFilletExecutive.h>
#include <SmAxis2Placement.h>
#include <SmTArray.h>
#include <SmNurbsCrv.h>
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
static SmVector3d r3() { double x = rd(), y = rd(), z = rd(); return SmVector3d(x, y, z); }
static void od(double value) {
    uint64_t bits;
    std::memcpy(&bits, &value, sizeof bits);
    std::printf(" d%016" PRIx64, bits);
}
static void oi(long value) { std::printf(" i%ld", value); }
static bool status(SmStatus s) { oi((long)s); return s == SM_SUCCESS; }
// A setup call whose failure is the case's answer; success prints nothing.
static bool setup(SmStatus s) { if (s != SM_SUCCESS) oi((long)s); return s == SM_SUCCESS; }

// The curve as stored: dimension, degree, counts, every knot, every control
// point's four stored coordinates.
static void raw_curve(const SmBSplineCurve &c) {
    gw_CURVE *g = c.GetGwNurbPointer();
    oi((long)c.GetDim());
    oi(g->p); oi(g->pol->n + 1); oi(g->knt->m + 1);
    for (long i = 0; i <= g->knt->m; ++i) od(g->knt->U[i]);
    for (long i = 0; i <= g->pol->n; ++i) {
        od(g->pol->Pw[i].x); od(g->pol->Pw[i].y); od(g->pol->Pw[i].z); od(g->pol->Pw[i].w);
    }
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

static void raw_surface(gw_SURFACE *s) {
    oi(s->p); oi(s->q); oi(s->net->n + 1); oi(s->net->m + 1);
    oi(s->knu->m + 1); oi(s->knv->m + 1);
    for (long i = 0; i <= s->knu->m; ++i) od(s->knu->U[i]);
    for (long i = 0; i <= s->knv->m; ++i) od(s->knv->U[i]);
    for (long i = 0; i <= s->net->n; ++i) for (long j = 0; j <= s->net->m; ++j) {
        gw_CPOINT p = s->net->Pw[i][j]; od(p.x); od(p.y); od(p.z); od(p.w);
    }
}
static gw_SURFACE *read_surface() {
    long p = ri(), q = ri(), n = ri(), m = ri();
    ri(); // rational flag, unused: weights are read as stored
    gw_SURFACE *raw = sm_AllocateNurbSurface(n - 1, m - 1, p, q, n + p, m + q);
    for (long i = 0; i < n + p + 1; ++i) raw->knu->U[i] = rd();
    for (long i = 0; i < m + q + 1; ++i) raw->knv->U[i] = rd();
    for (long i = 0; i < n; ++i) for (long j = 0; j < m; ++j) {
        raw->net->Pw[i][j].x = rd(); raw->net->Pw[i][j].y = rd();
        raw->net->Pw[i][j].z = rd(); raw->net->Pw[i][j].w = rd();
    }
    return raw;
}



static SmExtent2d read_domain() {
    double a = rd(), b = rd(), c = rd(), d = rd();
    return SmExtent2d(SmPoint2d(a, b), SmPoint2d(c, d));
}

static std::string g_dump;
static void dd(double value) { char b[40]; uint64_t bits; std::memcpy(&bits, &value, sizeof bits); std::snprintf(b, sizeof b, " d%016" PRIx64, bits); g_dump += b; }
static void di(long value) { char b[40]; std::snprintf(b, sizeof b, " i%ld", value); g_dump += b; }

// The fillet intersector with FlushCurve replaced: it prints the traced
// points the real FlushCurve would turn into a fillet surface.
struct DumpIntersector : public SmFilletIntersector {
    bool flushed = false;
    // a step a domain clipped without ending on its boundary, which SMLib re-solves
    // with LocalCurveIntersect on an offset surface
    bool saw_clipped = false;
    using SmFilletIntersector::SmFilletIntersector;
    SmStatus ComputeNextPoint(SmTsectPnt &t, double step, SmTsectPnt &next, SmBoolean &good, double &dev,
                              double &ang, SmBoolean &clipped, SmBoolean &hit) override {
        SmStatus st = SmFilletIntersector::ComputeNextPoint(t, step, next, good, dev, ang, clipped, hit);
        if (st == SM_SUCCESS && clipped && !hit) saw_clipped = true;
        return st;
    }
    SmStatus FlushCurve(SmTArray<SmCurve*> &, SmTArray<SmCurve*> &, SmTArray<SmCurve*> &,
                        SmTArray<SmTsectCurveType> &, SmTArray<double> &) override {
        flushed = true;
        SmTsectPnt *data[400];
        SmTArray<SmTsectPnt*> points(400, data);
        m_vCurvePoints.GetAllNodes(points);
        di((long)points.GetSize());
        for (ULONG i = 0; i < points.GetSize(); ++i) {
            SmTsectPnt &t = *points[i];
            di((long)t.m_ePointType);
            for (int s = 0; s < 2; ++s) {
                dd(t.UVPos(s).x); dd(t.UVPos(s).y); dd(t.UVDeriv(s).x); dd(t.UVDeriv(s).y);
            }
            dd(t.CrvPos().x); dd(t.CrvPos().y); dd(t.CrvPos().z);
            dd(t.CrvDeriv().x); dd(t.CrvDeriv().y); dd(t.CrvDeriv().z);
            dd(t.m_dCurveParameter); dd(t.m_dDeviation); dd(t.m_dTangentPlaneAngleRad);
        }
        di(m_bCurveIsClosed);
        m_vCurvePoints.Init();
        return SM_SUCCESS;
    }
};

int main() {
    SmContext context;
    std::string line, op;
    std::getline(std::cin, line);
    in.str(line);
    in >> op;
    std::printf("%s", op.c_str());
    // head: r1 r2 outward1 outward2 solver-tolerance surface1 surface2 domain1 domain2
    double r1 = rd(), r2 = rd();
    long o1 = ri(), o2 = ri();
    double tol = rd();
    gw_SURFACE *raw1 = read_surface();
    gw_SURFACE *raw2 = read_surface();
    SmBSplineSurface s1(raw1, TRUE, &context), s2(raw2, TRUE, &context);
    SmConstantRadiusFS fs(context, tol, 0.01, 0.001, r1, r2, s1, s2, o1 ? TRUE : FALSE, o2 ? TRUE : FALSE);
    SmExtent2d d1 = read_domain(), d2 = read_domain();
    DumpIntersector fi(*fs.GetSurface(0), d1, *fs.GetSurface(1), d2, fs);
    // request: extend-before extend-after has-start-dir [dir] has-end-dir [dir]
    //          has-approx [tol] has-angle [tol] extension-distance count (uv uv)...
    long extend_before = ri(), extend_after = ri();
    SmVector3d start_dir, end_dir;
    long has_start = ri(); if (has_start) start_dir = r3();
    long has_end = ri(); if (has_end) end_dir = r3();
    double approx = 0.0, angle = 0.0;
    long has_approx = ri(); if (has_approx) approx = rd();
    long has_angle = ri(); if (has_angle) angle = rd();
    fi.m_dExtensionDistance = rd();
    long count = ri();
    std::vector<SmTsectPnt> points(count);
    SmTArray<SmTsectPnt*> through;
    for (long i = 0; i < count; ++i) {
        double a = rd(), b = rd(), c = rd(), d = rd();
        points[i].UVPos(0) = SmPoint2d(a, b);
        points[i].UVPos(1) = SmPoint2d(c, d);
    }
    for (long i = 0; i < count; ++i) through.Add(&points[i]);
    SmBSplineCurve *c3 = NULL, *u1 = NULL, *u2 = NULL;
    SmTsectCurveType type;
    double deviation = 0.0;
    SmStatus st = fi.DoPointIntersection(context, through, extend_before ? TRUE : FALSE, extend_after ? TRUE : FALSE,
                                         has_start ? &start_dir : NULL, has_end ? &end_dir : NULL,
                                         has_approx ? SM_CAST_APPROXTOL3D_PTR(&approx) : NULL, has_angle ? &angle : NULL,
                                         c3, u1, u2, type, deviation);
    // the status first, then whether FlushCurve ran, then what it saw
    oi((long)st);
    oi(fi.saw_clipped ? 1 : 0);
    oi(fi.flushed ? 1 : 0);
    std::printf("%s", g_dump.c_str());
    for (long i = 0; i < count; ++i) { oi((long)points[i].m_ePointType); od(points[i].m_dCurveParameter); }
    std::printf("\n");
    return 0;
}
