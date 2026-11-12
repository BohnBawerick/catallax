// Test-only oracle driver for the catallax-trace fillet-intersector gate.
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


static void print_point(SmTsectPnt &t) {
    oi((long)t.m_ePointType);
    for (int s = 0; s < 2; ++s) {
        od(t.SrfPos(s).x); od(t.SrfPos(s).y); od(t.SrfPos(s).z);
        od(t.SrfDu(s).x); od(t.SrfDu(s).y); od(t.SrfDu(s).z);
        od(t.SrfDv(s).x); od(t.SrfDv(s).y); od(t.SrfDv(s).z);
        od(t.SrfNorm(s).x); od(t.SrfNorm(s).y); od(t.SrfNorm(s).z);
        od(t.UVPos(s).x); od(t.UVPos(s).y); od(t.UVDeriv(s).x); od(t.UVDeriv(s).y);
    }
    od(t.CrvPos().x); od(t.CrvPos().y); od(t.CrvPos().z);
    od(t.CrvDeriv().x); od(t.CrvDeriv().y); od(t.CrvDeriv().z);
    od(t.m_dTangentPlaneAngleRad);
    od(t.m_dCurveParameter);
}

static SmExtent2d read_domain() {
    double a = rd(), b = rd(), c = rd(), d = rd();
    return SmExtent2d(SmPoint2d(a, b), SmPoint2d(c, d));
}

int main() {
    SmContext context;
    std::string line, op;
    std::getline(std::cin, line);
    in.str(line);
    in >> op;
    std::printf("%s", op.c_str());
    // head: r1 r2 outward1 outward2 solver-tolerance surface1 surface2 domain1 domain2
    //       approx-tolerance angle-tolerance trace-direction through-point-flag
    double r1 = rd(), r2 = rd();
    long o1 = ri(), o2 = ri();
    double tol = rd();
    gw_SURFACE *raw1 = read_surface();
    gw_SURFACE *raw2 = read_surface();
    SmBSplineSurface s1(raw1, TRUE, &context), s2(raw2, TRUE, &context);
    SmConstantRadiusFS fs(context, tol, 0.01, 0.001, r1, r2, s1, s2, o1 ? TRUE : FALSE, o2 ? TRUE : FALSE);
    SmExtent2d d1 = read_domain(), d2 = read_domain();
    SmFilletIntersector fi(*fs.GetSurface(0), d1, *fs.GetSurface(1), d2, fs);
    if (op == "fi-ctor") {
        od(fi.m_dThisApproxTol3d); od(fi.m_dThisAngTolRad); od(fi.m_dCurveTraceDirection);
        oi(fi.m_bClippedByThroughPoint);
        std::printf("\n");
        return 0;
    }
    double approx = rd(), angle = rd(), direction = rd();
    long through = ri();
    fi.m_dThisApproxTol3d = approx;
    fi.m_dThisAngTolRad = angle;
    fi.m_dCurveTraceDirection = direction;
    fi.m_bClippedByThroughPoint = through ? TRUE : FALSE;
    auto read_uvs = [](SmPoint2d uv[2]) {
        for (int s = 0; s < 2; ++s) { double u = rd(), v = rd(); uv[s] = SmPoint2d(u, v); }
    };
    if (op == "fi-values") {
        SmPoint2d uv[2]; read_uvs(uv);
        long has_prev = ri();
        SmTsectPnt prev, p;
        if (has_prev) { prev.CrvDeriv() = r3(); }
        if (status(fi.ComputePointValues(uv, p, has_prev ? &prev : NULL, NULL))) print_point(p);
    } else if (op == "fi-delta") {
        SmPoint2d uv[2], delta[2]; read_uvs(uv); read_uvs(delta);
        double param = rd(), delta_param = rd();
        SmTsectPnt p, dp;
        p.m_dCurveParameter = param;
        dp.m_dCurveParameter = delta_param;
        if (status(fi.ComputePointValuesDelta(uv, delta, p, dp))) { print_point(p); print_point(dp); }
    } else if (op == "fi-law") {
        SmPoint2d uv[2]; read_uvs(uv);
        SmVector3d origin = r3(), normal = r3();
        SmFilletLaw law; // unused by the constant-radius solve
        SmExtent1d ivl(0.0, 1.0);
        double orient[2] = {1.0, 1.0};
        double param = -1.0;
        SmPoint2d out_uv[2] = {SmPoint2d(-7.0, -7.0), SmPoint2d(-7.0, -7.0)};
        if (status(fi.EvaluateLawPoint(0.0, uv[0], uv[1], law, FALSE, ivl, orient, origin, normal, param, out_uv))) {
            od(out_uv[0].x); od(out_uv[0].y); od(out_uv[1].x); od(out_uv[1].y); od(param);
        }
    } else if (op == "fi-span" || op == "fi-next") {
        // the current point comes from ComputePointValues at the given parameters
        SmPoint2d uv[2]; read_uvs(uv);
        SmTsectPnt cur;
        if (!setup(fi.ComputePointValues(uv, cur, NULL, NULL))) { std::printf("\n"); return 0; }
        cur.m_dCurveParameter = rd();
        if (op == "fi-span") {
            SmPoint2d nuv[2]; read_uvs(nuv);
            SmTsectPnt next;
            if (!setup(fi.ComputePointValues(nuv, next, &cur, NULL))) { std::printf("\n"); return 0; }
            next.m_dCurveParameter = rd();
            long has_mid = ri();
            SmTsectPnt mid;
            if (has_mid) {
                SmPoint2d muv[2]; read_uvs(muv);
                if (!setup(fi.ComputePointValues(muv, mid, &cur, NULL))) { std::printf("\n"); return 0; }
            }
            SmBoolean good = FALSE;
            double dev = 0.0, ang = 0.0;
            if (status(fi.TestSpanAccuracy(cur, next, has_mid ? &mid : NULL, good, dev, ang))) {
                oi(good); od(ang); oi((long)next.m_ePointType);
            }
        } else {
            double step = rd();
            SmTsectPnt next;
            SmBoolean good = FALSE, clipped = FALSE, hit = FALSE;
            double dev = 0.0, ang = 0.0;
            if (status(fi.ComputeNextPoint(cur, step, next, good, dev, ang, clipped, hit))) {
                oi(good); oi(clipped); oi(hit); od(dev); od(ang);
                oi((long)next.m_ePointType);
                od(next.UVPos(0).x); od(next.UVPos(0).y); od(next.UVPos(1).x); od(next.UVPos(1).y);
                od(next.UVDeriv(0).x); od(next.UVDeriv(0).y); od(next.UVDeriv(1).x); od(next.UVDeriv(1).y);
                od(next.CrvPos().x); od(next.CrvPos().y); od(next.CrvPos().z);
                od(next.CrvDeriv().x); od(next.CrvDeriv().y); od(next.CrvDeriv().z);
                od(next.m_dCurveParameter);
                // the current point, which the through-point branch may change
                od(cur.UVDeriv(0).x); od(cur.UVDeriv(0).y); od(cur.UVDeriv(1).x); od(cur.UVDeriv(1).y);
                od(cur.CrvDeriv().x); od(cur.CrvDeriv().y); od(cur.CrvDeriv().z);
            }
        }
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
