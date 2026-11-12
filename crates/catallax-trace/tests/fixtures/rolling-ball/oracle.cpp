// Test-only oracle driver for the catallax-trace rolling-ball gate.
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

int main() {
    SmContext context;
    std::string line, op;
    std::getline(std::cin, line);
    in.str(line);
    in >> op;
    std::printf("%s", op.c_str());
    if (op == "rb-surface" || op == "rb-point" || op == "rb-step") {
        double r1 = rd(), r2 = rd();
        long o1 = ri(), o2 = ri();
        double tol = rd();
        gw_SURFACE *raw1 = read_surface();
        gw_SURFACE *raw2 = read_surface();
        SmBSplineSurface s1(raw1, TRUE, &context), s2(raw2, TRUE, &context);
        SmConstantRadiusFS fs(context, tol, 0.01, 0.001, r1, r2, s1, s2, o1 ? TRUE : FALSE, o2 ? TRUE : FALSE);
        SmTsectPnt p;
        auto read_domain = []() { double a = rd(), b = rd(), c = rd(), d = rd(); return SmExtent2d(SmPoint2d(a, b), SmPoint2d(c, d)); };
        auto print_point = [&](SmTsectPnt &t) {
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
        };
        if (op == "rb-surface") {
            SmPoint2d uv[2];
            for (int s = 0; s < 2; ++s) { double u = rd(), v = rd(); uv[s] = SmPoint2d(u, v); }
            if (status(fs.ComputeSurfaceValues(uv, p))) print_point(p);
        } else if (op == "rb-point") {
            SmPoint2d uv[2];
            for (int s = 0; s < 2; ++s) { double u = rd(), v = rd(); uv[s] = SmPoint2d(u, v); }
            double dir = rd();
            SmExtent2d d1 = read_domain(), d2 = read_domain();
            long has_prev = ri();
            SmTsectPnt prev;
            if (has_prev) { prev.CrvDeriv() = r3(); }
            if (status(fs.ComputePointValues(uv, d1, d2, dir, p, has_prev ? &prev : NULL))) print_point(p);
        } else {
            for (int s = 0; s < 2; ++s) { double u = rd(), v = rd(); p.UVPos(s) = SmPoint2d(u, v); }
            for (int s = 0; s < 2; ++s) { double u = rd(), v = rd(); p.UVDeriv(s) = SmVector2d(u, v); }
            p.CrvPos() = r3();
            p.CrvDeriv() = r3();
            double step = rd();
            SmExtent2d d1 = read_domain(), d2 = read_domain();
            SmBoolean found, clipped, hit;
            SmTsectPnt next;
            if (status(fs.StepSolve(p, step, d1, d2, found, clipped, hit, next))) {
                oi(found); oi(clipped); oi(hit); od(step);
                if (found) { od(next.UVPos(0).x); od(next.UVPos(0).y); od(next.UVPos(1).x); od(next.UVPos(1).y); }
            }
        }
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
