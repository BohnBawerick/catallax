// Test-only oracle driver for the catallax-trace curve drop gate. Reads one
// request line on stdin and prints one response line in the token format of
// ../foundation/oracle.cpp. A surface is `s du dv nu nv rational knots_u...
// knots_v... (x y z w)...` as in ../intersector/oracle.cpp, or `p` and a
// plane; a curve is `degree count rational knots... (x y z w)...` as in
// ../curve-surface/oracle.cpp. Catallax never compiles or links this file.
#include "StdAfx.h"
#include <SmBSplineSurface.h>
#include <SmBSplineCurve.h>
#include <SmPlane.h>
#include <SmAxis2Placement.h>
#include <SmContext.h>
#include <SmSurfaceDropCurve.h>
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
static FILE *sink = stdout;
static void od(double value) {
    uint64_t bits;
    std::memcpy(&bits, &value, sizeof bits);
    std::fprintf(sink, " d%016" PRIx64, bits);
}
static void oi(long value) { std::fprintf(sink, " i%ld", value); }
static bool status(SmStatus s) { oi((long)s); return s == SM_SUCCESS; }

static void unique(const std::vector<double> &knots, SmTArray<double> &values, SmTArray<ULONG> &mults) {
    for (double k : knots) {
        if (values.GetSize() > 0 && values[values.GetSize() - 1] == k) mults[mults.GetSize() - 1]++;
        else { values.Add(k); mults.Add(1); }
    }
}

static SmBSplineSurface *rsurface(SmContext &context) {
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

static SmVector3d r3() { double x = rd(), y = rd(), z = rd(); return SmVector3d(x, y, z); }
static SmExtent2d rdomain() { double a = rd(), b = rd(), c = rd(), d = rd(); return SmExtent2d(SmPoint2d(a, b), SmPoint2d(c, d)); }

static SmPlane *rplane(SmContext &context) {
    long kind = ri();
    SmVector3d origin = r3(), x = r3(), y = r3();
    double su = rd(), sv = rd();
    SmExtent2d domain = rdomain();
    SmPlane *p = NULL;
    if (kind == 1) {
        SmAxis2Placement frame;
        frame.SetCanonical(origin, x, y);
        if (SmPlane::CreateCanonical(context, frame, p) != SM_SUCCESS || !p) { std::fprintf(stderr, "plane failed\n"); std::exit(4); }
    } else {
        p = new (context) SmPlane(origin, x, y, SmVector2d(su, sv), domain, &context);
        if (p->MakeNurb() != SM_SUCCESS) { std::fprintf(stderr, "make nurb failed\n"); std::exit(4); }
    }
    return p;
}

static SmSurface *rany(SmContext &context) {
    std::string kind;
    in >> kind;
    if (kind == "p") return rplane(context);
    return rsurface(context);
}

static SmBSplineCurve *rcurve(SmContext &context) {
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
    return new (context) SmBSplineCurve(3, count - 1, count + degree, (short)degree, points.data(), knots.data(),
                                        &context);
}

// A B-spline curve as stored: degree, control point and knot counts, every
// knot, then homogeneous control points.
static void raw(const SmCurve *curve) {
    const gw_CURVE *c = ((SmBSplineCurve *)curve)->GetOrCreateGwNurbPointer();
    oi(c->p); oi(c->pol->n + 1); oi(c->knt->m + 1);
    for (long i = 0; i <= c->knt->m; ++i) od(c->knt->U[i]);
    for (long i = 0; i <= c->pol->n; ++i) {
        od(c->pol->Pw[i].x); od(c->pol->Pw[i].y); od(c->pol->Pw[i].z); od(c->pol->Pw[i].w);
    }
}

int main() {
    SmContext context;
    std::string line;
    std::getline(std::cin, line);
    in.str(line);
    std::string op;
    in >> op;
    std::printf("%s", op.c_str());
    SmSurface *surface = rany(context);
    SmExtent2d domain = rdomain();
    SmBSplineCurve *curve = rcurve(context);
    double t0 = rd(), t1 = rd();
    SmExtent1d interval(t0, t1);
    double tolerance = rd();
    long keep_all = ri(), on_surface = ri();
    if (op == "drop" || op == "drefuse") {
        // SmSurface::DropAndTrimCurve: every parameter curve with the
        // largest drop distance and normal-line deviation.
        if (op == "drefuse") sink = stderr;
        double drop = -1, deviation = -1;
        SmTArray<SmBSplineCurve *> curves;
        if (status(surface->DropAndTrimCurve(context, domain, *curve, interval, tolerance, drop, deviation, curves,
                                             keep_all ? TRUE : FALSE, on_surface ? TRUE : FALSE))) {
            oi((long)curves.GetSize());
            for (ULONG i = 0; i < curves.GetSize(); ++i) raw(curves[i]);
            od(drop);
            od(deviation);
        }
        if (op == "drefuse") { sink = stdout; std::fprintf(stderr, "\n"); oi(-1); }
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
