// Test-only oracle driver for the catallax-trace offset-surface gate.
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
    if (op == "offset-point" || op == "offset-normal" || op == "offset-eval") {
        double distance = rd();
        gw_SURFACE *raw = read_surface();
        double u = rd(), v = rd();
        SmBSplineSurface base(raw, TRUE, &context);
        SmOffsetSurface offset(distance, base);
        SmPoint2d uv(u, v);
        if (op == "offset-point") {
            SmPoint3d p;
            if (status(offset.EvaluatePoint(uv, p))) { od(p.x); od(p.y); od(p.z); }
        } else if (op == "offset-normal") {
            SmVector3d n;
            long uleft = ri(), vleft = ri();
            if (status(offset.EvaluateNormal(uv, uleft ? TRUE : FALSE, vleft ? TRUE : FALSE, n))) { od(n.x); od(n.y); od(n.z); }
        } else {
            long simple = ri(), hu = ri(), hv = ri(), uleft = ri(), vleft = ri(), upper = ri(), nonzero = ri();
            SmVector3d cells[64];
            SmStatus s = simple
                ? offset.EvaluateSimple(uv, hu, hv, uleft ? TRUE : FALSE, vleft ? TRUE : FALSE, upper ? TRUE : FALSE, cells)
                : offset.Evaluate(uv, hu, hv, uleft ? TRUE : FALSE, vleft ? TRUE : FALSE, upper ? TRUE : FALSE, cells, nonzero ? TRUE : FALSE);
            if (status(s)) for (long i = 0; i < (hu + 1) * (hv + 1); ++i) { od(cells[i].x); od(cells[i].y); od(cells[i].z); }
        }
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
