// Test-only oracle driver for the catallax-trace curve-construction gate.
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

int main() {
    SmContext context;
    std::string line, op;
    std::getline(std::cin, line);
    in.str(line);
    in >> op;
    std::printf("%s", op.c_str());
    if (op == "line") {
        long dim = ri();
        SmPoint3d a = r3(), b = r3();
        SmBSplineCurve *c = NULL;
        if (status(SmBSplineCurve::CreateLineSegment(context, dim, a, b, c))) raw_curve(*c);
    } else if (op == "point") {
        long dim = ri();
        SmPoint3d a = r3();
        SmBSplineCurve *c = NULL;
        if (status(SmBSplineCurve::CreateDegenerateCurve(context, dim, a, c))) raw_curve(*c);
    } else if (op == "insert") {
        long dim;
        gw_CURVE *raw = read_curve(dim);
        double u = rd();
        long times = ri();
        SmBSplineCurve c(raw, dim, TRUE, &context);
        if (status(c.InsertOneKnot(u, times))) raw_curve(c);
    } else if (op == "refine") {
        long dim;
        gw_CURVE *raw = read_curve(dim);
        long count = ri();
        SmTArray<double> knots;
        for (long i = 0; i < count; ++i) knots.Add(rd());
        SmBSplineCurve c(raw, dim, TRUE, &context);
        if (status(c.InsertKnots(knots))) raw_curve(c);
    } else if (op == "reparam" || op == "reparam-nonincreasing") {
        long dim;
        gw_CURVE *raw = read_curve(dim);
        double lo = rd(), hi = rd();
        SmBSplineCurve c(raw, dim, TRUE, &context);
        if (status(c.EditParameterization(SmExtent1d(lo, hi), FALSE))) raw_curve(c);
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
