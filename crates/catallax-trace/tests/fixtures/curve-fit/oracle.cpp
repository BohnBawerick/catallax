// Test-only oracle driver for the catallax-trace curve-fit gate.
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

int main() {
    SmContext context;
    std::string line, op;
    std::getline(std::cin, line);
    in.str(line);
    in >> op;
    std::printf("%s", op.c_str());
    if (op == "interpolate") {
        long par = ri(), dim = ri(), degree = ri();
        long npoints = ri();
        SmTArray<SmPoint3d> points;
        for (long i = 0; i < npoints; ++i) points.Add(r3());
        long nvectors = ri();
        SmTArray<SmVector3d> vectors;
        for (long i = 0; i < nvectors; ++i) vectors.Add(r3());
        long ignore = ri();
        SmBSplineCurve *c = NULL;
        SmCurveParameterizationType type = par == 0 ? SM_CP_UNIFORM : par == 1 ? SM_CP_CHORDLENGTH : SM_CP_CENTRIPETAL;
        if (status(SmBSplineCurve::CreateInterpolatingCurve(context, type, dim, degree, points, vectors, NULL,
                                                            ignore ? TRUE : FALSE, c))) {
            raw_curve(*c);
            oi(vectors.GetSize());
            for (ULONG i = 0; i < vectors.GetSize(); ++i) { od(vectors[i].x); od(vectors[i].y); od(vectors[i].z); }
        }
    } else if (op == "fair") {
        SmTArray<SmPoint3d> points;
        SmTArray<SmVector3d> tangents, higher;
        for (int i = 0; i < 2; ++i) points.Add(r3());
        for (int i = 0; i < 2; ++i) tangents.Add(r3());
        long nhigher = ri();
        for (long i = 0; i < nhigher; ++i) higher.Add(r3());
        long prevent = ri();
        if (status(SmBSplineCurve::FairBlendCurveDerivatives(points, tangents, nhigher ? &higher : NULL, prevent ? TRUE : FALSE))) {
            for (int i = 0; i < 2; ++i) { od(tangents[i].x); od(tangents[i].y); od(tangents[i].z); }
            oi(higher.GetSize());
            for (ULONG i = 0; i < higher.GetSize(); ++i) { od(higher[i].x); od(higher[i].y); od(higher[i].z); }
        }
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
