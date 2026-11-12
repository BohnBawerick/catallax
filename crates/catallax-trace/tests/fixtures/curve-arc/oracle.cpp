// Test-only oracle driver for the catallax-trace circle-arc and end-point gate.
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
static bool nl_status(NL_FLAG flag) { oi(flag == NL_YES ? 1001 : 1000); return flag != NL_YES; }

// The curve as stored: dimension, degree, counts, every knot, every control
// point's four stored coordinates.
static void raw_gw(gw_CURVE *g, long dim) {
    oi(dim);
    oi(g->p); oi(g->pol->n + 1); oi(g->knt->m + 1);
    for (long i = 0; i <= g->knt->m; ++i) od(g->knt->U[i]);
    for (long i = 0; i <= g->pol->n; ++i) {
        od(g->pol->Pw[i].x); od(g->pol->Pw[i].y); od(g->pol->Pw[i].z); od(g->pol->Pw[i].w);
    }
}
static void raw_curve(const SmBSplineCurve &c) { raw_gw(c.GetGwNurbPointer(), (long)c.GetDim()); }

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
    if (op == "circle") {
        long dim = ri();
        SmPoint3d origin = r3();
        SmVector3d x = r3(), y = r3();
        double radius = rd(), start = rd(), end = rd();
        long kind = ri();
        SmAxis2Placement placement;
        if (placement.SetCanonical(origin, x, y) != SM_SUCCESS) { std::fprintf(stderr, "frame\n"); return 4; }
        SmBSplineCurve *c = NULL;
        if (status(SmBSplineCurve::CreateCircleSegment(context, dim, placement, radius, start, end,
                                                       kind == 0 ? SM_CO_QUADRATIC : SM_CO_QUINTIC, c))) raw_curve(*c);
    } else if (op == "derivs" || op == "closest" || op == "extract") {
        long dim;
        gw_CURVE *raw = read_curve(dim);
        if (op == "derivs") {
            double u = rd();
            long left = ri(), order = ri();
            NL_POINT d[20];
            if (nl_status(N_CrvDerivs((NL_CURVE *)raw, u, left ? NL_LEFT : NL_RIGHT, order, d)))
                for (long i = 0; i <= order; ++i) { od(d[i].x); od(d[i].y); od(d[i].z); }
        } else if (op == "closest") {
            NL_POINT p; p.x = rd(); p.y = rd(); p.z = rd();
            double u0 = rd(), top = rd(), toc = rd();
            NL_PARAMETER u; NL_POINT q;
            if (nl_status(N_CrvClosestPt((NL_CURVE *)raw, p, u0, top, toc, &u, &q))) { od(u); od(q.x); od(q.y); od(q.z); }
        } else {
            double ul = rd(), ur = rd();
            NL_STACKS stack; N_InitNurbs(&stack);
            NL_CURVE out; N_CrvInitArrays(&out);
            if (nl_status(N_CrvExtractCrvSeg((NL_CURVE *)raw, ul, ur, &out, &stack, &stack))) raw_gw((gw_CURVE *)&out, dim);
        }
    } else if (op == "edit-end") {
        long dim;
        gw_CURVE *raw = read_curve(dim);
        SmPoint3d point = r3();
        long start = ri(), has_weight = ri();
        double weight = rd();
        SmBSplineCurve c(raw, dim, TRUE, &context);
        if (status(c.EditEndPoint(point, start ? TRUE : FALSE, has_weight ? &weight : NULL,
                                  SM_EE_CONSTRUCTED_END_POINT, FALSE))) raw_curve(c);
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
