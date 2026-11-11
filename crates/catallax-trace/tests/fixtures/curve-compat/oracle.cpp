// Test-only oracle driver for the catallax-trace curve-compatibility gate.
// Reads one request line on stdin: an operation name and its inputs as
// 16-digit hexadecimal binary64 bit patterns or decimal integers. Prints one
// response line: tokens `d<hex>` for binary64 values and `i<decimal>` for
// integers, statuses, booleans and enumerations. After a failing status
// only the status is printed. Catallax never compiles or links this file.
#include "StdAfx.h"
#include <SmContext.h>
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
static void od(double value) {
    uint64_t bits;
    std::memcpy(&bits, &value, sizeof bits);
    std::printf(" d%016" PRIx64, bits);
}
static void oi(long value) { std::printf(" i%ld", value); }
// NLib's flag: 0 success, 1 error. SMLib reports an error as SM_ERR.
static bool nl_status(NL_FLAG flag) { oi(flag == NL_YES ? 1001 : 1000); return flag != NL_YES; }

static void raw_curve(gw_CURVE *g, long dim) {
    oi(dim);
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
    if (op == "elevate") {
        long dim;
        gw_CURVE *raw = read_curve(dim);
        long t = ri();
        NL_STACKS stack;
        N_InitNurbs(&stack);
        NL_CURVE out;
        N_CrvInitArrays(&out);
        if (nl_status(N_CrvElevateDegree((NL_CURVE *)raw, t, &out, &stack, &stack))) raw_curve((gw_CURVE *)&out, dim);
    } else if (op == "remove") {
        long dim;
        gw_CURVE *raw = read_curve(dim);
        double u = rd();
        long times = ri();
        double tol = rd();
        NL_STACKS stack;
        N_InitNurbs(&stack);
        NL_CURVE out;
        N_CrvInitArrays(&out);
        NL_INDEX removed = 0;
        if (nl_status(N_CrvRemoveKnot((NL_CURVE *)raw, u, times, tol, &removed, &out, &stack))) {
            oi(removed);
            raw_curve((gw_CURVE *)&out, dim);
        }
    } else if (op == "knots") {
        long count = ri();
        double tol = rd();
        std::vector<std::vector<double>> values(count);
        std::vector<NL_KNOTVECTOR> vectors(count);
        std::vector<NL_KNOTVECTOR *> pointers(count);
        for (long r = 0; r < count; ++r) {
            long length = ri();
            values[r].resize(length);
            for (auto &v : values[r]) v = rd();
            vectors[r].m = length - 1;
            vectors[r].U = values[r].data();
            pointers[r] = &vectors[r];
        }
        NL_STACKS stack;
        N_InitNurbs(&stack);
        NL_KNOTVECTOR **additions = NULL;
        if (nl_status(N_GetCompatibleKnotVectorToTol(pointers.data(), count - 1, tol, &additions, &stack))) {
            for (long r = 0; r < count; ++r) {
                oi(vectors[r].m + 1);
                for (long i = 0; i <= vectors[r].m; ++i) od(vectors[r].U[i]);
                oi(additions[r]->m + 1);
                for (long i = 0; i <= additions[r]->m; ++i) od(additions[r]->U[i]);
            }
        }
    } else if (op == "sync") {
        long count = ri();
        double tol = rd();
        SmTArray<SmBSplineCurve *> curves;
        std::vector<long> dims;
        for (long r = 0; r < count; ++r) {
            long dim;
            gw_CURVE *raw = read_curve(dim);
            dims.push_back(dim);
            curves.Add(new (context) SmBSplineCurve(raw, dim, TRUE, &context));
        }
        SmStatus s = SmBSplineCurve::SyncronizeKnotsOfCurves(curves, tol);
        oi((long)s);
        if (s == SM_SUCCESS)
            for (long r = 0; r < count; ++r) raw_curve(curves[r]->GetGwNurbPointer(), curves[r]->GetDim());
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
