// Test-only oracle driver for the catallax-trace offset curve-surface gate.
// Reads one request line on stdin: an operation name and its inputs as
// 16-digit hexadecimal binary64 bit patterns or decimal integers. Prints one
// response line: tokens `d<hex>` for binary64 values and `i<decimal>` for
// integers, statuses, booleans and enumerations. After a failing status
// only the status is printed. Catallax never compiles or links this file.
#include "StdAfx.h"
#include <SmContext.h>
#include <SmExtent1d.h>
#include <SmExtent2d.h>
#include <SmVector3d.h>
#include <SmBSplineSurface.h>
#include <SmOffsetSurface.h>
#include <SmIsoCurve.h>
#include <SmSolutionArray.h>
#include <SmTArray.h>
#include <SmNurbsSrf.h>
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
static bool status(SmStatus s) { oi((long)s); return s == SM_SUCCESS; }

// `du dv nu nv rational knots_u... knots_v... (x y z w)...`, every value as
// stored.
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
    if (op == "oc-sing") {
        // IsSingularity on an offset surface, then GetSingularities with its poles.
        double distance = rd();
        SmBSplineSurface base(read_surface(), TRUE, &context);
        SmOffsetSurface offset(distance, base);
        double u = rd(), v = rd(), tolerance = rd();
        long point_only = ri();
        SmSurfParamType direction = SM_SP_NEITHER;
        SmBoolean singular = offset.IsSingularity(SmPoint2d(u, v), direction, tolerance, point_only ? TRUE : FALSE);
        oi(singular ? 1 : 0); oi((long)direction);
        SmTArray<SmPoint3d> poles;
        // SmOffsetSurface hides the virtual with a cached flag; call SmSurface's.
        ULONG flags = static_cast<SmSurface &>(offset).GetSingularities(&poles);
        oi((long)flags);
        for (ULONG i = 0; i < 4; ++i) {
            if (poles[i].IsInitialized()) { od(poles[i].x); od(poles[i].y); od(poles[i].z); }
            else oi(0);
        }
    } else if (op == "oc-intersect" || op == "oc-solve") {
        // The iso-curve of offset surface A meets offset surface B, as
        // SmFilletSolver::StepSolve asks it for a clipped step.
        double distance_a = rd();
        SmBSplineSurface base_a(read_surface(), TRUE, &context);
        double distance_b = rd();
        SmBSplineSurface base_b(read_surface(), TRUE, &context);
        SmOffsetSurface offset_a(distance_a, base_a), offset_b(distance_b, base_b);
        long axis = ri();
        double iso = rd();
        long clamp = ri();
        SmIsoCurve curve(offset_a, axis == 0 ? SM_SP_U : SM_SP_V, iso, TRUE, clamp ? TRUE : FALSE);
        curve.SetContext(NULL);
        double u0 = rd(), v0 = rd(), u1 = rd(), v1 = rd();
        SmExtent2d domain(SmPoint2d(u0, v0), SmPoint2d(u1, v1));
        SmExtent1d interval = curve.GetNaturalInterval();
        if (op == "oc-intersect") {
            double tolerance = rd(), t = rd(), gu = rd(), gv = rd();
            SmBoolean found = FALSE;
            SmPoint2d uv;
            double tf = 0, deviation = 0;
            if (status(offset_b.LocalCurveIntersect(domain, curve, interval, tolerance, SmPoint2d(gu, gv), t, found,
                                                    uv, tf, deviation))) {
                oi(found ? 1 : 0);
                if (found) { od(tf); od(uv.x); od(uv.y); od(deviation); }
            }
        } else {
            long operation = ri();
            double tolerance = rd(), t = rd(), gu = rd(), gv = rd();
            SmBoolean found = FALSE;
            SmSolution solution;
            if (status(offset_b.LocalCurveSolve(domain, curve, interval, (SmSolverOperationType)operation, tolerance,
                                                NULL, NULL, SmPoint2d(gu, gv), t, found, solution))) {
                oi(found ? 1 : 0);
                if (found) {
                    od(solution.m_vStart.m_adParameters[0]); od(solution.m_vStart.m_adParameters[1]);
                    od(solution.m_vStart.m_adParameters[2]); od(solution.m_vStart.m_dSolutionValue);
                }
            }
        }
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
