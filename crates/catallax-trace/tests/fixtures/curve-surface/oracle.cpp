// Test-only oracle driver for the catallax-trace curve-surface gate.
// Reads one request line on stdin and prints one response line in the token
// format of ../foundation/oracle.cpp. A surface is `du dv nu nv rational
// knots_u... knots_v... (x y z w)...` with every knot listed, Euclidean
// coordinates, weights, and v varying fastest; a curve is `degree count
// rational knots... (x y z w)...`. Catallax never compiles or links this file.
#include "StdAfx.h"
#include <SmBSplineSurface.h>
#include <SmBSplineCurve.h>
#include <SmContext.h>
#include <SmSolutionArray.h>
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
static void o3(const SmVector3d &v) { od(v.x); od(v.y); od(v.z); }
static bool status(SmStatus s) { oi((long)s); return s == SM_SUCCESS; }

// Unique knots and multiplicities of a full knot vector.
static void unique(const std::vector<double> &knots, SmTArray<double> &values, SmTArray<ULONG> &mults) {
    for (double k : knots) {
        if (values.GetSize() > 0 && values[values.GetSize() - 1] == k) mults[mults.GetSize() - 1]++;
        else { values.Add(k); mults.Add(1); }
    }
}

SmBSplineSurface *rsurface(SmContext &context) {
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

int main() {
    SmContext context;
    std::string line;
    std::getline(std::cin, line);
    in.str(line);
    std::string op;
    in >> op;
    std::printf("%s", op.c_str());
    SmBSplineSurface *surface = rsurface(context);
    if (op == "sing") {
        // IsSingularity at a point, then GetSingularities with its poles.
        double u = rd(), v = rd(), tolerance = rd();
        long point_only = ri();
        SmSurfParamType direction = SM_SP_NEITHER;
        SmBoolean singular = surface->IsSingularity(SmPoint2d(u, v), direction, tolerance, point_only);
        oi(singular ? 1 : 0); oi((long)direction);
        SmTArray<SmPoint3d> poles;
        ULONG flags = surface->GetSingularities(&poles);
        oi((long)flags);
        for (ULONG i = 0; i < 4; ++i)
            if (poles[i].IsInitialized()) o3(poles[i]); else oi(0);
        std::printf("\n");
        return 0;
    }
    SmBSplineCurve *curve = rcurve(context);
    double u0 = rd(), v0 = rd(), u1 = rd(), v1 = rd();
    SmExtent2d domain(SmPoint2d(u0, v0), SmPoint2d(u1, v1));
    double t0 = rd(), t1 = rd();
    SmExtent1d interval(t0, t1);
    if (op == "intersect") {
        double tolerance = rd(), t = rd(), gu = rd(), gv = rd();
        SmBoolean found = FALSE;
        SmPoint2d uv;
        double tf = 0, deviation = 0;
        if (status(surface->LocalCurveIntersect(domain, *curve, interval, tolerance, SmPoint2d(gu, gv), t, found, uv,
                                                tf, deviation))) {
            oi(found ? 1 : 0);
            if (found) { od(tf); od(uv.x); od(uv.y); od(deviation); }
        }
    } else if (op == "global" || op == "grefuse") {
        double tolerance = rd();
        SmSolutionArray solutions;
        SmStatus st = surface->GlobalCurveIntersect(domain, *curve, interval, tolerance, solutions);
        if (op == "grefuse") {
            // The port refuses this case; SMLib's answer stays in the raw evidence.
            std::fprintf(stderr, "status %ld count %lu\n", (long)st, (unsigned long)solutions.GetSize());
            for (ULONG i = 0; i < solutions.GetSize(); ++i)
                std::fprintf(stderr, "  type %d t %.17g u %.17g v %.17g value %.17g\n",
                             (int)solutions[i].m_eSolutionType, solutions[i].m_vStart[0], solutions[i].m_vStart[1],
                             solutions[i].m_vStart[2], solutions[i].m_vStart.m_dSolutionValue);
            oi(-1);
        } else if (status(st)) {
            oi((long)solutions.GetSize());
            for (ULONG i = 0; i < solutions.GetSize(); ++i) {
                oi((long)solutions[i].m_eSolutionType);
                for (int k = 0; k < 3; ++k) od(solutions[i].m_vStart.m_adParameters[k]);
                od(solutions[i].m_vStart.m_dSolutionValue);
                if (solutions[i].m_eSolutionType == SM_ST_RANGE_OF_VALUES) {
                    for (int k = 0; k < 3; ++k) od(solutions[i].m_vEnd.m_adParameters[k]);
                    od(solutions[i].m_vEnd.m_dSolutionValue);
                }
            }
        }
    } else if (op == "coincide") {
        double tolerance = rd();
        long partial = ri();
        SmBoolean found = FALSE;
        SmSolution solution;
        if (status(surface->SimpleCoincidenceChecker(domain, *curve, interval, tolerance, found, solution,
                                                     partial ? TRUE : FALSE))) {
            oi(found ? 1 : 0);
            if (found) {
                for (int k = 0; k < 3; ++k) od(solution.m_vStart.m_adParameters[k]);
                od(solution.m_vStart.m_dSolutionValue);
                for (int k = 0; k < 3; ++k) od(solution.m_vEnd.m_adParameters[k]);
                od(solution.m_vEnd.m_dSolutionValue);
            }
        }
    } else if (op == "solve") {
        long operation = ri();
        double tolerance = rd(), t = rd(), gu = rd(), gv = rd();
        SmBoolean found = FALSE;
        SmSolution solution;
        if (status(surface->LocalCurveSolve(domain, *curve, interval, (SmSolverOperationType)operation, tolerance,
                                            NULL, NULL, SmPoint2d(gu, gv), t, found, solution))) {
            oi(found ? 1 : 0);
            if (found) {
                od(solution.m_vStart.m_adParameters[0]); od(solution.m_vStart.m_adParameters[1]);
                od(solution.m_vStart.m_adParameters[2]); od(solution.m_vStart.m_dSolutionValue);
            }
        }
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
