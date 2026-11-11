// Test-only oracle driver for the catallax-trace surface point solver gate.
// Reads one request line on stdin and prints one response line in the token
// format of ../foundation/oracle.cpp. A surface is `du dv nu nv rational
// knots_u... knots_v... (x y z w)...` with every knot listed, Euclidean
// coordinates, weights, and v varying fastest. Catallax never compiles or
// links this file.
#include "StdAfx.h"
#include <SmBSplineSurface.h>
#include <SmContext.h>
#include <SmSolutionArray.h>
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

int main() {
    SmContext context;
    std::string line;
    std::getline(std::cin, line);
    in.str(line);
    std::string op;
    in >> op;
    std::printf("%s", op.c_str());
    SmBSplineSurface *surface = rsurface(context);
    double u0 = rd(), v0 = rd(), u1 = rd(), v1 = rd();
    SmExtent2d domain(SmPoint2d(u0, v0), SmPoint2d(u1, v1));
    long operation = ri();
    double x = rd(), y = rd(), z = rd();
    if (op == "global") {
        // SmSurface::GlobalPointSolve_1 with a distance tolerance.
        double tolerance = rd();
        SmSolutionArray solutions;
        SmStatus st = surface->GlobalPointSolve_1(domain, (SmSolverOperationType)operation, SmPoint3d(x, y, z),
                                                  tolerance, NULL, SM_SR_ALL, solutions);
        if (status(st)) {
            oi((long)solutions.GetSize());
            for (ULONG i = 0; i < solutions.GetSize(); ++i) {
                od(solutions[i].m_vStart.m_adParameters[0]);
                od(solutions[i].m_vStart.m_adParameters[1]);
                od(solutions[i].m_vStart.m_dSolutionValue);
            }
        }
        std::printf("\n");
        return 0;
    }
    double gu = rd(), gv = rd();
    if (op == "solve" || op == "refuse") {
        SmBoolean found = FALSE;
        SmSolution solution;
        SmTerminationReasonType reason = SM_TR_UNDEFINED;
        int iterations = -1;
        SmStatus st = surface->LocalPointSolve_1(domain, (SmSolverOperationType)operation, SmPoint3d(x, y, z),
                                                 SmPoint2d(gu, gv), found, solution, reason, iterations);
        if (op == "refuse") {
            // The port refuses this case; SMLib's answer stays in the raw evidence.
            std::fprintf(stderr, "status %ld found %d reason %d iterations %d\n", (long)st, (int)found,
                         (int)reason, iterations);
            oi(-1);
        } else if (status(st)) {
            oi(found ? 1 : 0);
            od(solution.m_vStart.m_adParameters[0]);
            od(solution.m_vStart.m_adParameters[1]);
            od(solution.m_vStart.m_dSolutionValue);
            oi((long)reason);
            oi(iterations);
        }
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
