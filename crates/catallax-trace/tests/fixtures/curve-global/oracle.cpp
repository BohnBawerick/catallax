// Test-only oracle driver for the catallax-trace global point-curve gate.
// Reads one request line on stdin and prints one response line in the token
// format of ../foundation/oracle.cpp. A curve is `degree count rational
// knots... (x y z w)...` with Euclidean coordinates and weights. Catallax
// never compiles or links this file.
#include "StdAfx.h"
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
static bool status(SmStatus s) { oi((long)s); return s == SM_SUCCESS; }

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
    SmBSplineCurve *curve = rcurve(context);
    double t0 = rd(), t1 = rd();
    SmExtent1d interval(t0, t1);
    long operation = ri();
    double x = rd(), y = rd(), z = rd(), tolerance = rd();
    SmPoint3d test(x, y, z);
    if (op == "global" || op == "grefuse") {
        long requested = ri();
        SmSolutionArray solutions;
        SmStatus st = curve->GlobalPointSolve(interval, (SmSolverOperationType)operation, test, tolerance, NULL, NULL,
                                              (SmSolutionRequestedType)requested, solutions);
        if (op == "grefuse") {
            // The port refuses this case; SMLib's answer stays in the raw evidence.
            std::fprintf(stderr, "status %ld count %lu\n", (long)st, (unsigned long)solutions.GetSize());
            for (ULONG i = 0; i < solutions.GetSize(); ++i)
                std::fprintf(stderr, "  type %d t %.17g value %.17g\n", (int)solutions[i].m_eSolutionType,
                             solutions[i].m_vStart[0], solutions[i].m_vStart.m_dSolutionValue);
            oi(-1);
        } else if (status(st)) {
            oi((long)solutions.GetSize());
            for (ULONG i = 0; i < solutions.GetSize(); ++i) {
                oi((long)solutions[i].m_eSolutionType);
                od(solutions[i].m_vStart.m_adParameters[0]);
                od(solutions[i].m_vStart.m_dSolutionValue);
                if (solutions[i].m_eSolutionType == SM_ST_RANGE_OF_VALUES) {
                    od(solutions[i].m_vEnd.m_adParameters[0]);
                    od(solutions[i].m_vEnd.m_dSolutionValue);
                }
            }
        }
    } else if (op == "drop" || op == "drefuse") {
        long has_inward = ri();
        double ix = rd(), iy = rd(), iz = rd();
        long has_guess = ri();
        double guess = rd();
        SmVector3d inward(ix, iy, iz);
        SmBoolean success = FALSE;
        double parameter = 0, distance = 0;
        SmStatus st = curve->DropPoint(interval, test, has_inward ? &inward : NULL, tolerance,
                                       has_guess ? &guess : NULL, success, parameter, distance,
                                       (SmSolverOperationType)operation);
        if (op == "drefuse") {
            std::fprintf(stderr, "status %ld success %d t %.17g distance %.17g\n", (long)st, (int)success, parameter,
                         distance);
            oi(-1);
        } else if (status(st)) {
            oi(success ? 1 : 0);
            od(parameter);
            od(distance);
        }
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
