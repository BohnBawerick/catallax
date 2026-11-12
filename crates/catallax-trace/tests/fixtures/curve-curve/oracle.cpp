// Test-only oracle driver for the catallax-trace curve-curve gate. Reads one
// request line on stdin and prints one response line in the token format of
// ../foundation/oracle.cpp. A curve is `degree count rational knots... (x y z
// w)...` with Euclidean coordinates and weights. Catallax never compiles or
// links this file.
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

static void solutions(SmSolutionArray &s) {
    oi((long)s.GetSize());
    for (ULONG i = 0; i < s.GetSize(); ++i) {
        oi((long)s[i].m_eSolutionType);
        od(s[i].m_vStart.m_adParameters[0]); od(s[i].m_vStart.m_adParameters[1]); od(s[i].m_vStart.m_dSolutionValue);
        if (s[i].m_eSolutionType == SM_ST_RANGE_OF_VALUES) {
            od(s[i].m_vEnd.m_adParameters[0]); od(s[i].m_vEnd.m_adParameters[1]); od(s[i].m_vEnd.m_dSolutionValue);
        }
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
    SmBSplineCurve *a = rcurve(context);
    double a0 = rd(), a1 = rd();
    SmBSplineCurve *b = rcurve(context);
    double b0 = rd(), b1 = rd();
    SmExtent1d ia(a0, a1), ib(b0, b1);
    double tolerance = rd();
    long flag = ri();
    SmSolutionArray found;
    if (op == "global" || op == "grefuse") {
        // SmCurve::GlobalCurveIntersect, flag = skip the coincidence check.
        SmStatus st = a->GlobalCurveIntersect(ia, *b, ib, tolerance, found, flag ? TRUE : FALSE);
        if (op == "grefuse") {
            // The port refuses this case; SMLib's answer stays in the raw evidence.
            std::fprintf(stderr, "status %ld count %lu\n", (long)st, (unsigned long)found.GetSize());
            oi(-1);
        } else if (status(st)) {
            solutions(found);
        }
    } else if (op == "coincide") {
        // SmCurve::GlobalCoincidenceChecker, flag = snap to ends.
        SmBoolean more = FALSE;
        if (status(a->GlobalCoincidenceChecker(ia, *b, ib, tolerance, more, found, flag ? TRUE : FALSE))) {
            oi(more ? 1 : 0);
            solutions(found);
        }
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
