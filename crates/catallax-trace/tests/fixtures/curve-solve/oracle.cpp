// Test-only oracle driver for the catallax-trace curve solver gate.
// Reads one request line on stdin and prints one response line in the token
// format of ../foundation/oracle.cpp. A curve is `degree count rational
// knots... (x y z w)...` with Euclidean coordinates and weights. Catallax
// never compiles or links this file.
#include "StdAfx.h"
#include <SmBSplineCurve.h>
#include <nurbs.h>
#include <SmContext.h>
#include <SmCurve.h>
#include <SmPseudoBox.h>
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
static SmContext *g_context;

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
static SmExtent1d rx() { double a = rd(), b = rd(); return SmExtent1d(a, b); }
static void od(double value) {
    uint64_t bits;
    std::memcpy(&bits, &value, sizeof bits);
    std::printf(" d%016" PRIx64, bits);
}
static void oi(long value) { std::printf(" i%ld", value); }
static bool status(SmStatus s) { oi((long)s); return s == SM_SUCCESS; }

static SmBSplineCurve *rcurve() {
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
    // Objects are freed with the process.
    return new (*g_context) SmBSplineCurve(3, count - 1, count + degree, (short)degree, points.data(), knots.data(),
                                           g_context);
}

static void solution(const SmSolution &s) {
    oi((long)s.m_eSolutionType);
    od(s.m_vStart.m_adParameters[0]); od(s.m_vStart.m_adParameters[1]); od(s.m_vStart.m_dSolutionValue);
    if (s.m_eSolutionType == SM_ST_RANGE_OF_VALUES) {
        od(s.m_vEnd.m_adParameters[0]); od(s.m_vEnd.m_adParameters[1]); od(s.m_vEnd.m_dSolutionValue);
    }
}

int main() {
    SmContext context;
    g_context = &context;
    std::string line;
    std::getline(std::cin, line);
    in.str(line);
    std::string op;
    in >> op;
    std::printf("%s", op.c_str());
    if (op == "point") {
        auto curve = rcurve();
        SmExtent1d interval = rx();
        long operation = ri();
        SmVector3d test = r3();
        double guess = rd();
        SmBoolean found = FALSE;
        SmSolution s;
        if (status(curve->LocalPointSolve(interval, (SmSolverOperationType)operation, test, NULL, NULL, NULL, guess,
                                          found, s))) {
            oi(found);
            if (found) { od(s.m_vStart.m_adParameters[0]); od(s.m_vStart.m_dSolutionValue); }
        }
    } else if (op == "intersect") {
        auto a = rcurve();
        auto b = rcurve();
        SmExtent1d ia = rx(), ib = rx();
        double tol = rd(), ga = rd(), gb = rd();
        SmBoolean found = FALSE;
        double ta = 0.0, tb = 0.0, dev = 0.0;
        if (status(a->LocalCurveIntersect(ia, *b, ib, tol, ga, gb, found, ta, tb, dev))) {
            oi(found);
            if (found) { od(ta); od(tb); od(dev); }
        }
    } else if (op == "solve") {
        auto a = rcurve();
        auto b = rcurve();
        SmExtent1d ia = rx(), ib = rx();
        long operation = ri();
        double tol = rd(), ga = rd(), gb = rd();
        SmBoolean found = FALSE;
        SmSolution s;
        if (status(a->LocalCurveSolve(ia, *b, ib, (SmSolverOperationType)operation, tol, NULL, NULL, ga, gb, found,
                                      s))) {
            oi(found);
            if (found) { od(s.m_vStart.m_adParameters[0]); od(s.m_vStart.m_adParameters[1]); od(s.m_vStart.m_dSolutionValue); }
        }
    } else if (op == "maxdist") {
        auto a = rcurve();
        auto b = rcurve();
        SmExtent1d ia = rx();
        double start = rd(), end = rd();
        long samples = ri(), has_needed = ri();
        double needed = rd();
        double distance = -1.0, ta = -1.0, tb = -1.0;
        SmBoolean wrong = 2;
        if (status(a->CurveMaxDistanceBetween(ia, *b, start, end, (ULONG)samples, has_needed ? &needed : NULL,
                                              distance, &ta, &tb, &wrong))) {
            od(distance); od(ta); od(tb); oi(wrong);
        }
    } else if (op == "coincide") {
        auto a = rcurve();
        auto b = rcurve();
        SmExtent1d ia = rx(), ib = rx();
        double tol = rd();
        long boxes = ri();
        SmPseudoBox box_a, box_b;
        if (boxes) {
            a->CalculateBoundingBox(ia, NULL, &box_a, NULL, FALSE);
            b->CalculateBoundingBox(ib, NULL, &box_b, NULL, FALSE);
        }
        ULONG count = 0;
        SmSolution s[2];
        if (status(a->SimpleCoincidenceChecker(ia, *b, ib, tol, boxes ? &box_a : NULL, boxes ? &box_b : NULL, count,
                                               s))) {
            oi((long)count);
            for (ULONG i = 0; i < count; ++i) solution(s[i]);
        }
    } else if (op == "closed") {
        auto a = rcurve();
        SmExtent1d ia = rx();
        double tol = rd();
        oi(a->IsClosed(ia, tol)); oi(a->IsPeriodic(ia));
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
