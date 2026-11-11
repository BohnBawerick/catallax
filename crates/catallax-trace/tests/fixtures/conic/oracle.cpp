// Test-only oracle driver for the catallax-trace conic gate. Reads one
// request line on stdin and prints one response line in the token format of
// ../foundation/oracle.cpp. Catallax never compiles or links this file.
#include "StdAfx.h"
#include <SmPolynomial.h>
#include <SmConic.h>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>

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

// `kind a b count (angle x y)...`: a conic moved by each transform in turn.
static SmConic rconic() {
    long kind = ri();
    double a = rd(), b = rd();
    SmConic conic((SmConicType)kind, a, b);
    long count = ri();
    for (long i = 0; i < count; ++i) {
        double angle = rd(), x = rd(), y = rd();
        conic.Transform2D(angle, SmPoint2d(x, y));
    }
    return conic;
}

int main() {
    std::string line;
    std::getline(std::cin, line);
    in.str(line);
    std::string op;
    in >> op;
    std::printf("%s", op.c_str());
    if (op == "cubic") {
        double c[4];
        for (double &x : c) x = rd();
        double zero = rd();
        ULONG n = 0;
        double s[3] = {0, 0, 0};
        oi((long)SmPolynomial::SolveCubicEqn(c, zero, n, s));
        oi((long)n);
        for (ULONG i = 0; i < n; ++i) od(s[i]);
    } else if (op == "quartic") {
        double c[5];
        for (double &x : c) x = rd();
        double zero = rd();
        ULONG n = 0;
        double s[4] = {0, 0, 0, 0};
        oi((long)SmPolynomial::SolveQuarticEqn(c, zero, n, s));
        oi((long)n);
        for (ULONG i = 0; i < n; ++i) od(s[i]);
    } else if (op == "conic" || op == "conicx") {
        SmConic one = rconic();
        SmConic two = rconic();
        double tolerance = rd();
        // Points along the first conic and their inversion on it.
        for (double t : {-0.75, -0.25, 0.3, 0.9, 1.5, 2.6}) {
            SmPoint3d p;
            one.EvaluatePoint(t, p);
            od(p.x); od(p.y);
            SmBoolean ok = FALSE;
            double back = 0.0;
            one.InvertPointOnConic(p, ok, back);
            oi(ok ? 1 : 0);
            if (ok) od(back);
        }
        double coeff[6];
        two.ComputeImplicitEquation(coeff);
        for (double c : coeff) od(c);
        SmSolutionArray solutions;
        oi((long)one.IntersectConic(one.GetNaturalInterval(), two, two.GetNaturalInterval(), tolerance, solutions));
        oi((long)solutions.GetSize());
        for (ULONG i = 0; i < solutions.GetSize(); ++i) {
            od(solutions[i].m_vStart[0]);
            od(solutions[i].m_vStart[1]);
            // The residual distance is rounding noise no caller reads; only
            // its bound (the tolerance) is checked.
            oi(solutions[i].m_vStart.m_dSolutionValue <= tolerance ? 1 : 0);
        }
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
