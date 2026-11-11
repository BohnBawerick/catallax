// Test-only oracle driver for the catallax-trace plane gate.
// Reads one request line on stdin and prints one response line in the token
// format of ../foundation/oracle.cpp. A plane is `kind origin x_axis y_axis
// u_scale v_scale domain`, kind 1
// for SmPlane::CreateCanonical (which ignores scale and domain) and 0 for the
// bounded constructor followed by MakeNurb; a curve is `degree count rational
// knots... (x y z w)...`. Catallax never compiles or links this file.
#include "StdAfx.h"
#include <SmBSplineSurface.h>
#include <SmPlane.h>
#include <SmAxis2Placement.h>
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

static SmVector3d r3() { double x = rd(), y = rd(), z = rd(); return SmVector3d(x, y, z); }
static SmExtent2d rdomain() { double a = rd(), b = rd(), c = rd(), d = rd(); return SmExtent2d(SmPoint2d(a, b), SmPoint2d(c, d)); }

static SmPlane *rplane(SmContext &context) {
    long kind = ri();
    SmVector3d origin = r3(), x = r3(), y = r3();
    double su = rd(), sv = rd();
    SmExtent2d domain = rdomain();
    SmPlane *p = NULL;
    if (kind == 1) {
        SmAxis2Placement frame;
        frame.SetCanonical(origin, x, y);
        if (SmPlane::CreateCanonical(context, frame, p) != SM_SUCCESS || !p) { std::fprintf(stderr, "plane failed\n"); std::exit(4); }
    } else {
        p = new (context) SmPlane(origin, x, y, SmVector2d(su, sv), domain, &context);
        if (p->MakeNurb() != SM_SUCCESS) { std::fprintf(stderr, "make nurb failed\n"); std::exit(4); }
    }
    return p;
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

static void solutions(SmSolutionArray &s, int variables) {
    oi((long)s.GetSize());
    for (ULONG i = 0; i < s.GetSize(); ++i) {
        oi((long)s[i].m_eSolutionType);
        for (int k = 0; k < variables; ++k) od(s[i].m_vStart.m_adParameters[k]);
        od(s[i].m_vStart.m_dSolutionValue);
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
    SmPlane *plane = rplane(context);
    if (op == "net") {
        SmExtent2d d = plane->GetNaturalUVDomain();
        oi(1000); od(d.GetMin().x); od(d.GetMin().y); od(d.GetMax().x); od(d.GetMax().y);
        for (ULONG i = 0; i < 2; ++i)
            for (ULONG j = 0; j < 2; ++j) {
                SmPoint3d c; double w = 0;
                if (status(plane->GetControlPoint(SM_CP_NON_RATIONAL, i, j, c, w))) o3(c);
            }
    } else if (op == "point") {
        double u = rd(), v = rd();
        SmPoint3d pt;
        if (status(plane->EvaluatePoint(SmPoint2d(u, v), pt))) o3(pt);
    } else if (op == "eval") {
        double u = rd(), v = rd();
        long hu = ri(), hv = ri(), upper = ri(), simple = ri();
        std::vector<SmVector3d> cells((hu + 1) * (hv + 1));
        SmStatus st = simple ? plane->EvaluateSimple(SmPoint2d(u, v), hu, hv, TRUE, TRUE, upper, cells.data())
                             : plane->Evaluate(SmPoint2d(u, v), hu, hv, TRUE, TRUE, upper, cells.data());
        if (status(st))
            for (long k = 0; k <= hu; ++k)
                for (long l = 0; l <= hv; ++l)
                    if (!upper || k + l <= hu) o3(cells[k * (hv + 1) + l]);
    } else if (op == "normal") {
        double u = rd(), v = rd();
        SmVector3d n;
        if (status(plane->EvaluateNormal(SmPoint2d(u, v), TRUE, TRUE, n))) o3(n);
    } else if (op == "flags") {
        double u = rd(), v = rd();
        SmSurfParamType dir = SM_SP_BOTH;
        oi(1000);
        oi(plane->IsSingularity(SmPoint2d(u, v), dir, 1e-5, FALSE) ? 1 : 0); oi((long)dir);
        oi(plane->IsClosed(plane->GetNaturalUVDomain(), SM_SP_U) ? 1 : 0);
        oi(plane->IsOnSeam(SmPoint2d(u, v)) ? 1 : 0);
    } else if (op == "drop") {
        long operation = ri();
        SmExtent2d domain = rdomain();
        SmVector3d test = r3();
        double tolerance = rd();
        long has_target = ri();
        double target = rd();
        SmSolutionArray found;
        if (status(plane->DropPointFast(domain, (SmSolverOperationType)operation, test, NULL, tolerance,
                                        has_target ? &target : NULL, SM_SR_SINGLE, found)))
            solutions(found, 2);
    } else if (op == "global") {
        long operation = ri();
        SmExtent2d domain = rdomain();
        SmVector3d test = r3();
        double tolerance = rd();
        SmSolutionArray found;
        if (status(plane->GlobalPointSolve(domain, (SmSolverOperationType)operation, test, tolerance, NULL,
                                           SM_SR_SINGLE, found)))
            solutions(found, 2);
    } else if (op == "cs") {
        SmBSplineCurve *curve = rcurve(context);
        SmExtent2d domain = rdomain();
        double t0 = rd(), t1 = rd(), tolerance = rd();
        SmSolutionArray found;
        if (status(plane->GlobalCurveIntersect(domain, *curve, SmExtent1d(t0, t1), tolerance, found)))
            solutions(found, 3);
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
