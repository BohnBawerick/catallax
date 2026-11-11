// Test-only oracle driver for the catallax-trace surface point drop gate.
// Reads one request line on stdin and prints one response line in the token
// format of ../foundation/oracle.cpp. A surface is `s du dv nu nv rational
// knots_u... knots_v... (x y z w)...` with every knot listed, Euclidean
// coordinates, weights, and v varying fastest, or `p` and a plane as in
// ../plane/oracle.cpp. Catallax never compiles or links this file.
#include "StdAfx.h"
#include <SmBSplineSurface.h>
#include <SmPlane.h>
#include <SmAxis2Placement.h>
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


static SmSurface *rany(SmContext &context) {
    std::string kind;
    in >> kind;
    if (kind == "p") return rplane(context);
    return rsurface(context);
}

int main() {
    SmContext context;
    std::string line;
    std::getline(std::cin, line);
    in.str(line);
    std::string op;
    in >> op;
    std::printf("%s", op.c_str());
    SmSurface *surface = rany(context);
    if (op == "drop") {
        SmExtent2d domain = rdomain();
        SmVector3d point = r3();
        long has_guess = ri();
        double gu = rd(), gv = rd();
        long operation = ri();
        SmPoint2d guess(gu, gv);
        SmBoolean success = FALSE, multi = FALSE;
        SmPoint2d uv(0, 0);
        double gap = 0;
        if (status(surface->DropPoint(point, domain, has_guess ? &guess : NULL, success, uv, gap, multi,
                                      (SmSolverOperationType)operation))) {
            oi(success ? 1 : 0);
            if (success) { od(uv.x); od(uv.y); od(gap); oi(multi ? 1 : 0); }
        }
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
