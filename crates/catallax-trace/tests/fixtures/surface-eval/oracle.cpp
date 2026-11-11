// Test-only oracle driver for the catallax-trace surface evaluation gate.
// Reads one request line on stdin and prints one response line in the token
// format of ../foundation/oracle.cpp. A surface is `du dv nu nv rational
// knots_u... knots_v... (x y z w)...` with every knot listed, Euclidean
// coordinates, weights, and v varying fastest. Catallax never compiles or
// links this file.
#include "StdAfx.h"
#include <SmBSplineSurface.h>
#include <SmContext.h>
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
    double u = rd(), v = rd();
    SmPoint2d uv(u, v);
    if (op == "point") {
        SmVector3d point;
        if (status(surface->EvaluatePoint(uv, point))) o3(point);
    } else if (op == "normal") {
        long left_u = ri(), left_v = ri();
        SmVector3d normal;
        if (status(surface->EvaluateNormal(uv, left_u, left_v, normal))) o3(normal);
    } else if (op == "eval" || op == "simple") {
        long hu = ri(), hv = ri(), left_u = ri(), left_v = ri(), upper = ri(), nonzero = ri();
        std::vector<SmVector3d> cells((hu + 1) * (hv + 1), SmVector3d(0, 0, 0));
        SmStatus st = op == "eval"
            ? surface->Evaluate(uv, hu, hv, left_u, left_v, upper, cells.data(), nonzero)
            : surface->EvaluateSimple(uv, hu, hv, left_u, left_v, upper, cells.data());
        if (status(st))
            for (long k = 0; k <= hu; ++k)
                for (long l = 0; l <= hv; ++l)
                    if (!upper || k + l <= hu) o3(cells[k * (hv + 1) + l]);
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
