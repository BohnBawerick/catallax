// Test-only oracle driver for the catallax-trace foundation gate.
// Reads one request line on stdin: an operation name and its inputs as
// 16-digit hexadecimal binary64 bit patterns or decimal integers. Prints one
// response line: tokens `d<hex>` for binary64 values and `i<decimal>` for
// integers, statuses, booleans and enumerations. After a failing status
// only the status is printed. Catallax never compiles or links this file.
#include "StdAfx.h"
#include <SmContext.h>
#include <SmExtent1d.h>
#include <SmExtent2d.h>
#include <SmExtent3d.h>
#include <SmGeomUtility.h>
#include <SmVector2d.h>
#include <SmVector3d.h>
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
static SmVector3d r3() { double x = rd(), y = rd(), z = rd(); return SmVector3d(x, y, z); }
static SmVector2d r2() { double x = rd(), y = rd(); return SmVector2d(x, y); }

static void od(double value) {
    uint64_t bits;
    std::memcpy(&bits, &value, sizeof bits);
    std::printf(" d%016" PRIx64, bits);
}
static void oi(long value) { std::printf(" i%ld", value); }
static void o3(const SmVector3d &v) { od(v.x); od(v.y); od(v.z); }
static void o2(const SmVector2d &v) { od(v.x); od(v.y); }
static void ob2(const SmExtent2d &e) { o2(e.GetMin()); o2(e.GetMax()); }
static void ob3(const SmExtent3d &e) { o3(e.GetMin()); o3(e.GetMax()); }
static bool status(SmStatus s) { oi((long)s); return s == SM_SUCCESS; }

int main() {
    SmContext context;
    std::string line;
    std::getline(std::cin, line);
    in.str(line);
    std::string op;
    in >> op;
    std::printf("%s", op.c_str());
    if (op == "angle") {
        SmVector3d a = r3(), b = r3();
        double angle = 0.0;
        if (status(a.AngleBetween(b, angle))) od(angle);
    } else if (op == "unitize") {
        SmVector3d a = r3();
        SmStatus s = a.Unitize(FALSE);
        oi((long)s); o3(a);
    } else if (op == "parallel") {
        SmVector3d a = r3(), b = r3();
        double degrees = rd();
        SmBoolean anti = 2;
        oi(a.IsParallelTo(b, degrees, &anti)); oi(anti);
    } else if (op == "perp") {
        SmVector3d a = r3(), b = r3();
        oi(a.IsPerpendicularTo(b, rd()));
    } else if (op == "ortho") {
        SmVector3d a = r3();
        long has = ri();
        SmVector3d reference = r3(), x, y, z;
        if (status(a.MakeUnitOrthoVectors(has ? &reference : NULL, x, y, z))) { o3(x); o3(y); o3(z); }
    } else if (op == "uderiv") {
        SmVector3d a = r3(), d = r3();
        o3(a.UnitizedDerivative(d));
    } else if (op == "divassign3") {
        SmVector3d a = r3();
        a /= rd();
        o3(a);
    } else if (op == "div2") {
        SmVector2d a = r2();
        o2(a / rd());
    } else if (op == "closer") {
        SmVector3d a = r3(), b = r3();
        oi(a.CloserThan(rd(), b));
    } else if (op == "differs2") {
        SmVector2d a = r2(), b = r2();
        oi(a != b);
    } else if (op == "ext1") {
        double lo = rd(), hi = rd(), x = rd(), tol = rd(), olo = rd(), ohi = rd(), u = rd();
        SmExtent1d e(lo, hi), o(olo, ohi);
        od(e.ClampValue(x)); od(e.Evaluate(u)); oi((long)e.ClassifyPoint(x, tol));
        double other = 0.0;
        oi(e.IsValueOnBoundary(x, tol, &other)); od(other);
        oi(e.AreDisjoint(o, tol)); oi(e.IsContainedBy(o, tol)); oi(e.ContainsValue(x, tol));
        oi(e.AreEqual(o, tol)); oi(e.AreEqual(olo, ohi, tol)); oi(e.IsBounded());
        od(e.GetMid()); od(e.GetLength());
        SmExtent1d grown = e;
        grown.ExpandAbsolute(tol);
        od(grown.GetMin()); od(grown.GetMax());
        SmExtent1d set;
        oi((long)set.SetMinMax(olo, ohi));
    } else if (op == "ext2") {
        SmVector2d lo = r2(), hi = r2(), p = r2();
        double tol = rd();
        SmVector2d olo = r2(), ohi = r2();
        double u = rd(), v = rd();
        SmExtent2d e(lo, hi), o(olo, ohi);
        oi(e.ContainsPoint2d(p, tol));
        double distance = -1.0;
        oi(e.ContainsPoint2d(p, tol, &distance)); od(distance);
        SmVector2d normal(7.0, 7.0);
        oi(e.IsPoint2dOnBoundary(p, tol, &normal)); o2(normal);
        oi((long)e.GetPoint2dBoundaries(p, tol));
        o2(e.ClampPoint2d(p)); o2(e.Evaluate(u, v));
        SmExtent2d joined;
        e.Union(o, joined); ob2(joined);
        oi(e.AreDisjoint(o, tol)); oi(e.IsContainedBy(o, tol));
        o2(e.GetSize()); od(e.XLength()); od(e.YLength());
        od(e.GetUInterval().GetMin()); od(e.GetUInterval().GetMax());
        od(e.GetVInterval().GetMin()); od(e.GetVInterval().GetMax());
        SmExtent2d common;
        if (status(e.Intersect(o, common))) ob2(common);
    } else if (op == "clip") {
        SmVector2d lo = r2(), hi = r2(), point = r2(), vector = r2();
        SmExtent2d e(lo, hi);
        od(e.ClipLine2d(point, vector, rd()));
    } else if (op == "ext3") {
        long raw = ri();
        SmVector3d lo = r3(), hi = r3(), p = r3();
        double tol = rd();
        SmVector3d olo = r3(), ohi = r3(), n = r3();
        SmExtent3d e = raw ? SmExtent3d(lo.x, lo.y, lo.z, hi.x, hi.y, hi.z) : SmExtent3d(lo, hi);
        SmExtent3d o(olo, ohi);
        oi(e.ContainsPoint3d(p, tol)); oi(e.AreDisjoint(o, tol));
        SmExtent3d joined;
        e.Union(o, joined); ob3(joined);
        oi(e.IsPointSized(tol)); oi(e.HasNegativeVolume());
        SmVector3d center;
        double radius = 0.0, inner = 0.0;
        e.ComputeSphereBound(center, radius, &inner); o3(center); od(radius); od(inner);
        SmExtent3d grown = e;
        grown.ExpandAbsolute(tol); ob3(grown);
        o3(e.Evaluate(n.x, n.y, n.z)); o3(e.GetSize()); od(e.GetMaxDimension());
        SmVector3d axis;
        double angle = 0.0;
        if (status(e.ComputeConeBoundTo(o, axis, angle))) { o3(axis); od(angle); }
        SmExtent3d common;
        if (status(e.Intersect(o, common))) ob3(common);
    } else if (op == "triarea") {
        SmVector3d a = r3(), b = r3(), c = r3();
        double area = 0.0;
        if (status(smgu_TriangleAreaSquared(a, b, c, area))) od(area);
    } else if (op == "tripoint") {
        SmVector3d a = r3(), b = r3(), c = r3(), t = r3(), closest;
        double distance = 0.0, u = 0.0, v = 0.0;
        if (status(smgu_TrianglePointDistance(a, b, c, t, distance, closest, u, v))) {
            od(distance); o3(closest); od(u); od(v);
        }
    } else if (op == "linepoint") {
        SmVector3d p = r3(), d = r3(), t = r3();
        double distance = 0.0, parameter = 0.0;
        if (status(smgu_LinePointDistance(p, d, t, distance, &parameter))) { od(distance); od(parameter); }
    } else if (op == "lineclosest") {
        SmVector3d p = r3(), d = r3(), t = r3();
        double parameter = 0.0;
        if (status(smgu_LineClosestPoint(p, d, t, parameter))) od(parameter);
    } else if (op == "segpoint") {
        SmVector3d s = r3(), e = r3(), t = r3();
        double distance = 0.0, parameter = 0.0;
        if (status(smgu_SegmentPointDistance(s, e, t, distance, parameter))) { od(distance); od(parameter); }
    } else if (op == "segangle") {
        SmVector3d s = r3(), e = r3(), t = r3();
        double parameter = 0.0;
        double angle = smgu_SegmentPointAngleDegToPerp(s, e, t, parameter);
        od(angle); od(parameter);
    } else if (op == "linecircle") {
        SmVector3d p = r3(), d = r3(), c = r3(), n = r3();
        double radius = rd(), tol = rd();
        ULONG count = 0;
        double parameters[2] = {0.0, 0.0};
        if (status(smgu_LineCoPlanarCircleIntersect(p, d, c, n, radius, tol, count, parameters))) {
            oi((long)count);
            for (ULONG i = 0; i < count; ++i) od(parameters[i]);
        }
    } else if (op == "twoplanes") {
        SmVector3d p1 = r3(), n1 = r3(), p2 = r3(), n2 = r3(), point, direction;
        SmBoolean coincident = 2;
        SmStatus s = smgu_IntersectTwoPlanes(p1, n1, p2, n2, point, direction, &coincident);
        oi((long)s); oi(coincident);
        if (s == SM_SUCCESS) { o3(point); o3(direction); }
    } else if (op == "threeplanes") {
        SmVector3d p1 = r3(), n1 = r3(), p2 = r3(), n2 = r3(), p3 = r3(), n3 = r3(), point;
        if (status(smgu_IntersectThreePlanes(p1, n1, p2, n2, p3, n3, point))) o3(point);
    } else if (op == "planedist") {
        SmVector3d p = r3(), n = r3(), t = r3();
        double distance = 0.0;
        if (status(smgu_PlanePointDistance(p, n, t, distance))) od(distance);
    } else if (op == "ptproj") {
        SmVector3d p = r3(), q = r3(), n = r3(), out;
        if (status(smgu_PointProjectToPlane(p, q, n, out))) o3(out);
    } else if (op == "vecproj") {
        SmVector3d v = r3(), n = r3(), out;
        if (status(smgu_VectorProjectToPlane(v, n, out))) o3(out);
    } else if (op == "lincomb") {
        SmVector3d a = r3(), b = r3(), v = r3();
        SmVector2d uv;
        if (status(smgu_VecLinCombTwoVectors(a, b, v, uv))) o2(uv);
    } else if (op == "det3") {
        SmVector3d a = r3(), b = r3(), c = r3();
        od(smgu_Determinant3Vectors(a, b, c));
    } else if (op == "frenet") {
        long count = ri();
        SmVector3d values[4] = {r3(), r3(), r3(), r3()};
        SmVector3d frame[4];
        for (auto &f : frame) f.Set(5.0, 5.0, 5.0);
        oi((long)smgu_FrenetFromDerivatives(values, (ULONG)count, frame));
        for (auto &f : frame) o3(f);
    } else if (op == "continuity") {
        SmVector3d p1 = r3(), d1 = r3(), s1 = r3(), p2 = r3(), d2 = r3(), s2 = r3();
        double degrees = rd();
        long has = ri();
        SmVector3d t1 = r3(), t2 = r3();
        SmContinuityType type = SM_CT_UNDEFINED;
        if (status(smgu_EvaluateCurveContinuity(p1, d1, s1, p2, d2, s2, type, degrees,
                                               has ? &t1 : NULL, has ? &t2 : NULL)))
            oi((long)type);
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
