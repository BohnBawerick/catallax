// Test-only oracle driver for the catallax-trace output-curve gate.
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
#include <SmHermiteCurve.h>
#include <SmLine.h>
#include <SmPseudoBox.h>
#include <SmTol.h>
#include <SmVector2d.h>
#include <SmVector3d.h>
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
static SmVector3d r3() { double x = rd(), y = rd(), z = rd(); return SmVector3d(x, y, z); }
static SmVector2d r2() { double x = rd(), y = rd(); return SmVector2d(x, y); }

static void od(double value) {
    uint64_t bits;
    std::memcpy(&bits, &value, sizeof bits);
    std::printf(" d%016" PRIx64, bits);
}
static void oi(long value) { std::printf(" i%ld", value); }
static void o3(const SmVector3d &v) { od(v.x); od(v.y); od(v.z); }
static void ob3(const SmExtent3d &e) { o3(e.GetMin()); o3(e.GetMax()); }
static void opseudo(const SmPseudoBox &b) {
    for (ULONG i = 0; i < 3; ++i) o3(b.GetBasis(i));
    SmExtent1d a, c, d;
    b.GetIntervals(a, c, d);
    od(a.GetMin()); od(a.GetMax()); od(c.GetMin()); od(c.GetMax()); od(d.GetMin()); od(d.GetMax());
}
static bool status(SmStatus s) { oi((long)s); return s == SM_SUCCESS; }

#include <SmBSplineCurve.h>
#include <SmNurbsCrv.h>
#include <nurbs.h>

static void raw_curve(const gw_CURVE *c) {
    oi(c->p); oi(c->pol->n + 1); oi(c->knt->m + 1);
    for (long i = 0; i <= c->knt->m; ++i) od(c->knt->U[i]);
    for (long i = 0; i <= c->pol->n; ++i) {
        od(c->pol->Pw[i].x); od(c->pol->Pw[i].y); od(c->pol->Pw[i].z); od(c->pol->Pw[i].w);
    }
}

// `degree count knots... (x y z w)...`: the stored coordinates, weighted when
// the curve is rational, with `NL_NOZ` and `NL_NOW` markers as SMLib keeps them.
static gw_CURVE *read_curve() {
    long p = ri(), n = ri();
    gw_CURVE *raw = sm_AllocateNurbCurve(n - 1, p, n + p);
    for (long i = 0; i < n + p + 1; ++i) raw->knt->U[i] = rd();
    for (long i = 0; i < n; ++i) {
        raw->pol->Pw[i].x = rd(); raw->pol->Pw[i].y = rd(); raw->pol->Pw[i].z = rd(); raw->pol->Pw[i].w = rd();
    }
    return raw;
}

static void dump_line(SmLine &line, const SmContext &ctx) {
    double tol = rd(), t = rd();
    long der = ri();
    double step_t = rd(), ia = rd(), ib = rd();
    long reparam = ri();
    double rmin = rd(), rmax = rd();
    SmExtent1d iv = line.GetNaturalInterval();
    od(iv.GetMin()); od(iv.GetMax());
    oi(line.IsBounded());
    oi(line.IsDegenerate(tol, NULL));
    SmExtent1d part(ia, ib);
    oi(line.IsDegenerate(tol, &part));
    for (int pass = 0; pass < 2; ++pass) {
        if (pass == 1) {
            if (!reparam) break;
            SmExtent1d next(rmin, rmax);
            if (!status(line.EditParameterization(next))) return;
            SmExtent1d now = line.GetNaturalInterval();
            od(now.GetMin()); od(now.GetMax());
            t = now.Evaluate(0.3);
            step_t = now.Evaluate(0.7);
            part = SmExtent1d(now.Evaluate(0.2), now.Evaluate(0.9));
        }
        raw_curve(line.GetGwNurbPointer());
        double step;
        if (status(line.ConvertTFromNURBSToSTEP(t, step))) od(step);
        SmVector3d values[5];
        if (status(line.Evaluate(t, der, TRUE, values, TRUE))) for (long i = 0; i <= der; ++i) o3(values[i]);
        SmVector3d step_values[5];
        if (status(line.EvaluateSTEP(step_t, der, TRUE, step_values, 0.0))) for (long i = 0; i <= der; ++i) o3(step_values[i]);
        SmPoint3d point;
        if (status(line.EvaluatePoint(t, point))) o3(point);
        SmExtent3d box;
        SmPseudoBox pseudo;
        if (status(line.CalculateBoundingBox(part, &box, &pseudo, NULL, FALSE))) { ob3(box); opseudo(pseudo); }
    }
}

int main() {
    SmContext ctx;
    std::string line, op;
    std::getline(std::cin, line); in.str(line); in >> op;
    std::printf("%s", op.c_str());
    if (op == "consts") {
        od(SmTol::GetZoneTol3d());
    } else if (op == "canon") {
        long dim = ri(), degree = ri(), n = ri(), rational = ri(), has_trim = ri();
        SmExtent1d trim(0.0, 1.0);
        if (has_trim) { double a = rd(), b = rd(); trim = SmExtent1d(a, b); }
        long unique = ri();
        SmTArray<SmPoint3d> points;
        for (long i = 0; i < n; ++i) points.Add(r3());
        SmTArray<double> weights;
        if (rational) for (long i = 0; i < n; ++i) weights.Add(rd());
        SmTArray<double> knots;
        SmTArray<ULONG> mults;
        for (long i = 0; i < unique; ++i) { knots.Add(rd()); mults.Add((ULONG)ri()); }
        SmBSplineCurve *curve = NULL;
        SmStatus s = SmBSplineCurve::CreateCanonical(ctx, dim, degree, points, SM_CF_UNSPECIFIED, mults, knots,
                                                     SM_KT_UNSPECIFIED, rational ? &weights : NULL,
                                                     has_trim ? &trim : NULL, curve);
        if (status(s)) raw_curve(curve->GetGwNurbPointer());
    } else if (op == "reparam") {
        gw_CURVE *raw = read_curve();
        double a = rd(), b = rd();
        SmBSplineCurve c(raw, 3, TRUE, &ctx);
        SmExtent1d next(a, b);
        if (status(c.EditParameterization(next))) raw_curve(raw);
        sm_FreeNurbCurve(raw);
    } else if (op == "trim") {
        gw_CURVE *raw = read_curve();
        double a = rd(), b = rd();
        long mode = ri();
        SmBSplineCurve c(raw, 3, TRUE, &ctx);
        SmExtent1d got(a, b);
        SmTArray<double> unique_knots;
        sm_GetKnots(raw->knt, unique_knots, NULL, NULL);
        // Zone tolerance in parameter space at the requested ends, before any snapping.
        SmZoneTol3d zone = SmTol::GetZoneTol3d(&c);
        SmTol1d zone_min = SmTol::MapTo1d(zone, a, c);
        SmTol1d zone_max = SmTol::MapTo1d(zone, b, c);
        SmStatus s = c.Trim(got, mode ? UNSURE : TRUE);
        if (status(s)) {
            od(got.GetMin()); od(got.GetMax());
            // Whether SMLib's trim cut the curve and then entered its
            // degenerate-segment branch.
            long entered = 0;
            if (c.GetGwNurbPointer() != raw) {
                gw_CURVE *plain = NULL;
                if (sm_CrvTrim(raw, got.GetMin(), got.GetMax(), plain) == 0) {
                    double *U = plain->knt->U;
                    long m = plain->knt->m, p = plain->p;
                    if (((U[p + 1] - U[0]) < zone_min && unique_knots.GetSize() > 2) ||
                        ((U[m] - U[m - p - 1]) < zone_max && unique_knots.GetSize() > 2)) entered = 1;
                }
                sm_FreeNurbCurve(plain);
            }
            oi(entered);
            if (!mode) raw_curve(c.GetGwNurbPointer());
        }
        sm_FreeNurbCurve(raw);
    } else if (op == "line") {
        long mode = ri();
        if (mode == 0) {
            SmPoint3d p = r3();
            SmVector3d v = r3();
            double a = rd(), b = rd(), scale = rd();
            long dim = ri(), inside = ri();
            SmLine l(p, v, SmExtent1d(a, b), scale, dim, &ctx, NULL, inside);
            oi((long)SM_SUCCESS);
            dump_line(l, ctx);
        } else {
            SmPoint3d s = r3(), e = r3();
            long dim = ri(), has = ri();
            double a = rd(), b = rd();
            SmExtent1d ivl(a, b);
            SmLine *l = NULL;
            if (status(SmLine::CreateLineSegment(ctx, dim, s, e, l, has ? &ivl : NULL))) dump_line(*l, ctx);
        }
    } else if (op == "hermite") {
        long two_d = ri();
        SmVector3d p1 = r3(), d1 = r3(), p2 = r3(), d2 = r3();
        double t = rd();
        long der = ri(), nonzero = ri();
        SmHermiteCurve *h = two_d
            ? new (ctx) SmHermiteCurve(SmPoint2d(p1.x, p1.y), SmVector2d(d1.x, d1.y), SmPoint2d(p2.x, p2.y), SmVector2d(d2.x, d2.y), 2)
            : new (ctx) SmHermiteCurve(p1, d1, p2, d2, 3);
        SmExtent1d nat = h->GetNaturalInterval();
        od(nat.GetMin()); od(nat.GetMax());
        SmVector3d values[5];
        if (status(h->Evaluate(t, der, TRUE, values, nonzero))) for (long i = 0; i <= der; ++i) o3(values[i]);
        SmPoint3d point;
        if (status(h->EvaluatePoint(t, point))) o3(point);
        SmExtent3d box;
        SmPseudoBox pseudo;
        if (status(h->CalculateBoundingBox(nat, &box, &pseudo, NULL, FALSE))) { ob3(box); opseudo(pseudo); }
        // The process exits after one request, so the object is not freed.
    } else if (op == "clip") {
        SmVector2d lo = r2(), hi = r2(), point = r2(), vec = r2();
        double param = rd();
        SmExtent2d box(lo, hi);
        od(box.ClipLine2d(point, vec, param));
    } else {
        std::fprintf(stderr, "unknown op %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
