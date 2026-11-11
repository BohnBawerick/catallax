// Test-only oracle driver for the catallax-trace surface intersector gate.
// Reads one request line on stdin and prints one response line in the token
// format of ../foundation/oracle.cpp. A surface is `s du dv nu nv rational
// knots_u... knots_v... (x y z w)...` with every knot listed, Euclidean
// coordinates, weights, and v varying fastest, or `p` and a plane as in
// ../plane/oracle.cpp. Catallax never compiles or links this file.
#include "StdAfx.h"
#include <SmBSplineSurface.h>
#include <SmPlane.h>
#include <SmAdvSurfaceIntersector.h>
#include <SmSurfaceCache.h>
#include <SmTree.h>
#include <SmAxis2Placement.h>
#include <SmContext.h>
#include <SmSolutionArray.h>
#include <nurbs.h>
#include <SmBSplineCurve.h>
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
// Refusal cases send SMLib's answer to the error stream (kept in the raw
// evidence) and print only the refusal marker.
static FILE *sink = stdout;
static void od(double value) {
    uint64_t bits;
    std::memcpy(&bits, &value, sizeof bits);
    std::fprintf(sink, " d%016" PRIx64, bits);
}
static void oi(long value) { std::fprintf(sink, " i%ld", value); }
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

// Reads the protected state of the intersector after each stage.
class Probe : public SmAdvSurfaceIntersector {
public:
    Probe(const SmSurface &a, const SmExtent2d &da, const SmSurface &b, const SmExtent2d &db)
        : SmAdvSurfaceIntersector(a, da, b, db) {}
    void search(const double *tolerance) {
        if (tolerance) m_dThisApproxTol3d = *tolerance;
        oi(1000);
        od(m_dThisApproxTol3d);
        for (int i = 0; i < 2; ++i) {
            oi(m_bClosedU[i] ? 1 : 0); oi(m_bClosedV[i] ? 1 : 0);
            od(m_vUVDomain[i].GetMin().x); od(m_vUVDomain[i].GetMin().y);
            od(m_vUVDomain[i].GetMax().x); od(m_vUVDomain[i].GetMax().y);
        }
        SmSurfaceCache *c1 = smsurf_GetSurfaceCache(m_cpSurface[0]);
        SmSurfaceCache *c2 = smsurf_GetSurfaceCache(m_cpSurface[1]);
        m_apTrees[0] = c1->GetTree();
        m_apTrees[1] = c2->GetTree();
        SmExtent3d box = m_apTrees[0]->GetTopNode()->m_sBBox;
        box.ExpandAbsolute(m_dThisApproxTol3d);
        if (box.AreDisjoint(m_apTrees[1]->GetTopNode()->m_sBBox)) { oi(-1); return; }
        SmSolutionArray solutions;
        SmStatus st = SolveIt(SM_SO_INTERSECTION_TEST, SM_SR_ALL, m_dThisApproxTol3d, SM_BIG_DOUBLE, NULL, solutions);
        oi((long)st);
        oi((long)solutions.GetSize());
        if (solutions.GetSize() > 0) {
            for (int k = 0; k < 7; ++k) od(solutions[0].m_vStart.m_adParameters[k]);
            for (int k = 0; k < 8; ++k) od(solutions[0].m_vEnd.m_adParameters[k]);
        }
        if (!m_planes || solutions.GetSize() == 0) return;
        SmSolution &s = solutions[0];
        SmExtent3d bbox;
        bbox.SetMinMax(SmPoint3d(s.m_vStart[0], s.m_vStart[1], s.m_vStart[2]),
                       SmPoint3d(s.m_vStart[3], s.m_vStart[4], s.m_vStart[5]));
        SmTArray<SmCurve *> curves;
        SmTArray<SmTsectCurveType> types;
        SmTArray<double> deviations;
        SmBoundaryPlaneTestResult result;
        SmStatus pst = DoBoundaryPlaneTest(c1->GetPatchBoundaryBoundingPlanes(), c2->GetPatchBoundaryBoundingPlanes(),
                                           m_dThisApproxTol3d, bbox, curves, types, deviations, result);
        if (status(pst)) {
            if (result == SM_BR_GRAZE) {
                std::fprintf(stderr, "graze with %lu curves\n", (unsigned long)curves.GetSize());
                oi(-1);
            } else {
                oi(result == SM_BR_DISJOINT ? 1 : 0);
            }
        }
    }
    bool m_planes = false;
    // The whole DoIntersection: each curve's type and deviation, then the
    // stored 3D curve and the two parameter curves.
    static void raw(SmCurve *curve) {
        const gw_CURVE *c = ((SmBSplineCurve *)curve)->GetOrCreateGwNurbPointer();
        oi(c->p); oi(c->pol->n + 1); oi(c->knt->m + 1);
        for (long i = 0; i <= c->knt->m; ++i) od(c->knt->U[i]);
        for (long i = 0; i <= c->pol->n; ++i) {
            od(c->pol->Pw[i].x); od(c->pol->Pw[i].y); od(c->pol->Pw[i].z); od(c->pol->Pw[i].w);
        }
    }
    void intersect(const SmContext &context, const double *tolerance) {
        if (tolerance) m_dThisApproxTol3d = *tolerance;
        SmTArray<SmCurve *> curves, uv1, uv2;
        SmTArray<SmTsectCurveType> types;
        SmTArray<double> deviations;
        if (!status(DoIntersection(context, NULL, NULL, &curves, &uv1, &uv2, &types, &deviations))) return;
        oi((long)curves.GetSize());
        for (ULONG i = 0; i < curves.GetSize(); ++i) {
            oi((long)types[i]);
            od(deviations[i]);
            raw(curves[i]);
            raw(uv1[i]);
            raw(uv2[i]);
        }
    }
    // DoIntersection up to FindBoundaryStartPoints: the overlap box, the
    // side-plane test, then the start points, through points and touching
    // curves found on the sides.
    void boundary(const SmContext &context, const double *tolerance) {
        if (tolerance) m_dThisApproxTol3d = *tolerance;
        m_cpContext = &context;
        SmSurfaceCache *c1 = smsurf_GetSurfaceCache(m_cpSurface[0]);
        SmSurfaceCache *c2 = smsurf_GetSurfaceCache(m_cpSurface[1]);
        m_apTrees[0] = c1->GetTree();
        m_apTrees[1] = c2->GetTree();
        SmExtent3d box = m_apTrees[0]->GetTopNode()->m_sBBox;
        box.ExpandAbsolute(m_dThisApproxTol3d);
        if (box.AreDisjoint(m_apTrees[1]->GetTopNode()->m_sBBox)) { oi(-1); return; }
        SmSolutionArray solutions;
        SmStatus st = SolveIt(SM_SO_INTERSECTION_TEST, SM_SR_ALL, m_dThisApproxTol3d, SM_BIG_DOUBLE, NULL, solutions);
        if (st != SM_SUCCESS) { oi((long)st); return; }
        if (solutions.GetSize() == 0) { oi(-1); return; }
        SmSolution &s = solutions[0];
        SmExtent3d bbox;
        bbox.SetMinMax(SmPoint3d(s.m_vStart[0], s.m_vStart[1], s.m_vStart[2]),
                       SmPoint3d(s.m_vStart[3], s.m_vStart[4], s.m_vStart[5]));
        if (s.m_vStart[6] < 2.0) {
            SmVector3d size = bbox.GetSize();
            double limit = 2.0 * m_dThisApproxTol3d;
            if (size.x < limit || size.y < limit || size.z < limit) { oi(-2); return; }
        }
        oi(0);
        SmTArray<SmCurve *> curves, uv1, uv2;
        SmTArray<SmTsectCurveType> types;
        SmTArray<double> deviations;
        SmBoundaryPlaneTestResult result;
        if (!status(DoBoundaryPlaneTest(c1->GetPatchBoundaryBoundingPlanes(), c2->GetPatchBoundaryBoundingPlanes(),
                                        m_dThisApproxTol3d, bbox, curves, types, deviations, result))) return;
        if (result == SM_BR_GRAZE) { oi(-3); return; }
        oi(result == SM_BR_DISJOINT ? 1 : 0);
        if (result == SM_BR_DISJOINT) return;
        if (!status(FindBoundaryStartPoints(bbox, curves, uv1, uv2, types, deviations))) return;
        SmTArray<SmTsectPnt *> starts, throughs;
        m_vStartPoints.GetAllNodes(starts);
        m_vThroughPoints.GetAllNodes(throughs);
        oi((long)starts.GetSize());
        for (ULONG i = 0; i < starts.GetSize(); ++i) point(*starts[i]);
        oi((long)throughs.GetSize());
        oi((long)curves.GetSize());
        for (ULONG i = 0; i < curves.GetSize(); ++i) {
            oi((long)types[i]);
            od(deviations[i]);
            SmPoint3d p;
            curves[i]->EvaluatePoint(curves[i]->GetNaturalInterval().GetMin(), p);
            o3(p);
            for (SmCurve *uv : {uv1[i], uv2[i]}) {
                uv->EvaluatePoint(uv->GetNaturalInterval().GetMin(), p);
                od(p.x); od(p.y);
            }
        }
    }
    static void point(SmTsectPnt &t) {
        oi((long)t.m_ePointType);
        for (int i = 0; i < 2; ++i) { oi((long)t.UParamType(i)); oi((long)t.VParamType(i)); }
        for (int i = 0; i < 2; ++i) { od(t.UVPos(i).x); od(t.UVPos(i).y); od(t.UVDeriv(i).x); od(t.UVDeriv(i).y); }
        for (int i = 0; i < 2; ++i) {
            SmVector3d v[7] = {t.SrfPos(i), t.SrfDu(i), t.SrfDv(i), t.SrfDuu(i), t.SrfDuv(i), t.SrfDvv(i), t.SrfNorm(i)};
            for (auto &w : v) { od(w.x); od(w.y); od(w.z); }
        }
        od(t.CrvPos().x); od(t.CrvPos().y); od(t.CrvPos().z);
        od(t.CrvDeriv().x); od(t.CrvDeriv().y); od(t.CrvDeriv().z);
        od(t.m_dTangentPlaneAngleRad); od(t.m_dDeviation);
    }
    // ClassifyBoundaryStartPointIntersection at the meeting point near a
    // pair.
    void classify(const double *tolerance, SmPoint2d first[2]) {
        if (tolerance) m_dThisApproxTol3d = *tolerance;
        m_dCurveTraceDirection = 1.0;
        SmBoolean met = FALSE;
        SmSolution solution;
        if (!status(m_cpSurface[0]->LocalSurfaceIntersect(m_vUVDomain[0], *m_cpSurface[1], m_vUVDomain[1],
                                                          m_dThisApproxTol3d, first[0], first[1], NULL, NULL, TRUE,
                                                          met, solution))) return;
        oi(met ? 1 : 0);
        if (!met) return;
        first[0] = SmPoint2d(solution.m_vStart[0], solution.m_vStart[1]);
        first[1] = SmPoint2d(solution.m_vStart[2], solution.m_vStart[3]);
        SmTsectPnt start;
        if (!status(ComputePointValues(first, start, NULL))) return;
        oi((long)start.m_ePointType);
        od(first[0].x); od(first[0].y); od(first[1].x); od(first[1].y);
        SmBoolean degenerate = FALSE;
        if (status(ClassifyBoundaryStartPointIntersection(start, degenerate))) oi(degenerate ? 1 : 0);
    }
    // TraceCurve from the meeting point near a pair, with no start, end or
    // through points and no earlier curves.
    void trace(const double *tolerance, SmPoint2d first[2]) {
        if (tolerance) m_dThisApproxTol3d = *tolerance;
        m_dCurveTraceDirection = 1.0;
        SmBoolean met = FALSE;
        SmSolution solution;
        if (!status(m_cpSurface[0]->LocalSurfaceIntersect(m_vUVDomain[0], *m_cpSurface[1], m_vUVDomain[1],
                                                          m_dThisApproxTol3d, first[0], first[1], NULL, NULL, TRUE,
                                                          met, solution))) return;
        oi(met ? 1 : 0);
        if (!met) return;
        first[0] = SmPoint2d(solution.m_vStart[0], solution.m_vStart[1]);
        first[1] = SmPoint2d(solution.m_vStart[2], solution.m_vStart[3]);
        SmTsectPnt start;
        if (!status(ComputePointValues(first, start, NULL))) return;
        m_pStartPoint = NULL;
        m_pEndPoint = NULL;
        SmTArray<SmCurve *> curves;
        if (!status(TraceCurve(start, curves))) return;
        SmTArray<SmTsectPnt *> points;
        m_vCurvePoints.GetAllNodes(points);
        oi((long)points.GetSize());
        oi(m_bCurveIsClosed ? 1 : 0);
        for (ULONG i = 0; i < points.GetSize(); ++i) {
            SmTsectPnt &t = *points[i];
            oi((long)t.m_ePointType);
            od(t.UVPos(0).x); od(t.UVPos(0).y); od(t.UVPos(1).x); od(t.UVPos(1).y);
            od(t.CrvPos().x); od(t.CrvPos().y); od(t.CrvPos().z);
            od(t.CrvDeriv().x); od(t.CrvDeriv().y); od(t.CrvDeriv().z);
            od(t.m_dCurveParameter); od(t.m_dDeviation); od(t.m_dTangentPlaneAngleRad);
        }
    }
    // ComputeNextPoint from the point at a pair, with a given step.
    void step(const double *tolerance, SmPoint2d first[2], double size) {
        if (tolerance) m_dThisApproxTol3d = *tolerance;
        m_dCurveTraceDirection = 1.0;
        // Start on the meeting curve: the local intersection from the pair.
        SmBoolean met = FALSE;
        SmSolution solution;
        if (!status(m_cpSurface[0]->LocalSurfaceIntersect(m_vUVDomain[0], *m_cpSurface[1], m_vUVDomain[1],
                                                          m_dThisApproxTol3d, first[0], first[1], NULL, NULL, TRUE,
                                                          met, solution))) return;
        oi(met ? 1 : 0);
        if (!met) return;
        first[0] = SmPoint2d(solution.m_vStart[0], solution.m_vStart[1]);
        first[1] = SmPoint2d(solution.m_vStart[2], solution.m_vStart[3]);
        SmTsectPnt a, b;
        if (!status(ComputePointValues(first, a, NULL))) return;
        a.m_dCurveParameter = 0.0;
        double next_size = 0.0;
        ComputeStepSize(a, size, next_size);
        od(next_size);
        SmBoolean good = FALSE, clipped = FALSE, hit = FALSE;
        double deviation = 0.0, angle = 0.0;
        if (!status(ComputeNextPoint(a, next_size, b, good, deviation, angle, clipped, hit))) return;
        oi(good ? 1 : 0); oi(clipped ? 1 : 0); oi(hit ? 1 : 0); od(deviation); od(angle);
        if (good) { point(b); od(b.m_dCurveParameter); }
    }
    // ComputePointValues at one pair, then at a second with the first as
    // the previous point.
    void points(const double *tolerance, SmPoint2d first[2], SmPoint2d second[2], bool two) {
        if (tolerance) m_dThisApproxTol3d = *tolerance;
        SmTsectPnt a, b;
        if (!status(ComputePointValues(first, a, NULL))) return;
        point(a);
        if (!two) return;
        if (!status(ComputePointValues(second, b, &a))) return;
        point(b);
    }
};

int main() {
    SmContext context;
    std::string line;
    std::getline(std::cin, line);
    in.str(line);
    std::string op;
    in >> op;
    std::printf("%s", op.c_str());
    SmSurface *a = rany(context);
    SmExtent2d da = rdomain();
    SmSurface *b = rany(context);
    SmExtent2d db = rdomain();
    long has_tolerance = ri();
    double tolerance = rd();
    Probe probe(*a, da, *b, db);
    if (op == "intersect" || op == "irefuse") {
        if (op == "irefuse") sink = stderr;
        probe.intersect(context, has_tolerance ? &tolerance : NULL);
        if (op == "irefuse") { sink = stdout; std::fprintf(stderr, "\n"); oi(-1); }
    } else if (op == "boundary" || op == "brefuse") {
        if (op == "brefuse") sink = stderr;
        probe.boundary(context, has_tolerance ? &tolerance : NULL);
        if (op == "brefuse") { sink = stdout; std::fprintf(stderr, "\n"); oi(-1); }
    } else if (op == "classify") {
        SmPoint2d first[2];
        for (auto &q : first) { double u = rd(), v = rd(); q = SmPoint2d(u, v); }
        probe.classify(has_tolerance ? &tolerance : NULL, first);
    } else if (op == "trace") {
        SmPoint2d first[2];
        for (auto &q : first) { double u = rd(), v = rd(); q = SmPoint2d(u, v); }
        probe.trace(has_tolerance ? &tolerance : NULL, first);
    } else if (op == "step") {
        SmPoint2d first[2];
        for (auto &q : first) { double u = rd(), v = rd(); q = SmPoint2d(u, v); }
        double size = rd();
        probe.step(has_tolerance ? &tolerance : NULL, first, size);
    } else if (op == "point") {
        long two = ri();
        SmPoint2d first[2], second[2];
        for (auto &q : first) { double u = rd(), v = rd(); q = SmPoint2d(u, v); }
        for (auto &q : second) { double u = rd(), v = rd(); q = SmPoint2d(u, v); }
        probe.points(has_tolerance ? &tolerance : NULL, first, second, two != 0);
    } else if (op == "search" || op == "planes") {
        probe.m_planes = op == "planes";
        probe.search(has_tolerance ? &tolerance : NULL);
    } else {
        std::fprintf(stderr, "unknown operation %s\n", op.c_str());
        return 2;
    }
    std::printf("\n");
    return 0;
}
