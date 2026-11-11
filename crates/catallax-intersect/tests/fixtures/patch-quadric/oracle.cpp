// Fixture producer for SMLib c9c5979f04cea8c38213004d4f04df0dcfb40d05.
// Catallax does not compile or link this driver.
#include "StdAfx.h"
#include <SmCylinder.h>
#include <SmSphere.h>
#include <SmBSplineSurface.h>
#include <SmBSplineCurve.h>
#include <SmContext.h>
#include <nurbs.h>
#include <cstdio>
#include <string>

// Case names are "<cylinder|sphere>-<loft|lofttight|hill|tiny|tangent>".
int main(int argc, char **argv) {
    const std::string which = argc > 1 ? argv[1] : "cylinder-lofttight";
    const bool cylinder = which.rfind("cylinder-", 0) == 0;
    const std::string cut = which.substr(which.find('-') + 1);
    const bool loft = cut == "loft" || cut == "lofttight";
    SmContext ctx;
    SmTArray<SmPoint3d> points;
    SmTArray<double> weights, ku, kv;
    SmTArray<ULONG> mu, mv;
    if (loft) {
        const double bulge[4] = {1.0, 2.0, 0.5, 1.5};
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 4; ++j) {
            points.Add(SmPoint3d(i, i == 1 ? bulge[j] : 0.0, j));
            weights.Add(i == 1 ? 0.75 : 1.0);
        }
        mu.Add(3); mu.Add(3); mv.Add(4); mv.Add(4);
        ku.Add(0); ku.Add(1); kv.Add(0); kv.Add(1);
    } else {
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j)
            points.Add(SmPoint3d(i, j, (i >= 2 && i <= 3 && j >= 2 && j <= 3) ? 1.0 : 0.0));
        const ULONG multiplicities[4] = {4, 1, 1, 4};
        for (int i = 0; i < 4; ++i) {
            mu.Add(multiplicities[i]); mv.Add(multiplicities[i]); ku.Add(i); kv.Add(i);
        }
    }
    SmBSplineSurface *surface = NULL;
    if (SmBSplineSurface::CreateCanonical(ctx, loft ? 2 : 3, 3, points, SM_SF_UNSPECIFIED,
        mu, mv, ku, kv, SM_KT_UNSPECIFIED, loft ? &weights : NULL, NULL, surface) != SM_SUCCESS || !surface) return 2;
    // Hill cuts: radius 1, lowest point 0.3, 1e-6 below, or exactly at the top.
    const double peak = 225.0 / 256;
    const double low = cut == "hill" ? 0.3 : cut == "tiny" ? peak - 1e-6 : peak;
    SmSurface *other = NULL;
    SmExtent2d domain;
    if (cylinder) {
        SmAxis2Placement frame;
        // Loft: axis along y through the origin. Hill: axis along x.
        if (loft) frame.SetCanonical(SmPoint3d(0, 0, 0), SmVector3d(0, 0, 1), SmVector3d(1, 0, 0));
        else frame.SetCanonical(SmPoint3d(0, 2.5, low + 1.0), SmVector3d(0, 1, 0), SmVector3d(0, 0, 1));
        SmCylinder *c = NULL;
        if (SmCylinder::CreateCanonical(ctx, frame, loft ? 2.5 : 1.0, c) != SM_SUCCESS || !c) return 2;
        other = c;
        domain = SmExtent2d(0, -20, 360, 20);
    } else {
        SmAxis2Placement frame;
        if (loft) frame.SetCanonical(SmPoint3d(0, 0, 0), SmVector3d(1, 0, 0), SmVector3d(0, 1, 0));
        else frame.SetCanonical(SmPoint3d(2.5, 2.5, low + 1.0), SmVector3d(1, 0, 0), SmVector3d(0, 1, 0));
        SmSphere *s = NULL;
        if (SmSphere::CreateCanonical(ctx, frame, loft ? 2.5 : 1.0, s) != SM_SUCCESS || !s) return 2;
        other = s;
        domain = SmExtent2d(0, -90, 360, 90);
    }
    SmBoolean edges[2] = {TRUE, TRUE};
    SmTArray<SmCurve *> c3, uv1, uv2;
    SmTArray<SmTsectCurveType> types;
    SmTArray<double> deviation;
    SmApproxTol3d tolerance(1e-5);
    const SmStatus status = surface->GlobalSurfaceIntersect(ctx, surface->GetNaturalUVDomain(), *other,
        domain, edges, cut == "lofttight" ? &tolerance : NULL,
        NULL, &c3, &uv1, &uv2, &types, &deviation);
    std::printf("STATUS %ld CURVES %lu\n", (long)status, (unsigned long)uv1.GetSize());
    if (status != SM_SUCCESS) return 3;
    for (ULONG i = 0; i < uv1.GetSize(); ++i) {
        if (!uv1[i]->IsKindOf(SmBSplineCurve_TYPE)) return 4;
        SmBSplineCurve *curve = (SmBSplineCurve *)uv1[i];
        if (curve->GetDim() != 2 || curve->IsRational()) return 5;
        const gw_CURVE *g = curve->GetOrCreateGwNurbPointer();
        std::printf("CURVE %lu DEGREE %d\nKNOTS", (unsigned long)i, (int)g->p);
        for (long k = 0; k <= g->knt->m; ++k) std::printf(" %.17g", g->knt->U[k]);
        std::printf("\nPOINTS");
        for (long k = 0; k <= g->pol->n; ++k)
            std::printf(" %.17g %.17g 1", g->pol->Pw[k].x, g->pol->Pw[k].y);
        std::printf("\n");
    }
    return 0;
}
