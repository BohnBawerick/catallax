// Fixture producer for SMLib c9c5979f04cea8c38213004d4f04df0dcfb40d05.
// Catallax does not compile or link this driver.
#include "StdAfx.h"
#include <SmPlane.h>
#include <SmBSplineSurface.h>
#include <SmBSplineCurve.h>
#include <SmContext.h>
#include <nurbs.h>
#include <cmath>
#include <cstdio>
#include <string>

int main(int argc, char **argv) {
    const std::string which = argc > 1 ? argv[1] : "lofttight";
    const bool loft = which == "loft" || which == "lofttight";
    const bool hill = which == "hill" || which == "tiny" || which == "tangent";
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
        std::printf("HEIGHTS");
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
            const double z = hill ? ((i >= 2 && i <= 3 && j >= 2 && j <= 3) ? 1.0 : 0.0)
                : std::sin(0.9 * i + 0.3) * std::cos(0.7 * j - 0.3) + 0.05 * i;
            points.Add(SmPoint3d(i, j, z));
            std::printf(" %.17g", z);
        }
        std::printf("\n");
        const ULONG multiplicities[4] = {4, 1, 1, 4};
        for (int i = 0; i < 4; ++i) {
            mu.Add(multiplicities[i]); mv.Add(multiplicities[i]); ku.Add(i); kv.Add(i);
        }
    }
    SmBSplineSurface *surface = NULL;
    if (SmBSplineSurface::CreateCanonical(ctx, loft ? 2 : 3, 3, points, SM_SF_UNSPECIFIED,
        mu, mv, ku, kv, SM_KT_UNSPECIFIED, loft ? &weights : NULL, NULL, surface) != SM_SUCCESS || !surface) return 2;
    SmAxis2Placement frame;
    if (loft) frame.SetCanonical(SmPoint3d(2.5, 0, 0), SmVector3d(1, -1, 0), SmVector3d(1, 1, -2));
    else {
        const double height = which == "hill" ? 0.3 : which == "tiny" ? 225.0/256 - 1e-6
            : which == "tangent" ? 225.0/256 : 0.35;
        frame.SetCanonical(SmPoint3d(0, 0, height), SmVector3d(1, 0, 0), SmVector3d(0, 1, 0));
    }
    SmPlane *plane = NULL;
    if (SmPlane::CreateCanonical(ctx, frame, plane) != SM_SUCCESS || !plane) return 2;
    SmBoolean edges[2] = {TRUE, TRUE};
    SmTArray<SmCurve *> c3, uv1, uv2;
    SmTArray<SmTsectCurveType> types;
    SmTArray<double> deviation;
    SmApproxTol3d tolerance(1e-5);
    const SmStatus status = surface->GlobalSurfaceIntersect(ctx, surface->GetNaturalUVDomain(), *plane,
        SmExtent2d(-20, -20, 20, 20), edges, which == "lofttight" ? &tolerance : NULL,
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
