// Test-only probe of ComputePointValues with initialized working bounds.
#include "StdAfx.h"
#define private public
#include <SmSurfaceDropCurve.h>
#undef private
#define main curve_drop_oracle_main
#include "../drop-curve/oracle.cpp"
#undef main

int main() {
    SmContext context;
    std::string line, op;
    std::getline(std::cin, line);
    in.str(line);
    in >> op;
    std::printf("%s", op.c_str());
    SmSurface *surface = rany(context);
    SmExtent2d domain = rdomain();
    SmBSplineCurve *curve = rcurve(context);
    double a = rd(), b = rd();
    SmExtent1d interval(a, b);
    double u = rd(), v = rd(), t = rd();
    SmSurfaceDropCurve tracer(surface, domain, *curve, interval);
    tracer.m_vWorkingInterval = interval;
    tracer.m_vWorkingUVDomain = domain;
    SmTracePnt point;
    point.m_dCurveParameter = t;
    if (status(tracer.ComputePointValues(SmPoint2d(u, v), point))) {
        for (const SmVector3d &p : {point.m_vSurfacePV[0][0], point.m_vSurfacePV[1][0],
                                  point.m_vSurfacePV[0][1], point.m_vSurfacePV[2][2],
                                  point.m_v3DCurvePV[1]}) { od(p.x); od(p.y); od(p.z); }
        od(point.m_vUVCurvePV[0].x); od(point.m_vUVCurvePV[0].y);
        od(point.m_vUVCurvePV[1].x); od(point.m_vUVCurvePV[1].y);
        od(point.m_dSurfDropDist); od(point.m_dSurfNormLineDist);
    }
    std::printf("\n");
    return 0;
}
