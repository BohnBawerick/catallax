// Test-only oracle driver for the catallax-trace spline object gate.
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

#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmNurbsCrv.h>
#include <SmNurbsSrf.h>
#include <nurbs.h>
#include <vector>

extern gw_FLAG sm_SurfaceGetENet(gw_SURFACE*,gw_INDEX,gw_INDEX,gw_INDEX,gw_INDEX,NL_ENET*);

static gw_CPOINT cp() { gw_CPOINT p; p.x=rd(); p.y=rd(); p.z=rd(); p.w=rd(); return p; }
static void raw_curve(gw_CURVE *c) {
    oi(c->p);oi(c->pol->n+1);oi(c->knt->m+1);
    for(long i=0;i<=c->knt->m;++i) od(c->knt->U[i]);
    for(long i=0;i<=c->pol->n;++i) {od(c->pol->Pw[i].x);od(c->pol->Pw[i].y);od(c->pol->Pw[i].z);od(c->pol->Pw[i].w);}
}
static void raw_surface(gw_SURFACE *s) {
    oi(s->p);oi(s->q);oi(s->net->n+1);oi(s->net->m+1);
    oi(s->knu->m+1);oi(s->knv->m+1);
    for(long i=0;i<=s->knu->m;++i) od(s->knu->U[i]);
    for(long i=0;i<=s->knv->m;++i) od(s->knv->U[i]);
    for(long i=0;i<=s->net->n;++i) for(long j=0;j<=s->net->m;++j) {
        auto p=s->net->Pw[i][j];od(p.x);od(p.y);od(p.z);od(p.w);
    }
}
static void knots(gw_KNOTVECTOR *k) {
    SmTArray<double> unique; SmTArray<ULONG> mult;
    sm_GetKnots(k,unique,&mult,NULL);
    oi(unique.GetSize());
    for (ULONG i=0;i<unique.GetSize();++i) { od(unique[i]); oi(mult[i]); }
}
int main() {
    SmContext ctx;
    std::string line, op;
    std::getline(std::cin,line); in.str(line); in >> op;
    std::printf("%s",op.c_str());
    if(op=="curve") {
        long p=ri(), n=ri(), dim=ri();
        std::vector<double> u(n+p+1);
        for(auto &v:u) v=rd();
        std::vector<gw_CPOINT> pts(n);
        for(auto &v:pts) v=cp();
        gw_CURVE *raw=sm_AllocateNurbCurve(n-1,p,u.size()-1);
        for(size_t i=0;i<u.size();++i) raw->knt->U[i]=u[i];
        for(long i=0;i<n;++i) raw->pol->Pw[i]=pts[i];
        SmBSplineCurve c(raw,dim,TRUE,&ctx);
        oi(c.GetDegree()); oi(c.GetNumberControlPoints()); oi(c.GetNumberNaturalKnots());
        oi(c.IsRational()); oi(c.GetDim());
        SmExtent1d iv=c.GetNaturalInterval(); od(iv.GetMin()); od(iv.GetMax());
        SmPoint3d a,b; double wa,wb; c.GetEnds(a,b,&wa,&wb); o3(a); od(wa); o3(b); od(wb);
        for(long i=0;i<n+1;++i) od(c.GetGrevilleAbscissa(i));
        knots(raw->knt);
        for(int form=0;form<3;++form) for(long i=0;i<n;++i) {
            SmPoint3d pt; double w;
            if(status(c.GetControlPoint((SmControlPointFormType)form,i,pt,w))) { o3(pt); od(w); }
        }
        double tol=rd(); oi(c.IsDegenerate(tol)); oi(c.HasKnotMultiplicityGreaterThanDegree());
        oi(c.IsPiecewiseC0()); oi(c.HasRepeatedEndControlPoints());
        SmExtent3d box; if(status(c.CalculateBoundingBox(iv,&box,NULL,NULL,FALSE))) ob3(box);
        SmExtent1d partial(iv.Evaluate(.2),iv.Evaluate(.8));
        if(status(c.CalculateBoundingBox(partial,&box,NULL,NULL,FALSE))) ob3(box);
        double t=rd(); long der=ri(), left=ri(), nonzero=ri();
        gw_INDEX span; N_BasisFindSpan(raw->knt,p,t,left?NL_LEFT:NL_RIGHT,&span); oi(span);
        SmVector3d values[5];
        if(status(c.Evaluate(t,der,left,values,nonzero))) for(long i=0;i<=der;++i) o3(values[i]);
        gw_CURVE *l=NULL,*r=NULL;
        long split=sm_CrvSplit(raw,iv.GetMid(),l,r);oi(split);
        if(!split) {raw_curve(l);raw_curve(r);}
        sm_FreeNurbCurve(l);sm_FreeNurbCurve(r);
        sm_FreeNurbCurve(raw);
    } else if(op=="surface") {
        long p=ri(),q=ri(),n=ri(),m=ri();
        gw_SURFACE *raw=sm_AllocateNurbSurface(n-1,m-1,p,q,n+p,m+q);
        for(long i=0;i<n+p+1;++i) raw->knu->U[i]=rd();
        for(long i=0;i<m+q+1;++i) raw->knv->U[i]=rd();
        for(long i=0;i<n;++i) for(long j=0;j<m;++j) raw->net->Pw[i][j]=cp();
        SmBSplineSurface s(raw,TRUE,&ctx);
        for(auto axis : {SM_SP_U,SM_SP_V}) {
            oi(s.GetDegree(axis)); oi(s.GetNumberControlPoints(axis)); oi(s.GetNumberNaturalKnots(axis));
        }
        oi(s.IsRational()); oi(s.IsBilinear()); oi(s.HasKnotMultiplicityGreaterThanDegree()); oi(s.HasRepeatedEndControlPoints());
        SmExtent2d domain=s.GetNaturalUVDomain(); ob2(domain);
        for(auto axis : {SM_SP_U,SM_SP_V}) for(long i=0;i<(axis==SM_SP_U?n:m)+1;++i) od(s.GetGrevilleAbscissa(axis,i));
        knots(raw->knu); knots(raw->knv);
        for(int form=0;form<3;++form) for(long i=0;i<n;++i) for(long j=0;j<m;++j) {
            SmPoint3d pt; double w;
            if(status(s.GetControlPoint((SmControlPointFormType)form,i,j,pt,w))) { o3(pt); od(w); }
        }
        std::vector<NL_POINT> epoints(n*m); std::vector<NL_POINT*> erows(n);
        for(long i=0;i<n;++i) erows[i]=epoints.data()+i*m;
        NL_ENET enet={}; enet.P=erows.data();
        long enet_flag=sm_SurfaceGetENet(raw,0,n-1,0,m-1,&enet);oi(enet_flag);
        if(!enet_flag) for(auto p:epoints) {od(p.x);od(p.y);od(p.z);}
        double tol=rd(); oi(sm_IsNSrfDegenerate(raw,tol));
        SmPoint3d test=r3(); long count=ri();
        double lo[11]={},hi[11]={}; ULONG li[11][2]={},hi_idx[11][2]={},nl=0,nh=0;
        if(status(sm_FindNetExtrema(raw,test,count,lo,li,hi,hi_idx,nl,nh))) {
            oi(nl); oi(nh);
            for(ULONG i=0;i<nl;++i) { od(lo[i]); oi(li[i][0]); oi(li[i][1]); }
            for(ULONG i=0;i<nh;++i) { od(hi[i]); oi(hi_idx[i][0]); oi(hi_idx[i][1]); }
        }
        SmPoint2d uv0=r2(),uv1=r2(); ULONG u0,v0,u1,v1;
        s.FindSpans(uv0,FALSE,FALSE,u0,v0); s.FindSpans(uv1,TRUE,TRUE,u1,v1);
        oi(u0);oi(v0);oi(u1);oi(v1);
        SmExtent3d box; SmPoint3d c00,c10,c01,c11; double ch[2],ang[2];
        if(status(sm_ComputePartialNetConstants(raw,u0,v0,u1,v1,&box,NULL,c00,c10,c01,c11,&ch[0],&ch[1],&ang[0],&ang[1]))) {
            o3(c00);o3(c01);o3(c10);o3(c11);ob3(box);od(ch[0]);od(ch[1]);od(ang[0]);od(ang[1]);
        }
        for(auto axis : {SM_SP_U,SM_SP_V}) {
            SmExtent1d iv=axis==SM_SP_U?domain.GetUInterval():domain.GetVInterval();
            for(double t : {iv.GetMin()-5e-7,iv.GetMin(),iv.GetMid(),iv.GetMax(),iv.GetMax()+5e-7,iv.GetMax()+2e-6}) {
            SmBSplineCurve *iso=NULL;
            if(status(s.CreateIsoParametricCurve(ctx,axis,t,SmApproxTol3d(0.),iso))) {
                oi(iso->GetDegree()); oi(iso->GetNumberControlPoints());
                for(ULONG i=0;i<iso->GetNumberControlPoints();++i) {
                    SmPoint3d pt;double w;iso->GetControlPoint(SM_CP_HOMOGENEOUS_RATIONAL,i,pt,w);o3(pt);od(w);
                }
            }
            delete iso;
            }
        }
        for(int axis=0;axis<2;++axis) {
            double t=(axis==0?domain.GetUInterval():domain.GetVInterval()).GetMid();
            gw_SURFACE *l=NULL,*r=NULL;
            long split=sm_SplitSrf(raw,t,axis==0?NL_UDIR:NL_VDIR,l,r);oi(split);
            if(!split) {raw_surface(l);raw_surface(r);}
            smos_Free(l);smos_Free(r);
        }
        smos_Free(raw);
    } else return 4;
    std::puts("");
}
