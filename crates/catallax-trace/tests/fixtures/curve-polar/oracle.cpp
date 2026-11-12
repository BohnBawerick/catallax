// Test-only driver for the polar parts of CalculateBoundingBox.
// Reuse the output-curve gate's binary64 protocol and constructors.
#define main output_curves_main
#include "../output-curves/oracle.cpp"
#undef main
#include <SmPolarBox.h>

static void polar(const SmCurve &curve, const SmExtent1d &part) {
    SmPolarBox box;
    if (!status(curve.CalculateBoundingBox(part, NULL, NULL, &box, FALSE))) return;
    oi(box.GetPolarBoxState());
    // SMLib leaves the basis unset until a nonzero vector is added.
    if (box.GetPolarBoxState() == SM_PS_UNINITIALIZED) return;
    for (ULONG i = 0; i < 3; ++i) o3(*box.GetBasis(i));
    SmExtent2d domain = box.GetDomain();
    od(domain.GetMin().x); od(domain.GetMin().y);
    od(domain.GetMax().x); od(domain.GetMax().y);
}

static void line_polar(SmLine &curve) {
    rd(); rd(); ri(); rd();
    double a = rd(), b = rd();
    long reparam = ri();
    double lo = rd(), hi = rd();
    polar(curve, SmExtent1d(a, b));
    if (reparam) {
        if (!status(curve.EditParameterization(SmExtent1d(lo, hi)))) return;
        SmExtent1d now = curve.GetNaturalInterval();
        polar(curve, SmExtent1d(now.Evaluate(0.2), now.Evaluate(0.9)));
    }
}

int main() {
    SmContext ctx;
    std::string request, op;
    std::getline(std::cin, request); in.str(request); in >> op;
    std::printf("%s", op.c_str());
    if (op == "line") {
        long mode = ri();
        SmPoint3d p = r3();
        SmVector3d v = r3();
        if (mode == 0) {
            double a = rd(), b = rd(), scale = rd();
            long dim = ri(), inside = ri();
            SmLine curve(p, v, SmExtent1d(a, b), scale, dim, &ctx, NULL, inside);
            oi(SM_SUCCESS);
            line_polar(curve);
        } else {
            long dim = ri(), has = ri();
            double a = rd(), b = rd();
            SmExtent1d range(a, b);
            SmLine *curve = NULL;
            if (status(SmLine::CreateLineSegment(ctx, dim, p, v, curve, has ? &range : NULL)))
                line_polar(*curve);
        }
    } else if (op == "hermite") {
        long planar = ri();
        SmVector3d p1 = r3(), d1 = r3(), p2 = r3(), d2 = r3();
        rd(); ri(); ri();
        SmHermiteCurve *curve = planar
            ? new (ctx) SmHermiteCurve(SmPoint2d(p1.x, p1.y), SmVector2d(d1.x, d1.y), SmPoint2d(p2.x, p2.y), SmVector2d(d2.x, d2.y), 2)
            : new (ctx) SmHermiteCurve(p1, d1, p2, d2, 3);
        for (SmExtent1d part : {SmExtent1d(0., 1.), SmExtent1d(.2, .7), SmExtent1d(.5, .5)})
            polar(*curve, part);
    } else {
        return 2;
    }
    std::printf("\n");
}
