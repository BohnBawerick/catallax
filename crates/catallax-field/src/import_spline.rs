//! Imported spline faces that the field can hold exactly.
//!
//! The shared importer admits a spline face only as a complete embedded
//! rational patch whose boundary rows its edges cover. This module reads such a
//! patch back from the imported model when its definition is a bilinear patch
//! with equal weights. Such a patch is the image of the unit square under
//! `a + u b + v c + u v d`. Its three vectors either lie in one plane, in which
//! case the face is flat and the planar import route holds it as that plane, or
//! they span space, in which case the patch lies on the quadric `gamma = alpha
//! beta` of its own coordinates (see `quadric`).
//!
//! Anything else (a degree above one, unequal weights, a patch whose controls
//! are not rational) stays a typed refusal that names the face and the reason.

use crate::import_step::{rational, unsupported};
use catallax_numbers::{Rational, Refusal, Sign};
use catallax_shapes::{Curve, RPoint3, RVector3, Surface};

fn rpoint(p: &catallax_shapes::Point3) -> Result<RPoint3, Refusal> {
    p.as_rational()
        .ok_or_else(|| unsupported("a spline control point is not rational after placement"))
}

/// The four controls of a bilinear equal-weight patch, `[P00, P01, P10, P11]`
/// with `u` the first index.
pub(crate) struct Bilinear {
    pub corners: [RPoint3; 4],
}

impl Bilinear {
    /// `a`, `b`, `c` and `d` of `a + u b + v c + u v d`.
    pub(crate) fn basis(&self) -> (RPoint3, RVector3, RVector3, RVector3) {
        let [p00, p01, p10, p11] = &self.corners;
        let b = p10.minus(p00);
        let c = p01.minus(p00);
        let d = p11.minus(p10).minus(&p01.minus(p00));
        (p00.clone(), b, c, d)
    }

    /// Whether the four controls lie in one plane, so the patch is flat.
    pub(crate) fn is_flat(&self) -> bool {
        let (_, b, c, d) = self.basis();
        b.cross(&c).dot(&d).is_zero()
    }

    /// The patch normal `du x dv` at the centre of the unit square. For an
    /// embedded patch the sign is the same everywhere.
    pub(crate) fn natural_normal(&self) -> RVector3 {
        let (_, b, c, d) = self.basis();
        let half = Rational::from_ratio(1, 2).expect("half");
        b.cross(&c)
            .plus(&b.cross(&d).minus(&c.cross(&d)).scaled(&half))
    }
}

/// Reads a bilinear equal-weight patch from a model surface, `None` for any
/// other surface.
pub(crate) fn bilinear(surface: &Surface) -> Result<Option<Bilinear>, Refusal> {
    use catallax_shapes::SplineSupport;
    let Surface::BSpline(support) = surface else {
        return Ok(None);
    };
    let SplineSupport::RationalPatch(patch) = support else {
        return Err(unsupported(
            "a spline face is not a single rational patch; only bilinear patches build in the field",
        ));
    };
    let carrier = patch.carrier();
    if carrier.degrees() != [1, 1] || carrier.control_count() != [2, 2] {
        return Err(unsupported(
            "a spline face has a degree above one; only bilinear patches build in the field, the rest wait for M10",
        ));
    }
    let weights = carrier
        .weights()
        .iter()
        .map(|w| rational(w, "a spline weight"))
        .collect::<Result<Vec<_>, _>>()?;
    if weights.iter().any(|w| *w != weights[0]) {
        return Err(unsupported(
            "a spline face has unequal weights; only unweighted bilinear patches build in the field",
        ));
    }
    let controls = carrier
        .controls()
        .iter()
        .map(rpoint)
        .collect::<Result<Vec<_>, _>>()?;
    let [p00, p01, p10, p11] = <[RPoint3; 4]>::try_from(controls)
        .map_err(|_| unsupported("a bilinear patch does not have four controls"))?;
    Ok(Some(Bilinear {
        corners: [p00, p01, p10, p11],
    }))
}

/// Whether a spline edge is exactly a straight segment: every control lies on
/// the line through its first and last.
pub(crate) fn straight(curve: &Curve) -> bool {
    let Curve::BSpline(spline) = curve else {
        return false;
    };
    let points: Option<Vec<RPoint3>> = spline.controls().map(|p| p.as_rational()).collect();
    let Some(points) = points else {
        return false;
    };
    let (Some(first), Some(last)) = (points.first(), points.last()) else {
        return false;
    };
    let axis = last.minus(first);
    !axis.is_zero()
        && points.iter().all(|p| p.minus(first).cross(&axis).is_zero())
        && points
            .iter()
            .all(|p| p.minus(first).dot(&axis).sign() != Sign::Negative)
}

/// The plane that carries a flat bilinear spline patch, oriented by the patch's
/// own normal, or `None` for every other surface. Comparison code uses it to
/// give a flat spline face of a B-rep the same support equation as the plane
/// the field holds in its place.
pub fn flat_support(surface: &Surface) -> Option<Surface> {
    let patch = bilinear(surface).ok()??;
    if !patch.is_flat() || patch.natural_normal().is_zero() {
        return None;
    }
    let (origin, b, c, _) = patch.basis();
    Some(Surface::Plane(catallax_shapes::Plane {
        origin: origin.to_real(),
        normal: patch.natural_normal().to_real(),
        ref_dir: if b.is_zero() { c } else { b }.to_real(),
    }))
}
