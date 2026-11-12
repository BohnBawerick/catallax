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

use crate::import_step::rational;

/// A spline face the field does not hold, with its reason.
fn spline(detail: impl Into<String>) -> Refusal {
    Refusal::new("field.import-spline", detail)
}
use catallax_numbers::{Rational, Refusal, Sign};
use catallax_shapes::{Curve, RPoint3, RVector3, Surface};

fn rpoint(p: &catallax_shapes::Point3) -> Result<RPoint3, Refusal> {
    p.as_rational()
        .ok_or_else(|| spline("a spline control point is not rational after placement"))
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
        return Err(spline(
            "a spline face is not a single rational patch; only bilinear patches build in the field",
        ));
    };
    let carrier = patch.carrier();
    if carrier.degrees() != [1, 1] || carrier.control_count() != [2, 2] {
        return Err(spline(
            "a spline face has a degree above one; only bilinear patches build in the field, the rest wait for M10",
        ));
    }
    let weights = carrier
        .weights()
        .iter()
        .map(|w| rational(w, "a spline weight").map_err(|e| spline(e.detail)))
        .collect::<Result<Vec<_>, _>>()?;
    if weights.iter().any(|w| *w != weights[0]) {
        return Err(spline(
            "a spline face has unequal weights; only unweighted bilinear patches build in the field",
        ));
    }
    let controls = carrier
        .controls()
        .iter()
        .map(rpoint)
        .collect::<Result<Vec<_>, _>>()?;
    let [p00, p01, p10, p11] = <[RPoint3; 4]>::try_from(controls)
        .map_err(|_| spline("a bilinear patch does not have four controls"))?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{Decider, Real};
    use catallax_shapes::{BSplineSurface, Point3, RationalPatch, SplineSupport};

    fn patch(corners: [[i64; 3]; 4], weights: [i64; 4]) -> Surface {
        // `u` is the first index: P00, P01, P10, P11.
        let mut dec = Decider::default();
        let controls = corners
            .iter()
            .map(|c| {
                RPoint3::new(
                    Rational::from_i64(c[0]),
                    Rational::from_i64(c[1]),
                    Rational::from_i64(c[2]),
                )
                .to_real()
            })
            .collect::<Vec<Point3>>();
        let knots = vec![
            Rational::zero(),
            Rational::zero(),
            Rational::one(),
            Rational::one(),
        ];
        let carrier = BSplineSurface::new(
            [1, 1],
            [knots.clone(), knots],
            [2, 2],
            controls,
            weights.map(Real::from_i64).to_vec(),
            &mut dec,
        )
        .unwrap();
        Surface::BSpline(SplineSupport::RationalPatch(
            RationalPatch::new(carrier, &mut dec).unwrap(),
        ))
    }

    #[test]
    fn a_flat_trapezoid_patch_with_a_twist_vector_has_the_plane_and_normal_of_its_corners() {
        // The base is a trapezoid, so d = P11 - P10 - P01 + P00 is not zero.
        let surface = patch([[0, 0, 4], [0, 3, 4], [2, 0, 4], [3, 3, 4]], [1; 4]);
        let flat = bilinear(&surface).unwrap().unwrap();
        let (_, b, c, d) = flat.basis();
        assert!(!d.is_zero());
        assert!(flat.is_flat());
        // du x dv at the centre points along +z when u runs along x and v along y.
        let normal = flat.natural_normal();
        assert!(normal.x.is_zero() && normal.y.is_zero());
        assert_eq!(normal.z.sign(), Sign::Positive);
        assert_eq!(b.cross(&c).z.sign(), Sign::Positive);
        // Swapping the parameters turns the normal over.
        let turned = patch([[0, 0, 4], [2, 0, 4], [0, 3, 4], [3, 3, 4]], [1; 4]);
        let turned = bilinear(&turned).unwrap().unwrap();
        assert_eq!(turned.natural_normal().z.sign(), Sign::Negative);
        // The plane's support equation vanishes on the corners and not off them.
        let equation = flat_support(&surface).unwrap().implicit().unwrap();
        let value = |x: i64, y: i64, z: i64| equation.eval_real(&[x, y, z].map(Real::from_i64));
        assert!(value(3, 3, 4).rational_ref().unwrap().is_zero());
        assert!(!value(3, 3, 5).rational_ref().unwrap().is_zero());
    }

    #[test]
    fn spline_helper_refusals_carry_the_spline_code() {
        let weighted = patch([[0, 0, 4], [0, 3, 4], [2, 0, 4], [2, 3, 4]], [1, 2, 1, 1]);
        let refusal = bilinear(&weighted).err().expect("unequal weights refuse");
        assert_eq!(refusal.code, "field.import-spline");
        assert!(
            refusal.detail.contains("unequal weights"),
            "{}",
            refusal.detail
        );
    }
}
