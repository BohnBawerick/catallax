//! Exact contour generators of an extruded Bezier profile.
//!
//! The face is `S(u, h) = C(u) + h D` for a quadratic or cubic Bezier `C` and
//! the extrusion vector `D`, between two exact levels `h0 < h1`. Its normal
//! `C'(u) x D` does not depend on `h`, so the generators are whole rulings
//! `u = u*` where `f(u) = (C'(u) x D) . d` vanishes, as on a cylinder
//! (do Carmo, "Differential Geometry of Curves and Surfaces", section 3-5,
//! ruled surfaces; OpenCASCADE's `Contap_ContAna` gives the cylinder's two
//! rulings the same way). `f` has Bernstein coefficients
//! `n ((P_(i+1) - P_i) x D) . d` (Farin, "Curves and Surfaces for CAGD",
//! section 4.2, the hodograph), so it is linear or quadratic and its roots
//! are exact: one division, or the quadratic formula with one square root.
//! Each ruling is an exact line between the face's two rails.

use catallax_numbers::{Decider, Real, Refusal};
use catallax_shapes::{Curve, Line, Point3, PointDef, Vector3};
use catallax_topology::{BezierExtrusionBand, EdgeId, FaceId, Model};

use super::super::{ContourArc, lies_along};
use crate::queries::sign;

/// `C(u)` by the Bernstein form of degree `n`.
fn point_at(controls: &[Point3], u: &Real) -> Point3 {
    let n = controls.len() - 1;
    let v = &Real::one() - u;
    let mut binomial = 1_i64;
    let mut sum = Vector3::new(Real::zero(), Real::zero(), Real::zero());
    for (i, control) in controls.iter().enumerate() {
        let weight =
            &Real::from_i64(binomial) * &(&u.pow_u32(i as u32) * &v.pow_u32((n - i) as u32));
        sum = sum.plus(&control.minus(&Point3::origin()).scaled(&weight));
        binomial = binomial * (n - i) as i64 / (i + 1) as i64;
    }
    Point3::origin().plus(&sum)
}

/// Roots of `f` in `[0, 1]` from its Bernstein coefficients, ascending.
/// `None` when `f` vanishes everywhere.
fn roots(c: &[Real], dec: &mut Decider) -> Result<Option<Vec<Real>>, Refusal> {
    let mut zero = true;
    for x in c {
        zero &= sign(x, dec)? == 0;
    }
    if zero {
        return Ok(None);
    }
    let inside = |u: &Real, dec: &mut Decider| -> Result<bool, Refusal> {
        Ok(sign(u, dec)? >= 0 && sign(&(&Real::one() - u), dec)? >= 0)
    };
    // Power form a u^2 + b u + c; a is zero for a quadratic profile.
    let (a, b, k) = match c {
        [c0, c1] => (Real::zero(), c1 - c0, c0.clone()),
        [c0, c1, c2] => (&(c0 - &(c1 + c1)) + c2, &(c1 - c0) + &(c1 - c0), c0.clone()),
        _ => {
            return Err(Refusal::new(
                "app.contour-trim",
                "an extruded profile span is neither quadratic nor cubic",
            ));
        }
    };
    let mut out = Vec::new();
    if sign(&a, dec)? == 0 {
        if sign(&b, dec)? != 0 {
            let u = -(&k / &b);
            if inside(&u, dec)? {
                out.push(u);
            }
        }
        return Ok(Some(out));
    }
    let discriminant = &b.square() - &(&Real::from_i64(4) * &(&a * &k));
    let twice = &a + &a;
    match sign(&discriminant, dec)? {
        s if s < 0 => {}
        0 => {
            let u = -(&b / &twice);
            if inside(&u, dec)? {
                out.push(u);
            }
        }
        _ => {
            let root = discriminant.sqrt();
            let minus_b = -&b;
            for u in [&(&minus_b - &root) / &twice, &(&minus_b + &root) / &twice] {
                if inside(&u, dec)? {
                    out.push(u);
                }
            }
            if out.len() == 2 && sign(&(&out[0] - &out[1]), dec)? > 0 {
                out.swap(0, 1);
            }
        }
    }
    Ok(Some(out))
}

/// The exact ruling generators of one extruded Bezier face.
pub(super) fn generators(
    scratch: &mut Model,
    body_faces: &[FaceId],
    face: FaceId,
    view: &Vector3,
    arcs: &mut Vec<ContourArc>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let (support, [low, high], rulings) = {
        let band = BezierExtrusionBand::of_face(scratch, face, dec)?.ok_or_else(|| {
            Refusal::new(
                "app.contour-trim",
                "an extruded Bezier face needs its complete band",
            )
        })?;
        let rulings: Vec<EdgeId> = band
            .rows()
            .iter()
            .map(|row| row.0)
            .filter(|&edge| {
                matches!(
                    scratch.geom().curve(scratch.edge(edge).curve),
                    Curve::Line(_)
                )
            })
            .collect();
        let [low, high] = band.axial();
        (band.support().clone(), [low.clone(), high.clone()], rulings)
    };
    let controls = support.controls();
    let direction = support.direction();
    let coefficients: Vec<Real> = controls
        .windows(2)
        .map(|w| w[1].minus(&w[0]).cross(direction).dot(view))
        .collect();
    // Seen along the extrusion, or a flat band: no generator.
    let Some(roots) = roots(&coefficients, dec)? else {
        return Ok(());
    };
    for u in roots {
        dec.check_deadline()?;
        let base = point_at(controls, &u);
        let start = base.plus(&direction.scaled(&low));
        let end = base.plus(&direction.scaled(&high));
        // A ruling at either end of the profile is the band's side edge.
        let on_faces = if sign(&u, dec)? == 0 || sign(&(&Real::one() - &u), dec)? == 0 {
            let middle = base.plus(&direction.scaled(&(&(&low + &high) / &Real::from_i64(2))));
            lies_along(scratch, body_faces, face, &middle, &rulings, dec)?.faces
        } else {
            vec![face]
        };
        let geom = scratch.geom_mut();
        let curve = geom.add_curve(Curve::Line(Line {
            origin: start.clone(),
            dir: direction.clone(),
        }));
        let start = geom.add_point(PointDef::Explicit(start));
        let end = geom.add_point(PointDef::Explicit(end));
        arcs.push(ContourArc {
            face,
            curve,
            start: Some(start),
            end: Some(end),
            on_faces,
            fit: None,
        });
    }
    Ok(())
}
