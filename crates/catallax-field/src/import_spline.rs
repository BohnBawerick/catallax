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
//! The curved case is a saddle roof: a closed solid of six faces, a curved
//! bilinear patch, a floor plane, and four wall planes that contain the patch's
//! twist direction `d`. Over the unit square of its coordinates the solid is
//! `{ floor <= gamma <= alpha beta }` or its mirror, and `saddle` holds it as the
//! intersection of the six supports after checking, exactly, that the imported
//! faces are those supports and meet at those eight corners.
//!
//! Anything else (a degree above one, unequal weights, a patch whose controls
//! are not rational, a curved patch in any other solid) stays a typed refusal
//! that names the face and the reason.

use crate::Face;
use crate::geometry::{self, Cell, Leaf};
use crate::import_step::{outward, rational};

/// A spline face the field does not hold, with its reason.
fn spline(detail: impl Into<String>) -> Refusal {
    Refusal::new("field.import-spline", detail)
}
use catallax_numbers::{Decider, MPoly, Rational, Real, RealPoly, Refusal, Sign};
use catallax_shapes::{Curve, RPoint3, RVector3, Surface};
use catallax_topology::{BodyId, Loop, Model, Sense};
use std::collections::BTreeSet;

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

/// The quadric `T^2 (gamma - alpha beta)` of a non-flat bilinear patch in world
/// coordinates, signed so that its gradient points along the patch normal. Here
/// `alpha`, `beta` and `gamma` are the coordinates of `x - a` in the basis
/// `b`, `c`, `d`, and `T` is their triple product. Every patch point has
/// `gamma = alpha beta`, because it is `a + u b + v c + u v d`.
pub(crate) fn quadric(patch: &Bilinear) -> Option<MPoly> {
    if patch.is_flat() {
        return None;
    }
    let (a, b, c, d) = patch.basis();
    let triple = b.cross(&c).dot(&d);
    let w = [0, 1, 2].map(|i| {
        let origin = [&a.x, &a.y, &a.z][i];
        MPoly::var(3, i) - MPoly::constant(3, origin.clone())
    });
    let linear = |v: &RVector3| w[0].scale(&v.x) + w[1].scale(&v.y) + w[2].scale(&v.z);
    let gamma = linear(&b.cross(&c));
    let alpha = linear(&c.cross(&d));
    let beta = linear(&d.cross(&b));
    let f = gamma.scale(&triple) - alpha * beta;
    // Orient along the normal at the patch centre, which is on the surface.
    let half = Rational::from_ratio(1, 2).expect("half");
    let quarter = Rational::from_ratio(1, 4).expect("quarter");
    let centre = a
        .plus(&b.scaled(&half))
        .plus(&c.scaled(&half))
        .plus(&d.scaled(&quarter));
    let at = [centre.x.clone(), centre.y.clone(), centre.z.clone()];
    let gradient = RVector3::new(
        f.partial(0).eval(&at),
        f.partial(1).eval(&at),
        f.partial(2).eval(&at),
    );
    match gradient.dot(&patch.natural_normal()).sign() {
        Sign::Positive => Some(f),
        Sign::Negative => Some(f.scale(&-Rational::one())),
        Sign::Zero => None,
    }
}

/// The support polynomial of a bilinear spline patch, in the patch's own
/// orientation: the plane of a flat patch, or the quadric of a curved one.
/// `None` for every other surface. Comparison code uses it to give a B-rep
/// spline face the same equation as the support the field holds.
pub fn support_polynomial(surface: &Surface) -> Option<RealPoly> {
    let patch = bilinear(surface).ok()??;
    if patch.is_flat() {
        return flat_support(surface)?.implicit();
    }
    quadric(&patch).map(RealPoly::from)
}

/// The ten coefficients of `support_polynomial`, ordered `xx yy zz xy xz yz x y z 1`.
pub(crate) fn quadric_coefficients(surface: &Surface) -> Option<[Real; 10]> {
    let polynomial = support_polynomial(surface)?;
    let mut out: [Real; 10] = std::array::from_fn(|_| Real::zero());
    for (monomial, coefficient) in polynomial.terms() {
        let slot = match monomial.exponents() {
            [2, 0, 0] => 0,
            [0, 2, 0] => 1,
            [0, 0, 2] => 2,
            [1, 1, 0] => 3,
            [1, 0, 1] => 4,
            [0, 1, 1] => 5,
            [1, 0, 0] => 6,
            [0, 1, 0] => 7,
            [0, 0, 1] => 8,
            [0, 0, 0] => 9,
            _ => return None,
        };
        out[slot] = coefficient;
    }
    Some(out)
}

/// A point of the quadric of a curved bilinear patch near `p`: the point on the
/// line through `p` along the patch's twist direction. It lies on the untrimmed
/// quadric exactly, and is as far from `p` as `p` is from the quadric along that
/// direction. When `p` projects outside the unit parameter square, the point
/// lies outside the finite patch, so this is neither a point of the face nor a
/// nearest point: callers keep their own trim checks. `None` for a flat patch,
/// whose plane gives its projection, and for every other surface.
pub fn patch_witness(surface: &Surface, p: &RPoint3) -> Option<RPoint3> {
    let patch = bilinear(surface).ok()??;
    if patch.is_flat() {
        return None;
    }
    let (a, b, c, d) = patch.basis();
    let triple = b.cross(&c).dot(&d);
    let w = p.minus(&a);
    let alpha = c.cross(&d).dot(&w) / triple.clone();
    let beta = d.cross(&b).dot(&w) / triple;
    Some(
        a.plus(&b.scaled(&alpha))
            .plus(&c.scaled(&beta))
            .plus(&d.scaled(&(&alpha * &beta))),
    )
}

/// The patch with its two parameters swapped, which reverses its normal.
fn transposed(surface: &Surface, dec: &mut Decider) -> Result<Surface, Refusal> {
    use catallax_shapes::{BSplineSurface, RationalPatch, SplineSupport};
    let Surface::BSpline(SplineSupport::RationalPatch(patch)) = surface else {
        return Err(spline("only a rational patch can be turned over"));
    };
    let carrier = patch.carrier();
    let [u, v] = carrier.knots();
    let mut controls = carrier.controls().to_vec();
    controls.swap(1, 2);
    let mut weights = carrier.weights().to_vec();
    weights.swap(1, 2);
    let carrier = BSplineSurface::new(
        [1, 1],
        [v.to_vec(), u.to_vec()],
        [2, 2],
        controls,
        weights,
        dec,
    )?;
    Ok(Surface::BSpline(SplineSupport::RationalPatch(
        RationalPatch::new(carrier, dec)?,
    )))
}

/// The corners of the patch face's loop and of every other face's loop, as
/// points, or `None` unless the face has one loop of four straight edges.
fn corners(model: &Model, id: catallax_topology::FaceId) -> Option<BTreeSet<[Rational; 3]>> {
    let face = model.face(id);
    let [lp] = face.loops.as_slice() else {
        return None;
    };
    let Loop::Coedges(cs) = model.loop_(*lp) else {
        return None;
    };
    if cs.len() != 4 {
        return None;
    }
    let mut points = BTreeSet::new();
    for c in cs {
        let coedge = model.coedge(*c);
        let edge = model.edge(coedge.edge);
        let curve = model.geom().curve(edge.curve);
        if !matches!(curve, catallax_shapes::Curve::Line(_)) && !straight(curve) {
            return None;
        }
        let vertex = if coedge.sense == Sense::Same {
            edge.start
        } else {
            edge.end
        }?;
        let p = model
            .geom()
            .position(model.vertex(vertex).point)
            .as_rational()?;
        points.insert([p.x, p.y, p.z]);
    }
    (points.len() == 4).then_some(points)
}

/// Holds an imported saddle roof (a curved bilinear patch, a floor and four
/// walls) as the intersection of its six supports, with its exact volume.
///
/// With `a`, `b`, `c`, `d` the patch basis and `alpha`, `beta`, `gamma` the
/// coordinates in it, the patch is `gamma = alpha beta`. The four walls must be
/// the planes `alpha = 0, 1` and `beta = 0, 1` (they contain `d`), the floor a
/// plane `gamma = p + q alpha + r beta` strictly on one side of the patch at all
/// four corners (a bilinear difference takes its extremes at the corners, so it
/// never crosses), and the faces must meet at exactly the eight corners that
/// these supports give. Then the imported solid and the intersection of the six
/// halfspaces are bounded by the same faces on the same supports with the same
/// outward normals, and are the same solid; a point strictly inside every
/// halfspace confirms the orientation.
pub(crate) fn saddle(
    model: &Model,
    body: BodyId,
    dec: &mut Decider,
) -> Result<(Leaf, Vec<Face>), Refusal> {
    let ids = model.faces_of_body(body);
    let not_a_roof = |detail: &str| {
        spline(format!(
            "a curved spline face builds only in a saddle roof of six faces: {detail}"
        ))
    };
    if ids.len() != 6 {
        return Err(not_a_roof(&format!("this part has {} faces", ids.len())));
    }
    let mut roof = None;
    let mut planes = Vec::new();
    for id in &ids {
        let surface = model.geom().surface(model.face(*id).surface);
        match bilinear(surface)? {
            Some(patch) if !patch.is_flat() => {
                if roof.replace((*id, patch)).is_some() {
                    return Err(not_a_roof("there is more than one curved patch"));
                }
            }
            _ => planes.push((*id, outward(model, *id)?)),
        }
    }
    let Some((roof_id, patch)) = roof else {
        return Err(not_a_roof("there is no curved patch"));
    };
    // The roof as the field holds it: its own normal points out of the solid.
    let reversed = model.face(roof_id).sense == Sense::Reversed;
    let roof_surface = if reversed {
        transposed(model.geom().surface(model.face(roof_id).surface), dec)?
    } else {
        model.geom().surface(model.face(roof_id).surface).clone()
    };
    let [p00, p01, p10, p11] = patch.corners;
    let oriented = Bilinear {
        corners: if reversed {
            [p00, p10, p01, p11]
        } else {
            [p00, p01, p10, p11]
        },
    };
    let (a, b, c, d) = oriented.basis();
    let t = b.cross(&c).dot(&d);
    let (cd, db) = (c.cross(&d), d.cross(&b));
    let rational_point = |i: i64, j: i64, gamma: &Rational| {
        a.plus(&b.scaled(&Rational::from_i64(i)))
            .plus(&c.scaled(&Rational::from_i64(j)))
            .plus(&d.scaled(gamma))
    };
    let negative_t = t.sign() == Sign::Negative;
    // Walls by their planes; the remaining plane is the floor.
    let walls = [
        (&cd, a.clone(), true),
        (&cd, a.plus(&b), false),
        (&db, a.clone(), true),
        (&db, a.plus(&c), false),
    ];
    let mut wall_faces: [Option<usize>; 4] = [None; 4];
    let mut floor = None;
    for (index, (_, plane)) in planes.iter().enumerate() {
        let n = plane
            .normal
            .as_rational()
            .ok_or_else(|| not_a_roof("a wall is not rational"))?;
        let o = rpoint(&plane.origin)?;
        if n.dot(&d).is_zero() {
            let mut matched = None;
            for (k, (direction, base, low)) in walls.iter().enumerate() {
                let outward_sign = n.dot(direction).sign();
                let expected = if *low != negative_t {
                    Sign::Negative
                } else {
                    Sign::Positive
                };
                if n.cross(direction).is_zero()
                    && o.minus(base).dot(direction).is_zero()
                    && outward_sign == expected
                {
                    matched = Some(k);
                }
            }
            let k =
                matched.ok_or_else(|| not_a_roof("a wall is not a plane of the patch's square"))?;
            if wall_faces[k].replace(index).is_some() {
                return Err(not_a_roof("two faces lie on one wall plane"));
            }
        } else if floor.replace(index).is_some() {
            return Err(not_a_roof("there is more than one floor"));
        }
    }
    let (Some(floor), true) = (floor, wall_faces.iter().all(Option::is_some)) else {
        return Err(not_a_roof(
            "the walls and floor are not four walls and one floor",
        ));
    };
    let floor_plane = &planes[floor].1;
    let n = floor_plane
        .normal
        .as_rational()
        .ok_or_else(|| not_a_roof("the floor is not rational"))?;
    let o = rpoint(&floor_plane.origin)?;
    // gamma = p + q alpha + r beta on the floor.
    let slope = n.dot(&d);
    let p = -(n.dot(&a.minus(&o)) / slope.clone());
    let q = -(n.dot(&b) / slope.clone());
    let r = -(n.dot(&c) / slope);
    let floor_at =
        |i: i64, j: i64| &(&p + &(&q * &Rational::from_i64(i))) + &(&r * &Rational::from_i64(j));
    let above = |i: i64, j: i64| Rational::from_i64(i * j) - floor_at(i, j);
    let gaps = [above(0, 0), above(1, 0), above(0, 1), above(1, 1)];
    let sign = gaps[0].sign();
    if sign == Sign::Zero || gaps.iter().any(|g| g.sign() != sign) {
        return Err(not_a_roof("the floor meets or crosses the patch"));
    }
    // Corners the faces must meet at.
    let roof_corner = |i, j| {
        let p = rational_point(i, j, &Rational::from_i64(i * j));
        [p.x, p.y, p.z]
    };
    let floor_corner = |i, j| {
        let p = rational_point(i, j, &floor_at(i, j));
        [p.x, p.y, p.z]
    };
    let expected = |cells: [(i64, i64, bool); 4]| -> BTreeSet<[Rational; 3]> {
        cells
            .into_iter()
            .map(|(i, j, on_roof)| {
                if on_roof {
                    roof_corner(i, j)
                } else {
                    floor_corner(i, j)
                }
            })
            .collect()
    };
    let checks = [
        (
            roof_id,
            expected([(0, 0, true), (1, 0, true), (0, 1, true), (1, 1, true)]),
        ),
        (
            planes[floor].0,
            expected([(0, 0, false), (1, 0, false), (0, 1, false), (1, 1, false)]),
        ),
        (
            planes[wall_faces[0].unwrap()].0,
            expected([(0, 0, true), (0, 1, true), (0, 0, false), (0, 1, false)]),
        ),
        (
            planes[wall_faces[1].unwrap()].0,
            expected([(1, 0, true), (1, 1, true), (1, 0, false), (1, 1, false)]),
        ),
        (
            planes[wall_faces[2].unwrap()].0,
            expected([(0, 0, true), (1, 0, true), (0, 0, false), (1, 0, false)]),
        ),
        (
            planes[wall_faces[3].unwrap()].0,
            expected([(0, 1, true), (1, 1, true), (0, 1, false), (1, 1, false)]),
        ),
    ];
    for (id, want) in &checks {
        if corners(model, *id).as_ref() != Some(want) {
            return Err(not_a_roof(
                "a face does not meet at the corners of its supports",
            ));
        }
    }
    // The six supports, roof first, each with its outward side.
    let mut surfaces = vec![roof_surface];
    let mut named = vec![roof_id];
    for (id, plane) in &planes {
        surfaces.push(Surface::Plane(plane.clone()));
        named.push(*id);
    }
    // Orientation: the point between the roof and the floor over the square's
    // centre lies strictly inside every halfspace.
    let half = Rational::from_ratio(1, 2).expect("half");
    let quarter = Rational::from_ratio(1, 4).expect("quarter");
    let floor_centre = &(&p + &(&q * &half)) + &(&r * &half);
    let centre_gamma = &(&floor_centre + &quarter) * &half;
    let inside = rational_point(0, 0, &Rational::zero())
        .plus(&b.scaled(&half))
        .plus(&c.scaled(&half))
        .plus(&d.scaled(&centre_gamma));
    let point = [inside.x.clone(), inside.y.clone(), inside.z.clone()];
    let polynomial =
        support_polynomial(&surfaces[0]).ok_or_else(|| not_a_roof("the roof has no quadric"))?;
    if geometry::sign(
        &polynomial.eval_real(&point.each_ref().map(|x| Real::from(x.clone()))),
        dec,
    )? != Sign::Negative
    {
        return Err(not_a_roof("the part is inside out at the roof"));
    }
    for (_, plane) in &planes {
        let n = plane.normal.as_rational().expect("checked rational");
        let o = rpoint(&plane.origin)?;
        if n.dot(&inside.minus(&o)).sign() != Sign::Negative {
            return Err(not_a_roof("the part is inside out at a plane"));
        }
    }
    // Volume: |T| times the mean height between the patch and the floor.
    let mean_gap = &quarter - &floor_centre;
    let volume = &t.abs() * &mean_gap.abs();
    let radius = checks
        .iter()
        .flat_map(|(_, set)| set.iter())
        .flat_map(|point| point.iter())
        .map(Rational::abs)
        .max()
        .expect("eight corners");
    let cell = Cell::new(surfaces.clone(), Some(radius))?.with_volume(volume);
    let mut faces = Vec::new();
    for (surface, id) in surfaces.iter().zip(&named) {
        for origin in model.face_origins(*id) {
            faces.push(geometry::face(surface.clone(), origin.clone(), dec)?);
        }
    }
    Ok((Leaf::Cell(cell), faces))
}

/// Points of a curved bilinear patch at the centres of a regular grid of its
/// parameters, `divisions` per direction, for candidate sampling. They lie in
/// the patch's interior, off its edges, where its face has one normal.
pub(crate) fn patch_points(surface: &Surface, divisions: u32) -> Option<Vec<RPoint3>> {
    let patch = bilinear(surface).ok()??;
    let (a, b, c, d) = patch.basis();
    let n = i64::from(divisions.max(1));
    let mut points = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let u = Rational::from_ratio(2 * i + 1, 2 * n)?;
            let v = Rational::from_ratio(2 * j + 1, 2 * n)?;
            points.push(
                a.plus(&b.scaled(&u))
                    .plus(&c.scaled(&v))
                    .plus(&d.scaled(&(&u * &v))),
            );
        }
    }
    Some(points)
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::Decider;
    use catallax_shapes::{BSplineSurface, Point3, RationalPatch, SplineSupport};

    fn patch(corners: [[i64; 3]; 4]) -> Surface {
        weighted_patch(corners, [1; 4])
    }

    fn weighted_patch(corners: [[i64; 3]; 4], weights: [i64; 4]) -> Surface {
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
        // P00 (0,0,4), P01 (0,3,4), P10 (2,0,4), P11 (3/2 * 2, 3, 4) in integers:
        // the base is a trapezoid, so d = P11 - P10 - P01 + P00 is not zero.
        let surface = patch([[0, 0, 4], [0, 3, 4], [2, 0, 4], [3, 3, 4]]);
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
        let turned = patch([[0, 0, 4], [2, 0, 4], [0, 3, 4], [3, 3, 4]]);
        let turned = bilinear(&turned).unwrap().unwrap();
        assert_eq!(turned.natural_normal().z.sign(), Sign::Negative);
        // The plane's support equation vanishes on the corners and not off them.
        let equation = support_polynomial(&surface).unwrap();
        let value = |x: i64, y: i64, z: i64| equation.eval_real(&[x, y, z].map(Real::from_i64));
        assert!(value(3, 3, 4).rational_ref().unwrap().is_zero());
        assert!(!value(3, 3, 5).rational_ref().unwrap().is_zero());
    }

    #[test]
    fn a_curved_patch_quadric_vanishes_on_its_surface_and_points_along_its_normal() {
        // z = 4 + u v / 2 over the 2 by 3 base: P11 lifted by 1/2 (here by 1).
        let surface = patch([[0, 0, 4], [0, 3, 4], [2, 0, 4], [2, 3, 5]]);
        let equation = support_polynomial(&surface).unwrap();
        let at = |x: &str, y: &str, z: &str| {
            equation.eval_real(&[x, y, z].map(|t| Real::from(Rational::parse(t).unwrap())))
        };
        // Surface points: u = x / 2, v = y / 3, z = 4 + u v.
        assert!(at("1", "3/2", "17/4").rational_ref().unwrap().is_zero());
        assert!(at("2", "3", "5").rational_ref().unwrap().is_zero());
        // Above the surface (along +z, the patch normal) the value is positive.
        assert_eq!(
            at("1", "3/2", "6").rational_ref().unwrap().sign(),
            Sign::Positive
        );
        assert_eq!(
            at("1", "3/2", "1").rational_ref().unwrap().sign(),
            Sign::Negative
        );
        let witness = patch_witness(
            &surface,
            &RPoint3::new(
                Rational::one(),
                Rational::from_ratio(3, 2).unwrap(),
                Rational::from_i64(6),
            ),
        )
        .unwrap();
        assert_eq!(witness.z, Rational::from_ratio(17, 4).unwrap());
    }

    #[test]
    fn spline_helper_refusals_carry_the_spline_code() {
        let weighted = weighted_patch([[0, 0, 4], [0, 3, 4], [2, 0, 4], [2, 3, 4]], [1, 2, 1, 1]);
        let refusal = bilinear(&weighted).err().expect("unequal weights refuse");
        assert_eq!(refusal.code, "field.import-spline");
        assert!(
            refusal.detail.contains("unequal weights"),
            "{}",
            refusal.detail
        );
    }
}
