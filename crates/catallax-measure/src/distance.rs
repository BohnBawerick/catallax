//! Certified trimmed-face distances, independent of tessellation.

use catallax_intersect::trace;
use catallax_numbers::{Ball, Decider, Rational, Refusal};
use catallax_shapes::{
    Curve, Point3, RPoint3, RVector3, Surface,
    distance::{self, Bounds, Domain, Finite, Membership, Options},
};
use catallax_topology::{FaceId, Loop, Model, Sense};
use std::cell::RefCell;

/// Dyadic probe levels of the parallel planar witness search: at most 85 points.
const PLANAR_LEVELS: u32 = 3;

/// Reusable analytic chart cover and exact boundary-loop traces for one face.
/// Valid only while the model's face definitions and boundary loops are unchanged.
pub struct PreparedFace {
    face: FaceId,
    extent: crate::BoundingBox,
    domain: Domain,
    boundary: Vec<[Ball; 3]>,
    finite: Option<Finite>,
}
fn marked() -> Refusal {
    Refusal::new(
        "measure.distance-marked",
        "distance requires unmarked face bounds and boundary traces",
    )
}
impl PreparedFace {
    /// Builds a certified compact cover using analytic extrema and exact trims.
    pub fn new(model: &Model, face: FaceId, dec: &mut Decider) -> Result<Self, Refusal> {
        let bounds = crate::face_bounding_box(
            model,
            face,
            &crate::Accuracy {
                relative: Rational::zero(),
                absolute: Rational::from_ratio(1, 1_000_000).unwrap(),
            },
            dec,
        )?;
        if !bounds.marks.is_empty() {
            return Err(marked());
        }
        let extent = bounds.value;
        let bounds = std::array::from_fn(|i| {
            Ball::from_rational(&extent.min[i], 128).hull(&Ball::from_rational(&extent.max[i], 128))
        });
        let domain = Domain::new(model.geom().surface(model.face(face).surface), bounds, dec)?;
        let mut boundary = Vec::new();
        for lp in &model.face(face).loops {
            match model.loop_(*lp) {
                Loop::Vertex(v) => {
                    let p = model.geom().position(model.vertex(*v).point);
                    boundary.push([
                        dec.enclose_budgeted_at(&p.x, 128)?,
                        dec.enclose_budgeted_at(&p.y, 128)?,
                        dec.enclose_budgeted_at(&p.z, 128)?,
                    ]);
                }
                Loop::Coedges(cs) => {
                    for c in cs {
                        let edge = model.edge(model.coedge(*c).edge);
                        let (a, b) = if edge.sense == Sense::Same {
                            (edge.start, edge.end)
                        } else {
                            (edge.end, edge.start)
                        };
                        let traced = trace(
                            model.geom(),
                            edge.curve,
                            a.map(|v| model.vertex(v).point),
                            b.map(|v| model.vertex(v).point),
                            dec,
                        )?;
                        if !traced.marks.is_empty() {
                            return Err(marked());
                        }
                        for piece in traced.value.pieces() {
                            for i in 0..8 {
                                let lo = Rational::from_ratio(i, 8).unwrap();
                                let hi = Rational::from_ratio(i + 1, 8).unwrap();
                                let t = Ball::from_rational(&lo, 128)
                                    .hull(&Ball::from_rational(&hi, 128));
                                boundary.push(piece.evaluate(&t, dec)?);
                            }
                        }
                    }
                }
            }
        }
        let finite = Self::recognize(model, face, dec)?;
        Ok(Self {
            face,
            extent,
            domain,
            boundary,
            finite,
        })
    }
    /// Recognizes complete circular cylinder/disc, oriented spherical or conical
    /// regions. A cone tip may have one collapsed vertex loop. Every other loop
    /// must be one closed circle edge. Other trims keep the complete search.
    fn recognize(
        model: &Model,
        face: FaceId,
        dec: &mut Decider,
    ) -> Result<Option<Finite>, Refusal> {
        let mut rims = Vec::new();
        let mut apex = None;
        for lp in &model.face(face).loops {
            let cs = match model.loop_(*lp) {
                Loop::Vertex(v) => {
                    if apex.is_some() {
                        return Ok(None);
                    }
                    apex = Some(model.geom().position(model.vertex(*v).point));
                    continue;
                }
                Loop::Coedges(cs) => cs,
            };
            let [c] = cs.as_slice() else {
                return Ok(None);
            };
            let coedge = model.coedge(*c);
            let edge = model.edge(coedge.edge);
            if edge.start.is_some() || edge.end.is_some() {
                return Ok(None);
            }
            let Curve::Circle(circle) = model.geom().curve(edge.curve) else {
                return Ok(None);
            };
            let forward = (coedge.sense == edge.sense) == (model.face(face).sense == Sense::Same);
            rims.push((circle, forward));
        }
        if let Some(apex) = &apex {
            let catallax_shapes::Surface::Cone(cone) =
                model.geom().surface(model.face(face).surface)
            else {
                return Ok(None);
            };
            return Finite::recognize_conical(cone, &rims, Some(apex), dec);
        }
        Finite::recognize_oriented(model.geom().surface(model.face(face).surface), &rims, dec)
    }
    fn classify(
        &self,
        model: &mut Model,
        p: &Point3,
        cover: &[Ball; 3],
        dec: &mut Decider,
    ) -> Result<Membership, Refusal> {
        let witness = catallax_boolean::face_contains(model, self.face, p, dec)?;
        let boundary_free = self
            .boundary
            .iter()
            .all(|b| (0..3).any(|i| cover[i].hi() < b[i].lo() || b[i].hi() < cover[i].lo()));
        // A connected chart cell cannot enter a closed trimmed region without
        // meeting its boundary. Traces enclose every boundary loop, including holes.
        Ok(Membership {
            witness,
            excluded: !witness && boundary_free,
        })
    }
    /// Bounds point-to-face distance. Upper witnesses satisfy exact closed-loop
    /// membership; lower bounds cover all unexcluded parameter cells. Complete
    /// cylinder bands, perpendicular circular discs and oriented spherical caps
    /// or complete coaxial bands and positive-nappe conical bands (including
    /// explicitly retained tips) are answered directly.
    pub fn point_to_face(
        &self,
        model: &mut Model,
        p: &Point3,
        options: &Options,
        dec: &mut Decider,
    ) -> Result<Bounds, Refusal> {
        if let Some(finite) = &self.finite
            && let Some((bounds, _)) = finite.point_distance(p, options, dec)?
        {
            return Ok(bounds);
        }
        distance::point_distance(
            &self.domain,
            p,
            options,
            &mut |p: &Point3, b: &[Ball; 3], d: &mut Decider| self.classify(model, p, b, d),
            dec,
        )
    }
    /// Bounds minimum distance to a second prepared face in the same model.
    pub fn face_to_face(
        &self,
        other: &Self,
        model: &mut Model,
        options: &Options,
        dec: &mut Decider,
    ) -> Result<Bounds, Refusal> {
        let model = RefCell::new(model);
        distance::face_distance(
            &self.domain,
            &other.domain,
            options,
            &mut |p: &Point3, b: &[Ball; 3], d: &mut Decider| {
                self.classify(&mut model.borrow_mut(), p, b, d)
            },
            &mut |p: &Point3, b: &[Ball; 3], d: &mut Decider| {
                other.classify(&mut model.borrow_mut(), p, b, d)
            },
            dec,
        )
    }
    /// Bounds minimum distance to a parallel planar face in the same model.
    /// Both supporting planes need exact rational frames; their gap is a global
    /// lower bound. A probed point admitted by this face's exact trim, whose
    /// perpendicular translate the other trim admits, attains the gap. Probes
    /// are dyadic grid centres over both faces' overlapping boxes, coarsest
    /// first, and count as visited cells. Without such a witness the general
    /// trimmed search answers. A curved face refuses `measure.nonplanar-face`,
    /// planes that are not parallel refuse `measure.nonparallel-faces`, and an
    /// inexact frame refuses `measure.irrational-plane`.
    pub fn parallel_face_to_face(
        &self,
        other: &Self,
        model: &mut Model,
        options: &Options,
        dec: &mut Decider,
    ) -> Result<Bounds, Refusal> {
        let (origin, normal) = plane(model, self.face)?;
        let (other_origin, other_normal) = plane(model, other.face)?;
        if !normal.cross(&other_normal).is_zero() {
            return Err(Refusal::new(
                "measure.nonparallel-faces",
                format!(
                    "faces {:?} and {:?} are not parallel",
                    self.face, other.face
                ),
            ));
        }
        let along = &normal.dot(&other_origin.minus(&origin)) / &normal.norm_squared();
        let shift = normal.scaled(&along);
        let gap = crate::pair::root(shift.norm_squared(), dec)?;
        let o = [&origin.x, &origin.y, &origin.z];
        let n = [&normal.x, &normal.y, &normal.z];
        let v = [&shift.x, &shift.y, &shift.z];
        let axis = (0..3).max_by_key(|k| n[*k].abs()).unwrap();
        let [i, j] = [(axis + 1) % 3, (axis + 2) % 3];
        let range = |k: usize| {
            let lo = (&other.extent.min[k] - v[k]).max(self.extent.min[k].clone());
            let hi = (&other.extent.max[k] - v[k]).min(self.extent.max[k].clone());
            (lo <= hi).then_some((lo, hi))
        };
        let (Some((lo_i, hi_i)), Some((lo_j, hi_j))) = (range(i), range(j)) else {
            return self.face_to_face(other, model, options, dec);
        };
        let mut probes = 0;
        for level in 0..=PLANAR_LEVELS {
            let count = 1i64 << level;
            let at = |lo: &Rational, hi: &Rational, m: i64| {
                lo + &(&(hi - lo) * &Rational::from_ratio(2 * m + 1, 2 * count).unwrap())
            };
            for u in 0..count {
                for w in 0..count {
                    dec.check_deadline()?;
                    let mut p = [Rational::zero(), Rational::zero(), Rational::zero()];
                    p[i] = at(&lo_i, &hi_i, u);
                    p[j] = at(&lo_j, &hi_j, w);
                    let offset = &(n[i] * &(&p[i] - o[i])) + &(n[j] * &(&p[j] - o[j]));
                    p[axis] = o[axis] - &(&offset / n[axis]);
                    let [x, y, z] = p;
                    let p = RPoint3::new(x, y, z);
                    let q = p.plus(&shift);
                    probes += 1;
                    if catallax_boolean::face_contains(model, self.face, &p, dec)?
                        && catallax_boolean::face_contains(model, other.face, &q, dec)?
                    {
                        return Ok(Bounds {
                            converged: gap.width() <= options.width,
                            lower: gap.lo,
                            upper: Some(gap.hi),
                            witnesses: Some([p.into(), q.into()]),
                            cells: probes,
                        });
                    }
                }
            }
        }
        self.face_to_face(other, model, options, dec)
    }
}

/// Exact rational origin and normal of a planar face's supporting plane.
fn plane(model: &Model, face: FaceId) -> Result<(RPoint3, RVector3), Refusal> {
    let Surface::Plane(plane) = model.geom().surface(model.face(face).surface) else {
        return Err(Refusal::new(
            "measure.nonplanar-face",
            format!("face {face:?} is not planar"),
        ));
    };
    plane
        .origin
        .as_rational()
        .zip(plane.normal.as_rational())
        .ok_or_else(|| {
            Refusal::new(
                "measure.irrational-plane",
                format!("face {face:?} has no exact rational plane frame"),
            )
        })
}

#[cfg(test)]
mod tests;
