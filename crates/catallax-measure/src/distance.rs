//! Certified trimmed-face distances, independent of tessellation.

use catallax_intersect::trace;
use catallax_numbers::{Ball, Decider, Rational, Refusal};
use catallax_shapes::{
    Point3,
    distance::{self, Bounds, Domain, Membership, Options},
};
use catallax_topology::{FaceId, Loop, Model, Sense};
use std::cell::RefCell;

/// Reusable analytic chart cover and exact boundary-loop traces for one face.
/// Valid only while the model's face definitions and boundary loops are unchanged.
pub struct PreparedFace {
    face: FaceId,
    domain: Domain,
    boundary: Vec<[Ball; 3]>,
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
        let bounds = bounds.value;
        let bounds = std::array::from_fn(|i| {
            Ball::from_rational(&bounds.min[i], 128).hull(&Ball::from_rational(&bounds.max[i], 128))
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
        Ok(Self {
            face,
            domain,
            boundary,
        })
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
    /// membership; lower bounds cover all unexcluded parameter cells.
    pub fn point_to_face(
        &self,
        model: &mut Model,
        p: &Point3,
        options: &Options,
        dec: &mut Decider,
    ) -> Result<Bounds, Refusal> {
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
}

#[cfg(test)]
mod tests;
