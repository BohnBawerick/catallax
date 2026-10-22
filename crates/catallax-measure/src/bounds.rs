//! Coordinate extrema from certified boundary traces and compact carrier supports.

use catallax_intersect::{TracePiece, trace};
use catallax_numbers::{Answer, Answered, Ball, Decider, MarkLog, Rational, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{BodyId, FaceId, Loop, Model, Sense};

use crate::{Accuracy, Enclosure};

/// A closed axis-aligned box containing the exact geometry, in model coordinates.
/// Compact curved carriers can give a conservative box for a trimmed face.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BoundingBox {
    /// Inclusive coordinate lower bounds.
    pub min: [Rational; 3],
    /// Inclusive coordinate upper bounds.
    pub max: [Rational; 3],
}

impl BoundingBox {
    fn union(&mut self, other: &Self) {
        for i in 0..3 {
            self.min[i] = self.min[i].clone().min(other.min[i].clone());
            self.max[i] = self.max[i].clone().max(other.max[i].clone());
        }
    }
}

pub(crate) fn body(
    model: &Model,
    body: BodyId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<Option<BoundingBox>> {
    accuracy.validate();
    let mut bounds: Option<BoundingBox> = None;
    let mut marks = MarkLog::new();
    for face_id in model.faces_of_body(body) {
        let answer = face(model, face_id, accuracy, dec)?;
        marks.extend(&answer.marks);
        if let Some(bounds) = &mut bounds {
            bounds.union(&answer.value);
        } else {
            bounds = Some(answer.value);
        }
    }
    Ok(Answer::with_marks(bounds, marks))
}

pub(crate) fn face(
    model: &Model,
    face_id: FaceId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<BoundingBox> {
    accuracy.validate();
    let policy = *dec.policy();
    let bits = policy.max_precision_bits.min(256);
    if bits == 0 || policy.max_work == 0 {
        return Err(Refusal::new(
            "measure.budget-exhausted",
            "bounds require coordinate enclosures",
        ));
    }
    let surface = model.geom().surface(model.face(face_id).surface);
    let compact = match surface {
        Surface::Sphere(s) => Some((&s.center, std::array::from_fn(|_| s.radius.clone()))),
        Surface::Torus(t) => {
            let axis = [&t.axis.x, &t.axis.y, &t.axis.z];
            let squared = t.axis.norm_squared();
            let radii: [catallax_numbers::Real; 3] = std::array::from_fn(|i| {
                &(&t.major
                    * &(catallax_numbers::Real::from(1) - &axis[i].square() / &squared).sqrt())
                    + &t.minor
            });
            Some((&t.center, radii))
        }
        Surface::Plane(_)
        | Surface::Cylinder(_)
        | Surface::Cone(_)
        | Surface::EllipticCylinder(_) => None,
        _ => {
            return Err(Refusal::new(
                "measure.unsupported-surface",
                "surface has no proven bounds definition",
            ));
        }
    };
    if let Some((center, radii)) = compact {
        let coordinates = [&center.x, &center.y, &center.z];
        let mut min = std::array::from_fn(|_| Rational::zero());
        let mut max = min.clone();
        let mut marks = MarkLog::new();
        for i in 0..3 {
            let low = Enclosure::from_ball(&dec.enclose_at(&(coordinates[i] - &radii[i]), bits)?);
            let high = Enclosure::from_ball(&dec.enclose_at(&(coordinates[i] + &radii[i]), bits)?);
            accuracy.mark(&mut marks, &[&low, &high], "bounding coordinate");
            min[i] = low.lo;
            max[i] = high.hi;
        }
        return Ok(Answer::with_marks(BoundingBox { min, max }, marks));
    }
    let mut pieces = Vec::new();
    let mut marks = MarkLog::new();
    let mut points = Vec::new();
    for lp in &model.face(face_id).loops {
        match model.loop_(*lp) {
            Loop::Vertex(v) => {
                let p = model.geom().position(model.vertex(*v).point);
                let xyz = [&p.x, &p.y, &p.z].map(|x| dec.enclose_at(x, bits));
                points.push([xyz[0].clone()?, xyz[1].clone()?, xyz[2].clone()?]);
            }
            Loop::Coedges(cs) => {
                for c in cs {
                    let edge = model.edge(model.coedge(*c).edge);
                    let (mut a, mut b) = (edge.start, edge.end);
                    if edge.sense == Sense::Reversed {
                        std::mem::swap(&mut a, &mut b);
                    }
                    let traced = trace(
                        model.geom(),
                        edge.curve,
                        a.map(|v| model.vertex(v).point),
                        b.map(|v| model.vertex(v).point),
                        dec,
                    )?;
                    marks.extend(&traced.marks);
                    pieces.extend(traced.value.pieces().iter().cloned());
                }
            }
        }
    }
    struct Panel {
        piece: usize,
        lo: Rational,
        hi: Rational,
        xyz: [Ball; 3],
    }
    let evaluate = |piece: &TracePiece, lo: &Rational, hi: &Rational| {
        piece.evaluate(
            &Ball::from_rational(lo, bits).hull(&Ball::from_rational(hi, bits)),
            dec,
        )
    };
    let mut panels = Vec::new();
    for (i, piece) in pieces.iter().enumerate() {
        let (lo, hi) = piece.range();
        points.push(evaluate(piece, &lo, &lo)?);
        points.push(evaluate(piece, &hi, &hi)?);
        panels.push(Panel {
            piece: i,
            xyz: evaluate(piece, &lo, &hi)?,
            lo,
            hi,
        });
    }
    if points.is_empty() {
        return Err(Refusal::new(
            "measure.unbounded-face",
            "noncompact carrier has no bounded face boundary",
        ));
    }
    let mut remaining = policy.max_work;
    loop {
        dec.check_deadline()?;
        let lower: [Enclosure; 3] = std::array::from_fn(|i| Enclosure {
            lo: panels
                .iter()
                .map(|p| p.xyz[i].lo())
                .chain(points.iter().map(|p| p[i].lo()))
                .min()
                .unwrap(),
            hi: points.iter().map(|p| p[i].hi()).min().unwrap(),
        });
        let upper: [Enclosure; 3] = std::array::from_fn(|i| Enclosure {
            lo: points.iter().map(|p| p[i].lo()).max().unwrap(),
            hi: panels
                .iter()
                .map(|p| p.xyz[i].hi())
                .chain(points.iter().map(|p| p[i].hi()))
                .max()
                .unwrap(),
        });
        let result = BoundingBox {
            min: lower.each_ref().map(|v| v.lo.clone()),
            max: upper.each_ref().map(|v| v.hi.clone()),
        };
        if lower.iter().chain(&upper).all(|v| accuracy.reached(v)) || remaining < 64 {
            accuracy.mark(&mut marks, &lower.each_ref(), "face bounding minimum");
            accuracy.mark(&mut marks, &upper.each_ref(), "face bounding maximum");
            return Ok(Answer::with_marks(result, marks));
        }
        let mut refined = Vec::new();
        for p in panels {
            let relevant = (0..3).any(|i| {
                (!accuracy.reached(&lower[i]) && p.xyz[i].lo() <= lower[i].hi)
                    || (!accuracy.reached(&upper[i]) && p.xyz[i].hi() >= upper[i].lo)
            });
            if relevant && remaining >= 64 {
                remaining -= 64;
                let mid = (&p.lo + &p.hi).mul_pow2(-1);
                points.push(evaluate(&pieces[p.piece], &mid, &mid)?);
                for (lo, hi) in [(p.lo, mid.clone()), (mid, p.hi)] {
                    refined.push(Panel {
                        piece: p.piece,
                        xyz: evaluate(&pieces[p.piece], &lo, &hi)?,
                        lo,
                        hi,
                    });
                }
            } else {
                refined.push(p);
            }
        }
        panels = refined;
    }
}
