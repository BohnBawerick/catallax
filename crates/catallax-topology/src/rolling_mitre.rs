//! The whole band a ball sweeps over a mitre, read from a face's edges.

use catallax_numbers::{Decider, NotExact, Outcome, Rational, Real, Refusal};
use catallax_shapes::{BoundaryRow, Curve, Point3, RollingMitre, SplineSupport, Surface, Vector3};

use crate::{BodyId, EdgeId, FaceId, Loop, Model, Sense};

// An edge, its coedge sense in the face, and whether it starts at its row's
// first chart point.
type Row = (EdgeId, Sense, bool);

/// A face that keeps its [`RollingMitre`] carrier's whole chart: one loop of
/// the end section's arc and both contact rails, meeting at the band's two
/// feet and its apex.
#[derive(Clone, Copy, Debug)]
pub struct RollingMitreBand<'a> {
    /// The band's face.
    pub face: FaceId,
    /// Its exact support.
    pub support: &'a RollingMitre,
    // The contact rails, then the end section.
    rows: [Row; 3],
}

fn identical([a, b]: [[&Real; 3]; 2]) -> bool {
    (0..3).all(|i| a[i].same_definition(b[i]))
}

/// Exact identity tests for recognising a band. A test that is not proven
/// counts as a mismatch, and the first undecided one is kept for the caller.
struct Tests<'d> {
    dec: &'d mut Decider,
    undecided: Option<NotExact>,
}

impl Tests<'_> {
    fn zero(&mut self, value: Real) -> bool {
        match self.dec.sign(&value) {
            Outcome::ProvenCoincident(_) => true,
            Outcome::NotExact(n) => {
                self.undecided.get_or_insert(n);
                false
            }
            _ => false,
        }
    }
    // Equal definitions prove equal values without any refinement; only
    // different constructions of one value need a proof.
    fn equal(&mut self, a: &Real, b: &Real) -> bool {
        a.same_definition(b) || self.zero(a - b)
    }
    fn same(&mut self, a: &Point3, b: &Point3) -> bool {
        identical([[&a.x, &a.y, &a.z], [&b.x, &b.y, &b.z]]) || self.zero(a.minus(b).norm_squared())
    }
    fn parallel(&mut self, a: &Vector3, b: &Vector3) -> bool {
        identical([[&a.x, &a.y, &a.z], [&b.x, &b.y, &b.z]]) || self.zero(a.cross(b).norm_squared())
    }
}

impl<'a> RollingMitreBand<'a> {
    /// Checks that the face is exactly the carrier's complete band. Other
    /// trims, a missing boundary or a moved end remain unproved.
    pub fn of_face(
        model: &'a Model,
        face: FaceId,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        Ok(Self::checked(model, face, dec)?.ok())
    }

    /// As [`Self::of_face`], but a face that is not the whole band returns the
    /// first exact identity test that stopped undecided, if any, so that a
    /// caller's refusal can keep it.
    pub fn checked(
        model: &'a Model,
        face: FaceId,
        dec: &mut Decider,
    ) -> Result<Result<Self, Option<NotExact>>, Refusal> {
        let mut tests = Tests {
            dec,
            undecided: None,
        };
        let band = Self::recognise(model, face, &mut tests)?;
        Ok(band.ok_or(tests.undecided))
    }

    fn recognise(
        model: &'a Model,
        face: FaceId,
        tests: &mut Tests<'_>,
    ) -> Result<Option<Self>, Refusal> {
        let f = model.face(face);
        let Surface::BSpline(SplineSupport::RollingMitre(support)) =
            model.geom().surface(f.surface)
        else {
            return Ok(None);
        };
        let [lp] = f.loops[..] else {
            return Ok(None);
        };
        let Loop::Coedges(coedges) = model.loop_(lp) else {
            return Ok(None);
        };
        if coedges.len() != 3 {
            return Ok(None);
        }
        let section = support.section(tests.dec)?;
        let rails = [support.rail(0, tests.dec)?, support.rail(1, tests.dec)?];
        let feet =
            [Rational::zero(), Rational::one()].map(|v| support.point(&Rational::zero(), &v));
        let apex = support.apex();
        let mut seen = [false; 3];
        let mut rows: [Option<Row>; 3] = [None; 3];
        for coedge in coedges {
            let use_ = model.coedge(*coedge);
            let edge = model.edge(use_.edge);
            let (Some(start), Some(end)) = (edge.start, edge.end) else {
                return Ok(None);
            };
            let ends = [start, end].map(|v| model.geom().position(model.vertex(v).point));
            let joins = |a: &Point3, b: &Point3, tests: &mut Tests<'_>| {
                (tests.same(&ends[0], a) && tests.same(&ends[1], b))
                    || (tests.same(&ends[0], b) && tests.same(&ends[1], a))
            };
            match model.geom().curve(edge.curve) {
                Curve::Circle(circle)
                    if !seen[2]
                        && tests.same(&circle.center, &section.center)
                        && tests.equal(&circle.radius_sq, &section.radius_sq)
                        && tests.parallel(&circle.normal, &section.normal)
                        && joins(&feet[0], &feet[1], tests) =>
                {
                    seen[2] = true;
                    rows[2] = Some((use_.edge, use_.sense, tests.same(&ends[0], &feet[0])));
                }
                Curve::Ellipse(ellipse) => {
                    let mut matched = false;
                    for (i, rail) in rails.iter().enumerate() {
                        if !seen[i]
                            && tests.same(&ellipse.center, &rail.center)
                            && tests.parallel(&ellipse.normal, &rail.normal)
                            && tests.equal(&ellipse.major_sq, &rail.major_sq)
                            && joins(&feet[i], &apex, tests)
                        {
                            seen[i] = true;
                            rows[i] = Some((use_.edge, use_.sense, tests.same(&ends[0], &feet[i])));
                            matched = true;
                            break;
                        }
                    }
                    if !matched {
                        return Ok(None);
                    }
                }
                _ => return Ok(None),
            }
        }
        let [Some(a), Some(b), Some(c)] = rows else {
            return Ok(None);
        };
        Ok((seen == [true; 3]).then_some(Self {
            face,
            support,
            rows: [a, b, c],
        }))
    }
    /// Whether a body holds a band over a mitre and all its spline faces are
    /// whole exact bands: such bands or rolling pipes. A band that is not
    /// its carrier's whole chart, or a body that fails its checks, refuses.
    pub fn in_body(model: &'a Model, body: BodyId, dec: &mut Decider) -> Result<bool, Refusal> {
        let mut found = false;
        let mut other_spline = false;
        for face in model.faces_of_body(body) {
            match model.geom().surface(model.face(face).surface) {
                Surface::BSpline(SplineSupport::RollingMitre(_)) => {
                    if let Err(undecided) = Self::checked(model, face, dec)? {
                        return Err(Refusal::new(
                            catallax_numbers::GEOMETRIC_GAP_UNPROVEN,
                            "a rolling mitre band is not its carrier's whole chart",
                        )
                        .keeping_any(undecided));
                    }
                    found = true;
                }
                Surface::BSpline(SplineSupport::RollingPipe(_)) => {}
                Surface::BSpline(_) => other_spline = true,
                _ => {}
            }
        }
        if found && !model.check(body).is_empty() {
            return Err(Refusal::new(
                catallax_numbers::GEOMETRIC_GAP_UNPROVEN,
                "a body with a rolling mitre band has an unproved boundary",
            ));
        }
        Ok(found && !other_spline)
    }
    /// The boundary rows in display order: the first rail along chart `u` at
    /// `v = 0`, the second at `v = 1`, and the end section along `v` at
    /// `u = 0`. The apex column has no edge. The sense is the face's coedge
    /// sense; `forward` starts at that row's first chart point.
    pub fn rows(&self) -> impl Iterator<Item = (EdgeId, BoundaryRow, Sense, bool)> + '_ {
        self.rows
            .iter()
            .zip([BoundaryRow::Rail, BoundaryRow::Mouth, BoundaryRow::Negative])
            .map(|((edge, sense, forward), row)| (*edge, row, *sense, *forward))
    }
}

/// A starved decision must reach the caller with its code and bound intact.
/// See docs/NOT-EXACT-INVENTORY.md.
#[cfg(test)]
mod undecided_evidence {
    use super::*;
    use catallax_numbers::{BudgetPolicy, NotExact};

    fn starved() -> Decider {
        Decider::new(BudgetPolicy {
            max_precision_bits: 64,
            max_work: 0,
        })
    }

    fn undecided(value: &Real) -> NotExact {
        match starved().sign(value) {
            Outcome::NotExact(decision) => decision,
            other => panic!("an empty allowance must stop undecided: {other:?}"),
        }
    }

    #[test]
    fn recognition_keeps_its_first_undecided_identity_test() {
        let first = Point3::new(Real::from(2).sqrt(), Real::from(0), Real::from(0));
        let second = Point3::new(Real::from(3).sqrt(), Real::from(0), Real::from(0));
        let origin = Point3::new(Real::from(1), Real::from(0), Real::from(0));
        let decision = undecided(&first.minus(&origin).norm_squared());
        let mut dec = starved();
        let mut tests = Tests {
            dec: &mut dec,
            undecided: None,
        };
        assert!(!tests.same(&first, &origin));
        assert!(!tests.same(&second, &origin));
        // An identical definition needs no decision and records none.
        assert!(tests.same(&first, &first.clone()));
        assert_eq!(tests.undecided, Some(decision));
    }
}
