//! Physical, oriented four-row domains of checked exact rational patches.

use std::collections::BTreeSet;

use catallax_numbers::{Decider, GEOMETRIC_GAP_UNPROVEN, Outcome, Real, Refusal, Sign, ZeroProof};
use catallax_shapes::{BoundaryRow, Point3, RationalPatch, SplineBoundary, SplineSupport, Surface};

use crate::{EdgeId, FaceId, Loop, Model, Sense};

fn bad(detail: &str) -> Refusal {
    Refusal::new(GEOMETRIC_GAP_UNPROVEN, detail)
}

fn same(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(&a.minus(b).norm_squared()) {
        Outcome::Proven(_) => Ok(false),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => Ok(true),
        Outcome::ProvenCoincident(_) => Err(bad(
            "a patch row endpoint needs exact correspondence, not bounded coincidence",
        )),
        Outcome::NotExact(error) => Err(Refusal::from(error)),
        Outcome::Refused(error) => Err(error),
    }
}

/// A complete finite patch face, proved from its physical boundary curves,
/// endpoint identity, sewn vertex cycle and material-relative loop direction.
/// This does not certify nonintersection of different faces in a body.
pub struct RationalPatchFace<'a> {
    face: FaceId,
    support: &'a RationalPatch,
    rows: [(EdgeId, BoundaryRow, Sense, bool); 4],
}

impl<'a> RationalPatchFace<'a> {
    /// Checks all four complete rows. Additional trims, missing or duplicate
    /// rows, forged endpoint agreement and wrong orientation retain a refusal.
    pub fn of_face(
        model: &'a Model,
        face: FaceId,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        let Surface::BSpline(SplineSupport::RationalPatch(support)) =
            model.geom().surface(model.face(face).surface)
        else {
            return Ok(None);
        };
        Self::check(model, face, support, dec).map(Some)
    }

    /// As [`Self::of_face`], and also for a face on an atlas of one chart,
    /// checked as that chart's patch.
    pub fn of_chart_face(
        model: &'a Model,
        face: FaceId,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        match model.geom().surface(model.face(face).surface) {
            Surface::BSpline(SplineSupport::RationalPatch(support)) => {
                Self::check(model, face, support, dec).map(Some)
            }
            Surface::BSpline(SplineSupport::RationalAtlas(atlas)) if atlas.charts().len() == 1 => {
                Self::check(model, face, atlas.charts()[0].patch(), dec).map(Some)
            }
            _ => Ok(None),
        }
    }

    fn check(
        model: &'a Model,
        face: FaceId,
        support: &'a RationalPatch,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let definition = model.face(face);
        let [loop_id] = definition.loops.as_slice() else {
            return Err(bad("a complete rational patch needs one boundary loop"));
        };
        let Loop::Coedges(coedges) = model.loop_(*loop_id) else {
            return Err(bad("a rational patch needs physical boundary edges"));
        };
        if coedges.len() != 4 {
            return Err(bad(
                "the rational patch has missing rows or additional trims",
            ));
        }
        let [nu, nv] = support.carrier().control_count();
        let controls = support.carrier().controls();
        let corners = [
            &controls[0],
            &controls[(nu - 1) * nv],
            &controls[nu * nv - 1],
            &controls[nv - 1],
        ];
        let sides = [
            SplineBoundary::VMin,
            SplineBoundary::VMax,
            SplineBoundary::UMin,
            SplineBoundary::UMax,
        ];
        let labels = [
            BoundaryRow::Rail,
            BoundaryRow::Mouth,
            BoundaryRow::Negative,
            BoundaryRow::Positive,
        ];
        let endpoints = [[0, 1], [3, 2], [0, 3], [1, 2]];
        let positive_walk = [true, false, false, true];
        let mut rows = [None; 4];
        let mut seen = BTreeSet::new();
        let mut cycle = Vec::with_capacity(4);
        for coedge in coedges {
            dec.check_deadline()?;
            let use_ = model.coedge(*coedge);
            let edge = model.edge(use_.edge);
            if !seen.insert(use_.edge) {
                return Err(bad("the rational patch repeats a physical edge"));
            }
            let (Some(a), Some(b)) = (edge.start, edge.end) else {
                return Err(bad("a rational patch row needs two finite vertices"));
            };
            let start = model.geom().position(model.vertex(a).point);
            let end = model.geom().position(model.vertex(b).point);
            cycle.push(if use_.sense == Sense::Same {
                (a, b)
            } else {
                (b, a)
            });
            let mut matched = false;
            for row in 0..4 {
                let [lo, hi] = endpoints[row].map(|i| corners[i]);
                let forward = if same(&start, lo, dec)? && same(&end, hi, dec)? {
                    true
                } else if same(&start, hi, dec)? && same(&end, lo, dec)? {
                    false
                } else {
                    continue;
                };
                let curve_forward = support
                    .boundary_orientation(sides[row], model.geom().curve(edge.curve), dec)?
                    .ok_or_else(|| {
                        bad("a patch row has matching endpoints but no whole-curve identity")
                    })?;
                if rows[row].is_some() {
                    return Err(bad("the rational patch repeats a parameter row"));
                }
                if (edge.sense == Sense::Same) != (curve_forward == forward) {
                    return Err(bad(
                        "a rational patch edge reverses its natural carrier parameter",
                    ));
                }
                let wanted = positive_walk[row] == (definition.sense == Sense::Same);
                if ((use_.sense == Sense::Same) == forward) != wanted {
                    return Err(bad(
                        "the rational patch loop disagrees with its chart orientation",
                    ));
                }
                rows[row] = Some((use_.edge, labels[row], use_.sense, forward));
                matched = true;
                break;
            }
            if !matched {
                return Err(bad(
                    "a physical edge does not cover a complete rational patch row",
                ));
            }
        }
        for i in 0..4 {
            if cycle[i].1 != cycle[(i + 1) % 4].0 {
                return Err(bad(
                    "rational patch rows do not share one sewn vertex cycle",
                ));
            }
        }
        let [Some(a), Some(b), Some(c), Some(d)] = rows else {
            return Err(bad("the rational patch is missing a boundary row"));
        };
        Ok(Self {
            face,
            support,
            rows: [a, b, c, d],
        })
    }

    /// The face whose complete physical parameter domain was checked.
    pub fn face(&self) -> FaceId {
        self.face
    }
    /// The immutable exact carrier and its embedding certificate.
    pub fn support(&self) -> &'a RationalPatch {
        self.support
    }
    /// VMin, VMax, UMin, UMax rows, with physical coedge use and whether the edge
    /// start-to-end direction agrees with the increasing chart parameter.
    pub fn rows(&self) -> [(EdgeId, BoundaryRow, Sense, bool); 4] {
        self.rows
    }

    /// Closed finite membership: negative in the chart interior, zero on an
    /// actual boundary row, positive off the support or outside its domain.
    pub fn locate(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        let Some(parameters) = self.support.parameters_at(point, dec)? else {
            return Ok(Sign::Positive);
        };
        for (parameter, domain) in parameters.iter().zip(self.support.domain()) {
            for bound in domain {
                match dec.sign(&(parameter - &Real::from_rational(bound.clone()))) {
                    Outcome::Proven(_) => (),
                    Outcome::ProvenCoincident(proof)
                        if proof.proof != ZeroProof::SeparationBound =>
                    {
                        return Ok(Sign::Zero);
                    }
                    Outcome::ProvenCoincident(_) => {
                        return Err(bad("finite patch membership needs exact row incidence"));
                    }
                    Outcome::NotExact(error) => return Err(Refusal::from(error)),
                    Outcome::Refused(error) => return Err(error),
                }
            }
        }
        Ok(Sign::Negative)
    }
}
