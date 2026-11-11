//! Checked rectangular trims of exact polynomial ruled patches.

use crate::{BodyId, Edge, EdgeId, FaceId, Loop, Model, Sense, VertexId};
use catallax_numbers::{
    Decider, GEOMETRIC_GAP_UNPROVEN, Outcome, Rational, Real, Refusal, Sign, StrictSign, ZeroProof,
};
use catallax_shapes::{BezierExtrusion, BoundaryRow, Curve, Point3, SplineSupport, Surface};

fn bad(detail: &str) -> Refusal {
    Refusal::new(GEOMETRIC_GAP_UNPROVEN, detail)
}
fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => {
            Ok(Sign::Zero)
        }
        Outcome::ProvenCoincident(_) => Err(bad("Bezier boundary incidence needs an exact proof")),
        Outcome::NotExact(n) => Err(Refusal::from(n)),
        Outcome::Refused(error) => Err(error),
    }
}
fn same(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(&a.minus(b).norm_squared(), dec)? == Sign::Zero)
}
fn ends(
    model: &Model,
    edge: &Edge,
    wanted: [&Point3; 2],
    dec: &mut Decider,
) -> Result<([VertexId; 2], bool), Refusal> {
    let (Some(a), Some(b)) = (edge.start, edge.end) else {
        return Err(bad("a Bezier band edge needs finite endpoints"));
    };
    let p = model.geom().position(model.vertex(a).point);
    let q = model.geom().position(model.vertex(b).point);
    if same(&p, wanted[0], dec)? && same(&q, wanted[1], dec)? {
        return Ok(([a, b], true));
    }
    if same(&p, wanted[1], dec)? && same(&q, wanted[0], dec)? {
        return Ok(([b, a], false));
    }
    Err(bad(
        "a Bezier band edge does not cover its complete boundary row",
    ))
}

struct Row {
    edge: EdgeId,
    coedge: Sense,
    forward: bool,
    vertices: [VertexId; 2],
    height: Real,
}

/// A complete polynomial span between two exact axial levels.
/// Physical vertices, curve coverage and loop orientation establish its domain.
pub struct BezierExtrusionBand<'a> {
    support: &'a BezierExtrusion,
    axial: [Real; 2],
    rows: [(EdgeId, BoundaryRow, Sense, bool); 4],
}

impl<'a> BezierExtrusionBand<'a> {
    /// Checks one four-edge loop with two exact rails and two endpoint generators.
    /// Additional trims, duplicate rows and unproved incidence keep a refusal.
    pub fn of_face(
        model: &'a Model,
        face: FaceId,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        let f = model.face(face);
        let Surface::BSpline(SplineSupport::BezierExtrusion(support)) =
            model.geom().surface(f.surface)
        else {
            return Ok(None);
        };
        let [lp] = f.loops.as_slice() else {
            return Err(bad("a Bezier band needs one boundary loop"));
        };
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            return Err(bad("a Bezier band needs physical boundary edges"));
        };
        if coedges.len() != 4 {
            return Err(bad("a Bezier band has extra or missing trims"));
        }
        let mut rails = Vec::new();
        let mut generators = Vec::new();
        let mut cycle = Vec::new();
        for coedge in coedges {
            dec.check_deadline()?;
            let c = model.coedge(*coedge);
            let edge = model.edge(c.edge);
            let (Some(a), Some(b)) = (edge.start, edge.end) else {
                return Err(bad("a Bezier band cannot have a closed seam edge"));
            };
            cycle.push(if c.sense == Sense::Same {
                (a, b)
            } else {
                (b, a)
            });
            match model.geom().curve(edge.curve) {
                Curve::BSpline(curve) => {
                    let (height, curve_forward) = support
                        .section_parameter(curve, dec)?
                        .ok_or_else(|| bad("a Bezier rail must be a complete exact section"))?;
                    let offset = support.direction().scaled(&height);
                    let wanted = [
                        support.controls()[0].plus(&offset),
                        support.controls()[support.controls().len() - 1].plus(&offset),
                    ];
                    let (vertices, forward) = ends(model, edge, wanted.each_ref(), dec)?;
                    if (edge.sense == Sense::Same) != (curve_forward == forward) {
                        return Err(bad("a Bezier rail reverses its parameter sense"));
                    }
                    rails.push(Row {
                        edge: c.edge,
                        coedge: c.sense,
                        forward,
                        vertices,
                        height,
                    });
                }
                Curve::Line(line) => generators.push((c, edge, line)),
                _ => return Err(bad("a Bezier band has an unproved boundary curve")),
            }
        }
        if rails.len() != 2 || generators.len() != 2 {
            return Err(bad("a Bezier band needs two distinct rails and generators"));
        }
        for i in 0..4 {
            if cycle[i].1 != cycle[(i + 1) % 4].0 {
                return Err(bad("Bezier band edges do not form one sewn vertex cycle"));
            }
        }
        match sign(&(&rails[0].height - &rails[1].height), dec)? {
            Sign::Positive => rails.swap(0, 1),
            Sign::Zero => return Err(bad("a Bezier band has zero axial width")),
            Sign::Negative => (),
        }
        let ccw = f.sense == Sense::Same;
        for (layer, rail) in rails.iter().enumerate() {
            if ((rail.coedge == Sense::Same) == rail.forward) != ((layer == 0) == ccw) {
                return Err(bad(
                    "a Bezier rail loop disagrees with its chart orientation",
                ));
            }
        }
        let mut found = [None; 2];
        for (c, edge, line) in generators {
            let mut matched = false;
            for (end, present) in found.iter_mut().enumerate() {
                let base = &support.controls()[end * (support.controls().len() - 1)];
                let wanted = rails
                    .iter()
                    .map(|r| base.plus(&support.direction().scaled(&r.height)))
                    .collect::<Vec<_>>();
                let physical = [rails[0].vertices[end], rails[1].vertices[end]];
                if [edge.start, edge.end] != physical.map(Some)
                    && [edge.end, edge.start] != physical.map(Some)
                {
                    continue;
                }
                let (vertices, forward) = ends(model, edge, [&wanted[0], &wanted[1]], dec)?;
                if vertices != physical || present.is_some() {
                    return Err(bad("Bezier endpoint generators are duplicated or unsewn"));
                }
                if sign(&line.dir.norm_squared(), dec)? != Sign::Positive
                    || sign(&line.dir.cross(support.direction()).norm_squared(), dec)? != Sign::Zero
                    || sign(
                        &line
                            .origin
                            .minus(&wanted[0])
                            .cross(support.direction())
                            .norm_squared(),
                        dec,
                    )? != Sign::Zero
                {
                    return Err(bad(
                        "a Bezier endpoint trim is not its exact generator line",
                    ));
                }
                let positive = sign(&line.dir.dot(support.direction()), dec)? == Sign::Positive;
                if (edge.sense == Sense::Same) != (positive == forward)
                    || ((c.sense == Sense::Same) == forward) != ((end == 1) == ccw)
                {
                    return Err(bad("a Bezier generator has an inconsistent orientation"));
                }
                *present = Some((c.edge, c.sense, forward));
                matched = true;
            }
            if !matched {
                return Err(bad("a Bezier generator is not sewn to both complete rails"));
            }
        }
        if !found.iter().all(Option::is_some) {
            return Err(bad("a Bezier endpoint generator is missing"));
        }
        Ok(Some(Self {
            support,
            axial: [rails[0].height.clone(), rails[1].height.clone()],
            rows: [
                (
                    rails[0].edge,
                    BoundaryRow::Rail,
                    rails[0].coedge,
                    rails[0].forward,
                ),
                (
                    rails[1].edge,
                    BoundaryRow::Mouth,
                    rails[1].coedge,
                    rails[1].forward,
                ),
                {
                    let (edge, sense, forward) = found[0].expect("checked endpoint row");
                    (edge, BoundaryRow::Negative, sense, forward)
                },
                {
                    let (edge, sense, forward) = found[1].expect("checked endpoint row");
                    (edge, BoundaryRow::Positive, sense, forward)
                },
            ],
        }))
    }

    /// The checked finite carrier.
    pub fn support(&self) -> &'a BezierExtrusion {
        self.support
    }
    /// Exact ordered axial parameters of the two physical rails.
    pub fn axial(&self) -> [&Real; 2] {
        self.axial.each_ref()
    }
    /// Lower and upper complete polynomial rail edges.
    pub fn rails(&self) -> [EdgeId; 2] {
        [self.rows[0].0, self.rows[1].0]
    }
    /// Complete chart rows, with coedge sense and increasing-parameter direction.
    pub fn rows(&self) -> [(EdgeId, BoundaryRow, Sense, bool); 4] {
        self.rows
    }
    /// A point strictly within the checked finite chart.
    pub fn interior_point(&self, dec: &mut Decider) -> Result<Point3, Refusal> {
        let half = Real::from(Rational::from_ratio(1, 2).expect("nonzero denominator"));
        let parameter = &(&self.axial[0] + &self.axial[1]) * &half;
        let base = self.support.lower_edge().point_at(&half, dec)?;
        Ok(base.plus(&self.support.direction().scaled(&parameter)))
    }

    /// Closed finite-domain membership, with zero on a physical boundary.
    pub fn locate(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        let Some((u, v)) = self.support.parameters_at(point, dec)? else {
            return Ok(Sign::Positive);
        };
        let lo = sign(&(&v - &self.axial[0]), dec)?;
        let hi = sign(&(&v - &self.axial[1]), dec)?;
        if lo == Sign::Negative || hi == Sign::Positive {
            return Ok(Sign::Positive);
        }
        if lo == Sign::Zero
            || hi == Sign::Zero
            || sign(&u, dec)? == Sign::Zero
            || sign(&(&u - &Real::one()), dec)? == Sign::Zero
        {
            Ok(Sign::Zero)
        } else {
            Ok(Sign::Negative)
        }
    }

    /// Checks every Bezier band and the body's complete structural boundary.
    /// Returns false when the body has no band from this exact family.
    pub fn in_body(model: &'a Model, body: BodyId, dec: &mut Decider) -> Result<bool, Refusal> {
        let mut found = false;
        let mut other_spline = false;
        for face in model.faces_of_body(body) {
            match model.geom().surface(model.face(face).surface) {
                Surface::BSpline(SplineSupport::BezierExtrusion(_)) => {
                    Self::of_face(model, face, dec)?;
                    found = true;
                }
                Surface::BSpline(_) => other_spline = true,
                _ => (),
            }
        }
        if found && (other_spline || !model.check(body).is_empty()) {
            return Err(bad("a Bezier band body has an unproved complete boundary"));
        }
        Ok(found)
    }
}
