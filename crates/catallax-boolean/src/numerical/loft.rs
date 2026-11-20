//! Numerical material selection after an exact finite intersection and imprint.
//! The selected faces are stored as a private surface collection, never an exact
//! Boolean answer. The embedding callback must return an admitted numerical body.

use catallax_intersect::numerical_ray::{NumericalRayHit, patch_meetings};
use catallax_numbers::{
    Ball, Decider, MarkLog, NumericalCheck, NumericalEvidence, Rational, Real, Refusal, Sign,
    current_operation_work,
};
use catallax_shapes::{Plane, Point3, RationalPatch, SplineSupport, Surface, Vector3};
use catallax_topology::{
    Body, BodyId, BodyKind, Coedge, Face, FaceId, Loop, Lump, Model, Sense, Shell,
    numerical::{NumericalBody, NumericalOptions},
};

use crate::BoolOp;
use std::collections::BTreeMap;

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("boolean.numerical-loft", detail)
}

fn reversed(sense: Sense) -> Sense {
    if sense == Sense::Same {
        Sense::Reversed
    } else {
        Sense::Same
    }
}

fn outward(plane: &Plane, sense: Sense) -> Vector3 {
    plane
        .normal
        .scaled(&Real::from_i64(if sense == Sense::Same { 1 } else { -1 }))
}

fn zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    let value = dec.enclose_budgeted_at(value, 128)?;
    Ok(value.lo().is_zero() && value.hi().is_zero())
}

#[derive(Clone, Debug)]
struct Loft {
    patches: Vec<RationalPatch>,
    sense: Sense,
    planes: Vec<(Plane, Sense)>,
    ray: Vector3,
}

impl Loft {
    fn original(model: &Model, body: BodyId, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut patches = Vec::new();
        let mut planes = Vec::new();
        let mut sense = None;
        for id in model.faces_of_body(body) {
            dec.charge_work(1)?;
            let face = model.face(id);
            if face.loops.len() != 1 {
                return Err(bad(
                    "the original loft family requires one physical loop per face",
                ));
            }
            match model.geom().surface(face.surface) {
                Surface::Plane(plane) => planes.push((plane.clone(), face.sense)),
                Surface::BSpline(support) if sense.is_none() => {
                    sense = Some(face.sense);
                    match support {
                        SplineSupport::RationalPatch(patch) => patches.push(patch.clone()),
                        SplineSupport::RationalAtlas(atlas) => {
                            dec.charge_work(
                                u64::try_from(atlas.charts().len())
                                    .map_err(|_| bad("too many source charts"))?,
                            )?;
                            patches
                                .extend(atlas.charts().iter().map(|chart| chart.patch().clone()));
                        }
                        _ => return Err(bad("the loft needs retained finite rational charts")),
                    }
                }
                _ => {
                    return Err(bad(
                        "the original loft family needs one spline wall and three planes",
                    ));
                }
            }
        }
        if planes.len() != 3 || patches.is_empty() {
            return Err(bad(
                "the original loft family needs one spline wall and three planes",
            ));
        }
        let mut floor = None;
        for i in 0..3 {
            for j in i + 1..3 {
                dec.charge_work(1)?;
                let a = outward(&planes[i].0, planes[i].1);
                let b = outward(&planes[j].0, planes[j].1);
                if zero(&a.cross(&b).norm_squared(), dec)?
                    && dec.enclose_budgeted_at(&a.dot(&b), 128)?.hi().sign() == Sign::Negative
                {
                    let k = 3 - i - j;
                    let normal = outward(&planes[k].0, planes[k].1);
                    if zero(&normal.dot(&a), dec)? {
                        floor = Some(k);
                    }
                }
            }
        }
        let floor = floor
            .ok_or_else(|| bad("the loft needs opposed caps perpendicular to its planar floor"))?;
        let ray = outward(&planes[floor].0, planes[floor].1).scaled(&Real::from_i64(-1));
        Ok(Self {
            patches,
            sense: sense.unwrap(),
            planes,
            ray,
        })
    }

    fn classify(
        &self,
        point: &Point3,
        tolerance: &Rational,
        dec: &mut Decider,
    ) -> Result<Option<(bool, Vec<NumericalEvidence>)>, Refusal> {
        let Some((inside, plane_evidence)) = planes(&self.planes, point, tolerance, dec)? else {
            return Ok(None);
        };
        if !inside {
            return Ok(Some((false, vec![plane_evidence])));
        }
        let mut reports = vec![plane_evidence];
        let mut hits = Vec::<NumericalRayHit>::new();
        for patch in &self.patches {
            dec.charge_work(1)?;
            let result = patch_meetings(patch, point, &self.ray, tolerance, dec)?;
            hits.extend(result.hits);
            reports.push(result.evidence);
        }
        let count = u64::try_from(hits.len()).map_err(|_| bad("too many numerical ray reports"))?;
        dec.charge_work(count.saturating_mul(count))?;
        hits.sort_by(|a, b| a.distance_mm.cmp(&b.distance_mm));
        let Some(first) = hits.first() else {
            return Ok(Some((false, reports)));
        };
        if first.distance_mm <= *tolerance {
            return Ok(None);
        }
        let direction = |hit: &NumericalRayHit| -> Option<bool> {
            if hit.normal_dot.lo().sign() == Sign::Positive {
                Some(self.sense == Sense::Same)
            } else if hit.normal_dot.hi().sign() == Sign::Negative {
                Some(self.sense == Sense::Reversed)
            } else {
                None
            }
        };
        let Some(inside) = direction(first) else {
            return Ok(None);
        };
        for hit in &hits[1..] {
            dec.charge_work(1)?;
            if &hit.distance_mm - &first.distance_mm > *tolerance {
                break;
            }
            if direction(hit) != Some(inside) {
                return Ok(None);
            }
        }
        Ok(Some((inside, reports)))
    }
}

fn planes(
    planes: &[(Plane, Sense)],
    point: &Point3,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<Option<(bool, NumericalEvidence)>, Refusal> {
    let mut border = false;
    for (plane, sense) in planes {
        dec.charge_work(1)?;
        let normal = outward(plane, *sense);
        let offset = point.minus(&plane.origin).dot(&normal) / normal.norm_squared().sqrt();
        let offset = dec.enclose_budgeted_at(&offset, 128)?;
        if offset.lo() > *tolerance {
            return Ok(Some((
                false,
                plane_evidence(tolerance, "exterior-halfspace")?,
            )));
        }
        if offset.hi() >= -tolerance {
            border = true;
        }
    }
    if border {
        Ok(None)
    } else {
        Ok(Some((
            true,
            plane_evidence(tolerance, "interior-halfspaces")?,
        )))
    }
}

fn plane_evidence(tolerance: &Rational, scope: &str) -> Result<NumericalEvidence, Refusal> {
    NumericalEvidence::new(
        tolerance.clone(),
        vec![NumericalCheck::Passed {
            name: "material.offset-clear-of-preset-tolerance-band".into(),
            scope: scope.into(),
        }],
        vec![],
    )
}

/// One retained numerical material witness, with all successful ray checks.
#[derive(Clone, Debug)]
pub struct NumericalLoftFaceClass {
    /// Original imprinted face identity.
    pub face: FaceId,
    /// True for a loft piece classified against the planar stock.
    pub from_loft: bool,
    /// Actual interior point used by the classifier.
    pub witness: Point3,
    /// Certified coordinate enclosures of that witness, for numerical reporting.
    pub witness_bounds: [Ball; 3],
    /// Numerical inside/outside decision against the other operand.
    pub inside: bool,
    /// Successful finite query and tolerance checks.
    pub evidence: Vec<NumericalEvidence>,
}

/// Numerical face selection with no exact Boolean or solid authority.
#[derive(Clone, Debug)]
pub struct NumericalLoftSelection {
    model: Model,
    collection: BodyId,
    /// Recorded material witnesses for every imprinted face, including discarded faces.
    pub classifications: Vec<NumericalLoftFaceClass>,
    /// Source admission, intersection, imprint, classification and collection work.
    pub work_stages: [u64; 5],
    tolerance_mm: Rational,
}

impl NumericalLoftSelection {
    /// Embeds the private face collection and admits it as a numerical body.
    /// The callback receives a surface collection with no closed-solid claim.
    pub fn admit_with(
        &self,
        dec: &mut Decider,
        embed: impl FnOnce(
            &Model,
            BodyId,
            NumericalOptions,
            &mut Decider,
        ) -> Result<NumericalBody, Refusal>,
    ) -> Result<NumericalBody, Refusal> {
        embed(
            &self.model,
            self.collection,
            NumericalOptions {
                tolerance_mm: self.tolerance_mm.clone(),
                vertex_tolerance_mm: None,
            },
            dec,
        )
    }
}

/// Chooses a top stock spanning the width, strictly between the caps and above
/// the floor, with a roof point strictly inside it. This routing witness gives
/// no Boolean material, incidence or complete overlap authority.
pub fn stock_overlaps_loft_top(
    model: &Model,
    loft: BodyId,
    stock: BodyId,
    tolerance_mm: &Rational,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    if tolerance_mm.sign() != Sign::Positive {
        return Err(bad("top-stock routing needs a positive tolerance"));
    }
    let source = Loft::original(model, loft, dec)?;
    let mut width = None;
    for (plane, sense) in &source.planes {
        let axis = outward(plane, *sense).cross(&source.ray);
        if dec
            .enclose_budgeted_at(&axis.norm_squared(), 128)?
            .lo()
            .sign()
            == Sign::Positive
        {
            width = Some(axis);
            break;
        }
    }
    let axis = width.ok_or_else(|| bad("the loft has no resolved width direction"))?;
    let mut loft_range: Option<Ball> = None;
    for patch in &source.patches {
        for point in patch.carrier().controls() {
            dec.charge_work(1)?;
            let projection = &point.x * &axis.x + &point.y * &axis.y + &point.z * &axis.z;
            let bounds = dec.enclose_budgeted_at(&projection, 128)?;
            loft_range =
                Some(loft_range.map_or_else(|| bounds.clone(), |range| range.hull(&bounds)));
        }
    }
    let mut stock_range: Option<Ball> = None;
    for id in model.vertices_of_body(stock) {
        dec.charge_work(1)?;
        let point = model.geom().position(model.vertex(id).point);
        let projection = &point.x * &axis.x + &point.y * &axis.y + &point.z * &axis.z;
        let bounds = dec.enclose_budgeted_at(&projection, 128)?;
        stock_range = Some(stock_range.map_or_else(|| bounds.clone(), |range| range.hull(&bounds)));
    }
    let spans_width = match (loft_range, stock_range) {
        (Some(loft), Some(stock)) => stock.lo() < loft.lo() && stock.hi() > loft.hi(),
        _ => false,
    };
    if !spans_width {
        return Ok(false);
    }
    for id in model.vertices_of_body(stock) {
        dec.charge_work(1)?;
        let point = model.geom().position(model.vertex(id).point);
        if !matches!(
            planes(&source.planes, &point, tolerance_mm, dec)?,
            Some((true, _))
        ) {
            return Ok(false);
        }
    }
    let mut stock_planes = Vec::new();
    for id in model.faces_of_body(stock) {
        dec.charge_work(1)?;
        let face = model.face(id);
        let Surface::Plane(plane) = model.geom().surface(face.surface) else {
            return Ok(false);
        };
        stock_planes.push((plane.clone(), face.sense));
    }
    if stock_planes.len() != 6 {
        return Ok(false);
    }
    for patch in &source.patches {
        let points = patch.domain().each_ref().map(|range| {
            [
                range[0].clone(),
                (&range[0] + &range[1]).mul_pow2(-1),
                range[1].clone(),
            ]
        });
        for u in &points[0] {
            for v in &points[1] {
                dec.charge_work(1)?;
                let point = patch.point_at(
                    &[
                        Real::from_rational(u.clone()),
                        Real::from_rational(v.clone()),
                    ],
                    dec,
                )?;
                if matches!(
                    planes(&stock_planes, &point, tolerance_mm, dec)?,
                    Some((true, _))
                ) {
                    return Ok(true);
                }
            }
        }
    }
    Ok(false)
}

/// Selects an original loft against an original convex planar stock body.
/// The caller must supply unchanged constructor outputs, before any trimmed-face
/// edits. Numerical labels never enter exact Classification or incidence facts.
pub fn loft_box_selection(
    original: &Model,
    loft: BodyId,
    stock: BodyId,
    operation: BoolOp,
    tolerance_mm: Rational,
    dec: &mut Decider,
) -> Result<NumericalLoftSelection, Refusal> {
    let started = current_operation_work().used;
    if tolerance_mm.sign() != Sign::Positive {
        return Err(bad(
            "the material classification tolerance must be positive",
        ));
    }
    let source = Loft::original(original, loft, dec)?;
    let mut stock_planes = Vec::new();
    for id in original.faces_of_body(stock) {
        dec.charge_work(1)?;
        let face = original.face(id);
        let Surface::Plane(plane) = original.geom().surface(face.surface) else {
            return Err(bad("the stock input needs planar constructor faces"));
        };
        stock_planes.push((plane.clone(), face.sense));
    }
    if stock_planes.len() != 6 {
        return Err(bad("the stock input needs six constructor planes"));
    }
    let mut model = original.clone();
    let start = current_operation_work().used;
    let graph = crate::intersect(&mut model, loft, stock, dec)?.value;
    let after_intersect = current_operation_work().used;
    let imprinted = crate::imprint(&mut model, &graph, dec)?.value;
    let after_imprint = current_operation_work().used;
    let mut classifications = Vec::new();
    let mut selected = Vec::new();
    for (first, pieces) in [
        (true, imprinted.pieces_of_first()),
        (false, imprinted.pieces_of_second()),
    ] {
        for face in pieces {
            dec.charge_work(1)?;
            let class = crate::points::find(&mut model, face, dec, |_, _, point, dec| {
                let class = if first {
                    planes(&stock_planes, point, &tolerance_mm, dec)?
                        .map(|(inside, evidence)| (inside, vec![evidence]))
                } else {
                    source.classify(point, &tolerance_mm, dec)?
                };
                let Some((inside, evidence)) = class else {
                    return Ok(None);
                };
                let witness_bounds = [&point.x, &point.y, &point.z]
                    .map(|value| dec.enclose_budgeted_at(value, 128))
                    .into_iter()
                    .collect::<Result<Vec<_>, _>>()?
                    .try_into()
                    .unwrap();
                Ok(Some(NumericalLoftFaceClass {
                    face,
                    from_loft: first,
                    witness: point.clone(),
                    witness_bounds,
                    inside,
                    evidence,
                }))
            })?;
            let keep = match (first, operation) {
                (true, BoolOp::Unite | BoolOp::Subtract) => !class.inside,
                (true, BoolOp::Intersect) => class.inside,
                (false, BoolOp::Unite) => !class.inside,
                (false, BoolOp::Subtract | BoolOp::Intersect) => class.inside,
            };
            if keep {
                selected.push((face, !first && operation == BoolOp::Subtract));
            }
            classifications.push(class);
        }
    }
    let after_classify = current_operation_work().used;
    let mut faces = Vec::new();
    for (id, reverse) in selected {
        dec.charge_work(1)?;
        if !reverse {
            faces.push(id);
            continue;
        }
        let face = model.face(id).clone();
        let mut loops = Vec::new();
        for lp in face.loops {
            let Loop::Coedges(uses) = model.loop_(lp).clone() else {
                return Err(bad("the selected face needs explicit finite boundary uses"));
            };
            let mut reversed_uses = Vec::new();
            for id in uses.into_iter().rev() {
                dec.charge_work(1)?;
                let coedge = model.coedge(id).clone();
                reversed_uses.push(model.add_coedge(Coedge {
                    edge: coedge.edge,
                    sense: reversed(coedge.sense),
                }));
            }
            loops.push(model.add_loop(Loop::Coedges(reversed_uses)));
        }
        faces.push(model.add_face(Face {
            surface: face.surface,
            sense: reversed(face.sense),
            loops,
        }));
    }
    let mut pairs = BTreeMap::<_, [u64; 2]>::new();
    for id in &faces {
        for lp in &model.face(*id).loops {
            let Loop::Coedges(uses) = model.loop_(*lp) else {
                return Err(bad(
                    "the selected collection needs explicit finite boundary uses",
                ));
            };
            for id in uses {
                dec.charge_work(1)?;
                let coedge = model.coedge(*id);
                pairs.entry(coedge.edge).or_default()
                    [usize::from(coedge.sense == Sense::Reversed)] += 1;
            }
        }
    }
    for counts in pairs.values() {
        dec.charge_work(1)?;
        if *counts != [1, 1] {
            return Err(bad(
                "the numerically selected boundary does not close in opposed edge pairs",
            ));
        }
    }
    if faces.is_empty() {
        return Err(bad(
            "the numerical loft family requires a nonempty selection",
        ));
    }
    dec.charge_work(3)?;
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let collection = model.add_body(Body {
        kind: BodyKind::Surface,
        lumps: vec![lump],
        wires: Vec::new(),
        marks: MarkLog::new(),
    });
    dec.check_deadline()?;
    let finished = current_operation_work().used;
    Ok(NumericalLoftSelection {
        model,
        collection,
        classifications,
        work_stages: [
            start - started,
            after_intersect - start,
            after_imprint - after_intersect,
            after_classify - after_imprint,
            finished - after_classify,
        ],
        tolerance_mm,
    })
}
