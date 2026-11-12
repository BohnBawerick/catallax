//! Certified helpers shared by the stages and point queries.

use catallax_intersect::{Cyclic, curve_is_closed, cyclic_order, on_curve, order_on_curve};
use catallax_numbers::{
    Decider, Fact, FactKind, FactOrigin, NotExact, Outcome, Rational, Real, RealPoly, Refusal,
    Sign, StrictOrder, StrictSign,
};
use catallax_shapes::{Curve, CurveId, GeomKey, PointId, RVector3, Surface, Torus, Vector3};
use catallax_topology::{EdgeId, Model, Sense};

// Geometric search counts are independent of arithmetic working precision.
// Exhausting these searches refuses; it never turns an unresolved sign into
// a geometric answer.
pub(crate) const DIRECTION_ATTEMPTS: u32 = 64;
pub(crate) const PROBE_DEPTHS: u32 = 256;

/// The chart carrier used for membership in a trimmed face. An outer torus
/// sheet has the ring torus's `(u, v)` chart; its trims keep the face where
/// `major + minor cos v > 0`, and a queried point is already on the sheet.
pub(crate) fn chart(surface: &Surface) -> std::borrow::Cow<'_, Surface> {
    match surface {
        Surface::TorusSheet(s) => std::borrow::Cow::Owned(Surface::Torus(Torus {
            center: s.center.clone(),
            axis: s.axis.clone(),
            major: s.major.clone(),
            minor: s.minor.clone(),
            ref_dir: s.ref_dir.clone(),
        })),
        _ => std::borrow::Cow::Borrowed(surface),
    }
}

pub(crate) fn refusal(code: &'static str, detail: impl Into<String>) -> Refusal {
    Refusal::new(code, detail)
}

/// Refuses with a boolean code, keeping the undecided decision and its detail.
pub(crate) fn undecided(code: &'static str, decision: NotExact) -> Refusal {
    Refusal::from(decision).with_code(code)
}

pub(crate) fn sign(x: &Real, dec: &mut Decider, code: &'static str) -> Result<Sign, Refusal> {
    match dec.sign(x) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(n) => Err(undecided(code, n)),
        Outcome::Refused(r) => Err(r),
    }
}

pub(crate) fn coincident<T>(outcome: Outcome<T>, code: &'static str) -> Result<bool, Refusal> {
    match outcome {
        Outcome::Proven(_) => Ok(false),
        Outcome::ProvenCoincident(_) => Ok(true),
        Outcome::NotExact(n) => Err(undecided(code, n)),
        Outcome::Refused(r) => Err(r),
    }
}

pub(crate) fn on_surface(
    model: &Model,
    surface: catallax_shapes::SurfaceId,
    point: PointId,
    dec: &mut Decider,
    code: &'static str,
) -> Result<bool, Refusal> {
    if let Surface::BSpline(catallax_shapes::SplineSupport::RationalAtlas(support)) =
        model.geom().surface(surface)
    {
        let position = model.geom().position(point);
        for chart in support.charts() {
            if chart.patch().parameters_at(&position, dec)?.is_some() {
                return Ok(true);
            }
        }
        return Ok(false);
    }
    if let Surface::BSpline(catallax_shapes::SplineSupport::RationalPatch(support)) =
        model.geom().surface(surface)
    {
        return Ok(support
            .parameters_at(&model.geom().position(point), dec)?
            .is_some());
    }
    if let Surface::BSpline(catallax_shapes::SplineSupport::BezierExtrusion(support)) =
        model.geom().surface(surface)
    {
        return Ok(support
            .parameters_at(&model.geom().position(point), dec)?
            .is_some());
    }
    coincident(model.geom().side(surface, point, dec), code)
}

pub(crate) fn reverse(sense: Sense) -> Sense {
    match sense {
        Sense::Same => Sense::Reversed,
        Sense::Reversed => Sense::Same,
    }
}

pub(crate) fn product(a: Sense, b: Sense) -> Sense {
    if a == b { Sense::Same } else { Sense::Reversed }
}

pub(crate) fn oriented(v: &Vector3, sense: Sense) -> Vector3 {
    if sense == Sense::Same {
        v.clone()
    } else {
        v.scaled(&Real::from_i64(-1))
    }
}

pub(crate) fn natural_ends(model: &Model, edge: EdgeId) -> (Option<PointId>, Option<PointId>) {
    let edge = model.edge(edge);
    let a = edge.start.map(|v| model.vertex(v).point);
    let b = edge.end.map(|v| model.vertex(v).point);
    if edge.sense == Sense::Same {
        (a, b)
    } else {
        (b, a)
    }
}

pub(crate) fn on_arc(
    model: &Model,
    curve: CurveId,
    from: Option<PointId>,
    to: Option<PointId>,
    point: PointId,
    dec: &mut Decider,
    code: &'static str,
) -> Result<bool, Refusal> {
    if !coincident(on_curve(model.geom(), curve, point, dec), code)? {
        return Ok(false);
    }
    let (Some(a), Some(b)) = (from, to) else {
        assert!(from.is_none() && to.is_none());
        return Ok(true);
    };
    if coincident(model.geom().same_point(a, point, dec), code)?
        || coincident(model.geom().same_point(b, point, dec), code)?
    {
        return Ok(true);
    }
    if coincident(model.geom().same_point(a, b, dec), code)? {
        return Ok(true);
    }
    if !curve_is_closed(model.geom(), curve, dec)? {
        let after = match order_on_curve(model.geom(), curve, a, point, dec) {
            Outcome::Proven(StrictOrder::Less) => true,
            Outcome::Proven(_) | Outcome::ProvenCoincident(_) => false,
            Outcome::NotExact(n) => return Err(undecided(code, n)),
            Outcome::Refused(r) => return Err(r),
        };
        let before = match order_on_curve(model.geom(), curve, point, b, dec) {
            Outcome::Proven(StrictOrder::Less) => true,
            Outcome::Proven(_) | Outcome::ProvenCoincident(_) => false,
            Outcome::NotExact(n) => return Err(undecided(code, n)),
            Outcome::Refused(r) => return Err(r),
        };
        Ok(after && before)
    } else {
        match cyclic_order(model.geom(), curve, [a, point, b], dec) {
            Outcome::Proven(Cyclic::Forward) => Ok(true),
            Outcome::Proven(Cyclic::Backward) => Ok(false),
            Outcome::ProvenCoincident(_) => Ok(true),
            Outcome::NotExact(n) => Err(undecided(code, n)),
            Outcome::Refused(r) => Err(r),
        }
    }
}

/// Runs `work` on a copy of the model, then does to the model, in order,
/// what the work did to its copy, and returns the work's value with its new
/// handles where the model now has them.
///
/// The model ends as if the work had run on it: every entity, fact and name
/// the work added is added again through the same methods, in the same
/// order. This is the form a loop item must have before another thread can
/// run it, and the request's own thread runs its items this way too, so that
/// one thread and many do the same thing to the model.
pub(crate) fn on_a_copy<T: catallax_topology::Relocate>(
    model: &mut Model,
    work: impl FnOnce(&mut Model) -> T,
) -> T {
    let mut copy = model.clone();
    copy.journal_begin();
    let value = work(&mut copy);
    let journal = copy.journal_end();
    // The copy goes first: while it lives it shares the model's storage, and
    // every table the replay adds to would be copied before it is written.
    drop(copy);
    let remap = model.replay(&journal);
    value.relocate(&remap)
}

/// What another thread hands back for a loop item that adds to the model:
/// the item's value, what the item added to its own copy of the model, and
/// what it asked that copy that later additions can change the answer to.
pub(crate) struct Planned<T> {
    value: T,
    journal: catallax_topology::Journal,
    reads: catallax_topology::Reads,
    // The model the copy was taken from.
    base: std::sync::Arc<Model>,
}

/// Runs a loop item on a copy of `base`, on a thread that is ahead of the
/// writer. The writer [`Planned::commit`]s the result when it reaches the
/// item, if it [`Planned::holds`] then.
pub(crate) fn plan<T>(
    base: &std::sync::Arc<Model>,
    work: impl FnOnce(&mut Model) -> T,
) -> Planned<T> {
    let mut copy = Model::clone(base);
    copy.journal_begin();
    copy.watch();
    let value = work(&mut copy);
    let reads = copy.watched();
    let journal = copy.journal_end();
    Planned {
        value,
        journal,
        reads,
        base: base.clone(),
    }
}

impl<T: catallax_topology::Relocate> Planned<T> {
    /// The item's value, with the handles its copy gave what it made.
    pub(crate) fn value(&self) -> &T {
        &self.value
    }

    /// Whether the item would have done on `model` what it did on its copy.
    /// `model` has grown from the item's base by what earlier items added.
    /// The item would have done the same when everything it asked its copy
    /// has the same answer on `model`.
    pub(crate) fn holds(&self, model: &Model) -> bool {
        model.unchanged_where_read(&self.base, &self.reads)
    }

    /// Does to `model`, in order, what the item did to its copy, and returns
    /// the item's value with its new handles where `model` now has them.
    /// Only for a result that [`Planned::holds`] on `model`.
    pub(crate) fn commit(self, model: &mut Model) -> T {
        let remap = model.replay(&self.journal);
        self.value.relocate(&remap)
    }
}

pub(crate) fn record_incident(model: &mut Model, a: GeomKey, b: GeomKey) {
    if !model.geom().facts().holds(FactKind::Incident, a, b) {
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Incident,
            a,
            b,
            at: None,
            origin: FactOrigin::Construction("boolean.sew"),
        });
    }
}

pub(crate) fn record_same_point(model: &mut Model, a: PointId, b: PointId) {
    if a != b {
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Coincident,
            a: GeomKey::Point(a),
            b: GeomKey::Point(b),
            at: None,
            origin: FactOrigin::Proof(catallax_numbers::ZeroProof::Algebraic),
        });
    }
}

/// Whether an existing exact relation proves this equation zero at the source.
/// This says nothing about a different selected component or the rest of a path.
pub(crate) fn equation_is_source_contact(
    model: &Model,
    source: PointId,
    source_curve: CurveId,
    barrier: CurveId,
    equation_index: usize,
    equation: &RealPoly,
    source_equations: &[RealPoly],
) -> bool {
    let facts = model.geom().facts();
    if facts.holds(
        FactKind::Coincident,
        GeomKey::Curve(source_curve),
        GeomKey::Curve(barrier),
    ) || facts.holds(
        FactKind::Incident,
        GeomKey::Point(source),
        GeomKey::Curve(barrier),
    ) || source_equations.iter().any(|e| e.same_definition(equation))
    {
        return true;
    }
    if let Curve::Intersection(branch) = model.geom().curve(barrier) {
        let surface = GeomKey::Surface(branch.surfaces[equation_index]);
        if facts.holds(FactKind::Incident, GeomKey::Point(source), surface)
            || facts.holds(FactKind::Incident, GeomKey::Curve(source_curve), surface)
        {
            return true;
        }
    }
    let opposite = -equation;
    source_equations
        .iter()
        .any(|e| e.same_definition(&opposite))
}

pub(crate) fn canonical_point(
    model: &mut Model,
    points: &mut Vec<PointId>,
    point: PointId,
    dec: &mut Decider,
    code: &'static str,
) -> Result<PointId, Refusal> {
    // The first equal point ends the search. A comparison is a value of what
    // every thread shares, so other threads make later ones whole and this
    // one adopts them in order. Merging a point into the shared list is an
    // ordered append.
    #[cfg(feature = "work-span")]
    let _append = catallax_numbers::work_span::append();
    #[cfg(feature = "work-span")]
    let span_tasks = catallax_numbers::work_span::search("point merge");
    let mut ahead = catallax_numbers::values_ahead!(dec, points.len());
    for (index, &existing) in points.iter().enumerate() {
        #[cfg(feature = "work-span")]
        let _span_task = span_tasks.task();
        let adopted = ahead.take(index, dec, || {
            let (model, points) = (model.clone(), points.clone());
            move |index: usize, dec: &mut Decider| {
                same_canonical_point(&model, points[index], point, dec, code)
            }
        });
        let same = match adopted {
            Some(same) => same,
            None => same_canonical_point(model, existing, point, dec, code),
        }?;
        if same {
            drop(ahead);
            record_same_point(model, existing, point);
            return Ok(existing);
        }
    }
    drop(ahead);
    points.push(point);
    Ok(point)
}

/// Whether an incoming point is the canonical point `existing`.
fn same_canonical_point(
    model: &Model,
    existing: PointId,
    point: PointId,
    dec: &mut Decider,
    code: &'static str,
) -> Result<bool, Refusal> {
    let identity = if model.geom().facts().holds(
        FactKind::Coincident,
        GeomKey::Point(existing),
        GeomKey::Point(point),
    ) {
        None
    } else {
        crate::point_identity::same_cylinder_point(model.geom(), existing, point, dec)?
    };
    match identity {
        Some(equal) => Ok(equal),
        None => coincident(model.geom().same_point(existing, point, dec), code),
    }
}

pub(crate) fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
pub(crate) fn unit(axis: usize) -> Vector3 {
    RVector3::new(
        r(i64::from(axis == 0)),
        r(i64::from(axis == 1)),
        r(i64::from(axis == 2)),
    )
    .to_real()
}

pub(crate) fn perpendicular(
    vector: &Vector3,
    dec: &mut Decider,
    code: &'static str,
) -> Result<Vector3, Refusal> {
    for axis in 0..3 {
        let direction = vector.cross(&unit(axis));
        if sign(&direction.norm_squared(), dec, code)? != Sign::Zero {
            return Ok(direction);
        }
    }
    Err(refusal(code, "zero direction has no transverse axis"))
}

/// A diagonal enumeration of two signed Calkin-Wilf sequences. Its rational
/// pairs are dense, unlike a single polynomial family of ray directions.
pub(crate) fn probe_parameters(index: u32) -> (Rational, Rational) {
    fn fraction(index: u64) -> Rational {
        if index == 0 {
            return Rational::zero();
        }
        let positive = index.div_ceil(2);
        let mut numerator = 1i64;
        let mut denominator = 1i64;
        let highest = 63 - positive.leading_zeros();
        for bit in (0..highest).rev() {
            if positive & (1u64 << bit) == 0 {
                denominator += numerator;
            } else {
                numerator += denominator;
            }
        }
        if index.is_multiple_of(2) {
            numerator = -numerator;
        }
        Rational::from_ratio(numerator, denominator).expect("positive Calkin-Wilf denominator")
    }
    let index = u64::from(index);
    let diagonal = ((index * 8 + 1).isqrt() - 1) / 2;
    let left = index - diagonal * (diagonal + 1) / 2;
    (fraction(left), fraction(diagonal - left))
}

#[cfg(test)]
mod contact_fact_tests {
    use super::*;
    use catallax_shapes::{Circle, Line, Plane, Point3, PointDef, Surface};

    fn fixture() -> (Model, PointId, CurveId, CurveId) {
        let mut model = Model::new();
        let mut dec = Decider::default();
        let source = model.geom_mut().add_curve(Curve::Line(
            Line::new(Point3::origin(), unit(0), &mut dec).unwrap(),
        ));
        let barrier = model.geom_mut().add_curve(Curve::Circle(
            Circle::new(Point3::origin(), unit(2), Real::one(), unit(0), &mut dec).unwrap(),
        ));
        let point = model.geom_mut().add_point(PointDef::Explicit(Point3::new(
            Real::one(),
            Real::zero(),
            Real::zero(),
        )));
        (model, point, source, barrier)
    }
    #[test]
    fn direct_point_curve_fact_retains_analytic_contact() {
        let (mut model, point, source, barrier) = fixture();
        let source_equations = model.geom().curve(source).equations(model.geom()).unwrap();
        let equations = model.geom().curve(barrier).equations(model.geom()).unwrap();
        assert!(!equation_is_source_contact(
            &model,
            point,
            source,
            barrier,
            1,
            &equations[1],
            &source_equations
        ));
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Point(point),
            b: GeomKey::Curve(barrier),
            at: None,
            origin: FactOrigin::Proof(catallax_numbers::ZeroProof::ExactArithmetic),
        });
        assert!(equation_is_source_contact(
            &model,
            point,
            source,
            barrier,
            1,
            &equations[1],
            &source_equations
        ));
    }
    #[test]
    fn a_shared_plane_does_not_prove_the_other_circle_equation() {
        let (mut model, _, source, barrier) = fixture();
        let point = model
            .geom_mut()
            .add_point(PointDef::Explicit(Point3::origin()));
        let plane = model.geom_mut().add_surface(Surface::Plane(
            Plane::new(Point3::origin(), unit(2), unit(0), &mut Decider::default()).unwrap(),
        ));
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Point(point),
            b: GeomKey::Surface(plane),
            at: None,
            origin: FactOrigin::Proof(catallax_numbers::ZeroProof::ExactArithmetic),
        });
        let source_equations = model.geom().curve(source).equations(model.geom()).unwrap();
        let equations = model.geom().curve(barrier).equations(model.geom()).unwrap();
        assert!(!equation_is_source_contact(
            &model,
            point,
            source,
            barrier,
            1,
            &equations[1],
            &source_equations
        ));
    }
    #[test]
    fn opposite_source_equation_retains_the_same_zero_contact() {
        let (model, point, source, barrier) = fixture();
        let equations = model.geom().curve(barrier).equations(model.geom()).unwrap();
        let opposite = -&equations[0];
        assert!(equation_is_source_contact(
            &model,
            point,
            source,
            barrier,
            0,
            &equations[0],
            &[opposite]
        ));
    }
    #[test]
    fn coincident_curve_fact_retains_contact_under_a_different_definition() {
        let (mut model, point, source, _) = fixture();
        let barrier = model.geom_mut().add_curve(Curve::Line(
            Line::new(
                Point3::new(Real::from_i64(2), Real::zero(), Real::zero()),
                unit(0).scaled(&Real::from_i64(2)),
                &mut Decider::default(),
            )
            .unwrap(),
        ));
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Coincident,
            a: GeomKey::Curve(source),
            b: GeomKey::Curve(barrier),
            at: None,
            origin: FactOrigin::Proof(catallax_numbers::ZeroProof::ExactArithmetic),
        });
        let source_equations = model.geom().curve(source).equations(model.geom()).unwrap();
        let equations = model.geom().curve(barrier).equations(model.geom()).unwrap();
        for (index, equation) in equations.iter().enumerate() {
            assert!(equation_is_source_contact(
                &model,
                point,
                source,
                barrier,
                index,
                equation,
                &source_equations
            ));
        }
    }
}

/// A starved decision must reach the caller with its code and bound intact,
/// whatever code this wrapper reports. See docs/NOT-EXACT-INVENTORY.md.
#[cfg(test)]
mod undecided_evidence {
    use super::*;
    use catallax_numbers::{BudgetPolicy, Decider, NotExact, Outcome, Real};

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
    fn boolean_sign_keeps_the_decision() {
        let value = Real::from(2).sqrt();
        let decision = undecided(&value);
        let refusal = sign(&value, &mut starved(), "boolean.classify-ambiguous").unwrap_err();
        assert_eq!(refusal.code, "boolean.classify-ambiguous");
        assert_eq!(refusal.undecided(), Some(&decision));
    }

    #[test]
    fn boolean_sameness_keeps_the_decision() {
        let value = Real::from(2).sqrt();
        let decision = undecided(&value);
        let refusal = coincident(
            Outcome::<()>::NotExact(decision.clone()),
            "boolean.merge-ambiguous",
        )
        .unwrap_err();
        assert_eq!(refusal.code, "boolean.merge-ambiguous");
        assert_eq!(refusal.undecided(), Some(&decision));
    }
}
