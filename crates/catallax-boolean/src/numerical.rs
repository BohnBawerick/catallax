//! Open numerical face partitions, separate from exact Boolean incidence.

use catallax_intersect::{NumericalNativeIntersection, NumericalNativeWork};
use catallax_numbers::{
    Ball, Decider, NumericalCheck, NumericalEvidence, NumericalSample, Rational, Real, Refusal,
    Sign,
};
use catallax_shapes::{
    BSplineCurve, BSplineSurface, Circle, Cone, Curve, Cylinder, Point3, RationalPatch,
    SplineSupport, Surface, Vector3,
};
use catallax_topology::{
    CoedgeId, EdgeId, FaceId, LoopId, Model, Sense,
    numerical::{NumericalDefinition, NumericalNativeTrim, NumericalOptions},
};

mod band;
mod checks;
mod events;
mod material;
mod partition;
#[cfg(test)]
mod tests;

const BITS: u32 = 128;

/// One original complete circular ring and its use on the input cone face.
#[derive(Clone, Debug)]
pub struct NumericalConeRing {
    /// Original input loop, coedge and edge identities.
    pub original_loop: LoopId,
    /// Original edge-use identity.
    pub coedge: CoedgeId,
    /// Original edge identity, for sharing with adjacent cap/wall assembly.
    pub edge: EdgeId,
    /// Original traversal relative to its circle carrier.
    pub sense: Sense,
    /// Original edge direction relative to the circle carrier.
    pub edge_sense: Sense,
    /// Original analytic definition, including its parameter-origin direction.
    pub circle: Circle,
    /// Height used by the retained finite cone source.
    pub height_mm: Rational,
    /// Enclosure of the original circle centre's axial height.
    pub original_height_mm: [Rational; 2],
}

/// A retained ring segment linked to the original analytic boundary.
#[derive(Clone, Debug)]
pub struct NumericalRingLink {
    /// Local edge index in the open definition.
    pub edge: usize,
    /// Index in the partition's original ring table.
    pub ring: usize,
    /// The retained cone source's native angular coordinates in edge order.
    pub source_u: [Rational; 2],
    /// Original circle angles in radians, enclosed separately at both ends.
    /// These are unwrapped in the source traversal, so a periodic seam stays explicit.
    pub original_parameters_rad: [[Rational; 2]; 2],
}

/// One occurrence of a finite meeting, with its original native addresses.
#[derive(Clone, Debug)]
pub struct NumericalFiniteHit {
    /// Index of the retained intersection piece.
    pub piece: usize,
    /// Original native parameter of that piece.
    pub parameter: Rational,
    /// Retained cone and cylinder coordinates at that parameter.
    pub uv: [[Rational; 2]; 2],
    /// Requested rectangles for the producing search, in cone/cylinder order.
    pub requested_domains: [[[Rational; 2]; 2]; 2],
}

/// One physical finite event, including all duplicate quarter reports.
#[derive(Clone, Debug)]
pub struct NumericalFiniteEvent {
    /// Stable identity within this input face, derived from boundary ownership.
    pub key: String,
    /// Local vertex to share during assembly, if the event enters a retained loop.
    pub vertex: Option<usize>,
    /// True for a point contact with no nonconstant fitted incident arc.
    pub point_only: bool,
    /// Every reported native address, without replacing the retained definitions.
    pub hits: Vec<NumericalFiniteHit>,
}

/// A retained tool-side UV definition for the caller's cylinder wall assembly.
#[derive(Clone, Debug)]
pub struct NumericalCylinderTrim {
    /// Retained intersection piece and local shared 3D edge.
    pub piece: usize,
    /// Index in the partition's edge table.
    pub edge: usize,
    /// Canonical cylinder quarter ID, independent of the cone face key.
    pub quarter: usize,
    /// The unchanged native UV definition, before a final wall-use address exists.
    pub curve: BSplineCurve,
}

/// Work measured by this call, excluding the caller's recipe prefix.
#[derive(Clone, Copy, Debug)]
pub struct NumericalFaceWork {
    /// Reading and checking original rings and rounding finite heights.
    pub boundaries: u64,
    /// Source construction, tracing and independent retained-fit checks.
    pub intersection: NumericalNativeWork,
    /// Retained cylinder source construction and cone/cylinder chart admission.
    pub supports: u64,
    /// Finite ownership, explicit endpoint joins and kept loop construction.
    pub topology: u64,
}

/// An open face partition. The caller must assemble and admit the complete shell.
/// All geometry and checks remain numerical; this grants no exact incidence or
/// STEP authority. Original sources, unused contacts and actual checks remain
/// available alongside the open topology for save/replay and final reporting.
#[derive(Clone, Debug)]
pub struct NumericalFacePartition {
    /// Vertices, edges and retained cone regions with the input face sense.
    pub definition: NumericalDefinition,
    /// Native UV curves bound to stable cone face/loop/use addresses.
    pub native_trims: Vec<NumericalNativeTrim>,
    /// Every original 3D/UV fit, source definition and independent producer check.
    pub intersection: NumericalNativeIntersection,
    /// Original ring ownership and native interval links for adjacent faces.
    pub rings: Vec<NumericalConeRing>,
    /// Every retained segment of an original ring.
    pub ring_links: Vec<NumericalRingLink>,
    /// Finite meetings, including duplicate seam contacts and explicit vertex ownership.
    pub events: Vec<NumericalFiniteEvent>,
    /// Canonical cylinder quarters for the caller's wall assembly.
    pub cylinder_source: NumericalCylinderSource,
    /// Unbound tool-side UV templates, retaining their original native parameter.
    pub cylinder_trims: Vec<NumericalCylinderTrim>,
    /// Preset tolerance, actual sampled checks and explicit endpoint repairs.
    pub evidence: NumericalEvidence,
    /// Actual work per stage, under the caller's allowance and deadline.
    pub work: NumericalFaceWork,
}

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("boolean.numerical-cone-trim", detail)
}
fn q(value: i64) -> Rational {
    Rational::from_i64(value)
}
fn interval(value: &Real, dec: &mut Decider) -> Result<Ball, Refusal> {
    dec.enclose_budgeted_at(value, BITS)
}
fn unit(value: &Vector3, dec: &mut Decider) -> Result<Vector3, Refusal> {
    let length = value.norm_squared().sqrt();
    if interval(&length, dec)?.lo().sign() != Sign::Positive {
        return Err(bad("a retained frame direction vanished"));
    }
    Ok(value.scaled(&(&Real::one() / &length)))
}
fn upper(value: &Ball) -> Rational {
    value.lo().abs().max(value.hi().abs())
}
fn distance(a: &[Ball; 3], b: &[Ball; 3]) -> Result<Rational, Refusal> {
    let d: [Ball; 3] = std::array::from_fn(|i| &a[i] - &b[i]);
    Ok((d[0].square() + d[1].square() + d[2].square())
        .sqrt()
        .ok_or_else(|| bad("endpoint distance is unresolved"))?
        .hi())
}
fn patch(source: &BSplineSurface, quarter: usize, dec: &mut Decider) -> Result<Surface, Refusal> {
    let knots = source.knots();
    let span = source
        .span_indices()
        .nth(quarter)
        .ok_or_else(|| bad("missing retained quarter"))?;
    let bounds =
        std::array::from_fn(|axis| [&knots[axis][span[axis]], &knots[axis][span[axis] + 1]]);
    let retained = source.bezier_patch(span, bounds, dec)?;
    Ok(Surface::BSpline(SplineSupport::RationalPatch(
        RationalPatch::new(retained, dec)?,
    )))
}

/// Clips one complete circular cone band against a finite cylindrical tool.
/// The tool interval is signed millimetres along its unit axis from its original
/// origin. Reuse the same tool and interval for both cone faces. `None` means an
/// unsupported input face family. A supported separated or point-only face
/// returns `Some` with its retained region and actual numerical evidence.
/// The result is open topology; final cap/wall assembly and admission are required.
pub fn clip_cone_cylinder_face(
    model: &Model,
    cone_face: FaceId,
    face_key: &str,
    tool: &Cylinder,
    tool_axis_parameters_mm: [Rational; 2],
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<Option<NumericalFacePartition>, Refusal> {
    let run = |dec: &mut Decider| {
        build_face(
            model,
            cone_face,
            face_key,
            tool,
            tool_axis_parameters_mm,
            options,
            dec,
        )
    };
    if catallax_numbers::operation_budget_active() {
        return run(dec);
    }
    let (result, work) = catallax_numbers::with_operation_budget(20_000_000, || run(dec));
    if work.exhausted {
        return Err(Refusal::new(
            "numbers.operation-budget",
            "finite numerical face partition exhausted its operation allowance",
        ));
    }
    result
}

fn build_face(
    model: &Model,
    face: FaceId,
    key: &str,
    tool: &Cylinder,
    tool_parameters: [Rational; 2],
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<Option<NumericalFacePartition>, Refusal> {
    let start = catallax_numbers::current_operation_work().used;
    let Surface::Cone(cone) = model.geom().surface(model.face(face).surface) else {
        return Ok(None);
    };
    if key.is_empty() || options.tolerance_mm.sign() != Sign::Positive {
        return Err(bad("face key and positive preset tolerance are required"));
    }
    let Some(band) = band::read(model, face, cone, &options.tolerance_mm, dec)? else {
        return Ok(None);
    };
    let rings = band.rings;
    let mut checks = band.checks;
    let heights = [rings[0].height_mm.clone(), rings[1].height_mm.clone()];
    let before_intersection = catallax_numbers::current_operation_work().used;
    let originals = [Surface::Cone(cone.clone()), Surface::Cylinder(tool.clone())];
    let intersection = catallax_intersect::numerical_native_finite_intersection(
        [&originals[0], &originals[1]],
        heights,
        tool_parameters.clone(),
        dec,
    )?
    .expect("cone/cylinder family checked above");
    checks.extend_from_slice(intersection.evidence().checks());
    let before_charts = catallax_numbers::current_operation_work().used;
    let cone_quarters = (0..4)
        .map(|i| patch(&intersection.windows().sources()[0], i, dec))
        .collect::<Result<Vec<_>, _>>()?;
    let cylinder_source =
        NumericalCylinderSource::new(tool, tool_parameters, options.clone(), dec)?;
    if !cylinder_source
        .source
        .same_definition(&intersection.windows().sources()[1])
    {
        return Err(bad(
            "canonical cylinder constructor changed the finite query source definition",
        ));
    }
    checks.extend(
        cylinder_source
            .evidence
            .checks()
            .iter()
            .cloned()
            .map(|mut check| {
                let (NumericalCheck::Passed { scope, .. } | NumericalCheck::Samples { scope, .. }) =
                    &mut check;
                *scope = format!("canonical-cylinder/{scope}");
                check
            }),
    );
    let before_topology = catallax_numbers::current_operation_work().used;
    let mut result = partition::Input {
        key,
        sense: model.face(face).sense,
        tool,
        intersection,
        rings,
        cone_quarters,
        cylinder_source,
        options,
        checks,
    }
    .build(dec)?;
    result.work.boundaries = before_intersection - start;
    result.work.supports = before_topology - before_charts;
    result.work.topology = catallax_numbers::current_operation_work().used - before_topology;
    Ok(Some(result))
}

/// A canonical cylinder chart, independent of a cone face or assembly key.
#[derive(Clone, Debug)]
pub struct NumericalCylinderQuarter {
    /// Native quadrant, 0 through 3.
    id: usize,
    /// Deterministic descriptor-local key. The caller supplies assembly identities.
    key: String,
    /// Retained embedded patch over native u=[id,id+1], v=[0,1].
    surface: Surface,
    /// Signed millimetres at native v=0 and v=1 from the original tool origin.
    axis_parameters_mm: [Rational; 2],
}

impl NumericalCylinderQuarter {
    /// Native quadrant, 0 through 3.
    pub fn id(&self) -> usize {
        self.id
    }
    /// Descriptor-local key; callers provide their own assembly identity.
    pub fn key(&self) -> &str {
        &self.key
    }
    /// Retained patch over native u=[id,id+1], v=[0,1].
    pub fn surface(&self) -> &Surface {
        &self.surface
    }
    /// Signed millimetres at native v=0 and v=1 about the original origin.
    pub fn axis_parameters_mm(&self) -> &[Rational; 2] {
        &self.axis_parameters_mm
    }
}

/// A retained cylinder source and all four native quarter patches.
/// Constructing it supplies numerical source/control checks and patch embedding
/// checks. It does not construct a Boolean, admitted body or STEP result.
/// Construction keeps the retained definition, interval, quarters and evidence
/// together. Read-only accessors preserve that association.
///
/// ```compile_fail
/// use catallax_boolean::numerical::NumericalCylinderSource;
/// use catallax_numbers::Rational;
/// fn change(mut source: NumericalCylinderSource) {
///     source.axis_parameters_mm[1] = Rational::one();
/// }
/// ```
/// ```compile_fail
/// use catallax_boolean::numerical::NumericalCylinderSource;
/// fn change(mut source: NumericalCylinderSource) { source.quarters.clear(); }
/// ```
/// ```compile_fail
/// use catallax_boolean::numerical::NumericalCylinderSource;
/// fn change(mut source: NumericalCylinderSource) {
///     source.evidence = source.evidence().clone();
/// }
/// ```
/// ```compile_fail
/// use catallax_boolean::numerical::NumericalCylinderQuarter;
/// fn change(mut quarter: NumericalCylinderQuarter) { quarter.key.clear(); }
/// ```
#[derive(Clone, Debug)]
pub struct NumericalCylinderSource {
    /// Original analytic carrier, with the caller's origin, axis and seam.
    original: Cylinder,
    /// Complete finite numerical definition, without native knot changes.
    source: BSplineSurface,
    /// Signed millimetres at native v=0 and v=1 from the original origin.
    axis_parameters_mm: [Rational; 2],
    /// Canonical quarters in increasing native-u order.
    quarters: Vec<NumericalCylinderQuarter>,
    /// Preset tolerance and actual source/control/chart checks.
    evidence: NumericalEvidence,
    /// Actual work for source construction and its checks.
    source_work: u64,
    /// Actual work for admitting every retained quarter patch.
    chart_work: u64,
}
impl NumericalCylinderSource {
    /// Creates the same cylinder source used by finite cone intersections.
    /// Coordinates are signed millimetres along the unit cylinder axis from its
    /// original origin. Reuse the same carrier and interval across all face calls.
    pub fn new(
        cylinder: &Cylinder,
        axis_parameters_mm: [Rational; 2],
        options: NumericalOptions,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let run = |dec: &mut Decider| build(cylinder, axis_parameters_mm, options, dec);
        if catallax_numbers::operation_budget_active() {
            return run(dec);
        }
        let (result, work) = catallax_numbers::with_operation_budget(20_000_000, || run(dec));
        if work.exhausted {
            return Err(Refusal::new(
                "numbers.operation-budget",
                "numerical cylinder quarter construction exhausted its operation allowance",
            ));
        }
        result
    }
    /// Original analytic cylinder, including its authored origin and seam.
    pub fn original(&self) -> &Cylinder {
        &self.original
    }
    /// Complete retained finite source, without changing native knots.
    pub fn source(&self) -> &BSplineSurface {
        &self.source
    }
    /// Signed millimetres at native v=0 and v=1 about the original origin.
    pub fn axis_parameters_mm(&self) -> &[Rational; 2] {
        &self.axis_parameters_mm
    }
    /// All four retained quarters in increasing native-u order.
    pub fn quarters(&self) -> &[NumericalCylinderQuarter] {
        &self.quarters
    }
    /// Actual source/control/chart checks at the preset tolerance.
    pub fn evidence(&self) -> &NumericalEvidence {
        &self.evidence
    }
    /// Measured source construction and check work.
    pub fn source_work(&self) -> u64 {
        self.source_work
    }
    /// Measured retained quarter admission and check work.
    pub fn chart_work(&self) -> u64 {
        self.chart_work
    }
    /// Maps a native v in [0,1] to signed millimetres about the original origin.
    pub fn signed_mm_at_v(&self, v: &Rational) -> Result<Rational, Refusal> {
        if v < &Rational::zero() || v > &Rational::one() {
            return Err(Refusal::new(
                "boolean.numerical-cylinder-source",
                "native cylinder v leaves [0,1]",
            ));
        }
        let parameters = &self.axis_parameters_mm;
        Ok(&parameters[0] + &((&parameters[1] - &parameters[0]) * v.clone()))
    }
}
fn build(
    cylinder: &Cylinder,
    parameters: [Rational; 2],
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<NumericalCylinderSource, Refusal> {
    let before = catallax_numbers::current_operation_work().used;
    let retained = catallax_intersect::numerical_native_cylinder_source(
        cylinder,
        parameters.clone(),
        options.tolerance_mm,
        dec,
    )?;
    let after_source = catallax_numbers::current_operation_work().used;
    let source = retained.source();
    let knots = source.knots();
    let mut quarters = Vec::new();
    let mut checks = retained.evidence().checks().to_vec();
    for (id, span) in source.span_indices().enumerate() {
        dec.charge_work(1)?;
        let domain =
            std::array::from_fn(|axis| [&knots[axis][span[axis]], &knots[axis][span[axis] + 1]]);
        let carrier = source.bezier_patch(span, domain, dec)?;
        let patch = RationalPatch::new(carrier.clone(), dec)?;
        if !patch.carrier().same_definition(&carrier) {
            return Err(Refusal::new(
                "boolean.numerical-cylinder-source",
                "quarter admission changed its retained native carrier",
            ));
        }
        let key = format!("cylinder/quarter-{id}");
        checks.push(NumericalCheck::Passed {
            name: "native.cylinder-retained-quarter-embedding".into(),
            scope: key.clone(),
        });
        quarters.push(NumericalCylinderQuarter {
            id,
            key,
            surface: Surface::BSpline(SplineSupport::RationalPatch(patch)),
            axis_parameters_mm: parameters.clone(),
        });
    }
    if quarters.len() != 4 {
        return Err(Refusal::new(
            "boolean.numerical-cylinder-source",
            "retained cylinder requires exactly four native quarters",
        ));
    }
    dec.charge_work(checks.len() as u64)?;
    let evidence =
        NumericalEvidence::new(retained.evidence().tolerance_mm().clone(), checks, vec![])?;
    dec.check_deadline()?;
    let finished = catallax_numbers::current_operation_work().used;
    Ok(NumericalCylinderSource {
        original: retained.original().clone(),
        source: source.clone(),
        axis_parameters_mm: parameters,
        quarters,
        evidence,
        source_work: after_source - before,
        chart_work: finished - after_source,
    })
}

#[cfg(test)]
mod cylinder_tests {
    use super::*;
    use catallax_numbers::{Real, with_operation_budget};
    use catallax_shapes::{Cone, Point3, Vector3};
    fn q(i: i64) -> Rational {
        Rational::from_i64(i)
    }
    fn options() -> NumericalOptions {
        NumericalOptions {
            tolerance_mm: Rational::from_ratio(1, 10_000_000).unwrap(),
        }
    }
    fn cylinder(offset: &str, skew: bool) -> Cylinder {
        let axis = if skew {
            Vector3::new(1.into(), (-1).into(), 1.into())
                .scaled(&(&Real::one() / &Real::from(3).sqrt()))
        } else {
            Vector3::new(0.into(), (-1).into(), 0.into())
        };
        Cylinder {
            origin: Point3::new(Rational::parse(offset).unwrap().into(), 12.into(), 0.into()),
            axis,
            radius_sq: 1.into(),
            ref_dir: Vector3::new(1.into(), 0.into(), 0.into()),
        }
    }
    #[test]
    fn standalone_cylinder_sources_match_finite_queries_in_native_parameters() {
        for (offset, skew, extent) in [
            ("10", false, "49"),
            ("10.5", false, "99/2"),
            ("10.5", true, "99/2"),
        ] {
            let cylinder = cylinder(offset, skew);
            let extent = Rational::parse(extent).unwrap();
            let parameters = [-&extent, extent];
            let ((), work) = with_operation_budget(20_000_000, || {
                let mut dec = Decider::default();
                let source = NumericalCylinderSource::new(
                    &cylinder,
                    parameters.clone(),
                    options(),
                    &mut dec,
                )
                .unwrap();
                assert_eq!(
                    source.source_work + source.chart_work,
                    catallax_numbers::current_operation_work().used
                );
                assert_eq!(source.signed_mm_at_v(&q(0)).unwrap(), parameters[0]);
                assert_eq!(source.signed_mm_at_v(&q(1)).unwrap(), parameters[1]);
                assert_eq!(
                    source
                        .signed_mm_at_v(&Rational::from_ratio(1, 2).unwrap())
                        .unwrap(),
                    q(0)
                );
                assert!(source.signed_mm_at_v(&q(-1)).is_err());
                assert!(source.signed_mm_at_v(&q(2)).is_err());
                for upper in [false, true] {
                    let cone = Surface::Cone(Cone {
                        apex: Point3::new(
                            0.into(),
                            if upper { 23.into() } else { (-11).into() },
                            0.into(),
                        ),
                        axis: Vector3::new(
                            0.into(),
                            if upper { (-1).into() } else { 1.into() },
                            0.into(),
                        ),
                        slope: 1.into(),
                        ref_dir: Vector3::new(1.into(), 0.into(), 0.into()),
                    });
                    let tool = Surface::Cylinder(cylinder.clone());
                    let query = catallax_intersect::numerical_native_finite_intersection(
                        [&cone, &tool],
                        [q(11), q(12)],
                        parameters.clone(),
                        &mut dec,
                    )
                    .unwrap()
                    .unwrap();
                    assert!(source.source.same_definition(&query.windows().sources()[1]));
                }
                for (id, quarter) in source.quarters.iter().enumerate() {
                    assert_eq!(quarter.id, id);
                    assert_eq!(quarter.key, format!("cylinder/quarter-{id}"));
                    assert_eq!(quarter.axis_parameters_mm, parameters);
                    let Surface::BSpline(SplineSupport::RationalPatch(patch)) = &quarter.surface
                    else {
                        unreachable!()
                    };
                    assert_eq!(
                        patch.carrier().domain(),
                        [[&q(id as i64), &q(id as i64 + 1)], [&q(0), &q(1)]]
                    );
                    assert!(patch.carrier().has_rational_definition());
                }
                eprintln!(
                    "canonical cylinder offset={offset} skew={skew}: {} work",
                    source.source_work + source.chart_work
                );
            });
            assert!(!work.exhausted, "{work:?}");
        }
    }
    #[test]
    fn existing_outer_wall_uses_its_original_origin_radius_and_axis_span() {
        let outer = Cylinder {
            origin: Point3::new(0.into(), 0.into(), 0.into()),
            axis: Vector3::new(0.into(), 1.into(), 0.into()),
            radius_sq: 144.into(),
            ref_dir: Vector3::new(1.into(), 0.into(), 0.into()),
        };
        let mut dec = Decider::default();
        let descriptor =
            NumericalCylinderSource::new(&outer, [q(0), q(12)], options(), &mut dec).unwrap();
        for (v, y) in [(q(0), q(0)), (q(1), q(12))] {
            let point = descriptor
                .source
                .evaluate_box([&q(0), &q(0)], [&v, &v], 128, &dec)
                .unwrap();
            for (actual, expected) in point.iter().zip([q(12), y, q(0)]) {
                assert!(actual.lo() <= expected && actual.hi() >= expected)
            }
        }
        assert_eq!(descriptor.quarters.len(), 4);
        assert_eq!(
            descriptor
                .signed_mm_at_v(&Rational::from_ratio(1, 2).unwrap())
                .unwrap(),
            q(6)
        );
    }

    #[test]
    fn standalone_descriptor_is_identical_at_one_two_and_four_workers() {
        let jobs = [
            ("10", false, "49"),
            ("10.5", false, "99/2"),
            ("10.5", true, "99/2"),
        ];
        let build = |(offset, skew, extent): (&str, bool, &str)| {
            let extent = Rational::parse(extent).unwrap();
            NumericalCylinderSource::new(
                &cylinder(offset, skew),
                [-&extent, extent],
                options(),
                &mut Decider::default(),
            )
            .unwrap()
        };
        let baseline: Vec<_> = jobs.iter().copied().map(build).collect();
        for workers in [1, 2, 4] {
            let rows = std::thread::scope(|scope| {
                let handles: Vec<_> = (0..workers)
                    .map(|rank| {
                        let jobs = &jobs;
                        let build = &build;
                        scope.spawn(move || {
                            jobs.iter()
                                .copied()
                                .enumerate()
                                .filter(|(i, _)| i % workers == rank)
                                .map(|(i, job)| (i, build(job)))
                                .collect::<Vec<_>>()
                        })
                    })
                    .collect();
                handles
                    .into_iter()
                    .flat_map(|handle| handle.join().unwrap())
                    .collect::<Vec<_>>()
            });
            for (i, actual) in rows {
                assert!(actual.source.same_definition(&baseline[i].source));
                assert_eq!(actual.evidence, baseline[i].evidence);
                assert_eq!(
                    (actual.source_work, actual.chart_work),
                    (baseline[i].source_work, baseline[i].chart_work)
                );
                for (a, b) in actual.quarters.iter().zip(&baseline[i].quarters) {
                    assert_eq!(a.key, b.key);
                    assert_eq!(a.axis_parameters_mm, b.axis_parameters_mm);
                    let Surface::BSpline(SplineSupport::RationalPatch(a)) = &a.surface else {
                        unreachable!()
                    };
                    let Surface::BSpline(SplineSupport::RationalPatch(b)) = &b.surface else {
                        unreachable!()
                    };
                    assert!(a.carrier().same_definition(b.carrier()));
                }
            }
        }
    }
    #[test]
    fn canonical_source_honors_an_expired_caller_deadline() {
        let (result, work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            dec.set_deadline(Some(std::time::Instant::now()));
            NumericalCylinderSource::new(
                &cylinder("10", false),
                [q(-49), q(49)],
                options(),
                &mut dec,
            )
        });
        assert_eq!(result.unwrap_err().code, "numbers.deadline");
        assert!(!work.exhausted);
    }
    #[test]
    fn canonical_source_rejects_invalid_domains_and_honors_the_caller_limits() {
        let cylinder = cylinder("10", false);
        for parameters in [[q(1), q(1)], [q(2), q(1)]] {
            assert!(
                NumericalCylinderSource::new(
                    &cylinder,
                    parameters,
                    options(),
                    &mut Decider::default()
                )
                .is_err()
            );
        }
        let (result, work) = with_operation_budget(0, || {
            NumericalCylinderSource::new(
                &cylinder,
                [q(-49), q(49)],
                options(),
                &mut Decider::default(),
            )
        });
        assert!(result.is_err());
        assert!(work.exhausted);
        let mut options = options();
        options.tolerance_mm = q(0);
        assert!(
            NumericalCylinderSource::new(
                &cylinder,
                [q(-49), q(49)],
                options,
                &mut Decider::default()
            )
            .is_err()
        );
    }
}
