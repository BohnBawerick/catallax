//! Tolerance rules T0 to T8 of roadmap item 3.0A (decision Q7, section 4.3
//! of the Checkpoint 2 decisions).
//!
//! The rules, as the decisions state them:
//!
//! | Rule | Content |
//! |---|---|
//! | T0 | Exact entities store no tolerance and have no size limit |
//! | T1 | Floor for f64 geometry: eps(R) = max(1e-7 mm, 1e-11 x R), R the coordinate size in the frame the operation runs in. Warn when one body exceeds 1 km. Operations run in local frames reached by exact rational transforms |
//! | T2 | Own fitted geometry: fit to max(1e-6 mm, 10 x eps), tightened to L/100 near small features; cap max(1e-5 mm, 100 x eps) |
//! | T3 | Imports: measure every entity; the file's declared value is only a hint. Automatic cap 0.01 mm for edges and faces, 0.02 mm for vertices. Above that, refit first, then the user may opt in up to 0.1 mm |
//! | T4 | Margin: 1.1 x (measured gap + measurement error), applied only to fresh measurements, never stacked |
//! | T5 | Pair rule: for two entities that share no topology, error if t_a + t_b > d/2, warn if t_a + t_b > d/5. Found with a bounding-box tree inflated by 5t |
//! | T6 | Short edges and slivers are invalid only through T5 |
//! | T7 | After each operation: refit the 2D copy, measure again, lower where possible. Warn if a tolerance rises more than 10x above its inputs |
//! | T8 | Export with a receiver profile: refit edges above the receiver tolerance, warn on features under 10x and 4x the receiver tolerance, write the largest entity tolerance as the file uncertainty, offer scaled export for micro parts |
//!
//! "The cap that applies is always the smaller of the class cap and the
//! local term." The class cap is the body's limit,
//! [`NumericalOptions::tolerance_mm`](super::NumericalOptions); the local
//! term is rule T5, which refuses a pair whose tolerances add up to more than
//! half their distance. This module holds the numbers and the checks that
//! read them. The milestone note `docs/milestones/3.0A-tolerance-rules.md`
//! says where each rule acts.

use catallax_numbers::{Decider, NumericalCheck, NumericalEvidence, Rational, Refusal};

use super::{NumericalDefinition, gap};

mod pairs;
pub use pairs::fitted_target_for_definition;
pub(super) use pairs::pair_rule_check;
pub use pairs::{FeaturePair, features_within};

fn mm(numerator: i64, denominator: i64) -> Rational {
    Rational::from_ratio(numerator, denominator).expect("nonzero")
}

// ---------------------------------------------------------------- T0 and T1

/// Rule T1's fixed part: 1e-7 mm.
/// Source: rule T1 of the Checkpoint 2 decisions, section 4.3.
pub fn floor_minimum_mm() -> Rational {
    mm(1, 10_000_000)
}

/// Rule T1's share of the coordinate size: 1e-11.
/// Source: rule T1 of the Checkpoint 2 decisions, section 4.3 (Parasolid's
/// resolution-to-size-box ratio, `tolerance-scenarios.md`, section 1).
pub fn floor_ratio() -> Rational {
    mm(1, 100_000_000_000)
}

/// Rule T1: `eps(R) = max(1e-7 mm, 1e-11 x R)`, the floor for binary64
/// geometry whose coordinates are `R` millimetres in size in the frame its
/// operation runs in. Rule T0: an exact entity stores no tolerance, so this
/// floor and the size warning below apply to binary64 entities only.
pub fn floor_mm(coordinate_size_mm: &Rational) -> Rational {
    (floor_ratio() * coordinate_size_mm.abs()).max(floor_minimum_mm())
}

/// Rule T1's size warning: one body larger than 1 km.
/// Source: rule T1 of the Checkpoint 2 decisions, section 4.3.
pub fn body_size_warning_mm() -> Rational {
    Rational::from_i64(1_000_000)
}

/// The frame rule T1 measures a body in: centred on the middle of the box of
/// its vertices and edge controls, reached by an exact rational shift.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BodyFrame {
    /// The exact centre of the box, in millimetres.
    pub centre: [Rational; 3],
    /// `R`: the largest coordinate size in that frame, in millimetres.
    pub coordinate_size_mm: Rational,
    /// The largest side of the box, in millimetres. Rule T1 warns above
    /// [`body_size_warning_mm`].
    pub extent_mm: Rational,
}

impl BodyFrame {
    /// The frame of a definition: the box of every vertex and of the points
    /// that define each edge's curve (a line's ends, a circle's centre with
    /// its radius, a spline's controls). A face lies between its edges or
    /// inside its surface's controls, so this box holds the body's
    /// coordinates.
    pub fn of(definition: &NumericalDefinition, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut low: Option<[Rational; 3]> = None;
        let mut high: Option<[Rational; 3]> = None;
        let mut take = |point: [Rational; 3], reach: &Rational| {
            let lo = point.clone().map(|c| c - reach.clone());
            let hi = point.map(|c| c + reach.clone());
            low = Some(match low.take() {
                None => lo,
                Some(old) => std::array::from_fn(|i| old[i].clone().min(lo[i].clone())),
            });
            high = Some(match high.take() {
                None => hi,
                Some(old) => std::array::from_fn(|i| old[i].clone().max(hi[i].clone())),
            });
        };
        let zero = Rational::zero();
        for vertex in &definition.vertices {
            dec.charge_work(1)?;
            take(vertex.position.clone(), &zero);
        }
        for edge in &definition.edges {
            for (point, reach) in
                gap::defining_points(&edge.curve, edge.parameters.each_ref(), dec)?
            {
                dec.charge_work(1)?;
                take(point, &reach);
            }
        }
        let (Some(low), Some(high)) = (low, high) else {
            return Ok(Self {
                centre: [zero.clone(), zero.clone(), zero.clone()],
                coordinate_size_mm: zero.clone(),
                extent_mm: zero,
            });
        };
        let half = mm(1, 2);
        let centre: [Rational; 3] = std::array::from_fn(|i| (&low[i] + &high[i]) * half.clone());
        let extent_mm = (0..3)
            .map(|i| &high[i] - &low[i])
            .fold(Rational::zero(), Rational::max);
        Ok(Self {
            centre,
            coordinate_size_mm: extent_mm.clone() * half,
            extent_mm,
        })
    }

    /// The frame of a set of exact points, for a producer that has no
    /// numerical definition yet: an exact model's vertices.
    pub fn of_points(points: &[[Rational; 3]]) -> Self {
        let zero = Rational::zero();
        let Some(first) = points.first() else {
            return Self {
                centre: [zero.clone(), zero.clone(), zero.clone()],
                coordinate_size_mm: zero.clone(),
                extent_mm: zero,
            };
        };
        let mut low = first.clone();
        let mut high = first.clone();
        for point in &points[1..] {
            for axis in 0..3 {
                low[axis] = low[axis].clone().min(point[axis].clone());
                high[axis] = high[axis].clone().max(point[axis].clone());
            }
        }
        let half = mm(1, 2);
        let extent_mm = (0..3)
            .map(|i| &high[i] - &low[i])
            .fold(Rational::zero(), Rational::max);
        Self {
            centre: std::array::from_fn(|i| (&low[i] + &high[i]) * half.clone()),
            coordinate_size_mm: extent_mm.clone() * half,
            extent_mm,
        }
    }

    /// Rule T1's floor in this frame.
    pub fn floor_mm(&self) -> Rational {
        floor_mm(&self.coordinate_size_mm)
    }
}

// ---------------------------------------------------------------------- T2

/// The fit target and cap of rule T2 for own fitted geometry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FittedTolerance {
    /// What a fitter aims at, in millimetres.
    pub target_mm: Rational,
    /// The largest tolerance an entity may store, in millimetres: the class
    /// cap a producer states as its body's limit.
    pub cap_mm: Rational,
}

/// Rule T2's fixed part of the fit target: 1e-6 mm.
/// Source: rule T2 of the Checkpoint 2 decisions, section 4.3.
pub fn fitted_target_minimum_mm() -> Rational {
    mm(1, 1_000_000)
}

/// Rule T2's fixed part of the cap: 1e-5 mm.
/// Source: rule T2 of the Checkpoint 2 decisions, section 4.3.
pub fn fitted_cap_minimum_mm() -> Rational {
    mm(1, 100_000)
}

/// Rule T2's multiple of the floor in the fit target, 10, and in the cap,
/// 100. Source: rule T2 of the Checkpoint 2 decisions, section 4.3.
const TARGET_FLOORS: i64 = 10;
/// See [`TARGET_FLOORS`].
const CAP_FLOORS: i64 = 100;
/// Rule T2's share of a small feature's size: the target is tightened to
/// `L/100`. Source: rule T2 of the Checkpoint 2 decisions, section 4.3.
const FEATURE_SHARE: i64 = 100;

/// Rule T2 for a floor `eps` (rule T1) and, where the producer knows one, the
/// size `L` of the smallest feature near what it fits.
///
/// The target is `max(1e-6 mm, 10 x eps)`, tightened to `L/100`. The
/// tightening never goes under `eps`, which rule T1 makes the floor of every
/// binary64 entity. The cap is `max(1e-5 mm, 100 x eps)`; a small feature
/// does not lower it, because rule T5 is the local term the cap is compared
/// with.
pub fn own_fitted(eps: &Rational, small_feature_mm: Option<&Rational>) -> FittedTolerance {
    let mut target_mm = fitted_target_minimum_mm().max(eps * &Rational::from_i64(TARGET_FLOORS));
    if let Some(size) = small_feature_mm {
        let tightened = size / &Rational::from_i64(FEATURE_SHARE);
        target_mm = target_mm.min(tightened).max(eps.clone());
    }
    FittedTolerance {
        target_mm,
        cap_mm: fitted_cap_minimum_mm().max(eps * &Rational::from_i64(CAP_FLOORS)),
    }
}

/// The distance under which a feature tightens rule T2's target: a feature
/// of size `L` tightens it only when `L/100` is under the untightened target.
pub fn fitted_feature_reach_mm(eps: &Rational) -> Rational {
    own_fitted(eps, None).target_mm * Rational::from_i64(FEATURE_SHARE)
}

// ---------------------------------------------------------------------- T3

/// The kinds of entity a tolerance belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, catallax_numbers::Transfer)]
pub enum EntityKind {
    /// A vertex.
    Vertex,
    /// An edge.
    Edge,
    /// A face.
    Face,
}

impl EntityKind {
    /// A stable lower-case name.
    pub fn name(self) -> &'static str {
        match self {
            Self::Vertex => "vertex",
            Self::Edge => "edge",
            Self::Face => "face",
        }
    }
}

/// Rule T3's automatic cap for imported edges and faces: 0.01 mm.
/// Source: rule T3 of the Checkpoint 2 decisions, section 4.3.
pub fn import_cap_mm() -> Rational {
    mm(1, 100)
}

/// Rule T3's automatic cap for imported vertices: 0.02 mm.
/// Source: rule T3 of the Checkpoint 2 decisions, section 4.3.
pub fn import_vertex_cap_mm() -> Rational {
    mm(2, 100)
}

/// Rule T3's largest value a user may opt in to: 0.1 mm.
/// Source: rule T3 of the Checkpoint 2 decisions, section 4.3.
pub fn import_opt_in_limit_mm() -> Rational {
    mm(1, 10)
}

/// Rule T3's cap for one imported entity. Without an opt-in it is the
/// automatic cap of the entity's kind. An opt-in raises the cap of every
/// kind to the stated value, which may not exceed 0.1 mm; an opt-in under an
/// automatic cap does not lower it.
pub fn import_cap_for(kind: EntityKind, opt_in_mm: Option<&Rational>) -> Result<Rational, Refusal> {
    let automatic = match kind {
        EntityKind::Vertex => import_vertex_cap_mm(),
        EntityKind::Edge | EntityKind::Face => import_cap_mm(),
    };
    match opt_in_mm {
        None => Ok(automatic),
        Some(value) => {
            if *value <= Rational::zero() || *value > import_opt_in_limit_mm() {
                return Err(Refusal::new(
                    "numerical.import-opt-in",
                    format!(
                        "an import tolerance opt-in must be positive and at most {} mm (rule T3); {value} mm was asked",
                        import_opt_in_limit_mm()
                    ),
                ));
            }
            Ok(automatic.max(value.clone()))
        }
    }
}

/// Admission options for an import under rule T3 with a stated cap: the
/// automatic 0.01 mm, an opt-in up to 0.1 mm, or a stricter caller's value.
/// Edges and faces take the cap; a vertex takes 0.02 mm at the automatic
/// cap, the opt-in where that is larger, and a stricter caller's cap below
/// the automatic one.
pub fn import_options(cap_mm: Rational) -> super::NumericalOptions {
    let automatic = import_cap_mm();
    let vertex_tolerance_mm = if cap_mm < automatic {
        None
    } else {
        let opt_in = (cap_mm > automatic).then_some(&cap_mm);
        import_cap_for(EntityKind::Vertex, opt_in).ok()
    };
    super::NumericalOptions {
        tolerance_mm: cap_mm,
        vertex_tolerance_mm,
    }
}

/// Rule T3's caps by class on admission's records: the stored allowance
/// for an edge end at its vertex against `vertex_mm`, every other sampled
/// record against `cap_mm`. Refuses as the evidence does, with
/// `numbers.numerical-tolerance-exceeded`.
pub(super) fn check_by_class(
    checks: &[NumericalCheck],
    cap_mm: &Rational,
    vertex_mm: &Rational,
) -> Result<(), Refusal> {
    for check in checks {
        if let NumericalCheck::Samples {
            name,
            scope,
            samples,
        } = check
        {
            for sample in samples {
                let (stored, limit) = if is_vertex_record(name) {
                    (vertex_allowance(&sample.distance_upper_mm), vertex_mm)
                } else {
                    (sample.distance_upper_mm.clone(), cap_mm)
                };
                if stored > *limit {
                    return Err(Refusal::new(
                        "numbers.numerical-tolerance-exceeded",
                        format!("{name} on {scope}: tolerance {stored} mm exceeds {limit} mm",),
                    ));
                }
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------- T5

/// Rule T5's error share: error when `t_a + t_b > d/2`.
/// Source: rule T5 of the Checkpoint 2 decisions, section 4.3.
pub const PAIR_ERROR_DIVISOR: f64 = 2.0;
/// Rule T5's warning share: warn when `t_a + t_b > d/5`.
/// Source: rule T5 of the Checkpoint 2 decisions, section 4.3.
pub const PAIR_WARNING_DIVISOR: f64 = 5.0;
/// Rule T5's box inflation: each box grows by 5 times its entity's
/// tolerance, so two boxes overlap whenever their entities could be closer
/// than the warning distance `5 (t_a + t_b)`.
/// Source: rule T5 of the Checkpoint 2 decisions, section 4.3.
pub const PAIR_INFLATION: f64 = 5.0;

/// What rule T5 says of one pair.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PairVerdict {
    /// `t_a + t_b <= d/5`.
    Clear,
    /// `d/5 < t_a + t_b <= d/2`: a warning.
    Warn,
    /// `t_a + t_b > d/2`: an error.
    Error,
}

/// Rule T5 for two entities with tolerances `t_a` and `t_b` at distance `d`.
pub fn pair_verdict(t_a: f64, t_b: f64, distance_mm: f64) -> PairVerdict {
    let sum = t_a + t_b;
    if sum > distance_mm / PAIR_ERROR_DIVISOR {
        PairVerdict::Error
    } else if sum > distance_mm / PAIR_WARNING_DIVISOR {
        PairVerdict::Warn
    } else {
        PairVerdict::Clear
    }
}

// ---------------------------------------------------------------------- T7

/// Rule T7's growth factor: warn when a tolerance rises more than 10 times
/// above its inputs. Source: rule T7 of the Checkpoint 2 decisions, section
/// 4.3. The decisions call this number their weakest evidence: no source
/// gives one.
pub const GROWTH_FACTOR: i64 = 10;

/// Rule T7's trigger: `output > 10 x input`. A warning, never a refusal.
pub fn grew(input_mm: &Rational, output_mm: &Rational) -> bool {
    *output_mm > input_mm * &Rational::from_i64(GROWTH_FACTOR)
}

// ---------------------------------------------------------------------- T8

/// A receiving system's modelling tolerance (rule T8).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReceiverProfile {
    /// The receiver's name, for reports.
    pub name: String,
    /// Its modelling tolerance, in millimetres.
    pub tolerance_mm: Rational,
}

impl ReceiverProfile {
    /// Siemens NX at 0.01 mm.
    /// Source: rule T8 of the Checkpoint 2 decisions, section 4.3 (NX
    /// modelling distance tolerance, `tolerance-numbers-research.md`,
    /// section 3).
    pub fn nx() -> Self {
        Self {
            name: "NX".into(),
            tolerance_mm: mm(1, 100),
        }
    }

    /// CATIA V5 at 0.001 mm.
    /// Source: rule T8 of the Checkpoint 2 decisions, section 4.3 (CATIA
    /// join and healing merging distance, `tolerance-numbers-research.md`,
    /// section 3).
    pub fn catia() -> Self {
        Self {
            name: "CATIA".into(),
            tolerance_mm: mm(1, 1000),
        }
    }

    /// Another receiver, by name and positive tolerance.
    pub fn custom(name: impl Into<String>, tolerance_mm: Rational) -> Result<Self, Refusal> {
        let name = name.into();
        if name.trim().is_empty() || tolerance_mm <= Rational::zero() {
            return Err(Refusal::new(
                "numerical.receiver-profile",
                "a receiver profile needs a name and a positive tolerance",
            ));
        }
        Ok(Self { name, tolerance_mm })
    }

    /// A named profile: `nx` or `catia`, in any case.
    pub fn named(name: &str) -> Option<Self> {
        match name.to_ascii_lowercase().as_str() {
            "nx" => Some(Self::nx()),
            "catia" => Some(Self::catia()),
            _ => None,
        }
    }
}

/// Rule T8's two feature bands: a feature under 10 times, and one under 4
/// times, the receiver's tolerance. Source: rule T8 of the Checkpoint 2
/// decisions, section 4.3.
pub const RECEIVER_WIDE_BAND: i64 = 10;
/// See [`RECEIVER_WIDE_BAND`].
pub const RECEIVER_NARROW_BAND: i64 = 4;

/// Which of rule T8's bands a feature falls in.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum FeatureBand {
    /// Under 4 times the receiver's tolerance.
    UnderFour,
    /// Under 10 times, and at least 4 times.
    UnderTen,
}

/// Rule T8's feature warning for a feature of size `size_mm`, or none.
pub fn feature_band(size_mm: &Rational, receiver: &ReceiverProfile) -> Option<FeatureBand> {
    let narrow = &receiver.tolerance_mm * &Rational::from_i64(RECEIVER_NARROW_BAND);
    let wide = &receiver.tolerance_mm * &Rational::from_i64(RECEIVER_WIDE_BAND);
    if *size_mm < narrow {
        Some(FeatureBand::UnderFour)
    } else if *size_mm < wide {
        Some(FeatureBand::UnderTen)
    } else {
        None
    }
}

// ----------------------------------------------------------- entity report

/// One entity's tolerance, as admission stored it.
#[derive(Clone, Debug, PartialEq, Eq, catallax_numbers::Transfer)]
pub struct EntityTolerance {
    /// Vertex or edge. A face of a numerical body stores no tolerance: its
    /// surface is the authoritative copy, and the gap lives on its edges.
    pub kind: EntityKind,
    /// The stable key.
    pub key: String,
    /// The tolerance, in millimetres. Never under rule T1's floor.
    pub tolerance_mm: Rational,
}

/// A warning of the tolerance rules. Warnings never refuse.
#[derive(Clone, Debug, PartialEq, Eq, catallax_numbers::Transfer)]
pub enum ToleranceWarning {
    /// Rule T1: the body is larger than 1 km.
    LargeBody {
        /// The largest side of the body's box, in millimetres.
        extent_mm: Rational,
    },
    /// Rule T5: two entities that share no topology, with
    /// `d/5 < t_a + t_b <= d/2`.
    Pair(FeaturePair),
    /// Rule T7: an entity's tolerance rose more than 10 times above its
    /// input. The absolute values come with it: the ratio alone says nothing
    /// of risk when the input is near zero.
    Growth {
        /// Vertex or edge.
        kind: EntityKind,
        /// The entity's key.
        key: String,
        /// The input tolerance it is compared with, in millimetres: the same
        /// entity's tolerance before the operation, or the largest
        /// tolerance of the operation's input where the entity is new.
        input_mm: Rational,
        /// Whether `input_mm` is the entity's own earlier tolerance.
        same_entity: bool,
        /// The tolerance after the operation, in millimetres.
        output_mm: Rational,
        /// The body's class cap, in millimetres.
        cap_mm: Rational,
        /// `cap_mm - output_mm`: what is left under the cap.
        headroom_mm: Rational,
    },
}

impl ToleranceWarning {
    /// A stable code.
    pub fn code(&self) -> &'static str {
        match self {
            Self::LargeBody { .. } => "tolerance.t1-body-over-1-km",
            Self::Pair(_) => "tolerance.t5-pair-close",
            Self::Growth { .. } => "tolerance.t7-growth-over-10x",
        }
    }
}

/// What the tolerance rules found on one admitted body.
#[derive(Clone, Debug, Default, PartialEq, Eq, catallax_numbers::Transfer)]
pub struct ToleranceReport {
    /// Rule T1's coordinate size `R` in the body's frame, in millimetres.
    pub coordinate_size_mm: Rational,
    /// Rule T1's floor in that frame.
    pub floor_mm: Rational,
    /// The body's class cap.
    pub cap_mm: Rational,
    /// The vertex cap where the class has its own (rule T3).
    pub vertex_cap_mm: Option<Rational>,
    /// Every vertex and edge with its tolerance, in definition order.
    pub entities: Vec<EntityTolerance>,
    /// Warnings of rules T1, T5 and T7.
    pub warnings: Vec<ToleranceWarning>,
}

impl ToleranceReport {
    /// The largest entity tolerance, or rule T1's floor for a body with none.
    pub fn largest_mm(&self) -> Rational {
        self.entities
            .iter()
            .map(|entity| entity.tolerance_mm.clone())
            .fold(self.floor_mm.clone(), Rational::max)
    }

    /// The tolerance of one entity.
    pub fn of(&self, kind: EntityKind, key: &str) -> Option<&Rational> {
        self.entities
            .iter()
            .find(|entity| entity.kind == kind && entity.key == key)
            .map(|entity| &entity.tolerance_mm)
    }

    /// The headroom under the cap of the largest entity tolerance.
    pub fn headroom_mm(&self) -> Rational {
        &self.cap_mm - &self.largest_mm()
    }
}

/// Source: rule T4: a fresh measurement is stored as 1.1 times its value.
fn margin() -> Rational {
    mm(11, 10)
}

fn is_edge_surface_record(name: &str) -> bool {
    matches!(
        name,
        "trim.edge-surface-distance"
            | "trim.numerical-pair-edge-surface-distance"
            | "trim.edge-surface-derived-bound"
            | "trim.numerical-pair-edge-surface-derived-bound"
    )
}

fn is_vertex_record(name: &str) -> bool {
    matches!(
        name,
        "edge.vertex-distance" | "edge.numerical-pair-vertex-distance"
    )
}

fn vertex_allowance(distance_upper_mm: &Rational) -> Rational {
    &margin() * distance_upper_mm
}

/// The stored record of every edge use, in definition order (face, loop,
/// use), and whether it is a derived copy's bound.
pub(super) fn use_records(evidence: &NumericalEvidence) -> Vec<(Rational, bool)> {
    let mut out = Vec::new();
    for check in evidence.checks() {
        if let NumericalCheck::Samples { name, samples, .. } = check
            && is_edge_surface_record(name)
        {
            let stored = samples
                .iter()
                .map(|sample| sample.distance_upper_mm.clone())
                .fold(Rational::zero(), Rational::max);
            out.push((stored, name.ends_with("derived-bound")));
        }
    }
    out
}

/// Every vertex and edge tolerance of an admitted body, read from its
/// evidence in the order admission wrote it.
///
/// - An edge's tolerance is the largest stored record of its uses: rule T4's
///   1.1 times the measured gap, or a derived copy's bound.
/// - A vertex's tolerance is at least each of its edges' (vertex >= edge, the
///   order of Jackson 1995 that section 4.5 of the decisions keeps), and at
///   least rule T4's margin on the measured distance from each edge end to
///   the vertex.
/// - Rule T1's floor is the least value of either.
pub(super) fn entity_tolerances(
    definition: &NumericalDefinition,
    evidence: &NumericalEvidence,
    floor: &Rational,
    vertex_cap: Option<&Rational>,
) -> Result<Vec<EntityTolerance>, Refusal> {
    let mut edges = vec![floor.clone(); definition.edges.len()];
    let mut vertices = vec![floor.clone(); definition.vertices.len()];
    let mut vertex_records = Vec::new();
    let mut use_records = Vec::new();
    for check in evidence.checks() {
        if let NumericalCheck::Samples { name, samples, .. } = check {
            if is_vertex_record(name) {
                vertex_records.push(samples);
            } else if is_edge_surface_record(name) {
                use_records.push(samples);
            }
        }
    }
    let mismatch = || {
        Refusal::new(
            "numerical.invalid-definition",
            "the evidence does not match the definition's edges and uses",
        )
    };
    if vertex_records.len() != definition.edges.len() {
        return Err(mismatch());
    }
    let uses: Vec<usize> = definition
        .faces
        .iter()
        .flat_map(|face| face.loops.iter().flatten().map(|trim| trim.edge))
        .collect();
    if uses.len() != use_records.len() {
        return Err(mismatch());
    }
    for (edge, samples) in uses.iter().zip(&use_records) {
        for sample in samples.iter() {
            edges[*edge] = edges[*edge].clone().max(sample.distance_upper_mm.clone());
        }
    }
    for (index, edge) in definition.edges.iter().enumerate() {
        for (end, sample) in edge.vertices.iter().zip(vertex_records[index].iter()) {
            let measured = vertex_allowance(&sample.distance_upper_mm);
            vertices[*end] = vertices[*end]
                .clone()
                .max(measured)
                .max(edges[index].clone());
        }
    }
    let mut entities = Vec::with_capacity(vertices.len() + edges.len());
    for (vertex, tolerance_mm) in definition.vertices.iter().zip(vertices) {
        if let Some(cap) = vertex_cap
            && tolerance_mm > *cap
        {
            return Err(Refusal::new(
                "numbers.numerical-tolerance-exceeded",
                format!(
                    "vertex {}: tolerance {tolerance_mm} mm exceeds {cap} mm",
                    vertex.key
                ),
            ));
        }
        entities.push(EntityTolerance {
            kind: EntityKind::Vertex,
            key: vertex.key.clone(),
            tolerance_mm,
        });
    }
    for (edge, tolerance_mm) in definition.edges.iter().zip(edges) {
        entities.push(EntityTolerance {
            kind: EntityKind::Edge,
            key: edge.key.clone(),
            tolerance_mm,
        });
    }
    Ok(entities)
}

/// The largest entity tolerance a check records: the largest stored record
/// of an edge use, and rule T4's margin on the largest distance from an edge
/// end to its vertex. The same reading as [`ToleranceReport::largest_mm`],
/// for a check with no body of its own, such as an export's check of its
/// written values.
pub fn largest_recorded(evidence: &NumericalEvidence) -> Rational {
    let mut largest = Rational::zero();
    for check in evidence.checks() {
        if let NumericalCheck::Samples { name, samples, .. } = check {
            for sample in samples {
                let value = if is_vertex_record(name) {
                    vertex_allowance(&sample.distance_upper_mm)
                } else if is_edge_surface_record(name) {
                    sample.distance_upper_mm.clone()
                } else {
                    continue;
                };
                largest = largest.max(value);
            }
        }
    }
    largest
}

/// A report's entity tolerances by kind and key, charging one unit each.
pub(super) fn by_key<'a>(
    report: &'a ToleranceReport,
    dec: &mut Decider,
) -> Result<std::collections::BTreeMap<(EntityKind, &'a str), &'a Rational>, Refusal> {
    dec.charge_work(report.entities.len() as u64)?;
    Ok(report
        .entities
        .iter()
        .map(|entity| ((entity.kind, entity.key.as_str()), &entity.tolerance_mm))
        .collect())
}

/// Rule T7's growth warnings of one operation: every entity of `output`
/// whose tolerance is more than 10 times its input. The input of an entity
/// that `input` also has, by key and kind, is its own earlier tolerance; the
/// input of a new entity is the largest tolerance of `input`. Each entity
/// compared is charged.
///
/// An input with no recorded entities (a body read back from a layout that
/// had no report) gives no warning: there is nothing to compare with.
pub fn growth_warnings(
    input: &ToleranceReport,
    output: &ToleranceReport,
    dec: &mut Decider,
) -> Result<Vec<ToleranceWarning>, Refusal> {
    if input.entities.is_empty() {
        return Ok(Vec::new());
    }
    let largest = input.largest_mm();
    let earlier = by_key(input, dec)?;
    let mut warnings = Vec::new();
    for entity in &output.entities {
        dec.charge_work(1)?;
        let (input_mm, same_entity) = match earlier.get(&(entity.kind, entity.key.as_str())) {
            Some(value) => ((*value).clone(), true),
            None => (largest.clone(), false),
        };
        if grew(&input_mm, &entity.tolerance_mm) {
            warnings.push(ToleranceWarning::Growth {
                kind: entity.kind,
                key: entity.key.clone(),
                input_mm,
                same_entity,
                output_mm: entity.tolerance_mm.clone(),
                cap_mm: output.cap_mm.clone(),
                headroom_mm: &output.cap_mm - &entity.tolerance_mm,
            });
        }
    }
    Ok(warnings)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(text: &str) -> Rational {
        Rational::parse(text).expect("number")
    }

    #[test]
    fn t1_floor_is_the_larger_of_its_two_terms() {
        assert_eq!(floor_mm(&q("0")), q("0.0000001"));
        assert_eq!(floor_mm(&q("10000")), q("0.0000001"));
        assert_eq!(floor_mm(&q("30000")), q("0.0000003"));
        // About 1e-5 mm at 1 km, as the decisions say.
        assert_eq!(floor_mm(&q("1000000")), q("0.00001"));
        assert_eq!(body_size_warning_mm(), q("1000000"));
    }

    #[test]
    fn t2_target_and_cap_follow_the_floor() {
        let at_ten_metres = own_fitted(&floor_mm(&q("10000")), None);
        assert_eq!(at_ten_metres.target_mm, q("0.000001"));
        assert_eq!(at_ten_metres.cap_mm, q("0.00001"));
        let at_a_kilometre = own_fitted(&floor_mm(&q("1000000")), None);
        assert_eq!(at_a_kilometre.target_mm, q("0.0001"));
        assert_eq!(at_a_kilometre.cap_mm, q("0.001"));
        // A feature of 1e-5 mm tightens the target to its hundredth, but not
        // under the floor; the cap stays.
        let eps = floor_mm(&q("100"));
        let near = own_fitted(&eps, Some(&q("0.00005")));
        assert_eq!(near.target_mm, q("0.0000005"));
        assert_eq!(near.cap_mm, q("0.00001"));
        assert_eq!(own_fitted(&eps, Some(&q("0.000001"))).target_mm, eps);
        // A feature larger than 100 times the target does not tighten it.
        assert_eq!(fitted_feature_reach_mm(&eps), q("0.0001"));
        assert_eq!(
            own_fitted(&eps, Some(&q("0.0001"))).target_mm,
            q("0.000001")
        );
    }

    #[test]
    fn t3_caps_and_opt_in() {
        assert_eq!(import_cap_for(EntityKind::Edge, None).unwrap(), q("0.01"));
        assert_eq!(import_cap_for(EntityKind::Face, None).unwrap(), q("0.01"));
        assert_eq!(import_cap_for(EntityKind::Vertex, None).unwrap(), q("0.02"));
        let opt = q("0.1");
        assert_eq!(import_cap_for(EntityKind::Edge, Some(&opt)).unwrap(), opt);
        assert_eq!(
            import_cap_for(EntityKind::Vertex, Some(&q("0.015"))).unwrap(),
            q("0.02")
        );
        assert_eq!(
            import_cap_for(EntityKind::Edge, Some(&q("0.11")))
                .unwrap_err()
                .code,
            "numerical.import-opt-in"
        );
    }

    #[test]
    fn t5_thresholds() {
        // t_a + t_b = 0.75: clear at d = 3.75, a warning just under it and
        // at d = 1.5, an error just under 1.5.
        assert_eq!(pair_verdict(0.25, 0.5, 3.75), PairVerdict::Clear);
        assert_eq!(pair_verdict(0.25, 0.5, 3.74), PairVerdict::Warn);
        assert_eq!(pair_verdict(0.25, 0.5, 1.5), PairVerdict::Warn);
        assert_eq!(pair_verdict(0.25, 0.5, 1.49), PairVerdict::Error);
        assert_eq!(PAIR_INFLATION, 5.0);
    }

    #[test]
    fn t7_trigger_is_more_than_ten_times() {
        assert!(!grew(&q("0.000001"), &q("0.00001")));
        assert!(grew(&q("0.000001"), &q("0.0000100001")));
    }

    #[test]
    fn t8_profiles_and_bands() {
        let nx = ReceiverProfile::nx();
        let catia = ReceiverProfile::catia();
        assert_eq!(nx.tolerance_mm, q("0.01"));
        assert_eq!(catia.tolerance_mm, q("0.001"));
        assert_eq!(ReceiverProfile::named("CATIA"), Some(catia.clone()));
        assert_eq!(feature_band(&q("0.039"), &nx), Some(FeatureBand::UnderFour));
        assert_eq!(feature_band(&q("0.04"), &nx), Some(FeatureBand::UnderTen));
        assert_eq!(feature_band(&q("0.1"), &nx), None);
        assert_eq!(feature_band(&q("0.1"), &catia), None);
        assert!(ReceiverProfile::custom("x", q("0")).is_err());
    }
}
