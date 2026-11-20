//! The gap rule for imports (decision Q9, and section 4.2 step 2 of the
//! Checkpoint 2 decisions).
//!
//! "Own builds: same-parameter gap. Imports: geometric gap first; if the
//! parametric gap is much larger, reparametrise instead of raising the
//! tolerance."
//!
//! Admission stores the same-parameter gap of every edge use, because the
//! kernel reads both copies of an edge at one parameter. An own build makes
//! its copies in one parameter, so that is the whole rule for it. An import
//! promises no shared parameter: its copies can lie on each other and still
//! run at different speeds. Storing the same-parameter gap would then raise
//! a tolerance to cover a difference of parameter, not of geometry.
//!
//! So an imported use whose same-parameter gap could be lowered is read
//! geometrically, both ways ([`gap::measure_geometric`]). Where the
//! same-parameter gap is much larger than the geometric gap, the face copy
//! is given the parameter of the 3D curve ([`gap::refit::reparametrise`]),
//! and admission continues with the same-parameter gap of the new copy.
//!
//! "Much larger" is read with the one number rule T4 already has, its margin
//! of 1.1: a same-parameter gap within 1.1 times the geometric gap stores no
//! more than rule T4's tolerance of a gap that size already covers. Above
//! that, the stored tolerance would be covering the parameter. The same test
//! says when the change of parameter is done, and whether it was worth
//! keeping. Neither kernel this was checked against has a ratio of its own.
//! OpenCASCADE's `BRepLib::SameParameter` reparametrises an edge whose
//! samples are off at the same parameter by more than its working tolerance,
//! and keeps the new face curve whenever its deviation is not larger than
//! the old one's. SMLib ties no parameter of a u/v trim to its edge:
//! `SmCrvCrvGapFunction` drops each sample of the edge on the lifted trim,
//! so its edge gaps are geometric throughout.
//!
//! Work: section 4.2 step 8 gives an edge one hard work cap that refuses by
//! name. The rule reads a use up to three times, and all of it draws on that
//! use's one allowance ([`read_use`]). The cap's refusal,
//! `numerical.gap-work-cap`, is returned under its own name wherever it
//! comes from, and never answered with the body as imported
//! ([`settled_body`]).

use super::gap::{self, EdgeCopies, FaceCopyRef, GapSurface, Hints, Limits};
use super::*;

/// What the gap rule did with one imported edge use.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ImportGapOutcome {
    /// The same-parameter gap is within rule T4's margin of the geometric
    /// gap: the gap is in the geometry, and the copy is kept as imported.
    Geometric,
    /// The same-parameter gap was much larger than the geometric gap. The
    /// face copy was given the parameter of the 3D curve, and the body was
    /// admitted with it.
    Reparametrised,
    /// The same-parameter gap was much larger than the geometric gap, and no
    /// change of parameter lowered it by more than rule T4's margin. The copy
    /// is kept as imported and its tolerance covers the same-parameter gap.
    NotLowered,
    /// A change of parameter lowered the gap, and the body was not admitted
    /// with the new copy. The copy is kept as imported.
    NotAdmitted,
}

/// One imported edge use the gap rule read geometrically.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ImportGapRow {
    /// Key of the face of the use.
    pub face: String,
    /// Key of the edge.
    pub edge: String,
    /// Loop index on the face.
    pub boundary: usize,
    /// Use index in the loop.
    pub use_index: usize,
    /// Rule T4's tolerance of the same-parameter gap of the copy as
    /// imported, in millimetres.
    pub parametric_mm: Rational,
    /// The geometric gap as reported, in millimetres: the larger of the two
    /// directions plus its rounding error. It carries no margin, because it
    /// is never stored.
    pub geometric_mm: Rational,
    /// Rule T4's tolerance of the same-parameter gap of the reparametrised
    /// copy, where one was made, in millimetres.
    pub reparametrised_mm: Option<Rational>,
    /// Spans of the reparametrised copy; zero where none was made.
    pub spans: usize,
    /// What the rule did.
    pub outcome: ImportGapOutcome,
}

/// Whether `a` is much larger than `b`: more than rule T4's margin times it.
/// Source: rule T4 of the tolerance rules (Checkpoint 2 decisions, section
/// 4.3), read as the module's header says.
fn much_larger(a: f64, b: f64) -> bool {
    a > gap::MARGIN * b
}

fn exact(value: f64) -> Result<Rational, Refusal> {
    Rational::from_f64_exact(value)
        .ok_or_else(|| bad("trim", "an import gap is not a finite number"))
}

/// What the rule read of one edge use.
#[derive(Clone, Debug)]
pub(super) struct UseReading {
    /// The same-parameter gap of the copy as imported.
    pub parametric: gap::GapMeasure,
    /// Its geometric gap.
    pub geometric: gap::GeometricMeasure,
    /// The copy in the parameter of the 3D curve, where the same-parameter
    /// gap was much larger than the geometric gap and the old copy has a
    /// change of parameter.
    pub reparametrised: Option<gap::refit::Refit>,
    /// What the rule did: `Geometric`, `NotLowered` or `Reparametrised`.
    pub outcome: ImportGapOutcome,
}

/// The rule for one edge use, under the use's one work allowance.
///
/// Section 4.2 step 8: "A hard per-edge work cap refuses by name." The rule
/// does up to three things with a use, and `limits` is the allowance for all
/// of them together:
///
/// 1. The same-parameter gap of the copy as imported. Admission read it
///    before and kept the stored tolerance only, so it is read again here
///    and counted.
/// 2. The geometric gap, both ways, under what the first reading left.
/// 3. The change of parameter, under what the two readings left: its
///    projections, its readings of the old copy and the checker's work on
///    every candidate.
///
/// Where the allowance runs out at any of these, the use refuses with
/// `numerical.gap-work-cap`. Nothing read before that stands in for the
/// answer: not the copy as imported, and not an earlier candidate.
///
/// `None` for a use the rule leaves as it is: a straight residual, a copy
/// that leaves its surface's domain, a use that stores no more than
/// `floor_mm`, and a copy with no u/v of its own.
pub(super) fn read_use(
    copies: &EdgeCopies<'_>,
    limits: Limits,
    floor_mm: f64,
    dec: &mut Decider,
) -> Result<Option<UseReading>, Refusal> {
    let parametric = gap::measure(copies, &Hints::default(), limits, dec)?;
    if parametric.method == gap::Method::StraightResidual
        || parametric.domain_exit.is_some()
        || parametric.stored_mm() <= floor_mm
    {
        return Ok(None);
    }
    let left = Limits {
        work: limits.work.saturating_sub(parametric.work),
    };
    let Some(geometric) = gap::measure_geometric(copies, &Hints::default(), left, dec)? else {
        return Ok(None);
    };
    let left = Limits {
        work: left.work.saturating_sub(geometric.work),
    };
    let mut reparametrised = None;
    let mut outcome = ImportGapOutcome::Geometric;
    if much_larger(parametric.reported_mm(), geometric.reported_mm()) {
        // Done when the new copy's same-parameter gap is no longer much
        // larger than the geometric gap, or stores the floor.
        let target = (gap::MARGIN * gap::MARGIN * geometric.reported_mm()).max(floor_mm);
        outcome = ImportGapOutcome::NotLowered;
        reparametrised = gap::refit::reparametrise(copies, target, left, dec)?;
        if let Some(next) = &reparametrised
            && much_larger(parametric.stored_mm(), next.stored_mm)
        {
            outcome = ImportGapOutcome::Reparametrised;
        }
    }
    Ok(Some(UseReading {
        parametric,
        geometric,
        reparametrised,
        outcome,
    }))
}

/// Reads every imported edge use the gap rule applies to and gives the uses
/// whose copies run at different speeds a face copy in the parameter of
/// their 3D curve.
///
/// `known` holds the stored record of every use, in definition order, where
/// the definition was admitted already. A use whose stored tolerance is at
/// rule T1's least floor is not read: no tolerance is stored under the
/// floor, so nothing could be lowered. Without `known` every use is measured
/// here first.
///
/// These uses are left as they are:
///
/// - A derived copy, which admission does not measure.
/// - A curve projected on a plane. It has no u/v of its own, and its gap is
///   the distance of the 3D curve from the plane.
/// - A straight residual. Both copies are straight at constant speed, so
///   their largest distance is at an end, and the ends are paired with each
///   other in the geometric gap too: the two gaps are equal.
///
/// `each_use` is the one work allowance of an edge use, for all the rule
/// does with it ([`read_use`]): the checker's cap. A use that spends it
/// refuses the admission with `numerical.gap-work-cap`.
///
/// Returns the retained u/v curves with the new copies in place, or `None`
/// where no copy changed, and one row per use that was read geometrically.
#[allow(clippy::type_complexity)]
pub(super) fn settle(
    definition: &NumericalDefinition,
    native_trims: &[NumericalNativeTrim],
    known: Option<&[(Rational, bool)]>,
    each_use: Limits,
    dec: &mut Decider,
) -> Result<(Option<Vec<NumericalNativeTrim>>, Vec<ImportGapRow>), Refusal> {
    let floor = tolerance::floor_minimum_mm();
    let floor_mm = floor.to_f64_nearest();
    let mut rows = Vec::new();
    let mut replaced: Vec<NumericalNativeTrim> = Vec::new();
    let mut index = 0;
    for face in &definition.faces {
        for (boundary, row) in face.loops.iter().enumerate() {
            for (use_index, trim) in row.iter().enumerate() {
                dec.charge_work(1)?;
                let record = known.and_then(|records| records.get(index));
                index += 1;
                let Some(edge) = definition.edges.get(trim.edge) else {
                    continue;
                };
                let binding = native_at(native_trims, &face.key, boundary, use_index);
                let tag = binding.map_or(&trim.tag, |binding| &binding.tag);
                if *tag != CopyTag::Independent
                    || record.is_some_and(|(stored, derived)| *derived || *stored <= floor)
                {
                    continue;
                }
                let copy = match binding {
                    Some(binding) => FaceCopyRef::Native(&binding.curve),
                    None if analytic::projected_edge(face, edge) => continue,
                    None => FaceCopyRef::Linear(&trim.uv),
                };
                let copies = EdgeCopies {
                    curve: &edge.curve,
                    parameters: edge_gap::traversal(trim, edge),
                    surface: GapSurface::Kernel(&face.surface),
                    face_copy: copy,
                };
                let Some(read) = read_use(&copies, each_use, floor_mm, dec)? else {
                    continue;
                };
                let mut found = ImportGapRow {
                    face: face.key.clone(),
                    edge: edge.key.clone(),
                    boundary,
                    use_index,
                    parametric_mm: exact(read.parametric.stored_mm())?,
                    geometric_mm: exact(read.geometric.reported_mm())?,
                    reparametrised_mm: None,
                    spans: 0,
                    outcome: read.outcome,
                };
                if let Some(next) = read.reparametrised {
                    found.reparametrised_mm = Some(exact(next.stored_mm)?);
                    found.spans = next.spans;
                    if read.outcome == ImportGapOutcome::Reparametrised {
                        replaced.push(NumericalNativeTrim {
                            face: face.key.clone(),
                            boundary,
                            use_index,
                            edge: edge.key.clone(),
                            curve: next.curve,
                            tag: CopyTag::Independent,
                        });
                    }
                }
                rows.push(found);
            }
        }
    }
    if replaced.is_empty() {
        return Ok((None, rows));
    }
    let mut next: Vec<NumericalNativeTrim> = native_trims
        .iter()
        .filter(|old| {
            !replaced.iter().any(|new| {
                new.face == old.face
                    && new.boundary == old.boundary
                    && new.use_index == old.use_index
            })
        })
        .cloned()
        .collect();
    next.extend(replaced);
    Ok((Some(next), rows))
}

/// Whether a refusal says that work or time ran out: the operation's
/// allowance, its deadline, or the one hard work cap of an edge (section 4.2
/// step 8: "A hard per-edge work cap refuses by name"). That is no verdict
/// on a copy. It is returned under its own name and is never answered with
/// another body.
pub(super) fn resource_limit(refusal: &Refusal) -> bool {
    refusal.code == gap::WORK_CAP
        || refusal.code.contains("budget")
        || refusal.code.contains("deadline")
}

/// The body an import admission returns once the gap rule has run.
///
/// `first` is the admission of the copies as imported, and `with_new_copies`
/// the admission with the reparametrised ones, where any were made. The new
/// copies win where they are admitted. Where they are refused by a resource
/// limit, the work cap of an edge among them, that refusal is the answer:
/// the body as imported was admitted before, and it is still not returned.
/// Where they are refused for any other reason, their rows say so and the
/// copies as imported stand: a body the rule cannot improve is not refused
/// for it. Where the copies as imported were over the cap, `refit` is rule
/// T3's refit of them.
pub(super) fn settled_body<B>(
    first: Result<B, Refusal>,
    with_new_copies: Option<Result<B, Refusal>>,
    rows: &mut [ImportGapRow],
    refit: impl FnOnce(Refusal) -> Result<B, Refusal>,
) -> Result<B, Refusal> {
    match with_new_copies {
        Some(Ok(body)) => return Ok(body),
        Some(Err(refusal)) if resource_limit(&refusal) => return Err(refusal),
        Some(Err(_)) => {
            for row in rows.iter_mut() {
                if row.outcome == ImportGapOutcome::Reparametrised {
                    row.outcome = ImportGapOutcome::NotAdmitted;
                }
            }
        }
        None => (),
    }
    first.or_else(refit)
}

/// The quarter cylinder of the tests under `tests/`.
#[cfg(test)]
#[path = "../../tests/common/analytic.rs"]
mod fixture;

#[cfg(test)]
mod tests {
    use super::*;
    use fixture::{q, quarter, stretched_copy};

    /// The small-allowance control through the admission itself (section
    /// 4.2 step 8: "A hard per-edge work cap refuses by name").
    ///
    /// One circle use of the quarter cylinder's wall has a copy at another
    /// speed. The body as imported is admitted, so an admission that fell
    /// back would have a body to return. The rule's three steps for that use
    /// are counted first. An admission that gives each use one unit less
    /// than their sum refuses with `numerical.gap-work-cap`, and so does
    /// every smaller allowance tried: no admission answers with the body as
    /// imported. An admission that gives each use exactly the sum admits,
    /// with the rows and the tolerances of the full cap.
    #[test]
    fn an_admission_never_answers_a_spent_allowance_with_the_body_as_imported() {
        let source = quarter();
        let (stretched, edge) = stretched_copy(&source, "0.001");
        let options = || NumericalOptions {
            tolerance_mm: q("0.01"),
            vertex_tolerance_mm: None,
        };
        let admit = |each_use: u64| {
            admit_import_within(
                source.clone(),
                vec![stretched.clone()],
                options(),
                Limits { work: each_use },
                &mut Decider::default(),
            )
        };
        build_with_native_trims(
            source.clone(),
            vec![stretched.clone()],
            tolerance::import_options(q("0.01")),
            &mut Decider::default(),
        )
        .expect("the body as imported is admitted");
        let full = admit(gap::EDGE_WORK).expect("the full cap");
        let [row] = full.import_gaps() else {
            panic!("{:?}", full.import_gaps())
        };
        assert_eq!(row.edge, edge);
        assert_eq!(row.outcome, ImportGapOutcome::Reparametrised);

        // The three steps of the rule for that use.
        let face = source
            .faces
            .iter()
            .find(|face| face.key == stretched.face)
            .expect("the wall");
        let trim = &face.loops[stretched.boundary][stretched.use_index];
        let used = &source.edges[trim.edge];
        let copies = EdgeCopies {
            curve: &used.curve,
            parameters: edge_gap::traversal(trim, used),
            surface: GapSurface::Kernel(&face.surface),
            face_copy: FaceCopyRef::Native(&stretched.curve),
        };
        let read = read_use(
            &copies,
            Limits::default(),
            tolerance::floor_minimum_mm().to_f64_nearest(),
            &mut Decider::default(),
        )
        .expect("the full cap")
        .expect("a use the rule reads");
        let next = read.reparametrised.as_ref().expect("a new copy");
        let steps = [read.parametric.work, read.geometric.work, next.work];
        let total: u64 = steps.iter().sum();

        let smaller = [
            0,
            steps[0],
            steps[0] + steps[1],
            steps[2],
            total / 2,
            total - 1,
        ];
        for each_use in smaller {
            match admit(each_use) {
                Err(refusal) => assert_eq!(
                    refusal.code, "numerical.gap-work-cap",
                    "{each_use} of {total} units"
                ),
                Ok(body) => panic!(
                    "{each_use} of {total} units admitted a body with {:?}",
                    body.import_gaps()
                ),
            }
        }
        let same = admit(total).expect("the sum of the steps is enough");
        assert_eq!(same.import_gaps(), full.import_gaps());
        assert_eq!(
            same.tolerances().largest_mm(),
            full.tolerances().largest_mm()
        );
        println!(
            "quarter cylinder, one copy at another speed: {steps:?} units, least allowance {total}"
        );
    }

    fn row(outcome: ImportGapOutcome) -> ImportGapRow {
        ImportGapRow {
            face: "f".into(),
            edge: "e".into(),
            boundary: 0,
            use_index: 0,
            parametric_mm: Rational::one(),
            geometric_mm: Rational::zero(),
            reparametrised_mm: None,
            spans: 0,
            outcome,
        }
    }

    fn refused(code: &'static str) -> Refusal {
        Refusal::new(code, "test")
    }

    fn never(_: Refusal) -> Result<&'static str, Refusal> {
        panic!("no refit")
    }

    /// The reparametrised copies win where they are admitted, whether or not
    /// the copies as imported were.
    #[test]
    fn admitted_new_copies_are_the_body() {
        for first in [Ok("imported"), Err(refused(TOLERANCE_EXCEEDED))] {
            let mut rows = [row(ImportGapOutcome::Reparametrised)];
            let body = settled_body(first, Some(Ok("new")), &mut rows, never);
            assert_eq!(body.unwrap(), "new");
            assert_eq!(rows[0].outcome, ImportGapOutcome::Reparametrised);
        }
    }

    /// New copies that are refused leave the body as imported, and their
    /// rows say so. Other rows keep what they said.
    #[test]
    fn refused_new_copies_leave_the_body_as_imported() {
        let mut rows = [
            row(ImportGapOutcome::Reparametrised),
            row(ImportGapOutcome::Geometric),
        ];
        let body = settled_body(
            Ok("imported"),
            Some(Err(refused("numerical.invalid-trim"))),
            &mut rows,
            never,
        );
        assert_eq!(body.unwrap(), "imported");
        assert_eq!(rows[0].outcome, ImportGapOutcome::NotAdmitted);
        assert_eq!(rows[1].outcome, ImportGapOutcome::Geometric);
    }

    /// Copies as imported that were over the cap go to rule T3's refit, with
    /// their own refusal, when no new copy was made or none was admitted.
    #[test]
    fn copies_over_the_cap_go_to_the_refit() {
        for with_new_copies in [None, Some(Err(refused("numerical.invalid-trim")))] {
            let mut rows = [row(ImportGapOutcome::NotLowered)];
            let body = settled_body(
                Err(refused(TOLERANCE_EXCEEDED)),
                with_new_copies,
                &mut rows,
                |refusal| {
                    assert_eq!(refusal.code, TOLERANCE_EXCEEDED);
                    Ok("refitted")
                },
            );
            assert_eq!(body.unwrap(), "refitted");
        }
    }

    /// A resource limit is no verdict on the new copies: it is the answer.
    #[test]
    fn a_resource_limit_on_the_new_copies_is_returned() {
        for code in ["numbers.operation-budget", "numbers.deadline"] {
            let mut rows = [row(ImportGapOutcome::Reparametrised)];
            let refusal = settled_body(Ok("imported"), Some(Err(refused(code))), &mut rows, never)
                .unwrap_err();
            assert_eq!(refusal.code, code);
            assert_eq!(rows[0].outcome, ImportGapOutcome::Reparametrised);
        }
    }

    /// Section 4.2 step 8: "A hard per-edge work cap refuses by name." Where
    /// the admission of the new copies hits the work cap of an edge, that
    /// refusal is returned by name. The body as imported was admitted, and
    /// it is not accepted in its place; no row is marked `NotAdmitted`; a
    /// body over the cap is not sent to rule T3's refit.
    #[test]
    fn the_work_cap_of_an_edge_is_returned_by_name_and_no_body_stands_in() {
        let capped = || Refusal::new(gap::WORK_CAP, "test");
        assert_eq!(capped().code, "numerical.gap-work-cap");
        assert!(resource_limit(&capped()));
        for first in [Ok("imported"), Err(refused(TOLERANCE_EXCEEDED))] {
            let mut rows = [row(ImportGapOutcome::Reparametrised)];
            let answer = settled_body(first, Some(Err(capped())), &mut rows, never);
            assert_eq!(answer.unwrap_err().code, "numerical.gap-work-cap");
            assert_eq!(rows[0].outcome, ImportGapOutcome::Reparametrised);
        }
        // A refusal that only shares words with the cap's name is a verdict
        // on the copy, as before.
        assert!(!resource_limit(&refused("numerical.gap-unsupported")));
    }
}
