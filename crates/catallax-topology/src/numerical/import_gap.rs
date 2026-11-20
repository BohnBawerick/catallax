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
//! name. An import admission has one account per edge use
//! ([`UseAccount`]). Admission's reading of the use draws on it first and
//! is kept there. The rule uses that reading and does not take it again;
//! its geometric reading and its change of parameter draw on what admission
//! left ([`read_use`]). The admission of the new copies uses the readings
//! already taken. Rule T3's refit, where a tolerance is still over the cap,
//! draws on the same accounts, and its admissions use the kept readings
//! too. The cap's refusal, `numerical.gap-work-cap`, is returned under its
//! own name wherever it comes from, and never answered with the body as
//! imported ([`settled_body`]).

use super::edge_gap::UseAccount;
use super::gap::{self, EdgeCopies, FaceCopyRef, GapSurface, Hints};
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
    /// Work counted against the one allowance of the use in this admission:
    /// admission's reading of the copy as imported, the geometric reading,
    /// the change of parameter and, where rule T3 refitted the use, that
    /// refit.
    pub work: u64,
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
    /// The same-parameter gap of the copy as imported: admission's reading.
    pub parametric: gap::GapMeasure,
    /// Its geometric gap.
    pub geometric: gap::GeometricMeasure,
    /// The copy in the parameter of the 3D curve, where the same-parameter
    /// gap was much larger than the geometric gap and the old copy has a
    /// change of parameter.
    pub reparametrised: Option<gap::refit::Refit>,
    /// What the rule did: `Geometric`, `NotLowered` or `Reparametrised`.
    pub outcome: ImportGapOutcome,
    /// All the work counted against the use's allowance so far: admission's
    /// reading, the geometric reading and the change of parameter.
    pub work: u64,
}

/// The rule for one edge use, on the use's account.
///
/// Section 4.2 step 8: "A hard per-edge work cap refuses by name." The
/// account holds the use's one allowance for the whole import admission,
/// and these draw on it in this order:
///
/// 1. The same-parameter gap of the copy as imported. Admission took this
///    reading, counted it and kept it in the account, and the rule uses it
///    as it stands. A use that admission did not read is read here, on the
///    same account.
/// 2. The geometric gap, both ways, under what the first reading left.
/// 3. The change of parameter, under what the two readings left: its
///    projections, its readings of the old copy and the checker's work on
///    every candidate.
///
/// Where the allowance runs out at any of these, the use refuses with
/// `numerical.gap-work-cap`. Nothing read before that stands in for the
/// answer: not the copy as imported, and not an earlier candidate.
///
/// Where the new copy is kept, the checker's reading of it is kept in the
/// account in place of the old copy's, so the admission of the new copy
/// does not read it again.
///
/// `None` for a use the rule leaves as it is: a straight residual, a copy
/// that leaves its surface's domain, a use that stores no more than
/// `floor_mm`, and a copy with no u/v of its own.
pub(super) fn read_use(
    copies: &EdgeCopies<'_>,
    account: &mut UseAccount,
    floor_mm: f64,
    dec: &mut Decider,
) -> Result<Option<UseReading>, Refusal> {
    let parametric = match account.reading() {
        Some(kept) => kept.clone(),
        None => {
            let found = gap::measure(copies, &Hints::default(), account.left(), dec)?;
            account.keep(std::slice::from_ref(&found));
            found
        }
    };
    if parametric.method == gap::Method::StraightResidual
        || parametric.domain_exit.is_some()
        || parametric.stored_mm() <= floor_mm
    {
        return Ok(None);
    }
    let Some(geometric) = gap::measure_geometric(copies, &Hints::default(), account.left(), dec)?
    else {
        return Ok(None);
    };
    account.spend(geometric.work);
    let mut reparametrised = None;
    let mut outcome = ImportGapOutcome::Geometric;
    if much_larger(parametric.reported_mm(), geometric.reported_mm()) {
        // Done when the new copy's same-parameter gap is no longer much
        // larger than the geometric gap, or stores the floor.
        let target = (gap::MARGIN * gap::MARGIN * geometric.reported_mm()).max(floor_mm);
        outcome = ImportGapOutcome::NotLowered;
        let change = gap::refit::reparametrise(copies, target, account.left(), dec)?;
        account.spend(change.work);
        reparametrised = change.copy;
        if let Some(next) = &reparametrised
            && much_larger(parametric.stored_mm(), next.stored_mm)
        {
            outcome = ImportGapOutcome::Reparametrised;
            account.replace(next.reading.clone());
        }
    }
    Ok(Some(UseReading {
        parametric,
        geometric,
        reparametrised,
        outcome,
        work: account.spent(),
    }))
}

/// What the rule did with the uses of a definition.
pub(super) struct Settled {
    /// The retained u/v curves with the new copies in place, or `None` where
    /// no copy changed.
    pub native_trims: Option<Vec<NumericalNativeTrim>>,
    /// One row per use that was read geometrically.
    pub rows: Vec<ImportGapRow>,
    /// The place of each row's use in definition order, where its account
    /// is.
    pub slots: Vec<usize>,
}

/// Reads every imported edge use the gap rule applies to and gives the uses
/// whose copies run at different speeds a face copy in the parameter of
/// their 3D curve.
///
/// `accounts` holds the account of every use, in definition order, with the
/// reading admission took of it. A use whose reading stores no more than
/// rule T1's least floor is not read geometrically: no tolerance is stored
/// under the floor, so nothing could be lowered.
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
/// The rule draws on each use's account, after admission ([`read_use`]). A
/// use that spends its allowance refuses the admission with
/// `numerical.gap-work-cap`.
///
pub(super) fn settle(
    definition: &NumericalDefinition,
    native_trims: &[NumericalNativeTrim],
    accounts: &mut [UseAccount],
    dec: &mut Decider,
) -> Result<Settled, Refusal> {
    let floor_mm = tolerance::floor_minimum_mm().to_f64_nearest();
    let mut rows = Vec::new();
    let mut slots = Vec::new();
    let mut replaced: Vec<NumericalNativeTrim> = Vec::new();
    let mut index = 0;
    for face in &definition.faces {
        for (boundary, row) in face.loops.iter().enumerate() {
            for (use_index, trim) in row.iter().enumerate() {
                dec.charge_work(1)?;
                let slot = index;
                index += 1;
                let Some(edge) = definition.edges.get(trim.edge) else {
                    continue;
                };
                let binding = native_at(native_trims, &face.key, boundary, use_index);
                let tag = binding.map_or(&trim.tag, |binding| &binding.tag);
                if *tag != CopyTag::Independent {
                    continue;
                }
                let account = accounts
                    .get_mut(slot)
                    .ok_or_else(|| bad("definition", "an edge use has no work account"))?;
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
                let Some(read) = read_use(&copies, account, floor_mm, dec)? else {
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
                    work: read.work,
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
                slots.push(slot);
            }
        }
    }
    if replaced.is_empty() {
        return Ok(Settled {
            native_trims: None,
            rows,
            slots,
        });
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
    Ok(Settled {
        native_trims: Some(next),
        rows,
        slots,
    })
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
    use crate::numerical::gap::Limits;
    use fixture::{bowed_copy, q, quarter, stretched_copy};

    fn accounts(source: &NumericalDefinition, each_use: u64) -> Vec<UseAccount> {
        UseAccount::for_uses(source, Limits { work: each_use })
    }

    /// The import admission of the quarter cylinder with one retained u/v
    /// curve, at the automatic cap, on accounts of `each_use` units, and the
    /// accounts it leaves.
    fn admit(
        source: &NumericalDefinition,
        copy: &NumericalNativeTrim,
        each_use: u64,
    ) -> (Result<NumericalBody, Refusal>, Vec<UseAccount>) {
        let mut accounts = accounts(source, each_use);
        let admitted = admit_import_on(
            source.clone(),
            vec![copy.clone()],
            NumericalOptions {
                tolerance_mm: q("0.01"),
                vertex_tolerance_mm: None,
            },
            &mut accounts,
            &mut Decider::default(),
        );
        (admitted, accounts)
    }

    /// The body as imported alone, on accounts of `each_use` units.
    fn as_imported(
        source: &NumericalDefinition,
        copy: &NumericalNativeTrim,
        each_use: u64,
    ) -> (Result<NumericalBody, Refusal>, Vec<UseAccount>) {
        let mut accounts = accounts(source, each_use);
        let admitted = admit_on(
            source.clone(),
            vec![copy.clone()],
            tolerance::import_options(q("0.01")),
            &mut accounts,
            &mut Decider::default(),
        );
        (admitted, accounts)
    }

    /// The use a retained u/v curve belongs to: its place in definition
    /// order and its two copies.
    fn use_of<'a>(
        source: &'a NumericalDefinition,
        copy: &'a NumericalNativeTrim,
    ) -> (usize, EdgeCopies<'a>) {
        let wall = source
            .faces
            .iter()
            .position(|face| face.key == copy.face)
            .expect("the wall");
        let face = &source.faces[wall];
        let trim = &face.loops[copy.boundary][copy.use_index];
        let used = &source.edges[trim.edge];
        (
            use_slot(source, wall, copy.boundary, copy.use_index),
            EdgeCopies {
                curve: &used.curve,
                parameters: edge_gap::traversal(trim, used),
                surface: GapSurface::Kernel(&face.surface),
                face_copy: FaceCopyRef::Native(&copy.curve),
            },
        )
    }

    fn refuses_by_name(admitted: Result<NumericalBody, Refusal>, at: &str) {
        match admitted {
            Err(refusal) => assert_eq!(refusal.code, "numerical.gap-work-cap", "{at}"),
            Ok(body) => panic!("{at}: admitted a body with {:?}", body.import_gaps()),
        }
    }

    /// The small-allowance control through the admission itself (section
    /// 4.2 step 8: "A hard per-edge work cap refuses by name").
    ///
    /// One circle use of the quarter cylinder's wall has a copy at another
    /// speed. The use has one account for the whole admission. Admission's
    /// own reading of the use is the first to draw on it, then the rule's
    /// geometric reading and its change of parameter.
    ///
    /// The three are counted apart from the admission first: the checker's
    /// reading of the copy as imported, which is the reading admission
    /// takes, the geometric reading, and the change of parameter. The
    /// admission at the full cap counts exactly their sum for the use, also
    /// after the body is admitted with the new copy. So admission's reading
    /// is in the account once, the rule did not take it again, and the new
    /// copy was not read again. Every other use is counted once.
    ///
    /// The boundary includes admission's reading. An allowance that covers
    /// the rule's two steps and not admission's reading before them refuses
    /// with `numerical.gap-work-cap`, and so does one unit less than the sum
    /// of all three. Under both, the body as imported is admitted on its
    /// own, so an admission that fell back would have had a body to return.
    /// Exactly the sum admits, with the rows and tolerances of the full cap.
    #[test]
    fn one_allowance_covers_admission_and_the_rule_and_no_body_stands_in() {
        let source = quarter();
        let (stretched, edge) = stretched_copy(&source, "0.001");
        let (slot, copies) = use_of(&source, &stretched);

        // The three, counted apart from the admission.
        let floor_mm = tolerance::floor_minimum_mm().to_f64_nearest();
        let mut dec = Decider::default();
        let first = gap::measure(&copies, &Hints::default(), Limits::default(), &mut dec)
            .expect("the copy as imported")
            .work;
        let geometric =
            gap::measure_geometric(&copies, &Hints::default(), Limits::default(), &mut dec)
                .expect("the geometric gap")
                .expect("a copy with its own u/v");
        let target = (gap::MARGIN * gap::MARGIN * geometric.reported_mm()).max(floor_mm);
        let change = gap::refit::reparametrise(&copies, target, Limits::default(), &mut dec)
            .expect("the change of parameter")
            .work;
        let rule_alone = geometric.work + change;
        let total = first + rule_alone;

        let (full, counted) = admit(&source, &stretched, gap::EDGE_WORK);
        let full = full.expect("the full cap");
        let [row] = full.import_gaps() else {
            panic!("{:?}", full.import_gaps())
        };
        assert_eq!(row.edge, edge);
        assert_eq!(row.outcome, ImportGapOutcome::Reparametrised);
        assert_eq!(row.work, total);
        assert_eq!(counted[slot].spent(), total);
        let (imported, read_once) = as_imported(&source, &stretched, gap::EDGE_WORK);
        imported.expect("the body as imported is admitted");
        assert_eq!(read_once[slot].spent(), first);
        for (index, (after, once)) in counted.iter().zip(&read_once).enumerate() {
            if index != slot {
                assert_eq!(after.spent(), once.spent(), "use {index}");
            }
        }
        // The body admitted on kept readings stores what a fresh admission
        // of the same copies stores.
        let fresh = build_with_native_trims(
            source.clone(),
            full.native_trims().to_vec(),
            tolerance::import_options(q("0.01")),
            &mut Decider::default(),
        )
        .expect("the new copies on their own");
        assert_eq!(fresh.evidence().checks(), full.evidence().checks());
        assert_eq!(fresh.tolerances(), full.tolerances());

        // The boundary. With these allowances the body as imported is
        // admitted, and the admission refuses by name all the same.
        for each_use in [rule_alone, total - 1] {
            as_imported(&source, &stretched, each_use)
                .0
                .expect("the body as imported is admitted under this allowance");
            refuses_by_name(
                admit(&source, &stretched, each_use).0,
                &format!("{each_use} of {total} units"),
            );
        }
        // With less than admission's own reading needs, admission itself
        // refuses by name.
        for each_use in [0, first - 1] {
            refuses_by_name(
                admit(&source, &stretched, each_use).0,
                &format!("{each_use} units"),
            );
        }
        let same = admit(&source, &stretched, total)
            .0
            .expect("the sum is enough");
        assert_eq!(same.import_gaps(), full.import_gaps());
        assert_eq!(same.tolerances(), full.tolerances());
        println!(
            "quarter cylinder, one copy at another speed: admission's reading {first}, geometric reading {}, change of parameter {change}, least allowance {total}",
            geometric.work
        );
    }

    /// The work boundary of rule T3's refit in an import (section 4.2 step
    /// 8: "A hard per-edge work cap refuses by name").
    ///
    /// One circle use of the wall has a copy bowed off its circle: a gap in
    /// the geometry of the copy, over the automatic cap. The rule reads it
    /// and leaves it, the body as imported is over the cap, and rule T3
    /// refits the copy from its 3D curve.
    ///
    /// The use's account covers all of it. Counted apart from the admission:
    /// the checker's reading of the copy as imported, the geometric reading,
    /// and the refit under the cap. The admission at the full cap counts
    /// exactly their sum for the use and counts every other use once. So
    /// the probe at the opt-in limit and the admissions after the refit
    /// read nothing again, and the refit did not start a new cap.
    ///
    /// With the sum the import is admitted, and the body stores what a
    /// fresh admission of the refitted copy stores. With one unit less it
    /// refuses with `numerical.gap-work-cap`. With what admission and the
    /// rule need and nothing for the refit it refuses by name as well: the
    /// answer is neither a body nor the tolerance refusal with its opt-in.
    #[test]
    fn rule_t3s_refit_draws_on_the_account_of_its_use() {
        let source = quarter();
        let (bowed, edge) = bowed_copy(&source, "0.03");
        let (slot, copies) = use_of(&source, &bowed);

        // The three, counted apart from the admission.
        let mut dec = Decider::default();
        let first = gap::measure(&copies, &Hints::default(), Limits::default(), &mut dec)
            .expect("the copy as imported")
            .work;
        let geometric =
            gap::measure_geometric(&copies, &Hints::default(), Limits::default(), &mut dec)
                .expect("the geometric gap")
                .expect("a copy with its own u/v")
                .work;
        let cap_mm = q("0.01").to_f64_nearest();
        let refit = gap::refit::refit_within(&copies, cap_mm, Limits::default(), &mut dec)
            .expect("the refit of rule T3");
        assert!(refit.copy.is_some());
        let before_refit = first + geometric;
        let total = before_refit + refit.work;

        // The body as imported is over the cap, and admission reads every
        // use once before it says so.
        let (imported, read_once) = as_imported(&source, &bowed, gap::EDGE_WORK);
        assert_eq!(imported.unwrap_err().code, TOLERANCE_EXCEEDED);
        assert_eq!(read_once[slot].spent(), first);

        let (full, counted) = admit(&source, &bowed, gap::EDGE_WORK);
        let full = full.expect("the full cap");
        let [row] = full.import_gaps() else {
            panic!("{:?}", full.import_gaps())
        };
        assert_eq!(row.edge, edge);
        assert_eq!(row.outcome, ImportGapOutcome::Geometric);
        assert_eq!(row.work, total);
        assert_eq!(counted[slot].spent(), total);
        for (index, (after, once)) in counted.iter().zip(&read_once).enumerate() {
            if index != slot {
                assert_eq!(after.spent(), once.spent(), "use {index}");
            }
        }
        assert!(full.tolerances().largest_mm() <= q("0.01"));
        let fresh = build_with_native_trims(
            source.clone(),
            full.native_trims().to_vec(),
            tolerance::import_options(q("0.01")),
            &mut Decider::default(),
        )
        .expect("the refitted copy on its own");
        assert_eq!(fresh.evidence().checks(), full.evidence().checks());
        assert_eq!(fresh.tolerances(), full.tolerances());

        // The boundary.
        for each_use in [before_refit, total - 1] {
            refuses_by_name(
                admit(&source, &bowed, each_use).0,
                &format!("{each_use} of {total} units"),
            );
        }
        let same = admit(&source, &bowed, total).0.expect("the sum is enough");
        assert_eq!(same.import_gaps(), full.import_gaps());
        assert_eq!(same.tolerances(), full.tolerances());
        println!(
            "quarter cylinder, one bowed copy: admission's reading {first}, geometric reading {geometric}, refit of rule T3 {}, least allowance {total}",
            refit.work
        );
    }

    /// A new copy that is refused sends its use back to the copy as
    /// imported, on the same account.
    ///
    /// The wall is 0.02 mm wider than its circles, over the automatic cap,
    /// and one copy also runs up to 2% of the arc ahead. The rule gives that
    /// copy the circle's parameter. The new copy still stores the 0.02 mm
    /// that is in the geometry, rule T3's refit cannot lower that, and the
    /// new copy is refused. Rule T3 then runs on the copies as imported and
    /// refuses in turn, with the refusal of the body as imported and the
    /// opt-in.
    ///
    /// The account of the use holds the reading of its copy as imported
    /// again, and everything that was counted for the use on the way is
    /// still counted: admission's reading, the geometric reading, the
    /// change of parameter and both refits.
    #[test]
    fn a_refused_new_copy_sends_its_use_back_to_the_copy_as_imported() {
        let mut source = quarter();
        let Surface::Cylinder(cylinder) = &mut source.faces[4].surface else {
            unreachable!()
        };
        cylinder.radius_sq = (&q("1.02") * &q("1.02")).into();
        let (stretched, _) = stretched_copy(&source, "0.02");
        let (slot, copies) = use_of(&source, &stretched);

        let mut dec = Decider::default();
        let floor_mm = tolerance::floor_minimum_mm().to_f64_nearest();
        let first = gap::measure(&copies, &Hints::default(), Limits::default(), &mut dec)
            .expect("the copy as imported");
        let geometric =
            gap::measure_geometric(&copies, &Hints::default(), Limits::default(), &mut dec)
                .expect("the geometric gap")
                .expect("a copy with its own u/v");
        assert!(much_larger(first.reported_mm(), geometric.reported_mm()));
        let target = (gap::MARGIN * gap::MARGIN * geometric.reported_mm()).max(floor_mm);
        let change = gap::refit::reparametrise(&copies, target, Limits::default(), &mut dec)
            .expect("the change of parameter");
        let next = change.copy.as_ref().expect("a new copy");
        assert!(much_larger(first.stored_mm(), next.stored_mm));
        let cap_mm = q("0.01").to_f64_nearest();
        assert!(next.stored_mm > cap_mm);
        // Rule T3's refit of the new copy, then of the copy as imported.
        let of_new = gap::refit::refit_within(
            &EdgeCopies {
                face_copy: FaceCopyRef::Native(&next.curve),
                ..copies
            },
            cap_mm,
            Limits::default(),
            &mut dec,
        )
        .expect("the refit of the new copy");
        let of_imported = gap::refit::refit_within(&copies, cap_mm, Limits::default(), &mut dec)
            .expect("the refit of the copy as imported");
        let total = first.work + geometric.work + change.work + of_new.work + of_imported.work;

        let (imported, _) = as_imported(&source, &stretched, gap::EDGE_WORK);
        let over_the_cap = imported.unwrap_err();
        assert_eq!(over_the_cap.code, TOLERANCE_EXCEEDED);

        let (admitted, counted) = admit(&source, &stretched, gap::EDGE_WORK);
        let refusal = admitted.map(|_| ()).expect_err("over the cap");
        assert_eq!(refusal.code, TOLERANCE_EXCEEDED);
        assert_eq!(refusal.detail, opt_in_refusal(over_the_cap).detail);
        let kept = counted[slot].reading().expect("a kept reading");
        assert_eq!(kept.raw_mm.to_bits(), first.raw_mm.to_bits());
        assert_eq!(kept.work, first.work);
        assert_eq!(counted[slot].spent(), total);

        // The other circle use of the wider wall is read and refitted too,
        // so the boundary is taken on this use's account alone: the other
        // uses keep the full cap. One unit less than all of it refuses by
        // name, at the last refit, and exactly all of it gives the answer
        // of the full cap.
        let with = |allowance: u64| {
            let mut accounts = accounts(&source, gap::EDGE_WORK);
            accounts[slot] = UseAccount::new(Limits { work: allowance });
            admit_import_on(
                source.clone(),
                vec![stretched.clone()],
                NumericalOptions {
                    tolerance_mm: q("0.01"),
                    vertex_tolerance_mm: None,
                },
                &mut accounts,
                &mut Decider::default(),
            )
        };
        refuses_by_name(with(total - 1), &format!("{} of {total} units", total - 1));
        assert_eq!(with(total).map(|_| ()).expect_err("over the cap"), refusal);
        println!(
            "wider wall, one copy at another speed: admission's reading {}, geometric reading {}, change of parameter {}, refit of the new copy {}, refit of the copy as imported {}, least allowance {total}",
            first.work, geometric.work, change.work, of_new.work, of_imported.work
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
            work: 0,
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
