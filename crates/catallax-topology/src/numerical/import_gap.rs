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
//! name. An import admission has one allowance per edge, shared by the
//! edge's uses ([`Accounts`]). Admission's reading of a use draws on it
//! first and is kept. The rule uses that reading and does not take it
//! again; its geometric reading and its change of parameter draw on what is
//! left of the edge's allowance ([`read_use`]). The admission of the new
//! copies uses the readings already taken. Rule T3's refit, where a
//! tolerance is still over the cap, draws on the same accounts, and its
//! admissions use the kept readings too. The cap's refusal,
//! `numerical.gap-work-cap`, is returned under its own name wherever it
//! comes from, and never answered with the body as imported
//! ([`settled_body`]).

use super::edge_gap::{Accounts, UseDraw};
use super::gap::{self, EdgeCopies, FaceCopyRef, GapCurve, GapSurface, Hints};
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
    /// Work this use counted against the one allowance of its edge in this
    /// admission: admission's reading of the copy as imported, the geometric
    /// reading, the change of parameter and, where rule T3 refitted the use,
    /// that refit. The edge's other uses draw on the same allowance.
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
    /// All the work this use has counted so far: admission's reading, the
    /// geometric reading and the change of parameter.
    pub work: u64,
}

/// The rule for one edge use, on the allowance of the use's edge.
///
/// Section 4.2 step 8: "A hard per-edge work cap refuses by name." The edge
/// has one allowance for the whole import admission, shared by its uses,
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
/// Each step's work is taken off the edge's allowance as the step counts
/// it, also where the step is refused ([`UseDraw::draw_on`]).
///
/// Where the allowance runs out at any of these, the edge refuses with
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
    account: &mut UseDraw<'_>,
    floor_mm: f64,
    dec: &mut Decider,
) -> Result<Option<UseReading>, Refusal> {
    let parametric = match account.reading() {
        Some(kept) => kept.clone(),
        None => {
            let found = account.draw_on(|limits, spent| {
                gap::measure_counted(copies, &Hints::default(), limits, spent, dec)
            })?;
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
    let Some(geometric) = account.draw_on(|limits, spent| {
        gap::measure_geometric_counted(copies, &Hints::default(), limits, spent, dec)
    })?
    else {
        return Ok(None);
    };
    let mut reparametrised = None;
    let mut outcome = ImportGapOutcome::Geometric;
    if much_larger(parametric.reported_mm(), geometric.reported_mm()) {
        // Done when the new copy's same-parameter gap is no longer much
        // larger than the geometric gap, or stores the floor.
        let target = (gap::MARGIN * gap::MARGIN * geometric.reported_mm()).max(floor_mm);
        outcome = ImportGapOutcome::NotLowered;
        let change = account.draw_on(|limits, spent| {
            gap::refit::reparametrise_counted(copies, target, limits, spent, dec)
        })?;
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

/// What the gap rule read of one imported edge use that a reader holds in
/// its own forms ([`read_import_use`]).
#[derive(Clone, Debug)]
pub struct ImportUseReading {
    /// The same-parameter gap of the copy as imported, read by the A+
    /// checker.
    pub parametric: gap::GapMeasure,
    /// Its geometric gap, both ways. `None` for a use the rule leaves as it
    /// is: a straight residual, a copy that leaves its surface's domain, a
    /// use that stores no more than rule T1's least floor, and a copy with
    /// no u/v of its own.
    pub geometric: Option<gap::GeometricMeasure>,
    /// The copy in the parameter of the 3D curve, where the same-parameter
    /// gap was much larger than the geometric gap and the old copy has a
    /// change of parameter. It is the use's copy only where `outcome` is
    /// [`ImportGapOutcome::Reparametrised`].
    pub reparametrised: Option<gap::refit::Refit>,
    /// What the rule did: `Geometric`, `NotLowered` or `Reparametrised`.
    /// `None` for a use the rule leaves as it is.
    pub outcome: Option<ImportGapOutcome>,
    /// All the work this use counted against the allowance it was given:
    /// the same-parameter reading, the geometric reading and the change of
    /// parameter.
    pub work: u64,
}

impl ImportUseReading {
    /// The copy the rule gave the use in place of the imported one.
    pub fn new_copy(&self) -> Option<&gap::refit::Refit> {
        self.reparametrised
            .as_ref()
            .filter(|_| self.outcome == Some(ImportGapOutcome::Reparametrised))
    }

    /// The tolerance rule T4 stores for the copy the rule keeps: 1.1 times
    /// its same-parameter gap, as the full checker read it.
    pub fn stored_mm(&self) -> f64 {
        self.new_copy()
            .map_or_else(|| self.parametric.stored_mm(), |copy| copy.stored_mm)
    }
}

/// The gap rule for one imported edge use (decision Q9, and section 4.2 step
/// 2 of the Checkpoint 2 decisions), for a reader that holds an import in
/// forms of its own and admits it use by use.
///
/// This is the rule [`admit_import`] runs on every use of a numerical
/// definition, on one use: the A+ checker reads the same-parameter gap of
/// the copy as imported; where that could be lowered, the geometric gap is
/// read both ways; where the same-parameter gap is much larger, the face
/// copy is given the parameter of its 3D curve, and the checker's reading of
/// the new copy is what the use continues with.
///
/// `limits` is what the use's edge has left of its one work allowance
/// (section 4.2 step 8). The allowance belongs to the edge, so a reader
/// gives the first use of an edge the checker's cap, takes the returned
/// `work` off it, and gives each later use of the same edge what is left.
/// All three steps draw on `limits` in turn, every unit is charged to the
/// operation, and a use that spends it refuses with
/// `numerical.gap-work-cap`.
pub fn read_import_use(
    copies: &EdgeCopies<'_>,
    limits: gap::Limits,
    dec: &mut Decider,
) -> Result<ImportUseReading, Refusal> {
    read_import_use_counted(copies, limits, &mut 0, dec)
}

/// [`read_import_use`], with all the work the use counted added to `spent`
/// whether the rule finished or was refused. On a finished reading that is
/// the reading's `work`. A reader that goes on after a refusal that is not a
/// resource limit, to another use of the same edge, takes `spent` off the
/// edge's allowance either way.
pub fn read_import_use_counted(
    copies: &EdgeCopies<'_>,
    limits: gap::Limits,
    spent: &mut u64,
    dec: &mut Decider,
) -> Result<ImportUseReading, Refusal> {
    let mut accounts = Accounts::of_one_use(limits);
    let read = read_one_use(copies, &mut accounts, dec);
    *spent = spent.saturating_add(accounts.spent_by_use(0));
    read
}

fn read_one_use(
    copies: &EdgeCopies<'_>,
    accounts: &mut Accounts,
    dec: &mut Decider,
) -> Result<ImportUseReading, Refusal> {
    let mut account = accounts.draw(0)?;
    let floor_mm = tolerance::floor_minimum_mm().to_f64_nearest();
    Ok(match read_use(copies, &mut account, floor_mm, dec)? {
        Some(read) => ImportUseReading {
            parametric: read.parametric,
            geometric: Some(read.geometric),
            reparametrised: read.reparametrised,
            outcome: Some(read.outcome),
            work: read.work,
        },
        None => ImportUseReading {
            parametric: account
                .reading()
                .cloned()
                .ok_or_else(|| bad("trim", "an edge use has no gap reading"))?,
            geometric: None,
            reparametrised: None,
            outcome: None,
            work: account.spent(),
        },
    })
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
/// `accounts` holds the allowance of every edge and the record of every
/// use, in definition order, with the reading admission took of it. A use whose reading stores no more than
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
/// The rule draws on the allowance of each use's edge, after admission
/// ([`read_use`]). An edge whose uses spend its allowance refuses the
/// admission with `numerical.gap-work-cap`.
pub(super) fn settle(
    definition: &NumericalDefinition,
    native_trims: &[NumericalNativeTrim],
    accounts: &mut Accounts,
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
                let mut account = accounts.draw(slot)?;
                let copy = match binding {
                    Some(binding) => FaceCopyRef::Native(&binding.curve),
                    None if analytic::projected_edge(face, edge) => continue,
                    None => FaceCopyRef::Linear(&trim.uv),
                };
                let copies = EdgeCopies {
                    curve: GapCurve::Kernel(&edge.curve),
                    parameters: edge_gap::traversal(trim, edge),
                    surface: GapSurface::Kernel(&face.surface),
                    face_copy: copy,
                };
                let Some(read) = read_use(&copies, &mut account, floor_mm, dec)? else {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numerical::gap::Limits;
    use crate::numerical::test_fixture::{bowed_copy, q, quarter, stretched_copy};

    /// The accounts of `source` with the checker's cap for every edge, and
    /// `work` units for the edge at `edge`.
    fn accounts(source: &NumericalDefinition, edge: usize, work: u64) -> Accounts {
        let mut accounts = Accounts::for_definition(source, Limits::default());
        accounts.limit_edge(edge, work);
        accounts
    }

    /// The edge of the use a retained u/v curve belongs to, as its place in
    /// the edge table.
    fn edge_of(source: &NumericalDefinition, copy: &NumericalNativeTrim) -> usize {
        source
            .edges
            .iter()
            .position(|edge| edge.key == copy.edge)
            .expect("the edge of the copy")
    }

    /// The uses of the edge at `edge`, as their places in definition order.
    fn uses_of(source: &NumericalDefinition, edge: usize) -> Vec<usize> {
        source
            .faces
            .iter()
            .flat_map(|face| face.loops.iter().flatten())
            .enumerate()
            .filter(|(_, trim)| trim.edge == edge)
            .map(|(slot, _)| slot)
            .collect()
    }

    /// The import admission of the quarter cylinder with one retained u/v
    /// curve, at the automatic cap, with `work` units for the edge at `edge`
    /// and the checker's cap for every other, and the accounts it leaves.
    fn admit(
        source: &NumericalDefinition,
        copy: &NumericalNativeTrim,
        edge: usize,
        work: u64,
    ) -> (Result<NumericalBody, Refusal>, Accounts) {
        let mut accounts = accounts(source, edge, work);
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

    /// The body as imported alone, with `work` units for the edge at `edge`.
    fn as_imported(
        source: &NumericalDefinition,
        copy: &NumericalNativeTrim,
        edge: usize,
        work: u64,
    ) -> (Result<NumericalBody, Refusal>, Accounts) {
        let mut accounts = accounts(source, edge, work);
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
                curve: GapCurve::Kernel(&used.curve),
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
    /// speed. Its edge has one allowance for the whole admission, which the
    /// edge's other use, on the cap, draws on as well. Admission's own
    /// reading of the use draws on it, then the rule's geometric reading and
    /// its change of parameter.
    ///
    /// The three are counted apart from the admission first: the checker's
    /// reading of the copy as imported, which is the reading admission
    /// takes, the geometric reading, and the change of parameter. The
    /// admission at the full cap counts exactly their sum for the use, also
    /// after the body is admitted with the new copy. So admission's reading
    /// is in the account once, the rule did not take it again, and the new
    /// copy was not read again. Every other use is counted once.
    ///
    /// The boundary includes admission's reading and the edge's other use.
    /// With what the other use counts set aside, an allowance that covers
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

        let at = edge_of(&source, &stretched);
        let (full, counted) = admit(&source, &stretched, at, gap::EDGE_WORK);
        let full = full.expect("the full cap");
        let [row] = full.import_gaps() else {
            panic!("{:?}", full.import_gaps())
        };
        assert_eq!(row.edge, edge);
        assert_eq!(row.outcome, ImportGapOutcome::Reparametrised);
        assert_eq!(row.work, total);
        assert_eq!(counted.spent_by_use(slot), total);
        // What the edge's other use counted on the same allowance.
        let others = counted.spent_by_edge(at) - total;
        assert!(others > 0);
        let (imported, read_once) = as_imported(&source, &stretched, at, gap::EDGE_WORK);
        imported.expect("the body as imported is admitted");
        assert_eq!(read_once.spent_by_use(slot), first);
        for index in (0..counted.uses()).filter(|index| *index != slot) {
            assert_eq!(
                counted.spent_by_use(index),
                read_once.spent_by_use(index),
                "use {index}"
            );
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
        for work in [others + rule_alone, others + total - 1] {
            as_imported(&source, &stretched, at, work)
                .0
                .expect("the body as imported is admitted under this allowance");
            refuses_by_name(
                admit(&source, &stretched, at, work).0,
                &format!("{work} of {} units", others + total),
            );
        }
        // With less than admission's own readings need, admission itself
        // refuses by name.
        for work in [0, others + first - 1] {
            refuses_by_name(
                admit(&source, &stretched, at, work).0,
                &format!("{work} units"),
            );
        }
        let same = admit(&source, &stretched, at, others + total)
            .0
            .expect("the sum is enough");
        assert_eq!(same.import_gaps(), full.import_gaps());
        assert_eq!(same.tolerances(), full.tolerances());
        println!(
            "quarter cylinder, one copy at another speed: admission's reading {first}, geometric reading {}, change of parameter {change}, the edge's other use {others}, least allowance of the edge {}",
            geometric.work,
            others + total
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
    /// The allowance of the use's edge covers all of it. Counted apart from
    /// the admission: the checker's reading of the copy as imported, the
    /// geometric reading, and the refit under the cap. The admission at the
    /// full cap counts exactly their sum for the use and counts every other
    /// use once. So the probe at the opt-in limit and the admissions after
    /// the refit read nothing again, and the refit did not start a new cap.
    ///
    /// With the sum, and what the edge's other use counts, the import is
    /// admitted, and the body stores what a fresh admission of the refitted
    /// copy stores. With one unit less it refuses with
    /// `numerical.gap-work-cap`. With what admission and the rule need and
    /// nothing for the refit it refuses by name as well: the answer is
    /// neither a body nor the tolerance refusal with its opt-in.
    #[test]
    fn rule_t3s_refit_draws_on_the_allowance_of_its_edge() {
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
        let at = edge_of(&source, &bowed);
        let (imported, read_once) = as_imported(&source, &bowed, at, gap::EDGE_WORK);
        assert_eq!(imported.unwrap_err().code, TOLERANCE_EXCEEDED);
        assert_eq!(read_once.spent_by_use(slot), first);

        let (full, counted) = admit(&source, &bowed, at, gap::EDGE_WORK);
        let full = full.expect("the full cap");
        let [row] = full.import_gaps() else {
            panic!("{:?}", full.import_gaps())
        };
        assert_eq!(row.edge, edge);
        assert_eq!(row.outcome, ImportGapOutcome::Geometric);
        assert_eq!(row.work, total);
        assert_eq!(counted.spent_by_use(slot), total);
        let others = counted.spent_by_edge(at) - total;
        for index in (0..counted.uses()).filter(|index| *index != slot) {
            assert_eq!(
                counted.spent_by_use(index),
                read_once.spent_by_use(index),
                "use {index}"
            );
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
        for work in [others + before_refit, others + total - 1] {
            refuses_by_name(
                admit(&source, &bowed, at, work).0,
                &format!("{work} of {} units", others + total),
            );
        }
        let same = admit(&source, &bowed, at, others + total)
            .0
            .expect("the sum is enough");
        assert_eq!(same.import_gaps(), full.import_gaps());
        assert_eq!(same.tolerances(), full.tolerances());
        println!(
            "quarter cylinder, one bowed copy: admission's reading {first}, geometric reading {geometric}, refit of rule T3 {}, the edge's other use {others}, least allowance of the edge {}",
            refit.work,
            others + total
        );
    }

    /// A new copy that is refused sends its use back to the copy as
    /// imported, on the same allowance.
    ///
    /// The wall is 0.02 mm wider than its circles, over the automatic cap,
    /// and one copy also runs up to 2% of the arc ahead. The rule gives that
    /// copy the circle's parameter. The new copy still stores the 0.02 mm
    /// that is in the geometry, rule T3's refit cannot lower that, and the
    /// new copy is refused. Rule T3 then runs on the copies as imported and
    /// refuses in turn, with the refusal of the body as imported and the
    /// opt-in.
    ///
    /// The record of the use holds the reading of its copy as imported
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

        let at = edge_of(&source, &stretched);
        let (imported, _) = as_imported(&source, &stretched, at, gap::EDGE_WORK);
        let over_the_cap = imported.unwrap_err();
        assert_eq!(over_the_cap.code, TOLERANCE_EXCEEDED);

        let (admitted, counted) = admit(&source, &stretched, at, gap::EDGE_WORK);
        let refusal = admitted.map(|_| ()).expect_err("over the cap");
        assert_eq!(refusal.code, TOLERANCE_EXCEEDED);
        assert_eq!(refusal.detail, opt_in_refusal(over_the_cap).detail);
        let kept = counted.reading_of(slot).expect("a kept reading");
        assert_eq!(kept.raw_mm.to_bits(), first.raw_mm.to_bits());
        assert_eq!(kept.work, first.work);
        assert_eq!(counted.spent_by_use(slot), total);
        let others = counted.spent_by_edge(at) - total;

        // The other circle edge of the wider wall is read and refitted too,
        // so the boundary is taken on this edge's allowance alone: the other
        // edges keep the full cap. One unit less than all of it refuses by
        // name, at the last refit, and exactly all of it gives the answer
        // of the full cap.
        let with = |work: u64| admit(&source, &stretched, at, work).0;
        refuses_by_name(
            with(others + total - 1),
            &format!("{} of {} units", others + total - 1, others + total),
        );
        assert_eq!(
            with(others + total).map(|_| ()).expect_err("over the cap"),
            refusal
        );
        println!(
            "wider wall, one copy at another speed: admission's reading {}, geometric reading {}, change of parameter {}, refit of the new copy {}, refit of the copy as imported {}, the edge's other use {others}, least allowance of the edge {}",
            first.work,
            geometric.work,
            change.work,
            of_new.work,
            of_imported.work,
            others + total
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

    /// The cap of section 4.2 step 8 is one allowance for the topological
    /// edge, shared by its uses: "A hard per-edge work cap refuses by name."
    ///
    /// The stretched circle of the quarter cylinder is an edge of two faces,
    /// the wall and a cap, so it has two uses. Each counts work on the one
    /// allowance of the edge. An allowance that covers the costlier use by
    /// itself, and not the two together, refuses with
    /// `numerical.gap-work-cap`. So does one unit less than the two
    /// together. Exactly the two together admits, with the answers of the
    /// full cap.
    #[test]
    fn the_uses_of_an_edge_share_its_one_work_cap() {
        let source = quarter();
        let (stretched, _) = stretched_copy(&source, "0.001");
        let (slot, _) = use_of(&source, &stretched);
        let at = edge_of(&source, &stretched);
        let uses = uses_of(&source, at);
        assert_eq!(uses.len(), 2, "an edge between two faces has two uses");
        assert!(uses.contains(&slot));

        let (full, counted) = admit(&source, &stretched, at, gap::EDGE_WORK);
        let full = full.expect("the full cap");
        let by_use: Vec<u64> = uses
            .iter()
            .map(|use_| counted.spent_by_use(*use_))
            .collect();
        let whole_edge = counted.spent_by_edge(at);
        assert_eq!(whole_edge, by_use.iter().sum::<u64>());
        assert!(by_use.iter().all(|work| *work > 0));
        let costliest = *by_use.iter().max().expect("two uses");
        assert!(costliest < whole_edge);

        for work in [costliest, whole_edge - 1] {
            refuses_by_name(
                admit(&source, &stretched, at, work).0,
                &format!("{work} of {whole_edge} units"),
            );
        }
        let (same, counted) = admit(&source, &stretched, at, whole_edge);
        let same = same.expect("the two uses together");
        assert_eq!(counted.spent_by_edge(at), whole_edge);
        assert_eq!(same.import_gaps(), full.import_gaps());
        assert_eq!(same.tolerances(), full.tolerances());
        println!(
            "quarter cylinder, one edge of two faces: its uses count {by_use:?}, the edge {whole_edge}"
        );
    }

    /// The same in an own build, which has no import rule: every admission
    /// gives the uses of one edge one allowance. A circle of the quarter
    /// cylinder is read on its cap and on the wall.
    #[test]
    fn the_uses_of_an_edge_share_its_one_work_cap_in_an_own_build() {
        let source = quarter();
        let at = source
            .edges
            .iter()
            .position(|edge| matches!(edge.curve, Curve::Circle(_)))
            .expect("a circle of the wall");
        let uses = uses_of(&source, at);
        assert_eq!(uses.len(), 2, "an edge between two faces has two uses");
        let build = |work: u64| {
            let mut accounts = accounts(&source, at, work);
            let built = admit_on(
                source.clone(),
                Vec::new(),
                NumericalOptions::default(),
                &mut accounts,
                &mut Decider::default(),
            );
            (built, accounts)
        };
        let (full, counted) = build(gap::EDGE_WORK);
        let full = full.expect("the full cap");
        let by_use: Vec<u64> = uses
            .iter()
            .map(|use_| counted.spent_by_use(*use_))
            .collect();
        let whole_edge = counted.spent_by_edge(at);
        assert_eq!(whole_edge, by_use.iter().sum::<u64>());
        assert!(by_use.iter().all(|work| *work > 0));
        let costliest = *by_use.iter().max().expect("two uses");

        for work in [costliest, whole_edge - 1] {
            refuses_by_name(build(work).0, &format!("{work} of {whole_edge} units"));
        }
        let same = build(whole_edge).0.expect("the two uses together");
        assert_eq!(same.evidence().checks(), full.evidence().checks());
        assert_eq!(same.tolerances(), full.tolerances());
        // The public build is this admission on accounts of its own.
        let public = build_with_native_trims(
            source.clone(),
            Vec::new(),
            NumericalOptions::default(),
            &mut Decider::default(),
        )
        .expect("the quarter cylinder");
        assert_eq!(public.evidence().checks(), full.evidence().checks());
    }
}
