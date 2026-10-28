//! Applying the repair policy: automatic closure, accepted proposals and declines.
//!
//! A repair moves whole vertices of one side of a gap onto their partners and
//! refits each planar face through its moved vertices. Both sides are tried and
//! the one that keeps every face exactly planar with the smaller movement wins.
//! The result is rebuilt and checked before anything is returned; a gap that no
//! side can close without changing a surface stays unrepaired.

use std::collections::{BTreeMap, BTreeSet};

use catallax_measure::{Accuracy, volume};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::StepId;
use catallax_topology::{BodyId, Model};

use crate::exact::{cross, dot, is_zero, sub};
use crate::gap::{Gap, ImportOptions, decimal, measure, proposals, sqrt_enclosure};
use crate::{
    Approval, Bounds, Change, Curve, Decision, DerivativeRef, Disposition, GapKind, ImportRefusal,
    Placement, REPAIR_FORMAT, RepairRecord, Solid, SourceRef, StepFile, Surface, Vec3, build_solid,
    read_step, sha256_hex,
};

/// What an import produced.
#[derive(Clone, Debug)]
pub struct Imported {
    /// Solid entity name and body, for every solid that is closed after the decisions.
    pub bodies: Vec<(String, BodyId)>,
    /// Solids still open because a gap was declined. They are not built.
    pub left_open: Vec<String>,
    /// The change record to persist beside the recipe.
    pub record: RepairRecord,
    /// Canonical bytes of the repaired derivative, when any gap was repaired.
    pub derivative: Option<Vec<u8>>,
}

fn invalid(id: &str, detail: impl Into<String>) -> ImportRefusal {
    ImportRefusal::InvalidDecision {
        id: id.to_string(),
        detail: detail.into(),
    }
}

fn failed(id: &str, detail: impl Into<String>) -> ImportRefusal {
    ImportRefusal::ValidationFailed {
        id: id.to_string(),
        detail: detail.into(),
    }
}

/// Reads the source, applies automatic repairs and the application's decisions,
/// and builds the closed solids into `model`.
///
/// Every `ask` proposal needs a decision. Only a decline is accepted for an
/// automatic proposal, and an opening can only be declined: closing one is never
/// done. The source bytes are never changed.
pub fn apply(
    model: &mut Model,
    source: &[u8],
    name: &str,
    options: &ImportOptions,
    decisions: &[Decision],
    step: Option<&StepId>,
) -> Result<Imported, ImportRefusal> {
    let file = read_step(source)?;
    let sha = sha256_hex(source);
    crate::check_self_intersection(&file)?;
    let gaps = measure(&file, &options.search_radius_mm)?;
    let props = proposals(&file, &sha, &gaps);

    let mut by_id: BTreeMap<&str, &Decision> = BTreeMap::new();
    for d in decisions {
        let id = match d {
            Decision::Accept { id, .. } | Decision::Decline { id } => id.as_str(),
        };
        if !props.iter().any(|p| p.id == id) {
            return Err(invalid(id, "no such proposal"));
        }
        if by_id.insert(id, d).is_some() {
            return Err(invalid(id, "decided twice"));
        }
    }

    // What each gap asks for: Some(tolerance, approver) to repair, None for a decline.
    struct Plan {
        repair: Option<(Rational, Option<String>)>,
    }
    let mut plans = Vec::new();
    for (p, gap) in props.iter().zip(&gaps) {
        let plan = match (by_id.get(p.id.as_str()), p.disposition) {
            (Some(Decision::Decline { .. }), _) => Plan { repair: None },
            (None, Disposition::Auto) => Plan {
                repair: Some((crate::gap::auto_cap(), None)),
            },
            (None, Disposition::Ask) => return Err(invalid(&p.id, "needs a decision")),
            (Some(Decision::Accept { .. }), Disposition::Auto) => {
                return Err(invalid(
                    &p.id,
                    "automatic repairs take no acceptance; only a decline",
                ));
            }
            (
                Some(Decision::Accept {
                    tolerance_mm,
                    approver,
                    ..
                }),
                Disposition::Ask,
            ) => {
                if gap.kind != GapKind::SeamClose {
                    return Err(invalid(
                        &p.id,
                        "closing an opening is not supported; decline",
                    ));
                }
                let tol = Rational::parse(tolerance_mm)
                    .map_err(|_| invalid(&p.id, "tolerance is not a decimal"))?;
                let (_, hi) = sqrt_enclosure(&gap.d2);
                if tol < hi {
                    return Err(invalid(&p.id, "tolerance is below the measured distance"));
                }
                if approver.is_empty() {
                    return Err(invalid(&p.id, "an approver is required"));
                }
                Plan {
                    repair: Some((tol, Some(approver.clone()))),
                }
            }
        };
        plans.push(plan);
    }

    // Repair per solid on a working copy.
    let mut work = file.clone();
    let mut changes = Vec::new();
    let mut open: BTreeSet<usize> = BTreeSet::new();
    let mut repaired: BTreeSet<usize> = BTreeSet::new();
    for ((p, gap), plan) in props.iter().zip(&gaps).zip(&plans) {
        let before = bounds(&file.solids[gap.solid]);
        let (movement, after, approval, approver, tol) = match &plan.repair {
            None => {
                open.insert(gap.solid);
                (
                    Rational::zero(),
                    before.clone(),
                    Approval::Declined,
                    None,
                    Rational::parse(&p.tolerance_mm).expect("decimal"),
                )
            }
            Some((tol, approver)) => {
                let (solid, move2) =
                    best_repair(&work.solids[gap.solid], gap).map_err(|detail| {
                        failed(
                            &p.id,
                            format!("gap {} ({}): {detail}", p.id, p.entities.join(", ")),
                        )
                    })?;
                if move2 > tol * tol {
                    return Err(failed(&p.id, "measured movement exceeds the tolerance"));
                }
                work.solids[gap.solid] = solid;
                repaired.insert(gap.solid);
                let after = bounds(&work.solids[gap.solid]);
                (
                    move2,
                    after,
                    if approver.is_some() {
                        Approval::User
                    } else {
                        Approval::Automatic
                    },
                    approver.clone(),
                    tol.clone(),
                )
            }
        };
        let (_, hi) = sqrt_enclosure(&movement);
        changes.push(Change {
            id: p.id.clone(),
            entities: p.entities.clone(),
            kind: p.kind,
            bounds_before: before,
            bounds_after: after,
            max_movement_mm: decimal(&hi),
            tolerance_mm: decimal(&tol),
            reason: match approval {
                Approval::Declined => format!("declined: {}", p.reason),
                Approval::User => format!("approved: {}", p.reason),
                Approval::Automatic => p.reason.clone(),
            },
            approval,
            approver,
        });
    }

    // Refit and validate every repaired solid, then build the closed ones.
    let mut bodies = Vec::new();
    let mut left_open = Vec::new();
    for (si, solid) in work.solids.iter().enumerate() {
        if open.contains(&si) {
            left_open.push(solid.name.clone());
            continue;
        }
        let first = props
            .iter()
            .zip(&gaps)
            .find(|(_, g)| g.solid == si)
            .map(|(p, _)| p.id.clone());
        if let Some(id) = &first {
            let mut scratch = Model::new();
            let body = build_solid(&mut scratch, solid, None).map_err(|e| match e {
                e @ ImportRefusal::SelfIntersecting { .. } => e,
                e => failed(id, format!("the repaired solid is not valid: {e:?}")),
            })?;
            let v = volume(
                &scratch,
                body,
                &Accuracy::default(),
                &mut Decider::default(),
            )
            .map_err(|_| failed(id, "the repaired volume could not be measured"))?;
            if v.value.lo.sign() != catallax_numbers::Sign::Positive {
                return Err(failed(id, "the repaired solid is inside out or empty"));
            }
        }
        bodies.push((solid.name.clone(), build_solid(model, solid, step)?));
    }

    let derivative = (!repaired.is_empty()).then(|| canonical_bytes(&work));
    let record = RepairRecord {
        format: REPAIR_FORMAT.to_string(),
        source: SourceRef {
            name: name.to_string(),
            sha256: sha,
        },
        derivative: derivative.as_ref().map(|b| DerivativeRef {
            sha256: sha256_hex(b),
        }),
        changes,
    };
    Ok(Imported {
        bodies,
        left_open,
        record,
        derivative,
    })
}

fn bounds(solid: &Solid) -> Bounds {
    let mut lo: Option<Vec3> = None;
    let mut hi: Option<Vec3> = None;
    for face in &solid.faces {
        for bound in &face.loops {
            for e in &bound.edges {
                for p in [&e.start.point, &e.end.point] {
                    let l = lo.get_or_insert_with(|| p.clone());
                    let h = hi.get_or_insert_with(|| p.clone());
                    for i in 0..3 {
                        if p[i] < l[i] {
                            l[i] = p[i].clone();
                        }
                        if p[i] > h[i] {
                            h[i] = p[i].clone();
                        }
                    }
                }
            }
        }
    }
    let text = |v: Option<Vec3>| -> [String; 3] {
        v.map(|v| v.map(|c| decimal(&c)))
            .unwrap_or_else(|| std::array::from_fn(|_| "0".to_string()))
    };
    [text(lo), text(hi)]
}

/// Tries moving each side of the gap and returns the valid repair that moves least.
///
/// Fewer changed surfaces wins first, so a seam is closed by extending faces along
/// their own surfaces where that works; smaller movement breaks ties.
fn best_repair(solid: &Solid, gap: &Gap) -> Result<(Solid, Rational), String> {
    let mut best: Option<(Solid, Rational, usize)> = None;
    let mut last = String::new();
    for mover in [0u8, 1u8] {
        match try_move(solid, gap, mover) {
            Ok((s, move2)) => {
                let changed = s
                    .faces
                    .iter()
                    .zip(&solid.faces)
                    .filter(|(a, b)| !same_surface(&a.surface, &b.surface))
                    .count();
                if best
                    .as_ref()
                    .is_none_or(|(_, m, c)| (changed, &move2) < (*c, m))
                {
                    best = Some((s, move2, changed));
                }
            }
            Err(e) => last = e,
        }
    }
    best.map(|(s, m, _)| (s, m)).ok_or(last)
}

fn same_surface(a: &Surface, b: &Surface) -> bool {
    match (a, b) {
        (Surface::Plane(a), Surface::Plane(b)) => {
            is_zero(&cross(&a.axis, &b.axis)) && dot(&sub(&a.origin, &b.origin), &b.axis).is_zero()
        }
        _ => a == b,
    }
}

fn try_move(solid: &Solid, gap: &Gap, mover: u8) -> Result<(Solid, Rational), String> {
    // Two-colour the faces joined by paired edges.
    let mut color: BTreeMap<usize, u8> = BTreeMap::new();
    let mut changed = true;
    if let Some(p) = gap.pairs.first() {
        color.insert(p.a.face, 0);
    }
    while changed {
        changed = false;
        for p in &gap.pairs {
            for (x, y) in [(p.a.face, p.b.face), (p.b.face, p.a.face)] {
                if let Some(&c) = color.get(&x) {
                    match color.get(&y) {
                        None => {
                            color.insert(y, 1 - c);
                            changed = true;
                        }
                        Some(&d) if d == c => {
                            return Err(
                                "the faces around the gap cannot be split into two sides".into()
                            );
                        }
                        _ => {}
                    }
                }
            }
        }
    }
    let mut moves: BTreeMap<Vec3, Vec3> = BTreeMap::new();
    for p in &gap.pairs {
        for (x, y) in [(&p.a.from, &p.b.to), (&p.a.to, &p.b.from)] {
            let (from, to) = if color[&p.a.face] == mover {
                (x, y)
            } else {
                (y, x)
            };
            if let Some(prev) = moves.insert(from.clone(), to.clone())
                && prev != *to
            {
                return Err("a vertex has two different partners".into());
            }
        }
    }
    let mut out = solid.clone();
    let mut move2 = Rational::zero();
    for (from, to) in &moves {
        let d = sub(from, to);
        move2 = move2.max(dot(&d, &d));
    }
    for (fi, face) in out.faces.iter_mut().enumerate() {
        let touched = face
            .loops
            .iter()
            .flat_map(|l| &l.edges)
            .any(|e| moves.contains_key(&e.start.point) || moves.contains_key(&e.end.point));
        if let Some(&c) = color.get(&fi)
            && c != mover
            && touched
        {
            return Err("a moved vertex belongs to a face on the other side".into());
        }
        if !touched {
            continue;
        }
        for bound in &mut face.loops {
            for e in &mut bound.edges {
                let moved = moves.contains_key(&e.start.point) || moves.contains_key(&e.end.point);
                for v in [&mut e.start, &mut e.end] {
                    if let Some(to) = moves.get(&v.point) {
                        v.point = to.clone();
                    }
                }
                if moved {
                    let Curve::Line(origin, dir) = &mut e.curve else {
                        return Err(format!(
                            "curved edge {} would have to move; curved gaps are outside the first path",
                            e.edge
                        ));
                    };
                    let d = sub(&e.end.point, &e.start.point);
                    if is_zero(&d) {
                        return Err("an edge would collapse".into());
                    }
                    *origin = e.start.point.clone();
                    *dir = if e.curve_same_sense { d } else { d.map(|c| -c) };
                }
            }
        }
        refit(face)?;
    }
    // Edges are shared by name across faces; keep every occurrence consistent.
    let mut edges: BTreeMap<String, (Vec3, Vec3)> = BTreeMap::new();
    for face in &out.faces {
        for e in face.loops.iter().flat_map(|l| &l.edges) {
            let ends = (e.start.point.clone(), e.end.point.clone());
            if let Some(prev) = edges.insert(e.edge.clone(), ends.clone())
                && prev != ends
            {
                return Err("an edge shared by faces would end in two places".into());
            }
        }
    }
    Ok((out, move2))
}

/// Refits a planar face through its current vertices, which must stay coplanar.
fn refit(face: &mut crate::Face) -> Result<(), String> {
    let Surface::Plane(old) = &face.surface else {
        return Err(format!(
            "curved face {} would have to move; curved gaps are outside the first path",
            face.name
        ));
    };
    let points: Vec<Vec3> = face
        .loops
        .iter()
        .flat_map(|l| &l.edges)
        .flat_map(|e| [e.start.point.clone(), e.end.point.clone()])
        .collect();
    let p0 = points[0].clone();
    let Some(p1) = points.iter().find(|p| **p != p0).cloned() else {
        return Err("a face would collapse".into());
    };
    let Some(normal) = points
        .iter()
        .map(|p| cross(&sub(&p1, &p0), &sub(p, &p0)))
        .find(|n| !is_zero(n))
    else {
        return Err("a face would collapse to a line".into());
    };
    if points.iter().any(|p| !dot(&sub(p, &p0), &normal).is_zero()) {
        return Err("the moved vertices of a face are no longer coplanar".into());
    }
    let sign = dot(&normal, &old.axis);
    let axis = if sign.sign() == catallax_numbers::Sign::Negative {
        normal.map(|c| -c)
    } else if sign.is_zero() {
        return Err("a face would turn through a right angle".into());
    } else {
        normal
    };
    let ref_dir = if dot(&old.ref_dir, &axis).is_zero() {
        old.ref_dir.clone()
    } else {
        sub(&p1, &p0)
    };
    face.surface = Surface::Plane(Placement {
        origin: p0,
        axis,
        ref_dir,
    });
    Ok(())
}

fn v3(v: &Vec3) -> String {
    format!("{} {} {}", v[0], v[1], v[2])
}

fn placement(p: &Placement) -> String {
    format!("{}|{}|{}", v3(&p.origin), v3(&p.axis), v3(&p.ref_dir))
}

/// Canonical text of recognised data: the bytes the derivative hash covers.
pub fn canonical_bytes(file: &StepFile) -> Vec<u8> {
    let mut out = String::from("catallax-import-derivative/1\n");
    out.push_str(&format!("unit {}\n", file.mm_per_unit));
    for solid in &file.solids {
        out.push_str(&format!("solid {}\n", solid.name));
        for face in &solid.faces {
            let surface = match &face.surface {
                Surface::Plane(p) => format!("plane {}", placement(p)),
                Surface::Cylinder(p, r) => format!("cylinder {} {r}", placement(p)),
                Surface::Cone45(p, r) => format!("cone45 {} {r}", placement(p)),
                Surface::Sphere(p, r) => format!("sphere {} {r}", placement(p)),
                Surface::Torus(p, a, b) => format!("torus {} {a} {b}", placement(p)),
            };
            out.push_str(&format!(
                "face {} {} {surface}\n",
                face.name, face.same_sense
            ));
            for bound in &face.loops {
                out.push_str(&format!("bound {} {}\n", bound.outer, bound.orientation));
                for e in &bound.edges {
                    let curve = match &e.curve {
                        Curve::Line(o, d) => format!("line {}|{}", v3(o), v3(d)),
                        Curve::Circle(p, r) => format!("circle {} {r}", placement(p)),
                        Curve::Ellipse(p, a, b) => format!("ellipse {} {a} {b}", placement(p)),
                    };
                    out.push_str(&format!(
                        "edge {} {} {} {}|{} {}|{} {curve}\n",
                        e.edge,
                        e.forward,
                        e.curve_same_sense,
                        e.start.name,
                        v3(&e.start.point),
                        e.end.name,
                        v3(&e.end.point)
                    ));
                }
            }
        }
    }
    out.into_bytes()
}

/// Rebuilds an import from the source bytes and a persisted record.
///
/// The record's approvals become the decisions, the repair runs again, and the
/// result must reproduce the record exactly, derivative hash included. Bytes that
/// are not the recorded source are a hash mismatch, as is a derivative that
/// differs from the recorded one.
pub fn rebuild(
    model: &mut Model,
    source: &[u8],
    options: &ImportOptions,
    record: &RepairRecord,
    step: Option<&StepId>,
) -> Result<Imported, ImportRefusal> {
    if record.format != REPAIR_FORMAT {
        return Err(ImportRefusal::Malformed {
            detail: format!("unknown repair record format `{}`", record.format),
        });
    }
    let actual = sha256_hex(source);
    if actual != record.source.sha256 {
        return Err(ImportRefusal::HashMismatch {
            expected: record.source.sha256.clone(),
            actual,
        });
    }
    let decisions: Vec<Decision> = record
        .changes
        .iter()
        .filter_map(|c| match c.approval {
            Approval::Automatic => None,
            Approval::Declined => Some(Decision::Decline { id: c.id.clone() }),
            Approval::User => Some(Decision::Accept {
                id: c.id.clone(),
                tolerance_mm: c.tolerance_mm.clone(),
                approver: c.approver.clone().unwrap_or_default(),
            }),
        })
        .collect();
    let imported = apply(
        model,
        source,
        &record.source.name,
        options,
        &decisions,
        step,
    )?;
    let hash = |r: &RepairRecord| r.derivative.as_ref().map(|d| d.sha256.clone());
    if hash(&imported.record) != hash(record) {
        return Err(ImportRefusal::HashMismatch {
            expected: hash(record).unwrap_or_else(|| "none".into()),
            actual: hash(&imported.record).unwrap_or_else(|| "none".into()),
        });
    }
    if let Some(change) = record
        .changes
        .iter()
        .zip(&imported.record.changes)
        .find(|(a, b)| a != b)
    {
        return Err(failed(&change.0.id, "the record does not match the repair"));
    }
    if record.changes.len() != imported.record.changes.len() {
        return Err(failed("", "the record lists a different set of gaps"));
    }
    Ok(imported)
}
