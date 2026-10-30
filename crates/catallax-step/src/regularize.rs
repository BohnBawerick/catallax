//! Export-only repairs of planar topology a STEP reader cannot resolve. The
//! exact model is never edited. A straight, parallel wire gap lost by decimal
//! rounding is reconnected: a contained overlapping segment joins the two
//! wires, and its two parallel neighbouring faces become one. A step narrower
//! than the shortest readable edge is removed by the `ribbon` module.

use std::collections::BTreeMap;

use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, Rational, Refusal};
use catallax_shapes::{Point3, Surface};
use catallax_topology::{
    BodyId, Coedge, Edge, EdgeId, Face, FaceId, Loop, Model, Sense, Shell, ShellId, VertexId,
};

use crate::{StepOptions, StepRegularization, number};

mod ribbon;

type P = [Rational; 3];
type Wires = Vec<Vec<Coedge>>;

pub(crate) struct Prepared {
    pub(crate) model: Model,
    pub(crate) bodies: Vec<BodyId>,
    pub(crate) records: Vec<StepRegularization>,
    pub(crate) face_errors: BTreeMap<FaceId, Rational>,
    pub(crate) edge_errors: BTreeMap<EdgeId, Rational>,
}

struct Contact {
    face: FaceId,
    a: (usize, usize),
    b: (usize, usize),
    neighbours: [FaceId; 2],
    bound: Rational,
}

enum Repair {
    Gap(Contact),
    Ribbon(ribbon::Ribbon),
}

fn sub(a: &P, b: &P) -> P {
    std::array::from_fn(|i| &a[i] - &b[i])
}
fn dot(a: &P, b: &P) -> Rational {
    (0..3).fold(Rational::zero(), |s, i| s + &a[i] * &b[i])
}
fn cross(a: &P, b: &P) -> P {
    std::array::from_fn(|i| &a[(i + 1) % 3] * &b[(i + 2) % 3] - &a[(i + 2) % 3] * &b[(i + 1) % 3])
}
fn zero(p: &P) -> bool {
    p.iter().all(Rational::is_zero)
}
fn point(p: &Point3) -> Option<P> {
    Some([p.x.as_rational()?, p.y.as_rational()?, p.z.as_rational()?])
}
fn rounded(p: &P, digits: u32) -> P {
    p.each_ref().map(|v| number::rational(v, digits).value)
}
fn flip(s: Sense) -> Sense {
    if s == Sense::Same {
        Sense::Reversed
    } else {
        Sense::Same
    }
}
fn ends(model: &Model, c: &Coedge) -> Option<[VertexId; 2]> {
    let e = model.edge(c.edge);
    let pair = [e.start?, e.end?];
    Some(if c.sense == Sense::Same {
        pair
    } else {
        [pair[1], pair[0]]
    })
}
fn positions(model: &Model, c: &Coedge) -> Option<[P; 2]> {
    let [a, b] = ends(model, c)?;
    Some([
        point(&model.geom().position(model.vertex(a).point))?,
        point(&model.geom().position(model.vertex(b).point))?,
    ])
}
fn coedges(model: &Model, face: FaceId) -> Option<Wires> {
    model
        .face(face)
        .loops
        .iter()
        .map(|l| match model.loop_(*l) {
            Loop::Coedges(cs) => Some(cs.iter().map(|c| model.coedge(*c).clone()).collect()),
            Loop::Vertex(_) => None,
        })
        .collect()
}

// Only an actual overlap of the written line segments qualifies. Nearness alone
// never changes a file: separated, valid STEP wires keep their original bytes.
fn overlap(a: &[P; 2], b: &[P; 2]) -> bool {
    let d = sub(&a[1], &a[0]);
    if zero(&d) || !zero(&cross(&d, &sub(&b[0], &a[0]))) || !zero(&cross(&d, &sub(&b[1], &a[0]))) {
        return false;
    }
    let length = dot(&d, &d);
    let high = dot(&sub(&b[0], &a[0]), &d);
    let low = dot(&sub(&b[1], &a[0]), &d);
    Rational::zero() < low && low < high && high < length
}

fn neighbour(model: &Model, shell: ShellId, face: FaceId, edge: EdgeId) -> Option<FaceId> {
    let mut found = None;
    for &f in &model.shell(shell).faces {
        for l in &model.face(f).loops {
            if let Loop::Coedges(cs) = model.loop_(*l) {
                for c in cs {
                    if model.coedge(*c).edge == edge && f != face {
                        if found.is_some() {
                            return None;
                        }
                        found = Some(f);
                    }
                }
            }
        }
    }
    found
}

fn certify(
    model: &Model,
    shell: ShellId,
    face: FaceId,
    edges: [&Coedge; 2],
    digits: u32,
    limit: &Rational,
    dec: &mut Decider,
) -> Result<Option<([FaceId; 2], Rational)>, Refusal> {
    let [a, b] = edges;
    let Some([pa, qa]) = positions(model, a) else {
        return Ok(None);
    };
    let Some([pb, qb]) = positions(model, b) else {
        return Ok(None);
    };
    let da = sub(&qa, &pa);
    let db = sub(&qb, &pb);
    if !zero(&cross(&da, &db)) {
        return Ok(None);
    }
    let (Some(fa), Some(fb)) = (
        neighbour(model, shell, face, a.edge),
        neighbour(model, shell, face, b.edge),
    ) else {
        return Ok(None);
    };
    if fa == fb {
        return Ok(None);
    }
    let (Surface::Plane(sa), Surface::Plane(sb)) = (
        model.geom().surface(model.face(fa).surface),
        model.geom().surface(model.face(fb).surface),
    ) else {
        return Ok(None);
    };
    let (Some(na), Some(nb)) = (
        point(&Point3::new(
            sa.normal.x.clone(),
            sa.normal.y.clone(),
            sa.normal.z.clone(),
        )),
        point(&Point3::new(
            sb.normal.x.clone(),
            sb.normal.y.clone(),
            sb.normal.z.clone(),
        )),
    ) else {
        return Ok(None);
    };
    if !zero(&cross(&na, &nb))
        || dot(&na, &nb).is_zero()
        || (dot(&na, &nb) > Rational::zero()) != (model.face(fa).sense == model.face(fb).sense)
    {
        return Ok(None);
    }
    // The entire ribbon is a constant translation between two parallel lines.
    // Its L1 width bounds Euclidean displacement, including the neighbouring
    // parallel face's change of support. No sampled distance decides this.
    let axis = (0..3).max_by_key(|&i| da[i].abs()).unwrap();
    let t = &(&pb[axis] - &pa[axis]) / &da[axis];
    let gap: P = std::array::from_fn(|i| &pb[i] - &(&pa[i] + &(&t * &da[i])));
    let bound = gap.iter().fold(Rational::zero(), |s, x| s + x.abs());
    if bound > *limit {
        return Ok(None);
    }
    let oa = point(&sa.origin);
    let ob = point(&sb.origin);
    let (Some(oa), Some(ob)) = (oa, ob) else {
        return Ok(None);
    };
    if !dot(&sub(&pa, &oa), &na).is_zero() || !dot(&sub(&pb, &ob), &nb).is_zero() {
        return Ok(None);
    }
    // Endpoints alone do not certify the written carriers. Check the actual
    // rounded LINE placements/directions and the two written PLANE supports.
    let mut line = |c: &Coedge| -> Result<Option<(P, P)>, Refusal> {
        let CurveForm::Line { origin, dir } =
            curve_form(model.geom(), model.edge(c.edge).curve, dec)?
        else {
            return Ok(None);
        };
        let origin = number::point(&origin, digits, dec)?.xyz.map(|n| n.value);
        let dir = number::direction(&dir, digits, dec)?.xyz.map(|n| n.value);
        Ok(Some((origin, dir)))
    };
    let (Some(la), Some(lb)) = (line(a)?, line(b)?) else {
        return Ok(None);
    };
    if !zero(&cross(&la.1, &lb.1)) || !zero(&cross(&la.1, &sub(&lb.0, &la.0))) {
        return Ok(None);
    }
    let wa = number::direction(&sa.normal, digits, dec)?
        .xyz
        .map(|n| n.value);
    let wb = number::direction(&sb.normal, digits, dec)?
        .xyz
        .map(|n| n.value);
    if !zero(&cross(&wa, &wb))
        || !dot(&sub(&rounded(&ob, digits), &rounded(&oa, digits)), &wa).is_zero()
    {
        return Ok(None);
    }
    // Both adjacent face wires must be ordinary coedge wires, with exactly one
    // use of the contact edge. This excludes singular/nonmanifold surgery.
    if [face, fa, fb].iter().any(|f| {
        model
            .face(*f)
            .loops
            .iter()
            .any(|l| !matches!(model.loop_(*l), Loop::Coedges(_)))
    }) {
        return Ok(None);
    }
    Ok(Some(([fa, fb], bound)))
}

fn find(
    model: &Model,
    shell: ShellId,
    options: &StepOptions,
    dec: &mut Decider,
) -> Result<Option<Repair>, Refusal> {
    let limit = options
        .curve_tolerance
        .clone()
        .min(Rational::from_ratio(1, 1_000_000).unwrap());
    if limit < Rational::zero() {
        return Ok(None);
    }
    if let Some(contact) = find_gap(model, shell, options, &limit, dec)? {
        return Ok(Some(Repair::Gap(contact)));
    }
    Ok(ribbon::find(model, shell, options, &limit, dec)?.map(Repair::Ribbon))
}

fn find_gap(
    model: &Model,
    shell: ShellId,
    options: &StepOptions,
    limit: &Rational,
    dec: &mut Decider,
) -> Result<Option<Contact>, Refusal> {
    let digits = options.significant_digits.max(1);
    for &face in &model.shell(shell).faces {
        let f = model.face(face);
        if f.loops.len() < 2 || !matches!(model.geom().surface(f.surface), Surface::Plane(_)) {
            continue;
        }
        for i in 0..f.loops.len() {
            let Loop::Coedges(a) = model.loop_(f.loops[i]) else {
                continue;
            };
            for j in i + 1..f.loops.len() {
                let Loop::Coedges(b) = model.loop_(f.loops[j]) else {
                    continue;
                };
                for (ai, ca) in a.iter().enumerate() {
                    let ca = model.coedge(*ca);
                    let Some(pa) = positions(model, ca) else {
                        continue;
                    };
                    let ra = pa.each_ref().map(|p| rounded(p, digits));
                    for (bi, cb) in b.iter().enumerate() {
                        let cb = model.coedge(*cb);
                        let Some(pb) = positions(model, cb) else {
                            continue;
                        };
                        let rb = pb.each_ref().map(|p| rounded(p, digits));
                        for (outer, inner, x, y, ia, ib) in [
                            (ca, cb, &ra, &rb, (i, ai), (j, bi)),
                            (cb, ca, &rb, &ra, (j, bi), (i, ai)),
                        ] {
                            if overlap(x, y)
                                && let Some((neighbours, bound)) =
                                    certify(model, shell, face, [outer, inner], digits, limit, dec)?
                            {
                                return Ok(Some(Contact {
                                    face,
                                    a: ia,
                                    b: ib,
                                    neighbours,
                                    bound,
                                }));
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(None)
}

fn store_face(model: &mut Model, old: FaceId, wires: Wires) -> FaceId {
    let f = model.face(old).clone();
    let loops = wires
        .into_iter()
        .map(|wire| {
            let cs = wire.into_iter().map(|c| model.add_coedge(c)).collect();
            model.add_loop(Loop::Coedges(cs))
        })
        .collect();
    model.add_face(Face {
        surface: f.surface,
        sense: f.sense,
        loops,
    })
}
fn after(wire: &[Coedge], edge: EdgeId) -> Vec<Coedge> {
    let at = wire.iter().position(|c| c.edge == edge).unwrap();
    wire[at + 1..].iter().chain(&wire[..at]).cloned().collect()
}

fn apply(prepared: &mut Prepared, shell: ShellId, contact: Contact) -> ShellId {
    let model = &mut prepared.model;
    let mut wires = coedges(model, contact.face).unwrap();
    let a = wires[contact.a.0][contact.a.1].clone();
    let b = wires[contact.b.0][contact.b.1].clone();
    let [start, end] = ends(model, &a).unwrap();
    let [high, low] = ends(model, &b).unwrap();
    let source = model.edge(a.edge).clone();
    let sense = if a.sense == Sense::Same {
        source.sense
    } else {
        flip(source.sense)
    };
    let [left, right] = [(start, low), (high, end)].map(|(start, end)| {
        model.add_edge(Edge {
            curve: source.curve,
            sense,
            start: Some(start),
            end: Some(end),
        })
    });
    let mut joined = vec![Coedge {
        edge: left,
        sense: Sense::Same,
    }];
    joined.extend(after(&wires[contact.b.0], b.edge));
    joined.push(Coedge {
        edge: right,
        sense: Sense::Same,
    });
    wires[contact.a.0].splice(contact.a.1..=contact.a.1, joined);
    wires.remove(contact.b.0);
    let h = store_face(model, contact.face, wires);

    let [fa, fb] = contact.neighbours;
    let mut wa = coedges(model, fa).unwrap();
    let mut wb = coedges(model, fb).unwrap();
    let at_b = wb
        .iter()
        .position(|w| w.iter().any(|c| c.edge == b.edge))
        .unwrap();
    let path = after(&wb.remove(at_b), b.edge);
    for wire in &mut wa {
        if let Some(i) = wire.iter().position(|c| c.edge == a.edge) {
            let mut replacement = vec![Coedge {
                edge: right,
                sense: Sense::Reversed,
            }];
            replacement.extend(path.iter().cloned());
            replacement.push(Coedge {
                edge: left,
                sense: Sense::Reversed,
            });
            wire.splice(i..=i, replacement);
            break;
        }
    }
    wa.extend(wb);
    let merged = store_face(model, fa, wa);
    let old_a = prepared.face_errors.get(&fa).cloned().unwrap_or_default();
    let old_b = prepared.face_errors.get(&fb).cloned().unwrap_or_default();
    prepared
        .face_errors
        .insert(merged, old_a.max(&old_b + &contact.bound));
    if let Some(error) = prepared.face_errors.get(&contact.face).cloned() {
        prepared.face_errors.insert(h, error);
    }
    let edge_error = prepared
        .edge_errors
        .get(&a.edge)
        .cloned()
        .unwrap_or_default()
        + contact.bound.clone();
    for edge in [left, right] {
        prepared.edge_errors.insert(edge, edge_error.clone());
    }
    prepared.records.push(StepRegularization {
        face: contact.face,
        edges: [a.edge, b.edge],
        merged_faces: [fa, fb],
        max_deviation: contact.bound,
        removed_faces: 1,
        removed_edges: 0,
    });
    let faces = model
        .shell(shell)
        .faces
        .iter()
        .filter_map(|f| {
            if *f == fb {
                None
            } else {
                Some(if *f == contact.face {
                    h
                } else if *f == fa {
                    merged
                } else {
                    *f
                })
            }
        })
        .collect();
    model.add_shell(Shell { faces })
}

pub(crate) fn prepare(
    model: &Model,
    bodies: &[BodyId],
    options: &StepOptions,
    dec: &mut Decider,
) -> Result<Option<Prepared>, Refusal> {
    let mut first = None;
    'search: for &body in bodies {
        for &lump in &model.body(body).lumps {
            for &shell in &model.lump(lump).shells {
                if let Some(repair) = find(model, shell, options, dec)? {
                    first = Some((shell, repair));
                    break 'search;
                }
            }
        }
    }
    let Some((first_shell, first_repair)) = first else {
        return Ok(None);
    };
    let mut first = Some(first_repair);
    let mut prepared = Prepared {
        model: model.clone(),
        bodies: Vec::with_capacity(bodies.len()),
        records: Vec::new(),
        face_errors: BTreeMap::new(),
        edge_errors: BTreeMap::new(),
    };
    for &body in bodies {
        let mut body = model.body(body).clone();
        for lump in &mut body.lumps {
            let mut value = model.lump(*lump).clone();
            for shell in &mut value.shells {
                let mut repair = if *shell == first_shell && first.is_some() {
                    first.take()
                } else {
                    find(&prepared.model, *shell, options, dec)?
                };
                while let Some(found) = repair {
                    *shell = match found {
                        Repair::Gap(contact) => apply(&mut prepared, *shell, contact),
                        Repair::Ribbon(ribbon) => {
                            ribbon::apply(&mut prepared, *shell, ribbon, dec)?
                        }
                    };
                    repair = find(&prepared.model, *shell, options, dec)?;
                }
            }
            *lump = prepared.model.add_lump(value);
        }
        prepared.bodies.push(prepared.model.add_body(body));
    }
    Ok(Some(prepared))
}
