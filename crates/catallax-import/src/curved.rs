//! Self-intersection of faces with a curved surface or a circular edge.
//!
//! Each such pair is decided with the certified intersection engine, never with
//! floating point. The carriers of the two faces meet in curves; each curve is cut
//! where it crosses either face's boundary, and one exact point inside every piece
//! is located in both closed faces. The faces intersect when a piece lies in both
//! and in the interior of at least one, the same rule as for planar faces. Where
//! the engine cannot decide, the pair is reported as not verified, never as clear.

use catallax_boolean::face_contains;
use catallax_intersect::{
    Component, CurveSurface, Cyclic, SurfacePair, arc_point, curve_is_closed, cyclic_order,
    intersect_curve_surface, intersect_surfaces, on_curve, order_on_curve, same_curve,
};
use catallax_numbers::{Decider, Outcome, Rational, Refusal, StrictOrder};
use catallax_shapes::{CurveId, PointDef, PointId};
use catallax_topology::{EdgeId, FaceId, Loop, Model, Sense};

use crate::exact::{cross, dot, is_zero, sub, unit};
use crate::{Curve, Face, Solid, Surface, Vec3};

/// Why a pair could not be decided.
pub(crate) type NotVerified = String;

fn why(r: Refusal) -> NotVerified {
    format!("{r:?}")
}

fn coincident<T>(o: Outcome<T>) -> Result<bool, NotVerified> {
    match o {
        Outcome::Proven(_) => Ok(false),
        Outcome::ProvenCoincident(_) => Ok(true),
        Outcome::NotExact(n) => Err(format!("{}: {}", n.code, n.detail)),
        Outcome::Refused(r) => Err(why(r)),
    }
}

fn face_edges(model: &Model, face: FaceId) -> Vec<EdgeId> {
    let mut out = Vec::new();
    for &l in &model.face(face).loops {
        if let Loop::Coedges(cs) = model.loop_(l) {
            for &c in cs {
                let e = model.coedge(c).edge;
                if !out.contains(&e) {
                    out.push(e);
                }
            }
        }
    }
    out
}

/// The edge's ends along its curve's natural direction; both absent for a ring.
fn arc_ends(model: &Model, edge: EdgeId) -> (Option<PointId>, Option<PointId>) {
    let e = model.edge(edge);
    let a = e.start.map(|v| model.vertex(v).point);
    let b = e.end.map(|v| model.vertex(v).point);
    if e.sense == Sense::Same {
        (a, b)
    } else {
        (b, a)
    }
}

fn on_arc(
    model: &Model,
    curve: CurveId,
    ends: (Option<PointId>, Option<PointId>),
    point: PointId,
    dec: &mut Decider,
) -> Result<bool, NotVerified> {
    if !coincident(on_curve(model.geom(), curve, point, dec))? {
        return Ok(false);
    }
    let (Some(a), Some(b)) = ends else {
        return Ok(true);
    };
    let geom = model.geom();
    if coincident(geom.same_point(a, point, dec))?
        || coincident(geom.same_point(b, point, dec))?
        || coincident(geom.same_point(a, b, dec))?
    {
        return Ok(true);
    }
    let closed = curve_is_closed(geom, curve, dec).map_err(why)?;
    if !closed {
        let after = matches!(
            order_on_curve(geom, curve, a, point, dec),
            Outcome::Proven(StrictOrder::Less)
        );
        let before = matches!(
            order_on_curve(geom, curve, point, b, dec),
            Outcome::Proven(StrictOrder::Less)
        );
        Ok(after && before)
    } else {
        match cyclic_order(geom, curve, [a, point, b], dec) {
            Outcome::Proven(Cyclic::Forward) | Outcome::ProvenCoincident(_) => Ok(true),
            Outcome::Proven(Cyclic::Backward) => Ok(false),
            Outcome::NotExact(n) => Err(format!("{}: {}", n.code, n.detail)),
            Outcome::Refused(r) => Err(why(r)),
        }
    }
}

fn sort_points(
    model: &Model,
    curve: CurveId,
    points: &mut [PointId],
    dec: &mut Decider,
) -> Result<(), NotVerified> {
    if points.len() < 2 {
        return Ok(());
    }
    let open = !curve_is_closed(model.geom(), curve, dec).map_err(why)?;
    let anchor = points[0];
    for i in 1..points.len() {
        let mut j = i;
        while j > 0 {
            let before = if open {
                match order_on_curve(model.geom(), curve, points[j], points[j - 1], dec) {
                    Outcome::Proven(StrictOrder::Less) => true,
                    Outcome::Proven(StrictOrder::Greater) => false,
                    Outcome::ProvenCoincident(_) => {
                        return Err("breakpoints on a curve were not distinct".into());
                    }
                    Outcome::NotExact(n) => return Err(format!("{}: {}", n.code, n.detail)),
                    Outcome::Refused(r) => return Err(why(r)),
                }
            } else if j == 1 {
                false
            } else {
                match cyclic_order(model.geom(), curve, [anchor, points[j], points[j - 1]], dec) {
                    Outcome::Proven(Cyclic::Forward) => true,
                    Outcome::Proven(Cyclic::Backward) => false,
                    Outcome::ProvenCoincident(_) => {
                        return Err("breakpoints on a curve were not distinct".into());
                    }
                    Outcome::NotExact(n) => return Err(format!("{}: {}", n.code, n.detail)),
                    Outcome::Refused(r) => return Err(why(r)),
                }
            };
            if !before {
                break;
            }
            points.swap(j, j - 1);
            j -= 1;
        }
    }
    Ok(())
}

fn exact<T>(answer: catallax_numbers::Answer<T>, what: &str) -> Result<T, NotVerified> {
    if answer.is_exact() {
        Ok(answer.value)
    } else {
        Err(format!("{what} is not exact"))
    }
}

/// Whether the two faces' interiors intersect, `Err` when it cannot be decided.
pub(crate) fn pair_intersects(
    model: &mut Model,
    fa: FaceId,
    fb: FaceId,
    dec: &mut Decider,
) -> Result<bool, NotVerified> {
    let (sa, sb) = (model.face(fa).surface, model.face(fb).surface);
    for f in [fa, fb] {
        if model.face(f).loops.is_empty()
            && !matches!(
                model.geom().surface(model.face(f).surface),
                catallax_shapes::Surface::Sphere(_) | catallax_shapes::Surface::Torus(_)
            )
        {
            return Err("an unbounded face has no boundary to trim against".into());
        }
    }
    let pair = intersect_surfaces(model.geom_mut(), sa, sb, dec).map_err(why)?;
    match exact(pair, "surface intersection")? {
        SurfacePair::Disjoint => Ok(false),
        SurfacePair::Coincident { same_sense } => {
            let same_way = same_sense == (model.face(fa).sense == model.face(fb).sense);
            if same_way {
                Err("the carriers coincide and the faces look the same way; overlap of their regions is not decided".into())
            } else {
                Ok(false)
            }
        }
        SurfacePair::Meets(components) => {
            for component in components {
                if let Component::Curve { curve, .. } = component
                    && curve_pieces_intersect(model, fa, fb, curve, dec)?
                {
                    return Ok(true);
                }
            }
            Ok(false)
        }
    }
}

fn curve_pieces_intersect(
    model: &mut Model,
    fa: FaceId,
    fb: FaceId,
    k: CurveId,
    dec: &mut Decider,
) -> Result<bool, NotVerified> {
    let mut cuts: Vec<PointId> = Vec::new();
    for (face, other) in [(fa, fb), (fb, fa)] {
        let other_surface = model.face(other).surface;
        for edge in face_edges(model, face) {
            let curve = model.edge(edge).curve;
            let ends = arc_ends(model, edge);
            let hits = intersect_curve_surface(model.geom_mut(), curve, other_surface, dec)
                .map_err(why)?;
            match exact(hits, "curve-surface intersection")? {
                CurveSurface::Disjoint => {}
                CurveSurface::Contained => {
                    // The edge lies in both carriers, so it is part of some meeting
                    // curve. Separate meeting curves do not touch, so an edge of
                    // another one cuts this curve nowhere.
                    let same = coincident(same_curve(model.geom(), k, curve, dec))?;
                    let mut on_k = false;
                    for v in [ends.0, ends.1].into_iter().flatten() {
                        if coincident(on_curve(model.geom(), k, v, dec))? {
                            cuts.push(v);
                            on_k = true;
                        }
                    }
                    if !same && !on_k && ends.0.is_none() {
                        // A ring has no end to test; test a point of it.
                        let sample = exact(
                            arc_point(model.geom(), curve, None, None, dec).map_err(why)?,
                            "ring point",
                        )?;
                        let id = model.geom_mut().add_point(PointDef::Explicit(sample));
                        if coincident(on_curve(model.geom(), k, id, dec))? {
                            return Err("a ring edge lies in both carriers and meets this curve without matching it".into());
                        }
                    }
                }
                CurveSurface::Meets(list) => {
                    for hit in list {
                        if on_arc(model, curve, ends, hit.point, dec)?
                            && coincident(on_curve(model.geom(), k, hit.point, dec))?
                        {
                            cuts.push(hit.point);
                        }
                    }
                }
            }
        }
    }
    let mut distinct: Vec<PointId> = Vec::new();
    for p in cuts {
        let mut seen = false;
        for &q in &distinct {
            if coincident(model.geom().same_point(p, q, dec))? {
                seen = true;
                break;
            }
        }
        if !seen {
            distinct.push(p);
        }
    }
    sort_points(model, k, &mut distinct, dec)?;
    let closed = curve_is_closed(model.geom(), k, dec).map_err(why)?;
    let mut arcs: Vec<(Option<PointId>, Option<PointId>)> = Vec::new();
    if closed {
        match distinct.len() {
            0 => arcs.push((None, None)),
            1 => arcs.push((Some(distinct[0]), Some(distinct[0]))),
            n => {
                for i in 0..n {
                    arcs.push((Some(distinct[i]), Some(distinct[(i + 1) % n])));
                }
            }
        }
    } else {
        for w in distinct.windows(2) {
            arcs.push((Some(w[0]), Some(w[1])));
        }
    }
    for (from, to) in arcs {
        let q = exact(
            arc_point(model.geom(), k, from, to, dec).map_err(why)?,
            "arc point",
        )?;
        let (in_a, in_b) = (
            face_contains(model, fa, &q, dec).map_err(why)?,
            face_contains(model, fb, &q, dec).map_err(why)?,
        );
        if !(in_a && in_b) {
            continue;
        }
        let id = model.geom_mut().add_point(PointDef::Explicit(q));
        let mut boundary = [false; 2];
        for (slot, face) in [fa, fb].into_iter().enumerate() {
            for edge in face_edges(model, face) {
                let curve = model.edge(edge).curve;
                let ends = arc_ends(model, edge);
                if on_arc(model, curve, ends, id, dec)? {
                    boundary[slot] = true;
                    break;
                }
            }
        }
        if !(boundary[0] && boundary[1]) {
            return Ok(true);
        }
    }
    Ok(false)
}

type Bounds = [Vec3; 2];

fn hull(boxes: impl IntoIterator<Item = Bounds>) -> Option<Bounds> {
    let mut out: Option<Bounds> = None;
    for [lo, hi] in boxes {
        match &mut out {
            None => out = Some([lo, hi]),
            Some([l, h]) => {
                for i in 0..3 {
                    if lo[i] < l[i] {
                        l[i] = lo[i].clone();
                    }
                    if hi[i] > h[i] {
                        h[i] = hi[i].clone();
                    }
                }
            }
        }
    }
    out
}

fn around(c: &Vec3, r: &Rational) -> Bounds {
    [
        std::array::from_fn(|i| &c[i] - r),
        std::array::from_fn(|i| &c[i] + r),
    ]
}

/// A sound box around the face, or `None` when none is proven.
pub(crate) fn face_bounds(face: &Face) -> Option<Bounds> {
    match &face.surface {
        Surface::Sphere(p, r) => return Some(around(&p.origin, r)),
        Surface::Torus(p, major, minor) => return Some(around(&p.origin, &(major + minor))),
        _ => {}
    }
    if face.loops.is_empty() {
        return None;
    }
    let axial = match &face.surface {
        Surface::Cylinder(p, _) | Surface::Cone45(p, _) => Some(p),
        _ => None,
    };
    let mut boxes = Vec::new();
    let mut ts: Vec<Rational> = Vec::new();
    for e in face.loops.iter().flat_map(|l| &l.edges) {
        match &e.curve {
            Curve::Line(..) => {
                for p in [&e.start.point, &e.end.point] {
                    boxes.push([p.clone(), p.clone()]);
                    if let Some(a) = axial {
                        ts.push(&dot(&sub(p, &a.origin), &a.axis) / &dot(&a.axis, &a.axis));
                    }
                }
            }
            Curve::Circle(c, r) | Curve::Ellipse(c, r, _) => {
                let radius = match &e.curve {
                    Curve::Ellipse(_, a, b) => a.clone().max(b.clone()),
                    _ => r.clone(),
                };
                boxes.push(around(&c.origin, &radius));
                if let Some(a) = axial {
                    // Only a section square to the axis has one axial position.
                    if !is_zero(&cross(&c.axis, &a.axis)) {
                        return None;
                    }
                    ts.push(&dot(&sub(&c.origin, &a.origin), &a.axis) / &dot(&a.axis, &a.axis));
                }
            }
        }
    }
    let planar = matches!(face.surface, Surface::Plane(_));
    if planar {
        return hull(boxes);
    }
    let a = axial?;
    let (tmin, tmax) = (ts.iter().min()?.clone(), ts.iter().max()?.clone());
    let radius = match &face.surface {
        Surface::Cylinder(_, r) => r.clone(),
        Surface::Cone45(_, r) => {
            // Radius grows with the unit-axis distance from the placement.
            let u = unit(&a.axis)?;
            let len = &dot(&a.axis, &u) / &dot(&u, &u);
            let (r0, r1) = (r + &(&tmin * &len), r + &(&tmax * &len));
            r0.abs().max(r1.abs())
        }
        _ => return None,
    };
    let at = |t: &Rational| -> Vec3 { std::array::from_fn(|i| &a.origin[i] + &(&a.axis[i] * t)) };
    let (p0, p1) = (at(&tmin), at(&tmax));
    hull([around(&p0, &radius), around(&p1, &radius)])
}

/// Whether two sound boxes are strictly apart.
pub(crate) fn apart(a: &Bounds, b: &Bounds) -> bool {
    (0..3).any(|i| a[1][i] < b[0][i] || b[1][i] < a[0][i])
}

/// Faces of a solid, scratch-built, ready for pair queries.
pub(crate) struct Scratch {
    pub model: Model,
    pub faces: Vec<FaceId>,
}

pub(crate) fn scratch(solid: &Solid) -> Result<Scratch, NotVerified> {
    let mut model = Model::new();
    let faces = crate::build::scratch_faces(&mut model, solid)
        .map_err(|e| format!("the faces could not be set up for the intersection engine: {e:?}"))?;
    Ok(Scratch { model, faces })
}
