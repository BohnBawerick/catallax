//! Exact circular spline repair on a complete capped cylinder.

use crate::incidence_arithmetic::{
    arithmetic as op, checkpoint, compare, cross, dot, equal, maximum, minimum, same,
    same_sequence, sqrt_enclosure, sub,
};
use crate::incidence_repair::{Net, Prepared};
use crate::{Curve, Face, OrientedEdge, Placement, Solid, Surface, Vec3};
use catallax_numbers::{BoundKind, Decider, Mark, MarkLog, Rational, Sign};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

struct Frame {
    origin: Vec3,
    axes: [Vec3; 3],
    radius: Rational,
}
fn scaled(a: &Vec3, b: &Rational, dec: &mut Decider) -> Result<Vec3, String> {
    Ok([
        op(&a[0], b, '*', dec)?,
        op(&a[1], b, '*', dec)?,
        op(&a[2], b, '*', dec)?,
    ])
}
fn unit(a: &Vec3, dec: &mut Decider) -> Result<Vec3, String> {
    let mut denominator = Rational::one();
    for x in a {
        denominator = op(
            &denominator,
            &Rational::from_integer(x.denominator().clone()),
            '*',
            dec,
        )?;
    }
    let a = scaled(a, &denominator, dec)?;
    let norm = dot(&a, &a, dec)?;
    let root = sqrt_enclosure(&norm, dec)?.0;
    if root.is_zero() || !equal(&op(&root, &root, '*', dec)?, &norm, dec)? {
        return Err("curved rational-frame certificate: direction norm has no exact supported rational square root".into());
    }
    Ok([
        op(&a[0], &root, '/', dec)?,
        op(&a[1], &root, '/', dec)?,
        op(&a[2], &root, '/', dec)?,
    ])
}
fn zero(a: &Vec3) -> bool {
    a.iter().all(Rational::is_zero)
}
impl Frame {
    fn new(p: &Placement, r: &Rational, dec: &mut Decider) -> Result<Self, String> {
        if r.sign() != Sign::Positive {
            return Err(
                "curved rational-frame certificate: cylinder radius is not positive".into(),
            );
        }
        let z = unit(&p.axis, dec)?;
        let x = unit(&p.ref_dir, dec)?;
        if !dot(&x, &z, dec)?.is_zero() {
            return Err("curved rational-frame certificate: nonorthogonal cylinder frame".into());
        }
        let y = cross(&z, &x, dec)?;
        Ok(Self {
            origin: p.origin.clone(),
            axes: [x, y, z],
            radius: r.clone(),
        })
    }
    fn local(&self, p: &Vec3, dec: &mut Decider) -> Result<Vec3, String> {
        let d = sub(p, &self.origin, dec)?;
        Ok([
            dot(&d, &self.axes[0], dec)?,
            dot(&d, &self.axes[1], dec)?,
            dot(&d, &self.axes[2], dec)?,
        ])
    }
    fn centre(&self, h: &Rational, dec: &mut Decider) -> Result<Vec3, String> {
        let v = scaled(&self.axes[2], h, dec)?;
        Ok([
            op(&self.origin[0], &v[0], '+', dec)?,
            op(&self.origin[1], &v[1], '+', dec)?,
            op(&self.origin[2], &v[2], '+', dec)?,
        ])
    }
    fn on_cylinder(&self, p: &Vec3, dec: &mut Decider) -> Result<bool, String> {
        let a = self.local(p, dec)?;
        let x = op(&a[0], &a[0], '*', dec)?;
        let y = op(&a[1], &a[1], '*', dec)?;
        equal(
            &op(&x, &y, '+', dec)?,
            &op(&self.radius, &self.radius, '*', dec)?,
            dec,
        )
    }
}
fn local_nets(
    edge: &OrientedEdge,
    frame: &Frame,
    h: &Rational,
    dec: &mut Decider,
) -> Result<Vec<Net>, String> {
    let mut nets = crate::incidence_circle::source_nets(edge, dec)?;
    let centre = frame.centre(h, dec)?;
    for net in &mut nets {
        for p in net {
            let offset = scaled(&centre, &p[3], dec)?;
            let a = sub(&[p[0].clone(), p[1].clone(), p[2].clone()], &offset, dec)?;
            *p = [
                dot(&a, &frame.axes[0], dec)?,
                dot(&a, &frame.axes[1], dec)?,
                dot(&a, &frame.axes[2], dec)?,
                p[3].clone(),
            ];
        }
    }
    Ok(nets)
}
fn rotated(nets: &[Net], rotation: usize, dec: &mut Decider) -> Result<Vec<Net>, String> {
    let negative = Rational::from_i64(-1);
    let mut result = nets.to_vec();
    for net in &mut result {
        for p in net {
            checkpoint(dec)?;
            let [x, y] = [p[0].clone(), p[1].clone()];
            let [x, y] = match rotation {
                0 => [x, y],
                1 => [y, op(&x, &negative, '*', dec)?],
                2 => [op(&x, &negative, '*', dec)?, op(&y, &negative, '*', dec)?],
                _ => [op(&y, &negative, '*', dec)?, x],
            };
            p[0] = x;
            p[1] = y;
        }
    }
    Ok(result)
}
fn replacement(
    edge: &OrientedEdge,
    frame: &Frame,
    dec: &mut Decider,
) -> Result<(OrientedEdge, Vec<Net>), String> {
    if !frame.on_cylinder(&edge.start.point, dec)? || !frame.on_cylinder(&edge.end.point, dec)? {
        return Err(format!(
            "curved endpoint incidence: edge {} needs exact shared vertices on the unchanged cylinder",
            edge.edge
        ));
    }
    let a = frame.local(&edge.start.point, dec)?;
    let b = frame.local(&edge.end.point, dec)?;
    let mut target = edge.clone();
    if equal(&a[2], &b[2], dec)? {
        let nets = local_nets(edge, frame, &a[2], dec)?;
        let mut chosen = None;
        for rotation in 0..4 {
            let local = rotated(&nets, rotation, dec)?;
            let mut valid = true;
            for net in &local {
                for p in net {
                    let rw = op(&frame.radius, &p[3], '*', dec)?;
                    if op(&rw, &p[0], '+', dec)?.sign() != Sign::Positive {
                        valid = false;
                    }
                }
            }
            if valid {
                chosen = Some((rotation, local));
                break;
            }
        }
        let (rotation, local) = chosen.ok_or(
            "circle chart certificate: source control hull crosses every supported chart pole",
        )?;
        let (sign, nets) = crate::incidence_circle::image(&local, &frame.radius, dec)?;
        let reference = match rotation {
            0 => frame.axes[0].clone(),
            1 => frame.axes[1].clone(),
            2 => scaled(&frame.axes[0], &Rational::from_i64(-1), dec)?,
            _ => scaled(&frame.axes[1], &Rational::from_i64(-1), dec)?,
        };
        target.curve = Curve::Circle(
            Placement {
                origin: frame.centre(&a[2], dec)?,
                axis: frame.axes[2].clone(),
                ref_dir: reference,
            },
            frame.radius.clone(),
        );
        if sign == Sign::Negative {
            target.curve_same_sense = !target.curve_same_sense;
        }
        Ok((target, nets))
    } else {
        Err("curved target certificate: physical generator edges need a separate complete-neighbour feature proof".into())
    }
}

fn physical_edges<'a>(face: &'a Face, dec: &mut Decider) -> Result<Vec<&'a OrientedEdge>, String> {
    let mut out = Vec::new();
    for bound in &face.loops {
        for (i, e) in bound.edges.iter().enumerate() {
            checkpoint(dec)?;
            let cancelled = bound
                .edges
                .iter()
                .enumerate()
                .any(|(j, other)| j != i && e.edge == other.edge && e.forward != other.forward);
            if !cancelled {
                out.push(e);
            }
        }
    }
    Ok(out)
}
fn angle(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<std::cmp::Ordering, String> {
    let half = |a: &Vec3| {
        a[1].sign() == Sign::Positive || (a[1].is_zero() && a[0].sign() != Sign::Negative)
    };
    let (ha, hb) = (half(a), half(b));
    if ha != hb {
        return Ok(if ha {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Greater
        });
    }
    let xy = op(&a[0], &b[1], '*', dec)?;
    let yx = op(&a[1], &b[0], '*', dec)?;
    Ok(match op(&xy, &yx, '-', dec)?.sign() {
        Sign::Positive => std::cmp::Ordering::Less,
        Sign::Negative => std::cmp::Ordering::Greater,
        Sign::Zero => std::cmp::Ordering::Equal,
    })
}
fn ring(
    edges: &[&OrientedEdge],
    frame: &Frame,
    h: &Rational,
    dec: &mut Decider,
) -> Result<Rational, String> {
    if edges.is_empty() {
        return Err("curved ring certificate: empty cap boundary".into());
    }
    let mut points: Vec<Vec3> = Vec::new();
    for e in edges {
        checkpoint(dec)?;
        let Curve::Circle(p, r) = &e.curve else {
            return Err(
                "curved ring certificate: remaining physical boundary is not a native circle"
                    .into(),
            );
        };
        let centre = frame.local(&p.origin, dec)?;
        if !centre[0].is_zero()
            || !centre[1].is_zero()
            || !equal(&centre[2], h, dec)?
            || !equal(r, &frame.radius, dec)?
            || !zero(&cross(&p.axis, &frame.axes[2], dec)?)
            || zero(&p.axis)
            || zero(&p.ref_dir)
            || !dot(&p.axis, &p.ref_dir, dec)?.is_zero()
        {
            return Err(
                "curved ring certificate: circle carrier differs from its exact cylinder cap"
                    .into(),
            );
        }
        for v in [&e.start.point, &e.end.point] {
            if !frame.on_cylinder(v, dec)? {
                return Err("curved ring certificate: vertex is off the exact cylinder".into());
            }
            let local = frame.local(v, dec)?;
            if !equal(&local[2], h, dec)? {
                return Err("curved ring certificate: vertex is off the exact cap".into());
            }
            let mut found = false;
            for p in &points {
                found |= same(&local, p, dec)?;
            }
            if !found {
                points.push(local);
            }
        }
    }
    for i in 1..points.len() {
        let mut j = i;
        while j > 0 && angle(&points[j], &points[j - 1], dec)?.is_lt() {
            points.swap(j, j - 1);
            j -= 1;
        }
    }
    let mut starts = BTreeSet::new();
    for e in edges {
        let mut indices = Vec::new();
        for vertex in [&e.start.point, &e.end.point] {
            let p = frame.local(vertex, dec)?;
            let mut found = None;
            for (i, q) in points.iter().enumerate() {
                if same(&p, q, dec)? {
                    found = Some(i);
                    break;
                }
            }
            indices.push(found.ok_or("ring endpoint index")?);
        }
        let Curve::Circle(p, _) = &e.curve else {
            unreachable!()
        };
        let positive = dot(&p.axis, &frame.axes[2], dec)?.sign() == Sign::Positive;
        if e.curve_same_sense != positive {
            indices.swap(0, 1);
        }
        if (indices[0] + 1) % points.len() != indices[1] || !starts.insert(indices[0]) {
            return Err("curved ring coverage certificate: circular arcs overlap, repeat or skip an angular sector".into());
        }
    }
    if starts.len() != points.len() {
        return Err("curved ring coverage certificate: cap has uncovered angular sectors".into());
    }
    let mut feature = op(&frame.radius, &frame.radius, '*', dec)?;
    for i in 0..points.len() {
        for j in i + 1..points.len() {
            let d = sub(&points[i], &points[j], dec)?;
            feature = minimum(feature, dot(&d, &d, dec)?, dec)?;
        }
    }
    Ok(feature)
}
fn feature(
    source: &Solid,
    frame: &Frame,
    heights: &[Rational; 2],
    dec: &mut Decider,
) -> Result<Rational, String> {
    let d = op(&heights[1], &heights[0], '-', dec)?;
    let mut feature = op(&d, &d, '*', dec)?;
    let mut cap_names = Vec::new();
    let mut cylinder_names = [BTreeSet::new(), BTreeSet::new()];
    for face in &source.faces {
        let edges = physical_edges(face, dec)?;
        if matches!(face.surface, Surface::Cylinder(..)) {
            let mut rings = [Vec::new(), Vec::new()];
            for e in edges {
                let h = frame.local(&e.start.point, dec)?[2].clone();
                let i = if equal(&h, &heights[0], dec)? {
                    0
                } else if equal(&h, &heights[1], dec)? {
                    1
                } else {
                    return Err(
                        "curved feature certificate: cylinder edge is outside both cap planes"
                            .into(),
                    );
                };
                cylinder_names[i].insert(e.edge.clone());
                rings[i].push(e);
            }
            for i in 0..2 {
                feature = minimum(feature, ring(&rings[i], frame, &heights[i], dec)?, dec)?;
            }
        } else {
            if face.loops.len() != 1 || !face.loops[0].outer {
                return Err(
                    "curved feature certificate: cap must have one outer ring and no holes".into(),
                );
            }
            let Surface::Plane(p) = &face.surface else {
                unreachable!()
            };
            let h = frame.local(&p.origin, dec)?[2].clone();
            let i = usize::from(!equal(&h, &heights[0], dec)?);
            feature = minimum(feature, ring(&edges, frame, &heights[i], dec)?, dec)?;
            let names = edges
                .iter()
                .map(|e| e.edge.clone())
                .collect::<BTreeSet<_>>();
            if names.len() != edges.len() {
                return Err("curved ring identity certificate: repeated source edge".into());
            }
            cap_names.push((i, names));
        }
    }
    for (i, names) in cap_names {
        if names != cylinder_names[i] {
            return Err("curved ring identity certificate: cap and cylinder boundaries use different source edges".into());
        }
    }
    Ok(feature)
}
#[derive(Serialize)]
struct Record<'a> {
    solid: &'a str,
    entity: &'a str,
    entity_kind: &'static str,
    source_geometry_sha256: &'a str,
    definition: &'static str,
    deviation_lower_mm: String,
    deviation_upper_mm: String,
    approved_tolerance_mm: String,
    feature_size_squared_mm2: String,
}

fn same_curve(a: &Curve, b: &Curve, dec: &mut Decider) -> Result<bool, String> {
    let placement = |a: &Placement, b: &Placement, dec: &mut Decider| -> Result<bool, String> {
        Ok(same(&a.origin, &b.origin, dec)?
            && same(&a.axis, &b.axis, dec)?
            && same(&a.ref_dir, &b.ref_dir, dec)?)
    };
    Ok(match (a, b) {
        (Curve::Line(a, d), Curve::Line(b, e)) => same(a, b, dec)? && same(d, e, dec)?,
        (Curve::Circle(a, r), Curve::Circle(b, s)) => placement(a, b, dec)? && equal(r, s, dec)?,
        (Curve::BSpline(a), Curve::BSpline(b)) => {
            let (a, b) = (&a.geometry, &b.geometry);
            if a.degree() != b.degree()
                || a.controls().len() != b.controls().len()
                || !same_sequence(a.knots(), b.knots(), dec)?
                || !same_sequence(a.weights(), b.weights(), dec)?
            {
                false
            } else {
                let mut same_controls = true;
                for (a, b) in a.controls().iter().zip(b.controls()) {
                    same_controls &= same(a, b, dec)?;
                }
                same_controls
            }
        }
        _ => false,
    })
}
fn identities(source: &Solid, dec: &mut Decider) -> Result<(), String> {
    let mut vertices: BTreeMap<&str, &Vec3> = BTreeMap::new();
    let mut edges: BTreeMap<&str, &OrientedEdge> = BTreeMap::new();
    for edge in source
        .faces
        .iter()
        .flat_map(|f| &f.loops)
        .flat_map(|b| &b.edges)
    {
        checkpoint(dec)?;
        for vertex in [&edge.start, &edge.end] {
            if let Some(point) = vertices.insert(&vertex.name, &vertex.point)
                && !same(point, &vertex.point, dec)?
            {
                return Err("curved source identity certificate: one vertex name has inconsistent coordinates".into());
            }
        }
        if let Some(old) = edges.insert(&edge.edge, edge)
            && (old.start.name != edge.start.name
                || old.end.name != edge.end.name
                || old.curve_same_sense != edge.curve_same_sense
                || !same_curve(&old.curve, &edge.curve, dec)?)
        {
            return Err("curved source identity certificate: one edge name has inconsistent endpoint or curve data".into());
        }
    }
    Ok(())
}

pub(crate) fn prepare(source: &Solid) -> Result<Prepared, String> {
    let mut dec = Decider::default();
    checkpoint(&mut dec)?;
    if source.faces.len() != 3 {
        return Err("curved feature certificate: this step requires one cylinder and two perpendicular planar caps".into());
    }
    let cylinders = source
        .faces
        .iter()
        .filter_map(|f| {
            if let Surface::Cylinder(p, r) = &f.surface {
                Some((p, r))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    if cylinders.len() != 1 {
        return Err(
            "curved support certificate: this step requires exactly one unchanged cylinder".into(),
        );
    }
    identities(source, &mut dec)?;
    let frame = Frame::new(cylinders[0].0, cylinders[0].1, &mut dec)?;
    let mut heights = Vec::new();
    for face in &source.faces {
        match &face.surface {
            Surface::Cylinder(..) => {},
            Surface::Plane(p) => {
                if !zero(&cross(&p.axis, &frame.axes[2], &mut dec)?) || zero(&p.axis) {
                    return Err("curved cap certificate: plane is not perpendicular to the cylinder axis".into());
                }
                heights.push(frame.local(&p.origin, &mut dec)?[2].clone());
            },
            _ => return Err("curved support certificate: spline and other analytic faces need a separate replacement and feature certificate".into()),
        }
    }
    if heights.len() != 2 || equal(&heights[0], &heights[1], &mut dec)? {
        return Err("curved cap certificate: two distinct cap heights are required".into());
    }
    let heights = [heights[0].clone(), heights[1].clone()];
    let cap = crate::gap::auto_cap();
    let mut physical = BTreeSet::new();
    for face in &source.faces {
        for edge in physical_edges(face, &mut dec)? {
            physical.insert(edge.edge.clone());
        }
    }
    let mut result = source.clone();
    let mut candidates = BTreeMap::new();
    for edge in result
        .faces
        .iter_mut()
        .flat_map(|f| &mut f.loops)
        .flat_map(|b| &mut b.edges)
    {
        checkpoint(&mut dec)?;
        if !physical.contains(&edge.edge) || !matches!(edge.curve, Curve::BSpline(_)) {
            continue;
        }
        let (target, nets) =
            replacement(edge, &frame, &mut dec).map_err(|e| format!("edge {}: {e}", edge.edge))?;
        candidates.insert(edge.edge.clone(), nets);
        *edge = target;
    }
    if candidates.is_empty() {
        return Err(
            "curved target certificate: no supported physical spline boundary to replace".into(),
        );
    }
    // The complete rings bound two disks and one cylindrical band. Radius,
    // height and distinct ring-vertex chords bound their widths and edge lengths.
    // Source polar order makes each repaired cap a radial graph. Moving its
    // boundary by at most delta reduces any neighbouring separation by at most
    // 2 delta, so 12 delta <= target lower bound proves the source tenth rule.
    let feature2 = feature(&result, &frame, &heights, &mut dec)?;
    let feature_lower = sqrt_enclosure(&feature2, &mut dec)?.0;
    let limit = minimum(
        cap.clone(),
        op(&feature_lower, &Rational::from_i64(12), '/', &mut dec)?,
        &mut dec,
    )?;
    let mut records = BTreeMap::new();
    let mut movement = Rational::zero();
    for (entity, nets) in candidates {
        let (lo, hi) = crate::incidence_circle::bound(nets, &limit, &mut dec)?;
        if compare(&hi, &limit, &mut dec)?.is_gt() {
            return Err(format!(
                "curved displacement or nearest-feature certificate: edge {entity} has parameter-correspondence displacement in [{}, {}] mm; permitted bound {} mm is unproved",
                crate::gap::decimal(&lo),
                crate::gap::decimal(&hi),
                limit
            ));
        }
        movement = maximum(movement, hi.clone(), &mut dec)?;
        if !hi.is_zero() {
            records.insert(entity, (lo, hi));
        }
    }
    let margin = op(&movement, &Rational::from_i64(2), '*', &mut dec)?;
    let feature = op(
        &sqrt_enclosure(&feature2, &mut dec)?.0,
        &margin,
        '-',
        &mut dec,
    )?;
    let required = op(&movement, &Rational::from_i64(10), '*', &mut dec)?;
    if feature.sign() != Sign::Positive || compare(&required, &feature, &mut dec)?.is_gt() {
        return Err("curved nearest-feature certificate: displacement exceeds a tenth of the source feature lower bound".into());
    }
    let feature2 = op(&feature, &feature, '*', &mut dec)?.to_string();
    let hash = crate::sha256_hex(&crate::canonical_bytes(&crate::StepFile {
        solids: vec![source.clone()],
        mm_per_unit: Rational::one(),
    }));
    let mut marks = MarkLog::new();
    for (entity, (lo, hi)) in records {
        let record = Record {
            solid: &source.name,
            entity: &entity,
            entity_kind: "edge",
            source_geometry_sha256: &hash,
            definition: "maximum pointwise parameter-correspondence displacement; upper bound also bounds both finite-curve Hausdorff directions; unchanged shared vertices and supports",
            deviation_lower_mm: crate::gap::decimal(&lo),
            deviation_upper_mm: crate::gap::decimal(&hi),
            approved_tolerance_mm: crate::gap::decimal(&cap),
            feature_size_squared_mm2: feature2.clone(),
        };
        marks.push(Mark::new(
            "import.incidence-curved-repair",
            BoundKind::Distance,
            hi,
            serde_json::to_string(&record).expect("serializable curved tolerance"),
        ));
    }
    Ok(Prepared {
        solid: result,
        marks,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rational_frames_accept_nonterminating_lengths_and_refuse_irrational_norms() {
        let mut dec = Decider::default();
        let axis = [
            Rational::zero(),
            Rational::zero(),
            Rational::from_ratio(1, 3).unwrap(),
        ];
        assert_eq!(
            unit(&axis, &mut dec).unwrap(),
            [0, 0, 1].map(Rational::from_i64)
        );
        assert!(unit(&[1, 1, 0].map(Rational::from_i64), &mut dec).is_err());
    }
}
