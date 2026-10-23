//! Planar trims need exact predicates before positions are rounded for display.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Decider, Outcome, Rational, Real, Refusal, Sign, StrictSign};
use catallax_shapes::Plane;
use catallax_topology::{FaceId, Loop, Model, Sense};

use crate::{
    FaceMesh, Mesh, MeshControl, MeshOptions, Triangle, UnavailableControl, chart,
    edges::Boundaries, failure,
};

type Point = [Real; 2];

fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        Outcome::Refused(r) => Err(r),
    }
}

fn cross(a: &Point, b: &Point, c: &Point) -> Real {
    (&b[0] - &a[0]) * (&c[1] - &a[1]) - (&b[1] - &a[1]) * (&c[0] - &a[0])
}

struct Polygon<'a> {
    points: BTreeMap<usize, Point>,
    dec: &'a mut Decider,
}

impl Polygon<'_> {
    fn turn(&mut self, a: usize, b: usize, c: usize) -> Result<Sign, Refusal> {
        sign(
            &cross(&self.points[&a], &self.points[&b], &self.points[&c]),
            self.dec,
        )
    }

    fn area(&self, ring: &[usize]) -> Real {
        let origin = &self.points[&ring[0]];
        ring.windows(2).skip(1).fold(Real::from(0), |s, edge| {
            s + cross(origin, &self.points[&edge[0]], &self.points[&edge[1]])
        })
    }

    fn contains(&mut self, ring: &[usize], p: &Point) -> Result<bool, Refusal> {
        let mut winding = 0i32;
        for i in 0..ring.len() {
            let a = &self.points[&ring[i]];
            let b = &self.points[&ring[(i + 1) % ring.len()]];
            let ay = sign(&(&a[1] - &p[1]), self.dec)?;
            let by = sign(&(&b[1] - &p[1]), self.dec)?;
            if ay != Sign::Positive && by == Sign::Positive {
                if sign(&cross(a, b, p), self.dec)? == Sign::Positive {
                    winding += 1;
                }
            } else if ay == Sign::Positive
                && by != Sign::Positive
                && sign(&cross(a, b, p), self.dec)? == Sign::Negative
            {
                winding -= 1;
            }
        }
        Ok(winding != 0)
    }

    fn blocks(&mut self, a: usize, b: usize, c: usize, d: usize) -> Result<bool, Refusal> {
        if a == c || a == d || b == c || b == d {
            return Ok(false);
        }
        let ab_c = self.turn(a, b, c)?;
        let ab_d = self.turn(a, b, d)?;
        if ab_c == ab_d && ab_c != Sign::Zero {
            return Ok(false);
        }
        let cd_a = self.turn(c, d, a)?;
        let cd_b = self.turn(c, d, b)?;
        if cd_a == cd_b && cd_a != Sign::Zero {
            return Ok(false);
        }
        // Collinear disjoint segments do not obstruct a bridge.
        if [ab_c, ab_d, cd_a, cd_b].iter().all(|s| *s == Sign::Zero) {
            for axis in 0..2 {
                let ac = sign(&(&self.points[&a][axis] - &self.points[&c][axis]), self.dec)?;
                let ad = sign(&(&self.points[&a][axis] - &self.points[&d][axis]), self.dec)?;
                let bc = sign(&(&self.points[&b][axis] - &self.points[&c][axis]), self.dec)?;
                let bd = sign(&(&self.points[&b][axis] - &self.points[&d][axis]), self.dec)?;
                if ac != Sign::Zero && ac == ad && ac == bc && ac == bd {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    fn in_cone(&mut self, ring: &[usize], i: usize, b: usize) -> Result<bool, Refusal> {
        let a = ring[i];
        let prev = ring[(i + ring.len() - 1) % ring.len()];
        let next = ring[(i + 1) % ring.len()];
        let left = self.turn(prev, a, b)? == Sign::Positive;
        let right = self.turn(a, next, b)? == Sign::Positive;
        Ok(if self.turn(prev, a, next)? != Sign::Negative {
            left && right
        } else {
            left || right
        })
    }

    fn bridge(
        &mut self,
        ring: &mut Vec<usize>,
        hole: &[usize],
        loops: &[Vec<usize>],
    ) -> Result<(), Refusal> {
        for i in 0..ring.len() {
            for j in 0..hole.len() {
                let (a, b) = (ring[i], hole[j]);
                // Repeated bridge endpoints have different interior sectors.
                // Crossing tests alone cannot choose the correct occurrence.
                if !self.in_cone(ring, i, b)? || !self.in_cone(hole, j, a)? {
                    continue;
                }
                let mut clear = true;
                for boundary in loops
                    .iter()
                    .map(Vec::as_slice)
                    .chain(std::iter::once(ring.as_slice()))
                {
                    for k in 0..boundary.len() {
                        if self.blocks(a, b, boundary[k], boundary[(k + 1) % boundary.len()])? {
                            clear = false;
                            break;
                        }
                    }
                    if !clear {
                        break;
                    }
                }
                if !clear {
                    continue;
                }
                let midpoint = std::array::from_fn(|k| {
                    (&self.points[&a][k] + &self.points[&b][k]) / Real::from(2)
                });
                let mut inside = false;
                for boundary in loops {
                    inside ^= self.contains(boundary, &midpoint)?;
                }
                if !inside {
                    continue;
                }
                let mut joined = Vec::with_capacity(ring.len() + hole.len() + 2);
                joined.extend_from_slice(&ring[..=i]);
                joined.extend(hole[j..].iter().chain(&hole[..=j]).copied());
                joined.extend_from_slice(&ring[i..]);
                *ring = joined;
                return Ok(());
            }
        }
        Err(failure("no certified interior bridge for a planar hole"))
    }

    fn clip(
        &mut self,
        mut ring: Vec<usize>,
        face: FaceId,
        triangles: &mut Vec<Triangle>,
    ) -> Result<(), Refusal> {
        triangles.reserve(ring.len().saturating_sub(2));
        while ring.len() > 3 {
            let mut ear = None;
            for i in 0..ring.len() {
                let a = ring[(i + ring.len() - 1) % ring.len()];
                let b = ring[i];
                let c = ring[(i + 1) % ring.len()];
                if self.turn(a, b, c)? != Sign::Positive {
                    continue;
                }
                let mut clear = true;
                for &p in &ring {
                    if p == a || p == b || p == c {
                        continue;
                    }
                    if self.turn(a, b, p)? != Sign::Negative
                        && self.turn(b, c, p)? != Sign::Negative
                        && self.turn(c, a, p)? != Sign::Negative
                    {
                        clear = false;
                        break;
                    }
                }
                if clear {
                    ear = Some((i, [a, b, c]));
                    break;
                }
            }
            let Some((i, vertices)) = ear else {
                return Err(failure("no certified ear in a planar boundary"));
            };
            triangles.push(Triangle { vertices, face });
            ring.remove(i);
        }
        if self.turn(ring[0], ring[1], ring[2])? != Sign::Positive {
            return Err(failure("degenerate final planar triangle"));
        }
        triangles.push(Triangle {
            vertices: [ring[0], ring[1], ring[2]],
            face,
        });
        Ok(())
    }
}

fn refine(
    triangles: &mut Vec<Triangle>,
    points: &mut BTreeMap<usize, [Real; 3]>,
    loops: &[Vec<usize>],
    mesh: &mut Mesh,
    options: &MeshOptions,
    dec: &Decider,
) -> Result<Rational, Refusal> {
    let Some(limit) = &options.max_edge_length else {
        return Ok(Rational::zero());
    };
    let limit2 = limit * limit;
    let boundary: BTreeSet<_> = loops
        .iter()
        .flat_map(|ring| {
            (0..ring.len()).map(|i| {
                let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
                (a.min(b), a.max(b))
            })
        })
        .collect();
    let half = Real::from(Rational::from_ratio(1, 2).unwrap());
    let mut error = Rational::zero();
    for _ in 0..48 {
        dec.check_deadline()?;
        let mut pending = BTreeSet::new();
        for triangle in triangles.iter() {
            for i in 0..3 {
                let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                if crate::controls::length_squared(mesh.positions[a], mesh.positions[b]) > limit2 {
                    let key = (a.min(b), a.max(b));
                    if boundary.contains(&key) {
                        return Err(failure(
                            "sampled planar boundary exceeds maximum edge length",
                        ));
                    }
                    pending.insert(key);
                }
            }
        }
        if pending.is_empty() {
            return Ok(error);
        }
        if pending.len() > options.max_vertices.saturating_sub(mesh.positions.len()) {
            return Err(failure("planar refinement reached max_vertices"));
        }
        let mut midpoints = BTreeMap::new();
        for (a, b) in pending {
            let p = std::array::from_fn::<_, 3, _>(|i| &(&points[&a][i] + &points[&b][i]) * &half);
            let balls = [
                dec.enclose_at(&p[0], chart::BITS)?,
                dec.enclose_at(&p[1], chart::BITS)?,
                dec.enclose_at(&p[2], chart::BITS)?,
            ];
            let position = chart::position(&balls)?;
            let bits = position.map(f64::to_bits);
            if [a, b]
                .iter()
                .any(|&i| mesh.positions[i].map(f64::to_bits) == bits)
            {
                return Err(failure(
                    "coordinate precision cannot shorten a planar mesh edge",
                ));
            }
            error = error.max(chart::residual(&position, &balls));
            if error > options.chord_tolerance {
                return Err(failure("planar midpoint rounding exceeds chord tolerance"));
            }
            let index = mesh.positions.len();
            mesh.positions.push(position);
            points.insert(index, p);
            midpoints.insert((a, b), index);
        }
        let mut refined = Vec::with_capacity(triangles.len() + 2 * midpoints.len());
        for triangle in triangles.iter() {
            let vertices = triangle.vertices;
            let mids = std::array::from_fn::<_, 3, _>(|i| {
                let (a, b) = (vertices[i], vertices[(i + 1) % 3]);
                midpoints.get(&(a.min(b), a.max(b))).copied()
            });
            let count = mids.iter().filter(|m| m.is_some()).count();
            let mut emit = |vertices| {
                refined.push(Triangle {
                    vertices,
                    face: triangle.face,
                })
            };
            match count {
                0 => emit(vertices),
                1 => {
                    let i = mids.iter().position(Option::is_some).unwrap();
                    let [a, b, c] = std::array::from_fn(|j| vertices[(i + j) % 3]);
                    let m = mids[i].unwrap();
                    emit([a, m, c]);
                    emit([m, b, c]);
                }
                2 => {
                    let i = (0..3)
                        .find(|&i| mids[i].is_some() && mids[(i + 1) % 3].is_some())
                        .unwrap();
                    let [a, b, c] = std::array::from_fn(|j| vertices[(i + j) % 3]);
                    let ab = mids[i].unwrap();
                    let bc = mids[(i + 1) % 3].unwrap();
                    emit([a, ab, c]);
                    emit([ab, bc, c]);
                    emit([ab, b, bc]);
                }
                3 => {
                    let [a, b, c] = vertices;
                    let [ab, bc, ca] = mids.map(Option::unwrap);
                    emit([a, ab, ca]);
                    emit([ab, b, bc]);
                    emit([ca, bc, c]);
                    emit([ab, bc, ca]);
                }
                _ => unreachable!(),
            }
        }
        *triangles = refined;
    }
    Err(failure(
        "planar refinement reached its deterministic subdivision limit",
    ))
}

pub(crate) fn build(
    model: &Model,
    face: FaceId,
    plane: &Plane,
    boundaries: &mut Boundaries,
    options: &MeshOptions,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z];
    let mut candidates = [0, 1, 2];
    candidates.sort_by(|&a, &b| {
        normal[b]
            .to_f64()
            .abs()
            .total_cmp(&normal[a].to_f64().abs())
    });
    let mut projection = None;
    for axis in candidates {
        let direction = sign(normal[axis], dec)?;
        if direction != Sign::Zero {
            projection = Some((axis, direction));
            break;
        }
    }
    let (drop, direction) = projection.ok_or_else(|| failure("zero planar normal"))?;
    let mut axes = [(drop + 1) % 3, (drop + 2) % 3];
    if (direction == Sign::Negative) != (model.face(face).sense == Sense::Reversed) {
        axes.swap(0, 1);
    }
    let mut points = BTreeMap::<usize, [Real; 3]>::new();
    let mut loops = Vec::new();
    for &l in &model.face(face).loops {
        let Loop::Coedges(coedges) = model.loop_(l) else {
            return Err(failure("a planar face has a vertex-only boundary"));
        };
        let mut ring = Vec::new();
        for &coedge in coedges {
            let coedge = model.coedge(coedge);
            let edge = model.edge(coedge.edge);
            let sampled = &boundaries.mesh.edges[boundaries.edges[&coedge.edge]].vertices;
            for v in [edge.start, edge.end].into_iter().flatten() {
                points.entry(boundaries.vertices[&v]).or_insert_with(|| {
                    let p = model.geom().position(model.vertex(v).point);
                    [p.x.clone(), p.y.clone(), p.z.clone()]
                });
            }
            for &index in sampled {
                if let std::collections::btree_map::Entry::Vacant(entry) = points.entry(index) {
                    let p = boundaries
                        .exact_positions
                        .get(&index)
                        .ok_or_else(|| failure("planar boundary lacks its exact trace sample"))?;
                    entry.insert([p.x.clone(), p.y.clone(), p.z.clone()]);
                }
            }
            if coedge.sense == Sense::Same {
                ring.extend_from_slice(&sampled[..sampled.len() - 1]);
            } else {
                ring.extend(sampled[1..].iter().rev().copied());
            }
        }
        if ring.len() < 3 {
            return Err(failure("planar boundary has fewer than three vertices"));
        }
        loops.push(ring);
    }
    if boundaries.mesh.positions.len() > options.max_vertices {
        return Err(failure("planar triangulation reached max_vertices"));
    }
    // A planar support is affine. Convex interpolation of endpoint residuals
    // bounds every point of a triangle, including strips collapsed by binary64.
    let mut error = Rational::zero();
    for (&index, p) in &points {
        let balls = [
            dec.enclose_at(&p[0], chart::BITS)?,
            dec.enclose_at(&p[1], chart::BITS)?,
            dec.enclose_at(&p[2], chart::BITS)?,
        ];
        error = error.max(chart::residual(&boundaries.mesh.positions[index], &balls));
    }
    if error > options.chord_tolerance {
        return Err(failure("planar position rounding exceeds chord tolerance"));
    }
    let mut polygon = Polygon {
        points: points
            .iter()
            .map(|(&i, p)| (i, axes.map(|a| p[a].clone())))
            .collect(),
        dec,
    };
    let mut outer = Vec::new();
    let mut holes = Vec::new();
    let areas: Vec<_> = loops.iter().map(|ring| polygon.area(ring)).collect();
    for (i, area) in areas.iter().enumerate() {
        match sign(area, polygon.dec)? {
            Sign::Positive => outer.push(i),
            Sign::Negative => holes.push(i),
            Sign::Zero => return Err(failure("zero-area planar boundary")),
        }
    }
    let mut assigned = BTreeMap::<usize, Vec<usize>>::new();
    for hole in holes {
        let p = polygon.points[&loops[hole][0]].clone();
        let mut owner: Option<usize> = None;
        for &candidate in &outer {
            if polygon.contains(&loops[candidate], &p)?
                && (owner.is_none()
                    || sign(&(&areas[candidate] - &areas[owner.unwrap()]), polygon.dec)?
                        == Sign::Negative)
            {
                owner = Some(candidate);
            }
        }
        assigned
            .entry(owner.ok_or_else(|| failure("planar hole has no outer boundary"))?)
            .or_default()
            .push(hole);
    }
    let mut triangles = Vec::new();
    for i in outer {
        let mut ring = loops[i].clone();
        if let Some(holes) = assigned.get(&i) {
            for &hole in holes {
                polygon.bridge(&mut ring, &loops[hole], &loops)?;
            }
        }
        polygon.clip(ring, face, &mut triangles)?;
    }
    error = error.max(refine(
        &mut triangles,
        &mut points,
        &loops,
        &mut boundaries.mesh,
        options,
        polygon.dec,
    )?);
    let chart = chart::Chart::new(&catallax_shapes::Surface::Plane(plane.clone()), polygon.dec)?;
    let mut aspect = Some(Rational::zero());
    let mut angular_reached = true;
    for triangle in &triangles {
        let mut positions = triangle.vertices.map(|i| boundaries.mesh.positions[i]);
        // Quality is compared with the carrier's normal, before face reversal.
        if model.face(face).sense == Sense::Reversed {
            positions.swap(1, 2);
        }
        let (ratio, reached) = crate::controls::quality(
            &chart,
            &[[0.0; 2]; 3],
            positions,
            options.angular_deviation.as_ref(),
            polygon.dec,
        )?;
        aspect = aspect.zip(ratio).map(|(a, b)| a.max(b));
        angular_reached &= reached;
    }
    if options
        .max_aspect_ratio
        .as_ref()
        .is_some_and(|limit| aspect.as_ref().is_none_or(|value| value > limit))
    {
        boundaries.mesh.unavailable_controls.push(UnavailableControl {
            face,
            control: MeshControl::AspectRatio,
            reason: "rounded planar triangles have no aspect-ratio bound within the requested limit",
        });
    }
    if options.angular_deviation.is_some() && !angular_reached {
        boundaries
            .mesh
            .unavailable_controls
            .push(UnavailableControl {
                face,
                control: MeshControl::AngularDeviation,
                reason: "rounded planar triangles have no certified facet-normal bound",
            });
    }
    for t in &mut triangles {
        let min = (0..3).min_by_key(|&i| t.vertices[i]).unwrap();
        t.vertices.rotate_left(min);
    }
    triangles.sort_unstable_by_key(|t| t.vertices);
    let start = boundaries.mesh.triangles.len();
    boundaries.mesh.triangles.extend(triangles);
    boundaries.mesh.faces.push(FaceMesh {
        face,
        triangles: start..boundaries.mesh.triangles.len(),
        max_chord_error: error,
        max_aspect_ratio: aspect,
        max_angular_deviation: options
            .angular_deviation
            .clone()
            .filter(|_| angular_reached),
    });
    Ok(())
}
