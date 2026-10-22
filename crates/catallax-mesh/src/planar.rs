//! Planar trims need exact predicates before positions are rounded for display.

use std::collections::BTreeMap;

use catallax_numbers::{Decider, Outcome, Rational, Real, Refusal, Sign, StrictSign};
use catallax_shapes::Plane;
use catallax_topology::{FaceId, Loop, Model, Sense};

use crate::{FaceMesh, MeshOptions, Triangle, chart, edges::Boundaries, failure};

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

pub(crate) fn build(
    model: &Model,
    face: FaceId,
    plane: &Plane,
    boundaries: &mut Boundaries,
    options: &MeshOptions,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z];
    let origin = [&plane.origin.x, &plane.origin.y, &plane.origin.z];
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
                points.entry(index).or_insert_with(|| {
                    let mut p = boundaries.mesh.positions[index].map(|x| Real::from(chart::q(x)));
                    let offset = (0..3)
                        .filter(|&i| i != drop)
                        .fold(Real::from(0), |s, i| s + normal[i] * &(&p[i] - origin[i]));
                    p[drop] = origin[drop] - &(&offset / normal[drop]);
                    p
                });
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
            .into_iter()
            .map(|(i, p)| (i, axes.map(|a| p[a].clone())))
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
    });
    Ok(())
}
