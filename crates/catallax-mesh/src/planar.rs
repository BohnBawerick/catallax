//! Planar trims need exact predicates before positions are rounded for display.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use catallax_numbers::{Decider, Outcome, Rational, Real, Refusal, Sign, StrictSign};
use catallax_shapes::Plane;
use catallax_topology::{FaceId, Loop, Model, Sense};

use crate::{
    FaceMesh, Mesh, MeshControl, MeshOptions, Triangle, UnavailableControl, chart,
    edges::{Boundaries, Position},
    failure,
};

#[derive(Clone)]
struct Point {
    position: Position,
    axes: [usize; 2],
    bounds: [[f64; 2]; 2],
}

impl Point {
    fn exact(&self, dec: &Decider) -> Result<[Real; 2], Refusal> {
        let p = self.position.exact(dec)?;
        Ok(self.axes.map(|axis| p[axis].clone()))
    }
}

fn difference(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [(a[0] - b[1]).next_down(), (a[1] - b[0]).next_up()]
}

fn product(a: [f64; 2], b: [f64; 2]) -> Option<[f64; 2]> {
    let values = [a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]];
    if values.iter().any(|v| v.is_nan()) {
        return None;
    }
    Some([
        values
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min)
            .next_down(),
        values
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
            .next_up(),
    ])
}

fn filtered_turn(a: &Point, b: &Point, c: &Point) -> Option<Sign> {
    let ab = std::array::from_fn::<_, 2, _>(|i| difference(b.bounds[i], a.bounds[i]));
    let ac = std::array::from_fn::<_, 2, _>(|i| difference(c.bounds[i], a.bounds[i]));
    let determinant = difference(product(ab[0], ac[1])?, product(ab[1], ac[0])?);
    if determinant[0] > 0.0 {
        Some(Sign::Positive)
    } else if determinant[1] < 0.0 {
        Some(Sign::Negative)
    } else {
        None
    }
}

fn dec_bounds(value: &Real, dec: &Decider) -> Result<[f64; 2], Refusal> {
    let ball = dec.enclose_at(value, chart::BITS)?;
    Ok([
        ball.lo().to_f64_nearest().next_down(),
        ball.hi().to_f64_nearest().next_up(),
    ])
}

#[derive(Clone, Copy)]
struct CandidateVertex {
    position: spade::Point2<f64>,
    index: usize,
}

impl spade::HasPosition for CandidateVertex {
    type Scalar = f64;
    fn position(&self) -> spade::Point2<f64> {
        self.position
    }
}

fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        Outcome::Refused(r) => Err(r),
    }
}

fn cross(a: &Point, b: &Point, c: &Point, dec: &Decider) -> Result<Real, Refusal> {
    let a = a.exact(dec)?;
    let b = b.exact(dec)?;
    let c = c.exact(dec)?;
    Ok((&b[0] - &a[0]) * (&c[1] - &a[1]) - (&b[1] - &a[1]) * (&c[0] - &a[0]))
}

struct Polygon<'a> {
    points: BTreeMap<usize, Point>,
    dec: &'a mut Decider,
}

impl Polygon<'_> {
    fn turn(&mut self, a: usize, b: usize, c: usize) -> Result<Sign, Refusal> {
        if a == b || b == c || a == c {
            return Ok(Sign::Zero);
        }
        if let Some(value) = filtered_turn(&self.points[&a], &self.points[&b], &self.points[&c]) {
            return Ok(value);
        }
        sign(
            &cross(
                &self.points[&a],
                &self.points[&b],
                &self.points[&c],
                self.dec,
            )?,
            self.dec,
        )
    }

    fn area(&self, ring: &[usize]) -> Result<Real, Refusal> {
        let origin = &self.points[&ring[0]];
        ring.windows(2).skip(1).try_fold(Real::from(0), |s, edge| {
            Ok(s + cross(
                origin,
                &self.points[&edge[0]],
                &self.points[&edge[1]],
                self.dec,
            )?)
        })
    }

    fn contains(&mut self, ring: &[usize], p: &Point) -> Result<bool, Refusal> {
        let mut winding = 0i32;
        let p_exact = p.exact(self.dec)?;
        for i in 0..ring.len() {
            let a = &self.points[&ring[i]];
            let b = &self.points[&ring[(i + 1) % ring.len()]];
            let ay = sign(&(&a.exact(self.dec)?[1] - &p_exact[1]), self.dec)?;
            let by = sign(&(&b.exact(self.dec)?[1] - &p_exact[1]), self.dec)?;
            if ay != Sign::Positive && by == Sign::Positive {
                if sign(&cross(a, b, p, self.dec)?, self.dec)? == Sign::Positive {
                    winding += 1;
                }
            } else if ay == Sign::Positive
                && by != Sign::Positive
                && sign(&cross(a, b, p, self.dec)?, self.dec)? == Sign::Negative
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
            let a = self.points[&a].exact(self.dec)?;
            let b = self.points[&b].exact(self.dec)?;
            let c = self.points[&c].exact(self.dec)?;
            let d = self.points[&d].exact(self.dec)?;
            for axis in 0..2 {
                let ac = sign(&(&a[axis] - &c[axis]), self.dec)?;
                let ad = sign(&(&a[axis] - &d[axis]), self.dec)?;
                let bc = sign(&(&b[axis] - &c[axis]), self.dec)?;
                let bd = sign(&(&b[axis] - &d[axis]), self.dec)?;
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

    fn coordinate_order(&mut self, a: usize, b: usize, axis: usize) -> Result<Sign, Refusal> {
        let (a, b) = (&self.points[&a], &self.points[&b]);
        if a.bounds[axis][1] < b.bounds[axis][0] {
            Ok(Sign::Negative)
        } else if a.bounds[axis][0] > b.bounds[axis][1] {
            Ok(Sign::Positive)
        } else {
            sign(
                &(&a.exact(self.dec)?[axis] - &b.exact(self.dec)?[axis]),
                self.dec,
            )
        }
    }

    fn hull(
        &mut self,
        positions: &[[f64; 3]],
        axes: [usize; 2],
    ) -> Result<Option<BTreeSet<(usize, usize)>>, Refusal> {
        let mut order: Vec<_> = self.points.keys().copied().collect();
        order.sort_unstable_by(|&a, &b| {
            positions[a][axes[0]]
                .total_cmp(&positions[b][axes[0]])
                .then_with(|| positions[a][axes[1]].total_cmp(&positions[b][axes[1]]))
                .then_with(|| a.cmp(&b))
        });
        // Floating sorting is only a proposal; certify its exact lexicographic order.
        for pair in order.windows(2) {
            let x = self.coordinate_order(pair[0], pair[1], 0)?;
            if x == Sign::Positive
                || (x == Sign::Zero
                    && self.coordinate_order(pair[0], pair[1], 1)? != Sign::Negative)
            {
                return Ok(None);
            }
        }
        let mut hull = Vec::with_capacity(order.len() * 2);
        for &index in &order {
            while hull.len() >= 2
                && self.turn(hull[hull.len() - 2], hull[hull.len() - 1], index)? == Sign::Negative
            {
                hull.pop();
            }
            hull.push(index);
        }
        let lower = hull.len();
        for &index in order.iter().rev().skip(1) {
            while hull.len() > lower
                && self.turn(hull[hull.len() - 2], hull[hull.len() - 1], index)? == Sign::Negative
            {
                hull.pop();
            }
            hull.push(index);
        }
        hull.pop();
        let edges: BTreeSet<_> = (0..hull.len())
            .map(|i| {
                let (a, b) = (hull[i], hull[(i + 1) % hull.len()]);
                (a.min(b), a.max(b))
            })
            .collect();
        Ok((hull.len() >= 3 && edges.len() == hull.len()).then_some(edges))
    }

    fn candidate(
        &mut self,
        loops: &[Vec<usize>],
        face: FaceId,
        positions: &[[f64; 3]],
        axes: [usize; 2],
    ) -> Result<Option<Vec<Triangle>>, Refusal> {
        use spade::Triangulation;
        let mut cdt = spade::ConstrainedDelaunayTriangulation::<CandidateVertex>::new();
        let mut handles = BTreeMap::new();
        for &index in self.points.keys() {
            self.dec.check_deadline()?;
            let p = positions[index];
            let previous_count = cdt.num_vertices();
            let Ok(handle) = cdt.insert(CandidateVertex {
                position: spade::Point2::new(p[axes[0]], p[axes[1]]),
                index,
            }) else {
                return Ok(None);
            };
            if cdt.num_vertices() == previous_count {
                return Ok(None);
            }
            handles.insert(index, handle);
        }
        let mut boundary = BTreeSet::new();
        for ring in loops {
            for i in 0..ring.len() {
                self.dec.check_deadline()?;
                let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
                if a == b || boundary.contains(&(b, a)) || !boundary.insert((a, b)) {
                    return Ok(None);
                }
                if !cdt.can_add_constraint(handles[&a], handles[&b]) {
                    return Ok(None);
                }
                cdt.add_constraint(handles[&a], handles[&b]);
            }
        }
        let Some(expected_hull) = self.hull(positions, axes)? else {
            return Ok(None);
        };
        let mut labels = vec![None; cdt.num_all_faces()];
        let mut queue = VecDeque::new();
        let mut group = Vec::new();
        for start in cdt.inner_faces() {
            if labels[start.index()].is_some() {
                continue;
            }
            group.clear();
            queue.push_back(start.fix());
            labels[start.index()] = Some(false);
            let mut seed = None;
            while let Some(current) = queue.pop_front() {
                group.push(current.index());
                for edge in cdt.face(current).adjacent_edges() {
                    if edge.is_constraint_edge() {
                        let pair = (edge.from().data().index, edge.to().data().index);
                        let left = if boundary.contains(&pair) {
                            true
                        } else if boundary.contains(&(pair.1, pair.0)) {
                            false
                        } else {
                            return Ok(None);
                        };
                        if seed.is_some_and(|value| value != left) {
                            return Ok(None);
                        }
                        seed = Some(left);
                    } else if let Some(other) = edge.rev().face().as_inner()
                        && labels[other.index()].is_none()
                    {
                        labels[other.index()] = Some(false);
                        queue.push_back(other.fix());
                    }
                }
            }
            for index in &group {
                labels[*index] = Some(seed.unwrap_or(false));
            }
        }
        let mut triangles = Vec::new();
        let mut incidence = BTreeMap::<(usize, usize), (usize, isize)>::new();
        let mut actual_hull = BTreeSet::new();
        for triangle in cdt.inner_faces() {
            let vertices = triangle.vertices().map(|vertex| vertex.data().index);
            if self.turn(vertices[0], vertices[1], vertices[2])? != Sign::Positive {
                return Ok(None);
            }
            for edge in triangle.adjacent_edges() {
                if edge.rev().face().as_inner().is_none() {
                    let (a, b) = (edge.from().data().index, edge.to().data().index);
                    actual_hull.insert((a.min(b), a.max(b)));
                }
            }
            if labels[triangle.index()] != Some(true) {
                continue;
            }
            for i in 0..3 {
                let (a, b) = (vertices[i], vertices[(i + 1) % 3]);
                let uses = incidence.entry((a.min(b), a.max(b))).or_default();
                uses.0 += 1;
                uses.1 += if a < b { 1 } else { -1 };
            }
            triangles.push(Triangle { vertices, face });
        }
        // Certify the entire disk embedding, including triangles outside the
        // trims: a convex exact hull and positive triangles rule out overlaps.
        if actual_hull != expected_hull {
            return Ok(None);
        }
        // Positive exact triangles with precisely the oriented input boundary
        // have the same winding function as the trims. Float CDT coordinates
        // propose connectivity only; neither orientation nor incidence is trusted.
        for (&(a, b), &(count, winding)) in &incidence {
            let expected = if boundary.contains(&(a, b)) {
                (1, 1)
            } else if boundary.contains(&(b, a)) {
                (1, -1)
            } else {
                (2, 0)
            };
            if (count, winding) != expected {
                return Ok(None);
            }
        }
        if boundary
            .iter()
            .any(|&(a, b)| !incidence.contains_key(&(a.min(b), a.max(b))))
        {
            return Ok(None);
        }
        Ok(Some(triangles))
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
                let a = self.points[&a].exact(self.dec)?;
                let b = self.points[&b].exact(self.dec)?;
                let coordinates = std::array::from_fn(|k| (&a[k] + &b[k]) / Real::from(2));
                let bounds = coordinates
                    .each_ref()
                    .map(|value| dec_bounds(value, self.dec));
                let [x, y] = bounds;
                let midpoint = Point {
                    position: Position::Exact([
                        coordinates[0].clone(),
                        coordinates[1].clone(),
                        Real::zero(),
                    ]),
                    axes: [0, 1],
                    bounds: [x?, y?],
                };
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
    points: &mut BTreeMap<usize, Position>,
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
            let pa = points[&a].exact(dec)?;
            let pb = points[&b].exact(dec)?;
            let p = std::array::from_fn::<_, 3, _>(|i| &(&pa[i] + &pb[i]) * &half);
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
            points.insert(index, Position::Exact(p));
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

fn split_rail_ears(
    triangles: &mut Vec<Triangle>,
    points: &mut BTreeMap<usize, Position>,
    enclosures: &BTreeMap<usize, [catallax_numbers::Ball; 3]>,
    rails: &crate::edges::Rails,
    mesh: &mut Mesh,
    options: &MeshOptions,
    dec: &Decider,
) -> Result<Rational, Refusal> {
    let third = Rational::from_ratio(1, 3).unwrap();
    let factor = catallax_numbers::Ball::from_rational(&third, chart::BITS);
    let mut error = Rational::zero();
    let original = triangles.len();
    for i in 0..original {
        dec.check_deadline()?;
        let vertices = triangles[i].vertices;
        if !rails.contains_triangle(vertices) {
            continue;
        }
        if mesh.positions.len() >= options.max_vertices {
            return Err(failure("planar rail-ear refinement reached max_vertices"));
        }
        let [a, b, c] = vertices;
        // The exact affine centroid lies strictly inside the already-certified
        // planar triangle. Enclose that same centroid directly from the cached
        // endpoints, without expanding three trace samples into exact rationals.
        // Length refinement still retains the exact construction it consumes.
        let mut exact = None;
        let balls = if options.max_edge_length.is_some() {
            let [pa, pb, pc] = [a, b, c].map(|v| points[&v].exact(dec));
            let [pa, pb, pc] = [pa?, pb?, pc?];
            let third = Real::from(third.clone());
            let p = std::array::from_fn::<_, 3, _>(|axis| {
                &(&(&pa[axis] + &pb[axis]) + &pc[axis]) * &third
            });
            let balls = [
                dec.enclose_at(&p[0], chart::BITS)?,
                dec.enclose_at(&p[1], chart::BITS)?,
                dec.enclose_at(&p[2], chart::BITS)?,
            ];
            exact = Some(p);
            balls
        } else {
            std::array::from_fn(|axis| {
                &(&(&enclosures[&a][axis] + &enclosures[&b][axis]) + &enclosures[&c][axis])
                    * &factor
            })
        };
        let position = chart::position(&balls)?;
        error = error.max(chart::residual(&position, &balls));
        if error > options.chord_tolerance {
            return Err(failure("planar rail-ear rounding exceeds chord tolerance"));
        }
        let index = mesh.positions.len();
        mesh.positions.push(position);
        if let Some(p) = exact {
            points.insert(index, Position::Exact(p));
        }
        let face = triangles[i].face;
        triangles[i].vertices = [a, b, index];
        triangles.push(Triangle {
            vertices: [b, c, index],
            face,
        });
        triangles.push(Triangle {
            vertices: [c, a, index],
            face,
        });
    }
    Ok(error)
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
    let mut points = BTreeMap::<usize, Position>::new();
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
                    Position::Exact([p.x.clone(), p.y.clone(), p.z.clone()])
                });
            }
            for &index in sampled {
                if let std::collections::btree_map::Entry::Vacant(entry) = points.entry(index) {
                    let p = boundaries
                        .exact_positions
                        .get(&index)
                        .ok_or_else(|| failure("planar boundary lacks its exact trace sample"))?;
                    entry.insert(p.clone());
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
    let mut projected = BTreeMap::new();
    for (&index, p) in &points {
        let balls = match boundaries.enclosures.entry(index) {
            std::collections::btree_map::Entry::Occupied(entry) => entry.into_mut(),
            std::collections::btree_map::Entry::Vacant(entry) => entry.insert(p.enclose(dec)?),
        };
        error = error.max(chart::residual(&boundaries.mesh.positions[index], balls));
        projected.insert(
            index,
            Point {
                position: p.clone(),
                axes,
                bounds: axes.map(|a| {
                    [
                        balls[a].lo().to_f64_nearest().next_down(),
                        balls[a].hi().to_f64_nearest().next_up(),
                    ]
                }),
            },
        );
    }
    if error > options.chord_tolerance {
        return Err(failure("planar position rounding exceeds chord tolerance"));
    }
    let mut polygon = Polygon {
        points: projected,
        dec,
    };
    let mut triangles = if let Some(triangles) =
        polygon.candidate(&loops, face, &boundaries.mesh.positions, axes)?
    {
        triangles
    } else {
        let mut outer = Vec::new();
        let mut holes = Vec::new();
        let areas: Vec<_> = loops
            .iter()
            .map(|ring| polygon.area(ring))
            .collect::<Result<_, _>>()?;
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
        triangles
    };
    error = error.max(split_rail_ears(
        &mut triangles,
        &mut points,
        &boundaries.enclosures,
        &boundaries.rails,
        &mut boundaries.mesh,
        options,
        polygon.dec,
    )?);
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
