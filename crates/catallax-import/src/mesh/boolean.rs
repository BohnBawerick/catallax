//! Booleans of two closed triangle meshes, RECIPE.md section 4.21.
//!
//! The pipeline is SMLib's `SmPolyMerge::ManifoldBoolean`: find candidate
//! triangle pairs with a box tree, intersect each pair, insert the
//! intersection edges into both meshes, sort the pieces into inside and
//! outside of the other solid by a flood fill that stops at intersection
//! edges, then keep the pieces the operation asks for, the tool's reversed for
//! a difference.
//!
//! SMLib runs every step in binary64 and joins points within its tolerances.
//! A mesh import joins nothing by distance (section 4.19), so here every
//! decision is an exact orientation sign and every new vertex is the exact
//! rational point where an edge crosses a plane. That is the corefinement of
//! CGAL's Polygon Mesh Processing package under its exact-construction kernel.
//! Each cut triangle is retriangulated in its own plane: its new points go in
//! by splitting, and each intersection segment is recovered by edge flips
//! (Sloan, "A fast algorithm for generating constrained Delaunay
//! triangulations", 1993). A piece is classified by the exact winding number of
//! one interior point, counted along the read's symbolically perturbed ray.
//!
//! Every input has passed the closed-mesh proof of section 4.19, and the result
//! passes it again before it is returned, so the check runs once per result.

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet, VecDeque};

use catallax_numbers::{Decider, Rational, Refusal, Sign};

use super::predicates::{self, P, in_triangle_2d, orient2d, orient3d, projection, record_exact};
use super::{ClosedMesh, Meter, outward_box};

/// Which material the result keeps.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind {
    /// Material in either mesh.
    Union,
    /// Material in the target and not in the tool.
    Subtract,
    /// Material in both meshes.
    Intersect,
}

/// The validated result and where each of its triangles came from.
#[derive(Clone, Debug)]
pub struct Outcome {
    /// The result, in millimetres.
    pub mesh: ClosedMesh,
    /// For each result triangle, the operand (0 target, 1 tool) and that
    /// operand's triangle it lies in.
    pub sources: Vec<(u8, u32)>,
}

fn refuse(code: &'static str, detail: impl Into<String>) -> Refusal {
    Refusal::new(code, detail)
}

fn internal(detail: impl Into<String>) -> Refusal {
    refuse("mesh-boolean.internal", detail)
}

fn flip(s: Sign) -> Sign {
    match s {
        Sign::Positive => Sign::Negative,
        Sign::Negative => Sign::Positive,
        Sign::Zero => Sign::Zero,
    }
}

/// Exact points, each stored once, so equal points are one vertex.
#[derive(Default)]
struct Store {
    exact: Vec<[Rational; 3]>,
    approx: Vec<[f64; 3]>,
    close: Vec<bool>,
    index: HashMap<[Rational; 3], u32>,
}

impl Store {
    fn add(&mut self, p: [Rational; 3]) -> u32 {
        if let Some(&i) = self.index.get(&p) {
            return i;
        }
        let i = self.exact.len() as u32;
        let f = p.each_ref().map(Rational::to_f64_nearest);
        self.close.push(predicates::close(&p, &f));
        self.approx.push(f);
        self.index.insert(p.clone(), i);
        self.exact.push(p);
        i
    }

    fn p(&self, i: u32) -> P<'_> {
        P {
            e: &self.exact[i as usize],
            f: self.approx[i as usize],
            close: self.close[i as usize],
        }
    }

    /// Lexicographic order of the exact coordinates. On a line it is
    /// monotone, so it orders points along any segment.
    fn cmp(&self, a: u32, b: u32) -> Ordering {
        if a == b {
            return Ordering::Equal;
        }
        let (x, y) = (&self.exact[a as usize], &self.exact[b as usize]);
        record_exact(3, x.iter().chain(y));
        x.cmp(y)
    }

    fn describe(&self, i: u32) -> String {
        let [x, y, z] = self.approx[i as usize];
        format!("({x}, {y}, {z})")
    }
}

/// The exact `det[b - a, c - a, d - a]`.
fn orient_value(a: P, b: P, c: P, d: P) -> Rational {
    record_exact(23, a.e.iter().chain(b.e).chain(c.e).chain(d.e));
    let diff = |p: P| std::array::from_fn::<Rational, 3, _>(|i| &p.e[i] - &a.e[i]);
    let (u, v, w) = (diff(b), diff(c), diff(d));
    let det = &(&u[0] * &(&(&v[1] * &w[2]) - &(&v[2] * &w[1])))
        - &(&u[1] * &(&(&v[0] * &w[2]) - &(&v[2] * &w[0])));
    &det + &(&u[2] * &(&(&v[0] * &w[1]) - &(&v[1] * &w[0])))
}

/// What two triangles share.
enum Meet {
    Empty,
    Point(u32),
    Segment(u32, u32),
    /// Coplanar with overlapping interiors.
    Overlap,
}

/// Signs of `points`' corners against the plane of `plane`.
fn sides(store: &Store, plane: [u32; 3], points: [u32; 3]) -> [Sign; 3] {
    let [a, b, c] = plane.map(|v| store.p(v));
    points.map(|v| orient3d(a, b, c, store.p(v)))
}

fn one_sided(s: [Sign; 3]) -> bool {
    s[0] != Sign::Zero && s.iter().all(|&x| x == s[0])
}

/// The points of triangle `tri` on the plane of `plane`: corners on it and
/// the exact crossings of edges whose ends lie strictly on opposite sides.
fn cut(store: &mut Store, tri: [u32; 3], plane: [u32; 3], signs: [Sign; 3]) -> Vec<u32> {
    let mut out: Vec<u32> = (0..3)
        .filter(|&i| signs[i] == Sign::Zero)
        .map(|i| tri[i])
        .collect();
    for i in 0..3 {
        let j = (i + 1) % 3;
        if signs[i] == Sign::Zero || signs[j] == Sign::Zero || signs[i] == signs[j] {
            continue;
        }
        let point = {
            let [a, b, c] = plane.map(|v| store.p(v));
            let (p, q) = (store.p(tri[i]), store.p(tri[j]));
            let (op, oq) = (orient_value(a, b, c, p), orient_value(a, b, c, q));
            record_exact(8, p.e.iter().chain(q.e).chain([&op, &oq]));
            let t = op
                .checked_div(&(&op - &oq))
                .expect("strictly opposite sides");
            std::array::from_fn(|k| &p.e[k] + &(&(&q.e[k] - &p.e[k]) * &t))
        };
        out.push(store.add(point));
    }
    out
}

/// The first and last of collinear points.
fn span(store: &Store, points: &[u32]) -> (u32, u32) {
    let mut lo = points[0];
    let mut hi = points[0];
    for &p in &points[1..] {
        if store.cmp(p, lo) == Ordering::Less {
            lo = p;
        }
        if store.cmp(p, hi) == Ordering::Greater {
            hi = p;
        }
    }
    (lo, hi)
}

fn interval(store: &Store, lo: u32, hi: u32) -> Meet {
    match store.cmp(lo, hi) {
        Ordering::Greater => Meet::Empty,
        Ordering::Equal => Meet::Point(lo),
        Ordering::Less => Meet::Segment(lo, hi),
    }
}

/// Exactly what triangles `a` and `b` share (Guigue and Devillers, "Fast and
/// robust triangle-triangle overlap test using orientation predicates", 2003,
/// for the plane tests; the shared part is the overlap of each triangle's cut
/// by the other's plane, both on the planes' common line).
fn meet(store: &mut Store, a: [u32; 3], b: [u32; 3]) -> Meet {
    let sa = sides(store, b, a);
    if one_sided(sa) {
        return Meet::Empty;
    }
    let sb = sides(store, a, b);
    if one_sided(sb) {
        return Meet::Empty;
    }
    if sa.iter().all(|&s| s == Sign::Zero) {
        return coplanar(store, a, b);
    }
    let ca = cut(store, a, b, sa);
    let cb = cut(store, b, a, sb);
    let (alo, ahi) = span(store, &ca);
    let (blo, bhi) = span(store, &cb);
    let lo = if store.cmp(alo, blo) == Ordering::Less {
        blo
    } else {
        alo
    };
    let hi = if store.cmp(ahi, bhi) == Ordering::Greater {
        bhi
    } else {
        ahi
    };
    interval(store, lo, hi)
}

/// Two triangles in one plane: an overlap of interiors, or the segment or
/// point where they touch. Convex sets with disjoint interiors are separated
/// by an edge line of one of them.
fn coplanar(store: &Store, a: [u32; 3], b: [u32; 3]) -> Meet {
    let pa = a.map(|v| store.p(v));
    let pb = b.map(|v| store.p(v));
    let (k, s_a) = projection(pa[0], pa[1], pa[2]).expect("validated triangles");
    let s_b = orient2d(k, pb[0], pb[1], pb[2]);
    let separated = |t: [P; 3], s: Sign, o: [P; 3]| {
        (0..3).any(|i| o.iter().all(|&q| orient2d(k, t[i], t[(i + 1) % 3], q) != s))
    };
    if !separated(pa, s_a, pb) && !separated(pb, s_b, pa) {
        return Meet::Overlap;
    }
    let mut touch: Vec<u32> = (0..3)
        .filter(|&i| in_triangle_2d(k, s_b, pb, pa[i]))
        .map(|i| a[i])
        .chain(
            (0..3)
                .filter(|&i| in_triangle_2d(k, s_a, pa, pb[i]))
                .map(|i| b[i]),
        )
        .collect();
    touch.sort_unstable();
    touch.dedup();
    if touch.is_empty() {
        return Meet::Empty;
    }
    let (lo, hi) = span(store, &touch);
    interval(store, lo, hi)
}

/// Where a point of a triangle's plane, known to lie in the closed triangle, sits.
enum Loc {
    Corner,
    Edge(usize),
    Inside,
}

fn locate(store: &Store, t: [u32; 3], k: usize, q: u32) -> Loc {
    if t.contains(&q) {
        return Loc::Corner;
    }
    for i in 0..3 {
        if orient2d(k, store.p(t[i]), store.p(t[(i + 1) % 3]), store.p(q)) == Sign::Zero {
            return Loc::Edge(i);
        }
    }
    Loc::Inside
}

fn key(a: u32, b: u32) -> (u32, u32) {
    (a.min(b), a.max(b))
}

/// A triangulation of one input triangle with points added and segments
/// enforced, in the plane of the triangle projected along axis `k`.
struct Cdt<'s> {
    store: &'s Store,
    k: usize,
    /// The projected orientation of the input triangle; every triangle here keeps it.
    s: Sign,
    points: Vec<u32>,
    triangles: Vec<[usize; 3]>,
    alive: Vec<bool>,
    /// Directed edge to the live triangle that has it.
    edges: HashMap<(usize, usize), usize>,
}

impl Cdt<'_> {
    fn orient(&self, a: usize, b: usize, c: usize) -> Sign {
        let p = |i: usize| self.store.p(self.points[i]);
        orient2d(self.k, p(a), p(b), p(c))
    }

    fn add(&mut self, t: [usize; 3]) {
        let index = self.triangles.len();
        self.triangles.push(t);
        self.alive.push(true);
        for i in 0..3 {
            self.edges.insert((t[i], t[(i + 1) % 3]), index);
        }
    }

    fn kill(&mut self, index: usize) {
        self.alive[index] = false;
        let t = self.triangles[index];
        for i in 0..3 {
            if self.edges.get(&(t[i], t[(i + 1) % 3])) == Some(&index) {
                self.edges.remove(&(t[i], t[(i + 1) % 3]));
            }
        }
    }

    /// The corner of triangle `index` after the directed edge `(x, y)`.
    fn third(&self, index: usize, x: usize, y: usize) -> usize {
        *self.triangles[index]
            .iter()
            .find(|&&v| v != x && v != y)
            .expect("a triangle has three corners")
    }

    /// Puts `p`, strictly inside the segment `ab`, into both triangles on it.
    fn split_edge(&mut self, a: usize, b: usize, p: usize) {
        for (x, y) in [(a, b), (b, a)] {
            if let Some(&t) = self.edges.get(&(x, y)) {
                let c = self.third(t, x, y);
                self.kill(t);
                self.add([x, p, c]);
                self.add([p, y, c]);
            }
        }
    }

    fn insert_inside(&mut self, p: usize, meter: &mut Meter) -> Result<(), Refusal> {
        for t in 0..self.triangles.len() {
            meter.tick(1)?;
            if !self.alive[t] {
                continue;
            }
            let [a, b, c] = self.triangles[t];
            let o = [
                self.orient(a, b, p),
                self.orient(b, c, p),
                self.orient(c, a, p),
            ];
            if o.contains(&flip(self.s)) {
                continue;
            }
            let zeros: Vec<usize> = (0..3).filter(|&i| o[i] == Sign::Zero).collect();
            match zeros.as_slice() {
                [] => {
                    self.kill(t);
                    self.add([a, b, p]);
                    self.add([b, c, p]);
                    self.add([c, a, p]);
                }
                &[i] => {
                    let corners = [a, b, c];
                    self.split_edge(corners[i], corners[(i + 1) % 3], p);
                }
                _ => return Err(internal("an inserted point repeats a corner")),
            }
            return Ok(());
        }
        Err(internal("an intersection point lies outside its triangle"))
    }

    /// Whether segments `uv` and `cd` cross at a point inside both.
    fn crosses(&self, u: usize, v: usize, c: usize, d: usize) -> bool {
        let (x, y) = (self.orient(u, v, c), self.orient(u, v, d));
        if x == Sign::Zero || y == Sign::Zero || x == y {
            return false;
        }
        let (x, y) = (self.orient(c, d, u), self.orient(c, d, v));
        x != Sign::Zero && y != Sign::Zero && x != y
    }

    fn between(&self, u: usize, v: usize, w: usize) -> bool {
        let (pu, pv, pw) = (self.points[u], self.points[v], self.points[w]);
        let (lo, hi) = if self.store.cmp(pu, pv) == Ordering::Less {
            (pu, pv)
        } else {
            (pv, pu)
        };
        self.store.cmp(lo, pw) == Ordering::Less && self.store.cmp(pw, hi) == Ordering::Less
    }

    /// Makes `uv` an edge. A point on the segment splits it first; then every
    /// edge crossing it is flipped, a non-convex one after the others, which
    /// ends because each flip removes a crossing (Sloan 1993).
    fn enforce(
        &mut self,
        u: usize,
        v: usize,
        done: &mut HashSet<(usize, usize)>,
        meter: &mut Meter,
    ) -> Result<(), Refusal> {
        meter.tick(self.points.len() as u64)?;
        let mut on: Vec<usize> = (0..self.points.len())
            .filter(|&w| {
                w != u && w != v && self.orient(u, v, w) == Sign::Zero && self.between(u, v, w)
            })
            .collect();
        if !on.is_empty() {
            let forward = self.store.cmp(self.points[u], self.points[v]) == Ordering::Less;
            on.sort_by(|&x, &y| {
                let o = self.store.cmp(self.points[x], self.points[y]);
                if forward { o } else { o.reverse() }
            });
            let chain: Vec<usize> = std::iter::once(u).chain(on).chain([v]).collect();
            for pair in chain.windows(2) {
                self.enforce(pair[0], pair[1], done, meter)?;
            }
            return Ok(());
        }
        let undirected = (u.min(v), u.max(v));
        if self.edges.contains_key(&(u, v)) || self.edges.contains_key(&(v, u)) {
            done.insert(undirected);
            return Ok(());
        }
        let mut queue = VecDeque::new();
        let mut seen = HashSet::new();
        for t in 0..self.triangles.len() {
            if !self.alive[t] {
                continue;
            }
            let tri = self.triangles[t];
            for i in 0..3 {
                let (a, b) = (tri[i], tri[(i + 1) % 3]);
                if self.crosses(u, v, a, b) && seen.insert((a.min(b), a.max(b))) {
                    queue.push_back((a, b));
                }
            }
        }
        meter.tick(self.triangles.len() as u64)?;
        while let Some((a, b)) = queue.pop_front() {
            meter.tick(1)?;
            if done.contains(&(a.min(b), a.max(b))) {
                return Err(internal(
                    "two intersection segments cross inside one triangle",
                ));
            }
            let (Some(&t1), Some(&t2)) = (self.edges.get(&(a, b)), self.edges.get(&(b, a))) else {
                return Err(internal("an intersection segment leaves its triangle"));
            };
            let (c, d) = (self.third(t1, a, b), self.third(t2, b, a));
            let (oa, ob) = (self.orient(c, d, a), self.orient(c, d, b));
            if oa != Sign::Zero && ob != Sign::Zero && oa != ob {
                self.kill(t1);
                self.kill(t2);
                self.add([a, d, c]);
                self.add([d, b, c]);
                if self.crosses(u, v, c, d) {
                    queue.push_back((c, d));
                }
            } else {
                queue.push_back((a, b));
            }
        }
        if !(self.edges.contains_key(&(u, v)) || self.edges.contains_key(&(v, u))) {
            return Err(internal("an intersection segment was not recovered"));
        }
        done.insert(undirected);
        Ok(())
    }
}

/// The pieces of one input triangle: `edge_points[i]` lie strictly inside the
/// edge from corner `i` to corner `i + 1`, in that order; `inside` lie in its
/// interior; `segments` must become edges. Returns the pieces, oriented as
/// the input triangle, and the segment pieces that became edges.
#[allow(clippy::type_complexity)]
fn pieces(
    store: &Store,
    tri: [u32; 3],
    edge_points: &[Vec<u32>; 3],
    inside: &[u32],
    segments: &[(u32, u32)],
    meter: &mut Meter,
) -> Result<(Vec<[u32; 3]>, Vec<(u32, u32)>), Refusal> {
    let corners = tri.map(|v| store.p(v));
    let (k, s) = projection(corners[0], corners[1], corners[2]).expect("validated triangles");
    let mut cdt = Cdt {
        store,
        k,
        s,
        points: tri.to_vec(),
        triangles: Vec::new(),
        alive: Vec::new(),
        edges: HashMap::new(),
    };
    let mut local: HashMap<u32, usize> = tri.iter().enumerate().map(|(i, &v)| (v, i)).collect();
    cdt.add([0, 1, 2]);
    for (i, list) in edge_points.iter().enumerate() {
        let end = (i + 1) % 3;
        let mut previous = i;
        for &p in list {
            let index = cdt.points.len();
            cdt.points.push(p);
            local.insert(p, index);
            cdt.split_edge(previous, end, index);
            previous = index;
        }
    }
    for &p in inside {
        if local.contains_key(&p) {
            continue;
        }
        let index = cdt.points.len();
        cdt.points.push(p);
        local.insert(p, index);
        cdt.insert_inside(index, meter)?;
    }
    let mut done = HashSet::new();
    for &(a, b) in segments {
        let (Some(&u), Some(&v)) = (local.get(&a), local.get(&b)) else {
            return Err(internal("a segment end is not a point of its triangle"));
        };
        cdt.enforce(u, v, &mut done, meter)?;
    }
    let triangles = (0..cdt.triangles.len())
        .filter(|&t| cdt.alive[t])
        .map(|t| cdt.triangles[t].map(|i| cdt.points[i]))
        .collect();
    let mut cuts: Vec<(u32, u32)> = done
        .into_iter()
        .map(|(a, b)| key(cdt.points[a], cdt.points[b]))
        .collect();
    cuts.sort_unstable();
    Ok((triangles, cuts))
}

/// The Boolean of two placed meshes (scale one, in millimetres). Work is
/// charged to `dec`.
pub fn boolean(
    target: &ClosedMesh,
    tool: &ClosedMesh,
    kind: Kind,
    dec: &mut Decider,
) -> Result<Outcome, Refusal> {
    if target.scale != Rational::one() || tool.scale != Rational::one() {
        return Err(internal("mesh Booleans take placed meshes"));
    }
    let operands = [target, tool];
    let mut meter = Meter::new(dec);
    let mut store = Store::default();
    let mut corners: [Vec<[u32; 3]>; 2] = [Vec::new(), Vec::new()];
    for (side, mesh) in operands.iter().enumerate() {
        meter.tick(mesh.exact.len() as u64)?;
        let ids: Vec<u32> = mesh.exact.iter().map(|p| store.add(p.clone())).collect();
        corners[side] = mesh
            .triangles
            .iter()
            .map(|t| t.map(|v| ids[v as usize]))
            .collect();
    }

    // Every pair of triangles that meet, found through the tool's box tree.
    let mut contacts: [Vec<Vec<(u32, u32)>>; 2] = [
        vec![Vec::new(); target.triangles.len()],
        vec![Vec::new(); tool.triangles.len()],
    ];
    for (a, tri) in target.triangles.iter().enumerate() {
        let query = outward_box(tri.map(|v| &target.approx[v as usize]));
        let mut candidates = Vec::new();
        tool.tree.query(&query, |b| {
            candidates.push(b);
            true
        });
        meter.tick(1 + candidates.len() as u64)?;
        candidates.sort_unstable();
        for b in candidates {
            meter.tick(1)?;
            let (ta, tb) = (corners[0][a], corners[1][b as usize]);
            let found = match meet(&mut store, ta, tb) {
                Meet::Empty => continue,
                Meet::Point(p) => (p, p),
                Meet::Segment(p, q) => (p, q),
                Meet::Overlap => {
                    return Err(refuse(
                        "mesh-boolean.coplanar",
                        format!(
                            "triangle {a} of the target and triangle {b} of the tool overlap in one \
                             plane near {}; coplanar contact is not in this first build",
                            store.describe(ta[0])
                        ),
                    ));
                }
            };
            contacts[0][a].push(found);
            contacts[1][b as usize].push(found);
        }
    }

    // Each operand's triangles, cut along the intersection.
    let mut all_pieces: [Vec<([u32; 3], u32)>; 2] = [Vec::new(), Vec::new()];
    let mut cuts: [HashSet<(u32, u32)>; 2] = [HashSet::new(), HashSet::new()];
    for side in 0..2 {
        // Intersection points on each edge, shared by both triangles on it.
        let mut on_edge: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
        for (t, list) in contacts[side].iter().enumerate() {
            if list.is_empty() {
                continue;
            }
            let tri = corners[side][t];
            let c = tri.map(|v| store.p(v));
            let (k, _) = projection(c[0], c[1], c[2]).expect("validated triangles");
            for &(p, q) in list {
                for x in [p, q] {
                    meter.tick(1)?;
                    if let Loc::Edge(i) = locate(&store, tri, k, x) {
                        on_edge
                            .entry(key(tri[i], tri[(i + 1) % 3]))
                            .or_default()
                            .push(x);
                    }
                }
            }
        }
        for list in on_edge.values_mut() {
            list.sort_unstable();
            list.dedup();
        }
        for t in 0..corners[side].len() {
            let tri = corners[side][t];
            let mut edge_points: [Vec<u32>; 3] = Default::default();
            for (i, points) in edge_points.iter_mut().enumerate() {
                let (from, to) = (tri[i], tri[(i + 1) % 3]);
                if let Some(list) = on_edge.get(&key(from, to)) {
                    let mut list = list.clone();
                    let forward = store.cmp(from, to) == Ordering::Less;
                    list.sort_by(|&x, &y| {
                        let o = store.cmp(x, y);
                        if forward { o } else { o.reverse() }
                    });
                    *points = list;
                }
            }
            let list = &contacts[side][t];
            if list.is_empty() && edge_points.iter().all(Vec::is_empty) {
                all_pieces[side].push((tri, t as u32));
                continue;
            }
            let c = tri.map(|v| store.p(v));
            let (k, _) = projection(c[0], c[1], c[2]).expect("validated triangles");
            let mut inside = Vec::new();
            let mut segments = Vec::new();
            for &(p, q) in list {
                for x in [p, q] {
                    if matches!(locate(&store, tri, k, x), Loc::Inside) {
                        inside.push(x);
                    }
                }
                if p != q {
                    segments.push(key(p, q));
                }
            }
            inside.sort_unstable();
            inside.dedup();
            segments.sort_unstable();
            segments.dedup();
            let (triangles, cut) =
                pieces(&store, tri, &edge_points, &inside, &segments, &mut meter)?;
            all_pieces[side].extend(triangles.into_iter().map(|p| (p, t as u32)));
            cuts[side].extend(cut);
        }
    }

    // Pieces joined by an edge that is not an intersection edge lie on the
    // same side of the other solid; one exact winding number decides each set.
    let mut keep: Vec<([u32; 3], (u8, u32))> = Vec::new();
    for side in 0..2 {
        let other = operands[1 - side];
        let pieces = &all_pieces[side];
        let mut by_edge: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
        for (i, (tri, _)) in pieces.iter().enumerate() {
            for j in 0..3 {
                by_edge
                    .entry(key(tri[j], tri[(j + 1) % 3]))
                    .or_default()
                    .push(i as u32);
            }
        }
        meter.tick(pieces.len() as u64)?;
        let mut patch = vec![u32::MAX; pieces.len()];
        let mut inside: Vec<bool> = Vec::new();
        for start in 0..pieces.len() {
            if patch[start] != u32::MAX {
                continue;
            }
            let id = inside.len() as u32;
            patch[start] = id;
            let mut stack = vec![start as u32];
            while let Some(i) = stack.pop() {
                meter.tick(1)?;
                let tri = pieces[i as usize].0;
                for j in 0..3 {
                    let edge = key(tri[j], tri[(j + 1) % 3]);
                    if cuts[side].contains(&edge) {
                        continue;
                    }
                    for &n in &by_edge[&edge] {
                        if patch[n as usize] == u32::MAX {
                            patch[n as usize] = id;
                            stack.push(n);
                        }
                    }
                }
            }
            let [a, b, c] = pieces[start].0.map(|v| &store.exact[v as usize]);
            record_exact(9, a.iter().chain(b).chain(c));
            let three = Rational::from_i64(3);
            let centroid: [Rational; 3] = std::array::from_fn(|i| {
                (&(&a[i] + &b[i]) + &c[i])
                    .checked_div(&three)
                    .expect("three is not zero")
            });
            meter.flush()?;
            inside.push(match other.winding(&centroid, meter.dec)? {
                Some(0) => false,
                Some(1) => true,
                found => {
                    return Err(internal(format!(
                        "a piece near {} has winding number {found:?} in the other mesh",
                        store.describe(pieces[start].0[0])
                    )));
                }
            });
        }
        for (i, &(tri, source)) in pieces.iter().enumerate() {
            let within = inside[patch[i] as usize];
            let reversed = match (kind, side, within) {
                (Kind::Union, _, false)
                | (Kind::Intersect, _, true)
                | (Kind::Subtract, 0, false) => false,
                (Kind::Subtract, 1, true) => true,
                _ => continue,
            };
            let tri = if reversed {
                [tri[0], tri[2], tri[1]]
            } else {
                tri
            };
            keep.push((tri, (side as u8, source)));
        }
    }

    // The kept pieces, on their own vertices, proved closed and non-crossing.
    let mut index: HashMap<u32, u32> = HashMap::new();
    let mut points = Vec::new();
    let mut triangles = Vec::with_capacity(keep.len());
    let mut sources = Vec::with_capacity(keep.len());
    for (tri, source) in keep {
        triangles.push(tri.map(|v| {
            *index.entry(v).or_insert_with(|| {
                points.push(store.exact[v as usize].clone());
                points.len() as u32 - 1
            })
        }));
        sources.push(source);
    }
    meter.flush()?;
    let mesh = ClosedMesh::from_parts(points, triangles, meter.dec).map_err(|refusal| {
        match refusal.code {
            "import.mesh-nonmanifold-edge" | "import.mesh-pinched-vertex" => refuse(
                "mesh-boolean.non-manifold",
                format!(
                    "the meshes touch along an edge or at a point, so the result is not a manifold solid: {}",
                    refusal.detail
                ),
            ),
            code if code.starts_with("import.mesh-") => internal(format!(
                "the result failed the closed-mesh check: {code}: {}",
                refusal.detail
            )),
            _ => refusal,
        }
    })?;
    Ok(Outcome { mesh, sources })
}
