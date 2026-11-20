//! Decimation of a closed triangle mesh, RECIPE.md section 4.22.
//!
//! The method is SMLib's `SmPolyDecimate::DoQuadricDecimation`, which is the
//! quadric error metric of Garland and Heckbert ("Surface simplification using
//! quadric error metrics", SIGGRAPH 1997). Each vertex carries the sum of the
//! squared-distance quadrics of the planes of its triangles, each weighted by
//! a third of the triangle's area as SMLib weights them, and edges collapse
//! cheapest first. Three choices differ from SMLib, each for a reason:
//!
//! - The surviving vertex goes to whichever of the two endpoints and their
//!   midpoint has the least quadric error. That is Garland and Heckbert's own
//!   fallback set. Their optimal point `-Q^-1 b` would need a rank test with a
//!   threshold and is a fitted binary64 point; endpoints and midpoints are
//!   exact rationals, so the result stays an exact mesh like its input.
//! - A collapse must keep a 2-manifold: the link condition of Dey,
//!   Edelsbrunner, Guha and Nekhayev ("Topology preserving edge contraction",
//!   1999). SMLib checks no topology.
//! - A collapse may not turn a triangle over or flatten it: each moved
//!   triangle's exact normal keeps a positive dot product with its old normal,
//!   the mesh-inversion guard Garland and Heckbert describe.
//!
//! The result is then proved closed, oriented and free of crossings, exactly,
//! as an imported mesh is.

use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeSet, BinaryHeap};

use catallax_numbers::{Decider, Rational, Refusal, Sign};

use super::predicates::record_exact;
use super::{ClosedMesh, Meter};

/// The validated result and where each of its triangles came from.
#[derive(Clone, Debug)]
pub struct Outcome {
    /// The result, in millimetres.
    pub mesh: ClosedMesh,
    /// For each result triangle, the input triangle it was; a collapse moves
    /// corners of the triangles it keeps and creates none.
    pub sources: Vec<u32>,
}

fn refuse(code: &'static str, detail: impl Into<String>) -> Refusal {
    Refusal::new(code, detail)
}

type V = [f64; 3];

fn sub(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] - b[i])
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The quadric `x^T A x + 2 b.x + c` of Garland and Heckbert.
#[derive(Clone, Copy, Default)]
struct Quadric {
    /// The symmetric `A` as xx, xy, xz, yy, yz, zz.
    a: [f64; 6],
    b: V,
    c: f64,
}

impl Quadric {
    /// The squared distance to the plane `n.x + d = 0`, `n` a unit normal,
    /// times `weight`.
    fn plane(n: V, d: f64, weight: f64) -> Self {
        Self {
            a: [
                n[0] * n[0],
                n[0] * n[1],
                n[0] * n[2],
                n[1] * n[1],
                n[1] * n[2],
                n[2] * n[2],
            ]
            .map(|x| x * weight),
            b: n.map(|x| x * d * weight),
            c: d * d * weight,
        }
    }

    fn plus(&self, o: &Self) -> Self {
        Self {
            a: std::array::from_fn(|i| self.a[i] + o.a[i]),
            b: std::array::from_fn(|i| self.b[i] + o.b[i]),
            c: self.c + o.c,
        }
    }

    fn error(&self, x: V) -> f64 {
        let [xx, xy, xz, yy, yz, zz] = self.a;
        let ax = [
            xx * x[0] + xy * x[1] + xz * x[2],
            xy * x[0] + yy * x[1] + yz * x[2],
            xz * x[0] + yz * x[1] + zz * x[2],
        ];
        dot(x, ax) + 2.0 * dot(self.b, x) + self.c
    }
}

/// A candidate collapse, valid while both vertices keep their stamps.
#[derive(Clone, Copy, PartialEq)]
struct Candidate {
    cost: f64,
    keep: u32,
    gone: u32,
    stamps: (u32, u32),
    /// 0: the kept vertex's place, 1: the removed vertex's, 2: their midpoint.
    place: u8,
}

impl Eq for Candidate {}
impl PartialOrd for Candidate {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Candidate {
    fn cmp(&self, other: &Self) -> Ordering {
        self.cost
            .total_cmp(&other.cost)
            .then(self.keep.cmp(&other.keep))
            .then(self.gone.cmp(&other.gone))
            .then(self.stamps.cmp(&other.stamps))
            .then(self.place.cmp(&other.place))
    }
}

struct State<'m> {
    exact: Vec<[Rational; 3]>,
    approx: Vec<V>,
    triangles: Vec<[u32; 3]>,
    alive: Vec<bool>,
    /// Live triangles around each vertex.
    around: Vec<Vec<u32>>,
    quadrics: Vec<Quadric>,
    stamps: Vec<u32>,
    meter: Meter<'m>,
}

impl State<'_> {
    fn neighbours(&self, v: u32) -> BTreeSet<u32> {
        self.around[v as usize]
            .iter()
            .flat_map(|&t| self.triangles[t as usize])
            .filter(|&w| w != v)
            .collect()
    }

    fn midpoint(&self, u: u32, v: u32) -> [Rational; 3] {
        let (a, b) = (&self.exact[u as usize], &self.exact[v as usize]);
        record_exact(6, a.iter().chain(b));
        let two = Rational::from_i64(2);
        std::array::from_fn(|i| (&a[i] + &b[i]).checked_div(&two).expect("two is not zero"))
    }

    fn place(&self, keep: u32, gone: u32, place: u8) -> [Rational; 3] {
        match place {
            0 => self.exact[keep as usize].clone(),
            1 => self.exact[gone as usize].clone(),
            _ => self.midpoint(keep, gone),
        }
    }

    /// The cheapest placement for collapsing `gone` into `keep`.
    fn candidate(&self, keep: u32, gone: u32) -> Candidate {
        let q = self.quadrics[keep as usize].plus(&self.quadrics[gone as usize]);
        let (a, b) = (self.approx[keep as usize], self.approx[gone as usize]);
        // Choosing a placement is a binary64 decision; only the chosen place
        // is made exactly.
        let middle = std::array::from_fn(|i| a[i] * 0.5 + b[i] * 0.5);
        let mut best = (q.error(a), 0u8);
        for (place, x) in [(1u8, b), (2u8, middle)] {
            let e = q.error(x);
            if e < best.0 {
                best = (e, place);
            }
        }
        Candidate {
            cost: best.0,
            keep,
            gone,
            stamps: (self.stamps[keep as usize], self.stamps[gone as usize]),
            place: best.1,
        }
    }

    /// The exact normal `(b - a) x (c - a)`.
    fn normal(points: [&[Rational; 3]; 3]) -> [Rational; 3] {
        let [a, b, c] = points;
        record_exact(9, a.iter().chain(b).chain(c));
        let u: [Rational; 3] = std::array::from_fn(|i| &b[i] - &a[i]);
        let v: [Rational; 3] = std::array::from_fn(|i| &c[i] - &a[i]);
        [
            &(&u[1] * &v[2]) - &(&u[2] * &v[1]),
            &(&u[2] * &v[0]) - &(&u[0] * &v[2]),
            &(&u[0] * &v[1]) - &(&u[1] * &v[0]),
        ]
    }

    /// Whether collapsing `gone` into `keep` at `x` keeps a 2-manifold and
    /// turns no triangle over.
    fn allowed(&mut self, keep: u32, gone: u32, x: &[Rational; 3]) -> Result<bool, Refusal> {
        // The link condition: the common neighbours are exactly the two
        // corners opposite the edge, and no triangle joins those two corners
        // to either end, which would close a tetrahedron.
        let shared: Vec<u32> = self.around[keep as usize]
            .iter()
            .copied()
            .filter(|&t| self.triangles[t as usize].contains(&gone))
            .collect();
        if shared.len() != 2 {
            return Ok(false);
        }
        let opposite: BTreeSet<u32> = shared
            .iter()
            .flat_map(|&t| self.triangles[t as usize])
            .filter(|&w| w != keep && w != gone)
            .collect();
        let common: BTreeSet<u32> = self
            .neighbours(keep)
            .intersection(&self.neighbours(gone))
            .copied()
            .collect();
        self.meter.tick(common.len() as u64 + 1)?;
        if common != opposite || opposite.len() != 2 {
            return Ok(false);
        }
        let ends: Vec<u32> = opposite.iter().copied().collect();
        for v in [keep, gone] {
            if self.around[v as usize].iter().any(|&t| {
                let tri = self.triangles[t as usize];
                tri.contains(&ends[0]) && tri.contains(&ends[1])
            }) {
                return Ok(false);
            }
        }
        // No moved triangle turns over or becomes flat.
        let moved: Vec<u32> = self.around[keep as usize]
            .iter()
            .chain(&self.around[gone as usize])
            .copied()
            .filter(|t| !shared.contains(t))
            .collect();
        for t in moved {
            self.meter.tick(1)?;
            let tri = self.triangles[t as usize];
            let old = Self::normal(tri.map(|v| &self.exact[v as usize]));
            let new = Self::normal(tri.map(|v| {
                if v == keep || v == gone {
                    x
                } else {
                    &self.exact[v as usize]
                }
            }));
            record_exact(5, old.iter().chain(&new));
            let d = (0..3).fold(Rational::zero(), |sum, i| &sum + &(&old[i] * &new[i]));
            if d.sign() != Sign::Positive {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn collapse(&mut self, keep: u32, gone: u32, x: [Rational; 3]) {
        for &t in &self.around[gone as usize].clone() {
            if self.triangles[t as usize].contains(&keep) {
                self.alive[t as usize] = false;
                for v in self.triangles[t as usize] {
                    self.around[v as usize].retain(|&s| s != t);
                }
            } else {
                for v in &mut self.triangles[t as usize] {
                    if *v == gone {
                        *v = keep;
                    }
                }
                self.around[keep as usize].push(t);
            }
        }
        self.around[gone as usize].clear();
        self.around[keep as usize].sort_unstable();
        self.approx[keep as usize] = x.each_ref().map(Rational::to_f64_nearest);
        self.exact[keep as usize] = x;
        self.quadrics[keep as usize] =
            self.quadrics[keep as usize].plus(&self.quadrics[gone as usize]);
        self.stamps[keep as usize] += 1;
        self.stamps[gone as usize] += 1;
    }
}

/// Collapses edges of `mesh` until at most `target` triangles remain or no
/// allowed collapse is left. Work is charged to `dec`.
pub fn decimate(mesh: &ClosedMesh, target: usize, dec: &mut Decider) -> Result<Outcome, Refusal> {
    let exact: Vec<[Rational; 3]> = (0..mesh.vertex_count() as u32)
        .map(|v| mesh.vertex(v))
        .collect();
    let approx: Vec<V> = exact
        .iter()
        .map(|p| p.each_ref().map(Rational::to_f64_nearest))
        .collect();
    let mut state = State {
        around: vec![Vec::new(); exact.len()],
        quadrics: vec![Quadric::default(); exact.len()],
        stamps: vec![0; exact.len()],
        triangles: mesh.triangles.clone(),
        alive: vec![true; mesh.triangles.len()],
        exact,
        approx,
        meter: Meter::new(dec),
    };
    state.meter.tick(mesh.triangles.len() as u64 + 1)?;
    for (t, tri) in mesh.triangles.iter().enumerate() {
        let [a, b, c] = tri.map(|v| state.approx[v as usize]);
        let n = cross(sub(b, a), sub(c, a));
        let length = dot(n, n).sqrt();
        // A validated triangle has positive area; a binary64 normal can still
        // round to zero on a sliver, and then adds no plane.
        if length > 0.0 && length.is_finite() {
            let unit = n.map(|x| x / length);
            let weight = length / 2.0 / 3.0;
            let plane = Quadric::plane(unit, -dot(unit, a), weight);
            for &v in tri {
                state.quadrics[v as usize] = state.quadrics[v as usize].plus(&plane);
            }
        }
        for &v in tri {
            state.around[v as usize].push(t as u32);
        }
    }
    let mut heap = BinaryHeap::new();
    for tri in &mesh.triangles {
        for i in 0..3 {
            let (u, v) = (tri[i], tri[(i + 1) % 3]);
            // Each edge once, from the triangle that runs it upward.
            if u < v {
                heap.push(Reverse(state.candidate(u, v)));
            }
        }
    }
    state.meter.tick(heap.len() as u64)?;
    let mut count = mesh.triangles.len();
    while count > target {
        let Some(Reverse(c)) = heap.pop() else {
            break;
        };
        state.meter.tick(1)?;
        if c.stamps != (state.stamps[c.keep as usize], state.stamps[c.gone as usize]) {
            continue;
        }
        let x = state.place(c.keep, c.gone, c.place);
        if !state.allowed(c.keep, c.gone, &x)? {
            continue;
        }
        state.collapse(c.keep, c.gone, x);
        count -= 2;
        for w in state.neighbours(c.keep) {
            state.meter.tick(1)?;
            heap.push(Reverse(state.candidate(c.keep, w)));
        }
    }

    let mut index = vec![u32::MAX; state.exact.len()];
    let mut points = Vec::new();
    let mut triangles = Vec::with_capacity(count);
    let mut sources = Vec::with_capacity(count);
    for t in 0..state.triangles.len() {
        if !state.alive[t] {
            continue;
        }
        triangles.push(state.triangles[t].map(|v| {
            if index[v as usize] == u32::MAX {
                index[v as usize] = points.len() as u32;
                points.push(state.exact[v as usize].clone());
            }
            index[v as usize]
        }));
        sources.push(t as u32);
    }
    let State { mut meter, .. } = state;
    meter.flush()?;
    let mesh =
        ClosedMesh::from_parts(points, triangles, &mut *meter.dec).map_err(
            |refusal| match refusal.code {
                "import.mesh-crossing" => refuse(
                    "mesh-decimate.crossing",
                    format!(
                        "the decimated surface would cross itself, so it bounds no solid: {}",
                        refusal.detail
                    ),
                ),
                code if code.starts_with("import.mesh-") => refuse(
                    "mesh-decimate.internal",
                    format!(
                        "the result failed the closed-mesh check: {code}: {}",
                        refusal.detail
                    ),
                ),
                _ => refusal,
            },
        )?;
    Ok(Outcome { mesh, sources })
}
