//! Closed triangle mesh import, RECIPE.md section 4.19.
//!
//! [`ClosedMesh::read`] reads binary STL, ASCII STL or Wavefront OBJ bytes and
//! proves, exactly and with no tolerance, that the triangles bound a solid:
//! every triangle has positive area, every edge is used once in each direction
//! by exactly two triangles, every vertex has one fan of triangles, no two
//! triangles meet except in their shared edge or vertex, and the integer
//! winding number is one inside and zero outside. Inside and outside are then
//! decided by counting windings along a symbolically perturbed ray, never by a
//! threshold.
//!
//! All predicates work on the file's own coordinates. A positive scale does
//! not change any sign, so scaling is applied only to the points handed out.

mod bvh;
mod parse;
mod predicates;

use catallax_numbers::{Decider, Interval, Rational, Refusal, Sign};

use bvh::Bvh;
use predicates::{P, in_corner_2d, in_triangle_2d, orient2d, orient3d, projection};

/// The stated size limit: files with more triangles are refused with
/// `import.mesh-too-large` before any geometry is built. The B-rep bounds it:
/// at this size the B-rep build takes a few seconds and a few gigabytes, while
/// the field alone handles five times more. RECIPE.md section 4.19 records the
/// measurements.
pub const MAX_TRIANGLES: usize = 200_000;

/// One connected closed surface of the mesh.
#[derive(Clone, Debug)]
pub struct MeshShell {
    /// Its triangles, in file order.
    pub triangles: Vec<u32>,
    /// Whether its triangles face away from the region it encloses. An outer
    /// boundary faces outward; a cavity's triangles face into the cavity.
    pub outward: bool,
    /// The innermost other shell that encloses this one.
    pub parent: Option<usize>,
}

/// A validated closed, oriented, non-crossing triangle mesh.
#[derive(Clone)]
pub struct ClosedMesh {
    exact: Vec<[Rational; 3]>,
    approx: Vec<[f64; 3]>,
    /// Whether each vertex's binary64 copy is close enough for the filter.
    close: Vec<bool>,
    scale: Rational,
    triangles: Vec<[u32; 3]>,
    twins: Vec<[u32; 3]>,
    shells: Vec<MeshShell>,
    tree: Bvh,
}

impl std::fmt::Debug for ClosedMesh {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ClosedMesh")
            .field("vertices", &self.exact.len())
            .field("triangles", &self.triangles.len())
            .field("shells", &self.shells.len())
            .field("scale", &self.scale)
            .finish()
    }
}

fn refuse(code: &'static str, detail: impl Into<String>) -> Refusal {
    Refusal::new(code, detail)
}

fn outward_box(points: [&[f64; 3]; 3]) -> ([f64; 3], [f64; 3]) {
    let lo = std::array::from_fn(|i| {
        points
            .iter()
            .map(|p| p[i])
            .fold(f64::INFINITY, f64::min)
            .next_down()
    });
    let hi = std::array::from_fn(|i| {
        points
            .iter()
            .map(|p| p[i])
            .fold(f64::NEG_INFINITY, f64::max)
            .next_up()
    });
    (lo, hi)
}

/// Exact squared distance from `p` to the closed triangle `abc`.
fn triangle_distance_sq(
    p: &[Rational; 3],
    a: &[Rational; 3],
    b: &[Rational; 3],
    c: &[Rational; 3],
) -> Rational {
    let sub = |x: &[Rational; 3], y: &[Rational; 3]| -> [Rational; 3] {
        std::array::from_fn(|i| &x[i] - &y[i])
    };
    let dot = |x: &[Rational; 3], y: &[Rational; 3]| -> Rational {
        &(&(&x[0] * &y[0]) + &(&x[1] * &y[1])) + &(&x[2] * &y[2])
    };
    let zero = Rational::zero();
    let length = |v: [Rational; 3]| dot(&v, &v);
    let along = |origin: &[Rational; 3], v: &[Rational; 3], t: &Rational| -> [Rational; 3] {
        std::array::from_fn(|i| &origin[i] + &(&v[i] * t))
    };
    let (ab, ac, ap) = (sub(b, a), sub(c, a), sub(p, a));
    let (d1, d2) = (dot(&ab, &ap), dot(&ac, &ap));
    if d1 <= zero && d2 <= zero {
        return length(ap);
    }
    let bp = sub(p, b);
    let (d3, d4) = (dot(&ab, &bp), dot(&ac, &bp));
    if d3 >= zero && d4 <= d3 {
        return length(bp);
    }
    let vc = &(&d1 * &d4) - &(&d3 * &d2);
    if vc <= zero && d1 >= zero && d3 <= zero {
        let t = d1.checked_div(&(&d1 - &d3)).expect("positive edge length");
        return length(sub(p, &along(a, &ab, &t)));
    }
    let cp = sub(p, c);
    let (d5, d6) = (dot(&ab, &cp), dot(&ac, &cp));
    if d6 >= zero && d5 <= d6 {
        return length(cp);
    }
    let vb = &(&d5 * &d2) - &(&d1 * &d6);
    if vb <= zero && d2 >= zero && d6 <= zero {
        let t = d2.checked_div(&(&d2 - &d6)).expect("positive edge length");
        return length(sub(p, &along(a, &ac, &t)));
    }
    let va = &(&d3 * &d6) - &(&d5 * &d4);
    let (e1, e2) = (&d4 - &d3, &d5 - &d6);
    if va <= zero && e1 >= zero && e2 >= zero {
        let t = e1.checked_div(&(&e1 + &e2)).expect("positive edge length");
        return length(sub(p, &along(b, &sub(c, b), &t)));
    }
    let total = &(&va + &vb) + &vc;
    let v = vb.checked_div(&total).expect("positive triangle area");
    let w = vc.checked_div(&total).expect("positive triangle area");
    let foot: [Rational; 3] = std::array::from_fn(|i| &(&a[i] + &(&ab[i] * &v)) + &(&ac[i] * &w));
    length(sub(p, &foot))
}

/// Approximate squared distance from `p` to the triangle `abc`, for display.
fn fast_triangle_distance_sq(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let sub = |x: [f64; 3], y: [f64; 3]| -> [f64; 3] { std::array::from_fn(|i| x[i] - y[i]) };
    let dot = |x: [f64; 3], y: [f64; 3]| x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    let along =
        |o: [f64; 3], v: [f64; 3], t: f64| -> [f64; 3] { std::array::from_fn(|i| o[i] + v[i] * t) };
    let length = |v: [f64; 3]| dot(v, v);
    let (ab, ac, ap) = (sub(b, a), sub(c, a), sub(p, a));
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return length(ap);
    }
    let bp = sub(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 {
        return length(bp);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return length(sub(p, along(a, ab, d1 / (d1 - d3))));
    }
    let cp = sub(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return length(cp);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return length(sub(p, along(a, ac, d2 / (d2 - d6))));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
        return length(sub(
            p,
            along(b, sub(c, b), (d4 - d3) / ((d4 - d3) + (d5 - d6))),
        ));
    }
    let total = va + vb + vc;
    let (v, w) = (vb / total, vc / total);
    length(sub(
        p,
        std::array::from_fn(|i| a[i] + ab[i] * v + ac[i] * w),
    ))
}

/// A rigorous lower bound on squared distance from the box to a point that
/// lies in the per-axis intervals `[lo, hi]`.
fn box_lower_bound(b: &([f64; 3], [f64; 3]), lo: &[f64; 3], hi: &[f64; 3]) -> f64 {
    let mut sum = 0.0f64;
    for i in 0..3 {
        let gap = (b.0[i] - hi[i]).max(lo[i] - b.1[i]);
        if gap > 0.0 {
            let gap = gap.next_down().max(0.0);
            sum = (sum + (gap * gap).next_down()).next_down();
        }
    }
    sum.max(0.0)
}

fn fast_box_distance(b: &([f64; 3], [f64; 3]), p: &[f64; 3]) -> f64 {
    (0..3)
        .map(|i| {
            let gap = (b.0[i] - p[i]).max(p[i] - b.1[i]).max(0.0);
            gap * gap
        })
        .sum()
}

/// Sign of the projected orientation of `ab` and `q` with the axis-0 drop,
/// after moving `q` by `(0, e, e^2)` for an infinitesimal `e > 0`.
fn perturbed_yz(a: P, b: P, q: P) -> Sign {
    match orient2d(0, a, b, q) {
        Sign::Zero => {}
        s => return s,
    }
    match (&a.e[2] - &b.e[2]).sign() {
        Sign::Zero => (&b.e[1] - &a.e[1]).sign(),
        s => s,
    }
}

/// The result of the boundary-or-crossing walk along one ray.
enum Ray {
    Boundary,
    Clear,
}

/// Charges the caller's allowance for mesh work: one unit per triangle step,
/// candidate pair or visited triangle, plus the exact rational work the
/// predicates record. Charges go out in batches, and each batch also checks
/// the deadline.
struct Meter<'d> {
    dec: &'d mut Decider,
    pending: u64,
}

impl<'d> Meter<'d> {
    const BATCH: u64 = 4096;

    fn new(dec: &'d mut Decider) -> Self {
        Self { dec, pending: 0 }
    }

    fn tick(&mut self, units: u64) -> Result<(), Refusal> {
        self.pending = self.pending.saturating_add(units);
        if self
            .pending
            .saturating_add(predicates::pending_exact_work())
            >= Self::BATCH
        {
            self.flush()
        } else {
            Ok(())
        }
    }

    fn flush(&mut self) -> Result<(), Refusal> {
        let units = std::mem::take(&mut self.pending).saturating_add(predicates::take_exact_work());
        self.dec.charge_work(units)
    }
}

impl ClosedMesh {
    /// Reads and validates a mesh under the stated size limit. All work is
    /// charged to `dec`'s operation allowance and checked against its deadline.
    pub fn read(bytes: &[u8], scale: &Rational, dec: &mut Decider) -> Result<Self, Refusal> {
        Self::read_with_limit(bytes, scale, MAX_TRIANGLES, dec)
    }

    /// Reads and validates a mesh under a smaller limit, for callers with a
    /// tighter allowance. A limit above [`MAX_TRIANGLES`] is lowered to it.
    pub fn read_with_limit(
        bytes: &[u8],
        scale: &Rational,
        limit: usize,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if scale.sign() != Sign::Positive {
            return Err(refuse(
                "import.mesh-scale",
                format!("scale {scale} is not positive; a mesh import needs scale > 0"),
            ));
        }
        let mut meter = Meter::new(dec);
        // Reading is linear in the bytes; it is paid for before it starts.
        meter.tick(bytes.len() as u64 / 16 + 1)?;
        meter.flush()?;
        let soup = parse::read(bytes, limit.min(MAX_TRIANGLES))?;
        let mesh = Self::validate(soup.points, soup.triangles, scale.clone(), &mut meter)?;
        meter.flush()?;
        Ok(mesh)
    }

    fn validate(
        exact: Vec<[Rational; 3]>,
        triangles: Vec<[u32; 3]>,
        scale: Rational,
        meter: &mut Meter,
    ) -> Result<Self, Refusal> {
        let approx: Vec<[f64; 3]> = exact
            .iter()
            .map(|p| p.each_ref().map(Rational::to_f64_nearest))
            .collect();
        let close = exact
            .iter()
            .zip(&approx)
            .map(|(e, f)| predicates::close(e, f))
            .collect();
        let boxes = triangles
            .iter()
            .map(|t| outward_box(t.map(|v| &approx[v as usize])))
            .collect();
        let mut mesh = Self {
            exact,
            approx,
            close,
            scale,
            triangles,
            twins: Vec::new(),
            shells: Vec::new(),
            tree: Bvh::default(),
        };
        mesh.check_degenerate(meter)?;
        mesh.check_edges(meter)?;
        mesh.check_vertices(meter)?;
        meter.tick(mesh.triangles.len() as u64)?;
        mesh.tree = Bvh::new(boxes);
        mesh.check_crossings(meter)?;
        meter.tick(mesh.triangles.len() as u64)?;
        mesh.find_shells();
        mesh.check_windings(meter)?;
        Ok(mesh)
    }

    fn p(&self, v: u32) -> P<'_> {
        P {
            e: &self.exact[v as usize],
            f: self.approx[v as usize],
            close: self.close[v as usize],
        }
    }

    fn corners(&self, t: u32) -> [P<'_>; 3] {
        self.triangles[t as usize].map(|v| self.p(v))
    }

    fn describe(&self, v: u32) -> String {
        let [x, y, z] = self.approx[v as usize];
        format!("({x}, {y}, {z})")
    }

    fn check_degenerate(&self, meter: &mut Meter) -> Result<(), Refusal> {
        for (t, &[a, b, c]) in self.triangles.iter().enumerate() {
            meter.tick(1)?;
            if a == b || b == c || c == a || projection(self.p(a), self.p(b), self.p(c)).is_none() {
                return Err(refuse(
                    "import.mesh-degenerate",
                    format!(
                        "triangle {t} at {} has zero area, so it has no orientation",
                        self.describe(a)
                    ),
                ));
            }
        }
        Ok(())
    }

    fn check_edges(&mut self, meter: &mut Meter) -> Result<(), Refusal> {
        // Each directed edge as (undirected key, reversed flag, triangle, slot).
        let mut uses: Vec<(u64, bool, u32, u8)> = Vec::with_capacity(self.triangles.len() * 3);
        for (t, tri) in self.triangles.iter().enumerate() {
            for slot in 0..3 {
                let (a, b) = (tri[slot], tri[(slot + 1) % 3]);
                let key = (u64::from(a.min(b)) << 32) | u64::from(a.max(b));
                uses.push((key, a > b, t as u32, slot as u8));
            }
        }
        // One unit per directed edge covers sorting them; a comparison of two
        // integer keys is far cheaper than one exact rational operation.
        meter.tick(uses.len() as u64)?;
        uses.sort_unstable();
        self.twins = vec![[u32::MAX; 3]; self.triangles.len()];
        for group in uses.chunk_by(|x, y| x.0 == y.0) {
            meter.tick(1)?;
            let forward = group.iter().filter(|u| !u.1).count();
            let backward = group.len() - forward;
            let (u, v) = ((group[0].0 >> 32) as u32, group[0].0 as u32);
            let edge = format!("the edge from {} to {}", self.describe(u), self.describe(v));
            if forward != backward {
                return Err(if group.len() == 1 {
                    refuse(
                        "import.mesh-open",
                        format!(
                            "{edge} bounds only triangle {}; the mesh is not closed",
                            group[0].2
                        ),
                    )
                } else {
                    refuse(
                        "import.mesh-orientation",
                        format!(
                            "{edge} is used {forward} times one way and {backward} times the other; \
                             the triangles are not consistently oriented"
                        ),
                    )
                });
            }
            if group.len() != 2 {
                return Err(refuse(
                    "import.mesh-nonmanifold-edge",
                    format!(
                        "{edge} is shared by {} triangles; only two per edge are supported",
                        group.len()
                    ),
                ));
            }
            let (first, second) = (group[0], group[1]);
            self.twins[first.2 as usize][first.3 as usize] = second.2;
            self.twins[second.2 as usize][second.3 as usize] = first.2;
        }
        Ok(())
    }

    /// The slot of vertex `v` in triangle `t`.
    fn slot(&self, t: u32, v: u32) -> usize {
        self.triangles[t as usize]
            .iter()
            .position(|&x| x == v)
            .expect("the vertex belongs to the triangle")
    }

    fn check_vertices(&self, meter: &mut Meter) -> Result<(), Refusal> {
        let mut degree = vec![0u32; self.exact.len()];
        let mut first = vec![u32::MAX; self.exact.len()];
        for (t, tri) in self.triangles.iter().enumerate() {
            for &v in tri {
                degree[v as usize] += 1;
                if first[v as usize] == u32::MAX {
                    first[v as usize] = t as u32;
                }
            }
        }
        for v in 0..self.exact.len() as u32 {
            let start = first[v as usize];
            if start == u32::MAX {
                continue;
            }
            // Walk the fan: leave v along the edge after it, cross to the twin.
            let mut t = start;
            let mut steps = 0u32;
            loop {
                meter.tick(1)?;
                let slot = self.slot(t, v);
                t = self.twins[t as usize][slot];
                steps += 1;
                if t == start || steps > degree[v as usize] {
                    break;
                }
            }
            if steps != degree[v as usize] {
                return Err(refuse(
                    "import.mesh-pinched-vertex",
                    format!(
                        "vertex {} joins more than one fan of triangles; such pinched vertices are not supported",
                        self.describe(v)
                    ),
                ));
            }
        }
        Ok(())
    }

    /// Whether two distinct triangles meet anywhere other than in what they share.
    fn cross(&self, s: u32, t: u32) -> bool {
        let (ts, tt) = (self.triangles[s as usize], self.triangles[t as usize]);
        let shared: Vec<u32> = ts.iter().copied().filter(|v| tt.contains(v)).collect();
        let (ps, pt) = (self.corners(s), self.corners(t));
        let sides = |plane: [P; 3], points: [P; 3]| {
            points.map(|q| orient3d(plane[0], plane[1], plane[2], q))
        };
        match shared.len() {
            0 => {
                let a = sides(ps, pt);
                if a.iter().all(|&x| x == a[0]) && a[0] != Sign::Zero {
                    return false;
                }
                let b = sides(pt, ps);
                if b.iter().all(|&x| x == b[0]) && b[0] != Sign::Zero {
                    return false;
                }
                (0..3).any(|i| predicates::segment_meets_triangle(ps[i], ps[(i + 1) % 3], pt))
                    || (0..3)
                        .any(|i| predicates::segment_meets_triangle(pt[i], pt[(i + 1) % 3], ps))
            }
            1 => {
                let v = shared[0];
                let rotate = |tri: [u32; 3]| {
                    let i = tri.iter().position(|&x| x == v).expect("shared vertex");
                    [tri[i], tri[(i + 1) % 3], tri[(i + 2) % 3]].map(|x| self.p(x))
                };
                let (a, b) = (rotate(ts), rotate(tt));
                let one_sided = |plane: [P; 3], q: [P; 3]| {
                    let x = orient3d(plane[0], plane[1], plane[2], q[1]);
                    let y = orient3d(plane[0], plane[1], plane[2], q[2]);
                    x == y && x != Sign::Zero
                };
                if one_sided(a, b) && one_sided(b, a) {
                    return false;
                }
                let opposite = predicates::segment_meets_triangle(a[1], a[2], b)
                    || predicates::segment_meets_triangle(b[1], b[2], a);
                if opposite {
                    return true;
                }
                // An edge from the shared vertex enters the other triangle only
                // if it lies in that triangle's plane, inside its corner.
                let enters = |tri: [P; 3], edge_end: P| {
                    orient3d(tri[0], tri[1], tri[2], edge_end) == Sign::Zero && {
                        let (k, s) = projection(tri[0], tri[1], tri[2]).expect("nondegenerate");
                        in_corner_2d(k, s, tri, edge_end)
                    }
                };
                enters(b, a[1]) || enters(b, a[2]) || enters(a, b[1]) || enters(a, b[2])
            }
            2 => {
                let other = |tri: [u32; 3]| {
                    *tri.iter()
                        .find(|v| !shared.contains(v))
                        .expect("third vertex")
                };
                let (u, v) = (self.p(shared[0]), self.p(shared[1]));
                let (a, b) = (self.p(other(ts)), self.p(other(tt)));
                if orient3d(u, v, a, b) != Sign::Zero {
                    return false;
                }
                let (k, _) = projection(u, v, a).expect("nondegenerate");
                orient2d(k, u, v, a) == orient2d(k, u, v, b)
            }
            _ => true,
        }
    }

    fn check_crossings(&self, meter: &mut Meter) -> Result<(), Refusal> {
        let found = self.tree.pairs(|s, t| match meter.tick(1) {
            Err(refusal) => Some(Err(refusal)),
            Ok(()) => self.cross(s, t).then_some(Ok((s, t))),
        });
        if let Some(found) = found {
            let (s, t) = found?;
            return Err(refuse(
                "import.mesh-crossing",
                format!(
                    "triangles {s} and {t} meet away from their shared edges and vertices, near {}",
                    self.describe(self.triangles[s as usize][0])
                ),
            ));
        }
        Ok(())
    }

    fn find_shells(&mut self) {
        let n = self.triangles.len();
        let mut shell_of = vec![u32::MAX; n];
        let mut shells = Vec::new();
        for start in 0..n {
            if shell_of[start] != u32::MAX {
                continue;
            }
            let id = shells.len() as u32;
            let mut members = Vec::new();
            let mut stack = vec![start as u32];
            shell_of[start] = id;
            while let Some(t) = stack.pop() {
                members.push(t);
                for &u in &self.twins[t as usize] {
                    if shell_of[u as usize] == u32::MAX {
                        shell_of[u as usize] = id;
                        stack.push(u);
                    }
                }
            }
            members.sort_unstable();
            shells.push(MeshShell {
                triangles: members,
                outward: true,
                parent: None,
            });
        }
        self.shells = shells;
    }

    /// Sign of six times the signed volume of one shell. A closed surface that
    /// does not cross itself encloses positive volume, so the sign is never zero.
    fn shell_orientation(&self, shell: &MeshShell, meter: &mut Meter) -> Result<Sign, Refusal> {
        meter.tick(shell.triangles.len() as u64)?;
        let mut sum = 0.0f64;
        let mut error = 0.0f64;
        let mut magnitude = 0.0f64;
        for &t in &shell.triangles {
            let [a, b, c] = self.triangles[t as usize].map(|v| self.approx[v as usize]);
            let det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]);
            let permanent = a[0].abs() * (b[1].abs() * c[2].abs() + b[2].abs() * c[1].abs())
                + a[1].abs() * (b[0].abs() * c[2].abs() + b[2].abs() * c[0].abs())
                + a[2].abs() * (b[0].abs() * c[1].abs() + b[1].abs() * c[0].abs());
            sum += det;
            error += permanent;
            magnitude += det.abs();
        }
        let n = shell.triangles.len() as f64;
        let bound = 1e-14 * error + 4.0 * n * f64::EPSILON * magnitude;
        // As in the predicates, the filter holds only while every copy is
        // close to its exact coordinate and no product can underflow or overflow.
        let in_range = shell.triangles.iter().all(|&t| {
            let corners = self.triangles[t as usize];
            corners.iter().all(|&v| self.close[v as usize])
                && predicates::in_range(corners.iter().flat_map(|&v| &self.approx[v as usize]))
        });
        if in_range && sum.is_finite() && bound.is_finite() {
            if sum > bound {
                return Ok(Sign::Positive);
            }
            if sum < -bound {
                return Ok(Sign::Negative);
            }
        }
        let mut exact = Rational::zero();
        for &t in &shell.triangles {
            let [a, b, c] = self.triangles[t as usize].map(|v| &self.exact[v as usize]);
            predicates::record_exact(15, a.iter().chain(b).chain(c));
            meter.tick(0)?;
            let det = &(&a[0] * &(&(&b[1] * &c[2]) - &(&b[2] * &c[1])))
                - &(&a[1] * &(&(&b[0] * &c[2]) - &(&b[2] * &c[0])));
            let det = &det + &(&a[2] * &(&(&b[0] * &c[1]) - &(&b[1] * &c[0])));
            exact = &exact + &det;
        }
        Ok(exact.sign())
    }

    /// Walks the `+x` ray from `q`, moved by `(0, e, e^2)`, through every
    /// triangle not skipped, reporting each crossing in front of `q` with the
    /// sign of the triangle's `x` normal component.
    fn ray(
        &self,
        q: P,
        skip: impl Fn(u32) -> bool,
        mut hit: impl FnMut(u32, i64),
        meter: &mut Meter,
    ) -> Result<Ray, Refusal> {
        let lo = q.f.map(f64::next_down);
        let hi = q.f.map(f64::next_up);
        let query = ([lo[0], lo[1], lo[2]], [f64::INFINITY, hi[1], hi[2]]);
        let mut boundary = false;
        let mut failure = None;
        self.tree.query(&query, |t| {
            if let Err(refusal) = meter.tick(1) {
                failure = Some(refusal);
                return false;
            }
            if skip(t) {
                return true;
            }
            let [a, b, c] = self.corners(t);
            let s = orient2d(0, a, b, c);
            if s == Sign::Zero
                || [(a, b), (b, c), (c, a)]
                    .iter()
                    .any(|&(x, y)| perturbed_yz(x, y, q) != s)
            {
                return true;
            }
            match orient3d(a, b, c, q) {
                Sign::Zero => boundary = true,
                side if side != s => hit(t, if s == Sign::Positive { 1 } else { -1 }),
                _ => {}
            }
            !boundary
        });
        match failure {
            Some(refusal) => Err(refusal),
            None if boundary => Ok(Ray::Boundary),
            None => Ok(Ray::Clear),
        }
    }

    /// Whether `q` lies on a triangle.
    fn on_surface(&self, q: P, meter: &mut Meter) -> Result<bool, Refusal> {
        let query = (q.f.map(f64::next_down), q.f.map(f64::next_up));
        let mut found = false;
        let mut failure = None;
        self.tree.query(&query, |t| {
            if let Err(refusal) = meter.tick(1) {
                failure = Some(refusal);
                return false;
            }
            let c = self.corners(t);
            if orient3d(c[0], c[1], c[2], q) == Sign::Zero {
                let (k, s) = projection(c[0], c[1], c[2]).expect("nondegenerate");
                found = in_triangle_2d(k, s, c, q);
            }
            !found
        });
        failure.map_or(Ok(found), Err)
    }

    fn check_windings(&mut self, meter: &mut Meter) -> Result<(), Refusal> {
        let orientation = self
            .shells
            .iter()
            .map(|s| {
                Ok(if self.shell_orientation(s, meter)? == Sign::Positive {
                    1
                } else {
                    -1
                })
            })
            .collect::<Result<Vec<i64>, Refusal>>()?;
        let count = self.shells.len();
        let mut shell_of = vec![0u32; self.triangles.len()];
        for (i, s) in self.shells.iter().enumerate() {
            for &t in &s.triangles {
                shell_of[t as usize] = i as u32;
            }
        }
        // containers[s]: shells whose winding contribution at s is nonzero.
        let mut containers: Vec<Vec<usize>> = vec![Vec::new(); count];
        if count > 1 {
            for (i, shell) in self.shells.iter().enumerate() {
                let q = self.p(self.triangles[shell.triangles[0] as usize][0]);
                let mut windings = vec![0i64; count];
                let walk = self.ray(
                    q,
                    |t| shell_of[t as usize] == i as u32,
                    |t, w| {
                        windings[shell_of[t as usize] as usize] += w;
                    },
                    meter,
                )?;
                debug_assert!(matches!(walk, Ray::Clear), "shells do not touch");
                containers[i] = (0..count).filter(|&j| windings[j] != 0).collect();
            }
        }
        for i in 0..count {
            let outside: i64 = containers[i].iter().map(|&j| orientation[j]).sum();
            let inside = outside + orientation[i];
            if !matches!((outside, inside), (0, 1) | (1, 0)) {
                return Err(refuse(
                    "import.mesh-winding",
                    format!(
                        "shell {i} has winding number {outside} outside and {inside} inside; \
                         a solid needs zero outside and one inside, so the triangles face the wrong way"
                    ),
                ));
            }
            self.shells[i].outward = orientation[i] == 1;
            self.shells[i].parent = containers[i]
                .iter()
                .copied()
                .max_by_key(|&j| containers[j].len());
        }
        Ok(())
    }

    fn file_point(&self, p: &[Rational; 3]) -> [Rational; 3] {
        if self.scale == Rational::one() {
            p.clone()
        } else {
            p.each_ref()
                .map(|x| x.checked_div(&self.scale).expect("positive scale"))
        }
    }

    /// The integer winding number at a point in scaled mesh coordinates, or
    /// `None` on the surface. It is one inside the solid and zero outside.
    /// Work is charged to `dec`.
    pub fn winding(&self, p: &[Rational; 3], dec: &mut Decider) -> Result<Option<i64>, Refusal> {
        if self.triangles.is_empty() {
            return Ok(Some(0));
        }
        let mut meter = Meter::new(dec);
        let file = self.file_point(p);
        let q = P::new(&file);
        let winding = if self.on_surface(q, &mut meter)? {
            None
        } else {
            let mut winding = 0;
            match self.ray(q, |_| false, |_, w| winding += w, &mut meter)? {
                Ray::Boundary => None,
                Ray::Clear => Some(winding),
            }
        };
        meter.flush()?;
        Ok(winding)
    }

    /// A tentative winding number for display, counted in binary64 with the
    /// same perturbed ray. It carries no proof and may be wrong near the surface.
    pub fn fast_winding(&self, p: [f64; 3]) -> i64 {
        let s = self.scale.to_f64_nearest();
        let q = p.map(|x| x / s);
        let query = ([q[0], q[1], q[2]], [f64::INFINITY, q[1], q[2]]);
        let sign = |x: f64| i32::from(x > 0.0) - i32::from(x < 0.0);
        let orient = |a: [f64; 3], b: [f64; 3]| {
            let base = sign((b[1] - a[1]) * (q[2] - a[2]) - (b[2] - a[2]) * (q[1] - a[1]));
            match (base, sign(a[2] - b[2])) {
                (0, 0) => sign(b[1] - a[1]),
                (0, tie) => tie,
                (base, _) => base,
            }
        };
        let mut winding = 0;
        self.tree.query(&query, |t| {
            let [a, b, c] = self.triangles[t as usize].map(|v| self.approx[v as usize]);
            let area = sign((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]));
            if area == 0
                || [(a, b), (b, c), (c, a)]
                    .iter()
                    .any(|&(x, y)| orient(x, y) != area)
            {
                return true;
            }
            let n: [f64; 3] = {
                let (u, v) = (
                    std::array::from_fn::<f64, 3, _>(|i| b[i] - a[i]),
                    std::array::from_fn::<f64, 3, _>(|i| c[i] - a[i]),
                );
                [
                    u[1] * v[2] - u[2] * v[1],
                    u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0],
                ]
            };
            let side = (q[0] - a[0]) * n[0] + (q[1] - a[1]) * n[1] + (q[2] - a[2]) * n[2];
            if sign(side) == -area {
                winding += i64::from(area);
            }
            true
        });
        winding
    }

    /// Exact squared distance from a point in scaled mesh coordinates to the
    /// surface. `None` for the empty mesh. Work is charged to `dec`.
    pub fn squared_distance(
        &self,
        p: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Option<Rational>, Refusal> {
        let mut meter = Meter::new(dec);
        let mut failure = None;
        let file = self.file_point(p);
        let q = P::new(&file);
        let lo = q.f.map(f64::next_down);
        let hi = q.f.map(f64::next_up);
        let point = file.each_ref().map(Interval::from_rational);
        let mut best: Option<Rational> = None;
        let mut best_upper = f64::INFINITY;
        self.tree.nearest(
            |b| box_lower_bound(b, &lo, &hi),
            |t| {
                if let Err(refusal) = meter.tick(1) {
                    failure = Some(refusal);
                    return None;
                }
                // Distance to the triangle is at least distance to its plane;
                // an enclosure of that rules most candidates out without
                // exact work when many triangles are nearly as close.
                if self.plane_lower_bound(t, &point) > best_upper {
                    return Some(f64::INFINITY);
                }
                let [a, b, c] = self.triangles[t as usize].map(|v| &self.exact[v as usize]);
                predicates::record_exact(60, a.iter().chain(b).chain(c).chain(&file));
                let d = triangle_distance_sq(&file, a, b, c);
                let upper = d.to_f64_nearest().next_up();
                best_upper = best_upper.min(upper);
                if best.as_ref().is_none_or(|x| &d < x) {
                    best = Some(d);
                }
                Some(upper)
            },
        );
        if let Some(refusal) = failure {
            return Err(refusal);
        }
        meter.flush()?;
        Ok(best.map(|d| &d * &(&self.scale * &self.scale)))
    }

    /// A rigorous lower bound on the squared distance from `p` to the plane of
    /// triangle `t`, in file units.
    fn plane_lower_bound(&self, t: u32, p: &[Interval; 3]) -> f64 {
        let enclose = |v: u32| {
            self.approx[v as usize]
                .map(|x| Interval::new(x.next_down(), x.next_up()).unwrap_or(Interval::entire()))
        };
        let [a, b, c] = self.triangles[t as usize].map(enclose);
        let u: [Interval; 3] = std::array::from_fn(|i| b[i] - a[i]);
        let v: [Interval; 3] = std::array::from_fn(|i| c[i] - a[i]);
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let side = (p[0] - a[0]) * n[0] + (p[1] - a[1]) * n[1] + (p[2] - a[2]) * n[2];
        let length = n[0].square() + n[1].square() + n[2].square();
        let bound = (side.square() / length).lo();
        if bound.is_nan() { 0.0 } else { bound }
    }

    /// Approximate squared distance in scaled units, for display only.
    pub fn fast_squared_distance(&self, p: [f64; 3]) -> f64 {
        let s = self.scale.to_f64_nearest();
        let q = p.map(|x| x / s);
        let mut best = f64::INFINITY;
        self.tree.nearest(
            |b| fast_box_distance(b, &q),
            |t| {
                let [a, b, c] = self.triangles[t as usize].map(|v| self.approx[v as usize]);
                let d = fast_triangle_distance_sq(q, a, b, c);
                best = best.min(d);
                Some(d)
            },
        );
        best * s * s
    }

    /// The scale applied to file coordinates, in millimetres per file unit.
    pub fn scale(&self) -> &Rational {
        &self.scale
    }

    /// Number of welded vertices.
    pub fn vertex_count(&self) -> usize {
        self.exact.len()
    }

    /// The exact vertex in file units.
    pub fn file_vertex(&self, v: u32) -> &[Rational; 3] {
        &self.exact[v as usize]
    }

    /// The exact vertex after scaling, in millimetres.
    pub fn vertex(&self, v: u32) -> [Rational; 3] {
        self.exact[v as usize].each_ref().map(|x| x * &self.scale)
    }

    /// Triangles as vertex indices, in file order, each facing outward from
    /// the material.
    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.triangles
    }

    /// The triangle across the edge from corner `slot` to the next corner.
    pub fn twin(&self, t: u32, slot: usize) -> u32 {
        self.twins[t as usize][slot]
    }

    /// Whether the triangle across edge `slot` lies in the same plane. Work is
    /// charged to `dec`.
    pub fn coplanar_across(&self, t: u32, slot: usize, dec: &mut Decider) -> Result<bool, Refusal> {
        let u = self.twin(t, slot);
        let far = *self.triangles[u as usize]
            .iter()
            .find(|v| !self.triangles[t as usize].contains(v))
            .expect("twins share exactly one edge");
        let [a, b, c] = self.corners(t);
        let coplanar = orient3d(a, b, c, self.p(far)) == Sign::Zero;
        dec.charge_work(1 + predicates::take_exact_work())?;
        Ok(coplanar)
    }

    /// A projection axis in which triangle `t` keeps its area, with the sign
    /// of the projected area.
    pub fn projection(&self, t: u32) -> (usize, Sign) {
        let [a, b, c] = self.corners(t);
        projection(a, b, c).expect("validated triangles are nondegenerate")
    }

    /// Sign of the exact projected area of a closed vertex loop, axis `k`
    /// dropped. Work is charged to `dec`.
    pub fn loop_area_sign(
        &self,
        k: usize,
        vertices: &[u32],
        dec: &mut Decider,
    ) -> Result<Sign, Refusal> {
        let (i, j) = predicates::kept(k);
        let mut meter = Meter::new(dec);
        let mut twice = Rational::zero();
        for (n, &v) in vertices.iter().enumerate() {
            let a = &self.exact[v as usize];
            let b = &self.exact[vertices[(n + 1) % vertices.len()] as usize];
            predicates::record_exact(4, [&a[i], &a[j], &b[i], &b[j], &twice]);
            meter.tick(1)?;
            twice = &twice + &(&(&a[i] * &b[j]) - &(&a[j] * &b[i]));
        }
        meter.flush()?;
        Ok(twice.sign())
    }

    /// The connected closed surfaces, with their orientation and nesting.
    pub fn shells(&self) -> &[MeshShell] {
        &self.shells
    }

    /// Solid pieces: each outer shell with the cavities directly inside it.
    pub fn lumps(&self) -> Vec<(usize, Vec<usize>)> {
        let mut lumps: Vec<(usize, Vec<usize>)> = self
            .shells
            .iter()
            .enumerate()
            .filter(|(_, s)| s.outward)
            .map(|(i, _)| (i, Vec::new()))
            .collect();
        for (i, s) in self.shells.iter().enumerate() {
            if !s.outward {
                let parent = s.parent.expect("a cavity lies inside an outer shell");
                let lump = lumps
                    .iter_mut()
                    .find(|(outer, _)| *outer == parent)
                    .expect("a cavity's innermost container is an outer shell");
                lump.1.push(i);
            }
        }
        lumps
    }

    /// Exact bounds of the scaled vertices, `None` for the empty mesh.
    pub fn bounds(&self) -> Option<[[Rational; 2]; 3]> {
        let first = self.triangles.first()?[0];
        let mut bounds: [[Rational; 2]; 3] = std::array::from_fn(|i| {
            let x = &self.exact[first as usize][i];
            [x.clone(), x.clone()]
        });
        for tri in &self.triangles {
            for &v in tri {
                for (i, b) in bounds.iter_mut().enumerate() {
                    let x = &self.exact[v as usize][i];
                    if x < &b[0] {
                        b[0] = x.clone();
                    }
                    if x > &b[1] {
                        b[1] = x.clone();
                    }
                }
            }
        }
        Some(bounds.map(|b| b.map(|x| &x * &self.scale)))
    }
}
