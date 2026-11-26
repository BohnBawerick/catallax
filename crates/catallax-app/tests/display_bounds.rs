//! Sample the recipe corpus without using mesher charts, sample locations or proof predicates.
use std::{collections::BTreeMap, fs, path::Path, time::Instant};

use catallax_app::{EvaluatedApplicationOutcome, build_incremental_with_files_without_display};
use catallax_features::ImportFiles;
use catallax_intersect::{TracePiece, trace};
use catallax_mesh::{Mesh, MeshOptions, mesh_body};
use catallax_numbers::{Ball, Decider, Rational, Real, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{BSplineSurface, Point3, SaddleRing, SplineSupport, Surface};
use catallax_topology::{FaceId, Loop, Model, Sense};
use serde_json::Value;

const GRID: usize = 19;
const EDGE_STEPS: usize = 31;
const TAU: f64 = std::f64::consts::TAU;
type Point = [f64; 3];
type Uv = [f64; 2];
type Triangle = [Point; 3];

fn dot(a: Point, b: Point) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn sub(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}
fn mix(a: Point, b: Point, t: f64) -> Point {
    std::array::from_fn(|i| a[i] + (b[i] - a[i]) * t)
}
fn cross(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn unit(a: Point) -> Point {
    let length = dot(a, a).sqrt();
    a.map(|x| x / length)
}
fn gap(a: Point, b: Point) -> f64 {
    let d = sub(a, b);
    dot(d, d)
}
fn segment_gap(p: Point, a: Point, b: Point) -> f64 {
    let edge = sub(b, a);
    let length = dot(edge, edge);
    gap(
        p,
        if length == 0.0 {
            a
        } else {
            mix(a, b, (dot(sub(p, a), edge) / length).clamp(0.0, 1.0))
        },
    )
}
fn triangle_gap(p: Point, [a, b, c]: Triangle) -> f64 {
    let u = sub(b, a);
    let v = sub(c, a);
    let d = sub(p, a);
    let n = cross(u, v);
    let det = dot(n, n);
    if det > 0.0 {
        let s = dot(cross(d, v), n) / det;
        let t = dot(cross(u, d), n) / det;
        if s >= 0.0 && t >= 0.0 && s + t <= 1.0 {
            return dot(d, n).powi(2) / det;
        }
    }
    segment_gap(p, a, b)
        .min(segment_gap(p, b, c))
        .min(segment_gap(p, c, a))
}
fn near_mesh(p: Point, triangles: &[Triangle], bound: f64) -> bool {
    triangles.iter().any(|t| {
        (0..3).all(|i| {
            !t.iter().all(|v| v[i] < p[i] - bound) && !t.iter().all(|v| v[i] > p[i] + bound)
        }) && triangle_gap(p, *t) <= bound * bound
    })
}
fn number(x: &Real) -> f64 {
    x.rational_ref()
        .map_or_else(|| x.to_f64(), Rational::to_f64_nearest)
}
fn xyz(p: &Point3) -> Point {
    [&p.x, &p.y, &p.z].map(number)
}

// Basis functions at one parameter; only the degree + 1 nonzero entries are
// evaluated. This uses raw spline coefficients, not an interval mesh chart.
struct Basis {
    degree: usize,
    count: usize,
    knots: Vec<f64>,
}
impl Basis {
    fn domain(&self) -> [f64; 2] {
        [self.knots[self.degree], self.knots[self.count]]
    }
    fn at(&self, t: f64, weights: &mut [f64]) -> usize {
        let span = (self.knots.partition_point(|&k| k <= t) - 1).min(self.count - 1);
        weights[0] = 1.0;
        for j in 1..=self.degree {
            let mut saved = 0.0;
            for i in 0..j {
                let a = self.knots[span + 1 - j + i];
                let b = self.knots[span + i + 1];
                let w = weights[i] / (b - a);
                weights[i] = saved + (b - t) * w;
                saved = (t - a) * w;
            }
            weights[j] = saved;
        }
        span - self.degree
    }
}
struct Spline {
    axis: [Basis; 2],
    controls: Vec<(Point, f64)>,
}
impl Spline {
    fn new(s: &BSplineSurface) -> Self {
        let axis = std::array::from_fn(|i| Basis {
            degree: s.degrees()[i],
            count: s.control_count()[i],
            knots: s.knots()[i].iter().map(Rational::to_f64_nearest).collect(),
        });
        assert!(
            axis.iter().all(|a| a.degree < 16),
            "sampler supports degrees below 16"
        );
        let controls = s
            .controls()
            .iter()
            .zip(s.weights())
            .map(|(p, w)| (xyz(p), number(w)))
            .collect();
        Self { axis, controls }
    }
    fn at(&self, [u, v]: Uv) -> Point {
        let [a, b] = &self.axis;
        let (mut bu, mut bv) = ([0.0; 16], [0.0; 16]);
        let (i0, j0) = (a.at(u, &mut bu), b.at(v, &mut bv));
        let (mut p, mut sum) = ([0.0; 3], 0.0);
        for i in 0..=a.degree {
            for j in 0..=b.degree {
                let (c, w) = self.controls[(i0 + i) * b.count + j0 + j];
                let w = w * bu[i] * bv[j];
                sum += w;
                for k in 0..3 {
                    p[k] += w * c[k];
                }
            }
        }
        p.map(|x| x / sum)
    }
    fn locate(&self, p: Point, seeds: &[(Uv, Point)], periodic: [bool; 2]) -> Uv {
        let mut uv = seeds
            .iter()
            .min_by(|a, b| gap(a.1, p).total_cmp(&gap(b.1, p)))
            .unwrap()
            .0;
        let domain = self.axis.each_ref().map(Basis::domain);
        for _ in 0..12 {
            let q = self.at(uv);
            let mut jac = [[0.0; 3]; 2];
            for i in 0..2 {
                let [lo, hi] = domain[i];
                let h = (hi - lo) * 1e-5;
                let (mut a, mut b) = (uv, uv);
                a[i] = (uv[i] - h).max(lo);
                b[i] = (uv[i] + h).min(hi);
                jac[i] = sub(self.at(b), self.at(a)).map(|x| x / (b[i] - a[i]));
            }
            let [a, b] = jac;
            let d = sub(p, q);
            let (aa, ab, bb, da, db) = (dot(a, a), dot(a, b), dot(b, b), dot(d, a), dot(d, b));
            let det = aa * bb - ab * ab;
            if det <= 0.0 {
                break;
            }
            let step = [(bb * da - ab * db) / det, (aa * db - ab * da) / det];
            let mut moved = false;
            for scale in [1.0, 0.5, 0.25, 0.125] {
                let next = std::array::from_fn(|i| {
                    let [lo, hi] = domain[i];
                    let t = uv[i] + scale * step[i];
                    if periodic[i] {
                        lo + (t - lo).rem_euclid(hi - lo)
                    } else {
                        t.clamp(lo, hi)
                    }
                });
                if gap(p, self.at(next)) < gap(p, q) {
                    uv = next;
                    moved = true;
                    break;
                }
            }
            if !moved {
                break;
            }
        }
        uv
    }
}

#[derive(Clone, Copy)]
enum Shape {
    Plane,
    Cylinder(f64, f64),
    Cone(f64),
    Sphere(f64),
    Torus(f64, f64),
}
struct Analytic {
    shape: Shape,
    origin: Point,
    frame: [Point; 3],
}
impl Analytic {
    fn new(surface: &Surface) -> Self {
        let (shape, origin, axis, reference) = match surface {
            Surface::Plane(s) => (Shape::Plane, &s.origin, &s.normal, &s.ref_dir),
            Surface::Cylinder(s) => {
                let r = number(&s.radius());
                (Shape::Cylinder(r, r), &s.origin, &s.axis, &s.ref_dir)
            }
            Surface::EllipticCylinder(s) => (
                Shape::Cylinder(number(&s.major_sq).sqrt(), number(&s.minor_sq).sqrt()),
                &s.origin,
                &s.axis,
                &s.major_dir,
            ),
            Surface::Cone(s) => (Shape::Cone(number(&s.slope)), &s.apex, &s.axis, &s.ref_dir),
            Surface::Sphere(s) => (
                Shape::Sphere(number(&s.radius)),
                &s.center,
                &s.axis,
                &s.ref_dir,
            ),
            Surface::Torus(s) => (
                Shape::Torus(number(&s.major), number(&s.minor)),
                &s.center,
                &s.axis,
                &s.ref_dir,
            ),
            Surface::TorusSheet(s) => (
                Shape::Torus(number(&s.major), number(&s.minor)),
                &s.center,
                &s.axis,
                &s.ref_dir,
            ),
            _ => panic!("no independent chart for {:?}", surface.kind()),
        };
        let x = unit([&reference.x, &reference.y, &reference.z].map(number));
        let z = unit([&axis.x, &axis.y, &axis.z].map(number));
        Self {
            shape,
            origin: xyz(origin),
            frame: [x, cross(z, x), z],
        }
    }
    fn at(&self, [u, v]: Uv) -> Point {
        let (s, c) = u.sin_cos();
        let p = match self.shape {
            Shape::Plane => [u, v, 0.0],
            Shape::Cylinder(a, b) => [a * c, b * s, v],
            Shape::Cone(k) => [k * v * c, k * v * s, v],
            Shape::Sphere(r) => [r * v.cos() * c, r * v.cos() * s, r * v.sin()],
            Shape::Torus(r, a) => [(r + a * v.cos()) * c, (r + a * v.cos()) * s, a * v.sin()],
        };
        std::array::from_fn(|j| {
            self.origin[j] + (0..3).map(|i| p[i] * self.frame[i][j]).sum::<f64>()
        })
    }
    fn uv(&self, p: Point, project: bool) -> Uv {
        let d = sub(p, self.origin);
        let [x, y, z] = self.frame.map(|a| dot(a, d));
        let u = y.atan2(x);
        let r = x.hypot(y);
        match self.shape {
            Shape::Plane => [x, y],
            Shape::Cylinder(a, b) => [(y / b).atan2(x / a), z],
            Shape::Cone(k) => [
                u,
                if project {
                    ((z + k * r) / (1.0 + k * k)).max(0.0)
                } else {
                    z
                },
            ],
            Shape::Sphere(_) => [u, z.atan2(r)],
            Shape::Torus(a, _) => [u, z.atan2(r - a)],
        }
    }
    fn periods(&self) -> Uv {
        match self.shape {
            Shape::Plane => [0.0, 0.0],
            Shape::Torus(..) => [TAU, TAU],
            _ => [TAU, 0.0],
        }
    }
}

fn mesh_triangles(mesh: &Mesh, range: std::ops::Range<usize>) -> Vec<Triangle> {
    mesh.triangles[range]
        .iter()
        .map(|t| t.vertices.map(|i| mesh.positions[i]))
        .collect()
}

struct TrimArc {
    piece: TracePiece,
    points: [Point; EDGE_STEPS + 1],
}

impl TrimArc {
    fn point(&self, t: f64, dec: &Decider) -> Result<Point, catallax_numbers::Refusal> {
        let t = Ball::from_rational(&Rational::from_f64_exact(t).unwrap(), 80);
        Ok(self
            .piece
            .evaluate(&t, dec)?
            .map(|x| x.mid().to_f64_nearest()))
    }
}

struct Trims {
    loops: Vec<Vec<Point>>,
    arcs: Vec<TrimArc>,
}

impl Trims {
    fn nearest_boundary_sample(&self, p: Point, bound: f64, dec: &Decider) -> f64 {
        let mut best = f64::INFINITY;
        let mut nearest = None::<(f64, &TrimArc, usize, f64)>;
        for arc in &self.arcs {
            for (i, pair) in arc.points.windows(2).enumerate() {
                let [a, b] = [pair[0], pair[1]];
                best = best.min(gap(p, a)).min(gap(p, b));
                let ab = sub(b, a);
                let length = dot(ab, ab);
                if length == 0.0 {
                    continue;
                }
                let t = (dot(sub(p, a), ab) / length).clamp(0.0, 1.0);
                let d = gap(p, mix(a, b, t));
                if nearest.is_none_or(|(old, ..)| d < old) {
                    nearest = Some((d, arc, i, t));
                }
            }
        }
        if best <= bound * bound {
            return best;
        }
        if let Some((_, arc, i, t)) = nearest {
            // Chords of the independent trim samples seed a curve search but
            // are not witnesses themselves. Evaluate the source arc again.
            let eval = |t| arc.point(t, dec).map_or(f64::INFINITY, |q| gap(p, q));
            let n = EDGE_STEPS as f64;
            best = best.min(eval((i as f64 + t) / n));
            let mut lo = (i as f64 - 1.0).max(0.0) / n;
            let mut hi = ((i + 2) as f64 / n).min(1.0);
            for _ in 0..20 {
                if best <= bound * bound {
                    break;
                }
                let a = (2.0 * lo + hi) / 3.0;
                let b = (lo + 2.0 * hi) / 3.0;
                let (da, db) = (eval(a), eval(b));
                best = best.min(da).min(db);
                if da < db {
                    hi = b;
                } else {
                    lo = a;
                }
            }
        }
        best
    }
}

fn trim_loops(
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Trims, catallax_numbers::Refusal> {
    let mut arcs = Vec::new();
    let loops = model
        .face(face)
        .loops
        .iter()
        .map(|&ring| match model.loop_(ring) {
            Loop::Vertex(v) => Ok(vec![xyz(&model.geom().position(model.vertex(*v).point))]),
            Loop::Coedges(edges) => {
                let mut points = Vec::new();
                for &id in edges {
                    let coedge = model.coedge(id);
                    let edge = model.edge(coedge.edge);
                    let endpoints =
                        [edge.start, edge.end].map(|v| v.map(|v| model.vertex(v).point));
                    let forward = coedge.sense == edge.sense;
                    // Trace in the carrier's natural direction, then reverse
                    // traversal. Reversing the endpoints could select the other arc.
                    let [a, b] = if edge.sense == Sense::Same {
                        endpoints
                    } else {
                        [endpoints[1], endpoints[0]]
                    };
                    let path = trace(model.geom(), edge.curve, a, b, dec)?;
                    let pieces = path.value.pieces();
                    for j in 0..pieces.len() {
                        let piece = pieces[if forward { j } else { pieces.len() - 1 - j }].clone();
                        let mut arc = TrimArc {
                            piece,
                            points: [[0.0; 3]; EDGE_STEPS + 1],
                        };
                        for i in 0..=EDGE_STEPS {
                            arc.points[i] = arc.point(i as f64 / EDGE_STEPS as f64, dec)?;
                        }
                        if forward {
                            points.extend(arc.points[..EDGE_STEPS].iter().copied());
                        } else {
                            points.extend(arc.points[1..].iter().rev().copied());
                        }
                        arcs.push(arc);
                    }
                }
                Ok(points)
            }
        })
        .collect::<Result<_, catallax_numbers::Refusal>>()?;
    Ok(Trims { loops, arcs })
}

fn delta(x: f64, origin: f64, period: f64) -> f64 {
    let d = x - origin;
    if period > 0.0 {
        d - (d / period).round() * period
    } else {
        d
    }
}

// The closest directed crossing on a chart line determines which side of the
// trim is material. Work in periodic coordinates so noncontractible rings do
// not gain fictitious edges at the seam. Ignore probes near a sampled trim.
fn inside(uv: Uv, loops: &[Vec<Uv>], periods: Uv, interior_left: bool, guard: Uv) -> bool {
    let mut closest = None::<(f64, bool, f64)>;
    for i in 0..2 {
        let j = 1 - i;
        for ring in loops.iter().filter(|r| r.len() > 1) {
            for (&a, &b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
                let mut x = delta(a[i], uv[i], periods[i]);
                let dx = delta(b[i], a[i], periods[i]);
                if periods[i] > 0.0 {
                    x = delta(x + 0.5 * dx, 0.0, periods[i]) - 0.5 * dx;
                }
                let y = x + dx;
                if (x <= 0.0) == (y <= 0.0) {
                    continue;
                }
                let at = a[j] + delta(b[j], a[j], periods[j]) * (-x / dx);
                let offset = delta(uv[j], at, periods[j]);
                let distance = offset.abs();
                let signed = (offset > 0.0) == (dx > 0.0);
                if closest.is_none_or(|(d, _, _)| distance < d) {
                    closest = Some((distance, signed == (interior_left == (i == 0)), guard[j]));
                }
            }
        }
        // The scan axes have different units on some charts. Do not compare
        // distances from different axes; use the second axis only as fallback.
        if closest.is_some() {
            break;
        }
    }
    closest.map_or_else(
        || loops.iter().all(|r| r.len() <= 1),
        |(d, left, guard)| left && d > guard,
    )
}

fn sampled_range(domain: [[f64; 2]; 2]) -> impl Iterator<Item = Uv> {
    (0..GRID).flat_map(move |u| {
        (0..GRID).map(move |v| {
            std::array::from_fn(|i| {
                let t = ([u, v][i] as f64 + 0.5) / GRID as f64;
                let [a, b] = domain[i];
                a + (b - a) * t
            })
        })
    })
}

fn check_face(
    case: &str,
    model: &Model,
    mesh: &Mesh,
    face: &catallax_mesh::FaceMesh,
    dec: &mut Decider,
    findings: &mut Vec<String>,
) {
    let input = model.face(face.face);
    let surface = model.geom().surface(input.surface);
    // Recognize a whole band by its retained rail definitions and actual
    // coedges, not the fact that a producer supplied a spline carrier.
    // A crossing's two ends collapse to poles; UV inversion at those points
    // cannot recover the side of the band, and they are not periodic seams.
    let whole_band = match surface {
        Surface::BSpline(SplineSupport::RollingSaddle(_) | SplineSupport::RollingSphere(_)) => {
            catallax_topology::saddle_rows(model, face.face).is_some()
        }
        Surface::BSpline(SplineSupport::RollingPipe(_)) => {
            catallax_topology::pipe_rows(model, face.face, dec).is_ok_and(|rails| rails.is_some())
        }
        _ => false,
    };
    let atlas_domains: Vec<_> = match surface {
        Surface::BSpline(SplineSupport::RationalAtlas(atlas)) => atlas
            .charts()
            .iter()
            .map(|chart| {
                chart
                    .patch()
                    .domain()
                    .each_ref()
                    .map(|axis| axis.each_ref().map(Rational::to_f64_nearest))
            })
            .collect(),
        _ => Vec::new(),
    };
    let in_atlas = |uv: Uv| {
        atlas_domains.is_empty()
            || atlas_domains
                .iter()
                .any(|cell| (0..2).all(|i| (cell[i][0]..=cell[i][1]).contains(&uv[i])))
    };
    let triangles = mesh_triangles(mesh, face.triangles.clone());
    let label = format!(
        "{case} face {} {:?} {:?}",
        face.face.index(),
        surface.kind(),
        model.face_origins(face.face)
    );
    assert!(!triangles.is_empty(), "{label}: no triangles");
    let scale = triangles[0]
        .iter()
        .flatten()
        .map(|x| x.abs())
        .fold(1.0, f64::max);
    let face_bound = face.max_chord_error.to_f64_nearest() + scale * 2e-12;
    // Face reports cover the carrier-to-facet approximation; incident edge
    // reports separately cover trim displacement (notably on affine faces).
    let mut bound = face_bound;
    for &ring in &input.loops {
        if let Loop::Coedges(edges) = model.loop_(ring) {
            for &id in edges {
                let edge = model.coedge(id).edge;
                if let Some(line) = mesh.edges.iter().find(|e| e.edge == edge) {
                    bound = bound.max(line.max_chord_error.to_f64_nearest() + scale * 2e-12);
                }
            }
        }
    }
    let trims = match trim_loops(model, face.face, dec) {
        Ok(trims) => trims,
        Err(error) => {
            findings.push(format!("{label}: independent trim trace refused: {error}"));
            return;
        }
    };
    let loops = &trims.loops;
    let mut domain;
    let analytic;
    let spline;
    let seeds;
    let sample_spline;
    let invert_spline;
    let sample_analytic;
    let invert_analytic;
    let (at, uv, periods): (&dyn Fn(Uv) -> Point, &dyn Fn(Point) -> Uv, Uv) = match surface {
        Surface::BSpline(s) => {
            spline = Spline::new(s.carrier());
            domain = spline.axis.each_ref().map(Basis::domain);
            seeds = if atlas_domains.is_empty() {
                sampled_range(domain)
                    .map(|q| (q, spline.at(q)))
                    .collect::<Vec<_>>()
            } else {
                atlas_domains
                    .iter()
                    .flat_map(|&cell| sampled_range(cell))
                    .map(|q| (q, spline.at(q)))
                    .collect()
            };
            let mut periodic = s
                .param_kinds()
                .map(|k| matches!(k, catallax_shapes::ParamKind::PeriodicFinite));
            if matches!(s, SplineSupport::RollingSaddle(saddle)
                if matches!(saddle.definition().law().ring, SaddleRing::Crossing { .. }))
            {
                periodic[0] = false;
            }
            sample_spline = |q| spline.at(q);
            let spline_ref = &spline;
            let seed_ref = &seeds;
            invert_spline = move |p| spline_ref.locate(p, seed_ref, periodic);
            (
                &sample_spline,
                &invert_spline,
                std::array::from_fn(|i| {
                    if periodic[i] {
                        domain[i][1] - domain[i][0]
                    } else {
                        0.0
                    }
                }),
            )
        }
        _ => {
            analytic = Analytic::new(surface);
            domain = analytic.periods().map(|period| {
                if period > 0.0 {
                    [-period / 2.0, period / 2.0]
                } else {
                    [f64::INFINITY, f64::NEG_INFINITY]
                }
            });
            for &p in loops.iter().flatten() {
                let q = analytic.uv(p, false);
                for i in 0..2 {
                    if analytic.periods()[i] == 0.0 {
                        domain[i][0] = domain[i][0].min(q[i]);
                        domain[i][1] = domain[i][1].max(q[i]);
                    }
                }
            }
            if matches!(analytic.shape, Shape::Sphere(_)) {
                domain[1] = [-TAU / 4.0, TAU / 4.0];
            }
            sample_analytic = |q| analytic.at(q);
            invert_analytic = |p| analytic.uv(p, true);
            (&sample_analytic, &invert_analytic, analytic.periods())
        }
    };
    let trim: Vec<Vec<Uv>> = if whole_band {
        Vec::new()
    } else {
        loops
            .iter()
            .map(|ring| ring.iter().map(|&p| uv(p)).collect())
            .collect()
    };
    let guard = domain.map(|[a, b]| (b - a) * 1e-5);
    // Face sense reverses the output winding, not the carrier's trim domain.
    // Orient the independent chart's mask using one topological interior witness.
    let witness = if trim.is_empty() || !atlas_domains.is_empty() {
        None
    } else {
        catallax_boolean::face_point(&mut model.clone(), face.face, dec)
            .ok()
            .map(|w| uv(xyz(&w.value)))
    };
    let interior_left = witness.map_or(input.sense == Sense::Same, |q| {
        inside(q, &trim, periods, true, [0.0; 2])
    });
    let mut interior = 0;
    let mut worst = (0.0, [0.0; 3]);
    let mut check = |p| {
        if !near_mesh(p, &triangles, bound) {
            let d = triangles
                .iter()
                .map(|&t| triangle_gap(p, t))
                .fold(f64::INFINITY, f64::min)
                .sqrt();
            if d > worst.0 {
                worst = (d, p);
            }
        }
    };
    for &p in loops.iter().flatten() {
        check(p);
    }
    let cells = if atlas_domains.is_empty() {
        std::slice::from_ref(&domain)
    } else {
        &atlas_domains
    };
    for q in cells.iter().flat_map(|&cell| sampled_range(cell)) {
        if inside(q, &trim, periods, interior_left, guard) {
            interior += 1;
            check(at(q));
        }
    }
    // A narrow band can fall entirely between longitude stations on the full
    // cylinder. Add a grid around its trim samples, not around mesh vertices.
    if interior < 16
        && atlas_domains.is_empty()
        && let Some(witness) = witness
    {
        let local = std::array::from_fn(|i| {
            if periods[i] == 0.0 {
                return domain[i];
            }
            let [lo, hi] = trim
                .iter()
                .flatten()
                .map(|q| delta(q[i], witness[i], periods[i]))
                .fold([0.0_f64; 2], |[lo, hi], x| [lo.min(x), hi.max(x)]);
            [witness[i] + lo, witness[i] + hi]
        });
        for mut q in sampled_range(local) {
            for i in 0..2 {
                if periods[i] > 0.0 {
                    q[i] = domain[i][0] + (q[i] - domain[i][0]).rem_euclid(periods[i]);
                }
            }
            if inside(q, &trim, periods, interior_left, guard) {
                interior += 1;
                check(at(q));
            }
        }
    }
    if interior == 0 {
        findings.push(format!("{label}: no independent interior samples"));
    }
    if worst.0 > bound {
        findings.push(format!(
            "{label}: face -> mesh {} mm > {bound} mm at {:?}",
            worst.0, worst.1
        ));
    }

    let mut worst = (0.0, [0.0; 3]);
    let mut model = model.clone();
    let mut membership = catallax_boolean::PreparedFaceQueries::new(&mut model);
    let width = Rational::from_f64_exact((bound * 1e-5).max(1e-13)).unwrap();
    for &[a, b, c] in triangles.iter().step_by((triangles.len() / 128).max(1)) {
        for (s, t) in [(0.25, 0.25), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5)] {
            let p = std::array::from_fn(|i| a[i] + s * (b[i] - a[i]) + t * (c[i] - a[i]));
            let q = uv(p);
            let mut d = if in_atlas(q) && inside(q, &trim, periods, interior_left, [0.0; 2]) {
                gap(p, at(q))
            } else {
                f64::INFINITY
            };
            if d > bound * bound {
                d = d.min(trims.nearest_boundary_sample(p, bound, dec));
            }
            if d > bound * bound && !matches!(surface, Surface::BSpline(_) | Surface::TorusSheet(_))
            {
                // A coarse UV mask can reject a point inside a curved trim.
                // Recheck with an independently constructed exact witness.
                let target = catallax_shapes::RPoint3::new(
                    Rational::from_f64_exact(p[0]).unwrap(),
                    Rational::from_f64_exact(p[1]).unwrap(),
                    Rational::from_f64_exact(p[2]).unwrap(),
                )
                .to_real();
                if let Ok(Some(witness)) =
                    catallax_boolean::surface_witness(surface, &target, &width, dec)
                    && matches!(membership.contains(face.face, &witness, dec), Ok(true))
                {
                    d = d.min(gap(p, xyz(&witness.to_real())));
                }
            }
            if d > worst.0 {
                worst = (d, p);
            }
        }
    }
    if worst.0 > bound * bound {
        findings.push(format!(
            "{label}: mesh -> face {} mm > {bound} mm at {:?}",
            worst.0.sqrt(),
            worst.1
        ));
    }
}

// bad45f069: intersect the two offset cylinders, then sample the rolling ball
// between the two contact directions. The ideal law is not the fitted carrier.
fn check_rolling_ball(case: &str, model: &Model, mesh: &Mesh) -> usize {
    let mut saddles = 0;
    for face in &mesh.faces {
        let Surface::BSpline(SplineSupport::RollingSaddle(saddle)) =
            model.geom().surface(model.face(face.face).surface)
        else {
            continue;
        };
        let law = saddle.definition().law();
        let branches: &[i64] = match law.ring {
            SaddleRing::Through => std::slice::from_ref(&law.branch),
            SaddleRing::Overlap => &[-1, 1],
            SaddleRing::Crossing { .. } => continue,
        };
        saddles += 1;
        let [axis, side, across, along] =
            [law.narrow_axis(), law.side_axis(), law.across, law.along]
                .map(|v| [v.x, v.y, v.z].each_ref().map(number));
        let foot = xyz(&law.foot);
        let [wide, narrow, radius, gap, cosine, sine] = [
            law.wide, law.narrow, law.radius, law.gap, law.cosine, law.sine,
        ]
        .each_ref()
        .map(number);
        let [wide_side, narrow_side] = law.sides.map(|s| s as f64);
        let big = wide + wide_side * radius;
        let small = narrow + narrow_side * radius;
        let triangles = mesh_triangles(mesh, face.triangles.clone());
        let limit = (&face.max_chord_error + &saddle.bound().upper).to_f64_nearest() + 1e-9;
        let mut samples = 0;
        for station in 0..48 {
            let (sin, cos) = (TAU * (f64::from(station) + 0.5) / 48.0).sin_cos();
            let radial: Point = std::array::from_fn(|i| cos * across[i] + sin * side[i]);
            let x = gap + small * cos;
            if x.abs() >= big {
                continue;
            }
            for &branch in branches {
                let y = branch as f64 * (big * big - x * x).sqrt();
                let height = (y - cosine * small * sin) / sine;
                let centre: Point =
                    std::array::from_fn(|i| foot[i] + small * radial[i] + height * axis[i]);
                let wide_ray: Point =
                    std::array::from_fn(|i| -wide_side * (x * across[i] + y * along[i]) / big);
                for t in [0.125, 0.5, 0.875] {
                    let ray = unit(std::array::from_fn(|i| {
                        (1.0 - t) * wide_ray[i] - t * narrow_side * radial[i]
                    }));
                    let p = std::array::from_fn(|i| centre[i] + radius * ray[i]);
                    assert!(
                        near_mesh(p, &triangles, limit),
                        "{case}, saddle face {}, station {station}, branch {branch}, section {t}: exceeds {limit} mm",
                        face.face.index()
                    );
                    samples += 1;
                }
            }
        }
        assert!(samples >= 48, "{case}: too few offset-field samples");
    }
    saddles
}

#[test]
fn loop_boss_stays_near_the_rolling_ball_at_viewport_levels() {
    use serde_json::json;

    let recipe = Recipe::parse(
        &json!({
            "format":"catallax-recipe/3", "name":"loop boss",
            "steps":[
                {"id":"wide","op":"cylinder","radius":"10","height":"24",
                 "pose":{"rotation":["1","1","0","0"],"translation":["0","20","0"]}},
                {"id":"wall","op":"cylinder","radius":"4","height":"40",
                 "pose":{"rotation":["1","0","0","0"],"translation":["8","8","-20"]}},
                {"id":"joined","op":"boolean","kind":"union","target":"wide","tool":"wall"},
                {"id":"round","op":"fillet","body":"joined","radius":"1",
                 "edges":[{"between":[{"step":"wide","role":{"side":{"segment":"wall"}}},
                                      {"step":"wall","role":{"side":{"segment":"wall"}}}]}]}
            ],
            "result":"round"
        })
        .to_string(),
    )
    .unwrap();
    let mut dec = Decider::default();
    let base = catallax_app::build(&recipe, &MeshOptions::default(), &mut dec).unwrap();
    let positions = &base.display.value.positions;
    let diagonal = (0..3)
        .map(|i| {
            let lo = positions.iter().map(|p| p[i]).fold(f64::INFINITY, f64::min);
            let hi = positions
                .iter()
                .map(|p| p[i])
                .fold(f64::NEG_INFINITY, f64::max);
            (hi - lo).powi(2)
        })
        .sum::<f64>()
        .sqrt();
    for level in [4.0, 1.0, 0.25] {
        let options = MeshOptions {
            chord_tolerance: Rational::parse(&format!("{:.12}", 0.002 * diagonal * level)).unwrap()
                / Rational::from_i64(2),
            angular_deviation: Some(Rational::parse("0.2617993877991494").unwrap()),
            ..MeshOptions::default()
        };
        let mesh = mesh_body(&base.built.model, base.built.result, &options, &mut dec).unwrap();
        let case = format!("loop boss at viewport level {level}");
        assert!(
            check_rolling_ball(&case, &base.built.model, &mesh.value) > 0,
            "{case}: no rolling saddle faces"
        );
    }
}

#[test]
fn field_sphere_and_torus_charts_stay_within_display_bounds() {
    use catallax_field::{
        Plan,
        display::{Options, mesh_plan},
    };
    use serde_json::json;

    let fixtures = [
        (
            Shape::Sphere(2.0),
            json!({"op":"sphere","radius":"2"}),
            [-TAU / 4.0, TAU / 4.0],
        ),
        (
            Shape::Torus(3.0, 1.0),
            json!({"op":"torus","major":"3","minor":"1"}),
            [-TAU / 2.0, TAU / 2.0],
        ),
    ];
    for (shape, mut operation, v) in fixtures {
        let kind = operation["op"].as_str().unwrap().to_owned();
        operation["id"] = json!("stock");
        operation["pose"] = json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]});
        let recipe = Recipe::parse(
            &json!({"format":"catallax-recipe/3","name":kind,
            "steps":[operation],"result":"stock"})
            .to_string(),
        )
        .unwrap();
        let plan = Plan::compile(&recipe).unwrap();
        let mesh = mesh_plan(
            &plan,
            &Options {
                tolerance: Rational::from_ratio(1, 20).unwrap(),
                ..Options::default()
            },
            &mut Decider::default(),
        )
        .unwrap();
        let triangles: Vec<_> = mesh
            .triangles
            .iter()
            .map(|t| t.vertices.map(|i| mesh.positions[i]))
            .collect();
        let bound = mesh.max_deviation.to_f64_nearest() + 1e-12;
        let chart = Analytic {
            shape,
            origin: [0.0; 3],
            frame: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        };
        for q in sampled_range([[-TAU / 2.0, TAU / 2.0], v]) {
            let p = chart.at(q);
            assert!(
                near_mesh(p, &triangles, bound),
                "{kind}: source sample {p:?} exceeds {bound} mm"
            );
        }
        for [a, b, c] in triangles {
            for (u, v) in [(0.25, 0.25), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5)] {
                let [x, y, z] =
                    std::array::from_fn(|i| a[i] + u * (b[i] - a[i]) + v * (c[i] - a[i]));
                let distance = match shape {
                    Shape::Sphere(r) => (x.hypot(y).hypot(z) - r).abs(),
                    Shape::Torus(major, minor) => ((x.hypot(y) - major).hypot(z) - minor).abs(),
                    _ => unreachable!(),
                };
                assert!(
                    distance <= bound,
                    "{kind}: facet sample exceeds {bound} mm: {distance} mm"
                );
            }
        }
    }
}

fn recipe_of(
    root: &Path,
    row: &Value,
    case_files: &mut BTreeMap<String, BTreeMap<String, Value>>,
    sources: Option<&BTreeMap<String, std::path::PathBuf>>,
) -> Option<(Recipe, Option<(String, Vec<u8>)>)> {
    let mut source = None;
    let mut recipe: Value = if let Some(path) = row["recipe"].as_str() {
        serde_json::from_slice(&fs::read(root.join(path)).unwrap()).unwrap()
    } else if row["inline"].is_object() {
        row["inline"].clone()
    } else if let Some(path) = row["cases"].as_str() {
        case_files.entry(path.into()).or_insert_with(|| {
            fs::read_to_string(root.join(path))
                .unwrap()
                .lines()
                .filter(|l| !l.trim().is_empty())
                .map(|line| {
                    let mut row: Value = serde_json::from_str(line).unwrap();
                    (row["id"].as_str().unwrap().to_owned(), row["recipe"].take())
                })
                .collect()
        })[row["case"].as_str().unwrap()]
        .clone()
    } else {
        let model = row["import"]["model"].as_str().unwrap();
        let path = sources?
            .get(model)
            .unwrap_or_else(|| panic!("{model}: no selected STEP fixture"));
        let hash = row["import"]["sha256"].as_str().unwrap();
        source = Some((hash.to_owned(), fs::read(path).unwrap()));
        serde_json::json!({"format":"catallax-recipe/3", "name":format!("ABC {model}"),
            "steps":[{"id":"imported","op":"import_step",
                "file":{"name":"part.step","sha256":hash},
                "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],
            "result":"imported"})
    };
    if let Some(steps) = row["append"].as_array() {
        recipe["steps"]
            .as_array_mut()
            .unwrap()
            .extend(steps.iter().cloned());
        recipe["result"] = row["result"].clone();
    }
    if row["format"].is_string() {
        recipe["format"] = row["format"].clone();
    }
    Some((Recipe::parse(&recipe.to_string()).unwrap(), source))
}

#[test]
fn checked_corpus_faces_stay_within_their_reported_bounds() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let corpus: Value =
        serde_json::from_str(include_str!("../../../tools/predictor/corpus.json")).unwrap();
    let mut case_files = BTreeMap::new();
    let mut findings = Vec::new();
    let mut faces = 0;
    let mut private = 0;
    let sources = match (
        std::env::var_os("P5_SOURCES"),
        std::env::var_os("P5_SELECTION"),
    ) {
        (Some(root), Some(path)) => {
            let selection: Vec<(String, String)> =
                serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
            Some(
                selection
                    .into_iter()
                    .map(|(model, file)| (model, Path::new(&root).join(file)))
                    .collect(),
            )
        }
        (None, None) => None,
        _ => panic!("set both P5_SOURCES and P5_SELECTION for the private STEP corpus"),
    };
    let selected = std::env::var("CATALLAX_DISPLAY_CASES").ok();
    for row in corpus["rows"].as_array().unwrap() {
        let case = row["id"].as_str().unwrap();
        if selected
            .as_ref()
            .is_some_and(|names| !names.split(',').any(|name| name == case))
        {
            continue;
        }
        let Some((recipe, source)) = recipe_of(&root, row, &mut case_files, sources.as_ref())
        else {
            private += 1;
            continue;
        };
        let mut files = ImportFiles::new();
        if let Some((hash, bytes)) = &source {
            files.insert(hash.clone(), bytes);
        }
        let (body, _) = with_operation_budget(
            row["max_work"]
                .as_u64()
                .unwrap_or(catallax_app::DEFAULT_BUILD_WORK),
            || {
                build_incremental_with_files_without_display(
                    &recipe,
                    None,
                    &files,
                    &mut Decider::default(),
                )
            },
        );
        let body = match body {
            Ok(EvaluatedApplicationOutcome::Exact(body)) => body,
            Ok(EvaluatedApplicationOutcome::Numerical(_)) => {
                eprintln!("{case}: numerical outcome has no per-face exact display bound");
                continue;
            }
            Err(err) => {
                eprintln!("{case}: build refused: {err}");
                continue;
            }
        };
        let built = body.checkpoint.built();
        let started = Instant::now();
        let (mesh, work) = with_operation_budget(u64::MAX, || {
            mesh_body(
                &built.model,
                built.result,
                &MeshOptions::default(),
                &mut Decider::default(),
            )
        });
        let mesh = match mesh {
            Ok(mesh) => mesh.value,
            Err(err) => {
                match (case, err.code) {
                    ("m10-top-cut" | "fillet-vertex-4-two", "mesh.precision-or-budget") => {
                        eprintln!("{case}: {err} (10.6)")
                    }
                    ("fillet-partial-oblique", "mesh.precision-or-budget") => {
                        eprintln!("{case}: {err} (10.2)")
                    }
                    _ => findings.push(format!("{case}: display refused: {err}")),
                }
                continue;
            }
        };
        eprintln!(
            "{case}: display_ms={:.3} display_work={work:?} triangles={}",
            started.elapsed().as_secs_f64() * 1000.0,
            mesh.triangles.len()
        );
        assert_eq!(
            mesh.faces.iter().map(|f| f.face).collect::<Vec<_>>(),
            built.model.faces_of_body(built.result),
            "{case}"
        );
        let sampled_at = Instant::now();
        let mut dec = Decider::default();
        for face in &mesh.faces {
            check_face(case, &built.model, &mesh, face, &mut dec, &mut findings);
            faces += 1;
        }
        check_rolling_ball(case, &built.model, &mesh);
        eprintln!(
            "{case}: sampling_ms={:.3}",
            sampled_at.elapsed().as_secs_f64() * 1000.0
        );
    }
    eprintln!("{faces} checked faces; {private} private STEP rows not supplied");
    assert!(faces > 0, "corpus did not exercise any faces");
    assert!(
        findings.is_empty(),
        "{}\n{}",
        findings.len(),
        findings.join("\n")
    );
}
