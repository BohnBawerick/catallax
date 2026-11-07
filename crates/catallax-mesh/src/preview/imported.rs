// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: immutable exact source data, explicit rectangle stack,
// caller-owned work/cancellation and vertex limits, existing interval evaluator,
// shared STEP edge samples, and PreviewMesh only. No geometric proof is claimed.
// Rectangle subdivision adapted from rust-brep
// crates/kernel/src/tess_subdivide_mesh.rs, ST_TessSubdivideSrf,
// commit 4c83df59d777ccfcb14f441b1c8d984d1b5d93ab.
//! Unproved preview of imported source geometry. STEP entity names remain pick
//! keys; no fictitious certified Model or SplineSupport is constructed.

use crate::{PreviewEdge, PreviewFace, PreviewMesh, PreviewOptions, Triangle, chart};
use catallax_import::{
    Curve, Face, OrientedEdge, Placement, StepFile, Surface, read_step_with_decider, sha256_hex,
};
use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::{
    BSplineCurve, BSplineSurface, Cone, Cylinder, Plane, Point3, Sphere, Surface as GeomSurface,
    Torus, Vector3,
};
use spade::{ConstrainedDelaunayTriangulation, HasPosition, Point2, Triangulation};
use std::collections::{BTreeMap, BTreeSet};

type P = [f64; 3];
type Uv = [f64; 2];
const BITS: u32 = 128;

/// A source-file preview, not a validated solid or certified display mesh.
#[derive(Clone, Debug)]
pub struct StepPreview {
    /// SHA-256 of the original file, not of rounded preview coordinates.
    pub source_sha256: String,
    /// Retained exact source definitions and file entity identities.
    pub source: StepFile,
    /// Unproved triangles and shared edges keyed by STEP entity name.
    pub mesh: PreviewMesh<String, String>,
}

fn failure(detail: impl Into<String>) -> Refusal {
    Refusal::new("mesh.import-preview-unresolved", detail)
}
fn dot(a: P, b: P) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn sub(a: P, b: P) -> P {
    std::array::from_fn(|i| a[i] - b[i])
}
fn cross(a: P, b: P) -> P {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn norm(a: P) -> f64 {
    dot(a, a).sqrt()
}
fn real_point(p: &[Rational; 3]) -> Point3 {
    Point3::new(
        p[0].clone().into(),
        p[1].clone().into(),
        p[2].clone().into(),
    )
}
fn vector(p: &[Rational; 3]) -> Vector3 {
    Vector3::new(
        p[0].clone().into(),
        p[1].clone().into(),
        p[2].clone().into(),
    )
}
fn point(p: &[Rational; 3]) -> Result<P, Refusal> {
    let p = p.each_ref().map(Rational::to_f64_nearest);
    if p.iter().all(|v| v.is_finite()) {
        Ok(p)
    } else {
        Err(failure("source coordinates exceed binary64 preview range"))
    }
}
fn q(t: f64) -> Result<Rational, Refusal> {
    Rational::from_f64_exact(t).ok_or_else(|| failure("nonfinite preview parameter"))
}
fn parameter(t: f64, domain: [&Rational; 2]) -> Result<Rational, Refusal> {
    let fraction = q(t)?;
    Ok(domain[0] + &(&fraction * &(domain[1] - domain[0])))
}
fn normalized(p: &Placement) -> Result<[P; 3], Refusal> {
    let z = point(&p.axis)?;
    let x = point(&p.ref_dir)?;
    let (nz, nx) = (norm(z), norm(x));
    if nz == 0.0 || nx == 0.0 || !nz.is_finite() || !nx.is_finite() {
        return Err(failure("source frame is not finite"));
    }
    let (x, z) = (x.map(|v| v / nx), z.map(|v| v / nz));
    Ok([x, cross(z, x), z])
}

// Spline charts use normalized parameters, transformed back to exact source
// knots before evaluation. Non-binary decimal endpoints therefore stay inside
// their exact finite domains even when their binary64 rounding points outward.
enum SourceChart {
    Analytic(Box<chart::Chart>),
    Spline(BSplineSurface),
}
impl SourceChart {
    fn new(surface: &Surface, dec: &mut Decider) -> Result<Self, Refusal> {
        let geometry = match surface {
            Surface::BSpline(s) => return Ok(Self::Spline(s.geometry.to_surface(dec)?)),
            Surface::Plane(p) => GeomSurface::Plane(Plane {
                origin: real_point(&p.origin),
                normal: vector(&p.axis),
                ref_dir: vector(&p.ref_dir),
            }),
            Surface::Cylinder(p, r) => GeomSurface::Cylinder(Cylinder {
                origin: real_point(&p.origin),
                axis: vector(&p.axis),
                radius_sq: (r * r).into(),
                ref_dir: vector(&p.ref_dir),
            }),
            Surface::Sphere(p, r) => GeomSurface::Sphere(Sphere {
                center: real_point(&p.origin),
                axis: vector(&p.axis),
                radius: r.clone().into(),
                ref_dir: vector(&p.ref_dir),
            }),
            Surface::Torus(p, a, b) => GeomSurface::Torus(Torus {
                center: real_point(&p.origin),
                axis: vector(&p.axis),
                major: a.clone().into(),
                minor: b.clone().into(),
                ref_dir: vector(&p.ref_dir),
            }),
            Surface::Cone45(p, r) => {
                let axis = vector(&p.axis);
                let distance = &Real::from_rational(r.clone()) / &axis.norm_squared().sqrt();
                GeomSurface::Cone(Cone {
                    apex: real_point(&p.origin).plus(&axis.scaled(&(-distance))),
                    axis,
                    slope: Real::one(),
                    ref_dir: vector(&p.ref_dir),
                })
            }
        };
        Ok(Self::Analytic(Box::new(chart::Chart::new(&geometry, dec)?)))
    }
    fn periodic(&self) -> [bool; 2] {
        match self {
            Self::Analytic(c) => c.periodic,
            Self::Spline(_) => [false; 2],
        }
    }
    fn position(&self, uv: Uv, dec: &Decider) -> Result<P, Refusal> {
        dec.check_deadline()?;
        match self {
            Self::Analytic(c) => c.preview_position(c.canonical(uv), dec),
            Self::Spline(s) => {
                let d = s.domain();
                let u = parameter(uv[0], d[0])?;
                let v = parameter(uv[1], d[1])?;
                chart::position(&s.evaluate_box([&u, &u], [&v, &v], BITS, dec)?)
            }
        }
    }
    fn normal(&self, uv: Uv, dec: &Decider) -> Result<P, Refusal> {
        match self {
            Self::Analytic(c) => c
                .filtered_normal(c.canonical(uv))
                .map(|(n, _)| n.map(|v| v.mid()))
                .ok_or_else(|| failure("preview normal unavailable")),
            Self::Spline(s) => {
                let d = s.domain();
                let u = parameter(uv[0], d[0])?;
                let v = parameter(uv[1], d[1])?;
                let jet = s.jet_box([&u, &u], [&v, &v], BITS, dec)?;
                Ok(cross(
                    chart::position(&jet.first[0])?,
                    chart::position(&jet.first[1])?,
                ))
            }
        }
    }
    fn project(
        &self,
        p: P,
        previous: Option<Uv>,
        tolerance: f64,
        dec: &Decider,
    ) -> Result<Uv, Refusal> {
        let Self::Spline(s) = self else {
            let Self::Analytic(c) = self else {
                unreachable!()
            };
            let uv = c.project(p);
            return if uv.iter().all(|v| v.is_finite()) {
                Ok(uv)
            } else {
                Err(failure("source chart projection is not finite"))
            };
        };
        let domain = s.domain();
        let scale = domain.map(|d| (d[1] - d[0]).to_f64_nearest());
        let seeds = previous
            .into_iter()
            .chain((0..3).flat_map(|i| (0..3).map(move |j| [i as f64 * 0.5, j as f64 * 0.5])));
        for mut uv in seeds {
            for _ in 0..24 {
                dec.check_deadline()?;
                let u = parameter(uv[0], domain[0])?;
                let v = parameter(uv[1], domain[1])?;
                let jet = s.jet_box([&u, &u], [&v, &v], BITS, dec)?;
                let residual = sub(p, chart::position(&jet.value)?);
                if norm(residual) <= tolerance {
                    return Ok(uv);
                }
                let a = chart::position(&jet.first[0])?.map(|x| x * scale[0]);
                let b = chart::position(&jet.first[1])?.map(|x| x * scale[1]);
                let (aa, ab, bb) = (dot(a, a), dot(a, b), dot(b, b));
                let det = aa * bb - ab * ab;
                if !det.is_finite() || det <= 0.0 {
                    break;
                }
                let (ra, rb) = (dot(residual, a), dot(residual, b));
                let step = [(ra * bb - rb * ab) / det, (rb * aa - ra * ab) / det];
                let next = std::array::from_fn(|i| (uv[i] + step[i]).clamp(0.0, 1.0));
                // Parameter identity detects stagnation, not geometric equality.
                if chart::key(next) == chart::key(uv) {
                    break;
                }
                uv = next;
            }
        }
        Err(failure(
            "spline trim projection did not converge within the preview work limit",
        ))
    }
}

enum EdgeCurve {
    Line(P, P),
    Circle {
        center: P,
        axes: [P; 3],
        radii: [f64; 2],
    },
    Spline(BSplineCurve),
}
impl EdgeCurve {
    fn new(edge: &OrientedEdge, dec: &mut Decider) -> Result<Self, Refusal> {
        Ok(match &edge.curve {
            Curve::Line(origin, direction) => {
                let origin = point(origin)?;
                let direction = point(direction)?;
                let length = dot(direction, direction);
                if length == 0.0 || !length.is_finite() {
                    return Err(failure("source line direction is outside preview range"));
                }
                let project = |p: &[Rational; 3]| -> Result<P, Refusal> {
                    let parameter = dot(sub(point(p)?, origin), direction) / length;
                    Ok(std::array::from_fn(|i| {
                        origin[i] + parameter * direction[i]
                    }))
                };
                Self::Line(project(&edge.start.point)?, project(&edge.end.point)?)
            }
            Curve::Circle(p, r) => Self::Circle {
                center: point(&p.origin)?,
                axes: normalized(p)?,
                radii: [r.to_f64_nearest(); 2],
            },
            Curve::Ellipse(p, a, b) => Self::Circle {
                center: point(&p.origin)?,
                axes: normalized(p)?,
                radii: [a.to_f64_nearest(), b.to_f64_nearest()],
            },
            Curve::BSpline(c) => Self::Spline(c.geometry.to_curve(dec)?),
        })
    }
    fn closed(&self) -> bool {
        match self {
            Self::Circle { .. } => true,
            Self::Line(..) => false,
            Self::Spline(s) => {
                s.control(0).rational_components()
                    == s.control(s.control_count() - 1).rational_components()
            }
        }
    }
    fn value(&self, t: f64, dec: &Decider) -> Result<P, Refusal> {
        dec.check_deadline()?;
        let t = if self.closed() { t - t.floor() } else { t };
        match self {
            Self::Line(a, b) => Ok(std::array::from_fn(|i| a[i] + (b[i] - a[i]) * t)),
            Self::Circle {
                center,
                axes,
                radii,
            } => {
                let (s, c) = (t * std::f64::consts::TAU).sin_cos();
                Ok(std::array::from_fn(|i| {
                    center[i] + axes[0][i] * radii[0] * c + axes[1][i] * radii[1] * s
                }))
            }
            Self::Spline(s) => {
                let t = parameter(t, s.domain())?;
                chart::position(&s.evaluate_box([&t, &t], BITS, dec)?)
            }
        }
    }
    fn locate(&self, p: P, tolerance: f64, dec: &Decider) -> Result<f64, Refusal> {
        if let Self::Circle {
            center,
            axes,
            radii,
        } = self
        {
            let d = sub(p, *center);
            return Ok(
                (dot(d, axes[1]) / radii[1]).atan2(dot(d, axes[0]) / radii[0])
                    / std::f64::consts::TAU,
            );
        }
        let Self::Spline(s) = self else {
            return Err(failure("line endpoint parameters are explicit"));
        };
        let mut best: Option<(f64, f64)> = None;
        for i in 0..=16 {
            let t = i as f64 / 16.0;
            let error = norm(sub(self.value(t, dec)?, p));
            if !error.is_finite() {
                return Err(failure("edge endpoint seed is outside preview range"));
            }
            if best.is_none_or(|(e, _)| error < e) {
                best = Some((error, t));
            }
        }
        let mut t = best.expect("seventeen endpoint seeds").1;
        let domain = s.domain();
        let scale = (domain[1] - domain[0]).to_f64_nearest();
        for _ in 0..32 {
            let u = parameter(t, domain)?;
            let jet = s.jet_box([&u, &u], BITS, dec)?;
            let residual = sub(chart::position(&jet.value)?, p);
            if norm(residual) <= tolerance {
                return Ok(t);
            }
            let tangent = chart::position(&jet.first)?.map(|v| v * scale);
            let length = dot(tangent, tangent);
            if length == 0.0 || !length.is_finite() {
                break;
            }
            let next = (t - dot(residual, tangent) / length).clamp(0.0, 1.0);
            if next.to_bits() == t.to_bits() {
                break;
            }
            t = next;
        }
        Err(failure("edge endpoint projection did not converge"))
    }
}

fn push(
    mesh: &mut PreviewMesh<String, String>,
    p: P,
    options: &PreviewOptions,
) -> Result<usize, Refusal> {
    if mesh.positions.len() >= options.max_vertices {
        return Err(failure("import preview reached max_vertices"));
    }
    if !p.iter().all(|x| x.is_finite()) {
        return Err(failure("nonfinite preview sample"));
    }
    let index = mesh.positions.len();
    mesh.positions.push(p);
    Ok(index)
}

struct Segment {
    range: [f64; 2],
    ends: [usize; 2],
    depth: u32,
}

fn edge_segment(
    curve: &EdgeCurve,
    segment: Segment,
    indices: &mut Vec<usize>,
    mesh: &mut PreviewMesh<String, String>,
    options: &PreviewOptions,
    dec: &Decider,
) -> Result<(), Refusal> {
    let Segment { range, ends, depth } = segment;
    dec.check_deadline()?;
    let p = [0.25, 0.5, 0.75].map(|f| curve.value(range[0] + f * (range[1] - range[0]), dec));
    let [a, b, c] = p;
    let p = [a?, b?, c?];
    let a = mesh.positions[ends[0]];
    let b = mesh.positions[ends[1]];
    let target = options.chord_tolerance.to_f64_nearest() * 0.125;
    let small = p.iter().enumerate().all(|(j, p)| {
        let fraction = (j + 1) as f64 * 0.25;
        let chord = std::array::from_fn(|i| a[i] + fraction * (b[i] - a[i]));
        norm(sub(*p, chord)) <= target
    });
    if small {
        indices.push(ends[1]);
        return Ok(());
    }
    if depth >= 24 {
        return Err(failure("edge preview subdivision limit reached"));
    }
    let middle = push(mesh, p[1], options)?;
    let mid = (range[0] + range[1]) * 0.5;
    edge_segment(
        curve,
        Segment {
            range: [range[0], mid],
            ends: [ends[0], middle],
            depth: depth + 1,
        },
        indices,
        mesh,
        options,
        dec,
    )?;
    edge_segment(
        curve,
        Segment {
            range: [mid, range[1]],
            ends: [middle, ends[1]],
            depth: depth + 1,
        },
        indices,
        mesh,
        options,
        dec,
    )
}

#[derive(Clone, Copy)]
struct Node {
    uv: Uv,
    mesh: usize,
}
impl HasPosition for Node {
    type Scalar = f64;
    fn position(&self) -> Point2<f64> {
        Point2::new(self.uv[0], self.uv[1])
    }
}
type Cdt = ConstrainedDelaunayTriangulation<Node>;

fn inside(uv: Uv, loops: &[Vec<Node>]) -> bool {
    let mut inside = false;
    for lp in loops {
        for i in 0..lp.len() {
            let (a, b) = (lp[i].uv, lp[(i + 1) % lp.len()].uv);
            if (a[1] > uv[1]) != (b[1] > uv[1])
                && uv[0] < a[0] + (uv[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
            {
                inside = !inside;
            }
        }
    }
    inside
}
fn insert(cdt: &mut Cdt, node: Node) -> Result<spade::handles::FixedVertexHandle, Refusal> {
    cdt.insert(node)
        .map_err(|e| failure(format!("source chart triangulation: {e}")))
}

// Donor rectangle stack and midpoint bisection, without mutable knot rescaling
// or legacy status/sentinel arrays. Flatness samples are preview heuristics.
fn seed_rectangles(
    chart: &SourceChart,
    bounds: [Uv; 2],
    loops: &[Vec<Node>],
    cdt: &mut Cdt,
    mesh: &mut PreviewMesh<String, String>,
    options: &PreviewOptions,
    dec: &Decider,
) -> Result<(), Refusal> {
    let mut stack = vec![(bounds, 0u32)];
    while let Some(([lo, hi], depth)) = stack.pop() {
        dec.check_deadline()?;
        let mid: Uv = std::array::from_fn(|i| (lo[i] + hi[i]) * 0.5);
        let corners = [lo, [hi[0], lo[1]], hi, [lo[0], hi[1]]].map(|uv| chart.position(uv, dec));
        let [a, b, c, d] = corners;
        let corners = [a?, b?, c?, d?];
        let actual = chart.position(mid, dec)?;
        let linear = std::array::from_fn(|i| corners.iter().map(|p| p[i] * 0.25).sum());
        let mut error = norm(sub(actual, linear));
        for (uv, ends) in [
            ([mid[0], lo[1]], [0, 1]),
            ([hi[0], mid[1]], [1, 2]),
            ([mid[0], hi[1]], [3, 2]),
            ([lo[0], mid[1]], [0, 3]),
        ] {
            let linear = std::array::from_fn(|i| (corners[ends[0]][i] + corners[ends[1]][i]) * 0.5);
            error = error.max(norm(sub(chart.position(uv, dec)?, linear)));
        }
        if inside(mid, loops) {
            let index = push(mesh, actual, options)?;
            insert(
                cdt,
                Node {
                    uv: mid,
                    mesh: index,
                },
            )?;
        }
        if error <= options.chord_tolerance.to_f64_nearest() * 0.5 {
            continue;
        }
        if depth >= 24 || !(lo[0] < mid[0] && mid[0] < hi[0] && lo[1] < mid[1] && mid[1] < hi[1]) {
            return Err(failure("surface preview subdivision limit reached"));
        }
        let axis = usize::from((hi[1] - lo[1]) > (hi[0] - lo[0]));
        let mut a = hi;
        let mut b = lo;
        a[axis] = mid[axis];
        b[axis] = mid[axis];
        stack.push(([b, hi], depth + 1));
        stack.push(([lo, a], depth + 1));
    }
    Ok(())
}

fn face_mesh(
    face: &Face,
    shared: &BTreeMap<&str, usize>,
    mesh: &mut PreviewMesh<String, String>,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let chart = SourceChart::new(&face.surface, dec)?;
    let mut loops = Vec::new();
    let projection_target = options.chord_tolerance.to_f64_nearest() * 0.01;
    for lp in &face.loops {
        let mut nodes = Vec::new();
        for index in 0..lp.edges.len() {
            let edge = &lp.edges[if lp.orientation {
                index
            } else {
                lp.edges.len() - 1 - index
            }];
            let indices = &mesh.edges[shared[edge.edge.as_str()]].vertices;
            let forward = edge.forward == lp.orientation;
            let first = if forward {
                indices[0]
            } else {
                indices[indices.len() - 1]
            };
            if let Some(previous) = nodes.last() {
                let previous: &Node = previous;
                if previous.mesh != first {
                    return Err(failure("source trim edges do not form a connected loop"));
                }
                nodes.pop();
            }
            for offset in 0..indices.len() {
                let index = indices[if forward {
                    offset
                } else {
                    indices.len() - 1 - offset
                }];
                let previous = nodes.last().map(|n: &Node| n.uv);
                let mut uv =
                    chart.project(mesh.positions[index], previous, projection_target, dec)?;
                if let Some(before) = previous {
                    uv = chart::align(uv, before, chart.periodic());
                }
                nodes.push(Node { uv, mesh: index });
            }
        }
        if nodes.len() < 4 || nodes.first().map(|n| n.mesh) != nodes.last().map(|n| n.mesh) {
            return Err(failure("source trim loop is not closed"));
        }
        let last = nodes.pop().unwrap();
        if norm([
            last.uv[0] - nodes[0].uv[0],
            last.uv[1] - nodes[0].uv[1],
            0.0,
        ]) > 1e-8
        {
            return Err(failure(
                "periodic source trim requires an explicit seam chart",
            ));
        }
        loops.push(nodes);
    }
    if loops.is_empty() {
        return Err(failure("source face has no finite trimming loop"));
    }
    let mut bounds = [loops[0][0].uv; 2];
    for node in loops.iter().flatten() {
        let [lo, hi] = &mut bounds;
        for ((lo, hi), p) in lo.iter_mut().zip(hi).zip(node.uv) {
            *lo = (*lo).min(p);
            *hi = (*hi).max(p);
        }
    }
    let mut cdt = Cdt::new();
    for lp in &loops {
        for i in 0..lp.len() {
            dec.check_deadline()?;
            let a = insert(&mut cdt, lp[i])?;
            let b = insert(&mut cdt, lp[(i + 1) % lp.len()])?;
            if a == b || !cdt.can_add_constraint(a, b) {
                return Err(failure("crossing or collapsed source trim constraints"));
            }
            cdt.add_constraint(a, b);
        }
    }
    seed_rectangles(&chart, bounds, &loops, &mut cdt, mesh, options, dec)?;
    let start = mesh.triangles.len();
    for _ in 0..32 {
        let mut pending = BTreeMap::new();
        let mut triangles = Vec::new();
        for triangle in cdt.inner_faces() {
            dec.check_deadline()?;
            let nodes = triangle.vertices().map(|v| *v.data());
            let center = std::array::from_fn(|i| nodes.iter().map(|n| n.uv[i] / 3.0).sum());
            if !inside(center, &loops) {
                continue;
            }
            let p = nodes.map(|n| mesh.positions[n.mesh]);
            let normal = cross(sub(p[1], p[0]), sub(p[2], p[0]));
            if norm(normal) == 0.0 {
                return Err(failure(
                    "source preview triangle collapsed in spatial coordinates",
                ));
            }
            let mut refine = false;
            for weights in [
                [0.5, 0.5, 0.0],
                [0.0, 0.5, 0.5],
                [0.5, 0.0, 0.5],
                [1.0 / 3.0; 3],
            ] {
                let uv = std::array::from_fn(|i| (0..3).map(|j| weights[j] * nodes[j].uv[i]).sum());
                let linear = std::array::from_fn(|i| (0..3).map(|j| weights[j] * p[j][i]).sum());
                refine |= norm(sub(chart.position(uv, dec)?, linear))
                    > options.chord_tolerance.to_f64_nearest();
                if let Some(angle) = &options.angular_deviation {
                    let actual = chart.normal(uv, dec)?;
                    let cosine = dot(normal, actual) / (norm(normal) * norm(actual));
                    refine |= !cosine.is_finite()
                        || cosine.clamp(-1.0, 1.0).acos() > angle.to_f64_nearest();
                }
            }
            if refine {
                pending.insert(chart::key(center), center);
                continue;
            }
            let mut vertices = nodes.map(|n| n.mesh);
            if !face.same_sense {
                vertices.swap(1, 2);
            }
            let first = (0..3).min_by_key(|&i| vertices[i]).unwrap();
            vertices.rotate_left(first);
            triangles.push(Triangle {
                vertices,
                face: face.name.clone(),
            });
        }
        if pending.is_empty() {
            if triangles.is_empty() {
                return Err(failure("source face has no previewable interior"));
            }
            triangles.sort_unstable_by_key(|t| t.vertices);
            mesh.triangles.extend(triangles);
            mesh.faces.push(PreviewFace {
                face: face.name.clone(),
                triangles: start..mesh.triangles.len(),
            });
            return Ok(());
        }
        for uv in pending.into_values() {
            let index = push(mesh, chart.position(uv, dec)?, options)?;
            insert(&mut cdt, Node { uv, mesh: index })?;
        }
    }
    Err(failure("source triangle preview refinement limit reached"))
}

/// Reads exact STEP source data and previews every face with one shared sample
/// table per EDGE_CURVE. Output is explicitly unproved, including projection,
/// trimming, orientation and sampled flatness. Failure returns no partial mesh.
/// It never invokes solid certification or converts to a certified Mesh.
pub fn preview_step(
    bytes: &[u8],
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<StepPreview, Refusal> {
    super::validate(options)?;
    let target = options.chord_tolerance.to_f64_nearest();
    if !target.is_finite() || target <= 0.0 {
        return Err(failure("visual target is outside binary64 preview range"));
    }
    let source = read_step_with_decider(bytes, dec)
        .map_err(|e| failure(format!("STEP source refused: {e:?}")))?;
    let mut mesh = PreviewMesh::default();
    let mut shared = BTreeMap::new();
    let mut vertices = BTreeMap::new();
    let mut names = BTreeSet::new();
    for face in source.solids.iter().flat_map(|s| &s.faces) {
        if !names.insert(&face.name) {
            continue;
        }
        for edge in face.loops.iter().flat_map(|lp| &lp.edges) {
            if shared.contains_key(edge.edge.as_str()) {
                continue;
            }
            let mut ends = [0; 2];
            for (i, v) in [&edge.start, &edge.end].into_iter().enumerate() {
                ends[i] = if let Some(&index) = vertices.get(v.name.as_str()) {
                    index
                } else {
                    let index = push(&mut mesh, point(&v.point)?, options)?;
                    vertices.insert(v.name.as_str(), index);
                    index
                };
            }
            let curve = EdgeCurve::new(edge, dec)?;
            let range = if matches!(curve, EdgeCurve::Line(..)) {
                [0.0, 1.0]
            } else {
                let tolerance = options.chord_tolerance.to_f64_nearest() * 0.01;
                let a = curve.locate(mesh.positions[ends[0]], tolerance, dec)?;
                let mut b = curve.locate(mesh.positions[ends[1]], tolerance, dec)?;
                if curve.closed() {
                    if edge.curve_same_sense {
                        while b <= a {
                            b += 1.0;
                        }
                    } else {
                        while b >= a {
                            b -= 1.0;
                        }
                    }
                } else if edge.curve_same_sense != (b > a) {
                    return Err(failure(
                        "source spline edge sense contradicts its finite domain",
                    ));
                }
                [a, b]
            };
            let mut indices = vec![ends[0]];
            edge_segment(
                &curve,
                Segment {
                    range,
                    ends,
                    depth: 0,
                },
                &mut indices,
                &mut mesh,
                options,
                dec,
            )?;
            shared.insert(edge.edge.as_str(), mesh.edges.len());
            mesh.edges.push(PreviewEdge {
                edge: edge.edge.clone(),
                vertices: indices,
            });
        }
    }
    names.clear();
    for face in source.solids.iter().flat_map(|s| &s.faces) {
        if names.insert(&face.name) {
            face_mesh(face, &shared, &mut mesh, options, dec)?;
        }
    }
    dec.check_deadline()?;
    Ok(StepPreview {
        source_sha256: sha256_hex(bytes),
        source,
        mesh,
    })
}
