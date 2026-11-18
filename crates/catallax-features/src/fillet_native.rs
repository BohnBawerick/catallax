//! Rolling-ball fillets between planes and cylinders, built as native
//! cylinders and tori.
//!
//! A ball rolling against two planes, or against a plane and a cylinder it
//! is square or parallel to, sweeps a cylinder or a torus. This module takes
//! those carriers as they are instead of fitting a spline to them, so a strip
//! costs the same at every radius. Where two strips meet at a vertex it
//! computes the curve they share and trims both there, as the reference
//! kernels do after walking the ball: a common cross-section where the edge
//! runs on tangentially, or the intersection of the two strips at a sharp
//! corner, with the piece of the neighbouring wall the deeper strip exposes.
//!
//! Every position here is a binary64 proposal. The body is assembled beside
//! the exact model, embedded and admitted by the ordinary numerical checks at
//! the stated tolerance; nothing in this file is authority for the result.

use super::*;
use catallax_recipe::{EdgeName, Origin, Role, StepId};
use catallax_shapes::{Circle, Cylinder, Line, PointDef, Torus, Vector3};
use catallax_topology::{
    BodyId, Coedge, Edge, EdgeId, Face, FaceId, Loop, Lump, Model, Shell, Vertex, VertexId,
    numerical::NumericalOptions,
};

type V = [f64; 3];

/// A gap between two constructions of the same point, in millimetres.
const SAME: f64 = 1e-9;
/// Fewest and most spans of a fitted junction curve.
const SPANS: usize = 4;
const MAX_SPANS: usize = 1024;
/// How far a fitted junction curve may stray from the curve it samples, in
/// millimetres: well inside what admission then measures against the faces.
const FIT: f64 = 2e-10;

fn add(a: V, b: V) -> V {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn sub(a: V, b: V) -> V {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn mul(a: V, s: f64) -> V {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: V) -> f64 {
    dot(a, a).sqrt()
}
fn unit(a: V) -> Option<V> {
    let length = norm(a);
    (length.is_finite() && length > 1e-12).then(|| mul(a, 1.0 / length))
}
fn distance(a: V, b: V) -> f64 {
    norm(sub(a, b))
}
/// Whether two constructions of a point agree. A coordinate that is not a
/// number agrees with nothing.
fn same(a: V, b: V) -> bool {
    distance(a, b) < SAME
}

fn number(value: &Real, dec: &Decider) -> Result<f64, Refusal> {
    Ok(scalar(value, dec)?.to_f64_nearest())
}
fn point64(point: &Point3, dec: &Decider) -> Result<V, Refusal> {
    Ok([
        number(&point.x, dec)?,
        number(&point.y, dec)?,
        number(&point.z, dec)?,
    ])
}
fn vector64(vector: &Vector3, dec: &Decider) -> Result<V, Refusal> {
    Ok([
        number(&vector.x, dec)?,
        number(&vector.y, dec)?,
        number(&vector.z, dec)?,
    ])
}

/// A coordinate that is an integer to within rounding is that integer.
fn tidy(value: f64) -> f64 {
    let nearest = value.round();
    if (value - nearest).abs() < 1e-13 {
        nearest
    } else {
        value
    }
}
fn exact_point(point: V) -> Result<Point3, Refusal> {
    Ok(Point3::new(
        Real::from_rational(rational(point[0])?),
        Real::from_rational(rational(point[1])?),
        Real::from_rational(rational(point[2])?),
    ))
}
fn exact_vector(vector: V) -> Result<[Rational; 3], Refusal> {
    Ok([
        rational(tidy(vector[0]))?,
        rational(tidy(vector[1]))?,
        rational(tidy(vector[2]))?,
    ])
}
fn real_vector(vector: [Rational; 3]) -> Vector3 {
    let [x, y, z] = vector;
    Vector3::new(
        Real::from_rational(x),
        Real::from_rational(y),
        Real::from_rational(z),
    )
}
/// `reference` less its part along `axis`, exactly: a seam direction has to
/// be perpendicular to its axis, and rounded binary64 directions are not.
fn perpendicular(axis: &[Rational; 3], reference: V) -> Result<Vector3, Refusal> {
    let reference = exact_vector(reference)?;
    let along = (0..3).fold(q(0), |sum, i| sum + &reference[i] * &axis[i]);
    let length = (0..3).fold(q(0), |sum, i| sum + &axis[i] * &axis[i]);
    let share = along / length;
    Ok(real_vector(std::array::from_fn(|i| {
        &reference[i] - &(&share * &axis[i])
    })))
}

/// An input face as the ball sees it: its surface and the side material is on.
#[derive(Clone, Copy, Debug)]
enum Support {
    /// `normal` is the unit normal pointing out of the material.
    Plane { point: V, normal: V },
    /// `outward` is 1 where material is inside the cylinder and -1 in a bore.
    Cylinder {
        point: V,
        axis: V,
        radius: f64,
        outward: f64,
    },
}

impl Support {
    fn read(model: &Model, face: FaceId, dec: &Decider) -> Result<Option<Self>, Refusal> {
        let definition = model.face(face);
        let side = if definition.sense == Sense::Same {
            1.0
        } else {
            -1.0
        };
        Ok(match model.geom().surface(definition.surface) {
            Surface::Plane(plane) => {
                let point = point64(&plane.origin, dec)?;
                unit(vector64(&plane.normal, dec)?).map(|normal| Self::Plane {
                    point,
                    normal: mul(normal, side),
                })
            }
            Surface::Cylinder(cylinder) => {
                let point = point64(&cylinder.origin, dec)?;
                let radius = number(&cylinder.radius_sq, dec)?.sqrt();
                unit(vector64(&cylinder.axis, dec)?).map(|axis| Self::Cylinder {
                    point,
                    axis,
                    radius,
                    outward: side,
                })
            }
            _ => None,
        })
    }

    /// Signed distance of `p` from the surface, positive outside the material.
    fn signed(&self, p: V) -> f64 {
        match *self {
            Self::Plane { point, normal } => dot(sub(p, point), normal),
            Self::Cylinder {
                point,
                axis,
                radius,
                outward,
            } => {
                let offset = sub(p, point);
                let radial = sub(offset, mul(axis, dot(offset, axis)));
                outward * (norm(radial) - radius)
            }
        }
    }

    /// Unit normal out of the material at the surface point nearest `p`.
    fn normal(&self, p: V) -> Option<V> {
        match *self {
            Self::Plane { normal, .. } => Some(normal),
            Self::Cylinder {
                point,
                axis,
                outward,
                ..
            } => {
                let offset = sub(p, point);
                unit(sub(offset, mul(axis, dot(offset, axis)))).map(|radial| mul(radial, outward))
            }
        }
    }

    /// Normal curvature along the unit tangent `along`, positive where the
    /// surface bends toward its outward normal.
    fn curvature(&self, along: V) -> f64 {
        match *self {
            Self::Plane { .. } => 0.0,
            Self::Cylinder {
                axis,
                radius,
                outward,
                ..
            } => {
                let lengthwise = dot(along, axis);
                -outward * (1.0 - lengthwise * lengthwise) / radius
            }
        }
    }

    /// The surface point nearest `p`.
    fn foot(&self, p: V) -> Option<V> {
        Some(sub(p, mul(self.normal(p)?, self.signed(p))))
    }
}

/// The path of the ball's centre.
#[derive(Clone, Copy, Debug)]
enum Spine {
    /// `origin + t dir`, with a unit `dir`.
    Line { origin: V, dir: V },
    /// `centre + radius (cos t x + sin t y)`, with `x`, `y`, `axis` orthonormal
    /// and right-handed.
    Circle {
        centre: V,
        axis: V,
        radius: f64,
        x: V,
        y: V,
    },
}

impl Spine {
    fn point(&self, t: f64) -> V {
        match *self {
            Self::Line { origin, dir } => add(origin, mul(dir, t)),
            Self::Circle {
                centre,
                radius,
                x,
                y,
                ..
            } => add(
                centre,
                add(mul(x, radius * t.cos()), mul(y, radius * t.sin())),
            ),
        }
    }

    /// Unit direction of increasing parameter.
    fn tangent(&self, t: f64) -> V {
        match *self {
            Self::Line { dir, .. } => dir,
            Self::Circle { x, y, .. } => add(mul(x, -t.sin()), mul(y, t.cos())),
        }
    }

    /// Millimetres travelled per unit of parameter.
    fn speed(&self) -> f64 {
        match *self {
            Self::Line { .. } => 1.0,
            Self::Circle { radius, .. } => radius,
        }
    }

    /// Parameter of the spine point nearest `p`; on a circle, the turn
    /// nearest `near`.
    fn parameter(&self, p: V, near: f64) -> f64 {
        match *self {
            Self::Line { origin, dir } => dot(sub(p, origin), dir),
            Self::Circle { centre, x, y, .. } => {
                let offset = sub(p, centre);
                let angle = dot(offset, y).atan2(dot(offset, x));
                angle + ((near - angle) / std::f64::consts::TAU).round() * std::f64::consts::TAU
            }
        }
    }

    /// The spine point nearest `p`.
    fn nearest(&self, p: V) -> Option<V> {
        match *self {
            Self::Line { origin, dir } => Some(add(origin, mul(dir, dot(sub(p, origin), dir)))),
            Self::Circle {
                centre,
                axis,
                radius,
                ..
            } => {
                let offset = sub(p, centre);
                let radial = unit(sub(offset, mul(axis, dot(offset, axis))))?;
                Some(add(centre, mul(radial, radius)))
            }
        }
    }

    /// How a point carried rigidly along the path moves per unit of parameter.
    fn carried(&self, p: V) -> V {
        match *self {
            Self::Line { dir, .. } => dir,
            Self::Circle { centre, axis, .. } => cross(axis, sub(p, centre)),
        }
    }

    fn distance(&self, p: V) -> f64 {
        match *self {
            Self::Line { origin, dir } => {
                let offset = sub(p, origin);
                norm(sub(offset, mul(dir, dot(offset, dir))))
            }
            Self::Circle {
                centre,
                axis,
                radius,
                ..
            } => {
                let offset = sub(p, centre);
                let height = dot(offset, axis);
                let radial = norm(sub(offset, mul(axis, height)));
                (radial - radius).hypot(height)
            }
        }
    }
}

/// One rolled edge: its supports in input-face order, the centre path and the
/// path parameters under the edge's start and end vertices.
struct Strip {
    edge: EdgeId,
    faces: [FaceId; 2],
    supports: [Support; 2],
    spine: Spine,
    radius: f64,
    convex: bool,
    range: [f64; 2],
    name: EdgeName,
}

impl Strip {
    fn contact(&self, side: usize, t: f64) -> Option<V> {
        self.supports[side].foot(self.spine.point(t))
    }

    /// Centre, the unit direction to the first contact, the unit direction
    /// square to it in the section plane toward the second, and the angle
    /// between the two contacts.
    fn section(&self, t: f64) -> Option<(V, V, V, f64)> {
        let centre = self.spine.point(t);
        let first = unit(sub(self.contact(0, t)?, centre))?;
        let second = unit(sub(self.contact(1, t)?, centre))?;
        let cosine = dot(first, second).clamp(-1.0, 1.0);
        let across = unit(sub(second, mul(first, cosine)))?;
        Some((centre, first, across, cosine.acos()))
    }

    fn point(&self, t: f64, angle: f64) -> Option<V> {
        let (centre, first, across, _) = self.section(t)?;
        Some(add(
            centre,
            add(
                mul(first, self.radius * angle.cos()),
                mul(across, self.radius * angle.sin()),
            ),
        ))
    }

    /// How the strip point at `angle` moves per radian round its section.
    fn around(&self, t: f64, angle: f64) -> Option<V> {
        let (_, first, across, _) = self.section(t)?;
        Some(add(
            mul(first, -self.radius * angle.sin()),
            mul(across, self.radius * angle.cos()),
        ))
    }

    /// The motion of the strip point at `(t, angle)` per radian of section
    /// angle, as it stays on the zero set of a function whose gradient there
    /// is `gradient`.
    fn sliding(&self, t: f64, angle: f64, gradient: V) -> Option<V> {
        let around = self.around(t, angle)?;
        let along = self.spine.carried(self.point(t, angle)?);
        let rate = -dot(gradient, around) / dot(gradient, along);
        rate.is_finite().then(|| add(around, mul(along, rate)))
    }

    /// Path parameter and section angle of a point on the strip.
    fn locate(&self, p: V, near: f64) -> Option<(f64, f64)> {
        let t = self.spine.parameter(p, near);
        let (centre, first, across, _) = self.section(t)?;
        let offset = sub(p, centre);
        Some((t, dot(offset, across).atan2(dot(offset, first))))
    }

    /// Section angle of the contact with support `side`.
    fn side_angle(&self, side: usize, t: f64) -> Option<f64> {
        Some(if side == 0 { 0.0 } else { self.section(t)?.3 })
    }
}

struct Fan {
    faces: BTreeMap<EdgeId, Vec<(FaceId, Sense)>>,
    edges: BTreeMap<VertexId, Vec<EdgeId>>,
    positions: BTreeMap<VertexId, V>,
}

impl Fan {
    fn read(model: &Model, body: BodyId, dec: &Decider) -> Result<Option<Self>, Refusal> {
        let mut faces: BTreeMap<EdgeId, Vec<(FaceId, Sense)>> = BTreeMap::new();
        for face in model.faces_of_body(body) {
            for boundary in &model.face(face).loops {
                let Loop::Coedges(uses) = model.loop_(*boundary) else {
                    return Ok(None);
                };
                for coedge in uses {
                    let coedge = model.coedge(*coedge);
                    faces
                        .entry(coedge.edge)
                        .or_default()
                        .push((face, coedge.sense));
                }
            }
        }
        let mut edges: BTreeMap<VertexId, Vec<EdgeId>> = BTreeMap::new();
        let mut positions = BTreeMap::new();
        for edge in model.edges_of_body(body) {
            let definition = model.edge(edge);
            for vertex in [definition.start, definition.end].into_iter().flatten() {
                edges.entry(vertex).or_default().push(edge);
                if let std::collections::btree_map::Entry::Vacant(entry) = positions.entry(vertex) {
                    entry.insert(point64(
                        &model.geom().position(model.vertex(vertex).point),
                        dec,
                    )?);
                }
            }
        }
        Ok(Some(Self {
            faces,
            edges,
            positions,
        }))
    }
}

/// Whether the two supports meet in a straight edge.
fn straight(supports: &[Support; 2]) -> bool {
    match supports {
        [Support::Plane { .. }, Support::Plane { .. }] => true,
        [
            Support::Plane { normal, .. },
            Support::Cylinder { axis, .. },
        ]
        | [
            Support::Cylinder { axis, .. },
            Support::Plane { normal, .. },
        ] => dot(*normal, *axis).abs() < 1e-12,
        _ => false,
    }
}

/// Unit direction of travel from an edge's start to its end, at one of them.
fn travel(
    model: &Model,
    fan: &Fan,
    edge: EdgeId,
    supports: &[Support; 2],
    at_end: bool,
    dec: &Decider,
) -> Result<Option<V>, Refusal> {
    let definition = model.edge(edge);
    let (Some(start), Some(end)) = (definition.start, definition.end) else {
        return Ok(None);
    };
    let ends = [fan.positions[&start], fan.positions[&end]];
    if straight(supports) || matches!(model.geom().curve(definition.curve), Curve::Line(_)) {
        return Ok(unit(sub(ends[1], ends[0])));
    }
    let Curve::Circle(circle) = model.geom().curve(definition.curve) else {
        return Ok(None);
    };
    let Some(normal) = unit(vector64(&circle.normal, dec)?) else {
        return Ok(None);
    };
    let radial = sub(ends[usize::from(at_end)], point64(&circle.center, dec)?);
    let turn = if definition.sense == Sense::Same {
        1.0
    } else {
        -1.0
    };
    Ok(unit(mul(cross(normal, radial), turn)))
}

/// Whether the edge between two faces is convex, from their outward normals
/// and the direction the first face's boundary takes along it. `None` for a
/// smooth edge.
fn convexity(normals: [V; 2], travel: V, first_use: Sense) -> Option<bool> {
    let turn = dot(cross(normals[0], normals[1]), travel);
    (turn.abs() > 1e-9).then_some((turn > 0.0) == (first_use == Sense::Same))
}

fn strip(
    model: &Model,
    body: BodyId,
    fan: &Fan,
    edge: EdgeId,
    radius: f64,
    dec: &Decider,
) -> Result<Option<Strip>, Refusal> {
    let Some([(first, first_use), (second, _)]) = fan.faces.get(&edge).map(Vec::as_slice) else {
        return Ok(None);
    };
    let faces = [*first, *second];
    let (Some(a), Some(b)) = (
        Support::read(model, faces[0], dec)?,
        Support::read(model, faces[1], dec)?,
    ) else {
        return Ok(None);
    };
    let supports = [a, b];
    let definition = model.edge(edge);
    let (Some(start), Some(end)) = (definition.start, definition.end) else {
        return Ok(None);
    };
    let ends = [fan.positions[&start], fan.positions[&end]];
    let Some(along) = travel(model, fan, edge, &supports, false, dec)? else {
        return Ok(None);
    };
    let (Some(n0), Some(n1)) = (supports[0].normal(ends[0]), supports[1].normal(ends[0])) else {
        return Ok(None);
    };
    let Some(convex) = convexity([n0, n1], along, *first_use) else {
        return Ok(None);
    };
    // The centre sits one radius inside the material of a convex edge and one
    // radius outside it at a concave one.
    let depth = if convex { -radius } else { radius };
    let middle = mul(add(ends[0], ends[1]), 0.5);
    let spine = match supports {
        [Support::Plane { .. }, Support::Plane { .. }] => {
            let Some(dir) = unit(cross(n0, n1)) else {
                return Ok(None);
            };
            let share = depth / (1.0 + dot(n0, n1));
            Spine::Line {
                origin: add(middle, mul(add(n0, n1), share)),
                dir,
            }
        }
        [
            Support::Plane { point, normal },
            Support::Cylinder {
                point: on_axis,
                axis,
                radius: cylinder,
                outward,
            },
        ]
        | [
            Support::Cylinder {
                point: on_axis,
                axis,
                radius: cylinder,
                outward,
            },
            Support::Plane { point, normal },
        ] => {
            let offset = cylinder + outward * depth;
            if offset <= 1e-9 {
                return Ok(None);
            }
            let lean = dot(normal, axis);
            if lean.abs() > 1.0 - 1e-12 {
                let height = (depth - dot(sub(on_axis, point), normal)) / lean;
                let centre = add(on_axis, mul(axis, height));
                let Some(x) = unit(sub(
                    sub(ends[0], centre),
                    mul(axis, dot(sub(ends[0], centre), axis)),
                )) else {
                    return Ok(None);
                };
                Spine::Circle {
                    centre,
                    axis,
                    radius: offset,
                    x,
                    y: cross(axis, x),
                }
            } else if lean.abs() < 1e-12 {
                // The axis stands this far outside the plane the centre moves in.
                let height = dot(sub(on_axis, point), normal) - depth;
                let reach = offset * offset - height * height;
                if reach <= 0.0 {
                    return Ok(None);
                }
                let side = cross(axis, normal);
                let foot = sub(on_axis, mul(normal, height));
                let candidates = [1.0, -1.0].map(|sign| add(foot, mul(side, sign * reach.sqrt())));
                let gap = |origin: V| Spine::Line { origin, dir: axis }.distance(middle);
                let origin = if gap(candidates[0]) <= gap(candidates[1]) {
                    candidates[0]
                } else {
                    candidates[1]
                };
                Spine::Line { origin, dir: axis }
            } else {
                return Ok(None);
            }
        }
        _ => return Ok(None),
    };
    let first_parameter = spine.parameter(ends[0], 0.0);
    let range = match spine {
        Spine::Line { .. } => [first_parameter, spine.parameter(ends[1], 0.0)],
        Spine::Circle { axis, .. } => {
            // The edge's own circle says which way round it runs.
            let Curve::Circle(circle) = model.geom().curve(definition.curve) else {
                return Ok(None);
            };
            let Some(normal) = unit(vector64(&circle.normal, dec)?) else {
                return Ok(None);
            };
            let turn = dot(normal, axis).signum()
                * if definition.sense == Sense::Same {
                    1.0
                } else {
                    -1.0
                };
            let raw = spine.parameter(ends[1], first_parameter) - first_parameter;
            let swept = if turn > 0.0 {
                raw.rem_euclid(std::f64::consts::TAU)
            } else {
                -(-raw).rem_euclid(std::f64::consts::TAU)
            };
            [first_parameter, first_parameter + swept]
        }
    };
    if (range[1] - range[0]).abs() * spine.speed() < 1e-6 {
        return Ok(None);
    }
    let Some(name) = model.edge_name(body, edge) else {
        return Ok(None);
    };
    let strip = Strip {
        edge,
        faces,
        supports,
        spine,
        radius,
        convex,
        range,
        name,
    };
    // Both contacts and a quarter turn or so between them, or no strip.
    for t in [range[0], 0.5 * (range[0] + range[1]), range[1]] {
        let Some((centre, _, _, angle)) = strip.section(t) else {
            return Ok(None);
        };
        let touching = (0..2).all(|side| {
            strip
                .contact(side, t)
                .is_some_and(|contact| (distance(contact, centre) - radius).abs() < SAME)
        });
        if !touching || !(1e-3..std::f64::consts::PI - 1e-3).contains(&angle) {
            return Ok(None);
        }
    }
    Ok(Some(strip))
}

/// A root of `f` between `inside`, where it is not positive, and the first
/// sample beyond it, stepping by `step`, where it is.
fn first_exit(f: impl Fn(f64) -> Option<f64>, inside: f64, step: f64) -> Option<f64> {
    let mut low = inside;
    let mut high = None;
    for k in 1..=160 {
        let t = inside + step * k as f64;
        if f(t)? > 0.0 {
            high = Some(t);
            break;
        }
        low = t;
    }
    let mut high = high?;
    for _ in 0..80 {
        let middle = 0.5 * (low + high);
        if f(middle)? > 0.0 {
            high = middle;
        } else {
            low = middle;
        }
    }
    Some(0.5 * (low + high))
}

/// A root of `f` near `guess`, bracketed by widening steps.
fn root_near(f: impl Fn(f64) -> Option<f64>, guess: f64, step: f64) -> Option<f64> {
    let at = f(guess)?;
    if at.abs() <= 0.0 {
        return Some(guess);
    }
    let mut bracket = None;
    'search: for k in 1..=160 {
        for sign in [1.0, -1.0] {
            let t = guess + sign * step * k as f64;
            if let Some(value) = f(t)
                && (value > 0.0) != (at > 0.0)
            {
                bracket = Some((guess + sign * step * (k - 1) as f64, t));
                break 'search;
            }
        }
    }
    let (mut same, mut other) = bracket?;
    for _ in 0..80 {
        let middle = 0.5 * (same + other);
        if (f(middle)? > 0.0) == (at > 0.0) {
            same = middle;
        } else {
            other = middle;
        }
    }
    Some(0.5 * (same + other))
}

/// Position, first and second derivative of a junction curve at its knots.
type Samples = (Vec<V>, Vec<V>, Vec<V>);

/// Samples a curve over 0..1, given with its derivative, for a C2 quintic.
/// Spans are halved until the quintic is within [`FIT`] of the curve at the
/// middle of every span, where a quintic through exact ends, slopes and
/// curvatures is furthest off.
fn sampled(f: impl Fn(f64) -> Option<(V, V)>) -> Option<Samples> {
    // Second derivatives are fourth-order differences of the given first.
    // The step is wide because the first derivative of a curve that starts
    // where its two surfaces touch is only known roughly close to that start.
    let h = 1.0 / 256.0;
    let one_sided = |s: f64, here: V, direction: f64| -> Option<V> {
        let mut sum = mul(here, -25.0);
        for (k, weight) in [48.0, -36.0, 16.0, -3.0].into_iter().enumerate() {
            sum = add(sum, mul(f(s + direction * h * (k + 1) as f64)?.1, weight));
        }
        Some(mul(sum, direction / (12.0 * h)))
    };
    let bend = |s: f64, slope: V| -> Option<V> {
        if s < 2.0 * h {
            one_sided(s, slope, 1.0)
        } else if s > 1.0 - 2.0 * h {
            one_sided(s, slope, -1.0)
        } else {
            let outer = sub(f(s - 2.0 * h)?.1, f(s + 2.0 * h)?.1);
            let inner = sub(f(s + h)?.1, f(s - h)?.1);
            Some(mul(add(outer, mul(inner, 8.0)), 1.0 / (12.0 * h)))
        }
    };
    let mut spans = SPANS;
    loop {
        let width = 1.0 / spans as f64;
        let mut points = Vec::with_capacity(spans + 1);
        let mut slopes = Vec::with_capacity(spans + 1);
        let mut bends = Vec::with_capacity(spans + 1);
        for i in 0..=spans {
            let s = i as f64 * width;
            let (point, slope) = f(s)?;
            bends.push(bend(s, slope)?);
            points.push(point);
            slopes.push(slope);
        }
        let mut worst = 0.0f64;
        for i in 0..spans {
            let middle = add(
                add(
                    mul(add(points[i], points[i + 1]), 0.5),
                    mul(sub(slopes[i], slopes[i + 1]), 5.0 * width / 32.0),
                ),
                mul(add(bends[i], bends[i + 1]), width * width / 64.0),
            );
            worst = worst.max(distance(middle, f((i as f64 + 0.5) * width)?.0));
        }
        if worst < FIT {
            return Some((points, slopes, bends));
        }
        if spans >= MAX_SPANS || !worst.is_finite() {
            return None;
        }
        spans *= 2;
    }
}

/// The same curve over 0..1 by its share of arc length, at one speed, so
/// uniform spans of a fitted curve, and the uniform samples the embedding
/// takes of it, land where the curve is however fast its own parameter runs.
fn by_length(f: impl Fn(f64) -> Option<(V, V)>) -> Option<impl Fn(f64) -> Option<(V, V)>> {
    // Three-point Gauss-Legendre length of the piece from `a` to `b`.
    fn piece(f: &impl Fn(f64) -> Option<(V, V)>, a: f64, b: f64) -> Option<f64> {
        let (middle, half) = (0.5 * (a + b), 0.5 * (b - a));
        let offset = half * 0.6f64.sqrt();
        let speed = |s: f64| f(s).map(|(_, slope)| norm(slope));
        Some(
            half * (5.0 * (speed(middle - offset)? + speed(middle + offset)?)
                + 8.0 * speed(middle)?)
                / 9.0,
        )
    }
    // Pieces are halved until halving no longer changes their length.
    const PANELS: usize = 16;
    const DEPTH: usize = 40;
    let mut pending: Vec<(f64, f64, f64, usize)> = Vec::new();
    for i in (0..PANELS).rev() {
        let (a, b) = (i as f64 / PANELS as f64, (i + 1) as f64 / PANELS as f64);
        pending.push((a, b, piece(&f, a, b)?, 0));
    }
    // The start of each piece and the length before it, with the end last.
    let mut table = vec![(0.0, 0.0)];
    let mut total = 0.0;
    while let Some((a, b, whole, depth)) = pending.pop() {
        let middle = 0.5 * (a + b);
        let (left, right) = (piece(&f, a, middle)?, piece(&f, middle, b)?);
        if !(left + right).is_finite() {
            return None;
        }
        if (left + right - whole).abs() <= 1e-13 * (left + right) || depth >= DEPTH {
            total += left + right;
            table.push((b, total));
        } else {
            pending.push((middle, b, right, depth + 1));
            pending.push((a, middle, left, depth + 1));
        }
    }
    if total.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
        return None;
    }
    Some(move |share: f64| -> Option<(V, V)> {
        let scaled = |s: f64| f(s).map(|(point, slope)| (point, mul(slope, total / norm(slope))));
        if share <= 0.0 {
            return scaled(0.0);
        }
        if share >= 1.0 {
            return scaled(1.0);
        }
        let target = share * total;
        let after = table
            .partition_point(|(_, before)| *before <= target)
            .clamp(1, table.len() - 1);
        let ((a, before), (b, _)) = (table[after - 1], table[after]);
        // Newton on the length from the piece's start, kept inside the piece.
        let mut s = 0.5 * (a + b);
        for _ in 0..12 {
            let speed = norm(f(s)?.1);
            let step = (piece(&f, a, s)? - (target - before)) / speed;
            s = (s - step).clamp(a, b);
            if step.abs() <= 1e-15 {
                break;
            }
        }
        scaled(s)
    })
}

/// The C2 quintic through the samples over uniform spans of 0..1, its ends
/// exactly the first and last point.
fn quintic((points, slopes, bends): &Samples, dec: &mut Decider) -> Result<BSplineCurve, Refusal> {
    let spans = points.len() - 1;
    let width = 1.0 / spans as f64;
    let (fifth, twentieth) = (width / 5.0, width * width / 20.0);
    // Bezier points two and three of each span, and at each inner knot the
    // one control the two spans beside it share.
    let leaving = |i: usize| {
        add(
            add(points[i], mul(slopes[i], 2.0 * fifth)),
            mul(bends[i], twentieth),
        )
    };
    let arriving = |i: usize| {
        add(
            sub(points[i], mul(slopes[i], 2.0 * fifth)),
            mul(bends[i], twentieth),
        )
    };
    let mut controls = Vec::with_capacity(3 * spans + 3);
    controls.push(exact_point(points[0])?);
    controls.push(exact_point(add(points[0], mul(slopes[0], fifth)))?);
    for i in 0..spans {
        controls.push(exact_point(leaving(i))?);
        controls.push(exact_point(arriving(i + 1))?);
        if i + 1 < spans {
            controls.push(exact_point(sub(
                points[i + 1],
                mul(bends[i + 1], twentieth),
            ))?);
        }
    }
    controls.push(exact_point(sub(points[spans], mul(slopes[spans], fifth)))?);
    controls.push(exact_point(points[spans])?);
    let mut knots = vec![q(0); 6];
    for i in 1..spans {
        let knot = Rational::from_ratio(i as i64, spans as i64).expect("positive span count");
        knots.extend([knot.clone(), knot.clone(), knot]);
    }
    knots.extend(vec![q(1); 6]);
    let weights = vec![Real::one(); controls.len()];
    BSplineCurve::new(5, knots, controls, weights, dec)
}

/// A curve two faces share at a junction.
enum Shared {
    /// A cross-section of the ball: its centre, the unit normal of its plane.
    Section { centre: V, normal: V },
    /// A sampled curve over 0..1.
    Fitted(Samples),
}

/// How two strip ends close at one vertex.
struct Junction {
    vertex: VertexId,
    /// The two strips and which end of each is here.
    ends: [(usize, usize); 2],
    /// The unselected third edge of the vertex.
    third: EdgeId,
    /// Where both strips touch their common support.
    common: V,
    /// Where the third edge now ends.
    on_third: V,
    /// The curve from `common` that both strips share.
    shared: Shared,
    /// At unequal reach: which of the two strips' own rail the shared curve
    /// ends on, that end point, and the curve from it to `on_third` across
    /// the other strip, in that strip's wall.
    step: Option<(usize, V, Shared)>,
}

fn other(strip: &Strip, face: FaceId) -> usize {
    usize::from(strip.faces[0] == face)
}

fn junction(
    model: &Model,
    fan: &Fan,
    strips: &[Strip],
    vertex: VertexId,
    ends: [(usize, usize); 2],
    selected: &BTreeSet<EdgeId>,
    dec: &Decider,
) -> Result<Option<Junction>, Refusal> {
    let incident = &fan.edges[&vertex];
    let [third] = incident
        .iter()
        .copied()
        .filter(|edge| !selected.contains(edge))
        .collect::<Vec<_>>()[..]
    else {
        return Ok(None);
    };
    if incident.len() != 3 {
        return Ok(None);
    }
    let [a, b] = ends.map(|(strip, _)| &strips[strip]);
    let Some(common) = a.faces.into_iter().find(|face| b.faces.contains(face)) else {
        return Ok(None);
    };
    // Sides of the common support and of each strip's own wall.
    let sides = [a, b].map(|strip| 1 - other(strip, common));
    let walls = [a, b].map(|strip| other(strip, common));
    let wall_faces = [a.faces[walls[0]], b.faces[walls[1]]];
    let Some([(third_first, third_use), (third_second, _)]) =
        fan.faces.get(&third).map(Vec::as_slice)
    else {
        return Ok(None);
    };
    if ![[*third_first, *third_second], [*third_second, *third_first]].contains(&wall_faces) {
        return Ok(None);
    }
    let at = fan.positions[&vertex];
    let parameters = [0, 1].map(|i| [a, b][i].range[ends[i].1]);
    // The direction along each path that leaves the edge past this vertex.
    let out = [0, 1].map(|i| {
        let strip = [a, b][i];
        (strip.range[ends[i].1] - strip.range[1 - ends[i].1]).signum()
    });
    let (Some(wall_a), Some(wall_b)) = (
        a.supports[walls[0]].normal(at),
        b.supports[walls[1]].normal(at),
    ) else {
        return Ok(None);
    };
    if a.convex != b.convex {
        return Ok(None);
    }
    if dot(wall_a, wall_b) > 1.0 - 1e-12 {
        // The edge runs on tangentially: one cross-section serves both strips.
        let centre = a.spine.point(parameters[0]);
        if !same(centre, b.spine.point(parameters[1])) {
            return Ok(None);
        }
        let (Some(on_common), Some(on_third)) = (
            a.contact(sides[0], parameters[0]),
            a.contact(walls[0], parameters[0]),
        ) else {
            return Ok(None);
        };
        let agree = b
            .contact(sides[1], parameters[1])
            .is_some_and(|p| distance(p, on_common) < SAME)
            && b.contact(walls[1], parameters[1])
                .is_some_and(|p| distance(p, on_third) < SAME);
        if !agree {
            return Ok(None);
        }
        return Ok(Some(Junction {
            vertex,
            ends,
            third,
            common: on_common,
            on_third,
            shared: Shared::Section {
                centre,
                normal: a.spine.tangent(parameters[0]),
            },
            step: None,
        }));
    }
    // A sharp third edge. Two convex strips run on past the vertex only where
    // that edge is concave, and meet in the curve their surfaces share.
    let third_supports = if *third_first == wall_faces[0] {
        [a.supports[walls[0]], b.supports[walls[1]]]
    } else {
        [b.supports[walls[1]], a.supports[walls[0]]]
    };
    let third_definition = model.edge(third);
    let Some(third_travel) = travel(
        model,
        fan,
        third,
        &third_supports,
        third_definition.end == Some(vertex),
        dec,
    )?
    else {
        return Ok(None);
    };
    let third_normals = if *third_first == wall_faces[0] {
        [wall_a, wall_b]
    } else {
        [wall_b, wall_a]
    };
    if !a.convex || convexity(third_normals, third_travel, *third_use) != Some(false) {
        return Ok(None);
    }
    // Where the two centre paths cross: the ball touching all three supports.
    let mut t = parameters;
    for _ in 0..60 {
        let gap = sub(a.spine.point(t[0]), b.spine.point(t[1]));
        let (ta, tb) = (
            mul(a.spine.tangent(t[0]), a.spine.speed()),
            mul(b.spine.tangent(t[1]), b.spine.speed()),
        );
        let (aa, ab, bb) = (dot(ta, ta), dot(ta, tb), dot(tb, tb));
        let determinant = aa * bb - ab * ab;
        if determinant.abs() < 1e-18 {
            return Ok(None);
        }
        let (ga, gb) = (dot(gap, ta), dot(gap, tb));
        t[0] -= (ga * bb - gb * ab) / determinant;
        t[1] += (gb * aa - ga * ab) / determinant;
    }
    let centre = a.spine.point(t[0]);
    // The crossing lies past the vertex along both paths. A comparison with
    // a parameter that failed to converge is false, and so refuses.
    let past = |i: usize| (t[i] - parameters[i]) * out[i] > 0.0;
    if !same(centre, b.spine.point(t[1])) || !past(0) || !past(1) {
        return Ok(None);
    }
    let Some(on_common) = a.contact(sides[0], t[0]) else {
        return Ok(None);
    };
    let step = [a, b].map(|strip| strip.radius / strip.spine.speed() / 16.0);
    let angles = |strip: &Strip, side: usize, wall: usize, t: f64| {
        Some([strip.side_angle(side, t)?, strip.side_angle(wall, t)?])
    };
    let (Some(from_a), Some(from_b)) = (
        angles(a, sides[0], walls[0], t[0]),
        angles(b, sides[1], walls[1], t[1]),
    ) else {
        return Ok(None);
    };
    // The shared curve by the share of strip a's section it has crossed:
    // the point, both path parameters, and the share of strip b's section.
    let meet = |share: f64| -> Option<(V, f64, f64, f64)> {
        if share <= 0.0 {
            return Some((on_common, t[0], t[1], 0.0));
        }
        let angle = from_a[0] + share * (from_a[1] - from_a[0]);
        let ta = first_exit(
            |ta| Some(b.spine.distance(a.point(ta, angle)?) - b.radius),
            t[0],
            -out[0] * step[0],
        )?;
        let point = a.point(ta, angle)?;
        let (tb, across) = b.locate(point, t[1])?;
        Some((
            point,
            ta,
            tb,
            (across - from_b[0]) / (from_b[1] - from_b[0]),
        ))
    };
    let Some((_, _, _, reached)) = meet(1.0) else {
        return Ok(None);
    };
    // The curve ends where it first leaves either strip.
    if !reached.is_finite() {
        return Ok(None);
    }
    let (last, deeper) = if reached > 1.0 + 1e-9 {
        let (mut low, mut high) = (0.0, 1.0);
        for _ in 0..60 {
            let middle = 0.5 * (low + high);
            match meet(middle) {
                Some((_, _, _, share)) if share > 1.0 => high = middle,
                Some(_) => low = middle,
                None => return Ok(None),
            }
        }
        (0.5 * (low + high), Some(1))
    } else if reached < 1.0 - 1e-9 {
        (1.0, Some(0))
    } else {
        (1.0, None)
    };
    // The curve's direction, per unit of share.
    let heading = |share: f64| -> Option<V> {
        let (point, ta, _, _) = meet(share)?;
        let angle = from_a[0] + share * (from_a[1] - from_a[0]);
        let gradient = unit(sub(point, b.spine.nearest(point)?))?;
        Some(mul(a.sliding(ta, angle, gradient)?, from_a[1] - from_a[0]))
    };
    // Where the curve starts both strips are tangent to the common support,
    // so their gradients agree there and say nothing of its direction. It
    // leaves along a direction in which the two strips bend alike: across its
    // section each bends as the ball does, and along its rail as the common
    // support does. Of the two such lines, the curve takes the one that runs
    // into both strips.
    let start = {
        let support = a.supports[sides[0]];
        let Some(normal) = support.normal(on_common) else {
            return Ok(None);
        };
        let across = if a.convex { -1.0 } else { 1.0 } / a.radius;
        let rails = [a.spine.tangent(t[0]), b.spine.tangent(t[1])];
        let sections = rails.map(|rail| cross(normal, rail));
        let difference = |v: V| {
            across * (dot(v, sections[0]).powi(2) - dot(v, sections[1]).powi(2))
                + support.curvature(rails[0]) * dot(v, rails[0]).powi(2)
                - support.curvature(rails[1]) * dot(v, rails[1]).powi(2)
        };
        let (e1, e2) = (rails[0], sections[0]);
        let (alpha, gamma) = (difference(e1), difference(e2));
        let beta = 0.5 * (difference(add(e1, e2)) - alpha - gamma);
        let root = beta * beta - alpha * gamma;
        if root.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
            return Ok(None);
        }
        let lines: Vec<V> = if gamma.abs() > 1e-12 {
            [1.0, -1.0]
                .into_iter()
                .map(|sign| add(e1, mul(e2, (-beta + sign * root.sqrt()) / gamma)))
                .collect()
        } else {
            vec![e2, sub(mul(e1, 2.0 * beta), mul(e2, alpha))]
        };
        let (Some(into_a), Some(into_b)) = (
            a.around(t[0], from_a[0])
                .and_then(|v| unit(mul(v, from_a[1] - from_a[0]))),
            b.around(t[1], from_b[0])
                .and_then(|v| unit(mul(v, from_b[1] - from_b[0]))),
        ) else {
            return Ok(None);
        };
        let Some(direction) = lines
            .into_iter()
            .flat_map(|line| [line, mul(line, -1.0)])
            .filter_map(unit)
            .find(|v| dot(*v, into_a) > 1e-9 && dot(*v, into_b) > 1e-9)
        else {
            return Ok(None);
        };
        // One unit of share crosses the whole of strip a's section.
        let length = a.radius * (from_a[1] - from_a[0]).abs() / dot(direction, into_a);
        let start = mul(direction, length);
        // The directions just after the start have to lead back to it.
        let Some(after) = heading(1.0 / 64.0) else {
            return Ok(None);
        };
        if distance(start, after).partial_cmp(&(0.25 * length)) != Some(std::cmp::Ordering::Less) {
            return Ok(None);
        }
        start
    };
    let curve = |s: f64| -> Option<(V, V)> {
        let share = s * last;
        let slope = if s <= 0.0 { start } else { heading(share)? };
        Some((meet(share)?.0, mul(slope, last)))
    };
    let Some(shared) = sampled(curve) else {
        return Ok(None);
    };
    let Some((end, end_a, end_b, _)) = meet(last) else {
        return Ok(None);
    };
    let (on_third, step) = match deeper {
        None => (end, None),
        Some(deeper) => {
            // The shallower strip runs on to the deeper strip's wall.
            let shallow = 1 - deeper;
            let strip = [a, b][shallow];
            let wall = [a, b][deeper].supports[walls[deeper]];
            let start = [end_a, end_b][shallow];
            let Some((_, start_angle)) = strip.locate(end, start) else {
                return Ok(None);
            };
            let finish = [from_a, from_b][shallow][1];
            let width = step[shallow];
            let cut = |s: f64| -> Option<(V, V)> {
                let angle = start_angle + s * (finish - start_angle);
                let t = root_near(|t| Some(wall.signed(strip.point(t, angle)?)), start, width)?;
                let point = strip.point(t, angle)?;
                let slope = strip.sliding(t, angle, wall.normal(point)?)?;
                Some((point, mul(slope, finish - start_angle)))
            };
            let Some(across) = sampled(cut) else {
                return Ok(None);
            };
            let Some((on_third, _)) = cut(1.0) else {
                return Ok(None);
            };
            (on_third, Some((deeper, end, Shared::Fitted(across))))
        }
    };
    Ok(Some(Junction {
        vertex,
        ends,
        third,
        common: on_common,
        on_third,
        shared: Shared::Fitted(shared),
        step,
    }))
}

/// How one strip end closes against the face it runs into.
struct OpenEnd {
    vertex: VertexId,
    /// The strip and which end of it is here.
    end: (usize, usize),
    /// The unselected edge of each support at this vertex, in strip side order.
    borders: [EdgeId; 2],
    /// The face both borders bound, which cuts the strip.
    face: FaceId,
    /// The cut across the strip, from its first rail to its second.
    cut: Samples,
}

fn open_end(
    model: &Model,
    fan: &Fan,
    strips: &[Strip],
    vertex: VertexId,
    end: (usize, usize),
    dec: &Decider,
) -> Result<Option<OpenEnd>, Refusal> {
    let strip = &strips[end.0];
    let incident = &fan.edges[&vertex];
    if incident.len() != 3 {
        return Ok(None);
    }
    let mut borders = [None; 2];
    let mut beyond = [None; 2];
    for edge in incident.iter().copied().filter(|edge| *edge != strip.edge) {
        let Some([(first, _), (second, _)]) = fan.faces.get(&edge).map(Vec::as_slice) else {
            return Ok(None);
        };
        for side in 0..2 {
            let other = if *first == strip.faces[side] {
                *second
            } else if *second == strip.faces[side] {
                *first
            } else {
                continue;
            };
            if borders[side].replace(edge).is_some() {
                return Ok(None);
            }
            beyond[side] = Some(other);
        }
    }
    let ([Some(first), Some(second)], [Some(face), Some(same)]) = (borders, beyond) else {
        return Ok(None);
    };
    if first == second || face != same || strip.faces.contains(&face) {
        return Ok(None);
    }
    let Some(wall) = Support::read(model, face, dec)? else {
        return Ok(None);
    };
    let at = strip.range[end.1];
    let Some((_, _, _, across)) = strip.section(at) else {
        return Ok(None);
    };
    let width = strip.radius / strip.spine.speed() / 16.0;
    let cut = |s: f64| -> Option<(V, V)> {
        let angle = s * across;
        let t = root_near(|t| Some(wall.signed(strip.point(t, angle)?)), at, width)?;
        let point = strip.point(t, angle)?;
        let slope = strip.sliding(t, angle, wall.normal(point)?)?;
        Some((point, mul(slope, across)))
    };
    // The section angle is the cut's natural parameter and usually its
    // smoothest. Close to the largest radius that fits, the cut runs almost
    // along the strip and passes through most of its length in a sliver of
    // that angle; there it is smoother by arc length. The fit takes whichever
    // needs fewer spans, and the angle when they need the same.
    let by_angle = sampled(cut);
    let by_arc = by_length(cut).and_then(sampled);
    let cut = match (by_angle, by_arc) {
        (Some(angle), Some(arc)) if arc.0.len() < angle.0.len() => arc,
        (Some(angle), _) => angle,
        (None, Some(arc)) => arc,
        (None, None) => return Ok(None),
    };
    Ok(Some(OpenEnd {
        vertex,
        end,
        borders: [first, second],
        face,
        cut,
    }))
}

struct Scratch<'a> {
    model: &'a mut Model,
    radius_sq: Real,
}

impl Scratch<'_> {
    fn vertex(&mut self, at: V) -> Result<VertexId, Refusal> {
        let point = self
            .model
            .geom_mut()
            .add_point(PointDef::Explicit(exact_point(at)?));
        Ok(self.model.add_vertex(Vertex { point }))
    }

    fn edge(&mut self, curve: Curve, sense: Sense, ends: [VertexId; 2]) -> EdgeId {
        let curve = self.model.geom_mut().add_curve(curve);
        self.model.add_edge(Edge {
            curve,
            sense,
            start: Some(ends[0]),
            end: Some(ends[1]),
        })
    }

    /// The arc of a circle that leaves `from` toward `through`, a point less
    /// than half a turn on and no further than the arc's middle.
    fn arc(
        &mut self,
        centre: V,
        normal: V,
        radius_sq: Real,
        [from, through]: [V; 2],
        ends: [VertexId; 2],
    ) -> Result<EdgeId, Refusal> {
        let axis = exact_vector(normal)?;
        let turn = dot(cross(sub(from, centre), sub(through, centre)), normal);
        // The seam faces away from the arc, so the arc never crosses it.
        let seam = perpendicular(&axis, sub(centre, through))?;
        Ok(self.edge(
            Curve::Circle(Circle {
                center: exact_point(centre)?,
                normal: real_vector(axis),
                radius_sq,
                ref_dir: seam,
            }),
            if turn > 0.0 {
                Sense::Same
            } else {
                Sense::Reversed
            },
            ends,
        ))
    }

    fn shared(
        &mut self,
        shared: &Shared,
        from: V,
        to: V,
        ends: [VertexId; 2],
        dec: &mut Decider,
    ) -> Result<EdgeId, Refusal> {
        match shared {
            Shared::Section { centre, normal } => {
                let Some(through) = unit(add(sub(from, *centre), sub(to, *centre))) else {
                    return Err(bad("a shared cross-section spans half a turn"));
                };
                let radius = norm(sub(from, *centre));
                self.arc(
                    *centre,
                    *normal,
                    self.radius_sq.clone(),
                    [from, add(*centre, mul(through, radius))],
                    ends,
                )
            }
            Shared::Fitted(samples) => {
                let mut samples = samples.clone();
                samples.0[0] = from;
                *samples.0.last_mut().expect("sampled curve") = to;
                let curve = quintic(&samples, dec)?;
                Ok(self.edge(Curve::BSpline(curve), Sense::Same, ends))
            }
        }
    }
}

/// Whether an edge whose ends moved along it still runs from its start to
/// its end the way it did, with something of it left.
///
/// A straight edge is judged by its chord. An arc is not: one of 300 degrees
/// cut back to its last 100 has a chord that points against the old one and
/// is still the same arc, shorter. `turning` is the centre of a circular
/// edge and the axis it turns about in its own direction. The ends are then
/// placed by angle from the old start: an end in the gap the arc leaves
/// counts as before the start or past the end, whichever it is nearer.
fn keeps_its_way(old: [V; 2], new: [V; 2], turning: Option<(V, V)>) -> bool {
    use std::f64::consts::TAU;
    let Some((centre, axis)) = turning else {
        return dot(sub(new[1], new[0]), sub(old[1], old[0])) > 0.0;
    };
    let from_start = |to: V| {
        let (a, b) = (sub(old[0], centre), sub(to, centre));
        let angle = dot(cross(a, b), axis).atan2(dot(a, b));
        if angle < 0.0 { angle + TAU } else { angle }
    };
    // A closed ring spans the whole turn.
    let span = match from_start(old[1]) {
        angle if angle > 1e-12 => angle,
        _ => TAU,
    };
    let placed = |to: V| {
        let angle = from_start(to);
        if angle <= span || angle < 0.5 * (span + TAU) {
            angle
        } else {
            angle - TAU
        }
    };
    let left = placed(new[1]) - placed(new[0]);
    left > 0.0 && left <= TAU
}

/// A walk along new edges from one vertex to another, in either direction.
fn walk(
    model: &Model,
    path: &[EdgeId],
    from: VertexId,
    to: VertexId,
) -> Result<Vec<(EdgeId, Sense)>, Refusal> {
    let ends = |edge: EdgeId| {
        let edge = model.edge(edge);
        [edge.start, edge.end]
    };
    let forward = path
        .first()
        .is_some_and(|edge| ends(*edge)[0] == Some(from))
        && path.last().is_some_and(|edge| ends(*edge)[1] == Some(to));
    let backward = path.first().is_some_and(|edge| ends(*edge)[0] == Some(to))
        && path.last().is_some_and(|edge| ends(*edge)[1] == Some(from));
    if forward {
        Ok(path.iter().map(|edge| (*edge, Sense::Same)).collect())
    } else if backward {
        Ok(path
            .iter()
            .rev()
            .map(|edge| (*edge, Sense::Reversed))
            .collect())
    } else {
        Err(bad("a junction path does not join the two rails it closes"))
    }
}

/// The admitted body and the origins of its faces by numerical face key.
pub(crate) type Built = (NumericalBody, BTreeMap<String, Vec<Origin>>);

/// Rounds the selected edges where every one of them rolls a native strip
/// and every strip end closes against another strip at its vertex.
///
/// `None` leaves the selection to the general construction: a support that is
/// not a plane or a cylinder square or parallel to its neighbour, an open
/// strip end, or a vertex this module has no junction for.
#[allow(clippy::too_many_arguments)]
pub(crate) fn build(
    model: &Model,
    body: BodyId,
    selected: &[EdgeId],
    radius: &Rational,
    step: &StepId,
    name: String,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<Option<Built>, Refusal> {
    let Some(fan) = Fan::read(model, body, dec)? else {
        return Ok(None);
    };
    let chosen: BTreeSet<EdgeId> = selected.iter().copied().collect();
    let mut strips = Vec::with_capacity(selected.len());
    for edge in selected {
        dec.charge_work(64)?;
        let Some(strip) = strip(model, body, &fan, *edge, radius.to_f64_nearest(), dec)? else {
            return Ok(None);
        };
        strips.push(strip);
    }
    // A strip end closes against the end of one other strip at its vertex, or
    // alone against the face it runs into.
    let mut at_vertex: BTreeMap<VertexId, Vec<(usize, usize)>> = BTreeMap::new();
    for (index, strip) in strips.iter().enumerate() {
        let edge = model.edge(strip.edge);
        for (end, vertex) in [edge.start, edge.end].into_iter().enumerate() {
            let Some(vertex) = vertex else {
                return Ok(None);
            };
            at_vertex.entry(vertex).or_default().push((index, end));
        }
    }
    let mut junctions = Vec::with_capacity(at_vertex.len());
    let mut open_ends = Vec::new();
    for (vertex, ends) in &at_vertex {
        dec.charge_work(4096)?;
        match ends[..] {
            [end] => {
                let Some(found) = open_end(model, &fan, &strips, *vertex, end, dec)? else {
                    return Ok(None);
                };
                open_ends.push(found);
            }
            [first, second] => {
                let Some(found) =
                    junction(model, &fan, &strips, *vertex, [first, second], &chosen, dec)?
                else {
                    return Ok(None);
                };
                junctions.push(found);
            }
            _ => return Ok(None),
        }
    }

    let mut scratch_model = model.clone();
    let mut scratch = Scratch {
        model: &mut scratch_model,
        radius_sq: Real::from_rational(radius * radius),
    };
    // New vertices: each rail end, each moved end of a third edge.
    let mut rail_ends: BTreeMap<(usize, usize, usize), (VertexId, V)> = BTreeMap::new();
    let mut moved: BTreeMap<(EdgeId, VertexId), VertexId> = BTreeMap::new();
    let mut closings: BTreeMap<(usize, usize), Vec<EdgeId>> = BTreeMap::new();
    let mut steps: BTreeMap<(FaceId, VertexId), Vec<EdgeId>> = BTreeMap::new();
    // Where each new end of a moved edge is.
    let mut placed: BTreeMap<VertexId, V> = BTreeMap::new();
    for found in &junctions {
        let [a, b] = found.ends.map(|(strip, _)| &strips[strip]);
        let common = a
            .faces
            .into_iter()
            .find(|face| b.faces.contains(face))
            .expect("checked junction");
        let sides = [a, b].map(|strip| 1 - other(strip, common));
        let walls = [a, b].map(|strip| other(strip, common));
        let on_common = scratch.vertex(found.common)?;
        let on_third = scratch.vertex(found.on_third)?;
        moved.insert((found.third, found.vertex), on_third);
        placed.insert(on_third, found.on_third);
        for i in 0..2 {
            let (strip, end) = found.ends[i];
            rail_ends.insert((strip, sides[i], end), (on_common, found.common));
        }
        match &found.step {
            None => {
                let edge = scratch.shared(
                    &found.shared,
                    found.common,
                    found.on_third,
                    [on_common, on_third],
                    dec,
                )?;
                for i in 0..2 {
                    let (strip, end) = found.ends[i];
                    rail_ends.insert((strip, walls[i], end), (on_third, found.on_third));
                    closings.insert((strip, end), vec![edge]);
                }
            }
            Some((deeper, corner, cut)) => {
                let shallow = 1 - deeper;
                let at_corner = scratch.vertex(*corner)?;
                let meet = scratch.shared(
                    &found.shared,
                    found.common,
                    *corner,
                    [on_common, at_corner],
                    dec,
                )?;
                let cut =
                    scratch.shared(cut, *corner, found.on_third, [at_corner, on_third], dec)?;
                let (deep_strip, deep_end) = found.ends[*deeper];
                let (shallow_strip, shallow_end) = found.ends[shallow];
                rail_ends.insert((deep_strip, walls[*deeper], deep_end), (at_corner, *corner));
                rail_ends.insert(
                    (shallow_strip, walls[shallow], shallow_end),
                    (on_third, found.on_third),
                );
                closings.insert((deep_strip, deep_end), vec![meet]);
                closings.insert((shallow_strip, shallow_end), vec![meet, cut]);
                steps.insert(
                    (strips[deep_strip].faces[walls[*deeper]], found.vertex),
                    vec![cut],
                );
            }
        }
    }
    for found in &open_ends {
        let from = found.cut.0[0];
        let to = *found.cut.0.last().expect("sampled cut");
        let ends = [scratch.vertex(from)?, scratch.vertex(to)?];
        let cut = scratch.shared(&Shared::Fitted(found.cut.clone()), from, to, ends, dec)?;
        let (strip, end) = found.end;
        for side in 0..2 {
            let at = [from, to][side];
            rail_ends.insert((strip, side, end), (ends[side], at));
            moved.insert((found.borders[side], found.vertex), ends[side]);
            placed.insert(ends[side], at);
        }
        closings.insert((strip, end), vec![cut]);
        steps.insert((found.face, found.vertex), vec![cut]);
    }
    // Rails, and the strip carriers.
    let mut rails: BTreeMap<(usize, usize), EdgeId> = BTreeMap::new();
    let mut carriers = Vec::with_capacity(strips.len());
    for (index, strip) in strips.iter().enumerate() {
        let middle = 0.5 * (strip.range[0] + strip.range[1]);
        let Some((centre, first, across, angle)) = strip.section(middle) else {
            return Ok(None);
        };
        let bisector = add(
            mul(first, (0.5 * angle).cos()),
            mul(across, (0.5 * angle).sin()),
        );
        for side in 0..2 {
            let (Some((from, from_at)), Some((to, to_at))) = (
                rail_ends.get(&(index, side, 0)).copied(),
                rail_ends.get(&(index, side, 1)).copied(),
            ) else {
                return Ok(None);
            };
            let edge = match strip.spine {
                Spine::Line { .. } => scratch.edge(
                    Curve::Line(Line {
                        origin: exact_point(from_at)?,
                        dir: real_vector(exact_vector(sub(to_at, from_at))?),
                    }),
                    Sense::Same,
                    [from, to],
                ),
                Spine::Circle {
                    centre: hub, axis, ..
                } => {
                    let Some(through) = strip.contact(side, middle) else {
                        return Ok(None);
                    };
                    let height = dot(sub(through, hub), axis);
                    let ring = add(hub, mul(axis, height));
                    // The exact writer wants a rational radius, not only its square.
                    let reach = rational(distance(through, ring))?;
                    scratch.arc(
                        ring,
                        axis,
                        Real::from_rational(&reach * &reach),
                        [from_at, through],
                        [from, to],
                    )?
                }
            };
            rails.insert((index, side), edge);
        }
        let surface = match strip.spine {
            Spine::Line { origin, dir } => {
                let axis = exact_vector(dir)?;
                Surface::Cylinder(Cylinder {
                    origin: exact_point(origin)?,
                    ref_dir: perpendicular(&axis, mul(bisector, -1.0))?,
                    axis: real_vector(axis),
                    radius_sq: scratch.radius_sq.clone(),
                })
            }
            Spine::Circle {
                centre: hub,
                axis,
                radius: major,
                ..
            } => {
                let exact_axis = exact_vector(axis)?;
                Surface::Torus(Torus {
                    center: exact_point(hub)?,
                    ref_dir: perpendicular(&exact_axis, sub(hub, centre))?,
                    axis: real_vector(exact_axis),
                    major: Real::from_rational(rational(major)?),
                    minor: Real::from_rational(radius.clone()),
                })
            }
        };
        carriers.push(scratch.model.geom_mut().add_surface(surface));
    }
    // Third edges with their ends moved.
    let mut replaced: BTreeMap<EdgeId, EdgeId> = BTreeMap::new();
    let ends_of = |model: &Model, edge: EdgeId| {
        let edge = model.edge(edge);
        [edge.start, edge.end]
    };
    for (edge, _) in moved.keys() {
        if replaced.contains_key(edge) {
            continue;
        }
        let mut definition = scratch.model.edge(*edge).clone();
        let before = [definition.start, definition.end];
        for end in [&mut definition.start, &mut definition.end] {
            if let Some(vertex) = *end
                && let Some(new) = moved.get(&(*edge, vertex))
            {
                *end = Some(*new);
            }
        }
        // A strip may shorten the edge it ends on, not use it up: an end
        // moved past the other end means the ball is too wide for that edge,
        // and no body with the edge turned round is the fillet asked for.
        let place = |vertex: Option<VertexId>| {
            vertex.and_then(|vertex| {
                placed
                    .get(&vertex)
                    .or_else(|| fan.positions.get(&vertex))
                    .copied()
            })
        };
        // A circle's edge turns about its normal, or against it when reversed.
        let turning = match scratch.model.geom().curve(definition.curve) {
            Curve::Circle(circle) => {
                let axis = vector64(&circle.normal, dec)?;
                Some((
                    point64(&circle.center, dec)?,
                    if definition.sense == Sense::Same {
                        axis
                    } else {
                        mul(axis, -1.0)
                    },
                ))
            }
            _ => None,
        };
        if let ([Some(old_start), Some(old_end)], Some(new_start), Some(new_end)) = (
            before.map(place),
            place(definition.start),
            place(definition.end),
        ) && !keeps_its_way([old_start, old_end], [new_start, new_end], turning)
        {
            return Err(Refusal::new(
                crate::fillet::RADIUS_TOO_LARGE,
                "the ball is wider than an edge the fillet ends on: its rail end passes that edge's far end",
            ));
        }
        // A line or a circle carries its edge wherever the ends go. A traced
        // intersection is only known between its old ends, so a straight one
        // becomes the line through the new ends and any other is left to the
        // general construction.
        if !matches!(
            scratch.model.geom().curve(definition.curve),
            Curve::Line(_) | Curve::Circle(_)
        ) {
            let Some([(first, _), (second, _)]) = fan.faces.get(edge).map(Vec::as_slice) else {
                return Ok(None);
            };
            let (Some(a), Some(b)) = (
                Support::read(model, *first, dec)?,
                Support::read(model, *second, dec)?,
            ) else {
                return Ok(None);
            };
            let (Some(start), Some(end)) = (definition.start, definition.end) else {
                return Ok(None);
            };
            let position = |vertex: VertexId| {
                placed
                    .get(&vertex)
                    .or_else(|| fan.positions.get(&vertex))
                    .copied()
            };
            let (true, Some(from), Some(to)) = (straight(&[a, b]), position(start), position(end))
            else {
                return Ok(None);
            };
            definition.curve = scratch.model.geom_mut().add_curve(Curve::Line(Line {
                origin: exact_point(from)?,
                dir: real_vector(exact_vector(sub(to, from))?),
            }));
            definition.sense = Sense::Same;
        }
        replaced.insert(*edge, scratch.model.add_edge(definition));
    }
    let strip_of: BTreeMap<EdgeId, usize> = strips
        .iter()
        .enumerate()
        .map(|(index, strip)| (strip.edge, index))
        .collect();
    // Every input face again, its boundary following the rails.
    let mut faces = Vec::new();
    let mut rail_uses: BTreeMap<(usize, usize), Sense> = BTreeMap::new();
    for old in model.faces_of_body(body) {
        let mut face = model.face(old).clone();
        let mut loops = Vec::with_capacity(face.loops.len());
        for boundary in &face.loops {
            let Loop::Coedges(uses) = model.loop_(*boundary).clone() else {
                return Ok(None);
            };
            // (new edge, sense, the old vertex this use ends at)
            let mut mapped = Vec::with_capacity(uses.len());
            for coedge in uses {
                let coedge = model.coedge(coedge).clone();
                let old_ends = ends_of(model, coedge.edge);
                let leaves = old_ends[usize::from(coedge.sense == Sense::Same)];
                let edge = if let Some(index) = strip_of.get(&coedge.edge) {
                    let side = usize::from(strips[*index].faces[1] == old);
                    rail_uses.insert((*index, side), coedge.sense);
                    rails[&(*index, side)]
                } else {
                    replaced.get(&coedge.edge).copied().unwrap_or(coedge.edge)
                };
                mapped.push((edge, coedge.sense, leaves));
            }
            let mut coedges = Vec::with_capacity(mapped.len() + 2);
            for (index, (edge, sense, leaves)) in mapped.iter().enumerate() {
                coedges.push(scratch.model.add_coedge(Coedge {
                    edge: *edge,
                    sense: *sense,
                }));
                let (next, next_sense, _) = mapped[(index + 1) % mapped.len()];
                let here = ends_of(scratch.model, *edge)[usize::from(*sense == Sense::Same)];
                let there = ends_of(scratch.model, next)[usize::from(next_sense != Sense::Same)];
                if here == there {
                    continue;
                }
                let (Some(here), Some(there), Some(leaves)) = (here, there, *leaves) else {
                    return Err(bad("a rebuilt boundary has an unbounded edge"));
                };
                let path = steps
                    .get(&(old, leaves))
                    .ok_or_else(|| bad("a rebuilt boundary has a gap no junction closes"))?;
                for (edge, sense) in walk(scratch.model, path, here, there)? {
                    coedges.push(scratch.model.add_coedge(Coedge { edge, sense }));
                }
            }
            loops.push(scratch.model.add_loop(Loop::Coedges(coedges)));
        }
        face.loops = loops;
        let origins = model.face_origins(old).to_vec();
        faces.push(scratch.model.add_named_face(face, origins));
    }
    // The strips themselves.
    for (index, strip) in strips.iter().enumerate() {
        let ends = [0, 1].map(|side| [0, 1].map(|end| rail_ends[&(index, side, end)].0));
        let mut boundary = vec![(rails[&(index, 0)], Sense::Same)];
        boundary.extend(walk(
            scratch.model,
            &closings[&(index, 1)],
            ends[0][1],
            ends[1][1],
        )?);
        boundary.push((rails[&(index, 1)], Sense::Reversed));
        boundary.extend(walk(
            scratch.model,
            &closings[&(index, 0)],
            ends[1][0],
            ends[0][0],
        )?);
        let (Some(first_use), Some(second_use)) = (
            rail_uses.get(&(index, 0)).copied(),
            rail_uses.get(&(index, 1)).copied(),
        ) else {
            return Err(bad("a rail has no use on its support"));
        };
        // Each rail runs against its support's use of it.
        if first_use == Sense::Same {
            boundary.reverse();
            for (_, sense) in &mut boundary {
                *sense = flip(*sense);
            }
        }
        let second = boundary
            .iter()
            .find(|(edge, _)| *edge == rails[&(index, 1)])
            .expect("second rail")
            .1;
        if second == second_use {
            return Err(bad("the two rails of a strip turn the same way round it"));
        }
        let coedges = boundary
            .into_iter()
            .map(|(edge, sense)| scratch.model.add_coedge(Coedge { edge, sense }))
            .collect();
        let boundary = scratch.model.add_loop(Loop::Coedges(coedges));
        faces.push(scratch.model.add_named_face(
            Face {
                surface: carriers[index],
                sense: if strip.convex {
                    Sense::Same
                } else {
                    Sense::Reversed
                },
                loops: vec![boundary],
            },
            vec![Origin {
                step: step.clone(),
                role: Role::Blend {
                    edge: strip.name.clone(),
                },
            }],
        ));
    }
    let shell = scratch.model.add_shell(Shell {
        faces: faces.clone(),
    });
    let lump = scratch.model.add_lump(Lump {
        shells: vec![shell],
    });
    let mut rounded = scratch.model.body(body).clone();
    rounded.lumps = vec![lump];
    let rounded = scratch.model.add_body(rounded);
    // The proposals carry no exact incidence facts, and claim none: admission
    // measures every incidence instead. The topology itself has to be sound.
    let faults: Vec<_> = scratch_model
        .check(rounded)
        .into_iter()
        .filter(|fault| fault.code != "topology.missing-incidence")
        .collect();
    if !faults.is_empty() {
        return Err(bad(format!(
            "the native strip network is not a closed body: {faults:?}"
        )));
    }
    let (definition, trims) = crate::numerical_embedding::definition_from_exact_model(
        &scratch_model,
        rounded,
        name,
        tolerance,
        dec,
    )
    .map_err(|error| at_stage("native strip embedding", error))?;
    let admitted = catallax_topology::numerical::build_with_native_trims(
        definition,
        trims,
        NumericalOptions {
            tolerance_mm: tolerance.clone(),
        },
        dec,
    )
    .map_err(|error| at_stage("native strip admission", error))?;
    let origins = faces
        .into_iter()
        .map(|face| {
            (
                format!("{face:?}"),
                scratch_model.face_origins(face).to_vec(),
            )
        })
        .collect();
    Ok(Some((admitted, origins)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    const TOP: &str = "top";
    const ROUND: &str = "right-round";

    #[test]
    fn a_long_arc_cut_back_keeps_its_way_and_one_used_up_does_not() {
        let at = |degrees: f64| {
            let angle = degrees.to_radians();
            [10.0 * angle.cos(), 10.0 * angle.sin(), 3.0]
        };
        let turning = Some(([0.0, 0.0, 3.0], [0.0, 0.0, 1.0]));
        let old = [at(0.0), at(300.0)];
        // Cut back to its last 100 degrees: the chord turns against the old
        // one, which a chord test would call reversed.
        let cut = [at(200.0), at(300.0)];
        assert!(dot(sub(cut[1], cut[0]), sub(old[1], old[0])) < 0.0);
        assert!(keeps_its_way(old, cut, turning));
        // Shortened at both ends, and run on a little past either end.
        assert!(keeps_its_way(old, [at(40.0), at(250.0)], turning));
        assert!(keeps_its_way(old, [at(350.0), at(300.0)], turning));
        assert!(keeps_its_way(old, [at(0.0), at(320.0)], turning));
        // The start moved past the end, or the end back past the start.
        assert!(!keeps_its_way(old, [at(310.0), at(300.0)], turning));
        assert!(!keeps_its_way(old, [at(0.0), at(350.0)], turning));
        assert!(!keeps_its_way(old, [at(120.0), at(120.0)], turning));
        // The same arc run the other way turns about the opposite axis.
        let reversed = Some(([0.0, 0.0, 3.0], [0.0, 0.0, -1.0]));
        assert!(keeps_its_way(
            [at(300.0), at(0.0)],
            [at(300.0), at(200.0)],
            reversed
        ));
        assert!(!keeps_its_way(
            [at(300.0), at(0.0)],
            [at(300.0), at(310.0)],
            reversed
        ));
        // A straight edge is its chord.
        let line = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        assert!(keeps_its_way(
            line,
            [[2.0, 0.0, 0.0], [9.0, 0.0, 0.0]],
            None
        ));
        assert!(!keeps_its_way(
            line,
            [[9.5, 0.0, 0.0], [9.0, 0.0, 0.0]],
            None
        ));
    }

    #[test]
    fn a_strip_shortens_a_long_border_arc_at_both_ends() {
        // A disc of radius 10 with a flat 8 from its centre: the top's border
        // is one arc of 286 degrees from one end of the flat's top edge round
        // to the other. Rounding that edge cuts the arc back at both ends,
        // further the larger the radius: at R12 to 132 degrees.
        let recipe = |radius: &str| {
            catallax_recipe::Recipe::parse(&format!(
                r#"{{"format":"catallax-recipe/5","name":"flatted disc","steps":[
                {{"id":"profile","op":"sketch","plane":{{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
                  "loops":[[{{"kind":"line","id":"flat","to":["8","6"]}},
                            {{"kind":"arc","id":"round","centre":["0","0"],"to":["8","-6"],"turn":"ccw"}}]]}},
                {{"id":"disc","op":"extrude","profile":"profile","distance":"20"}},
                {{"id":"round","op":"fillet","body":"disc","radius":"{radius}","edges":[
                  {{"between":[{{"step":"disc","role":{{"side":{{"segment":"flat"}}}}}},{{"step":"disc","role":"end"}}]}}]}}],
                "result":"round"}}"#
            ))
            .unwrap()
        };
        for radius in ["1", "6", "12"] {
            let recipe = recipe(radius);
            let prefix = catallax_recipe::Recipe::parse(
                &serde_json::json!({
                    "format": "catallax-recipe/5", "name": "flatted disc",
                    "steps": &serde_json::to_value(&recipe).unwrap()["steps"].as_array().unwrap()[..2],
                    "result": "disc"
                })
                .to_string(),
            )
            .unwrap();
            let catallax_recipe::Operation::Fillet { edges, .. } = &recipe.steps[2].op else {
                panic!("the last step is the fillet");
            };
            let (built, work) = with_operation_budget(100_000_000, || {
                let mut dec = Decider::default();
                let source = crate::evaluate(&prefix, &mut dec).unwrap();
                let selected = crate::fillet_build::select(
                    &source.model,
                    source.result,
                    edges,
                    true,
                    &[],
                    &mut dec,
                )
                .unwrap();
                build(
                    &source.model,
                    source.result,
                    &selected,
                    &Rational::parse(radius).unwrap(),
                    &StepId::new("round").unwrap(),
                    "flatted disc".into(),
                    &Rational::parse("0.0000001").unwrap(),
                    &mut dec,
                )
                .unwrap()
            });
            assert!(!work.exhausted, "R{radius}: {work:?}");
            let (body, _) = built.unwrap_or_else(|| panic!("R{radius}: native strips apply"));
            // Top, bottom, flat, round wall and the strip.
            assert_eq!(body.definition().faces.len(), 5, "R{radius}");
        }
    }

    /// Rounds rim edges of the owner's pocket, each named by the outer face
    /// and the pocket wall it lies between.
    fn rim(radius: &str, edges: &[(&str, &str)]) -> (Option<Built>, u64) {
        let recipe = catallax_recipe::Recipe::parse(include_str!(
            "../tests/fixtures/owner-curved-pocket.json"
        ))
        .unwrap();
        let refs: Vec<String> = edges
            .iter()
            .map(|(outer, wall)| {
                format!(
                    r#"{{"between":[{{"step":"blank","role":{{"side":{{"segment":"{outer}"}}}}}},{{"step":"pocket-tool","role":{{"side":{{"segment":"{wall}"}}}}}}]}}"#
                )
            })
            .collect();
        let selection: catallax_recipe::EdgeSelection =
            serde_json::from_str(&format!("[{}]", refs.join(","))).unwrap();
        let (built, work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::default();
            let source = crate::evaluate(&recipe, &mut dec).unwrap();
            let selected = crate::fillet_build::select(
                &source.model,
                source.result,
                &selection,
                true,
                &[],
                &mut dec,
            )
            .unwrap();
            build(
                &source.model,
                source.result,
                &selected,
                &Rational::parse(radius).unwrap(),
                &StepId::new("round").unwrap(),
                "rim".into(),
                &Rational::parse("0.0000001").unwrap(),
                &mut dec,
            )
            .unwrap()
        });
        assert!(!work.exhausted, "{work:?}");
        (built, work.used)
    }

    const WHOLE_RIM: [(&str, &str); 6] = [
        (TOP, "near"),
        (TOP, "left"),
        (TOP, "far"),
        (ROUND, "near"),
        (ROUND, "right"),
        (ROUND, "far"),
    ];

    fn blends(built: &Built) -> Vec<&NumericalFace> {
        built
            .0
            .definition()
            .faces
            .iter()
            .filter(|face| {
                built.1[&face.key]
                    .iter()
                    .any(|origin| matches!(origin.role, Role::Blend { .. }))
            })
            .collect()
    }

    #[test]
    fn the_whole_pocket_rim_rolls_at_every_radius_under_its_wall_limit() {
        // The rounded side leaves 8 mm of wall beside the pocket, so the ball
        // along that wall fits only below a radius of 4.
        for radius in ["0.5", "1", "2", "3.99"] {
            let (built, work) = rim(radius, &WHOLE_RIM);
            let built = built.expect("planes and a cylinder square to them roll native strips");
            println!("pocket rim R{radius}: {work} work");
            // Thirteen input faces and one strip an edge: the corners are
            // curves the strips share, not faces.
            assert_eq!(built.0.definition().faces.len(), 19, "R{radius}");
            let strips = blends(&built);
            assert_eq!(strips.len(), 6, "R{radius}");
            let expected = Rational::parse(radius).unwrap();
            let mut tori = 0;
            for strip in strips {
                match &strip.surface {
                    Surface::Cylinder(cylinder) => assert_eq!(
                        cylinder.radius_sq.as_rational(),
                        Some(&expected * &expected),
                        "R{radius}"
                    ),
                    Surface::Torus(torus) => {
                        tori += 1;
                        assert_eq!(torus.minor.as_rational(), Some(expected.clone()));
                        // The ball's centre runs one radius inside the R25
                        // side, to the last digit of a binary64 proposal.
                        let major = torus.major.as_rational().expect("rational major radius");
                        let gap = (&major - &(&Rational::from_i64(25) - &expected)).abs();
                        assert!(
                            gap < Rational::parse("0.000000000001").unwrap(),
                            "R{radius}"
                        );
                    }
                    other => panic!("R{radius}: a strip is carried by {other:?}"),
                }
            }
            // Two arcs round the side, four straight edges.
            assert_eq!(tori, 2, "R{radius}");
            assert!(work < 20_000_000, "R{radius} took {work} work");
        }
    }

    #[test]
    fn an_open_strip_end_is_cut_by_the_wall_it_runs_into() {
        // The rim without its thin-wall edge: both arcs end on that wall.
        let five: Vec<_> = WHOLE_RIM
            .into_iter()
            .filter(|edge| *edge != (ROUND, "right"))
            .collect();
        let (built, _) = rim("1", &five);
        let built = built.expect("an open end against a plane is cut by it");
        assert_eq!(built.0.definition().faces.len(), 18);
        assert_eq!(blends(&built).len(), 5);
        // The thin-wall edge alone ends on the two walls beside it.
        let (built, _) = rim("1", &[(ROUND, "right")]);
        let built = built.expect("one straight strip with two open ends");
        assert_eq!(built.0.definition().faces.len(), 14);
        assert_eq!(blends(&built).len(), 1);
    }
}
