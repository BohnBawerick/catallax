//! Reader for SolidDraft's `soliddraft-recipe/1` wire format, per RECIPE.md §6.2.
//!
//! Input types stay private so the two repositories need not share a checkout or
//! introduce a dependency cycle. Numbers enter through their shortest binary64
//! round-trip decimal, then all geometric arithmetic uses exact rationals.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use catallax_numbers::Rational;
use catallax_recipe::{
    Axis, BooleanKind, EdgeSelection, ElementId, FaceRef, Loop, Number, Operation, OriginPattern,
    Point2, Point3, Pose, Recipe, Segment, Step, StepId, Turn, Vector3,
};
use serde::Deserialize;

/// SolidDraft's profile join tolerance in millimetres.
pub const JOIN_TOLERANCE: f64 = 1e-7;

/// A converted recipe and the tag correspondence required for application history.
#[derive(Clone, Debug)]
pub struct Converted {
    /// Exact recipe. Expected-failure annotations are not part of this definition.
    pub recipe: Recipe,
    /// Original tag to normalised, collision-free step id.
    pub steps: BTreeMap<String, StepId>,
}

/// A recipe error versus an unfinished conversion or operation, never conflated
/// with impossible geometry.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum ConversionError {
    /// Malformed JSON, structure, parameters or references.
    Invalid {
        /// Original application step tag, when known.
        step: Option<String>,
        /// Explanation of the invalid input.
        detail: String,
    },
    /// A valid application operation needs kernel or converter work.
    NotYet {
        /// Original application step tag.
        step: String,
        /// Stable operation or conversion code.
        code: &'static str,
        /// What is not implemented. This is not an impossible-input claim.
        detail: String,
    },
}
impl fmt::Display for ConversionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid { step, detail } => write!(f, "soliddraft.invalid at {step:?}: {detail}"),
            Self::NotYet { step, code, detail } => write!(f, "{code} at {step}: {detail}"),
        }
    }
}
impl std::error::Error for ConversionError {}

type Result<T, E = ConversionError> = std::result::Result<T, E>;
fn invalid(step: Option<&str>, detail: impl ToString) -> ConversionError {
    ConversionError::Invalid {
        step: step.map(str::to_owned),
        detail: detail.to_string(),
    }
}
fn not_yet(step: &str, code: &'static str, detail: &str) -> ConversionError {
    ConversionError::NotYet {
        step: step.into(),
        code,
        detail: detail.into(),
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Input {
    format: String,
    name: String,
    #[serde(default)]
    description: String,
    steps: Vec<serde_json::Value>,
    result: String,
}

#[derive(Deserialize)]
#[serde(tag = "op", rename_all = "snake_case", deny_unknown_fields)]
enum Command {
    Profile {
        plane: PlaneSpec,
        rectangle: Option<Rectangle>,
        circle: Option<Circle>,
        loops: Option<Loops>,
    },
    Extrude {
        profile: String,
        distance: f64,
        direction: Option<[f64; 3]>,
    },
    Revolve {
        profile: String,
        axis_origin: [f64; 3],
        axis_direction: [f64; 3],
        angle_deg: f64,
    },
    Boolean {
        kind: BooleanKind,
        target: String,
        tool: String,
    },
    Fillet {
        body: String,
        edges: Selector,
        radius: f64,
    },
    Chamfer {
        body: String,
        edges: Selector,
        distance: f64,
    },
    Shell {
        body: String,
        open_faces: Selector,
        thickness: f64,
    },
}
#[derive(Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
enum Selector {
    All,
    CreatedBy(String),
    Nearest([f64; 3]),
    Within { min: [f64; 3], max: [f64; 3] },
    AllOf(Vec<Selector>),
}
impl Selector {
    fn check(&self, known: &BTreeMap<String, StepId>, step: &str) -> Result<()> {
        match self {
            Self::All => Ok(()),
            Self::CreatedBy(tag) => known
                .get(tag)
                .map(|_| ())
                .ok_or_else(|| invalid(Some(step), format!("unknown selector tag {tag:?}"))),
            Self::Nearest(p) => vector(*p).map(|_| ()),
            Self::Within { min, max } => {
                vector(*min)?;
                vector(*max)?;
                if (0..3).any(|i| min[i] > max[i]) {
                    return Err(invalid(Some(step), "selector minimum exceeds maximum"));
                }
                Ok(())
            }
            Self::AllOf(selectors) => {
                for selector in selectors {
                    selector.check(known, step)?;
                }
                Ok(())
            }
        }
    }
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Rectangle {
    center: [f64; 2],
    width: f64,
    height: f64,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Circle {
    center: [f64; 2],
    radius: f64,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Loops {
    outer: Vec<InputSegment>,
    #[serde(default)]
    holes: Vec<Vec<InputSegment>>,
}
#[derive(Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
enum InputSegment {
    Line {
        from: [f64; 2],
        to: [f64; 2],
    },
    Arc {
        center: [f64; 2],
        from: [f64; 2],
        to: [f64; 2],
        #[serde(default)]
        clockwise: bool,
    },
    Circle {
        center: [f64; 2],
        radius: f64,
    },
}
#[derive(Clone, Copy, Deserialize)]
enum NamedPlane {
    XY,
    YZ,
    ZX,
}
#[derive(Deserialize)]
#[serde(untagged)]
enum PlaneSpec {
    Named(NamedPlane),
    Offset(OffsetPlane),
    Full(FullPlane),
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct OffsetPlane {
    base: NamedPlane,
    origin: [f64; 3],
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct FullPlane {
    origin: [f64; 3],
    normal: [f64; 3],
    x_dir: [f64; 3],
}

fn number(x: f64) -> Result<Number> {
    if !x.is_finite() {
        return Err(invalid(None, "non-finite number"));
    }
    Number::parse(&x.to_string()).map_err(|e| invalid(None, e))
}
fn point(p: [f64; 2]) -> Result<Point2> {
    Ok(Point2([number(p[0])?, number(p[1])?]))
}
fn vector(p: [f64; 3]) -> Result<Vector3> {
    Ok(Vector3([number(p[0])?, number(p[1])?, number(p[2])?]))
}
fn positive(x: f64, step: &str) -> Result<Number> {
    if x <= 0.0 {
        return Err(invalid(Some(step), "parameter must be positive"));
    }
    number(x)
}
fn norm(v: [f64; 3]) -> f64 {
    libm::sqrt(v.iter().map(|x| x * x).sum::<f64>())
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn named_plane(base: NamedPlane, origin: [f64; 3]) -> Result<Pose> {
    let q = match base {
        NamedPlane::XY => [1, 0, 0, 0],
        NamedPlane::YZ => [1, 1, 1, 1],
        // SolidDraft's ZX u-axis is +Z, v-axis +X, normal +Y.
        NamedPlane::ZX => [1, -1, -1, -1],
    };
    Ok(Pose {
        rotation: q.map(Number::from),
        translation: vector(origin)?,
    })
}
fn plane(spec: PlaneSpec, step: &str) -> Result<Pose> {
    let p = match spec {
        PlaneSpec::Named(base) => return named_plane(base, [0.0; 3]),
        PlaneSpec::Offset(p) => return named_plane(p.base, p.origin),
        PlaneSpec::Full(p) => p,
    };
    vector(p.origin)?;
    vector(p.normal)?;
    vector(p.x_dir)?;
    let (n, x) = (p.normal, p.x_dir);
    let n_length = norm(n);
    if (n_length - 1.0).abs() >= 1e-9 || (norm(x) - 1.0).abs() >= 1e-9 || dot(n, x).abs() >= 1e-9 {
        return Err(invalid(Some(step), "plane frame is not orthonormal"));
    }
    let z = n.map(|v| v / n_length);
    let projection = dot(x, z);
    let u = std::array::from_fn(|i| x[i] - projection * z[i]);
    let u_length = norm(u);
    let u = u.map(|v| v / u_length);
    let v = cross(z, u);
    let m = [[u[0], v[0], z[0]], [u[1], v[1], z[1]], [u[2], v[2], z[2]]];
    // Pick the largest diagonal quaternion component, avoiding cancellation at pi.
    let candidates = [
        1.0 + m[0][0] + m[1][1] + m[2][2],
        1.0 + m[0][0] - m[1][1] - m[2][2],
        1.0 - m[0][0] + m[1][1] - m[2][2],
        1.0 - m[0][0] - m[1][1] + m[2][2],
    ];
    let k = (0..4)
        .max_by(|a, b| candidates[*a].total_cmp(&candidates[*b]))
        .unwrap();
    let s = 2.0 * libm::sqrt(candidates[k]);
    let q = match k {
        0 => [
            s / 4.0,
            (m[2][1] - m[1][2]) / s,
            (m[0][2] - m[2][0]) / s,
            (m[1][0] - m[0][1]) / s,
        ],
        1 => [
            (m[2][1] - m[1][2]) / s,
            s / 4.0,
            (m[0][1] + m[1][0]) / s,
            (m[0][2] + m[2][0]) / s,
        ],
        2 => [
            (m[0][2] - m[2][0]) / s,
            (m[0][1] + m[1][0]) / s,
            s / 4.0,
            (m[1][2] + m[2][1]) / s,
        ],
        _ => [
            (m[1][0] - m[0][1]) / s,
            (m[0][2] + m[2][0]) / s,
            (m[1][2] + m[2][1]) / s,
            s / 4.0,
        ],
    };
    let q = [number(q[0])?, number(q[1])?, number(q[2])?, number(q[3])?];
    // Rational quaternions define exact rotations. Certify the two input axes are
    // within the source frame tolerance, rather than accepting a float estimate.
    let r = q.each_ref().map(|v| v.value());
    let two = Rational::from_i64(2);
    let d = r
        .iter()
        .map(|v| *v * *v)
        .fold(Rational::zero(), |a, b| a + b);
    let [w, a, b, c] = r;
    let ux = [
        &(w * w + a * a - b * b - c * c) / &d,
        &(&two * &((a * b) + (w * c))) / &d,
        &(&two * &((a * c) - (w * b))) / &d,
    ];
    let nz = [
        &(&two * &((a * c) + (w * b))) / &d,
        &(&two * &((b * c) - (w * a))) / &d,
        &(w * w - a * a - b * b + c * c) / &d,
    ];
    let tol = Rational::parse("1e-9").unwrap();
    for (actual, expected) in [(&ux, x), (&nz, n)] {
        let error = (0..3)
            .map(|i| {
                let e = &actual[i] - number(expected[i]).unwrap().value();
                &e * &e
            })
            .fold(Rational::zero(), |a, b| a + b);
        if error > &tol * &tol {
            return Err(not_yet(
                step,
                "soliddraft.plane-conversion",
                "could not certify frame within source tolerance",
            ));
        }
    }
    Ok(Pose {
        rotation: q,
        translation: vector(p.origin)?,
    })
}

fn end(s: &InputSegment) -> Result<Point2> {
    match s {
        InputSegment::Line { to, .. } | InputSegment::Arc { to, .. } => point(*to),
        InputSegment::Circle { .. } => Err(invalid(None, "circle must be alone")),
    }
}
fn squared(a: &Point2, b: &Point2) -> Rational {
    (0..2)
        .map(|i| {
            let d = a.0[i].value() - b.0[i].value();
            &d * &d
        })
        .fold(Rational::zero(), |a, b| a + b)
}
fn close(a: &Point2, b: &Point2) -> bool {
    squared(a, b) <= Rational::parse("1e-14").unwrap()
}
fn snap_arc(start: &Point2, centre: &Point2, target: &Point2, step: &str) -> Result<Point2> {
    let s: [Rational; 2] = std::array::from_fn(|i| start.0[i].value() - centre.0[i].value());
    let t: [Rational; 2] = std::array::from_fn(|i| target.0[i].value() - centre.0[i].value());
    let r = &s[0] * &s[0] + &s[1] * &s[1];
    if r.is_zero() {
        return Err(invalid(Some(step), "zero-radius arc"));
    }
    if r == &t[0] * &t[0] + &t[1] * &t[1] {
        return Ok(target.clone());
    }
    // Stereographic projection from the rational start. The exact tangent
    // parameter gives a rational point on the same circle, including near pi.
    let cross = &s[0] * &t[1] - &s[1] * &t[0];
    let dot = &s[0] * &t[0] + &s[1] * &t[1];
    let (u, v) = if cross.is_zero() && dot < Rational::zero() {
        (-Rational::one(), Rational::zero())
    } else {
        let len = catallax_numbers::Ball::from_rational(&(&t[0] * &t[0] + &t[1] * &t[1]), 128)
            .sqrt()
            .ok_or_else(|| invalid(Some(step), "zero arc endpoint"))?;
        let radius = catallax_numbers::Ball::from_rational(&r, 128)
            .sqrt()
            .unwrap();
        let d = ((len.lo() + len.hi()) * (radius.lo() + radius.hi())) / Rational::from_i64(4) + dot;
        if d.is_zero() {
            return Err(not_yet(
                step,
                "soliddraft.arc-conversion",
                "unresolved antipodal arc endpoint",
            ));
        }
        let h = cross / d;
        let h2 = &h * &h;
        let divisor = &Rational::one() + &h2;
        (
            &(Rational::one() - h2) / &divisor,
            Rational::from_i64(2) * h / divisor,
        )
    };
    let result = Point2([
        Number::new(centre.0[0].value() + &(&u * &s[0] - &v * &s[1])),
        Number::new(centre.0[1].value() + &(&v * &s[0] + &u * &s[1])),
    ]);
    if !close(&result, target) {
        return Err(invalid(
            Some(step),
            "arc endpoint is outside source join tolerance",
        ));
    }
    Ok(result)
}

fn convert_loop(mut input: Vec<InputSegment>, li: usize, step: &str) -> Result<Loop> {
    if input.is_empty() {
        return Err(invalid(Some(step), "empty profile loop"));
    }
    if let [InputSegment::Circle { center, radius }] = input.as_slice() {
        return Ok(Loop(vec![Segment::Circle {
            id: ElementId::new(format!("l{li}s0")).unwrap(),
            centre: point(*center)?,
            radius: positive(*radius, step)?,
        }]));
    }
    // End on a line where possible, so closure never moves an already snapped arc.
    let shift = input
        .iter()
        .rposition(|s| matches!(s, InputSegment::Line { .. }))
        .map_or(0, |i| (i + 1) % input.len());
    input.rotate_left(shift);
    let initial = end(input.last().unwrap())?;
    let mut previous = initial.clone();
    let mut output = Vec::with_capacity(input.len());
    for (i, s) in input.iter().enumerate() {
        let id = ElementId::new(format!("l{li}s{}", (i + shift) % input.len())).unwrap();
        let (from, to) = match s {
            InputSegment::Line { from, to } | InputSegment::Arc { from, to, .. } => {
                (point(*from)?, point(*to)?)
            }
            InputSegment::Circle { .. } => return Err(invalid(Some(step), "circle must be alone")),
        };
        let original_previous = end(&input[(i + input.len() - 1) % input.len()])?;
        if !close(&original_previous, &from) {
            return Err(invalid(Some(step), "open profile loop"));
        }
        if close(&from, &to) {
            return Err(invalid(Some(step), "degenerate profile segment"));
        }
        let segment = match s {
            InputSegment::Line { .. } => Segment::Line { id, to },
            InputSegment::Arc {
                center, clockwise, ..
            } => {
                let centre = point(*center)?;
                let snapped = snap_arc(&previous, &centre, &to, step)?;
                if i + 1 == input.len() && snapped != initial {
                    return Err(not_yet(
                        step,
                        "soliddraft.arc-closure",
                        "all-arc loop requires a shared exact closure construction",
                    ));
                }
                Segment::Arc {
                    id,
                    centre,
                    to: snapped,
                    turn: if *clockwise { Turn::Cw } else { Turn::Ccw },
                }
            }
            InputSegment::Circle { .. } => unreachable!(),
        };
        previous = match &segment {
            Segment::Line { to, .. } | Segment::Arc { to, .. } => to.clone(),
            _ => unreachable!(),
        };
        output.push(segment);
    }
    Ok(Loop(output))
}

/// Converts an application recipe without running unsupported features or choosing
/// positional references by guesswork. Not-yet operations return [`ConversionError::NotYet`].
pub fn convert(json: &str) -> Result<Converted> {
    let input: Input = serde_json::from_str(json).map_err(|e| invalid(None, e))?;
    if input.format != "soliddraft-recipe/1" {
        return Err(invalid(None, "expected soliddraft-recipe/1"));
    }
    let mut steps = BTreeMap::new();
    let mut used = BTreeSet::new();
    let mut output = Vec::with_capacity(input.steps.len());
    for value in input.steps {
        let serde_json::Value::Object(mut object) = value else {
            return Err(invalid(None, "step must be an object"));
        };
        let tag = object
            .remove("id")
            .and_then(|v| v.as_str().map(str::to_owned))
            .ok_or_else(|| invalid(None, "step needs a string id"))?;
        if steps.contains_key(&tag) {
            return Err(invalid(Some(&tag), "duplicate step tag"));
        }
        if let Some(expect) = object.remove("expect")
            && expect != "ok"
            && expect != "error"
        {
            return Err(invalid(Some(&tag), "invalid expected outcome"));
        }
        let command: Command =
            serde_json::from_value(object.into()).map_err(|e| invalid(Some(&tag), e))?;
        let lookup = |id: &str| {
            steps
                .get(id)
                .cloned()
                .ok_or_else(|| invalid(Some(&tag), format!("unknown or forward input {id:?}")))
        };
        let op = match command {
            Command::Profile {
                plane: p,
                rectangle,
                circle,
                loops,
            } => {
                let loops = match (rectangle, circle, loops) {
                    (Some(r), None, None) => {
                        let c = point(r.center)?;
                        let w = positive(r.width, &tag)?.value().clone() / Rational::from_i64(2);
                        let h = positive(r.height, &tag)?.value().clone() / Rational::from_i64(2);
                        let x = c.0[0].value();
                        let y = c.0[1].value();
                        let points = [
                            [x - &w, y - &h],
                            [x + &w, y - &h],
                            [x + &w, y + &h],
                            [x - &w, y + &h],
                        ];
                        vec![Loop(
                            (0..4)
                                .map(|i| Segment::Line {
                                    id: ElementId::new(format!("l0s{i}")).unwrap(),
                                    to: Point2(points[(i + 1) % 4].clone().map(Number::new)),
                                })
                                .collect(),
                        )]
                    }
                    (None, Some(c), None) => vec![Loop(vec![Segment::Circle {
                        id: ElementId::new("l0s0").unwrap(),
                        centre: point(c.center)?,
                        radius: positive(c.radius, &tag)?,
                    }])],
                    (None, None, Some(l)) => std::iter::once(l.outer)
                        .chain(l.holes)
                        .enumerate()
                        .map(|(i, l)| convert_loop(l, i, &tag))
                        .collect::<Result<Vec<_>>>()?,
                    _ => return Err(invalid(Some(&tag), "profile needs exactly one shape")),
                };
                Operation::Sketch {
                    plane: plane(p, &tag)?,
                    loops,
                }
            }
            Command::Extrude {
                profile,
                distance,
                direction,
            } => Operation::Extrude {
                profile: lookup(&profile)?,
                distance: positive(distance, &tag)?,
                direction: direction.map(vector).transpose()?,
            },
            Command::Revolve {
                profile,
                axis_origin,
                axis_direction,
                angle_deg,
            } => {
                if angle_deg > 360.0 {
                    return Err(invalid(Some(&tag), "revolve angle exceeds 360 degrees"));
                }
                Operation::Revolve {
                    profile: lookup(&profile)?,
                    axis: Axis {
                        point: Point3(vector(axis_origin)?.0),
                        direction: vector(axis_direction)?,
                    },
                    angle: positive(angle_deg, &tag)?,
                }
            }
            Command::Boolean { kind, target, tool } => Operation::Boolean {
                kind,
                target: lookup(&target)?,
                tool: lookup(&tool)?,
            },
            Command::Fillet {
                body,
                edges,
                radius,
            } => {
                let body = lookup(&body)?;
                edges.check(&steps, &tag)?;
                let radius = positive(radius, &tag)?;
                Operation::Fillet {
                    body,
                    edges: edge_selection(&edges, &tag, "fillet")?,
                    radius,
                }
            }
            Command::Chamfer {
                body,
                edges,
                distance,
            } => {
                let body = lookup(&body)?;
                edges.check(&steps, &tag)?;
                let distance = positive(distance, &tag)?;
                Operation::Chamfer {
                    body,
                    edges: edge_selection(&edges, &tag, "chamfer")?,
                    distance,
                    other_distance: None,
                }
            }
            Command::Shell {
                body,
                open_faces,
                thickness,
            } => {
                let body = lookup(&body)?;
                open_faces.check(&steps, &tag)?;
                let thickness = positive(thickness, &tag)?;
                let open = face_refs(&open_faces, &tag, &body, &steps, &output)?;
                Operation::Shell {
                    body,
                    open,
                    thickness,
                }
            }
        };
        let base: String = tag
            .to_ascii_lowercase()
            .chars()
            .map(|c| {
                if c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_' || c == '-' {
                    c
                } else {
                    '_'
                }
            })
            .take(64)
            .collect();
        let base = if base.is_empty() {
            "step".to_owned()
        } else {
            base
        };
        let mut candidate = base.clone();
        let mut suffix = 1;
        while used.contains(&candidate) {
            let tail = format!("_{suffix}");
            candidate = format!("{}{}", &base[..base.len().min(64 - tail.len())], tail);
            suffix += 1;
        }
        used.insert(candidate.clone());
        let id = StepId::new(candidate).unwrap();
        steps.insert(tag, id.clone());
        output.push(Step { id, op });
    }
    let result = steps
        .get(&input.result)
        .cloned()
        .ok_or_else(|| invalid(None, "unknown result"))?;
    let mut recipe = Recipe::new(input.name, output, result);
    recipe.description = input.description;
    recipe.check().map_err(|e| invalid(None, e))?;
    Ok(Converted { recipe, steps })
}

/// Maps a legacy edge selector that has an exact structural counterpart. Selectors
/// that depend on where edges lie in a built body stay unconverted.
fn edge_selection(selector: &Selector, tag: &str, operation: &str) -> Result<EdgeSelection> {
    match selector {
        Selector::All => Ok(EdgeSelection::All),
        _ => Err(not_yet(
            tag,
            "features.unsupported-operation",
            &format!("{operation} edge selector has no exact structural mapping yet"),
        )),
    }
}

/// A face reference for every step the body was built from, itself included, that
/// makes faces of its own, in recipe order. Every face of the body carries an origin
/// from one of them, so together they choose exactly its faces. References are
/// mandatory, so this holds only while each of those steps keeps a face: a boolean or
/// a shell can consume every face of an earlier step, and then no exact list exists
/// without evaluating the body. Such a history yields `None`.
fn ancestor_face_makers(body: &StepId, output: &[Step]) -> Option<Vec<FaceRef>> {
    let mut needed = BTreeSet::from([body.clone()]);
    let mut refs = Vec::new();
    for step in output.iter().rev() {
        if !needed.contains(&step.id) {
            continue;
        }
        needed.extend(step.op.inputs().into_iter().map(|(id, _)| id.clone()));
        match step.op {
            Operation::Boolean { .. } | Operation::Shell { .. } => return None,
            Operation::Sketch { .. } => {}
            _ => refs.push(FaceRef {
                face: OriginPattern {
                    step: step.id.clone(),
                    role: None,
                },
                touching: Vec::new(),
                near: None,
            }),
        }
    }
    refs.reverse();
    Some(refs)
}

/// Maps a legacy open-face selector to face references. The faces a step generated
/// become a reference to that step's origins, and every face of the body becomes a
/// reference to each face-making step it was built from. Positional forms are not
/// converted.
fn face_refs(
    selector: &Selector,
    tag: &str,
    body: &StepId,
    steps: &BTreeMap<String, StepId>,
    output: &[Step],
) -> Result<Vec<FaceRef>> {
    match selector {
        Selector::All => match ancestor_face_makers(body, output) {
            Some(refs) if !refs.is_empty() => Ok(refs),
            _ => Err(not_yet(
                tag,
                "features.unsupported-operation",
                "shell open faces of a body built through a boolean or shell have no exact structural mapping yet",
            )),
        },
        Selector::CreatedBy(source) => {
            let id = &steps[source];
            let makes_faces = output
                .iter()
                .find(|step| step.id == *id)
                .is_some_and(|step| {
                    !matches!(
                        step.op,
                        Operation::Sketch { .. } | Operation::Boolean { .. }
                    )
                });
            if !makes_faces {
                return Err(not_yet(
                    tag,
                    "features.unsupported-operation",
                    "shell open faces from a step that makes no faces of its own are not converted",
                ));
            }
            Ok(vec![FaceRef {
                face: OriginPattern {
                    step: id.clone(),
                    role: None,
                },
                touching: Vec::new(),
                near: None,
            }])
        }
        _ => Err(not_yet(
            tag,
            "features.unsupported-operation",
            "shell open-face selector has no exact structural mapping yet",
        )),
    }
}
