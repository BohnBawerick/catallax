//! A through hole in a body with tolerant entities.
//!
//! Checkpoint 2 decisions, section 4.5: "a hole after a numerical fillet
//! works", "Tolerant entities never reach the exact intersector" and "Any
//! decision with a tolerant participant is 'distance against the sum of
//! tolerances' through a filtered predicate".
//!
//! The exact hole subtracts a tool body with the Boolean. This one does not
//! intersect anything. It decides where the tool's cylinder is against every
//! entity of the body, each at a reach of the tool's radius plus the
//! entity's tolerance, and then writes the hole into the body's definition:
//!
//! - Every vertex must be farther from the axis than its reach
//!   ([`predicate::Line::point_beyond`]).
//! - Every edge must be farther than its reach along its whole length. A
//!   straight edge is one sign ([`predicate::Line::segment_beyond`]). Any
//!   other edge is cut into pieces, each enclosed in a box by the curve's
//!   own interval evaluation, until every box is beyond the reach. A piece
//!   whose box is no larger than the edge's tolerance and still not beyond
//!   the reach is within tolerance of the tool, and the hole refuses.
//! - A planar face is crossed where the axis passes inside its loops. The
//!   pieces of its edges are already known to be clear of the axis, so the
//!   winding of their end points round the axis is the winding of the loops
//!   themselves. A crossed face must be square to the axis.
//! - Every other face must be clear of the tool, by the same boxes over its
//!   parameter rectangle.
//!
//! OpenCASCADE's Boolean reaches the same decisions as interferences of
//! vertex, edge and face pairs, each a distance against a sum of tolerances
//! (`BOPAlgo_PaveFiller`, `IntTools_Context`), and then intersects the
//! faces. Here the tool is a cylinder square to the faces it crosses, so
//! every new edge is a circle or a line that is written down, not found.
//! Any other case refuses by name.
//!
//! The crossings come in pairs along the axis, into the material and out of
//! it. Each crossing gets four vertices and four quarter circles as a new
//! inner loop of its face, and each pair four lines and four quarter
//! cylinders between its two rings. Admission then measures the body,
//! re-tightens the new circles and lines, which lie exactly on their
//! surfaces, and runs `Model::check`.

use catallax_numbers::{Ball, Decider, Interval, Rational, Real, Refusal, Sign};
use catallax_shapes::{Circle, Curve, Cylinder, Line, Point3, Surface, Vector3};
use catallax_topology::{
    Sense,
    numerical::{
        CopyTag, NumericalBody, NumericalEdge, NumericalFace, NumericalTrim, NumericalVertex,
        predicate::{self, Stages},
        tolerance::{self, EntityKind},
    },
};

use super::HoleTool;

/// Working precision of the enclosures, in bits. Binary64 has 53; the
/// boxes are read as binary64 intervals, so more would be rounded away.
const ENCLOSURE_BITS: u32 = 64;

/// Precision of a quarter turn and of a ring point's parameters, in bits:
/// the precision numerical admission encloses at (`BITS` in
/// `catallax-topology/src/numerical.rs`), so what is written here is read
/// back no coarser than it is checked.
const ADMISSION_BITS: u32 = 128;

/// The body with the hole, and what the decisions cost.
#[derive(Clone, Debug)]
pub struct Drilled {
    /// The admitted body.
    pub body: NumericalBody,
    /// Keys of the new wall faces with the axial interval each one spans,
    /// in millimetres from the placement point.
    pub walls: Vec<(String, [Rational; 2])>,
    /// How the signs with a tolerant participant were settled.
    pub stages: Stages,
    /// Boxes read of curved edges and faces.
    pub boxes: u64,
}

fn refuse(code: &'static str, detail: impl Into<String>) -> Refusal {
    Refusal::new(code, detail)
}

/// The tool does not cross this way, or is too near an entity to tell.
const NOT_CLEAR: &str = "hole.numerical-not-clear";

/// The tool's frame in rationals: a point, the unit axis and four unit
/// directions across it a quarter turn apart.
struct Frame {
    line: predicate::Line,
    across: [[Rational; 3]; 4],
    radius: Rational,
}

fn rational(value: &Real, what: &str) -> Result<Rational, Refusal> {
    value.as_rational().ok_or_else(|| {
        refuse(
            "hole.numerical-frame",
            format!("the hole's {what} is not rational, so its ring has no rational points"),
        )
    })
}

fn rationals(x: &Real, y: &Real, z: &Real, what: &str) -> Result<[Rational; 3], Refusal> {
    Ok([rational(x, what)?, rational(y, what)?, rational(z, what)?])
}

fn triple(x: &Real, y: &Real, z: &Real) -> Option<[Rational; 3]> {
    Some([x.as_rational()?, y.as_rational()?, z.as_rational()?])
}

// Rational arithmetic of the construction. Each operation is reserved by the
// limb rule before it runs, as the decider's is, and counts itself.
fn times(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, true)?;
    Ok(a * b)
}

fn plus(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, false)?;
    Ok(a + b)
}

fn less(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, false)?;
    Ok(a - b)
}

/// Whether `a` is at most `b`. A comparison cross-multiplies and does not
/// count itself.
fn at_most(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<bool, Refusal> {
    dec.reserve_uncounted_arithmetic(a, b, true)?;
    Ok(a <= b)
}

fn dot(a: &[Rational; 3], b: &[Rational; 3], dec: &mut Decider) -> Result<Rational, Refusal> {
    let mut sum = times(&a[0], &b[0], dec)?;
    for i in 1..3 {
        sum = plus(&sum, &times(&a[i], &b[i], dec)?, dec)?;
    }
    Ok(sum)
}

fn cross(
    a: &[Rational; 3],
    b: &[Rational; 3],
    dec: &mut Decider,
) -> Result<[Rational; 3], Refusal> {
    let mut out = [Rational::zero(), Rational::zero(), Rational::zero()];
    for (k, (i, j)) in [(1, 2), (2, 0), (0, 1)].into_iter().enumerate() {
        out[k] = less(&times(&a[i], &b[j], dec)?, &times(&a[j], &b[i], dec)?, dec)?;
    }
    Ok(out)
}

fn minus(
    a: &[Rational; 3],
    b: &[Rational; 3],
    dec: &mut Decider,
) -> Result<[Rational; 3], Refusal> {
    Ok([
        less(&a[0], &b[0], dec)?,
        less(&a[1], &b[1], dec)?,
        less(&a[2], &b[2], dec)?,
    ])
}

fn along(
    origin: &[Rational; 3],
    direction: &[Rational; 3],
    by: &Rational,
    dec: &mut Decider,
) -> Result<[Rational; 3], Refusal> {
    let mut out = origin.clone();
    for i in 0..3 {
        out[i] = plus(&origin[i], &times(&direction[i], by, dec)?, dec)?;
    }
    Ok(out)
}

fn point(p: &[Rational; 3]) -> Point3 {
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

impl Frame {
    fn of(tool: &HoleTool, dec: &mut Decider) -> Result<Self, Refusal> {
        let span = tool.cylinder_span().ok_or_else(|| {
            refuse(
                "hole.numerical-tool",
                "only a hole of one radius has a construction on a body with tolerant entities",
            )
        })?;
        let origin = rationals(&tool.point.x, &tool.point.y, &tool.point.z, "point")?;
        let axis = rationals(&tool.axis.x, &tool.axis.y, &tool.axis.z, "direction")?;
        let reference = rationals(
            &tool.reference.x,
            &tool.reference.y,
            &tool.reference.z,
            "seam direction",
        )?;
        let length = crate::extrude::sqrt(&dot(&reference, &reference, dec)?)
            .filter(|length| length.sign() == Sign::Positive)
            .ok_or_else(|| {
                refuse(
                    "hole.numerical-frame",
                    "the hole's seam direction has no rational length, so its ring has no rational points",
                )
            })?;
        dec.reserve_uncounted_arithmetic(&Rational::one(), &length, true)?;
        let scale = Rational::one() / length;
        let first = along(
            &[Rational::zero(), Rational::zero(), Rational::zero()],
            &reference,
            &scale,
            dec,
        )?;
        let second = cross(&axis, &first, dec)?;
        let opposite =
            |v: &[Rational; 3]| -> [Rational; 3] { std::array::from_fn(|i| -v[i].clone()) };
        let radius = crate::extrude::sqrt(&rational(&span.cylinder.radius_sq, "radius")?)
            .ok_or_else(|| refuse("hole.numerical-frame", "the hole's radius is not rational"))?;
        Ok(Self {
            line: predicate::Line {
                origin,
                direction: axis,
            },
            across: [
                first.clone(),
                second.clone(),
                opposite(&first),
                opposite(&second),
            ],
            radius,
        })
    }

    fn axis(&self) -> &[Rational; 3] {
        &self.line.direction
    }

    /// Coordinates across the axis, from the tool's point.
    fn across_of(&self, p: &[Rational; 3], dec: &mut Decider) -> Result<[Rational; 2], Refusal> {
        let from = minus(p, &self.line.origin, dec)?;
        Ok([
            dot(&from, &self.across[0], dec)?,
            dot(&from, &self.across[1], dec)?,
        ])
    }
}

/// What is read and decided while the tool is placed.
struct Reading<'a> {
    frame: &'a Frame,
    stages: Stages,
    boxes: u64,
}

fn ball(lo: &Rational, hi: &Rational) -> Ball {
    Ball::from_rational(lo, ENCLOSURE_BITS).hull(&Ball::from_rational(hi, ENCLOSURE_BITS))
}

/// A box round the piece of `curve` between two parameters. A spline is read
/// over the exact interval, so an end of its domain stays inside it.
fn curve_box(
    curve: &Curve,
    lo: &Rational,
    hi: &Rational,
    dec: &Decider,
) -> Result<Option<[Ball; 3]>, Refusal> {
    match curve {
        Curve::BSpline(spline) => spline.evaluate_box([lo, hi], ENCLOSURE_BITS, dec).map(Some),
        other => other.evaluate(&ball(lo, hi), dec),
    }
}

/// A box round the piece of `surface` over a parameter rectangle.
fn surface_box(
    surface: &Surface,
    u: &[Rational; 2],
    v: &[Rational; 2],
    dec: &Decider,
) -> Result<[Ball; 3], Refusal> {
    match surface {
        Surface::BSpline(spline) => {
            spline
                .carrier()
                .evaluate_box([&u[0], &u[1]], [&v[0], &v[1]], ENCLOSURE_BITS, dec)
        }
        other => other.evaluate(&ball(&u[0], &u[1]), &ball(&v[0], &v[1]), dec),
    }
}

/// The largest side of a box.
fn size(enclosure: &[Ball; 3]) -> Rational {
    enclosure
        .iter()
        .map(|side| side.rad().mul_pow2(1))
        .max()
        .expect("three sides")
}

impl Reading<'_> {
    /// Whether every point of a box is farther from the axis than `reach`:
    /// the squared distance over the box, in outward-rounded binary64
    /// intervals, less the squared reach. The axis is a unit vector.
    fn box_beyond(&self, enclosure: &[Ball; 3], reach: &Rational) -> bool {
        let line = &self.frame.line;
        let from: [Interval; 3] = std::array::from_fn(|i| {
            enclosure[i].to_interval() - Interval::from_rational(&line.origin[i])
        });
        let d: [Interval; 3] = std::array::from_fn(|i| Interval::from_rational(&line.direction[i]));
        let w = [
            from[1] * d[2] - from[2] * d[1],
            from[2] * d[0] - from[0] * d[2],
            from[0] * d[1] - from[1] * d[0],
        ];
        let distance = w[0].square() + w[1].square() + w[2].square();
        let reach = Interval::from_rational(reach).square();
        (distance - reach).lo() > 0.0
    }

    /// The parameters that cut `curve` between `from` and `to` into pieces
    /// whose boxes are all beyond `reach`, in order. A box is one reading of
    /// the curve over a parameter interval, charged before it is made.
    fn clear_pieces(
        &mut self,
        curve: &Curve,
        [from, to]: [&Rational; 2],
        reach: &Rational,
        tolerance: &Rational,
        what: &dyn Fn() -> String,
        dec: &mut Decider,
    ) -> Result<Vec<Rational>, Refusal> {
        let mut cuts = vec![from.clone()];
        let mut pending = vec![[from.clone(), to.clone()]];
        while let Some([a, b]) = pending.pop() {
            dec.charge_work(1)?;
            self.boxes += 1;
            let (lo, hi) = if at_most(&a, &b, dec)? {
                (&a, &b)
            } else {
                (&b, &a)
            };
            let enclosure = curve_box(curve, lo, hi, dec)?.ok_or_else(|| {
                refuse(
                    NOT_CLEAR,
                    format!("{} has a curve with no enclosure", what()),
                )
            })?;
            if self.box_beyond(&enclosure, reach) {
                cuts.push(b);
                continue;
            }
            if at_most(&size(&enclosure), tolerance, dec)? {
                return Err(refuse(
                    NOT_CLEAR,
                    format!(
                        "{} is within its tolerance of the hole, or crosses it",
                        what()
                    ),
                ));
            }
            let middle = plus(&a, &b, dec)?.mul_pow2(-1);
            pending.push([middle.clone(), b]);
            pending.push([a, middle]);
        }
        Ok(cuts)
    }

    /// The point of `curve` at `t`, to the working precision.
    fn at(curve: &Curve, t: &Rational, dec: &mut Decider) -> Result<[Rational; 3], Refusal> {
        dec.charge_work(1)?;
        let enclosure = curve_box(curve, t, t, dec)?
            .ok_or_else(|| refuse(NOT_CLEAR, "an edge has a curve with no enclosure"))?;
        Ok(enclosure.map(|side| side.mid()))
    }

    /// Whether a parameter rectangle of `surface` is beyond `reach`, by
    /// halving its longer side.
    fn surface_clear(
        &mut self,
        surface: &Surface,
        rectangle: [[Rational; 2]; 2],
        reach: &Rational,
        tolerance: &Rational,
        what: &dyn Fn() -> String,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        // The two sides of a piece are compared by their shares of the
        // whole rectangle, so a long side is halved before a short one
        // whatever the units of the two parameters.
        let whole = [
            less(&rectangle[0][1], &rectangle[0][0], dec)?,
            less(&rectangle[1][1], &rectangle[1][0], dec)?,
        ];
        let mut pending = vec![rectangle];
        while let Some([u, v]) = pending.pop() {
            dec.charge_work(1)?;
            self.boxes += 1;
            let enclosure = surface_box(surface, &u, &v, dec)?;
            if self.box_beyond(&enclosure, reach) {
                continue;
            }
            if at_most(&size(&enclosure), tolerance, dec)? {
                return Err(refuse(
                    NOT_CLEAR,
                    format!(
                        "{} is not clear of the hole: only a planar face square to the hole can be drilled",
                        what()
                    ),
                ));
            }
            let mut share = |side: &[Rational; 2], whole: &Rational| -> Result<Rational, Refusal> {
                let width = less(&side[1], &side[0], dec)?;
                if whole.is_zero() {
                    return Ok(Rational::zero());
                }
                dec.reserve_uncounted_arithmetic(&width, whole, true)?;
                Ok(width / whole.clone())
            };
            let shares = [share(&u, &whole[0])?, share(&v, &whole[1])?];
            let halves = |side: &[Rational; 2], dec: &mut Decider| -> Result<_, Refusal> {
                let middle = plus(&side[0], &side[1], dec)?.mul_pow2(-1);
                Ok([[side[0].clone(), middle.clone()], [middle, side[1].clone()]])
            };
            if at_most(&shares[1], &shares[0], dec)? {
                for half in halves(&u, dec)? {
                    pending.push([half, v.clone()]);
                }
            } else {
                for half in halves(&v, dec)? {
                    pending.push([u.clone(), half]);
                }
            }
        }
        Ok(())
    }
}

/// One crossing of a planar face.
struct Crossing {
    face: usize,
    /// Millimetres along the axis from the tool's point.
    depth: Rational,
    /// Whether the tool goes into the material here.
    entry: bool,
}

/// The winding number of a closed polygon round the origin, by the turn of
/// each side that crosses the positive x ray's line (the crossing-number
/// rule with signs; Sunday, "Inclusion of a point in a polygon", 2001).
fn winding(
    polygon: &[[Rational; 2]],
    stages: &mut Stages,
    dec: &mut Decider,
) -> Result<i64, Refusal> {
    let mut turns = 0;
    for (i, p) in polygon.iter().enumerate() {
        let q = &polygon[(i + 1) % polygon.len()];
        let below = |point: &[Rational; 2]| point[1].sign() != Sign::Positive;
        if below(p) && !below(q) {
            if predicate::turn(p, q, stages, dec)? == Sign::Positive {
                turns += 1;
            }
        } else if !below(p) && below(q) && predicate::turn(p, q, stages, dec)? == Sign::Negative {
            turns -= 1;
        }
    }
    Ok(turns)
}

/// Drills `tool` through `source`. `key` names the new entities.
///
/// A tool that crosses no face leaves the body as it is, as the exact hole
/// does.
pub fn through(
    source: &NumericalBody,
    key: &str,
    tool: &HoleTool,
    dec: &mut Decider,
) -> Result<Drilled, Refusal> {
    let frame = Frame::of(tool, dec)?;
    let mut reading = Reading {
        frame: &frame,
        stages: Stages::default(),
        boxes: 0,
    };
    let definition = source.definition();
    // An exact entity has no tolerance. Its pieces are cut down to rule
    // T1's fixed floor, the smallest size the rules give any entity.
    let floor = tolerance::floor_minimum_mm();
    let stored = |kind, key: &str| {
        source
            .stored_tolerance(kind, key)
            .map(|tolerance| tolerance.mm().clone())
    };

    for vertex in &definition.vertices {
        let reach = match stored(EntityKind::Vertex, &vertex.key) {
            Some(tolerance) => plus(&frame.radius, &tolerance, dec)?,
            None => frame.radius.clone(),
        };
        let beyond = frame
            .line
            .point_beyond(&vertex.position, &reach, &mut reading.stages, dec)?;
        if beyond != Sign::Positive {
            return Err(refuse(
                NOT_CLEAR,
                format!(
                    "vertex {} is within its tolerance of the hole, or inside it",
                    vertex.key
                ),
            ));
        }
    }

    // The points of every edge between its proven pieces, from its first
    // vertex to its second.
    let mut nodes: Vec<Vec<[Rational; 3]>> = Vec::with_capacity(definition.edges.len());
    for edge in &definition.edges {
        let tolerance = stored(EntityKind::Edge, &edge.key);
        let reach = match &tolerance {
            Some(tolerance) => plus(&frame.radius, tolerance, dec)?,
            None => frame.radius.clone(),
        };
        let what = || format!("edge {}", edge.key);
        let [from, to] = &edge.parameters;
        let straight = match &edge.curve {
            Curve::Line(line) => triple(&line.origin.x, &line.origin.y, &line.origin.z)
                .zip(triple(&line.dir.x, &line.dir.y, &line.dir.z)),
            _ => None,
        };
        if let Some((origin, direction)) = straight {
            let ends = [
                along(&origin, &direction, from, dec)?,
                along(&origin, &direction, to, dec)?,
            ];
            let beyond =
                frame
                    .line
                    .segment_beyond(&ends[0], &ends[1], &reach, &mut reading.stages, dec)?;
            if beyond != Sign::Positive {
                return Err(refuse(
                    NOT_CLEAR,
                    format!(
                        "{} is within its tolerance of the hole, or crosses it",
                        what()
                    ),
                ));
            }
            nodes.push(ends.to_vec());
            continue;
        }
        let cuts = reading.clear_pieces(
            &edge.curve,
            [from, to],
            &reach,
            tolerance.as_ref().unwrap_or(&floor),
            &what,
            dec,
        )?;
        let mut points = Vec::with_capacity(cuts.len());
        for cut in &cuts {
            points.push(Reading::at(&edge.curve, cut, dec)?);
        }
        nodes.push(points);
    }

    let mut crossings = Vec::new();
    for (index, face) in definition.faces.iter().enumerate() {
        let what = || format!("face {}", face.key);
        if let Surface::Plane(plane) = &face.surface {
            let mut turns = 0;
            for boundary in &face.loops {
                let mut polygon = Vec::new();
                for trim in boundary {
                    let mut points = Vec::with_capacity(nodes[trim.edge].len());
                    for node in &nodes[trim.edge] {
                        points.push(frame.across_of(node, dec)?);
                    }
                    if trim.sense == Sense::Reversed {
                        points.reverse();
                    }
                    polygon.append(&mut points);
                }
                turns += winding(&polygon, &mut reading.stages, dec)?;
            }
            if turns == 0 {
                continue;
            }
            let normal = triple(&plane.normal.x, &plane.normal.y, &plane.normal.z);
            let origin = triple(&plane.origin.x, &plane.origin.y, &plane.origin.z);
            let square = match &normal {
                Some(normal) => cross(normal, frame.axis(), dec)?
                    .iter()
                    .all(Rational::is_zero),
                None => false,
            };
            let (Some(normal), Some(origin), true) = (normal, origin, square) else {
                return Err(refuse(
                    "hole.numerical-oblique-face",
                    format!(
                        "{} is crossed at an angle: only a planar face square to the hole can be drilled",
                        what()
                    ),
                ));
            };
            let from = minus(&origin, &frame.line.origin, dec)?;
            let outward = match face.sense {
                Sense::Same => dot(&normal, frame.axis(), dec)?,
                Sense::Reversed => -dot(&normal, frame.axis(), dec)?,
            };
            crossings.push(Crossing {
                face: index,
                depth: dot(&from, frame.axis(), dec)?,
                entry: outward.sign() == Sign::Negative,
            });
            continue;
        }
        // The parameter rectangle of the face: its trims' ends, and the
        // whole of every retained u/v curve.
        let mut corners: Vec<[Rational; 2]> = Vec::new();
        for trim in face.loops.iter().flatten() {
            dec.charge_work(1)?;
            corners.extend(trim.uv.iter().cloned());
        }
        for native in source.native_trims() {
            if native.face != face.key {
                continue;
            }
            dec.charge_work(1)?;
            let Some(edge) = definition.edges.iter().find(|edge| edge.key == native.edge) else {
                continue;
            };
            let [a, b] = &edge.parameters;
            let (lo, hi) = if at_most(a, b, dec)? { (a, b) } else { (b, a) };
            let [u, v, _] = native.curve.evaluate_box([lo, hi], ENCLOSURE_BITS, dec)?;
            corners.push([u.lo(), v.lo()]);
            corners.push([u.hi(), v.hi()]);
        }
        let mut rectangle: Option<[[Rational; 2]; 2]> = None;
        for [u, v] in corners {
            rectangle = Some(match rectangle.take() {
                None => [[u.clone(), u], [v.clone(), v]],
                Some([mut us, mut vs]) => {
                    for (range, value) in [(&mut us, u), (&mut vs, v)] {
                        if !at_most(&range[0], &value, dec)? {
                            range[0] = value;
                        } else if !at_most(&value, &range[1], dec)? {
                            range[1] = value;
                        }
                    }
                    [us, vs]
                }
            });
        }
        let Some(mut rectangle) = rectangle else {
            continue;
        };
        // The box of a retained u/v curve is rounded outward. The face
        // lies in its carrier's domain, so the rectangle is cut to it.
        if let Surface::BSpline(spline) = &face.surface {
            let domain = spline.carrier().domain();
            for (range, [lo, hi]) in rectangle.iter_mut().zip(domain) {
                if !at_most(lo, &range[0], dec)? {
                    range[0] = lo.clone();
                }
                if !at_most(&range[1], hi, dec)? {
                    range[1] = hi.clone();
                }
            }
        }
        // A face stores no tolerance: its surface is the copy the others
        // are measured against.
        reading.surface_clear(&face.surface, rectangle, &frame.radius, &floor, &what, dec)?;
    }

    if crossings.is_empty() {
        return Ok(Drilled {
            body: source.clone(),
            walls: Vec::new(),
            stages: reading.stages,
            boxes: reading.boxes,
        });
    }
    dec.charge_work(crossings.len() as u64)?;
    crossings.sort_by(|a, b| a.depth.cmp(&b.depth));
    let paired = crossings.len() % 2 == 0
        && crossings
            .chunks(2)
            .all(|pair| pair[0].entry && !pair[1].entry && pair[0].depth < pair[1].depth);
    if !paired {
        return Err(refuse(
            "hole.numerical-crossings",
            "the faces the hole crosses do not pair into the material and out of it",
        ));
    }

    let mut next = definition.clone();
    let quarter = Ball::pi(ADMISSION_BITS).mid().mul_pow2(-1);
    let zero = Rational::zero();
    let radius_sq: Real = times(&frame.radius, &frame.radius, dec)?.into();
    // Ring vertices and arcs of each crossing, by quarter.
    let mut rings: Vec<([usize; 4], [usize; 4])> = Vec::with_capacity(crossings.len());
    for (ring, crossing) in crossings.iter().enumerate() {
        dec.charge_work(8)?;
        let centre = along(&frame.line.origin, frame.axis(), &crossing.depth, dec)?;
        let mut positions = Vec::with_capacity(4);
        for direction in &frame.across {
            positions.push(along(&centre, direction, &frame.radius, dec)?);
        }
        let vertices: [usize; 4] = std::array::from_fn(|k| next.vertices.len() + k);
        for (k, position) in positions.iter().enumerate() {
            next.vertices.push(NumericalVertex {
                key: format!("{key}/ring/{ring}/v{k}"),
                position: position.clone(),
            });
        }
        let arcs: [usize; 4] = std::array::from_fn(|k| next.edges.len() + k);
        for k in 0..4 {
            next.edges.push(NumericalEdge {
                key: format!("{key}/ring/{ring}/arc{k}"),
                curve: Curve::Circle(Circle {
                    center: point(&centre),
                    normal: vector(frame.axis()),
                    radius_sq: radius_sq.clone(),
                    ref_dir: vector(&frame.across[k]),
                }),
                vertices: [vertices[k], vertices[(k + 1) % 4]],
                parameters: [zero.clone(), quarter.clone()],
            });
        }
        let face = &mut next.faces[crossing.face];
        let mut uv = Vec::with_capacity(4);
        for position in &positions {
            let enclosure = position
                .each_ref()
                .map(|value| Ball::from_rational(value, ADMISSION_BITS));
            let [u, v] = face
                .surface
                .parameters_of(&enclosure, dec)?
                .ok_or_else(|| {
                    refuse(
                        "hole.numerical-assembly",
                        format!("face {} has no parameters at a ring point", face.key),
                    )
                })?;
            uv.push([u.mid(), v.mid()]);
        }
        // An inner loop runs against the outer one: with the quarters
        // where the tool goes in, and against them where it comes out.
        let order = if crossing.entry {
            [0, 1, 2, 3]
        } else {
            [3, 2, 1, 0]
        };
        let row = order
            .into_iter()
            .map(|k| {
                let ends = [uv[k].clone(), uv[(k + 1) % 4].clone()];
                if crossing.entry {
                    NumericalTrim {
                        edge: arcs[k],
                        sense: Sense::Same,
                        uv: ends,
                        tag: CopyTag::Independent,
                    }
                } else {
                    let [a, b] = ends;
                    NumericalTrim {
                        edge: arcs[k],
                        sense: Sense::Reversed,
                        uv: [b, a],
                        tag: CopyTag::Independent,
                    }
                }
            })
            .collect();
        face.loops.push(row);
        rings.push((vertices, arcs));
    }

    let mut walls = Vec::new();
    for (pair, (rings, crossings)) in rings.chunks(2).zip(crossings.chunks(2)).enumerate() {
        dec.charge_work(8)?;
        let [(low_vertices, low_arcs), (high_vertices, high_arcs)] = [&rings[0], &rings[1]];
        let [low, high] = [&crossings[0].depth, &crossings[1].depth];
        let height = less(high, low, dec)?;
        let rise = along(
            &[Rational::zero(), Rational::zero(), Rational::zero()],
            frame.axis(),
            &height,
            dec,
        )?;
        let lines: [usize; 4] = std::array::from_fn(|k| next.edges.len() + k);
        for k in 0..4 {
            let start = &next.vertices[low_vertices[k]].position;
            next.edges.push(NumericalEdge {
                key: format!("{key}/wall/{pair}/line{k}"),
                curve: Curve::Line(Line {
                    origin: point(start),
                    dir: vector(&rise),
                }),
                vertices: [low_vertices[k], high_vertices[k]],
                parameters: [zero.clone(), Rational::one()],
            });
        }
        for k in 0..4 {
            let trim = |edge, sense, uv| NumericalTrim {
                edge,
                sense,
                uv,
                tag: CopyTag::Independent,
            };
            let at = |u: &Rational, v: &Rational| [u.clone(), v.clone()];
            let face_key = format!("{key}/wall/{pair}/quarter{k}");
            next.faces.push(NumericalFace {
                key: face_key.clone(),
                surface: Surface::Cylinder(Cylinder {
                    origin: point(&frame.line.origin),
                    axis: vector(frame.axis()),
                    radius_sq: radius_sq.clone(),
                    ref_dir: vector(&frame.across[k]),
                }),
                // The material is outside the tool's cylinder.
                sense: Sense::Reversed,
                loops: vec![vec![
                    trim(lines[k], Sense::Same, [at(&zero, low), at(&zero, high)]),
                    trim(
                        high_arcs[k],
                        Sense::Same,
                        [at(&zero, high), at(&quarter, high)],
                    ),
                    trim(
                        lines[(k + 1) % 4],
                        Sense::Reversed,
                        [at(&quarter, high), at(&quarter, low)],
                    ),
                    trim(
                        low_arcs[k],
                        Sense::Reversed,
                        [at(&quarter, low), at(&zero, low)],
                    ),
                ]],
            });
            walls.push((face_key, [low.clone(), high.clone()]));
        }
    }

    let body = source.rebuild(next, source.options().clone(), dec)?;
    Ok(Drilled {
        body,
        walls,
        stages: reading.stages,
        boxes: reading.boxes,
    })
}
