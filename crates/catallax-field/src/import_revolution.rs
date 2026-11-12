//! Imported bodies of revolution whose faces include trimmed spheres or tori.
//!
//! A STEP part is a body of revolution when every face turns about one common
//! axis (planes across it, and coaxial cylinders, cones, spheres and tori) and
//! every boundary is a whole circle about that axis, or a pole or apex on it.
//! Such a part is read back as its meridian profile and rebuilt as the field's
//! full revolution of that profile. The field then holds the part by exact
//! analytic definitions, never by the imported B-rep, and each revolved face
//! takes back the name or names of the imported face it came from.
//!
//! The reading is orientation-checked. Each face's outward normal and the
//! direction its boundary circle is travelled give the direction in which the
//! profile runs with material on its left. Every profile node must be left by
//! exactly one face and entered by exactly one, the axis intervals must be
//! entered at their top and left at their bottom, and each profile loop must
//! run counter-clockwise at even nesting depth and clockwise at odd depth. A
//! part with a reversed face or an inward-facing shell is refused, never
//! published as a filled solid.
//!
//! Anything that is not such a body (a sphere cut by an oblique plane, a
//! boundary split into arcs, faces about different axes) stays a typed refusal.

use crate::import_step::{rational, rational_sqrt, unsupported};
use crate::{Face, Location, features::Feature, profile};
use catallax_numbers::{Ball, Decider, Rational, Refusal, Sign};
use catallax_recipe::{
    Axis, ElementId, Loop as Sketch, Number, Point2, Point3, Pose, Role, Segment, StepId, Turn,
    Vector3,
};
use catallax_shapes::{Curve, RPoint3, RVector3, Surface};
use catallax_topology::{BodyId, FaceId, Loop, Model, Sense};
use std::collections::{BTreeMap, BTreeSet};

/// A point of the meridian half-plane: distance from the axis, then height.
type Pair = [Rational; 2];

fn dot(a: &Pair, b: &Pair) -> Rational {
    &a[0] * &b[0] + &a[1] * &b[1]
}
fn cross(a: &Pair, b: &Pair) -> Rational {
    &a[0] * &b[1] - &a[1] * &b[0]
}
fn sub(a: &Pair, b: &Pair) -> Pair {
    [&a[0] - &b[0], &a[1] - &b[1]]
}
/// A quarter turn counter-clockwise.
fn left(v: &Pair) -> Pair {
    [-&v[1], v[0].clone()]
}
fn int(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn number(q: &Rational) -> Number {
    Number::new(q.clone())
}

struct Ring {
    centre: RPoint3,
    normal: RVector3,
    reference: RVector3,
    radius: Rational,
    /// `1` when the boundary runs counter-clockwise about its own circle normal.
    travel: i64,
}
enum Bound {
    Ring(Box<Ring>),
    Point(RPoint3),
}
struct Piece {
    id: FaceId,
    surface: Surface,
    /// `1` when the face normal is the surface normal.
    sign: i64,
    bounds: Vec<Bound>,
}

fn rpoint(p: &catallax_shapes::Point3, what: &str) -> Result<RPoint3, Refusal> {
    p.as_rational()
        .ok_or_else(|| unsupported(format!("{what} is not rational after placement")))
}
fn rvector(v: &catallax_shapes::Vector3, what: &str) -> Result<RVector3, Refusal> {
    v.as_rational()
        .ok_or_else(|| unsupported(format!("{what} is not rational after placement")))
}
fn senses(sense: Sense) -> i64 {
    if sense == Sense::Same { 1 } else { -1 }
}

fn pieces(model: &Model, body: BodyId) -> Result<Vec<Piece>, Refusal> {
    let mut pieces = Vec::new();
    for id in model.faces_of_body(body) {
        let face = model.face(id);
        let mut bounds = Vec::new();
        for lp in &face.loops {
            bounds.push(match model.loop_(*lp) {
                Loop::Vertex(v) => Bound::Point(rpoint(
                    &model.geom().position(model.vertex(*v).point),
                    "a vertex loop",
                )?),
                Loop::Coedges(cs) => {
                    let [c] = cs.as_slice() else {
                        return Err(unsupported(
                            "a boundary is not one whole parallel circle about the axis",
                        ));
                    };
                    let coedge = model.coedge(*c);
                    let edge = model.edge(coedge.edge);
                    let Curve::Circle(circle) = model.geom().curve(edge.curve) else {
                        return Err(unsupported("a boundary is not a circle"));
                    };
                    if edge.start.is_some() || edge.end.is_some() {
                        return Err(unsupported(
                            "a boundary circle is split into arcs; only whole circles are read",
                        ));
                    }
                    let radius = rational_sqrt(&rational(&circle.radius_sq, "circle radius")?)
                        .ok_or_else(|| unsupported("a boundary circle has an irrational radius"))?;
                    Bound::Ring(Box::new(Ring {
                        centre: rpoint(&circle.center, "a circle centre")?,
                        normal: rvector(&circle.normal, "a circle axis")?,
                        reference: rvector(&circle.ref_dir, "a circle seam direction")?,
                        radius,
                        travel: senses(edge.sense) * senses(coedge.sense),
                    }))
                }
            });
        }
        pieces.push(Piece {
            id,
            surface: model.geom().surface(face.surface).clone(),
            sign: senses(face.sense),
            bounds,
        });
    }
    Ok(pieces)
}

struct Frame {
    origin: RPoint3,
    /// Unit axis direction, the local `y`.
    axis: RVector3,
    /// Unit radial reference, the local `x`.
    radial: RVector3,
}

fn unit(v: &RVector3) -> Option<RVector3> {
    let length = rational_sqrt(&v.norm_squared())?;
    Some(v.scaled(&(&Rational::one() / &length)))
}

impl Frame {
    fn height(&self, p: &RPoint3) -> Rational {
        p.minus(&self.origin).dot(&self.axis)
    }

    /// The rational quaternion `[w, x, y, z]` of the rotation carrying the local
    /// axes `x`, `y`, `z` onto the radial reference, the axis and their product.
    fn quaternion(&self) -> [Rational; 4] {
        let z = self.radial.cross(&self.axis);
        let columns = [&self.radial, &self.axis, &z];
        let m = |row: usize, column: usize| -> Rational {
            let v = columns[column];
            [&v.x, &v.y, &v.z][row].clone()
        };
        let trace = m(0, 0) + m(1, 1) + m(2, 2);
        if trace > int(-1) {
            return [
                &Rational::one() + &trace,
                m(2, 1) - m(1, 2),
                m(0, 2) - m(2, 0),
                m(1, 0) - m(0, 1),
            ];
        }
        // A half turn: `m + 1` is `2 v v^T`, so any column with a nonzero
        // diagonal entry is a multiple of the axis `v`.
        let i = (0..3)
            .max_by_key(|i| m(*i, *i))
            .expect("three diagonal entries");
        let entry = |row: usize| m(row, i) + if row == i { int(1) } else { int(0) };
        [int(0), entry(0), entry(1), entry(2)]
    }
}

fn frame(pieces: &[Piece]) -> Result<Frame, Refusal> {
    let mut directions = Vec::new();
    let mut points = Vec::new();
    let mut references = Vec::new();
    for piece in pieces {
        match &piece.surface {
            Surface::Plane(p) => directions.push(rvector(&p.normal, "a plane normal")?),
            Surface::Cylinder(c) => {
                directions.push(rvector(&c.axis, "a cylinder axis")?);
                points.push(rpoint(&c.origin, "a cylinder axis point")?);
            }
            Surface::Cone(c) => {
                directions.push(rvector(&c.axis, "a cone axis")?);
                points.push(rpoint(&c.apex, "a cone apex")?);
            }
            Surface::Sphere(s) => points.push(rpoint(&s.center, "a sphere centre")?),
            Surface::Torus(t) => {
                directions.push(rvector(&t.axis, "a torus axis")?);
                points.push(rpoint(&t.center, "a torus centre")?);
            }
            _ => {
                return Err(unsupported(
                    "a face is not planar, cylindrical, conical, spherical or toroidal",
                ));
            }
        }
        for bound in &piece.bounds {
            match bound {
                Bound::Ring(r) => {
                    directions.push(r.normal.clone());
                    points.push(r.centre.clone());
                    references.push(r.reference.clone());
                }
                Bound::Point(p) => points.push(p.clone()),
            }
        }
    }
    let axis = directions
        .iter()
        .find(|d| !d.is_zero())
        .ok_or_else(|| unsupported("the faces fix no axis direction"))?;
    let axis =
        unit(axis).ok_or_else(|| unsupported("the axis direction has no rational unit vector"))?;
    if directions.iter().any(|d| !d.cross(&axis).is_zero()) {
        return Err(unsupported(
            "the faces do not turn about one axis direction; only a body of revolution builds",
        ));
    }
    let origin = points
        .first()
        .cloned()
        .ok_or_else(|| unsupported("the faces fix no point on the axis"))?;
    if points
        .iter()
        .any(|p| !p.minus(&origin).cross(&axis).is_zero())
    {
        return Err(unsupported(
            "a sphere or torus centre, cylinder, cone apex, circle centre or pole is off the common axis",
        ));
    }
    let across = |e: [i64; 3]| axis.cross(&RVector3::new(int(e[0]), int(e[1]), int(e[2])));
    let radial = references
        .into_iter()
        .chain([across([1, 0, 0]), across([0, 1, 0]), across([0, 0, 1])])
        .filter(|v| !v.is_zero())
        .find_map(|v| unit(&v))
        .ok_or_else(|| unsupported("no rational unit vector across the axis"))?;
    Ok(Frame {
        origin,
        axis,
        radial,
    })
}

/// A face's carrier in the meridian half-plane.
enum Carrier {
    Plane {
        h: Rational,
        up: i64,
    },
    Cylinder {
        r: Rational,
    },
    Cone {
        apex: Rational,
        slope: Rational,
        dir: i64,
    },
    Sphere {
        c: Rational,
        r: Rational,
    },
    Torus {
        c: Rational,
        major: Rational,
        minor: Rational,
    },
}

impl Carrier {
    fn new(frame: &Frame, surface: &Surface) -> Result<Self, Refusal> {
        let along = |v: &catallax_shapes::Vector3, what: &str| -> Result<i64, Refusal> {
            Ok(match rvector(v, what)?.dot(&frame.axis).sign() {
                Sign::Negative => -1,
                _ => 1,
            })
        };
        Ok(match surface {
            Surface::Plane(p) => Carrier::Plane {
                h: frame.height(&rpoint(&p.origin, "a plane point")?),
                up: along(&p.normal, "a plane normal")?,
            },
            Surface::Cylinder(c) => Carrier::Cylinder {
                r: rational_sqrt(&rational(&c.radius_sq, "a cylinder radius")?)
                    .ok_or_else(|| unsupported("a cylinder radius is irrational"))?,
            },
            Surface::Cone(c) => Carrier::Cone {
                apex: frame.height(&rpoint(&c.apex, "a cone apex")?),
                slope: rational(&c.slope, "a cone slope")?,
                dir: along(&c.axis, "a cone axis")?,
            },
            Surface::Sphere(s) => Carrier::Sphere {
                c: frame.height(&rpoint(&s.center, "a sphere centre")?),
                r: rational(&s.radius, "a sphere radius")?,
            },
            Surface::Torus(t) => {
                let (major, minor) = (
                    rational(&t.major, "a torus major radius")?,
                    rational(&t.minor, "a torus minor radius")?,
                );
                if minor.sign() != Sign::Positive || major <= minor {
                    return Err(unsupported(
                        "a torus is a horn or spindle; only ring tori are read",
                    ));
                }
                Carrier::Torus {
                    c: frame.height(&rpoint(&t.center, "a torus centre")?),
                    major,
                    minor,
                }
            }
            _ => unreachable!("checked when the frame was fixed"),
        })
    }

    /// Whether the meridian point lies on the carrier, exactly.
    fn on(&self, n: &Pair) -> bool {
        let square = |x: Rational| &x * &x;
        match self {
            Carrier::Plane { h, .. } => n[1] == *h,
            Carrier::Cylinder { r } => n[0] == *r,
            Carrier::Cone { apex, slope, dir } => {
                let g = (&n[1] - apex) * int(*dir);
                g.sign() != Sign::Negative && n[0] == slope * &g
            }
            Carrier::Sphere { c, r } => {
                square(n[0].clone()) + square(&n[1] - c) == square(r.clone())
            }
            Carrier::Torus { c, major, minor } => {
                square(&n[0] - major) + square(&n[1] - c) == square(minor.clone())
            }
        }
    }

    /// The surface normal at a meridian point, pointing away from the solid the
    /// carrier bounds on its own: radially out of a cylinder, cone, sphere or
    /// tube, along the plane normal.
    fn normal(&self, n: &Pair) -> Pair {
        match self {
            Carrier::Plane { up, .. } => [int(0), int(*up)],
            Carrier::Cylinder { .. } => [int(1), int(0)],
            Carrier::Cone { slope, dir, .. } => [int(1), -(slope * &int(*dir))],
            Carrier::Sphere { c, .. } => [n[0].clone(), &n[1] - c],
            Carrier::Torus { c, major, .. } => [&n[0] - major, &n[1] - c],
        }
    }
}

struct Arc {
    centre: Pair,
    radius: Rational,
    ccw: bool,
}
/// One profile edge, directed so that material lies on its left.
struct Seg {
    piece: Option<usize>,
    from: Pair,
    to: Pair,
    arc: Option<Arc>,
}

fn segment(frame: &Frame, index: usize, piece: &Piece) -> Result<Seg, Refusal> {
    let carrier = Carrier::new(frame, &piece.surface)?;
    // Each whole circle, with the way its face departs from it. A boundary that
    // runs counter-clockwise about the axis leaves the face on its left, seen
    // from outside, so the face departs along the quarter turn of its normal.
    let mut rings: Vec<(Pair, Pair, i64)> = Vec::new();
    let mut points: Vec<Pair> = Vec::new();
    for bound in &piece.bounds {
        match bound {
            Bound::Ring(r) => {
                let node = [r.radius.clone(), frame.height(&r.centre)];
                if !carrier.on(&node) {
                    return Err(unsupported("a boundary circle is off its face's surface"));
                }
                let about_axis = match r.normal.dot(&frame.axis).sign() {
                    Sign::Negative => -1,
                    _ => 1,
                };
                let spin = r.travel * about_axis;
                let n = carrier.normal(&node);
                let n = [&n[0] * &int(piece.sign), &n[1] * &int(piece.sign)];
                let depart = [&int(spin) * &-&n[1], &int(spin) * &n[0]];
                rings.push((node, depart, spin));
            }
            Bound::Point(p) => {
                let node = [int(0), frame.height(p)];
                if !carrier.on(&node) {
                    return Err(unsupported("a vertex loop is off its face's surface"));
                }
                points.push(node);
            }
        }
    }
    let Some((a, depart, spin)) = rings.first().cloned() else {
        return Err(unsupported(
            "a face has no boundary circle; a complete sphere or torus builds only alone",
        ));
    };
    let end = match (rings.len(), points.len()) {
        (2, 0) => rings[1].0.clone(),
        (1, 1) => points[0].clone(),
        (1, 0) => match &carrier {
            Carrier::Plane { h, .. } => [int(0), h.clone()],
            Carrier::Cone { apex, .. } => [int(0), apex.clone()],
            Carrier::Sphere { c, r } => match depart[1].sign() {
                Sign::Positive => [int(0), c + r],
                Sign::Negative => [int(0), c - r],
                Sign::Zero => return Err(unsupported("a spherical boundary has no pole side")),
            },
            _ => {
                return Err(unsupported(
                    "a cylindrical or toroidal face needs a boundary circle at each end",
                ));
            }
        },
        (circles, vertices) => {
            return Err(unsupported(format!(
                "a face has {circles} boundary circles and {vertices} vertex loops; a face of revolution has one profile edge"
            )));
        }
    };
    if end == a {
        return Err(unsupported("a face starts and ends on the same parallel"));
    }
    if !carrier.on(&end) {
        return Err(unsupported("a face's end is off its surface"));
    }
    let round = match &carrier {
        Carrier::Sphere { c, r } => Some(([int(0), c.clone()], r.clone())),
        Carrier::Torus { c, major, minor } => Some(([major.clone(), c.clone()], minor.clone())),
        _ => None,
    };
    let arc = match round {
        None => {
            if dot(&depart, &sub(&end, &a)).sign() != Sign::Positive {
                return Err(unsupported(
                    "a boundary's direction contradicts its face's orientation",
                ));
            }
            if let Some((b, back, _)) = rings.get(1)
                && dot(back, &sub(&a, b)).sign() != Sign::Positive
            {
                return Err(unsupported(
                    "a boundary's direction contradicts its face's orientation",
                ));
            }
            None
        }
        Some((centre, radius)) => {
            let turn = |at: &Pair, depart: &Pair| match dot(depart, &left(&sub(at, &centre))).sign()
            {
                Sign::Positive => Some(true),
                Sign::Negative => Some(false),
                Sign::Zero => None,
            };
            let ccw = turn(&a, &depart)
                .ok_or_else(|| unsupported("a boundary leaves its face across the profile"))?;
            if matches!(carrier, Carrier::Sphere { .. })
                && depart[1].sign() != sub(&end, &a)[1].sign()
            {
                return Err(unsupported(
                    "a spherical face leaves its boundary away from its other end",
                ));
            }
            if let Some((b, back, _)) = rings.get(1)
                && turn(b, back) != Some(!ccw)
            {
                return Err(unsupported(
                    "a face's two boundary circles bound different arcs of its profile",
                ));
            }
            Some(Arc {
                centre,
                radius,
                ccw,
            })
        }
    };
    if let Some((_, _, other)) = rings.get(1)
        && spin + other != 0
    {
        return Err(unsupported(
            "a face's two boundary circles run the same way about it",
        ));
    }
    // Material lies on the left of the profile, so a face whose boundary runs
    // counter-clockwise about the axis is entered at that boundary's end.
    Ok(if spin > 0 {
        Seg {
            piece: Some(index),
            from: a,
            to: end,
            arc,
        }
    } else {
        Seg {
            piece: Some(index),
            from: end,
            to: a,
            arc: arc.map(|arc| Arc {
                ccw: !arc.ccw,
                ..arc
            }),
        }
    })
}

/// Whether `v` lies on the counter-clockwise arc from `u` to `w`, all about the
/// circle centre.
fn between(u: &Pair, v: &Pair, w: &Pair) -> bool {
    match cross(u, w).sign() {
        Sign::Positive => {
            cross(u, v).sign() != Sign::Negative && cross(v, w).sign() != Sign::Negative
        }
        Sign::Negative => {
            !(cross(w, v).sign() == Sign::Positive && cross(v, u).sign() == Sign::Positive)
        }
        Sign::Zero => cross(u, v).sign() != Sign::Negative,
    }
}

/// The tangent of a directed profile edge at one of its points.
fn tangent(seg: &Seg, at: &Pair) -> Pair {
    match &seg.arc {
        None => sub(&seg.to, &seg.from),
        Some(arc) => {
            let v = sub(at, &arc.centre);
            if arc.ccw {
                left(&v)
            } else {
                [v[1].clone(), -&v[0]]
            }
        }
    }
}

/// Whether a loop of directed profile edges runs counter-clockwise. The loop's
/// lowest point, leftmost if several, is a convex corner or lies inside an arc,
/// so the profile turns left there if and only if the loop runs counter-clockwise.
fn counter_clockwise(edges: &[&Seg]) -> Result<bool, Refusal> {
    let bottom = |arc: &Arc| [arc.centre[0].clone(), &arc.centre[1] - &arc.radius];
    let reaches_bottom = |seg: &Seg| {
        seg.arc.as_ref().is_some_and(|arc| {
            let (u, w) = (sub(&seg.from, &arc.centre), sub(&seg.to, &arc.centre));
            let v = sub(&bottom(arc), &arc.centre);
            if arc.ccw {
                between(&u, &v, &w)
            } else {
                between(&w, &v, &u)
            }
        })
    };
    // The lowest candidate by height, then by distance from the axis.
    let mut lowest: Option<Pair> = None;
    let mut consider = |point: Pair| {
        let better = lowest
            .as_ref()
            .is_none_or(|best| (&point[1], &point[0]) < (&best[1], &best[0]));
        if better {
            lowest = Some(point);
        }
    };
    for seg in edges {
        consider(seg.from.clone());
        consider(seg.to.clone());
        if let Some(arc) = &seg.arc
            && reaches_bottom(seg)
        {
            consider(bottom(arc));
        }
    }
    let lowest = lowest.ok_or_else(|| unsupported("an empty profile loop"))?;
    let leaves = edges.iter().find(|seg| seg.from == lowest);
    let enters = edges.iter().find(|seg| seg.to == lowest);
    let heading = if let (Some(leaves), Some(enters)) = (leaves, enters) {
        let (out, back) = (tangent(leaves, &lowest), tangent(enters, &lowest));
        match cross(&back, &out).sign() {
            Sign::Positive => return Ok(true),
            Sign::Negative => return Ok(false),
            // A straight or tangent continuation: the common heading decides.
            Sign::Zero => out[0].sign(),
        }
    } else {
        // Inside an arc, where the tangent is horizontal.
        let arc = edges
            .iter()
            .filter(|seg| reaches_bottom(seg))
            .find_map(|seg| seg.arc.as_ref().filter(|arc| bottom(arc) == lowest))
            .ok_or_else(|| unsupported("the profile's lowest point is on no edge"))?;
        if arc.ccw {
            Sign::Positive
        } else {
            Sign::Negative
        }
    };
    match heading {
        Sign::Positive => Ok(true),
        Sign::Negative => Ok(false),
        Sign::Zero => Err(unsupported(
            "the profile's lowest point does not fix its direction",
        )),
    }
}

/// The exact volume of a full revolution, `pi * (rational + sum weight * turn)`.
///
/// Green's theorem on the profile gives `V = pi * (loop integral of rho^2 dh)`
/// over loops that run counter-clockwise at even depth and clockwise at odd
/// depth, so holes subtract. A line contributes a rational number, and so does
/// every sphere arc. A torus arc adds its tube turn about the tube's axis,
/// weighted by `major * minor^2`, which is where `pi^2` comes in.
#[derive(Clone, Debug)]
pub(crate) struct Volume {
    rational: Rational,
    turns: Vec<TubeTurn>,
}
#[derive(Clone, Debug)]
struct TubeTurn {
    weight: Rational,
    from: Pair,
    to: Pair,
    ccw: bool,
}

impl Volume {
    fn new(segs: &[Seg]) -> Self {
        let mut rational = Rational::zero();
        let mut turns = Vec::new();
        for seg in segs {
            let (a, b) = (&seg.from, &seg.to);
            match &seg.arc {
                None => {
                    rational = rational
                        + (&b[1] - &a[1]) * (&a[0] * &a[0] + &a[0] * &b[0] + &b[0] * &b[0])
                            / int(3);
                }
                Some(arc) => {
                    let [cx, cy] = &arc.centre;
                    let square = &arc.radius * &arc.radius;
                    // The antiderivative of rho^2 dh, less its angle term, at a point.
                    let at = |p: &Pair| {
                        let (u, v) = (&p[1] - cy, &p[0] - cx);
                        let cube = &(&u * &u) * &u;
                        &(cx * cx) * &u + &(cx * &u) * &v + &square * &u - &cube / &int(3)
                    };
                    rational = rational + at(b) - at(a);
                    if !cx.is_zero() {
                        turns.push(TubeTurn {
                            weight: cx * &square,
                            from: sub(a, &arc.centre),
                            to: sub(b, &arc.centre),
                            ccw: arc.ccw,
                        });
                    }
                }
            }
        }
        Self { rational, turns }
    }

    /// The volume enclosed at `bits` bits of precision. Every pi and arctangent
    /// series runs under the caller's decider, and each tube turn is charged to
    /// the caller's operation allowance first, so an exhausted allowance or a
    /// passed deadline refuses instead of returning a volume.
    pub(crate) fn enclose(&self, bits: u32, dec: &mut Decider) -> Result<Ball, Refusal> {
        let failed = || {
            Refusal::new(
                "field.volume-definition",
                "a revolved tube turn has no enclosed angle",
            )
        };
        let pi = Ball::pi_checked(bits, dec)?;
        let half = &pi * &Ball::from_rational(&Rational::from_ratio(1, 2).expect("half"), bits);
        // The counter-clockwise angle from `a` to `b`, in `(0, 2 pi]`. The
        // rational signs of the cross and dot products place it in its quadrant;
        // only the arctangent of their ratio is computed, and checked.
        let ccw_angle = |a: &Pair, b: &Pair, dec: &mut Decider| -> Result<Ball, Refusal> {
            dec.charge_work(u64::from(bits))?;
            let (c, d) = (cross(a, b), dot(a, b));
            let theta = match (c.sign(), d.sign()) {
                (Sign::Zero, Sign::Negative) => return Ok(pi.clone()),
                (Sign::Zero, _) => return Err(failed()),
                (_, Sign::Zero) => {
                    if c.sign() == Sign::Positive {
                        half.clone()
                    } else {
                        -half.clone()
                    }
                }
                (sign, d_sign) => {
                    let ratio = c.checked_div(&d).ok_or_else(failed)?;
                    let atan = Ball::from_rational(&ratio, bits).atan_checked(dec)?;
                    match (d_sign, sign) {
                        (Sign::Positive, _) => atan,
                        (_, Sign::Positive) => &atan + &pi,
                        _ => &atan - &pi,
                    }
                }
            };
            dec.check_deadline()?;
            Ok(if c.sign() == Sign::Negative {
                &theta + &(&pi * &Ball::from_rational(&int(2), bits))
            } else {
                theta
            })
        };
        let mut total = Ball::from_rational(&self.rational, bits);
        for turn in &self.turns {
            let angle = if turn.ccw {
                ccw_angle(&turn.from, &turn.to, dec)?
            } else {
                &Ball::from_rational(&int(-1), bits) * &ccw_angle(&turn.to, &turn.from, dec)?
            };
            total = &total + &(&Ball::from_rational(&turn.weight, bits) * &angle);
        }
        dec.check_deadline()?;
        Ok(&total * &pi)
    }
}

/// The field revolution of an imported body of revolution, with every revolved
/// face named for the imported face it came from.
pub(crate) fn feature(
    model: &Model,
    body: BodyId,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(Feature, Vec<Face>, Vec<Role>), Refusal> {
    let pieces = pieces(model, body)?;
    let frame = frame(&pieces)?;
    let mut segs = Vec::new();
    for (index, piece) in pieces.iter().enumerate() {
        segs.push(segment(&frame, index, piece)?);
    }

    // Every circle is left by one face and entered by one. A pole or apex on the
    // axis is touched by one face; the axis between two such points carries the
    // profile's closing edge, which makes no face.
    let mut left_at: BTreeMap<Pair, usize> = BTreeMap::new();
    let mut entered_at: BTreeMap<Pair, usize> = BTreeMap::new();
    for seg in &segs {
        *left_at.entry(seg.from.clone()).or_default() += 1;
        *entered_at.entry(seg.to.clone()).or_default() += 1;
    }
    let mut on_axis: Vec<(Rational, bool)> = Vec::new();
    for node in left_at
        .keys()
        .chain(entered_at.keys())
        .collect::<BTreeSet<_>>()
    {
        let leaves = left_at.get(node).copied().unwrap_or(0);
        let enters = entered_at.get(node).copied().unwrap_or(0);
        if node[0].is_zero() {
            if leaves + enters != 1 {
                return Err(unsupported(
                    "two faces meet at one point of the axis; the part pinches there",
                ));
            }
            on_axis.push((node[1].clone(), leaves == 1));
        } else if leaves != 1 || enters != 1 {
            return Err(unsupported(
                "a boundary circle is not shared by exactly two faces, or the faces' orientations disagree",
            ));
        }
    }
    on_axis.sort();
    if !on_axis.len().is_multiple_of(2) {
        return Err(unsupported(
            "the profile meets the axis an odd number of times",
        ));
    }
    for pair in on_axis.chunks(2) {
        let [(low, leaves), (high, _)] = pair else {
            unreachable!("even count");
        };
        if !*leaves || pair[1].1 {
            return Err(unsupported(
                "the faces face into the material along the axis; the part is inside out",
            ));
        }
        segs.push(Seg {
            piece: None,
            from: [int(0), high.clone()],
            to: [int(0), low.clone()],
            arc: None,
        });
    }

    let next: BTreeMap<&Pair, usize> = segs
        .iter()
        .enumerate()
        .map(|(index, seg)| (&seg.from, index))
        .collect();
    let mut seen = vec![false; segs.len()];
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    for first in 0..segs.len() {
        if seen[first] {
            continue;
        }
        let mut cycle = Vec::new();
        let mut at = first;
        loop {
            seen[at] = true;
            cycle.push(at);
            at = *next
                .get(&segs[at].to)
                .ok_or_else(|| unsupported("a profile edge is not entered by another"))?;
            if at == first {
                break;
            }
            if seen[at] {
                return Err(unsupported("the profile does not close into loops"));
            }
        }
        cycles.push(cycle);
    }
    let mut counter = Vec::new();
    for cycle in &cycles {
        let edges: Vec<&Seg> = cycle.iter().map(|index| &segs[*index]).collect();
        counter.push(counter_clockwise(&edges)?);
    }

    let id = |index: usize, seg: &Seg| {
        let text = match seg.piece {
            Some(piece) => format!("f{piece}"),
            None => format!("x{index}"),
        };
        ElementId::new(text).map_err(|e| unsupported(e.to_string()))
    };
    let mut sketch = Vec::new();
    for cycle in &cycles {
        let mut loop_ = Vec::new();
        for index in cycle {
            let seg = &segs[*index];
            let (id, to) = (
                id(*index, seg)?,
                Point2([number(&seg.to[0]), number(&seg.to[1])]),
            );
            loop_.push(match &seg.arc {
                None => Segment::Line { id, to },
                Some(arc) => Segment::Arc {
                    id,
                    centre: Point2([number(&arc.centre[0]), number(&arc.centre[1])]),
                    to,
                    turn: if arc.ccw { Turn::Ccw } else { Turn::Cw },
                },
            });
        }
        sketch.push(Sketch(loop_));
    }
    let [w, x, y, z] = frame.quaternion().map(|q| number(&q));
    let world = |p: &RPoint3| [number(&p.x), number(&p.y), number(&p.z)];
    let pose = Pose {
        rotation: [w, x, y, z],
        translation: Vector3(world(&frame.origin)),
    };
    let profile = profile::Profile::new(&pose, &sketch, dec)?;
    // Even nesting depth must run counter-clockwise, odd depth clockwise.
    for (index, edges) in profile.loops.iter().enumerate() {
        let mut depth = 0;
        for (other, against) in profile.loops.iter().enumerate() {
            if other != index
                && profile::loop_location(against, &profile::real(&edges[0].start), dec)?
                    == Location::Inside
            {
                depth += 1;
            }
        }
        if counter[index] != (depth % 2 == 0) {
            return Err(unsupported(
                "the faces face into the material; the part is inside out",
            ));
        }
    }
    let axis = Axis {
        point: Point3(world(&frame.origin)),
        direction: Vector3([
            number(&frame.axis.x),
            number(&frame.axis.y),
            number(&frame.axis.z),
        ]),
    };
    let (feature, faces) = Feature::revolve(&profile, step, &axis, &int(360), dec)?
        .ok_or_else(|| unsupported("the profile of revolution is empty"))?;
    let feature = feature.with_volume(Volume::new(&segs));

    let mut named = Vec::new();
    let mut roles = Vec::new();
    let mut kept = BTreeSet::new();
    for face in faces {
        let piece = match &face.origin.role {
            Role::Side { segment, .. } => segment
                .as_str()
                .strip_prefix('f')
                .and_then(|n| n.parse::<usize>().ok())
                .filter(|n| *n < pieces.len()),
            _ => None,
        }
        .ok_or_else(|| unsupported("a revolved face has no imported face"))?;
        kept.insert(piece);
        for origin in model.face_origins(pieces[piece].id) {
            let mut renamed = face.clone();
            renamed.origin = origin.clone();
            named.push(renamed);
            // The revolution still samples and trims each face by its own role.
            roles.push(face.origin.role.clone());
        }
    }
    if kept.len() != pieces.len() {
        return Err(unsupported(format!(
            "the revolution keeps {} of {} imported face names; no partial import is published",
            kept.len(),
            pieces.len()
        )));
    }
    Ok((feature, named, roles))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    fn q(n: i64, d: i64) -> Rational {
        Rational::from_ratio(n, d).expect("ratio")
    }

    // The torus cap's turn: cross 24/25 and dot -7/25 take the arctangent
    // branch, once counter-clockwise and once clockwise.
    fn turns() -> Volume {
        let turn = |ccw| TubeTurn {
            weight: int(2),
            from: [q(4, 5), q(3, 5)],
            to: [q(-4, 5), q(3, 5)],
            ccw,
        };
        Volume {
            rational: q(1, 3),
            turns: vec![turn(true), turn(false)],
        }
    }

    #[test]
    fn an_allowance_that_covers_pi_but_not_the_tube_turn_angles_refuses() {
        const BITS: u32 = 96;
        let volume = turns();
        let (pi, pi_work) =
            with_operation_budget(u64::MAX / 2, || Ball::pi_checked(BITS, &Decider::default()));
        pi.unwrap();
        let (full, work) = with_operation_budget(u64::MAX / 2, || {
            volume.enclose(BITS, &mut Decider::default())
        });
        let full = full.unwrap();
        // The angles are charged on top of pi: with an unchecked arctangent
        // the whole volume cost only its pi.
        assert!(work.used > pi_work.used + 1, "{pi_work:?} {work:?}");
        let gap = work.used - pi_work.used;
        for limit in [
            pi_work.used + 1,
            pi_work.used + gap / 4,
            pi_work.used + gap / 2,
            work.used - 1,
        ] {
            // Pi alone fits this allowance; the angles do not.
            let (pi, _) =
                with_operation_budget(limit, || Ball::pi_checked(BITS, &Decider::default()));
            assert!(pi.is_ok(), "limit {limit}");
            let (answer, spent) =
                with_operation_budget(limit, || volume.enclose(BITS, &mut Decider::default()));
            assert!(spent.exhausted, "limit {limit}: {spent:?}");
            assert_eq!(
                answer.unwrap_err().code,
                "numbers.operation-budget",
                "limit {limit}"
            );
        }
        let (again, _) =
            with_operation_budget(work.used, || volume.enclose(BITS, &mut Decider::default()));
        let again = again.unwrap();
        assert_eq!((again.lo(), again.hi()), (full.lo(), full.hi()));
    }
}
