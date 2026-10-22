//! Definition-based profile sweeps, independent of B-rep construction.

use crate::{
    Face, Location, geometry,
    metric::Scalar,
    profile::{self, Edge, P, Profile},
    support::Support,
};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{Axis, Origin, Role, StepId};
use catallax_shapes::{Cone, Cylinder, Plane, Pose, RPoint3, RVector3, Sphere, Surface, Torus};

#[derive(Clone, Debug)]
enum Sweep {
    Extrude {
        vector: RVector3,
        through: bool,
    },
    Revolve {
        origin: RPoint3,
        axis: RVector3,
        angle: u16,
    },
}
#[derive(Clone, Debug)]
pub(crate) struct Feature {
    loops: Vec<Vec<Edge>>,
    inverse: Pose,
    sweep: Sweep,
    axis_norm: Rational,
    merge_caps: Vec<usize>,
    supports: Vec<Support>,
    constraints: Vec<geometry::Constraint>,
    fast: FastProfile,
}
fn bad(code: &'static str) -> Refusal {
    Refusal::new(
        code,
        "feature definition requires an unsupported exact surface or frame",
    )
}
fn vector(v: &catallax_recipe::Vector3) -> RVector3 {
    let [x, y, z] = v.0.each_ref().map(|n| n.value().clone());
    RVector3::new(x, y, z)
}
fn point(p: &P) -> RPoint3 {
    RPoint3::new(p[0].clone(), p[1].clone(), Rational::zero())
}
fn z() -> RVector3 {
    RVector3::new(Rational::zero(), Rational::zero(), Rational::one())
}
fn face(surface: Surface, step: &StepId, role: Role) -> Face {
    geometry::face(
        surface,
        Origin {
            step: step.clone(),
            role,
        },
    )
}
fn side(segment: &catallax_recipe::ElementId) -> Role {
    Role::Side {
        segment: segment.clone(),
        path: None,
    }
}

impl Feature {
    pub fn extrude(
        profile: &Profile,
        step: &StepId,
        distance: Option<&Rational>,
        direction: Option<&catallax_recipe::Vector3>,
        empty_target: bool,
    ) -> Result<Option<(Self, Vec<Face>)>, Refusal> {
        if distance.is_some_and(Rational::is_zero) {
            return Err(bad("extrude.zero-distance"));
        }
        let normal = profile.pose.apply_vector(&z());
        let direction = direction.map(vector).unwrap_or(normal.clone());
        if direction.is_zero() {
            return Err(bad("extrude.zero-direction"));
        }
        if direction.dot(&normal).is_zero() {
            return Err(bad("extrude.zero-thickness"));
        }
        if profile.loops.is_empty() || empty_target {
            return Ok(None);
        }
        let through = distance.is_none();
        let v = if let Some(distance) = distance {
            let length = profile::sqrt(&direction.norm_squared())
                .ok_or_else(|| bad("features.unsupported-algebraic-placement"))?;
            direction.scaled(&(distance / &length))
        } else {
            direction
        };
        let local = profile.pose.inverse().apply_vector(&v);
        let mut faces = Vec::new();
        for edge in profile.loops.iter().flatten() {
            let surface = if let Some((c, r, _)) = &edge.circle {
                if !v.cross(&normal).is_zero() {
                    return Err(bad("features.unsupported-elliptic-cylinder"));
                }
                Surface::Cylinder(Cylinder {
                    origin: profile.world(c),
                    axis: normal.clone(),
                    radius: profile::sqrt(r)
                        .ok_or_else(|| bad("features.unsupported-algebraic-cylinder-radius"))?,
                    ref_dir: profile.world(&edge.start).minus(&profile.world(c)),
                })
            } else {
                Surface::Plane(Plane {
                    origin: profile.world(&edge.start),
                    normal: profile
                        .world(&edge.end)
                        .minus(&profile.world(&edge.start))
                        .cross(&v),
                    ref_dir: v.clone(),
                })
            };
            for id in &edge.ids {
                faces.push(face(surface.clone(), step, side(id)));
            }
        }
        if !through {
            let origin = profile.pose.apply_point(&RPoint3::origin());
            let reference = profile.pose.apply_vector(&RVector3::new(
                Rational::one(),
                Rational::zero(),
                Rational::zero(),
            ));
            for (origin, role) in [(origin.clone(), Role::Start), (origin.plus(&v), Role::End)] {
                faces.push(face(
                    Surface::Plane(Plane {
                        origin,
                        normal: normal.clone(),
                        ref_dir: reference.clone(),
                    }),
                    step,
                    role,
                ));
            }
        }
        let sweep = Sweep::Extrude {
            vector: local,
            through,
        };
        Ok(Some((
            Self::new(profile, sweep, &faces, Rational::zero(), Vec::new()),
            faces,
        )))
    }
    pub fn revolve(
        profile: &Profile,
        step: &StepId,
        axis: &Axis,
        angle: &Rational,
    ) -> Result<Option<(Self, Vec<Face>)>, Refusal> {
        if angle.sign() != Sign::Positive || *angle > Rational::from_i64(360) {
            return Err(bad("revolve.angle"));
        }
        let world_axis = vector(&axis.direction);
        if world_axis.is_zero() {
            return Err(bad("revolve.zero-axis"));
        }
        let inverse = profile.pose.inverse();
        let origin = inverse.apply_point(&RPoint3::new(
            axis.point.0[0].value().clone(),
            axis.point.0[1].value().clone(),
            axis.point.0[2].value().clone(),
        ));
        let axis = inverse.apply_vector(&world_axis);
        if !origin.z.is_zero() || !axis.z.is_zero() {
            return Err(bad("revolve.axis-plane"));
        }
        if profile.loops.is_empty() {
            return Ok(None);
        }
        let norm = axis.norm_squared();
        let project =
            |p: &RPoint3| origin.plus(&axis.scaled(&(p.minus(&origin).dot(&axis) / norm.clone())));
        let signed = |p: &RPoint3| axis.cross(&p.minus(&origin)).z;
        let mut dec = Decider::default();
        let mut depths = Vec::new();
        for (i, edges) in profile.loops.iter().enumerate() {
            let mut depth = 0;
            for (j, other) in profile.loops.iter().enumerate() {
                if i != j
                    && profile::loop_location(other, &profile::real(&edges[0].start), &mut dec)?
                        == Location::Inside
                {
                    depth += 1;
                }
            }
            depths.push(depth);
        }
        let full = *angle == Rational::from_i64(360);
        let on_axis = |e: &Edge| {
            e.circle.is_none()
                && signed(&point(&e.start)).is_zero()
                && signed(&point(&e.end)).is_zero()
        };
        for (_, edges) in profile
            .loops
            .iter()
            .enumerate()
            .filter(|(i, _)| depths[*i] % 2 == 0)
        {
            let (mut positive, mut negative) = (false, false);
            let mut pinches_axis = false;
            for (j, e) in edges.iter().enumerate() {
                let mut values = vec![
                    Real::from(signed(&point(&e.start))),
                    Real::from(signed(&point(&e.end))),
                ];
                if let Some((c, r, _)) = &e.circle {
                    let scale = profile::positive_sqrt(&(r / &norm), &mut dec)?;
                    for flip in [-1, 1] {
                        let scale = &Real::from_i64(flip) * &scale;
                        let p = [
                            Real::from(c[0].clone()) - &Real::from(axis.y.clone()) * &scale,
                            Real::from(c[1].clone()) + &Real::from(axis.x.clone()) * &scale,
                        ];
                        if e.contains(&p, &mut dec)? {
                            values.push(
                                &Real::from(axis.x.clone())
                                    * &(&p[1] - &Real::from(origin.y.clone()))
                                    - &Real::from(axis.y.clone())
                                        * &(&p[0] - &Real::from(origin.x.clone())),
                            );
                        }
                    }
                    let distance = signed(&point(c));
                    if full && &distance * &distance == r * &norm {
                        let touch = project(&point(c));
                        let p = [touch.x, touch.y];
                        pinches_axis |= (e.full || (p != e.start && p != e.end))
                            && e.contains(&profile::real(&p), &mut dec)?;
                    }
                }
                for value in values {
                    match profile::sign(&value, &mut dec)? {
                        Sign::Positive => positive = true,
                        Sign::Negative => negative = true,
                        _ => (),
                    }
                }
                pinches_axis |= full
                    && !e.full
                    && signed(&point(&e.start)).is_zero()
                    && !on_axis(e)
                    && !on_axis(&edges[(j + edges.len() - 1) % edges.len()]);
            }
            if positive == negative {
                return Err(bad("revolve.axis-interior"));
            }
            if pinches_axis {
                return Err(bad("boolean.self-touching-shell"));
            }
        }
        let angle = if full {
            360
        } else if *angle == Rational::from_i64(90) {
            90
        } else if *angle == Rational::from_i64(180) {
            180
        } else if *angle == Rational::from_i64(270) {
            270
        } else {
            return Err(bad("features.unsupported-algebraic-revolve-frame"));
        };
        let length = if angle == 90 || angle == 270 {
            Some(
                profile::sqrt(&norm)
                    .ok_or_else(|| bad("features.unsupported-algebraic-placement"))?,
            )
        } else {
            None
        };
        let mut faces = Vec::new();
        let mut merge_caps = Vec::new();
        for (index, edges) in profile.loops.iter().enumerate() {
            for edge in edges {
                let start = point(&edge.start);
                let end = point(&edge.end);
                let radial = start.minus(&project(&start));
                let surface = if let Some((c, r, _)) = &edge.circle {
                    let centre = point(c);
                    let projection = project(&centre);
                    let offset = centre.minus(&projection);
                    let radius = profile::sqrt(r).ok_or_else(|| {
                        bad(if offset.is_zero() {
                            "features.unsupported-algebraic-sphere-radius"
                        } else {
                            "features.unsupported-algebraic-torus-radius"
                        })
                    })?;
                    if offset.is_zero() {
                        Surface::Sphere(Sphere {
                            center: centre,
                            radius,
                            axis: axis.clone(),
                            ref_dir: axis.cross(&z()),
                        })
                    } else {
                        let major = profile::sqrt(&offset.norm_squared())
                            .ok_or_else(|| bad("features.unsupported-algebraic-torus-radius"))?;
                        if major <= radius {
                            return Err(bad("features.unsupported-spindle-torus"));
                        }
                        Surface::Torus(Torus {
                            center: projection,
                            axis: axis.clone(),
                            major,
                            minor: radius,
                            ref_dir: offset,
                        })
                    }
                } else {
                    let tangent = end.minus(&start);
                    if on_axis(edge) {
                        if angle == 180 && depths[index] % 2 == 0 {
                            merge_caps.push(index);
                        }
                        continue;
                    }
                    let axial = tangent.dot(&axis);
                    if axial.is_zero() {
                        Surface::Plane(Plane {
                            origin: project(&start),
                            normal: axis.clone(),
                            ref_dir: if radial.is_zero() {
                                end.minus(&project(&end))
                            } else {
                                radial
                            },
                        })
                    } else if tangent.cross(&axis).is_zero() {
                        Surface::Cylinder(Cylinder {
                            origin: project(&start),
                            axis: axis.clone(),
                            radius: profile::sqrt(&radial.norm_squared()).ok_or_else(|| {
                                bad("features.unsupported-algebraic-cylinder-radius")
                            })?,
                            ref_dir: radial,
                        })
                    } else {
                        let change = axis.cross(&tangent).z;
                        let apex = start.plus(&tangent.scaled(&(-signed(&start) / change.clone())));
                        let direction = if start.minus(&apex).dot(&axis).sign() == Sign::Positive
                            || (start == apex
                                && end.minus(&apex).dot(&axis).sign() == Sign::Positive)
                        {
                            axis.clone()
                        } else {
                            axis.scaled(&Rational::from_i64(-1))
                        };
                        Surface::Cone(Cone {
                            apex,
                            axis: direction,
                            slope: (change / axial).abs(),
                            ref_dir: if radial.is_zero() {
                                end.minus(&project(&end))
                            } else {
                                radial
                            },
                        })
                    }
                };
                let surface = surface.transformed(&profile.pose);
                for id in &edge.ids {
                    faces.push(face(surface.clone(), step, side(id)));
                }
            }
        }
        if !full {
            let normal = profile.pose.apply_vector(&z());
            let local_end_normal = match angle {
                90 => axis
                    .cross(&z())
                    .scaled(&length.as_ref().unwrap().recip().unwrap()),
                180 => z().scaled(&Rational::from_i64(-1)),
                270 => axis
                    .cross(&z())
                    .scaled(&(-length.as_ref().unwrap().recip().unwrap())),
                _ => unreachable!(),
            };
            let end_normal = profile.pose.apply_vector(&local_end_normal);
            let reference = world_axis.clone();
            let world_origin = profile.pose.apply_point(&origin);
            faces.push(face(
                Surface::Plane(Plane {
                    origin: world_origin.clone(),
                    normal,
                    ref_dir: reference.clone(),
                }),
                step,
                Role::Start,
            ));
            faces.push(face(
                Surface::Plane(Plane {
                    origin: world_origin,
                    normal: end_normal,
                    ref_dir: reference,
                }),
                step,
                Role::End,
            ));
        }
        let sweep = Sweep::Revolve {
            origin,
            axis,
            angle,
        };
        Ok(Some((
            Self::new(profile, sweep, &faces, norm, merge_caps),
            faces,
        )))
    }
    fn new(
        profile: &Profile,
        sweep: Sweep,
        faces: &[Face],
        axis_norm: Rational,
        merge_caps: Vec<usize>,
    ) -> Self {
        let mut constraints = Vec::new();
        let mut supports = Vec::new();
        for face in faces {
            if !constraints
                .iter()
                .any(|c: &geometry::Constraint| c.polynomial == face.polynomial)
            {
                constraints.push(geometry::Constraint::new(face.polynomial.clone()));
                supports.push(Support::new(&face.surface));
            }
        }
        Self {
            loops: profile.loops.clone(),
            inverse: profile.pose.inverse(),
            fast: FastProfile::new(profile, &sweep, &axis_norm),
            axis_norm,
            merge_caps,
            sweep,
            constraints,
            supports,
        }
    }
    pub fn location(&self, p: &[Rational; 3]) -> Location {
        self.classify(p, &mut Decider::default())
            .unwrap_or(Location::Indeterminate)
    }
    fn classify(&self, p: &[Rational; 3], dec: &mut Decider) -> Result<Location, Refusal> {
        let local =
            self.inverse
                .apply_point(&RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone()));
        let norm = &self.axis_norm;
        match &self.sweep {
            Sweep::Extrude { vector, through } => {
                let t = &local.z / &vector.z;
                if !through && (t.sign() == Sign::Negative || t > Rational::one()) {
                    return Ok(Location::Outside);
                }
                let q = [&local.x - &(&t * &vector.x), &local.y - &(&t * &vector.y)];
                let loc = profile::location(&self.loops, &profile::real(&q), dec)?;
                Ok(
                    if !through && loc == Location::Inside && (t.is_zero() || t == Rational::one())
                    {
                        Location::Boundary
                    } else {
                        loc
                    },
                )
            }
            Sweep::Revolve {
                origin,
                axis,
                angle,
                ..
            } => {
                let v = local.minus(origin);
                let t = v.dot(axis) / norm.clone();
                let b = axis.cross(&v).z / norm.clone();
                let k = profile::positive_sqrt(&(&b * &b + (&v.z * &v.z) / norm.clone()), dec)?;
                let mut boundary = false;
                for flip in [-1, 1] {
                    let angular = angular(*angle, flip, b.sign(), v.z.sign());
                    if angular == Location::Outside {
                        continue;
                    }
                    let k = &Real::from_i64(flip) * &k;
                    let q = [
                        Real::from(&origin.x + &(&axis.x * &t)) - &Real::from(axis.y.clone()) * &k,
                        Real::from(&origin.y + &(&axis.y * &t)) + &Real::from(axis.x.clone()) * &k,
                    ];
                    let loc = profile::location(&self.loops, &q, dec)?;
                    if loc == Location::Inside && angular == Location::Inside {
                        return Ok(Location::Inside);
                    }
                    boundary |= loc != Location::Outside;
                }
                if boundary && *angle == 360 && b.is_zero() && v.z.is_zero() {
                    let q = [&local.x, &local.y].map(|q| q.clone());
                    for edge in self.loops.iter().flatten().filter(|e| e.circle.is_none()) {
                        if axis.cross(&point(&edge.start).minus(origin)).z.is_zero()
                            && axis.cross(&point(&edge.end).minus(origin)).z.is_zero()
                            && q != edge.start
                            && q != edge.end
                            && edge.contains(&profile::real(&q), dec)?
                        {
                            return Ok(Location::Inside);
                        }
                    }
                }
                Ok(if boundary {
                    Location::Boundary
                } else {
                    Location::Outside
                })
            }
        }
    }
    pub fn on_face(&self, role: &Role, p: &[Rational; 3]) -> Option<bool> {
        self.face_contains(role, p, &mut Decider::default()).ok()
    }
    fn face_contains(
        &self,
        role: &Role,
        p: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let local =
            self.inverse
                .apply_point(&RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone()));
        let norm = &self.axis_norm;
        let merge_caps = &self.merge_caps;
        match &self.sweep {
            Sweep::Extrude { vector, through } => {
                let t = &local.z / &vector.z;
                if !through && (t.sign() == Sign::Negative || t > Rational::one()) {
                    return Ok(false);
                }
                let q =
                    profile::real(&[&local.x - &(&t * &vector.x), &local.y - &(&t * &vector.y)]);
                match role {
                    Role::Side { segment, .. } => {
                        for e in self
                            .loops
                            .iter()
                            .flatten()
                            .filter(|e| e.ids.contains(segment))
                        {
                            if e.contains(&q, dec)? {
                                return Ok(true);
                            }
                        }
                        Ok(false)
                    }
                    Role::Start | Role::End => {
                        Ok(profile::location(&self.loops, &q, dec)? != Location::Outside)
                    }
                    _ => Ok(false),
                }
            }
            Sweep::Revolve {
                origin,
                axis,
                angle,
            } => {
                let v = local.minus(origin);
                let t = v.dot(axis) / norm.clone();
                let b = axis.cross(&v).z / norm.clone();
                if let Role::Side { segment, .. } = role {
                    let k = profile::positive_sqrt(&(&b * &b + (&v.z * &v.z) / norm.clone()), dec)?;
                    for flip in [-1, 1] {
                        if angular(*angle, flip, b.sign(), v.z.sign()) == Location::Outside {
                            continue;
                        }
                        let k = &Real::from_i64(flip) * &k;
                        let q = [
                            Real::from(&origin.x + &(&axis.x * &t))
                                - &Real::from(axis.y.clone()) * &k,
                            Real::from(&origin.y + &(&axis.y * &t))
                                + &Real::from(axis.x.clone()) * &k,
                        ];
                        for e in self
                            .loops
                            .iter()
                            .flatten()
                            .filter(|e| e.ids.contains(segment))
                        {
                            if e.contains(&q, dec)? {
                                return Ok(true);
                            }
                        }
                    }
                    return Ok(false);
                }
                let (cosine, sine) = if matches!(role, Role::Start) {
                    (1, 0)
                } else {
                    match angle {
                        90 => (0, 1),
                        180 => (-1, 0),
                        270 => (0, -1),
                        _ => (1, 0),
                    }
                };
                let radial = v.minus(&axis.scaled(&t));
                let rotated = if sine == 0 {
                    radial.scaled(&Rational::from_i64(cosine))
                } else {
                    axis.cross(&radial)
                        .scaled(&(Rational::from_i64(-sine) / profile::sqrt(norm).unwrap()))
                };
                if !rotated.z.is_zero() {
                    return Ok(false);
                }
                let q = origin.plus(&axis.scaled(&t)).plus(&rotated);
                let q = [q.x, q.y];
                if profile::location(&self.loops, &profile::real(&q), dec)? != Location::Outside {
                    return Ok(true);
                }
                if *angle == 180 {
                    let reflected = origin
                        .plus(&axis.scaled(&t))
                        .plus(&rotated.scaled(&Rational::from_i64(-1)));
                    let reflected = profile::real(&[reflected.x, reflected.y]);
                    if profile::location(&self.loops, &reflected, dec)? != Location::Outside {
                        for &i in merge_caps {
                            if profile::loop_location(&self.loops[i], &reflected, dec)?
                                != Location::Outside
                            {
                                return Ok(true);
                            }
                        }
                    }
                }
                Ok(false)
            }
        }
    }
    pub fn proven_distance(&self, p: &RPoint3, location: Location) -> Ball {
        if location == Location::Boundary {
            return Ball::integer(0);
        }
        let mut supports = self.supports.iter();
        let mut residual = supports.next().expect("nonempty sweep supports").proven(p);
        for support in supports {
            residual = residual.min(&support.proven(p));
        }
        match location {
            Location::Inside => -residual,
            Location::Indeterminate => residual.hull(&(-&residual)),
            _ => residual,
        }
    }
    pub fn ray(&self, p: &[Rational; 3], d: &[Rational; 3]) -> Location {
        // Pick an exact displacement within every support polynomial's initial
        // sign cell. This is a certified germ, not a geometric tolerance.
        let mut epsilon = Rational::one();
        for constraint in &self.constraints {
            let coefficients = geometry::ray_coefficients(constraint, p, d);
            if let Some((i, c)) = coefficients.iter().enumerate().find(|(_, c)| !c.is_zero()) {
                let tail = coefficients[i + 1..]
                    .iter()
                    .fold(Rational::zero(), |a, c| a + c.abs());
                let bound = c.abs() / (c.abs() + tail + Rational::one());
                epsilon = epsilon.min(bound);
            }
        }
        let sample = std::array::from_fn(|i| &p[i] + &(&epsilon * &d[i]));
        // A ray identically on a true boundary remains a boundary, never an interior witness.
        self.location(&sample)
    }
    pub fn fast_distance(&self, p: [f64; 3]) -> f64 {
        if !self.fast.valid {
            return f64::NAN;
        }
        let inside = self.fast.inside(p);
        let mut residual = f64::INFINITY;
        for s in &self.supports {
            let distance = s.fast(&p);
            if !distance.is_finite() {
                return f64::NAN;
            }
            residual = residual.min(distance);
        }
        if inside { -residual } else { residual }
    }
}

fn angular(angle: u16, flip: i64, b: Sign, z: Sign) -> Location {
    if angle == 360 {
        return Location::Inside;
    }
    let oriented = |s| {
        if flip == 1 {
            s
        } else {
            match s {
                Sign::Positive => Sign::Negative,
                Sign::Negative => Sign::Positive,
                Sign::Zero => Sign::Zero,
            }
        }
    };
    let (b, z) = (oriented(b), oriented(z));
    let outside = match angle {
        90 => b == Sign::Negative || z == Sign::Negative,
        180 => z == Sign::Negative,
        270 => b == Sign::Positive && z == Sign::Negative,
        _ => unreachable!(),
    };
    if outside {
        Location::Outside
    } else if z == Sign::Zero && (angle == 180 || b != Sign::Negative)
        || b == Sign::Zero
            && (angle != 180
                && z == if angle == 90 {
                    Sign::Positive
                } else {
                    Sign::Negative
                })
    {
        Location::Boundary
    } else {
        Location::Inside
    }
}

#[derive(Clone, Debug)]
struct FastEdge {
    start: [f64; 2],
    end: [f64; 2],
    circle: Option<([f64; 2], f64, bool)>,
    full: bool,
}
#[derive(Clone, Debug)]
struct FastProfile {
    matrix: [[f64; 3]; 3],
    translation: [f64; 3],
    edges: Vec<FastEdge>,
    sweep: [f64; 5],
    angle: u16,
    through: bool,
    valid: bool,
}
impl FastProfile {
    fn new(profile: &Profile, sweep: &Sweep, norm: &Rational) -> Self {
        let inverse = profile.pose.inverse();
        let matrix = std::array::from_fn(|i| {
            std::array::from_fn(|j| inverse.rotation.matrix()[i][j].to_f64_nearest())
        });
        let translation = [
            &inverse.translation.x,
            &inverse.translation.y,
            &inverse.translation.z,
        ]
        .map(Rational::to_f64_nearest);
        let edges = profile
            .loops
            .iter()
            .flatten()
            .map(|e| FastEdge {
                start: e.start.each_ref().map(Rational::to_f64_nearest),
                end: e.end.each_ref().map(Rational::to_f64_nearest),
                circle: e.circle.as_ref().map(|(c, r, t)| {
                    (
                        c.each_ref().map(Rational::to_f64_nearest),
                        r.to_f64_nearest(),
                        *t == catallax_recipe::Turn::Ccw,
                    )
                }),
                full: e.full,
            })
            .collect::<Vec<_>>();
        let (parameters, angle) = match sweep {
            Sweep::Extrude { vector, .. } => (
                [
                    vector.x.to_f64_nearest(),
                    vector.y.to_f64_nearest(),
                    vector.z.to_f64_nearest(),
                    0.0,
                    0.0,
                ],
                0,
            ),
            Sweep::Revolve {
                origin,
                axis,
                angle,
                ..
            } => (
                [
                    origin.x.to_f64_nearest(),
                    origin.y.to_f64_nearest(),
                    axis.x.to_f64_nearest(),
                    axis.y.to_f64_nearest(),
                    norm.to_f64_nearest(),
                ],
                *angle,
            ),
        };
        let valid = matrix
            .iter()
            .flatten()
            .chain(translation.iter())
            .chain(parameters.iter())
            .all(|v| v.is_finite())
            && edges.iter().all(|e| {
                e.start.iter().chain(e.end.iter()).all(|v| v.is_finite())
                    && e.circle
                        .as_ref()
                        .is_none_or(|(c, r, _)| c.iter().all(|v| v.is_finite()) && r.is_finite())
            });
        Self {
            matrix,
            translation,
            edges,
            sweep: parameters,
            angle,
            through: matches!(sweep, Sweep::Extrude { through: true, .. }),
            valid,
        }
    }
    fn inside(&self, p: [f64; 3]) -> bool {
        if !self.valid {
            return false;
        }
        let p: [f64; 3] = std::array::from_fn(|i| {
            ((self.matrix[i][0] * p[0] + self.matrix[i][1] * p[1]) + self.matrix[i][2] * p[2])
                + self.translation[i]
        });
        let a = self.sweep;
        if self.angle == 0 {
            let t = p[2] / a[2];
            return (self.through || (0.0..=1.0).contains(&t))
                && self.region([p[0] - t * a[0], p[1] - t * a[1]]);
        }
        let v = [p[0] - a[0], p[1] - a[1]];
        let t = (v[0] * a[2] + v[1] * a[3]) / a[4];
        let b = (a[2] * v[1] - a[3] * v[0]) / a[4];
        let k = Scalar::sqrt(&(b * b + p[2] * p[2] / a[4]));
        for flip in [-1.0, 1.0] {
            let angular = match self.angle {
                90 => b * flip >= 0.0 && p[2] * flip >= 0.0,
                180 => p[2] * flip >= 0.0,
                270 => b * flip <= 0.0 || p[2] * flip >= 0.0,
                _ => true,
            };
            if angular
                && self.region([
                    a[0] + a[2] * t - a[3] * k * flip,
                    a[1] + a[3] * t + a[2] * k * flip,
                ])
            {
                return true;
            }
        }
        false
    }
    fn region(&self, p: [f64; 2]) -> bool {
        let cross = |a: [f64; 2], b: [f64; 2]| a[0] * b[1] - a[1] * b[0];
        let sub = |a: [f64; 2], b: [f64; 2]| [a[0] - b[0], a[1] - b[1]];
        let mut crossings = 0;
        for e in &self.edges {
            if let Some((c, r, ccw)) = e.circle {
                let dy = p[1] - c[1];
                let discriminant = r - dy * dy;
                if discriminant <= 0.0 {
                    continue;
                }
                let root = Scalar::sqrt(&discriminant);
                for x in [c[0] - root, c[0] + root] {
                    if x <= p[0] {
                        continue;
                    }
                    let v = [x - c[0], dy];
                    let (start, end) = if ccw {
                        (e.start, e.end)
                    } else {
                        (e.end, e.start)
                    };
                    let a = sub(start, c);
                    let b = sub(end, c);
                    let av = cross(a, v);
                    let vb = cross(v, b);
                    let contains = e.full
                        || if cross(a, b) > 0.0 {
                            av >= 0.0 && vb >= 0.0
                        } else if cross(a, b) < 0.0 {
                            av >= 0.0 || vb >= 0.0
                        } else {
                            av >= 0.0
                        };
                    if contains {
                        // Half-open endpoints count a crossing once across a join.
                        let at_start = av == 0.0 && a[0] * v[0] + a[1] * v[1] > 0.0;
                        let at_end = vb == 0.0 && b[0] * v[0] + b[1] * v[1] > 0.0;
                        if e.full || (!at_start || a[0] > 0.0) && (!at_end || b[0] < 0.0) {
                            crossings += 1;
                        }
                    }
                }
            } else if (e.start[1] > p[1]) != (e.end[1] > p[1]) {
                let x = e.start[0]
                    + (p[1] - e.start[1]) * (e.end[0] - e.start[0]) / (e.end[1] - e.start[1]);
                if x > p[0] {
                    crossings += 1;
                }
            }
        }
        crossings % 2 == 1
    }
}
