//! Exact supporting surfaces, primitive trimming and separate rounded caches.

use crate::metric::{Kind, PRECISION, distance};
use crate::{Face, Location};
use catallax_numbers::{
    Ball, Decider, Outcome, Rational, Real, RealPoly, Refusal, Sign, StrictSign,
};
use catallax_recipe::{ElementId, Operation, Origin, Role, StepId};
use catallax_shapes::{
    Cone, Cylinder, Plane, Point3, Pose, RPoint3, RVector3, Rotation, Sphere, Surface, Torus,
    Vector3,
};

#[derive(Clone, Debug)]
pub(crate) struct Primitive {
    pub kind: Kind,
    pub parameters: [Rational; 3],
    inverse: Pose,
    fast_matrix: [[f64; 3]; 3],
    fast_translation: [f64; 3],
    fast_parameters: [f64; 5],
    balls: [Ball; 5],
    constraints: Vec<Constraint>,
}
#[derive(Clone, Debug)]
pub(crate) struct Constraint {
    surface: Surface,
    pub(crate) polynomial: RealPoly,
    terms: Terms,
}
#[derive(Clone, Debug)]
enum Terms {
    Rational(Vec<([u32; 3], Rational)>),
    Algebraic(Vec<([u32; 3], Real)>),
}

pub(crate) fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    if let Some(value) = value.rational_ref() {
        return Ok(value.sign());
    }
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        Outcome::Refused(r) => Err(r),
    }
}
pub(crate) fn enclose(value: &Real, precision: u32, dec: &mut Decider) -> Result<Ball, Refusal> {
    if let Some(value) = value.rational_ref() {
        Ok(Ball::from_rational(value, precision))
    } else {
        dec.enclose_budgeted_at(value, precision)
    }
}
fn unit(axis: usize, sign: i64) -> Vector3 {
    let [x, y, z] = std::array::from_fn(|i| Rational::from_i64(if i == axis { sign } else { 0 }));
    RVector3::new(x, y, z).to_real()
}
fn at_z(z: Rational) -> Point3 {
    RPoint3::new(Rational::zero(), Rational::zero(), z).to_real()
}
fn side(name: &str) -> Role {
    Role::Side {
        segment: ElementId::new(name).expect("fixed role id"),
        path: None,
    }
}
fn plane(axis: usize, coordinate: Rational, sign: i64) -> Surface {
    let [x, y, z] = std::array::from_fn(|i| {
        if i == axis {
            coordinate.clone()
        } else {
            Rational::zero()
        }
    });
    Surface::Plane(Plane {
        origin: RPoint3::new(x, y, z).to_real(),
        normal: unit(axis, sign),
        ref_dir: unit((axis + 1) % 3, 1),
    })
}
pub(crate) fn constraint(surface: Surface) -> Constraint {
    let polynomial = surface.implicit();
    let terms = if let Some(p) = polynomial.rational() {
        Terms::Rational(
            p.terms()
                .into_iter()
                .map(|(m, c)| (m.exponents().try_into().expect("three coordinates"), c))
                .collect(),
        )
    } else {
        Terms::Algebraic(
            polynomial
                .terms()
                .into_iter()
                .map(|(m, c)| (m.exponents().try_into().expect("three coordinates"), c))
                .collect(),
        )
    };
    Constraint {
        surface,
        polynomial,
        terms,
    }
}
impl Constraint {
    pub(crate) fn transform(&mut self, pose: &Pose) {
        *self = constraint(self.surface.transformed(pose));
    }
}
fn support_geometry(
    polynomial: &RealPoly,
    dec: &mut Decider,
) -> Result<(RealPoly, [RealPoly; 3]), Refusal> {
    let mut first = None;
    for (_, c) in polynomial.terms() {
        if sign(&c, dec)? != Sign::Zero {
            first = Some(c);
            break;
        }
    }
    let coefficient = first.ok_or_else(|| {
        Refusal::new(
            "numbers.zero-polynomial",
            "surface support polynomial vanishes identically",
        )
    })?;
    let support = polynomial.scale(&(&Real::one() / &coefficient));
    Ok((support, std::array::from_fn(|i| polynomial.partial(i))))
}
pub(crate) fn face_geometry(
    face: &mut Face,
    surface: Surface,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let polynomial = surface.implicit();
    let (support, gradient) = support_geometry(&polynomial, dec)?;
    face.support = support;
    face.gradient = gradient;
    face.polynomial = polynomial;
    face.surface = surface;
    Ok(())
}

pub(crate) fn compile(
    op: &Operation,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(Leaf, Vec<Face>), Refusal> {
    let q = |n: &catallax_recipe::Number| n.value().clone();
    let z = Rational::zero();
    let (kind, mut parameters, placement) = match op {
        Operation::Block { size, pose } => (
            Kind::Block,
            [q(&size[0]), q(&size[1]), q(&size[2]), z.clone(), z.clone()],
            pose,
        ),
        Operation::Cylinder {
            radius,
            height,
            pose,
        } => (
            Kind::Round,
            [q(radius), q(radius), q(height), z.clone(), z.clone()],
            pose,
        ),
        Operation::Cone {
            base_radius,
            top_radius,
            height,
            pose,
        } => (
            Kind::Round,
            [
                q(base_radius),
                q(top_radius),
                q(height),
                z.clone(),
                z.clone(),
            ],
            pose,
        ),
        Operation::Sphere { radius, pose } => (
            Kind::Sphere,
            [q(radius), z.clone(), z.clone(), z.clone(), z.clone()],
            pose,
        ),
        Operation::Torus { major, minor, pose } => (
            Kind::Torus,
            [q(major), q(minor), z.clone(), z.clone(), z.clone()],
            pose,
        ),
        _ => {
            return Err(Refusal::new(
                "field.unsupported-operation",
                format!("{step:?}: {}", op.name()),
            ));
        }
    };
    let positive = |r: &Rational| r.sign() == Sign::Positive;
    let valid = match kind {
        Kind::Block => parameters[..3].iter().all(positive),
        Kind::Round => {
            positive(&parameters[0])
                && parameters[1].sign() != Sign::Negative
                && positive(&parameters[2])
        }
        Kind::Sphere => positive(&parameters[0]),
        Kind::Torus => positive(&parameters[1]) && parameters[0] > parameters[1],
    };
    if !valid {
        return Err(Refusal::new(
            "field.impossible-primitive",
            format!("invalid dimensions at {step:?}"),
        ));
    }
    let [w, x, y, zz] = placement.rotation.clone().map(|n| q(&n));
    let rotation = Rotation::from_quaternion(w, x, y, zz)?;
    let [x, y, zz] = placement.translation.0.clone().map(|n| q(&n));
    let pose = Pose::new(rotation, RVector3::new(x, y, zz).to_real(), dec)?;
    let inverse = pose.inverse();
    let mut constraints = Vec::new();
    let mut faces = Vec::new();
    let mut add = |surface: Surface, role: Option<Role>| -> Result<(), Refusal> {
        let surface = surface.transformed(&pose);
        let constraint = constraint(surface.clone());
        if let Some(role) = role {
            let (support, gradient) = support_geometry(&constraint.polynomial, dec)?;
            faces.push(Face {
                origin: Origin {
                    step: step.clone(),
                    role,
                },
                surface,
                polynomial: constraint.polynomial.clone(),
                support,
                gradient,
                leaf: 0,
            });
        }
        constraints.push(constraint);
        Ok(())
    };
    match kind {
        Kind::Block => {
            for (axis, lo, hi) in [
                (0, side("x_min"), side("x_max")),
                (1, side("y_min"), side("y_max")),
                (2, Role::Start, Role::End),
            ] {
                add(plane(axis, z.clone(), -1), Some(lo))?;
                add(plane(axis, parameters[axis].clone(), 1), Some(hi))?;
            }
        }
        Kind::Round => {
            add(plane(2, z.clone(), -1), Some(Role::Start))?;
            add(
                plane(2, parameters[2].clone(), 1),
                positive(&parameters[1]).then_some(Role::End),
            )?;
            let difference = &parameters[0] - &parameters[1];
            let wall = if difference.is_zero() {
                Surface::Cylinder(Cylinder {
                    origin: Point3::origin(),
                    axis: unit(2, 1),
                    radius: parameters[0].clone().into(),
                    ref_dir: unit(0, 1),
                })
            } else {
                Surface::Cone(Cone {
                    apex: at_z(&(&parameters[2] * &parameters[0]) / &difference),
                    axis: unit(
                        2,
                        if difference.sign() == Sign::Positive {
                            -1
                        } else {
                            1
                        },
                    ),
                    slope: (&difference.abs() / &parameters[2]).into(),
                    ref_dir: unit(0, 1),
                })
            };
            add(wall, Some(side("wall")))?;
            parameters[3] = -difference;
            parameters[4] = &parameters[3] * &parameters[3] + &parameters[2] * &parameters[2];
        }
        Kind::Sphere => add(
            Surface::Sphere(Sphere {
                center: Point3::origin(),
                radius: parameters[0].clone().into(),
                axis: unit(2, 1),
                ref_dir: unit(0, 1),
            }),
            Some(side("wall")),
        )?,
        Kind::Torus => add(
            Surface::Torus(Torus {
                center: Point3::origin(),
                major: parameters[0].clone().into(),
                minor: parameters[1].clone().into(),
                axis: unit(2, 1),
                ref_dir: unit(0, 1),
            }),
            Some(side("wall")),
        )?,
    }
    let fast_parameters = parameters.each_ref().map(Rational::to_f64_nearest);
    let balls = parameters
        .each_ref()
        .map(|p| Ball::from_rational(p, PRECISION));
    let [a, b, c, _, _] = parameters;
    let leaf = Primitive {
        kind,
        fast_matrix: fast_matrix(&inverse),
        fast_translation: fast_translation(&inverse),
        fast_parameters,
        balls,
        parameters: [a, b, c],
        inverse,
        constraints,
    };
    Ok((Leaf::Primitive(leaf), faces))
}
pub(crate) fn fast_matrix(inverse: &Pose) -> [[f64; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| inverse.rotation.matrix()[i][j].to_f64()))
}
pub(crate) fn fast_translation(inverse: &Pose) -> [f64; 3] {
    [
        &inverse.translation.x,
        &inverse.translation.y,
        &inverse.translation.z,
    ]
    .map(Real::to_f64)
}

impl Primitive {
    pub fn same_geometry(&self, other: &Self, dec: &mut Decider) -> Result<bool, Refusal> {
        Ok(self.kind == other.kind
            && self.parameters == other.parameters
            && self.inverse.same_motion(&other.inverse, dec)?)
    }
    fn transform(&mut self, pose: &Pose) {
        self.inverse = pose.inverse().then(&self.inverse);
        self.fast_matrix = fast_matrix(&self.inverse);
        self.fast_translation = fast_translation(&self.inverse);
        for c in &mut self.constraints {
            *c = constraint(c.surface.transformed(pose));
        }
    }
    pub fn location(&self, point: &[Rational; 3], dec: &mut Decider) -> Result<Location, Refusal> {
        let mut boundary = false;
        for constraint in &self.constraints {
            let value = if let Some(p) = constraint.polynomial.rational() {
                p.eval(point).sign()
            } else {
                sign(
                    &constraint
                        .polynomial
                        .eval_real(&point.clone().map(Real::from)),
                    dec,
                )?
            };
            match value {
                Sign::Positive => return Ok(Location::Outside),
                Sign::Zero => boundary = true,
                Sign::Negative => (),
            }
        }
        Ok(if boundary {
            Location::Boundary
        } else {
            Location::Inside
        })
    }
    pub fn proven_distance(
        &self,
        point: &RPoint3,
        location: Location,
        dec: &mut Decider,
    ) -> Result<Ball, Refusal> {
        if location == Location::Boundary {
            return Ok(Ball::from_rational(&Rational::zero(), PRECISION));
        }
        let p = self.inverse.apply_point(&point.to_real());
        let p = [
            enclose(&p.x, PRECISION, dec)?,
            enclose(&p.y, PRECISION, dec)?,
            enclose(&p.z, PRECISION, dec)?,
        ];
        let d = distance(self.kind, &self.balls, &p);
        Ok(
            if self.kind == Kind::Round && location == Location::Inside {
                -d
            } else {
                d
            },
        )
    }
    pub fn fast_distance(&self, point: [f64; 3]) -> f64 {
        if self
            .fast_parameters
            .iter()
            .chain(self.fast_translation.iter())
            .any(|v| !v.is_finite())
            || (self.kind == Kind::Round && self.fast_parameters[4] <= 0.0)
        {
            return f64::NAN;
        }
        let p = std::array::from_fn(|i| {
            ((self.fast_matrix[i][0] * point[0] + self.fast_matrix[i][1] * point[1])
                + self.fast_matrix[i][2] * point[2])
                + self.fast_translation[i]
        });
        let d = distance(self.kind, &self.fast_parameters, &p);
        if self.kind != Kind::Round {
            return d;
        }
        let [rb, _, h, dr, _] = self.fast_parameters;
        let radius = rb + (p[2] / h) * dr;
        if p[2] >= 0.0 && p[2] <= h && p[0] * p[0] + p[1] * p[1] <= radius * radius {
            -d
        } else {
            d
        }
    }
}

/// Exact positive-infinitesimal ray classification, retaining coefficient proofs.
pub(crate) fn ray_location(
    leaf: &Primitive,
    point: &[Rational; 3],
    direction: &Vector3,
    dec: &mut Decider,
) -> Result<Location, Refusal> {
    let mut boundary = false;
    for poly in &leaf.constraints {
        let value = if let (Terms::Rational(terms), Some(direction)) =
            (&poly.terms, direction.rational_components())
        {
            rational_ray_coefficients(terms, point, direction)
                .into_iter()
                .find(|c| !c.is_zero())
                .map_or(Sign::Zero, |c| c.sign())
        } else {
            algebraic_ray_sign(poly, point, direction, dec)?
        };
        match value {
            Sign::Positive => return Ok(Location::Outside),
            Sign::Zero => boundary = true,
            Sign::Negative => (),
        }
    }
    Ok(if boundary {
        Location::Boundary
    } else {
        Location::Inside
    })
}
fn rational_ray_coefficients(
    terms: &[([u32; 3], Rational)],
    point: &[Rational; 3],
    direction: [&Rational; 3],
) -> [Rational; 5] {
    let mut coefficients: [Rational; 5] = std::array::from_fn(|_| Rational::zero());
    for (exponents, coefficient) in terms {
        let mut term: [Rational; 5] = std::array::from_fn(|_| Rational::zero());
        term[0] = coefficient.clone();
        let mut degree = 0;
        for axis in 0..3 {
            for _ in 0..exponents[axis] {
                for i in (1..=degree + 1).rev() {
                    term[i] = (&term[i] * &point[axis]) + (&term[i - 1] * direction[axis]);
                }
                term[0] = &term[0] * &point[axis];
                degree += 1;
            }
        }
        for i in 0..=degree {
            coefficients[i] = &coefficients[i] + &term[i];
        }
    }
    coefficients
}

pub(crate) fn face(surface: Surface, origin: Origin, dec: &mut Decider) -> Result<Face, Refusal> {
    let polynomial = surface.implicit();
    let (support, gradient) = support_geometry(&polynomial, dec)?;
    Ok(Face {
        origin,
        surface,
        polynomial,
        support,
        gradient,
        leaf: 0,
    })
}

#[derive(Clone, Debug)]
pub(crate) enum Leaf {
    Primitive(Primitive),
    Feature(crate::features::Feature),
}
impl Leaf {
    pub fn same_geometry(&self, other: &Self, dec: &mut Decider) -> Result<bool, Refusal> {
        match (self, other) {
            (Self::Primitive(a), Self::Primitive(b)) => a.same_geometry(b, dec),
            _ => Ok(false),
        }
    }
    pub fn transform(&mut self, pose: &Pose) {
        match self {
            Self::Primitive(a) => a.transform(pose),
            Self::Feature(a) => a.transform(pose),
        }
    }
    pub fn true_distance(&self) -> bool {
        matches!(self, Self::Primitive(_))
    }
    pub fn location(&self, p: &[Rational; 3], dec: &mut Decider) -> Result<Location, Refusal> {
        match self {
            Self::Primitive(a) => a.location(p, dec),
            Self::Feature(a) => a.location(p, dec),
        }
    }
    pub fn proven_distance(
        &self,
        p: &RPoint3,
        l: Location,
        dec: &mut Decider,
    ) -> Result<Ball, Refusal> {
        match self {
            Self::Primitive(a) => a.proven_distance(p, l, dec),
            Self::Feature(a) => a.proven_distance(p, l, dec),
        }
    }
    pub fn fast_distance(&self, p: [f64; 3]) -> f64 {
        match self {
            Self::Primitive(a) => a.fast_distance(p),
            Self::Feature(a) => a.fast_distance(p),
        }
    }
    pub fn ray(
        &self,
        p: &[Rational; 3],
        d: &Vector3,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        match self {
            Self::Primitive(a) => ray_location(a, p, d, dec),
            Self::Feature(a) => a.ray(p, d, dec),
        }
    }
    pub fn on_face(
        &self,
        role: &Role,
        p: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Option<bool>, Refusal> {
        match self {
            Self::Primitive(_) => Ok(Some(true)),
            Self::Feature(a) => a.on_face(role, p, dec).map(Some),
        }
    }
}

impl Primitive {
    pub(crate) fn sample_bounds(&self) -> Result<crate::samples::Bounds, Refusal> {
        let zero = Rational::zero();
        let p = &self.parameters;
        let local = match self.kind {
            Kind::Block => std::array::from_fn(|i| [zero.clone(), p[i].clone()]),
            Kind::Round => {
                let radius = p[0].clone().max(p[1].clone());
                [
                    [-radius.clone(), radius.clone()],
                    [-radius.clone(), radius],
                    [zero, p[2].clone()],
                ]
            }
            Kind::Sphere => std::array::from_fn(|_| [-p[0].clone(), p[0].clone()]),
            Kind::Torus => {
                let radius = &p[0] + &p[1];
                [
                    [-radius.clone(), radius.clone()],
                    [-radius.clone(), radius],
                    [-p[1].clone(), p[1].clone()],
                ]
            }
        };
        crate::samples::transform(&local, &self.inverse.inverse())
    }
}
fn algebraic_ray_sign(
    poly: &Constraint,
    point: &[Rational; 3],
    direction: &Vector3,
    dec: &mut Decider,
) -> Result<Sign, Refusal> {
    match ray_coefficients(poly, point, direction) {
        RayCoefficients::Rational(c) => Ok(c
            .into_iter()
            .find(|c| !c.is_zero())
            .map_or(Sign::Zero, |c| c.sign())),
        RayCoefficients::Algebraic(c) => {
            for coefficient in c {
                let s = sign(&coefficient, dec)?;
                if s != Sign::Zero {
                    return Ok(s);
                }
            }
            Ok(Sign::Zero)
        }
    }
}
pub(crate) enum RayCoefficients {
    Rational([Rational; 5]),
    Algebraic([Real; 5]),
}
pub(crate) fn ray_coefficients(
    poly: &Constraint,
    point: &[Rational; 3],
    direction: &Vector3,
) -> RayCoefficients {
    if let (Terms::Rational(terms), Some(d)) = (&poly.terms, direction.rational_components()) {
        return RayCoefficients::Rational(rational_ray_coefficients(terms, point, d));
    }
    let point = point.clone().map(Real::from);
    let direction = [&direction.x, &direction.y, &direction.z];
    let mut coefficients: [Real; 5] = std::array::from_fn(|_| Real::zero());
    let mut accumulate = |exponents: &[u32; 3], coefficient: Real| {
        let mut term: [Real; 5] = std::array::from_fn(|_| Real::zero());
        term[0] = coefficient;
        let mut degree = 0;
        for axis in 0..3 {
            for _ in 0..exponents[axis] {
                for i in (1..=degree + 1).rev() {
                    term[i] = (&term[i] * &point[axis]) + (&term[i - 1] * direction[axis]);
                }
                term[0] = &term[0] * &point[axis];
                degree += 1;
            }
        }
        for i in 0..=degree {
            coefficients[i] = &coefficients[i] + &term[i];
        }
    };
    match &poly.terms {
        Terms::Algebraic(terms) => {
            for (exponents, coefficient) in terms {
                accumulate(exponents, coefficient.clone());
            }
        }
        Terms::Rational(terms) => {
            for (exponents, coefficient) in terms {
                accumulate(exponents, coefficient.clone().into());
            }
        }
    }
    RayCoefficients::Algebraic(coefficients)
}
