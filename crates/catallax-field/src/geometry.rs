//! Exact supporting surfaces, primitive trimming and cached numeric constants.

use catallax_numbers::{Ball, MPoly, Rational, Refusal, Sign};
use catallax_recipe::{ElementId, Operation, Origin, Role, StepId};
use catallax_shapes::{
    Cone, Cylinder, Plane, Pose, RPoint3, RVector3, Rotation, Sphere, Surface, Torus,
};

use crate::metric::{Kind, PRECISION, distance};
use crate::{Face, Location};

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
    pub polynomial: MPoly,
    terms: Vec<([u32; 3], Rational)>,
}

fn unit(axis: usize, sign: i64) -> RVector3 {
    let [x, y, z] = std::array::from_fn(|i| Rational::from_i64(if i == axis { sign } else { 0 }));
    RVector3::new(x, y, z)
}
fn at_z(z: Rational) -> RPoint3 {
    RPoint3::new(Rational::zero(), Rational::zero(), z)
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
        origin: RPoint3::new(x, y, z),
        normal: unit(axis, sign),
        ref_dir: unit((axis + 1) % 3, 1),
    })
}

pub(crate) fn compile(op: &Operation, step: &StepId) -> Result<(Leaf, Vec<Face>), Refusal> {
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
    let pose = Pose::new(rotation, RVector3::new(x, y, zz));
    let inverse = pose.inverse();
    let mut constraints = Vec::new();
    let mut faces = Vec::new();
    let mut add = |surface: Surface, role: Option<Role>| {
        let surface = surface.transformed(&pose);
        let polynomial = surface.implicit();
        let terms = polynomial
            .terms()
            .into_iter()
            .map(|(m, c)| {
                let e: [u32; 3] = m.exponents().try_into().expect("three coordinates");
                assert!(
                    e.iter().sum::<u32>() <= 4,
                    "analytic constraints are at most quartic"
                );
                (e, c)
            })
            .collect();
        constraints.push(Constraint {
            polynomial: polynomial.clone(),
            terms,
        });
        if let Some(role) = role {
            let gradient = std::array::from_fn(|i| polynomial.partial(i));
            let coefficient = polynomial
                .terms()
                .into_iter()
                .next()
                .expect("surface polynomial")
                .1;
            let support = polynomial.scale(&coefficient.recip().expect("nonzero coefficient"));
            faces.push(Face {
                origin: Origin {
                    step: step.clone(),
                    role,
                },
                surface,
                polynomial,
                support,
                gradient,
                leaf: 0,
            });
        }
    };
    match kind {
        Kind::Block => {
            for (axis, lo, hi) in [
                (0, side("x_min"), side("x_max")),
                (1, side("y_min"), side("y_max")),
                (2, Role::Start, Role::End),
            ] {
                add(plane(axis, z.clone(), -1), Some(lo));
                add(plane(axis, parameters[axis].clone(), 1), Some(hi));
            }
        }
        Kind::Round => {
            add(plane(2, z.clone(), -1), Some(Role::Start));
            add(
                plane(2, parameters[2].clone(), 1),
                positive(&parameters[1]).then_some(Role::End),
            );
            let difference = &parameters[0] - &parameters[1];
            let wall = if difference.is_zero() {
                Surface::Cylinder(Cylinder {
                    origin: RPoint3::origin(),
                    axis: unit(2, 1),
                    radius: parameters[0].clone(),
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
                    slope: &difference.abs() / &parameters[2],
                    ref_dir: unit(0, 1),
                })
            };
            add(wall, Some(side("wall")));
            parameters[3] = -difference;
            parameters[4] = &parameters[3] * &parameters[3] + &parameters[2] * &parameters[2];
        }
        Kind::Sphere => add(
            Surface::Sphere(Sphere {
                center: RPoint3::origin(),
                radius: parameters[0].clone(),
                axis: unit(2, 1),
                ref_dir: unit(0, 1),
            }),
            Some(side("wall")),
        ),
        Kind::Torus => add(
            Surface::Torus(Torus {
                center: RPoint3::origin(),
                major: parameters[0].clone(),
                minor: parameters[1].clone(),
                axis: unit(2, 1),
                ref_dir: unit(0, 1),
            }),
            Some(side("wall")),
        ),
    }
    Ok((
        Leaf::Primitive(Primitive {
            kind,
            fast_matrix: std::array::from_fn(|i| {
                std::array::from_fn(|j| inverse.rotation.matrix()[i][j].to_f64_nearest())
            }),
            fast_translation: [
                &inverse.translation.x,
                &inverse.translation.y,
                &inverse.translation.z,
            ]
            .map(Rational::to_f64_nearest),
            fast_parameters: std::array::from_fn(|i| parameters[i].to_f64_nearest()),
            balls: std::array::from_fn(|i| Ball::from_rational(&parameters[i], PRECISION)),
            parameters: {
                let [a, b, c, _, _] = parameters;
                [a, b, c]
            },
            inverse,
            constraints,
        }),
        faces,
    ))
}

impl Primitive {
    pub fn location(&self, point: &[Rational; 3]) -> Location {
        let mut boundary = false;
        for poly in &self.constraints {
            match poly.polynomial.eval(point).sign() {
                Sign::Positive => return Location::Outside,
                Sign::Zero => boundary = true,
                Sign::Negative => (),
            }
        }
        if boundary {
            Location::Boundary
        } else {
            Location::Inside
        }
    }
    pub fn proven_distance(&self, point: &RPoint3, location: Location) -> Ball {
        if location == Location::Boundary {
            return Ball::from_rational(&Rational::zero(), PRECISION);
        }
        let p = self.inverse.apply_point(point);
        let p = [p.x, p.y, p.z].map(|r| Ball::from_rational(&r, PRECISION));
        let d = distance(self.kind, &self.balls, &p);
        if self.kind == Kind::Round && location == Location::Inside {
            -d
        } else {
            d
        }
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

/// Sign at `point + epsilon * direction` for all sufficiently small positive
/// epsilon. Each constraint is a polynomial, so its first nonzero coefficient
/// settles the ray without choosing a tolerance or sampling a finite distance.
pub(crate) fn ray_location(
    leaf: &Primitive,
    point: &[Rational; 3],
    direction: &[Rational; 3],
) -> Location {
    let mut boundary = false;
    for poly in &leaf.constraints {
        match ray_sign(poly, point, direction) {
            Sign::Positive => return Location::Outside,
            Sign::Zero => boundary = true,
            Sign::Negative => (),
        }
    }
    if boundary {
        Location::Boundary
    } else {
        Location::Inside
    }
}

fn ray_sign(poly: &Constraint, point: &[Rational; 3], direction: &[Rational; 3]) -> Sign {
    ray_coefficients(poly, point, direction)
        .iter()
        .find(|c| !c.is_zero())
        .map_or(Sign::Zero, Rational::sign)
}

pub(crate) fn ray_coefficients(
    poly: &Constraint,
    point: &[Rational; 3],
    direction: &[Rational; 3],
) -> [Rational; 5] {
    let mut coefficients: [Rational; 5] = std::array::from_fn(|_| Rational::zero());
    for (exponents, coefficient) in &poly.terms {
        let mut term: [Rational; 5] = std::array::from_fn(|_| Rational::zero());
        term[0] = coefficient.clone();
        let mut degree = 0;
        for axis in 0..3 {
            for _ in 0..exponents[axis] {
                for i in (1..=degree + 1).rev() {
                    term[i] = (&term[i] * &point[axis]) + (&term[i - 1] * &direction[axis]);
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

impl Constraint {
    pub fn new(polynomial: MPoly) -> Self {
        let terms = polynomial
            .terms()
            .into_iter()
            .map(|(m, c)| (m.exponents().try_into().expect("three coordinates"), c))
            .collect();
        Self { polynomial, terms }
    }
}

pub(crate) fn face(surface: Surface, origin: Origin) -> Face {
    let polynomial = surface.implicit();
    let gradient = std::array::from_fn(|i| polynomial.partial(i));
    let coefficient = polynomial
        .terms()
        .into_iter()
        .next()
        .expect("surface polynomial")
        .1;
    let support = polynomial.scale(&coefficient.recip().expect("nonzero coefficient"));
    Face {
        origin,
        surface,
        polynomial,
        support,
        gradient,
        leaf: 0,
    }
}

#[derive(Clone, Debug)]
pub(crate) enum Leaf {
    Primitive(Primitive),
    Feature(crate::features::Feature),
}
impl Leaf {
    pub fn same_geometry(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Primitive(a), Self::Primitive(b)) => {
                a.kind == b.kind && a.parameters == b.parameters && a.inverse == b.inverse
            }
            _ => false,
        }
    }
    pub fn true_distance(&self) -> bool {
        matches!(self, Self::Primitive(_))
    }
    pub fn location(&self, p: &[Rational; 3]) -> Location {
        match self {
            Self::Primitive(a) => a.location(p),
            Self::Feature(a) => a.location(p),
        }
    }
    pub fn proven_distance(&self, p: &RPoint3, l: Location) -> Ball {
        match self {
            Self::Primitive(a) => a.proven_distance(p, l),
            Self::Feature(a) => a.proven_distance(p, l),
        }
    }
    pub fn fast_distance(&self, p: [f64; 3]) -> f64 {
        match self {
            Self::Primitive(a) => a.fast_distance(p),
            Self::Feature(a) => a.fast_distance(p),
        }
    }
    pub fn ray(&self, p: &[Rational; 3], d: &[Rational; 3]) -> Location {
        match self {
            Self::Primitive(a) => ray_location(a, p, d),
            Self::Feature(a) => a.ray(p, d),
        }
    }
    pub fn on_face(&self, role: &Role, p: &[Rational; 3]) -> Option<bool> {
        match self {
            Self::Primitive(_) => Some(true),
            Self::Feature(a) => a.on_face(role, p),
        }
    }
}

impl Primitive {
    pub(crate) fn sample_bounds(&self) -> crate::samples::Bounds {
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
