use crate::math::*;
use crate::{Component, Contact, CurveForm, SurfacePair};
use catallax_numbers::{
    Answer, Answered, Decider, Fact, FactKind, FactOrigin, Rational, Real, Refusal, Sign,
};
use catallax_shapes::{
    Circle, Curve, Cylinder, Ellipse, EllipticCylinder, GeomKey, GeomStore, IntersectionCurve,
    Line, Meet, MeetOf, Plane, Point3, PointDef, Position, Sphere, Surface, SurfaceId, Vector3,
};
use std::{cell::RefCell, collections::VecDeque, sync::Arc};

mod conics;

#[derive(Clone, Debug)]
enum DividerValue {
    Plane {
        origin: Point3,
        vector: Vector3,
    },
    Radius {
        origin: Point3,
        axis: Vector3,
        inverse_norm_sq: Real,
        radius_sq: Real,
    },
}
#[derive(Clone, Debug)]
pub(super) struct Divider {
    value: DividerValue,
    pub length_lower: Rational,
    // The older quartic plans store a canonical direction and reverse it
    // across the divider. Conic plans store each component's own direction.
    pub canonical_direction: bool,
}
impl Divider {
    fn new(origin: Point3, vector: Vector3) -> Self {
        if let Some(v) = vector.as_rational() {
            let lower = v.x.abs().max(v.y.abs()).max(v.z.abs());
            return Self::linear(origin, vector, lower, true);
        }
        // Callers have proved this direction nonzero.
        let vector = vector.scaled(&(&Real::one() / &vector.norm_squared().sqrt()));
        Self::linear(origin, vector, r(1), true)
    }
    fn linear(
        origin: Point3,
        vector: Vector3,
        length_lower: Rational,
        canonical_direction: bool,
    ) -> Self {
        Self {
            value: DividerValue::Plane { origin, vector },
            length_lower,
            canonical_direction,
        }
    }
    fn conic(origin: Point3, vector: Vector3) -> Self {
        let mut divider = Self::new(origin, vector);
        divider.canonical_direction = false;
        divider
    }
    fn radial(origin: Point3, axis: Vector3, radius_sq: Real, radius_lower: Rational) -> Self {
        let inverse_norm_sq = s(1) / axis.norm_squared();
        Self {
            value: DividerValue::Radius {
                origin,
                axis,
                inverse_norm_sq,
                radius_sq,
            },
            length_lower: radius_lower,
            canonical_direction: false,
        }
    }
    pub fn value(&self, p: &Point3) -> Real {
        match &self.value {
            DividerValue::Plane { origin, vector } => p.minus(origin).dot(vector),
            DividerValue::Radius {
                origin,
                axis,
                inverse_norm_sq,
                radius_sq,
            } => {
                let d = p.minus(origin);
                let radial = d.norm_squared() - &d.dot(axis).square() * inverse_norm_sq;
                &radial - radius_sq
            }
        }
    }
}
#[derive(Clone, Debug)]
pub(super) enum Frame {
    CylinderSphere {
        origin: Point3,
        e1: Vector3,
        e2: Vector3,
        e3: Vector3,
        radius: Real,
        d: Real,
        k: Real,
        mode: u8,
        sense: i64,
    },
    Cylinders {
        origin: Point3,
        e1: Vector3,
        e2: Vector3,
        e3: Vector3,
        big: Real,
        small: Real,
        delta: Real,
        alpha: Real,
        gamma: Real,
        mode: u8,
        sense: i64,
    },
}
#[derive(Clone, Debug)]
pub(super) struct Spec {
    pub form: CurveForm,
    pub witness: Option<Point3>,
    pub contact: Contact,
    pub singular: Vec<usize>,
}
#[derive(Clone, Debug)]
pub(super) struct Plan {
    pub terminal: Option<SurfacePair>,
    pub points: Vec<Point3>,
    pub curves: Vec<Spec>,
    pub divider: Option<Divider>,
    pub frame: Option<Frame>,
}
impl Plan {
    fn empty() -> Self {
        Self {
            terminal: None,
            points: vec![],
            curves: vec![],
            divider: None,
            frame: None,
        }
    }
    fn disjoint() -> Self {
        Self {
            terminal: Some(SurfacePair::Disjoint),
            ..Self::empty()
        }
    }
    fn coincident(same_sense: bool) -> Self {
        Self {
            terminal: Some(SurfacePair::Coincident { same_sense }),
            ..Self::empty()
        }
    }
    fn point(p: Point3) -> Self {
        Self {
            points: vec![p],
            ..Self::empty()
        }
    }
    fn form(form: CurveForm, contact: Contact) -> Self {
        Self {
            curves: vec![Spec {
                form,
                witness: None,
                contact,
                singular: vec![],
            }],
            ..Self::empty()
        }
    }
    fn branches(
        witnesses: Vec<Point3>,
        node: Option<Point3>,
        divider: Option<Divider>,
        frame: Frame,
    ) -> Self {
        let singular = if node.is_some() { vec![0] } else { vec![] };
        Self {
            terminal: None,
            points: node.into_iter().collect(),
            curves: witnesses
                .into_iter()
                .map(|w| Spec {
                    form: CurveForm::Branch,
                    witness: Some(w),
                    contact: Contact::Crossing,
                    singular: singular.clone(),
                })
                .collect(),
            divider,
            frame: Some(frame),
        }
    }
}
fn s(value: i64) -> Real {
    Real::from_i64(value)
}
fn line(origin: Point3, dir: Vector3) -> CurveForm {
    CurveForm::Line { origin, dir }
}
fn circle(center: Point3, normal: Vector3, radius_sq: Real, ref_dir: Vector3) -> CurveForm {
    CurveForm::Circle {
        center,
        normal,
        radius_sq,
        ref_dir,
    }
}
fn is_zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(value, dec)? == Sign::Zero)
}
fn compare(a: &Real, b: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    Ok(sign(&(a - b), dec)?)
}
fn vector_zero(value: &Vector3, dec: &mut Decider) -> Result<bool, Refusal> {
    is_zero(&value.norm_squared(), dec)
}

fn two_lines(base: Point3, v: Vector3, h2: Real, dirs: [Vector3; 2], first_positive: bool) -> Plan {
    let root = root_real(&h2);
    let offset = v.scaled(&root);
    let origins = [base.plus(&offset), base.plus(&offset.scaled(&s(-1)))];
    let witness = root.rational_ref().is_none();
    Plan {
        curves: origins
            .into_iter()
            .zip(dirs)
            .map(|(origin, dir)| Spec {
                witness: witness.then(|| origin.clone()),
                form: line(origin, dir),
                contact: Contact::Crossing,
                singular: vec![],
            })
            .collect(),
        divider: Some(Divider::new(
            base,
            if first_positive { v } else { v.scaled(&s(-1)) },
        )),
        ..Plan::empty()
    }
}
fn planes(a: &Plane, b: &Plane, dec: &mut Decider) -> Result<Plan, Refusal> {
    let m = a.normal.cross(&b.normal);
    if vector_zero(&m, dec)? {
        return Ok(
            if is_zero(&b.origin.minus(&a.origin).dot(&a.normal), dec)? {
                Plan::coincident(sign(&a.normal.dot(&b.normal), dec)? == Sign::Positive)
            } else {
                Plan::disjoint()
            },
        );
    }
    let scale = b.origin.minus(&a.origin).dot(&b.normal) / m.norm_squared();
    let origin = a.origin.plus(&a.normal.cross(&m).scaled(&(-scale)));
    Ok(Plan::form(line(origin, m), Contact::Crossing))
}
fn plane_cylinder(
    p: &Plane,
    c: &Cylinder,
    sigma: Real,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let dot = p.normal.dot(&c.axis);
    let r2 = c.radius_sq.clone();
    if !is_zero(&dot, dec)? {
        let center = c.origin.plus(
            &c.axis
                .scaled(&(&p.origin.minus(&c.origin).dot(&p.normal) / &dot)),
        );
        let normal = p.normal.scaled(&sigma);
        if vector_zero(&c.axis.cross(&p.normal), dec)? {
            return Ok(Plan::form(
                circle(center, normal, r2, c.ref_dir.clone()),
                Contact::Crossing,
            ));
        }
        return Ok(Plan::form(
            CurveForm::Ellipse {
                center,
                normal,
                major_dir: c
                    .axis
                    .minus(&p.normal.scaled(&(&dot / &p.normal.norm_squared()))),
                major_sq: r2.clone() * p.normal.norm_squared() * c.axis.norm_squared()
                    / dot.square(),
                minor_sq: r2,
            },
            Contact::Crossing,
        ));
    }
    let k = c.origin.minus(&p.origin).dot(&p.normal);
    let e2 = k.square() / p.normal.norm_squared();
    let relation = compare(&e2, &r2, dec)?;
    if relation == Sign::Positive {
        return Ok(Plan::disjoint());
    }
    let foot = c
        .origin
        .plus(&p.normal.scaled(&(-k / p.normal.norm_squared())));
    if relation == Sign::Zero {
        return Ok(Plan::form(line(foot, c.axis.clone()), Contact::Tangent));
    }
    let v = p.normal.cross(&c.axis);
    let h2 = (r2 - e2) / v.norm_squared();
    Ok(two_lines(
        foot,
        v,
        h2,
        [c.axis.scaled(&(-&sigma)), c.axis.scaled(&sigma)],
        true,
    ))
}
fn plane_sphere(
    p: &Plane,
    sphere: &Sphere,
    sigma: Real,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let k = sphere.center.minus(&p.origin).dot(&p.normal);
    let e2 = k.square() / p.normal.norm_squared();
    let r2 = sphere.radius.square();
    let relation = compare(&e2, &r2, dec)?;
    if relation == Sign::Positive {
        return Ok(Plan::disjoint());
    }
    let foot = sphere
        .center
        .plus(&p.normal.scaled(&(-k / p.normal.norm_squared())));
    if relation == Sign::Zero {
        return Ok(Plan::point(foot));
    }
    Ok(Plan::form(
        circle(foot, p.normal.scaled(&sigma), r2 - e2, p.ref_dir.clone()),
        Contact::Crossing,
    ))
}
fn spheres(a: &Sphere, b: &Sphere, dec: &mut Decider) -> Result<Plan, Refusal> {
    let u = b.center.minus(&a.center);
    let d2 = u.norm_squared();
    if is_zero(&d2, dec)? {
        return Ok(if is_zero(&(&a.radius - &b.radius), dec)? {
            Plan::coincident(true)
        } else {
            Plan::disjoint()
        });
    }
    let sum = &a.radius + &b.radius;
    let diff = &a.radius - &b.radius;
    let outer = compare(&d2, &sum.square(), dec)?;
    let inner = compare(&d2, &diff.square(), dec)?;
    if outer == Sign::Positive || inner == Sign::Negative {
        return Ok(Plan::disjoint());
    }
    if outer == Sign::Zero {
        return Ok(Plan::point(a.center.plus(&u.scaled(&(&a.radius / &sum)))));
    }
    if inner == Sign::Zero {
        return Ok(Plan::point(a.center.plus(&u.scaled(&(&a.radius / &diff)))));
    }
    let t = (&d2 + &a.radius.square() - b.radius.square()) / (s(2) * d2.clone());
    Ok(Plan::form(
        circle(
            a.center.plus(&u.scaled(&t)),
            u.clone(),
            a.radius.square() - t.square() * d2,
            perp(&u, dec)?,
        ),
        Contact::Crossing,
    ))
}
fn parallel_cylinders(a: &Cylinder, b: &Cylinder, dec: &mut Decider) -> Result<Plan, Refusal> {
    let d = b.origin.minus(&a.origin);
    let u = d.minus(&a.axis.scaled(&(d.dot(&a.axis) / a.axis.norm_squared())));
    let d2 = u.norm_squared();
    if is_zero(&d2, dec)? {
        return Ok(if is_zero(&(&a.radius_sq - &b.radius_sq), dec)? {
            Plan::coincident(true)
        } else {
            Plan::disjoint()
        });
    }
    let a_radius = a.radius();
    let b_radius = b.radius();
    let sum = &a_radius + &b_radius;
    let diff = &a_radius - &b_radius;
    let outer = compare(&d2, &sum.square(), dec)?;
    let inner = compare(&d2, &diff.square(), dec)?;
    if outer == Sign::Positive || inner == Sign::Negative {
        return Ok(Plan::disjoint());
    }
    if outer == Sign::Zero {
        return Ok(Plan::form(
            line(
                a.origin.plus(&u.scaled(&(&a_radius / &sum))),
                a.axis.clone(),
            ),
            Contact::Tangent,
        ));
    }
    if inner == Sign::Zero {
        return Ok(Plan::form(
            line(
                a.origin.plus(&u.scaled(&(&a_radius / &diff))),
                a.axis.clone(),
            ),
            Contact::Tangent,
        ));
    }
    let t = (&d2 + &a.radius_sq - b.radius_sq.clone()) / (s(2) * d2.clone());
    let v = a.axis.cross(&u);
    let h2 = (a.radius_sq.clone() - t.square() * d2) / v.norm_squared();
    Ok(two_lines(
        a.origin.plus(&u.scaled(&t)),
        v,
        h2,
        [a.axis.clone(), a.axis.scaled(&s(-1))],
        true,
    ))
}
fn cylinders(
    a: &Cylinder,
    b: &Cylinder,
    rational: bool,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    if vector_zero(&a.axis.cross(&b.axis), dec)? {
        return parallel_cylinders(a, b, dec);
    }
    let a_is_big = compare(&a.radius_sq, &b.radius_sq, dec)? != Sign::Negative;
    let (big, small) = if a_is_big { (a, b) } else { (b, a) };
    let m = big.axis.cross(&small.axis);
    let d = small.origin.minus(&big.origin);
    let fb = big.origin.plus(
        &big.axis
            .scaled(&(d.cross(&small.axis).dot(&m) / m.norm_squared())),
    );
    let fs = small.origin.plus(
        &small
            .axis
            .scaled(&(d.cross(&big.axis).dot(&m) / m.norm_squared())),
    );
    let delta2 = d.dot(&m).square() / m.norm_squared();
    let big_radius = big.radius();
    let small_radius = small.radius();
    let sum = &big_radius + &small_radius;
    let diff = &big_radius - &small_radius;
    match compare(&delta2, &sum.square(), dec)? {
        Sign::Positive => return Ok(Plan::disjoint()),
        Sign::Zero => {
            return Ok(Plan::point(
                fb.plus(&fs.minus(&fb).scaled(&(&big_radius / &sum))),
            ));
        }
        Sign::Negative => {}
    }
    let zero_delta = is_zero(&delta2, dec)?;
    if zero_delta && is_zero(&diff, dec)? {
        let norm_product = big.axis.norm_squared() * small.axis.norm_squared();
        let k = if rational {
            // Preserve the existing rational CC6 coverage boundary and answers.
            rational_sqrt(
                norm_product
                    .rational_ref()
                    .expect("rational cylinder frame"),
            )
            .ok_or_else(refused_pair)?
            .into()
        } else {
            root_real(&norm_product)
        };
        let nm = big
            .axis
            .scaled(&small.axis.norm_squared())
            .minus(&small.axis.scaled(&k));
        let np = big
            .axis
            .scaled(&small.axis.norm_squared())
            .plus(&small.axis.scaled(&k));
        let cosine = big.axis.dot(&small.axis) / k;
        let offset = m.scaled(&root_real(&(&big.radius_sq / &m.norm_squared())));
        let points = vec![fb.plus(&offset), fb.plus(&offset.scaled(&s(-1)))];
        return Ok(Plan {
            points,
            curves: [
                (nm.clone(), np.clone(), s(1) - cosine.clone()),
                (np, nm, s(1) + cosine),
            ]
            .into_iter()
            .map(|(normal, major_dir, denom)| Spec {
                form: CurveForm::Ellipse {
                    center: fb.clone(),
                    normal,
                    major_dir,
                    major_sq: s(2) * big.radius_sq.clone() / denom,
                    minor_sq: big.radius_sq.clone(),
                },
                witness: None,
                contact: Contact::Crossing,
                singular: vec![0, 1],
            })
            .collect(),
            ..Plan::empty()
        });
    }
    let delta = root_real(&delta2);
    let radial = if zero_delta { m } else { fs.minus(&fb) };
    let transverse = radial.cross(&big.axis);
    let e3 = unit(&big.axis);
    let e2 = unit(&radial);
    let e1 = unit(&transverse);
    let alpha_numerator = small.axis.dot(&transverse);
    let alpha =
        &alpha_numerator / &root_real(&(small.axis.norm_squared() * transverse.norm_squared()));
    let sigma = if a_is_big { 1 } else { -1 };
    let sense = -sigma
        * if sign(&alpha_numerator, dec)? == Sign::Positive {
            1
        } else {
            -1
        };
    let gamma = small.axis.dot(&big.axis)
        / root_real(&(small.axis.norm_squared() * big.axis.norm_squared()));
    let mode = match compare(&delta2, &diff.square(), dec)? {
        Sign::Positive => 3,
        Sign::Zero => 4,
        Sign::Negative => 5,
    };
    let frame = Frame::Cylinders {
        origin: fb.clone(),
        e1: e1.clone(),
        e2: e2.clone(),
        e3: e3.clone(),
        big: big_radius.clone(),
        small: small_radius.clone(),
        delta: delta.clone(),
        alpha: alpha.clone(),
        gamma: gamma.clone(),
        mode,
        sense,
    };
    if mode == 3 {
        let z = (&small.radius_sq - &(&big_radius - &delta).square()).sqrt() / alpha;
        let witness = fb.plus(&e2.scaled(&big_radius)).plus(&e3.scaled(&z));
        return Ok(Plan::branches(vec![witness], None, None, frame));
    }
    let x = root_real(&(&big.radius_sq - &delta2));
    let wplus = fs
        .plus(&e1.scaled(&x))
        .plus(&e3.scaled(&((&gamma * &x + small_radius.clone()) / alpha.clone())));
    let wminus = fs
        .plus(&e1.scaled(&(-&x)))
        .plus(&e3.scaled(&((&small_radius - &(&gamma * &x)) / alpha)));
    let node = (mode == 4).then(|| fb.plus(&fs.minus(&fb).scaled(&(&big_radius / &diff))));
    Ok(Plan::branches(
        vec![wplus, wminus],
        node,
        Some(Divider::linear(fb, e1, r(1), true)),
        frame,
    ))
}
fn cylinder_sphere(
    c: &Cylinder,
    sphere: &Sphere,
    sigma: Real,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let foot = c.origin.plus(
        &c.axis
            .scaled(&(sphere.center.minus(&c.origin).dot(&c.axis) / c.axis.norm_squared())),
    );
    let offset = sphere.center.minus(&foot);
    let d2 = offset.norm_squared();
    let r2 = c.radius_sq.clone();
    let big_r2 = sphere.radius.square();
    if is_zero(&d2, dec)? {
        match compare(&big_r2, &r2, dec)? {
            Sign::Negative => return Ok(Plan::disjoint()),
            Sign::Zero => {
                return Ok(Plan::form(
                    circle(sphere.center.clone(), c.axis.clone(), r2, c.ref_dir.clone()),
                    Contact::Tangent,
                ));
            }
            Sign::Positive => {}
        }
        let k2 = (big_r2 - r2.clone()) / c.axis.norm_squared();
        let root = root_real(&k2);
        let v = c.axis.scaled(&root);
        let centers = [
            sphere.center.plus(&v),
            sphere.center.plus(&v.scaled(&s(-1))),
        ];
        let witness_offset = c
            .ref_dir
            .scaled(&root_real(&(&r2 / &c.ref_dir.norm_squared())));
        return Ok(Plan {
            curves: centers
                .into_iter()
                .zip([c.axis.scaled(&(-&sigma)), c.axis.scaled(&sigma)])
                .map(|(center, normal)| Spec {
                    witness: root
                        .rational_ref()
                        .is_none()
                        .then(|| center.plus(&witness_offset)),
                    form: circle(center, normal, r2.clone(), c.ref_dir.clone()),
                    contact: Contact::Crossing,
                    singular: vec![],
                })
                .collect(),
            divider: Some(Divider::new(sphere.center.clone(), c.axis.clone())),
            ..Plan::empty()
        });
    }
    let radius = c.radius();
    let sum = &radius + &sphere.radius;
    let diff = &radius - &sphere.radius;
    let positive_diff = sign(&diff, dec)? == Sign::Positive;
    let outer = compare(&d2, &sum.square(), dec)?;
    let inner = if positive_diff {
        compare(&d2, &diff.square(), dec)?
    } else {
        Sign::Positive
    };
    if outer == Sign::Positive || inner == Sign::Negative {
        return Ok(Plan::disjoint());
    }
    let d = root_real(&d2);
    let nearest = foot.plus(&offset.scaled(&(&radius / &d)));
    if outer == Sign::Zero || inner == Sign::Zero {
        return Ok(Plan::point(nearest));
    }
    let outer_diff = &sphere.radius - &radius;
    let mode = if sign(&outer_diff, dec)? == Sign::Positive {
        match compare(&outer_diff.square(), &d2, dec)? {
            Sign::Positive => 8,
            Sign::Zero => 7,
            Sign::Negative => 6,
        }
    } else {
        6
    };
    let z = ((&big_r2 - &(&radius - &d).square()) / c.axis.norm_squared()).sqrt();
    let upper = nearest.plus(&c.axis.scaled(&z));
    let lower = nearest.plus(&c.axis.scaled(&(-z)));
    let e3 = unit(&c.axis);
    let e1 = offset.scaled(&(&Real::one() / &d));
    let e2 = e3.cross(&e1);
    let k = big_r2 - r2 - d2;
    let frame = Frame::CylinderSphere {
        origin: foot.clone(),
        e1,
        e2,
        e3,
        radius: radius.clone(),
        d: d.clone(),
        k,
        mode,
        sense: if sign(&sigma, dec)? == Sign::Positive {
            1
        } else {
            -1
        },
    };
    let node = (mode == 7).then(|| foot.plus(&offset.scaled(&(&(-&radius) / &d))));
    let divider = (mode != 6).then(|| Divider::new(foot, c.axis.clone()));
    Ok(Plan::branches(
        if mode == 6 {
            vec![upper]
        } else {
            vec![upper, lower]
        },
        node,
        divider,
        frame,
    ))
}

fn plane_elliptic(
    p: &Plane,
    c: &EllipticCylinder,
    sigma: Real,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let e1 = unit(&c.major_dir);
    let e3 = unit(&c.axis);
    let e2 = e3.cross(&e1);
    let axial = p.normal.dot(&e3);
    if is_zero(&axial, dec)? {
        let n1 = p.normal.dot(&e1);
        let n2 = p.normal.dot(&e2);
        let h2 = &c.major_sq * &n1.square() + &c.minor_sq * &n2.square();
        let k = c.origin.minus(&p.origin).dot(&p.normal);
        let relation = compare(&k.square(), &h2, dec)?;
        if relation == Sign::Positive {
            return Ok(Plan::disjoint());
        }
        let shift = e1
            .scaled(&(&c.major_sq * &n1))
            .plus(&e2.scaled(&(&c.minor_sq * &n2)))
            .scaled(&(&(-&k) / &h2));
        let foot = c.origin.plus(&shift);
        if relation == Sign::Zero {
            return Ok(Plan::form(line(foot, c.axis.clone()), Contact::Tangent));
        }
        let direction = p.normal.cross(&c.axis);
        let a =
            &direction.dot(&e1).square() / &c.major_sq + &direction.dot(&e2).square() / &c.minor_sq;
        let height = (&Real::one() - &(&k.square() / &h2)) / a;
        return Ok(two_lines(
            foot,
            direction,
            height,
            [c.axis.scaled(&(-&sigma)), c.axis.scaled(&sigma)],
            true,
        ));
    }
    let center = c
        .origin
        .plus(&e3.scaled(&(&p.origin.minus(&c.origin).dot(&p.normal) / &axial)));
    let a = e1
        .minus(&e3.scaled(&(&p.normal.dot(&e1) / &axial)))
        .scaled(&root_real(&c.major_sq));
    let b = e2
        .minus(&e3.scaled(&(&p.normal.dot(&e2) / &axial)))
        .scaled(&root_real(&c.minor_sq));
    let aa = a.norm_squared();
    let bb = b.norm_squared();
    let ab = a.dot(&b);
    let gap = root_real(&((&aa - &bb).square() + s(4) * ab.square()));
    let major_sq = (&aa + &bb + gap.clone()) / s(2);
    let minor_sq = (&aa + &bb - gap.clone()) / s(2);
    let normal = p.normal.scaled(&sigma);
    if is_zero(&gap, dec)? {
        return Ok(Plan::form(
            circle(center, normal, major_sq, a),
            Contact::Crossing,
        ));
    }
    let major_dir = if is_zero(&ab, dec)? {
        if compare(&aa, &bb, dec)? == Sign::Positive {
            a
        } else {
            b
        }
    } else {
        a.scaled(&(&major_sq - &bb)).plus(&b.scaled(&ab))
    };
    Ok(Plan::form(
        CurveForm::Ellipse {
            center,
            normal,
            major_dir,
            major_sq,
            minor_sq,
        },
        Contact::Crossing,
    ))
}

struct CachedPlan {
    surfaces: [Surface; 2],
    plan: Arc<Plan>,
}
thread_local! {
    // Plans retain input Real nodes. Value-only keys would attach a fresh query
    // to an earlier query's DAG and change memoized proof work. Reuse only the
    // same input nodes; algebraic decision work is never cached here.
    static PLANS: RefCell<VecDeque<CachedPlan>> = const { RefCell::new(VecDeque::new()) };
}
pub(super) fn plan(
    geom: &GeomStore,
    ids: [SurfaceId; 2],
    dec: &mut Decider,
) -> Result<Arc<Plan>, Refusal> {
    let a = geom.surface(ids[0]);
    let b = geom.surface(ids[1]);
    let known = |kind| {
        geom.facts()
            .holds(kind, GeomKey::Surface(ids[0]), GeomKey::Surface(ids[1]))
    };
    if known(FactKind::Coincident) {
        dec.note_fact_decision();
        let same_sense = match (a, b) {
            (Surface::Plane(a), Surface::Plane(b)) => {
                sign(&a.normal.dot(&b.normal), dec)? == Sign::Positive
            }
            _ => true,
        };
        return Ok(Arc::new(Plan::coincident(same_sense)));
    }
    let radial_measures = match (a, b) {
        (Surface::Sphere(a), Surface::Sphere(b)) if known(FactKind::Concentric) => {
            Some((&a.radius, &b.radius))
        }
        (Surface::Cylinder(a), Surface::Cylinder(b)) if known(FactKind::Coaxial) => {
            Some((&a.radius_sq, &b.radius_sq))
        }
        _ => None,
    };
    if let Some((a, b)) = radial_measures {
        dec.note_fact_decision();
        return Ok(Arc::new(if is_zero(&(a - b), dec)? {
            Plan::coincident(true)
        } else {
            Plan::disjoint()
        }));
    }
    let rational = a.has_rational_definition() && b.has_rational_definition();
    if rational
        && let Some(plan) = PLANS.with(|plans| {
            plans
                .borrow()
                .iter()
                .find(|entry| {
                    entry.surfaces[0].shares_definition_nodes(a)
                        && entry.surfaces[1].shares_definition_nodes(b)
                })
                .map(|entry| entry.plan.clone())
        })
    {
        return Ok(plan);
    }
    let plan = Arc::new(uncached_plan(a, b, rational, dec)?);
    if rational {
        PLANS.with(|plans| {
            let mut plans = plans.borrow_mut();
            if plans.len() == 32 {
                plans.pop_front();
            }
            plans.push_back(CachedPlan {
                surfaces: [a.clone(), b.clone()],
                plan: plan.clone(),
            });
        });
    }
    Ok(plan)
}
fn uncached_plan(
    a: &Surface,
    b: &Surface,
    rational: bool,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    match (a, b) {
        (Surface::Plane(a), Surface::Plane(b)) => planes(a, b, dec),
        (Surface::Plane(p), Surface::Cylinder(c)) => plane_cylinder(p, c, s(1), dec),
        (Surface::Cylinder(c), Surface::Plane(p)) => plane_cylinder(p, c, s(-1), dec),
        (Surface::Plane(p), Surface::Sphere(c)) => plane_sphere(p, c, s(1), dec),
        (Surface::Sphere(c), Surface::Plane(p)) => plane_sphere(p, c, s(-1), dec),
        (Surface::Sphere(a), Surface::Sphere(b)) => spheres(a, b, dec),
        (Surface::Cylinder(a), Surface::Cylinder(b)) => cylinders(a, b, rational, dec),
        (Surface::Cylinder(c), Surface::Sphere(sphere)) => cylinder_sphere(c, sphere, s(1), dec),
        (Surface::Sphere(sphere), Surface::Cylinder(c)) => cylinder_sphere(c, sphere, s(-1), dec),
        (Surface::Plane(p), Surface::EllipticCylinder(c)) => plane_elliptic(p, c, s(1), dec),
        (Surface::EllipticCylinder(c), Surface::Plane(p)) => plane_elliptic(p, c, s(-1), dec),
        (Surface::Plane(p), Surface::Cone(c)) => conics::plane_cone(p, c, s(1), rational, dec),
        (Surface::Cone(c), Surface::Plane(p)) => conics::plane_cone(p, c, s(-1), rational, dec),
        (Surface::Cone(c), Surface::Cylinder(y)) => conics::cone_cylinder(c, y, s(1), dec),
        (Surface::Cylinder(y), Surface::Cone(c)) => conics::cone_cylinder(c, y, s(-1), dec),
        (Surface::Cone(c), Surface::Sphere(sphere)) => conics::cone_sphere(c, sphere, s(1), dec),
        (Surface::Sphere(sphere), Surface::Cone(c)) => conics::cone_sphere(c, sphere, s(-1), dec),
        (Surface::Cone(a), Surface::Cone(b)) => conics::cones(a, b, dec),
        (Surface::Plane(p), Surface::Torus(t)) => conics::plane_torus(p, t, s(1), dec),
        (Surface::Torus(t), Surface::Plane(p)) => conics::plane_torus(p, t, s(-1), dec),
        (Surface::Cylinder(c), Surface::Torus(t)) => conics::cylinder_torus(c, t, s(1), dec),
        (Surface::Torus(t), Surface::Cylinder(c)) => conics::cylinder_torus(c, t, s(-1), dec),
        (Surface::Sphere(sphere), Surface::Torus(t)) => conics::sphere_torus(sphere, t, s(1), dec),
        (Surface::Torus(t), Surface::Sphere(sphere)) => conics::sphere_torus(sphere, t, s(-1), dec),
        (Surface::Cone(c), Surface::Torus(t)) => conics::cone_torus(c, t, s(1), rational, dec),
        (Surface::Torus(t), Surface::Cone(c)) => conics::cone_torus(c, t, s(-1), rational, dec),
        (Surface::Torus(a), Surface::Torus(b)) => conics::tori(a, b, dec),
        _ => Err(refused_pair()),
    }
}
fn record(geom: &mut GeomStore, kind: FactKind, a: GeomKey, b: GeomKey, at: Option<GeomKey>) {
    geom.record_fact(Fact {
        kind,
        a,
        b,
        at,
        origin: FactOrigin::Construction("intersect.surfaces"),
    });
}
fn add_point(
    geom: &mut GeomStore,
    surfaces: [SurfaceId; 2],
    p: Point3,
) -> catallax_shapes::PointId {
    geom.add_point(PointDef::Meet(Meet {
        of: MeetOf::Surfaces(surfaces.to_vec()),
        position: Position::Closed(p),
    }))
}
pub(super) fn construct(
    geom: &mut GeomStore,
    a: SurfaceId,
    b: SurfaceId,
    dec: &mut Decider,
) -> Answered<SurfacePair> {
    let plan = plan(geom, [a, b], dec)?;
    if let Some(answer) = &plan.terminal {
        if matches!(answer, SurfacePair::Coincident { .. }) {
            record(
                geom,
                FactKind::Coincident,
                GeomKey::Surface(a),
                GeomKey::Surface(b),
                None,
            );
        }
        return Ok(Answer::exact(answer.clone()));
    }
    let points: Vec<_> = plan
        .points
        .iter()
        .cloned()
        .map(|p| add_point(geom, [a, b], p))
        .collect();
    let mut components = vec![];
    if plan.curves.is_empty() {
        for point in points {
            record(
                geom,
                FactKind::Tangent,
                GeomKey::Surface(a),
                GeomKey::Surface(b),
                Some(GeomKey::Point(point)),
            );
            components.push(Component::Point { point });
        }
    } else {
        for (index, &point) in points.iter().enumerate() {
            if !plan
                .curves
                .iter()
                .any(|spec| spec.singular.contains(&index))
            {
                record(
                    geom,
                    FactKind::Tangent,
                    GeomKey::Surface(a),
                    GeomKey::Surface(b),
                    Some(GeomKey::Point(point)),
                );
                components.push(Component::Point { point });
            }
        }
        for spec in &plan.curves {
            let witness = spec
                .witness
                .as_ref()
                .map(|p| add_point(geom, [a, b], p.clone()));
            let definition = if let Some(witness) = witness {
                Curve::Intersection(IntersectionCurve {
                    surfaces: [a, b],
                    witness,
                })
            } else {
                match &spec.form {
                    CurveForm::Line { origin, dir } => Curve::Line(Line {
                        origin: origin.clone(),
                        dir: dir.clone(),
                    }),
                    CurveForm::Circle {
                        center,
                        normal,
                        radius_sq,
                        ref_dir,
                    } => Curve::Circle(Circle {
                        center: center.clone(),
                        normal: normal.clone(),
                        radius_sq: radius_sq.clone(),
                        ref_dir: ref_dir.clone(),
                    }),
                    CurveForm::Ellipse {
                        center,
                        normal,
                        major_dir,
                        major_sq,
                        minor_sq,
                    } => Curve::Ellipse(Ellipse {
                        center: center.clone(),
                        normal: normal.clone(),
                        major_dir: major_dir.clone(),
                        major_sq: major_sq.clone(),
                        minor_sq: minor_sq.clone(),
                    }),
                    CurveForm::BSpline(curve) => Curve::BSpline(curve.clone()),
                    CurveForm::Branch
                    | CurveForm::Parabola { .. }
                    | CurveForm::Hyperbola { .. } => {
                        unreachable!("non-analytic form has a witness")
                    }
                }
            };
            let analytic = !matches!(definition, Curve::Intersection(_));
            let curve = geom.add_curve(definition);
            if analytic {
                for surface in [a, b] {
                    record(
                        geom,
                        FactKind::Incident,
                        GeomKey::Curve(curve),
                        GeomKey::Surface(surface),
                        None,
                    );
                }
            }
            // A circular component remains a circular section when its exact
            // definition needs an intersection witness. GeomStore can infer
            // this from Curve::Circle incidence, but cannot inspect this plan.
            if matches!(spec.form, CurveForm::Circle { .. }) {
                for surface in [a, b] {
                    if matches!(
                        geom.surface(surface),
                        Surface::Cylinder(_) | Surface::Cone(_)
                    ) {
                        record(
                            geom,
                            FactKind::Coaxial,
                            GeomKey::Curve(curve),
                            GeomKey::Surface(surface),
                            None,
                        );
                    }
                }
            }
            if let Some(witness) = witness {
                record(
                    geom,
                    FactKind::Incident,
                    GeomKey::Point(witness),
                    GeomKey::Curve(curve),
                    None,
                );
            }
            let singular: Vec<_> = spec.singular.iter().map(|&i| points[i]).collect();
            for point in &singular {
                record(
                    geom,
                    FactKind::Incident,
                    GeomKey::Point(*point),
                    GeomKey::Curve(curve),
                    None,
                );
            }
            if spec.contact == Contact::Tangent {
                record(
                    geom,
                    FactKind::Tangent,
                    GeomKey::Surface(a),
                    GeomKey::Surface(b),
                    Some(GeomKey::Curve(curve)),
                );
            }
            components.push(Component::Curve {
                curve,
                contact: spec.contact,
                singular,
            });
        }
    }
    Ok(Answer::exact(SurfacePair::Meets(components)))
}
