use crate::math::*;
use crate::{Component, Contact, CurveForm, SurfacePair};
use catallax_numbers::{
    Answer, Answered, Fact, FactKind, FactOrigin, Rational, Real, Refusal, Sign,
};
use catallax_shapes::{
    Circle, Curve, Cylinder, Ellipse, GeomKey, GeomStore, IntersectionCurve, Line, Meet, MeetOf,
    Plane, Point3, PointDef, Position, RPoint3, RVector3, Sphere, Surface, SurfaceId, Vector3,
};
use std::{cell::RefCell, collections::VecDeque, sync::Arc};

#[derive(Clone, Debug)]
pub(super) struct Divider {
    pub origin: Point3,
    pub vector: Vector3,
    pub length_lower: Rational,
}
impl Divider {
    fn rational(origin: RPoint3, vector: RVector3) -> Self {
        let lower = vector.x.abs().max(vector.y.abs()).max(vector.z.abs());
        Self {
            origin: origin.to_real(),
            vector: vector.to_real(),
            length_lower: lower,
        }
    }
    pub fn value(&self, p: &Point3) -> Real {
        p.minus(&self.origin).dot(&self.vector)
    }
}
#[derive(Clone, Debug)]
pub(super) enum Frame {
    CylinderSphere {
        origin: RPoint3,
        e1: Vector3,
        e2: Vector3,
        e3: Vector3,
        radius: Rational,
        d: Real,
        k: Real,
        mode: u8,
        sense: i64,
    },
    Cylinders {
        origin: RPoint3,
        e1: Vector3,
        e2: Vector3,
        e3: Vector3,
        big: Rational,
        small: Rational,
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
fn line(origin: Point3, dir: RVector3) -> CurveForm {
    CurveForm::Line { origin, dir }
}
fn circle(center: Point3, normal: RVector3, radius_sq: Rational, ref_dir: RVector3) -> CurveForm {
    CurveForm::Circle {
        center,
        normal,
        radius_sq,
        ref_dir,
    }
}
fn two_lines(
    base: RPoint3,
    v: RVector3,
    h2: Rational,
    dirs: [RVector3; 2],
    first_positive: bool,
) -> Plan {
    let offset = v.to_real().scaled(&root(&h2));
    let origins = [
        base.to_real().plus(&offset),
        base.to_real().plus(&offset.scaled(&Real::from_i64(-1))),
    ];
    let witness = rational_sqrt(&h2).is_none();
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
        divider: Some(Divider::rational(
            base,
            if first_positive { v } else { v.scaled(&r(-1)) },
        )),
        ..Plan::empty()
    }
}
fn planes(a: &Plane, b: &Plane) -> Plan {
    let m = a.normal.cross(&b.normal);
    if m.is_zero() {
        return if b.origin.minus(&a.origin).dot(&a.normal).is_zero() {
            Plan::coincident(a.normal.dot(&b.normal).sign() == Sign::Positive)
        } else {
            Plan::disjoint()
        };
    }
    let scale = b.origin.minus(&a.origin).dot(&b.normal) / m.norm_squared();
    let origin = a.origin.plus(&a.normal.cross(&m).scaled(&(-scale)));
    Plan::form(line(origin.to_real(), m), Contact::Crossing)
}
fn plane_cylinder(p: &Plane, c: &Cylinder, sigma: Rational) -> Plan {
    let s = p.normal.dot(&c.axis);
    let r2 = sq(&c.radius);
    if !s.is_zero() {
        let center = c.origin.plus(
            &c.axis
                .scaled(&(&p.origin.minus(&c.origin).dot(&p.normal) / &s)),
        );
        let normal = p.normal.scaled(&sigma);
        if c.axis.cross(&p.normal).is_zero() {
            return Plan::form(
                circle(center.to_real(), normal, r2, c.ref_dir.clone()),
                Contact::Crossing,
            );
        }
        return Plan::form(
            CurveForm::Ellipse {
                center: center.to_real(),
                normal,
                major_dir: c
                    .axis
                    .minus(&p.normal.scaled(&(&s / &p.normal.norm_squared()))),
                major_sq: r2.clone() * p.normal.norm_squared() * c.axis.norm_squared() / sq(&s),
                minor_sq: r2,
            },
            Contact::Crossing,
        );
    }
    let k = c.origin.minus(&p.origin).dot(&p.normal);
    let e2 = sq(&k) / p.normal.norm_squared();
    if e2 > r2 {
        return Plan::disjoint();
    }
    let foot = c
        .origin
        .plus(&p.normal.scaled(&(-k / p.normal.norm_squared())));
    if e2 == r2 {
        return Plan::form(line(foot.to_real(), c.axis.clone()), Contact::Tangent);
    }
    let v = p.normal.cross(&c.axis);
    let h2 = (r2 - e2) / v.norm_squared();
    two_lines(
        foot,
        v,
        h2,
        [c.axis.scaled(&(-&sigma)), c.axis.scaled(&sigma)],
        true,
    )
}
fn plane_sphere(p: &Plane, s: &Sphere, sigma: Rational) -> Plan {
    let k = s.center.minus(&p.origin).dot(&p.normal);
    let e2 = sq(&k) / p.normal.norm_squared();
    let r2 = sq(&s.radius);
    if e2 > r2 {
        return Plan::disjoint();
    }
    let foot = s
        .center
        .plus(&p.normal.scaled(&(-k / p.normal.norm_squared())));
    if e2 == r2 {
        return Plan::point(foot.to_real());
    }
    Plan::form(
        circle(
            foot.to_real(),
            p.normal.scaled(&sigma),
            r2 - e2,
            p.ref_dir.clone(),
        ),
        Contact::Crossing,
    )
}
fn spheres(a: &Sphere, b: &Sphere) -> Plan {
    let u = b.center.minus(&a.center);
    let d2 = u.norm_squared();
    if d2.is_zero() {
        return if a.radius == b.radius {
            Plan::coincident(true)
        } else {
            Plan::disjoint()
        };
    }
    let sum = &a.radius + &b.radius;
    let diff = &a.radius - &b.radius;
    if d2 > sq(&sum) || d2 < sq(&diff) {
        return Plan::disjoint();
    }
    if d2 == sq(&sum) {
        return Plan::point(a.center.plus(&u.scaled(&(&a.radius / &sum))).to_real());
    }
    if d2 == sq(&diff) {
        return Plan::point(a.center.plus(&u.scaled(&(&a.radius / &diff))).to_real());
    }
    let t = (d2.clone() + sq(&a.radius) - sq(&b.radius)) / (r(2) * d2.clone());
    Plan::form(
        circle(
            a.center.plus(&u.scaled(&t)).to_real(),
            u.clone(),
            sq(&a.radius) - sq(&t) * d2,
            perp(&u),
        ),
        Contact::Crossing,
    )
}
fn parallel_cylinders(a: &Cylinder, b: &Cylinder) -> Plan {
    let d = b.origin.minus(&a.origin);
    let u = d.minus(&a.axis.scaled(&(d.dot(&a.axis) / a.axis.norm_squared())));
    let d2 = u.norm_squared();
    if d2.is_zero() {
        return if a.radius == b.radius {
            Plan::coincident(true)
        } else {
            Plan::disjoint()
        };
    }
    let sum = &a.radius + &b.radius;
    let diff = &a.radius - &b.radius;
    if d2 > sq(&sum) || d2 < sq(&diff) {
        return Plan::disjoint();
    }
    if d2 == sq(&sum) {
        return Plan::form(
            line(
                a.origin.plus(&u.scaled(&(&a.radius / &sum))).to_real(),
                a.axis.clone(),
            ),
            Contact::Tangent,
        );
    }
    if d2 == sq(&diff) {
        return Plan::form(
            line(
                a.origin.plus(&u.scaled(&(&a.radius / &diff))).to_real(),
                a.axis.clone(),
            ),
            Contact::Tangent,
        );
    }
    let t = (d2.clone() + sq(&a.radius) - sq(&b.radius)) / (r(2) * d2.clone());
    let v = a.axis.cross(&u);
    let h2 = (sq(&a.radius) - sq(&t) * d2) / v.norm_squared();
    two_lines(
        a.origin.plus(&u.scaled(&t)),
        v,
        h2,
        [a.axis.clone(), a.axis.scaled(&r(-1))],
        true,
    )
}
fn cylinders(a: &Cylinder, b: &Cylinder) -> Result<Plan, Refusal> {
    if a.axis.cross(&b.axis).is_zero() {
        return Ok(parallel_cylinders(a, b));
    }
    let (big, small) = if a.radius >= b.radius { (a, b) } else { (b, a) };
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
    let delta2 = sq(&d.dot(&m)) / m.norm_squared();
    let sum = &big.radius + &small.radius;
    let diff = &big.radius - &small.radius;
    if delta2 > sq(&sum) {
        return Ok(Plan::disjoint());
    }
    if delta2 == sq(&sum) {
        return Ok(Plan::point(
            fb.plus(&fs.minus(&fb).scaled(&(&big.radius / &sum)))
                .to_real(),
        ));
    }
    if delta2.is_zero() && diff.is_zero() {
        let Some(k) = rational_sqrt(&(big.axis.norm_squared() * small.axis.norm_squared())) else {
            return Err(refused_pair());
        };
        let nm = big
            .axis
            .scaled(&small.axis.norm_squared())
            .minus(&small.axis.scaled(&k));
        let np = big
            .axis
            .scaled(&small.axis.norm_squared())
            .plus(&small.axis.scaled(&k));
        let cos = big.axis.dot(&small.axis) / k;
        let offset = m
            .to_real()
            .scaled(&root(&(sq(&big.radius) / m.norm_squared())));
        let points = vec![
            fb.to_real().plus(&offset),
            fb.to_real().plus(&offset.scaled(&Real::from_i64(-1))),
        ];
        return Ok(Plan {
            points,
            curves: [
                (nm.clone(), np.clone(), r(1) - cos.clone()),
                (np, nm, r(1) + cos),
            ]
            .into_iter()
            .map(|(normal, major_dir, denom)| Spec {
                form: CurveForm::Ellipse {
                    center: fb.to_real(),
                    normal,
                    major_dir,
                    major_sq: r(2) * sq(&big.radius) / denom,
                    minor_sq: sq(&big.radius),
                },
                witness: None,
                contact: Contact::Crossing,
                singular: vec![0, 1],
            })
            .collect(),
            ..Plan::empty()
        });
    }
    let delta = root(&delta2);
    let radial = if delta2.is_zero() { m } else { fs.minus(&fb) };
    let transverse = radial.cross(&big.axis);
    let e3 = unit(&big.axis);
    let e2 = unit(&radial);
    let e1 = unit(&transverse);
    // Cancel dot-product terms in the rational field before normalising.
    let alpha_numerator = small.axis.dot(&transverse);
    let alpha =
        real(&alpha_numerator) / root(&(small.axis.norm_squared() * transverse.norm_squared()));
    let sigma = if a.radius >= b.radius { 1 } else { -1 };
    let sense = -sigma
        * if alpha_numerator.sign() == Sign::Positive {
            1
        } else {
            -1
        };
    let gamma = real(&small.axis.dot(&big.axis))
        / root(&(small.axis.norm_squared() * big.axis.norm_squared()));
    let mode = if delta2 > sq(&diff) {
        3
    } else if delta2 == sq(&diff) {
        4
    } else {
        5
    };
    let frame = Frame::Cylinders {
        origin: fb.clone(),
        e1: e1.clone(),
        e2: e2.clone(),
        e3: e3.clone(),
        big: big.radius.clone(),
        small: small.radius.clone(),
        delta: delta.clone(),
        alpha: alpha.clone(),
        gamma: gamma.clone(),
        mode,
        sense,
    };
    if mode == 3 {
        let z = (real(&sq(&small.radius)) - (real(&big.radius) - delta).square()).sqrt() / alpha;
        let witness = fb
            .to_real()
            .plus(&e2.scaled(&real(&big.radius)))
            .plus(&e3.scaled(&z));
        return Ok(Plan::branches(vec![witness], None, None, frame));
    }
    let x = root(&(sq(&big.radius) - delta2));
    let base = fs.to_real();
    let wplus = base
        .plus(&e1.scaled(&x))
        .plus(&e3.scaled(&((&gamma * &x + real(&small.radius)) / alpha.clone())));
    let wminus = base
        .plus(&e1.scaled(&(-&x)))
        .plus(&e3.scaled(&((real(&small.radius) - &gamma * &x) / alpha)));
    let node = (mode == 4).then(|| {
        fb.plus(&fs.minus(&fb).scaled(&(&big.radius / &diff)))
            .to_real()
    });
    Ok(Plan::branches(
        vec![wplus, wminus],
        node,
        Some(Divider {
            origin: fb.to_real(),
            vector: e1,
            length_lower: r(1),
        }),
        frame,
    ))
}
fn cylinder_sphere(c: &Cylinder, s: &Sphere, sigma: Rational) -> Plan {
    let foot = c.origin.plus(
        &c.axis
            .scaled(&(s.center.minus(&c.origin).dot(&c.axis) / c.axis.norm_squared())),
    );
    let offset = s.center.minus(&foot);
    let d2 = offset.norm_squared();
    let r2 = sq(&c.radius);
    let big_r2 = sq(&s.radius);
    if d2.is_zero() {
        if s.radius < c.radius {
            return Plan::disjoint();
        }
        if s.radius == c.radius {
            return Plan::form(
                circle(s.center.to_real(), c.axis.clone(), r2, c.ref_dir.clone()),
                Contact::Tangent,
            );
        }
        let k2 = (big_r2 - r2.clone()) / c.axis.norm_squared();
        let v = c.axis.to_real().scaled(&root(&k2));
        let centers = [
            s.center.to_real().plus(&v),
            s.center.to_real().plus(&v.scaled(&Real::from_i64(-1))),
        ];
        let witness_offset = c
            .ref_dir
            .to_real()
            .scaled(&root(&(&r2 / &c.ref_dir.norm_squared())));
        return Plan {
            curves: centers
                .into_iter()
                .zip([c.axis.scaled(&(-&sigma)), c.axis.scaled(&sigma)])
                .map(|(center, normal)| Spec {
                    witness: rational_sqrt(&k2)
                        .is_none()
                        .then(|| center.plus(&witness_offset)),
                    form: circle(center, normal, r2.clone(), c.ref_dir.clone()),
                    contact: Contact::Crossing,
                    singular: vec![],
                })
                .collect(),
            divider: Some(Divider::rational(s.center.clone(), c.axis.clone())),
            ..Plan::empty()
        };
    }
    // Comparing d with r +/- R avoids comparing irrational distances.
    let sum = &c.radius + &s.radius;
    let diff = &c.radius - &s.radius;
    if d2 > sq(&sum) || (diff.sign() == Sign::Positive && d2 < sq(&diff)) {
        return Plan::disjoint();
    }
    let d = root(&d2);
    let nearest = foot
        .to_real()
        .plus(&offset.to_real().scaled(&(&real(&c.radius) / &d)));
    if d2 == sq(&sum) || (diff.sign() == Sign::Positive && d2 == sq(&diff)) {
        return Plan::point(nearest);
    }
    let outer_diff = &s.radius - &c.radius;
    let mode = if outer_diff.sign() == Sign::Positive && sq(&outer_diff) > d2 {
        8
    } else if outer_diff.sign() == Sign::Positive && sq(&outer_diff) == d2 {
        7
    } else {
        6
    };
    let z =
        ((real(&big_r2) - (&real(&c.radius) - &d).square()) / real(&c.axis.norm_squared())).sqrt();
    let upper = nearest.plus(&c.axis.to_real().scaled(&z));
    let lower = nearest.plus(&c.axis.to_real().scaled(&(-z)));
    let e3 = unit(&c.axis);
    let e1 = offset.to_real().scaled(&(&Real::one() / &d));
    let e2 = e3.cross(&e1);
    let k = real(&(big_r2 - r2 - d2.clone()));
    let frame = Frame::CylinderSphere {
        origin: foot.clone(),
        e1,
        e2,
        e3,
        radius: c.radius.clone(),
        d: d.clone(),
        k,
        mode,
        sense: if sigma.sign() == Sign::Positive {
            1
        } else {
            -1
        },
    };
    let node = (mode == 7).then(|| {
        foot.to_real()
            .plus(&offset.to_real().scaled(&(&-real(&c.radius) / &d)))
    });
    let divider = (mode != 6).then(|| Divider::rational(foot, c.axis.clone()));
    Plan::branches(
        if mode == 6 {
            vec![upper]
        } else {
            vec![upper, lower]
        },
        node,
        divider,
        frame,
    )
}
struct CachedPlan {
    surfaces: [Surface; 2],
    plan: Arc<Plan>,
}

thread_local! {
    // Keys are complete ordered definitions, never arena handles or addresses.
    // The cache cannot retain a whole model and has a fixed per-thread cap.
    static PLANS: RefCell<VecDeque<CachedPlan>> = const { RefCell::new(VecDeque::new()) };
}

pub(super) fn plan(a: &Surface, b: &Surface) -> Result<Arc<Plan>, Refusal> {
    if let Some(plan) = PLANS.with(|plans| {
        plans
            .borrow()
            .iter()
            .find(|entry| entry.surfaces[0] == *a && entry.surfaces[1] == *b)
            .map(|entry| entry.plan.clone())
    }) {
        return Ok(plan);
    }
    let plan = Arc::new(uncached_plan(a, b)?);
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
    Ok(plan)
}

fn uncached_plan(a: &Surface, b: &Surface) -> Result<Plan, Refusal> {
    Ok(match (a, b) {
        (Surface::Plane(a), Surface::Plane(b)) => planes(a, b),
        (Surface::Plane(p), Surface::Cylinder(c)) => plane_cylinder(p, c, r(1)),
        (Surface::Cylinder(c), Surface::Plane(p)) => plane_cylinder(p, c, r(-1)),
        (Surface::Plane(p), Surface::Sphere(s)) => plane_sphere(p, s, r(1)),
        (Surface::Sphere(s), Surface::Plane(p)) => plane_sphere(p, s, r(-1)),
        (Surface::Sphere(a), Surface::Sphere(b)) => spheres(a, b),
        (Surface::Cylinder(a), Surface::Cylinder(b)) => return cylinders(a, b),
        (Surface::Cylinder(c), Surface::Sphere(s)) => cylinder_sphere(c, s, r(1)),
        (Surface::Sphere(s), Surface::Cylinder(c)) => cylinder_sphere(c, s, r(-1)),
        _ => return Err(refused_pair()),
    })
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
pub(super) fn construct(geom: &mut GeomStore, a: SurfaceId, b: SurfaceId) -> Answered<SurfacePair> {
    let plan = plan(geom.surface(a), geom.surface(b))?;
    if let Some(answer) = &plan.terminal {
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
                        origin: point_rational(origin).expect("analytic rational origin"),
                        dir: dir.clone(),
                    }),
                    CurveForm::Circle {
                        center,
                        normal,
                        radius_sq,
                        ref_dir,
                    } => Curve::Circle(Circle {
                        center: point_rational(center).expect("analytic rational center"),
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
                        center: point_rational(center).expect("analytic rational center"),
                        normal: normal.clone(),
                        major_dir: major_dir.clone(),
                        major_sq: major_sq.clone(),
                        minor_sq: minor_sq.clone(),
                    }),
                    CurveForm::Branch => unreachable!("a branch has a witness"),
                }
            };
            let analytic = !matches!(definition, Curve::Intersection(_));
            let curve = geom.add_curve(definition);
            if analytic {
                for s in [a, b] {
                    record(
                        geom,
                        FactKind::Incident,
                        GeomKey::Curve(curve),
                        GeomKey::Surface(s),
                        None,
                    );
                }
            }
            if let Some(w) = witness {
                record(
                    geom,
                    FactKind::Incident,
                    GeomKey::Point(w),
                    GeomKey::Curve(curve),
                    None,
                );
            }
            let singular: Vec<_> = spec.singular.iter().map(|&i| points[i]).collect();
            for p in &singular {
                record(
                    geom,
                    FactKind::Incident,
                    GeomKey::Point(*p),
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
