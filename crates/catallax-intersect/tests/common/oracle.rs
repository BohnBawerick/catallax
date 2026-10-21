//! Expected answers are the contract's algebra, never calls into catallax-intersect.
use super::*;

#[derive(Clone, Debug)]
pub enum Expected {
    Disjoint,
    Coincident(bool),
    Refused,
    Meets(Vec<Part>),
}
#[derive(Clone, Debug)]
pub enum Shape {
    Analytic(Box<Curve>),
    Intersection {
        witness: Point3,
        form: Box<CurveForm>,
    },
}
#[derive(Clone, Debug)]
pub enum Part {
    Point(Point3),
    Curve {
        shape: Shape,
        contact: Contact,
        singular: Vec<Point3>,
    },
}
fn cv(shape: Shape, contact: Contact) -> Part {
    Part::Curve {
        shape,
        contact,
        singular: vec![],
    }
}
fn line(o: Point3, d: RVector3, contact: Contact) -> Part {
    cv(
        match rational_point(&o) {
            Some(origin) => Shape::Analytic(Box::new(Curve::Line(Line::new(origin, d).unwrap()))),
            None => Shape::Intersection {
                witness: o.clone(),
                form: Box::new(CurveForm::Line { origin: o, dir: d }),
            },
        },
        contact,
    )
}
fn circle(
    o: Point3,
    normal: RVector3,
    radius_sq: Rational,
    reference: RVector3,
    contact: Contact,
) -> Part {
    let shape = match rational_point(&o) {
        Some(center) => Shape::Analytic(Box::new(Curve::Circle(
            Circle::new(center, normal, radius_sq, reference).unwrap(),
        ))),
        None => {
            let witness = o.plus(
                &reference
                    .to_real()
                    .scaled(&root(&(&radius_sq / &reference.norm_squared()))),
            );
            Shape::Intersection {
                witness,
                form: Box::new(CurveForm::Circle {
                    center: o,
                    normal,
                    radius_sq: Real::from(radius_sq),
                    ref_dir: reference,
                }),
            }
        }
    };
    cv(shape, contact)
}
fn ellipse(o: RPoint3, normal: RVector3, major: RVector3, a: Rational, b: Rational) -> Part {
    cv(
        Shape::Analytic(Box::new(Curve::Ellipse(
            Ellipse::new(o, normal, major, a, b).unwrap(),
        ))),
        Contact::Crossing,
    )
}
fn branch(witness: Point3, singular: Vec<Point3>) -> Part {
    Part::Curve {
        shape: Shape::Intersection {
            witness,
            form: Box::new(CurveForm::Branch),
        },
        contact: Contact::Crossing,
        singular,
    }
}
fn meets(p: Part) -> Expected {
    Expected::Meets(vec![p])
}
fn point(p: Point3) -> Expected {
    meets(Part::Point(p))
}
fn sq(x: &Rational) -> Rational {
    x * x
}
fn neg(v: &RVector3) -> RVector3 {
    v.scaled(&n(-1))
}

pub fn surfaces(a: &Surface, b: &Surface) -> Expected {
    use Contact::*;
    match (a, b) {
        (Surface::Plane(a), Surface::Plane(b)) => {
            let m = a.normal.cross(&b.normal);
            let k = b.origin.minus(&a.origin).dot(&b.normal);
            if m.is_zero() {
                if a.origin.minus(&b.origin).dot(&a.normal).is_zero() {
                    Expected::Coincident(a.normal.dot(&b.normal) > n(0))
                } else {
                    Expected::Disjoint
                }
            } else {
                let o = a
                    .origin
                    .plus(&a.normal.cross(&m).scaled(&(-k / m.norm_squared())));
                meets(line(o.to_real(), m, Crossing))
            }
        }
        (Surface::Plane(p), Surface::Cylinder(c)) => pc(p, c, 1),
        (Surface::Cylinder(c), Surface::Plane(p)) => pc(p, c, -1),
        (Surface::Plane(p), Surface::Sphere(s)) => ps(p, s, 1),
        (Surface::Sphere(s), Surface::Plane(p)) => ps(p, s, -1),
        (Surface::Sphere(a), Surface::Sphere(b)) => {
            let u = b.center.minus(&a.center);
            let d2 = u.norm_squared();
            let ra = &a.radius;
            let rb = &b.radius;
            if d2.is_zero() {
                return if ra == rb {
                    Expected::Coincident(true)
                } else {
                    Expected::Disjoint
                };
            }
            if d2 > sq(&(ra + rb)) || d2 < sq(&(ra - rb)) {
                return Expected::Disjoint;
            }
            if d2 == sq(&(ra + rb)) {
                return point(a.center.plus(&u.scaled(&(ra / &(ra + rb)))).to_real());
            }
            if d2 == sq(&(ra - rb)) {
                return point(a.center.plus(&u.scaled(&(ra / &(ra - rb)))).to_real());
            }
            let t = (&d2 + &sq(ra) - sq(rb)) / (&n(2) * &d2);
            meets(circle(
                a.center.plus(&u.scaled(&t)).to_real(),
                u.clone(),
                sq(ra) - sq(&t) * d2,
                perp(&u),
                Crossing,
            ))
        }
        (Surface::Cylinder(a), Surface::Cylinder(b)) => cc(a, b),
        (Surface::Cylinder(c), Surface::Sphere(s)) => cs(c, s, 1),
        (Surface::Sphere(s), Surface::Cylinder(c)) => cs(c, s, -1),
        _ => Expected::Refused,
    }
}
fn pc(p: &Plane, c: &Cylinder, sigma: i64) -> Expected {
    use Contact::*;
    let s = p.normal.dot(&c.axis);
    let normal = p.normal.scaled(&n(sigma));
    if !s.is_zero() {
        let o = c.origin.plus(
            &c.axis
                .scaled(&(&p.origin.minus(&c.origin).dot(&p.normal) / &s)),
        );
        if c.axis.cross(&p.normal).is_zero() {
            return meets(circle(
                o.to_real(),
                normal,
                sq(&c.radius),
                c.ref_dir.clone(),
                Crossing,
            ));
        }
        let major = c
            .axis
            .minus(&p.normal.scaled(&(&s / &p.normal.norm_squared())));
        let a = sq(&c.radius) * p.normal.norm_squared() * c.axis.norm_squared() / sq(&s);
        return meets(ellipse(o, normal, major, a, sq(&c.radius)));
    }
    let k = c.origin.minus(&p.origin).dot(&p.normal);
    let f = c
        .origin
        .plus(&p.normal.scaled(&(-&k / p.normal.norm_squared())));
    let e2 = sq(&k) / p.normal.norm_squared();
    let r2 = sq(&c.radius);
    if e2 > r2 {
        return Expected::Disjoint;
    }
    if e2 == r2 {
        return meets(line(f.to_real(), c.axis.clone(), Tangent));
    }
    let v = p.normal.cross(&c.axis);
    let h = root(&((r2 - e2) / v.norm_squared()));
    Expected::Meets(vec![
        line(plus(&f, &v, &h), c.axis.scaled(&n(-sigma)), Crossing),
        line(plus(&f, &v, &(-h)), c.axis.scaled(&n(sigma)), Crossing),
    ])
}
fn ps(p: &Plane, s: &Sphere, sigma: i64) -> Expected {
    let k = s.center.minus(&p.origin).dot(&p.normal);
    let e2 = sq(&k) / p.normal.norm_squared();
    let r2 = sq(&s.radius);
    let f = s
        .center
        .plus(&p.normal.scaled(&(-k / p.normal.norm_squared())));
    if e2 > r2 {
        Expected::Disjoint
    } else if e2 == r2 {
        point(f.to_real())
    } else {
        meets(circle(
            f.to_real(),
            p.normal.scaled(&n(sigma)),
            r2 - e2,
            p.ref_dir.clone(),
            Contact::Crossing,
        ))
    }
}
fn cc(a: &Cylinder, b: &Cylinder) -> Expected {
    use Contact::*;
    if a.axis.cross(&b.axis).is_zero() {
        let d = b.origin.minus(&a.origin);
        let u = d.minus(&a.axis.scaled(&(d.dot(&a.axis) / a.axis.norm_squared())));
        let d2 = u.norm_squared();
        let ra = &a.radius;
        let rb = &b.radius;
        if d2.is_zero() {
            return if ra == rb {
                Expected::Coincident(true)
            } else {
                Expected::Disjoint
            };
        }
        if d2 > sq(&(ra + rb)) || d2 < sq(&(ra - rb)) {
            return Expected::Disjoint;
        }
        if d2 == sq(&(ra + rb)) {
            return meets(line(
                a.origin.plus(&u.scaled(&(ra / &(ra + rb)))).to_real(),
                a.axis.clone(),
                Tangent,
            ));
        }
        if d2 == sq(&(ra - rb)) {
            return meets(line(
                a.origin.plus(&u.scaled(&(ra / &(ra - rb)))).to_real(),
                a.axis.clone(),
                Tangent,
            ));
        }
        let t = (&d2 + &sq(ra) - sq(rb)) / (&n(2) * &d2);
        let v = a.axis.cross(&u);
        let f = a.origin.plus(&u.scaled(&t));
        let h = root(&((sq(ra) - sq(&t) * d2) / v.norm_squared()));
        return Expected::Meets(vec![
            line(plus(&f, &v, &h), a.axis.clone(), Crossing),
            line(plus(&f, &v, &(-h)), neg(&a.axis), Crossing),
        ]);
    }
    let (big, small) = if a.radius >= b.radius { (a, b) } else { (b, a) };
    let rb = &big.radius;
    let rs = &small.radius;
    let m = big.axis.cross(&small.axis);
    let d = small.origin.minus(&big.origin);
    let m2 = m.norm_squared();
    let fb = big
        .origin
        .plus(&big.axis.scaled(&(&d.cross(&small.axis).dot(&m) / &m2)));
    let fs = small
        .origin
        .plus(&small.axis.scaled(&(&d.cross(&big.axis).dot(&m) / &m2)));
    let delta2 = &sq(&d.dot(&m)) / &m2;
    if delta2 > sq(&(rb + rs)) {
        return Expected::Disjoint;
    }
    if delta2 == sq(&(rb + rs)) {
        return point(fb.plus(&fs.minus(&fb).scaled(&(rb / &(rb + rs)))).to_real());
    }
    if delta2.is_zero() && rb == rs {
        let k2 = big.axis.norm_squared() * small.axis.norm_squared();
        let Some(k) = square_root(&k2) else {
            return Expected::Refused;
        };
        let minus = big
            .axis
            .scaled(&small.axis.norm_squared())
            .minus(&small.axis.scaled(&k));
        let plus = big
            .axis
            .scaled(&small.axis.norm_squared())
            .plus(&small.axis.scaled(&k));
        let cos = big.axis.dot(&small.axis) / k;
        let h = root(&(&sq(rb) / &m2));
        let singular = vec![plus_point(&fb, &m, &h), plus_point(&fb, &m, &(-h))];
        let mut parts = vec![
            ellipse(
                fb.clone(),
                minus.clone(),
                plus.clone(),
                n(2) * sq(rb) / (&n(1) - &cos),
                sq(rb),
            ),
            ellipse(fb, plus, minus, n(2) * sq(rb) / (n(1) + cos), sq(rb)),
        ];
        for part in &mut parts {
            if let Part::Curve { singular: s, .. } = part {
                *s = singular.clone();
            }
        }
        return Expected::Meets(parts);
    }
    let delta = root(&delta2);
    let e3 = big
        .axis
        .to_real()
        .scaled(&(rr(1) / root(&big.axis.norm_squared())));
    let e2 = if delta2.is_zero() {
        m.to_real().scaled(&(rr(1) / root(&m2)))
    } else {
        fs.minus(&fb).to_real().scaled(&(&rr(1) / &delta))
    };
    let e1 = e2.cross(&e3);
    let ws = small
        .axis
        .to_real()
        .scaled(&(rr(1) / root(&small.axis.norm_squared())));
    let alpha = ws.dot(&e1);
    let gamma = ws.dot(&e3);
    let br = Real::from_rational(rb.clone());
    let sr = Real::from_rational(rs.clone());
    if delta2 > sq(&(rb - rs)) {
        let z = (sr.square() - (&br - &delta).square()).sqrt() / alpha;
        return meets(branch(
            fb.to_real().plus(&e2.scaled(&br)).plus(&e3.scaled(&z)),
            vec![],
        ));
    }
    let x = root(&(&sq(rb) - &delta2));
    let singular = if delta2 == sq(&(rb - rs)) {
        vec![fb.plus(&fs.minus(&fb).scaled(&(rb / &(rb - rs)))).to_real()]
    } else {
        vec![]
    };
    let wp = fs
        .to_real()
        .plus(&e1.scaled(&x))
        .plus(&e3.scaled(&(&(&(&gamma * &x) + &sr) / &alpha)));
    let wm = fs
        .to_real()
        .plus(&e1.scaled(&(-&x)))
        .plus(&e3.scaled(&((sr - gamma * x) / alpha)));
    Expected::Meets(vec![branch(wp, singular.clone()), branch(wm, singular)])
}
fn plus_point(o: &RPoint3, v: &RVector3, t: &Real) -> Point3 {
    plus(o, v, t)
}
fn cs(c: &Cylinder, s: &Sphere, sigma: i64) -> Expected {
    use Contact::*;
    let r = &c.radius;
    let big = &s.radius;
    let w2 = c.axis.norm_squared();
    let f = c.origin.plus(
        &c.axis
            .scaled(&(&s.center.minus(&c.origin).dot(&c.axis) / &w2)),
    );
    let offset = s.center.minus(&f);
    let d2 = offset.norm_squared();
    if d2.is_zero() {
        if big < r {
            return Expected::Disjoint;
        }
        if big == r {
            return meets(circle(
                s.center.to_real(),
                c.axis.clone(),
                sq(r),
                c.ref_dir.clone(),
                Tangent,
            ));
        }
        let h = root(&((sq(big) - sq(r)) / w2));
        return Expected::Meets(vec![
            circle(
                plus(&s.center, &c.axis, &h),
                c.axis.scaled(&n(-sigma)),
                sq(r),
                c.ref_dir.clone(),
                Crossing,
            ),
            circle(
                plus(&s.center, &c.axis, &(-h)),
                c.axis.scaled(&n(sigma)),
                sq(r),
                c.ref_dir.clone(),
                Crossing,
            ),
        ]);
    }
    if d2 > sq(&(r + big)) || (r > big && d2 < sq(&(r - big))) {
        return Expected::Disjoint;
    }
    let radial = offset
        .to_real()
        .scaled(&(Real::from_rational(r.clone()) / root(&d2)));
    let near = f.to_real().plus(&radial);
    if d2 == sq(&(r + big)) || (r > big && d2 == sq(&(r - big))) {
        return point(near);
    }
    let dz = (Real::from_rational(sq(big)) - (Real::from_rational(r.clone()) - root(&d2)).square())
        .sqrt()
        / root(&w2);
    let upper = near.plus(&c.axis.to_real().scaled(&dz));
    if big > r && d2 <= sq(&(big - r)) {
        let singular = if d2 == sq(&(big - r)) {
            vec![f.to_real().plus(&radial.scaled(&rr(-1)))]
        } else {
            vec![]
        };
        let lower = near.plus(&c.axis.to_real().scaled(&(-dz)));
        Expected::Meets(vec![
            branch(upper, singular.clone()),
            branch(lower, singular),
        ])
    } else {
        meets(branch(upper, vec![]))
    }
}
