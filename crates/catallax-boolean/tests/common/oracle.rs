use super::*;

pub fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
pub fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
pub fn p(x: &str, y: &str, z: &str) -> RPoint3 {
    RPoint3::new(q(x), q(y), q(z))
}
pub fn v(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(r(x), r(y), r(z))
}
pub fn sq(x: &Rational) -> Rational {
    x * x
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Location {
    Inside,
    Outside,
    On,
}

#[derive(Clone, Debug)]
pub enum Shape {
    Block([Rational; 3]),
    Cylinder(Rational, Rational),
    Ball(Rational),
}
#[derive(Clone, Debug)]
pub struct Primitive {
    pub shape: Shape,
    pub pose: Pose,
}
impl Primitive {
    pub fn build(&self, m: &mut Model) -> BodyId {
        match &self.shape {
            Shape::Block(s) => block(m, s.clone(), &self.pose),
            Shape::Cylinder(r, h) => cylinder(m, r.clone(), h.clone(), &self.pose),
            Shape::Ball(r) => sphere(m, r.clone(), &self.pose),
        }
        .unwrap()
    }
    pub fn transformed(&self, scale: &Rational, pose: &Pose) -> Self {
        let shape = match &self.shape {
            Shape::Block(s) => Shape::Block(s.each_ref().map(|x| x * scale)),
            Shape::Cylinder(r, h) => Shape::Cylinder(r * scale, h * scale),
            Shape::Ball(r) => Shape::Ball(r * scale),
        };
        Self {
            shape,
            pose: Pose::new(
                self.pose.rotation.clone(),
                self.pose.translation.scaled(scale),
            )
            .then(pose),
        }
    }
    // Every constraint is <= 0 in the closed primitive. These formulas are
    // independent of Surface::implicit, side(), and the boolean implementation.
    pub fn constraints(&self, world: &RPoint3) -> Vec<Rational> {
        let p = self.pose.inverse().apply_point(world);
        match &self.shape {
            Shape::Block(s) => vec![
                -&p.x,
                &p.x - &s[0],
                -&p.y,
                &p.y - &s[1],
                -&p.z,
                &p.z - &s[2],
            ],
            Shape::Cylinder(radius, height) => {
                vec![sq(&p.x) + sq(&p.y) - sq(radius), -&p.z, &p.z - height]
            }
            Shape::Ball(radius) => vec![sq(&p.x) + sq(&p.y) + sq(&p.z) - sq(radius)],
        }
    }
    pub fn locate(&self, p: &RPoint3) -> Location {
        let c = self.constraints(p);
        if c.iter().any(|x| x > &r(0)) {
            Location::Outside
        } else if c.iter().any(Rational::is_zero) {
            Location::On
        } else {
            Location::Inside
        }
    }
    pub fn outward(&self, world: &RPoint3) -> RVector3 {
        assert_eq!(self.locate(world), Location::On);
        let p = self.pose.inverse().apply_point(world);
        let n = match &self.shape {
            Shape::Block(s) => {
                let c = [&p.x, &p.y, &p.z];
                let mut found = Vec::new();
                for i in 0..3 {
                    if c[i].is_zero() {
                        found.push((i, -1));
                    }
                    if c[i] == &s[i] {
                        found.push((i, 1));
                    }
                }
                assert_eq!(
                    found.len(),
                    1,
                    "face interior must not lie on a primitive edge"
                );
                let mut n = [r(0), r(0), r(0)];
                n[found[0].0] = r(found[0].1);
                RVector3::new(n[0].clone(), n[1].clone(), n[2].clone())
            }
            Shape::Cylinder(_, h) => {
                if p.z.is_zero() {
                    v(0, 0, -1)
                } else if &p.z == h {
                    v(0, 0, 1)
                } else {
                    RVector3::new(p.x, p.y, r(0))
                }
            }
            Shape::Ball(_) => RVector3::new(p.x, p.y, p.z),
        };
        self.pose.apply_vector(&n)
    }
}
pub fn block_at(size: [i64; 3], t: [i64; 3]) -> Primitive {
    Primitive {
        shape: Shape::Block(size.map(r)),
        pose: translation(t),
    }
}
pub fn cyl_at(radius: i64, height: i64, t: [i64; 3], x: bool) -> Primitive {
    let rotation = if x {
        Rotation::from_quaternion(r(1), r(0), r(1), r(0)).unwrap()
    } else {
        Rotation::identity()
    };
    Primitive {
        shape: Shape::Cylinder(r(radius), r(height)),
        pose: Pose::new(rotation, v(t[0], t[1], t[2])),
    }
}
pub fn ball_at(radius: i64, t: [i64; 3]) -> Primitive {
    Primitive {
        shape: Shape::Ball(r(radius)),
        pose: translation(t),
    }
}
pub fn translation(t: [i64; 3]) -> Pose {
    Pose::new(Rotation::identity(), v(t[0], t[1], t[2]))
}
pub fn moved_point(p: &RPoint3, scale: &Rational, pose: &Pose) -> RPoint3 {
    pose.apply_point(&RPoint3::new(&p.x * scale, &p.y * scale, &p.z * scale))
}
pub fn set_rule(op: BoolOp, a: Location, b: Location) -> Option<Location> {
    if a == Location::On || b == Location::On {
        return None;
    }
    let (a, b) = (a == Location::Inside, b == Location::Inside);
    let inside = match op {
        BoolOp::Unite => a || b,
        BoolOp::Intersect => a && b,
        BoolOp::Subtract => a && !b,
        _ => panic!("new operation"),
    };
    Some(if inside {
        Location::Inside
    } else {
        Location::Outside
    })
}

// Normal steps can stay on discarded operand boundaries. Decide strict
// result membership from neighbourhoods, not from boundary-free set logic.
pub fn regularized_step(op: BoolOp, a: &Primitive, b: &Primitive, p: &RPoint3) -> Option<Location> {
    use Location::{Inside, On, Outside};
    let (la, lb) = (a.locate(p), b.locate(p));
    match (op, la, lb) {
        (BoolOp::Unite, Inside, _) | (BoolOp::Unite, _, Inside) => Some(Inside),
        (BoolOp::Intersect, Outside, _) | (BoolOp::Intersect, _, Outside) => Some(Outside),
        (BoolOp::Subtract, Outside, _) | (BoolOp::Subtract, _, Inside) => Some(Outside),
        (BoolOp::Unite, On, On)
            if matches!((&a.shape, &b.shape), (Shape::Block(_), Shape::Block(_)))
                && [a, b]
                    .iter()
                    .all(|s| s.constraints(p).iter().filter(|c| c.is_zero()).count() == 1) =>
        {
            // One active planar constraint on each box. Opposite normals
            // bound complementary half-neighbourhoods, so their union is open
            // at p even though both source predicates say On.
            let na = a.outward(p);
            let nb = b.outward(p);
            (na.cross(&nb).is_zero() && na.dot(&nb) < r(0)).then_some(Inside)
        }
        _ => set_rule(op, la, lb),
    }
}

// Evaluate a surface's rational implicit and outward gradient without calling
// geometric decision code. The positive side is the definition's normal side.
pub fn implicit_normal(s: &Surface, p: &RPoint3) -> (Rational, RVector3) {
    match s {
        Surface::Plane(s) => (p.minus(&s.origin).dot(&s.normal), s.normal.clone()),
        Surface::Sphere(s) => {
            let d = p.minus(&s.center);
            (d.norm_squared() - sq(&s.radius), d.scaled(&r(2)))
        }
        Surface::Cylinder(s) => {
            let d = p.minus(&s.origin);
            let aa = s.axis.norm_squared();
            let radial = d.scaled(&aa).minus(&s.axis.scaled(&d.dot(&s.axis)));
            (
                d.cross(&s.axis).norm_squared() - sq(&s.radius) * aa,
                radial.scaled(&r(2)),
            )
        }
        _ => panic!("unsupported oracle surface"),
    }
}

// For f(p+t*n)=c+l*t+a*t^2, t<=1 and
// t<|c|/(2*(|l|+|a|+1)) preserve every nonzero constraint's sign.
// For c=0,l!=0, t<|l|/(2*(|a|+1)) preserves the first-order sign.
// This is a certified clearance bound, not a geometric tolerance or sampling.
pub fn safe_step(a: &Primitive, b: &Primitive, p: &RPoint3, n: &RVector3) -> Rational {
    let mut t = r(1);
    for primitive in [a, b] {
        let c = primitive.constraints(p);
        let plus = primitive.constraints(&p.plus(n));
        let minus = primitive.constraints(&p.plus(&n.scaled(&r(-1))));
        for ((c, plus), minus) in c.into_iter().zip(plus).zip(minus) {
            let linear = (&plus - &minus) / r(2);
            let quadratic = (plus + minus - &r(2) * &c) / r(2);
            let abs = |x: &Rational| if x < &r(0) { -x } else { x.clone() };
            if !c.is_zero() {
                t = t.min(abs(&c) / (r(2) * (abs(&linear) + abs(&quadratic) + r(1))));
            } else if !linear.is_zero() {
                t = t.min(abs(&linear) / (r(2) * (abs(&quadratic) + r(1))));
            }
        }
    }
    t / r(2)
}

pub struct Rng(pub u64);
impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }
    pub fn point(&mut self) -> RPoint3 {
        let mut coord = || {
            Rational::from_ratio(
                (self.next() % 101) as i64 - 30,
                (self.next() % 11 + 1) as i64,
            )
            .unwrap()
        };
        RPoint3::new(coord(), coord(), coord())
    }
    pub fn pose(&mut self) -> Pose {
        let w = r((self.next() % 100 + 1) as i64);
        let x = r((self.next() % 201) as i64 - 100);
        let y = r((self.next() % 201) as i64 - 100);
        let z = r((self.next() % 201) as i64 - 100);
        let rotation = Rotation::from_quaternion(w, x, y, z).unwrap();
        let mut coord = || {
            Rational::from_ratio(
                (self.next() % 2_000_001) as i64 - 1_000_000,
                (self.next() % 17 + 1) as i64,
            )
            .unwrap()
        };
        Pose::new(rotation, RVector3::new(coord(), coord(), coord()))
    }
}
