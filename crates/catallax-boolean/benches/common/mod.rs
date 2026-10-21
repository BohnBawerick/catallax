//! Fixtures derived only from the numbered cases in the piece 5 contract.

use catallax_numbers::{Rational, Refusal};
use catallax_shapes::{Pose, RPoint3, Rotation};
use catallax_topology::{BodyId, Model, block, cone, cylinder, sphere};

pub(crate) fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn p(t: [i64; 3]) -> RPoint3 {
    RPoint3::new(r(t[0]), r(t[1]), r(t[2]))
}
#[derive(Clone)]
pub(crate) enum Kind {
    Block([Rational; 3]),
    Cylinder(Rational, Rational),
    Sphere(Rational),
    Cone(Rational, Rational, Rational),
}
#[derive(Clone)]
pub(crate) struct Primitive {
    pub(crate) kind: Kind,
    pub(crate) pose: Pose,
}
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Position {
    Inside,
    Outside,
    Boundary,
}
impl Primitive {
    pub(crate) fn build(&self, model: &mut Model) -> Result<BodyId, Refusal> {
        match &self.kind {
            Kind::Block(s) => block(model, s.clone(), &self.pose),
            Kind::Cylinder(r, h) => cylinder(model, r.clone(), h.clone(), &self.pose),
            Kind::Sphere(r) => sphere(model, r.clone(), &self.pose),
            Kind::Cone(base, top, height) => {
                cone(model, base.clone(), top.clone(), height.clone(), &self.pose)
            }
        }
    }
    pub(crate) fn classify(&self, point: &RPoint3) -> Position {
        let p = self.pose.inverse().apply_point(point);
        match &self.kind {
            Kind::Block(size) => {
                let x = [p.x, p.y, p.z];
                if x.iter().zip(size).any(|(x, s)| *x < r(0) || x > s) {
                    Position::Outside
                } else if x.iter().zip(size).any(|(x, s)| x.is_zero() || x == s) {
                    Position::Boundary
                } else {
                    Position::Inside
                }
            }
            Kind::Cylinder(radius, height) => {
                let value = &p.x * &p.x + &p.y * &p.y - radius * radius;
                if p.z < r(0) || p.z > *height || value > r(0) {
                    Position::Outside
                } else if p.z.is_zero() || p.z == *height || value.is_zero() {
                    Position::Boundary
                } else {
                    Position::Inside
                }
            }
            Kind::Sphere(radius) => {
                let value = &p.x * &p.x + &p.y * &p.y + &p.z * &p.z - radius * radius;
                if value < r(0) {
                    Position::Inside
                } else if value > r(0) {
                    Position::Outside
                } else {
                    Position::Boundary
                }
            }
            Kind::Cone(base, top, height) => {
                if p.z < r(0) || p.z > *height {
                    return Position::Outside;
                }
                let radius = (base * &(height - &p.z) + top * &p.z) / height.clone();
                let value = &p.x * &p.x + &p.y * &p.y - &radius * &radius;
                if value > r(0) {
                    Position::Outside
                } else if value.is_zero() || p.z.is_zero() || p.z == *height {
                    Position::Boundary
                } else {
                    Position::Inside
                }
            }
        }
    }
    pub(crate) fn bounds(&self) -> ([Rational; 3], [Rational; 3]) {
        let (lo, hi) = match &self.kind {
            Kind::Block(s) => ([r(0), r(0), r(0)], s.clone()),
            Kind::Cylinder(radius, h) => (
                [-radius, -radius, r(0)],
                [radius.clone(), radius.clone(), h.clone()],
            ),
            Kind::Sphere(radius) => (
                [-radius, -radius, -radius],
                [radius.clone(), radius.clone(), radius.clone()],
            ),
            Kind::Cone(base, top, height) => {
                let radius = base.max(top);
                (
                    [-radius, -radius, r(0)],
                    [radius.clone(), radius.clone(), height.clone()],
                )
            }
        };
        let points: Vec<_> = (0..8)
            .map(|i| {
                let [x, y, z] = std::array::from_fn(|axis| {
                    if i & (1 << axis) == 0 {
                        lo[axis].clone()
                    } else {
                        hi[axis].clone()
                    }
                });
                self.pose.apply_point(&RPoint3::new(x, y, z))
            })
            .collect();
        let lower = std::array::from_fn(|i| {
            points
                .iter()
                .map(|p| [&p.x, &p.y, &p.z][i].clone())
                .min()
                .unwrap()
        });
        let upper = std::array::from_fn(|i| {
            points
                .iter()
                .map(|p| [&p.x, &p.y, &p.z][i].clone())
                .max()
                .unwrap()
        });
        (lower, upper)
    }
    pub(crate) fn transformed(&self, scale: &Rational, placement: &Pose) -> Self {
        let kind = match &self.kind {
            Kind::Block(s) => Kind::Block(s.clone().map(|x| x * scale.clone())),
            Kind::Cylinder(r, h) => Kind::Cylinder(r * scale, h * scale),
            Kind::Sphere(r) => Kind::Sphere(r * scale),
            Kind::Cone(base, top, height) => Kind::Cone(base * scale, top * scale, height * scale),
        };
        let pose = Pose::new(
            self.pose.rotation.clone(),
            self.pose.translation.scaled(scale),
        )
        .then(placement);
        Self { kind, pose }
    }
}
fn at(kind: Kind, t: [i64; 3], x_axis: bool) -> Primitive {
    let rotation = if x_axis {
        Rotation::from_quaternion(r(1), r(0), r(1), r(0)).unwrap()
    } else {
        Rotation::identity()
    };
    Primitive {
        kind,
        pose: Pose::new(rotation, p(t).minus(&RPoint3::origin())),
    }
}
pub(crate) fn block_at(size: [i64; 3], t: [i64; 3]) -> Primitive {
    at(Kind::Block(size.map(r)), t, false)
}
pub(crate) fn cylinder_at(radius: i64, height: i64, t: [i64; 3], x_axis: bool) -> Primitive {
    at(Kind::Cylinder(r(radius), r(height)), t, x_axis)
}
pub(crate) fn ball_at(radius: i64, t: [i64; 3]) -> Primitive {
    at(Kind::Sphere(r(radius)), t, false)
}
pub(crate) type Counts = (usize, usize, usize, usize, usize, usize);
fn c(f: usize, e: usize, rings: usize, v: usize) -> Option<Counts> {
    Some((1, 1, f, e, rings, v))
}
fn lumps(f: usize, e: usize, rings: usize, v: usize) -> Option<Counts> {
    Some((2, 2, f, e, rings, v))
}
fn void(f: usize, e: usize, rings: usize, v: usize) -> Option<Counts> {
    Some((1, 2, f, e, rings, v))
}
const EMPTY: Option<Counts> = Some((0, 0, 0, 0, 0, 0));
pub(crate) struct Case {
    pub(crate) number: u8,
    pub(crate) variant: &'static str,
    pub(crate) a: Primitive,
    pub(crate) b: Primitive,
    pub(crate) expected: [Option<Counts>; 4],
}
fn case(
    number: u8,
    variant: &'static str,
    a: Primitive,
    b: Primitive,
    expected: [Option<Counts>; 4],
) -> Case {
    Case {
        number,
        variant,
        a,
        b,
        expected,
    }
}
pub(crate) fn cases() -> Vec<Case> {
    let cube = || block_at([10, 10, 10], [0, 0, 0]);
    let b2 = || block_at([2, 2, 2], [0, 0, 0]);
    let b1 = || block_at([1, 1, 1], [0, 0, 0]);
    let ball5 = || ball_at(5, [0, 0, 0]);
    let mut cases = vec![
        case(
            1,
            "",
            b2(),
            block_at([2, 2, 2], [1, 1, 1]),
            [
                c(12, 30, 0, 20),
                c(6, 12, 0, 8),
                c(9, 21, 0, 14),
                c(9, 21, 0, 14),
            ],
        ),
        case(
            2,
            "",
            b1(),
            block_at([1, 1, 1], [1, 0, 0]),
            [c(6, 12, 0, 8), EMPTY, c(6, 12, 0, 8), None],
        ),
        case(
            3,
            "",
            b2(),
            block_at([2, 2, 2], [2, 1, 0]),
            [c(10, 24, 0, 16), EMPTY, c(6, 12, 0, 8), None],
        ),
        case(
            4,
            "",
            b1(),
            block_at([1, 1, 1], [1, 1, 0]),
            [lumps(12, 24, 0, 16), EMPTY, c(6, 12, 0, 8), None],
        ),
        case(
            5,
            "",
            b1(),
            block_at([1, 1, 1], [1, 1, 1]),
            [lumps(12, 24, 0, 16), EMPTY, None, None],
        ),
        case(
            6,
            "",
            block_at([4, 4, 4], [0, 0, 0]),
            block_at([2, 2, 6], [1, 1, -1]),
            [c(16, 36, 0, 24), c(6, 12, 0, 8), c(10, 24, 0, 16), None],
        ),
        case(
            7,
            "",
            b2(),
            block_at([2, 2, 2], [1, 0, 0]),
            [c(6, 12, 0, 8), c(6, 12, 0, 8), c(6, 12, 0, 8), None],
        ),
        case(
            8,
            "",
            b2(),
            b2(),
            [c(6, 12, 0, 8), c(6, 12, 0, 8), EMPTY, None],
        ),
        case(
            10,
            "",
            cube(),
            cylinder_at(3, 12, [5, 5, -1], false),
            [c(10, 16, 4, 8), c(3, 2, 2, 0), c(7, 14, 2, 8), None],
        ),
        case(
            11,
            "",
            cube(),
            cylinder_at(3, 5, [5, 5, 6], false),
            [None, None, c(8, 14, 2, 8), None],
        ),
        case(
            12,
            "",
            cube(),
            cylinder_at(2, 20, [-5, 5, 12], true),
            [lumps(9, 14, 2, 8), EMPTY, c(6, 12, 0, 8), None],
        ),
        case(
            14,
            "",
            cube(),
            cylinder_at(5, 12, [5, 5, -1], false),
            [
                c(16, 30, 2, 16),
                c(3, 2, 2, 0),
                Some((4, 4, 20, 36, 0, 24)),
                None,
            ],
        ),
        case(
            15,
            "",
            cube(),
            cylinder_at(3, 12, [0, 5, -1], false),
            [None, c(4, 6, 0, 4), c(8, 18, 0, 12), None],
        ),
        case(
            16,
            "",
            cube(),
            cylinder_at(2, 4, [5, 5, 3], false),
            [c(6, 12, 0, 8), c(3, 2, 2, 0), void(9, 14, 2, 8), None],
        ),
        case(
            17,
            "",
            cube(),
            ball_at(3, [10, 10, 10]),
            [c(7, 15, 0, 10), c(4, 6, 0, 4), c(7, 15, 0, 10), None],
        ),
        case(
            18,
            "",
            cube(),
            ball_at(3, [5, 5, 10]),
            [c(7, 13, 1, 8), c(2, 1, 1, 0), c(7, 13, 1, 8), None],
        ),
        case(
            19,
            "tangent",
            cube(),
            ball_at(3, [5, 5, 13]),
            [lumps(7, 12, 0, 8), EMPTY, c(6, 12, 0, 8), None],
        ),
        case(
            20,
            "",
            cube(),
            ball_at(3, [5, 5, 7]),
            [c(6, 12, 0, 8), c(1, 0, 0, 0), void(7, 12, 0, 8), None],
        ),
        case(
            21,
            "",
            cube(),
            ball_at(5, [5, 5, 5]),
            [c(6, 12, 0, 8), c(1, 0, 0, 0), void(7, 12, 0, 8), None],
        ),
        case(
            22,
            "",
            cube(),
            ball_at(3, [10, 5, 5]),
            [None, c(2, 1, 1, 0), c(7, 13, 1, 8), None],
        ),
        case(
            23,
            "",
            cylinder_at(2, 10, [0, 0, -5], false),
            cylinder_at(1, 10, [-5, 0, 0], true),
            [
                c(7, 6, 6, 0),
                c(3, 2, 2, 0),
                c(4, 4, 4, 0),
                lumps(6, 4, 4, 0),
            ],
        ),
        case(
            24,
            "",
            cylinder_at(1, 4, [0, 0, -2], false),
            cylinder_at(1, 4, [-2, 0, 0], true),
            [c(8, 8, 4, 2), c(4, 4, 0, 2), lumps(6, 6, 2, 4), None],
        ),
        case(
            25,
            "",
            cylinder_at(5, 10, [0, 0, 0], false),
            cylinder_at(5, 10, [8, 0, 0], false),
            [c(4, 6, 0, 4), c(4, 6, 0, 4), c(4, 6, 0, 4), None],
        ),
        case(
            26,
            "",
            cylinder_at(5, 10, [0, 0, 0], false),
            cylinder_at(5, 10, [10, 0, 0], false),
            [lumps(6, 4, 4, 0), EMPTY, None, None],
        ),
        case(
            27,
            "stack",
            cylinder_at(2, 5, [0, 0, 0], false),
            cylinder_at(2, 5, [0, 0, 5], false),
            [c(3, 2, 2, 0), EMPTY, None, None],
        ),
        case(
            27,
            "tube",
            cylinder_at(3, 10, [0, 0, 0], false),
            cylinder_at(1, 12, [0, 0, -1], false),
            [c(7, 6, 6, 0), c(3, 2, 2, 0), c(4, 4, 4, 0), None],
        ),
        case(
            28,
            "",
            cylinder_at(2, 10, [0, 0, -5], false),
            cylinder_at(1, 10, [-5, 1, 0], true),
            [c(7, 6, 4, 1), c(3, 2, 0, 1), None, lumps(6, 4, 2, 2)],
        ),
        case(
            29,
            "tangent",
            cylinder_at(2, 10, [0, 0, -5], false),
            cylinder_at(1, 10, [-5, 3, 0], true),
            [lumps(6, 4, 4, 0), None, None, None],
        ),
        case(
            30,
            "radius5",
            cylinder_at(3, 20, [0, 0, -10], false),
            ball5(),
            [
                c(5, 4, 4, 0),
                c(3, 2, 2, 0),
                lumps(6, 4, 4, 0),
                c(2, 2, 2, 0),
            ],
        ),
        case(
            30,
            "radius6",
            cylinder_at(3, 20, [0, 0, -10], false),
            ball_at(6, [0, 0, 0]),
            [
                c(5, 4, 4, 0),
                c(3, 2, 2, 0),
                lumps(6, 4, 4, 0),
                c(2, 2, 2, 0),
            ],
        ),
        case(
            31,
            "",
            cylinder_at(5, 20, [0, 0, -10], false),
            ball5(),
            [c(3, 2, 2, 0), c(1, 0, 0, 0), lumps(6, 4, 4, 0), EMPTY],
        ),
        case(
            32,
            "",
            cylinder_at(1, 6, [0, 0, -3], false),
            ball_at(2, [1, 0, 0]),
            [c(5, 4, 2, 1), c(3, 2, 0, 1), lumps(6, 4, 2, 2), None],
        ),
        case(
            33,
            "",
            cylinder_at(1, 6, [0, 0, -3], false),
            ball_at(2, [2, 0, 0]),
            [c(4, 3, 3, 0), c(2, 1, 1, 0), c(4, 3, 3, 0), c(2, 1, 1, 0)],
        ),
        case(
            34,
            "",
            cylinder_at(1, 6, [0, 0, -3], false),
            ball_at(2, [3, 0, 0]),
            [lumps(4, 2, 2, 0), None, None, None],
        ),
        case(
            35,
            "",
            cylinder_at(1, 6, [0, 0, -3], false),
            ball_at(1, [1, 0, 0]),
            [c(4, 3, 3, 0), c(2, 1, 1, 0), c(4, 3, 3, 0), c(2, 1, 1, 0)],
        ),
        case(
            36,
            "",
            ball5(),
            ball_at(5, [6, 0, 0]),
            [c(2, 1, 1, 0), c(2, 1, 1, 0), c(2, 1, 1, 0), None],
        ),
        case(
            37,
            "tangent",
            ball5(),
            ball_at(5, [10, 0, 0]),
            [lumps(2, 0, 0, 0), EMPTY, None, None],
        ),
        case(
            38,
            "",
            ball5(),
            ball_at(2, [3, 0, 0]),
            [c(1, 0, 0, 0), c(1, 0, 0, 0), void(2, 0, 0, 0), EMPTY],
        ),
        case(
            39,
            "",
            ball5(),
            ball_at(3, [0, 0, 0]),
            [c(1, 0, 0, 0), c(1, 0, 0, 0), void(2, 0, 0, 0), None],
        ),
        case(
            40,
            "",
            ball5(),
            ball5(),
            [c(1, 0, 0, 0), c(1, 0, 0, 0), EMPTY, None],
        ),
        case(
            46,
            "cone",
            cube(),
            at(Kind::Cone(r(2), r(1), r(2)), [0, 0, 0], false),
            [None, None, None, None],
        ),
        case(
            47,
            "",
            cube(),
            cylinder_at(5, 6, [5, 5, 2], false),
            [None, None, None, None],
        ),
    ];
    let e = Rational::parse("1e-30").unwrap();
    for (name, sign) in [("gap", 1), ("overlap", -1)] {
        let mut b = block_at([1, 1, 1], [1, 0, 0]);
        b.pose.translation.x = &b.pose.translation.x + &(&e * &r(sign));
        cases.push(case(
            9,
            name,
            b1(),
            b,
            if sign == 1 {
                [lumps(12, 24, 0, 16), EMPTY, c(6, 12, 0, 8), None]
            } else {
                [c(6, 12, 0, 8), c(6, 12, 0, 8), c(6, 12, 0, 8), None]
            },
        ));
        let mut b = cylinder_at(2, 20, [-5, 5, 12], true);
        b.pose.translation.z = &b.pose.translation.z + &(&e * &r(sign));
        cases.push(case(
            13,
            name,
            cube(),
            b,
            if sign == 1 {
                [lumps(9, 14, 2, 8), EMPTY, c(6, 12, 0, 8), None]
            } else {
                [c(10, 20, 2, 12), c(4, 6, 0, 4), c(8, 18, 0, 12), None]
            },
        ));
        let mut b = ball_at(3, [5, 5, 13]);
        b.pose.translation.z = &b.pose.translation.z + &(&e * &r(sign));
        cases.push(case(
            19,
            name,
            cube(),
            b,
            if sign == 1 {
                [lumps(7, 12, 0, 8), EMPTY, c(6, 12, 0, 8), None]
            } else {
                [c(7, 13, 1, 8), c(2, 1, 1, 0), c(7, 13, 1, 8), None]
            },
        ));
        let mut b = ball_at(5, [10, 0, 0]);
        b.pose.translation.x = &b.pose.translation.x + &(&e * &r(sign));
        cases.push(case(
            37,
            name,
            ball5(),
            b,
            if sign == 1 {
                [lumps(2, 0, 0, 0), EMPTY, None, None]
            } else {
                [c(2, 1, 1, 0), None, None, None]
            },
        ));
    }
    let mut b = cylinder_at(1, 10, [-5, 3, 0], true);
    b.pose.translation.y = &b.pose.translation.y - &e;
    cases.push(case(
        29,
        "bite",
        cylinder_at(2, 10, [0, 0, -5], false),
        b,
        [c(6, 5, 5, 0), None, None, None],
    ));
    cases.sort_by_key(|c| (c.number, c.variant));
    cases
}
