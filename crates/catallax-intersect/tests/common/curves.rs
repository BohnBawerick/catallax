use super::*;
#[derive(Clone, Debug)]
pub enum Carrier {
    Analytic(Curve),
    Branch {
        a: Surface,
        b: Box<Surface>,
        witness: Point3,
    },
}
impl Carrier {
    pub fn add(&self, g: &mut GeomStore) -> CurveId {
        match self {
            Self::Analytic(c) => g.add_curve(c.clone()),
            Self::Branch { a, b, witness } => {
                let a = g.add_surface(a.clone());
                let b = g.add_surface((**b).clone());
                let witness = closed(g, witness.clone(), &[a, b]);
                g.add_curve(Curve::Intersection(IntersectionCurve {
                    surfaces: [a, b],
                    witness,
                }))
            }
        }
    }
    pub fn transformed(&self, scale: &Rational, pose: &Pose) -> Self {
        let k = Real::from_rational(scale.clone());
        let point =
            |p: &Point3| pose.apply_real_point(&Point3::new(&p.x * &k, &p.y * &k, &p.z * &k));
        match self {
            Self::Analytic(c) => Self::Analytic(scaled_curve(c, scale).transformed(pose).unwrap()),
            Self::Branch { a, b, witness } => Self::Branch {
                a: scaled_surface(a, scale).transformed(pose),
                b: Box::new(scaled_surface(b, scale).transformed(pose)),
                witness: point(witness),
            },
        }
    }
}
pub fn viviani(upper: bool) -> Carrier {
    Carrier::Branch {
        a: zc(1),
        b: Box::new(sphere(pi([1, 0, 0]), n(2))),
        witness: pi([1, 0, if upper { 2 } else { -2 }]).to_real(),
    }
}
pub fn loop25() -> Carrier {
    Carrier::Branch {
        a: zc(1),
        b: Box::new(sphere(pi([2, 0, 0]), n(2))),
        witness: Point3::new(rr(1), rr(0), rr(3).sqrt()),
    }
}
pub fn upper26() -> Carrier {
    Carrier::Branch {
        a: zc(1),
        b: Box::new(sphere(pi([1, 0, 0]), n(3))),
        witness: pi([1, 0, 3]).to_real(),
    }
}
#[derive(Clone, Debug)]
pub enum Hits {
    Disjoint,
    Contained,
    Meets(Vec<(Point3, Contact)>),
}
#[derive(Clone, Debug)]
pub struct HitCase {
    pub name: String,
    pub curve: Carrier,
    pub surface: Surface,
    pub expected: Hits,
}
fn points(p: Vec<Point3>, c: Contact) -> Hits {
    Hits::Meets(p.into_iter().map(|p| (p, c)).collect())
}
fn rational_hits(p: &[[i64; 3]], c: Contact) -> Hits {
    points(p.iter().map(|p| pi(*p).to_real()).collect(), c)
}
pub fn hits() -> Vec<HitCase> {
    use Contact::*;
    let mut out = vec![];
    let mut add = |name: &str, curve: Carrier, surface: Surface, expected: Hits| {
        out.push(HitCase {
            name: name.into(),
            curve,
            surface,
            expected,
        })
    };
    let l = || Carrier::Analytic(line28());
    let k = || Carrier::Analytic(k5());
    let e = || Carrier::Analytic(ellipse4());
    for s in [zc(5), sphere(pi([0, 0, 0]), n(5))] {
        add(
            "28 line secant",
            l(),
            s.clone(),
            rational_hits(&[[-5, 0, 0], [5, 0, 0]], Crossing),
        );
        let tangent = Curve::Line(
            Line::new(
                (pi([-10, 5, 0])).into(),
                (vi([1, 0, 0])).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        );
        add(
            "28 line tangent",
            Carrier::Analytic(tangent),
            s.clone(),
            rational_hits(&[[0, 5, 0]], Tangent),
        );
        let miss = Curve::Line(
            Line::new(
                (RPoint3::new(n(-10), n(5) + r("1e-30"), n(0))).into(),
                (vi([1, 0, 0])).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        );
        add(
            "28 line near miss",
            Carrier::Analytic(miss),
            s,
            Hits::Disjoint,
        );
    }
    add(
        "28 line plane",
        l(),
        axis_plane(0, n(3)),
        rational_hits(&[[3, 0, 0]], Crossing),
    );
    add(
        "28 line contained",
        l(),
        axis_plane(1, n(0)),
        Hits::Contained,
    );
    add(
        "29 circle plane",
        k(),
        axis_plane(0, n(3)),
        rational_hits(&[[3, 4, 0], [3, -4, 0]], Crossing),
    );
    add(
        "29 circle tangent origin",
        k(),
        axis_plane(0, n(5)),
        rational_hits(&[[5, 0, 0]], Tangent),
    );
    let y = rr(5) * rr(3).sqrt() / rr(2);
    add(
        "29 circle sphere",
        k(),
        sphere(pi([5, 0, 0]), n(5)),
        points(
            vec![
                Point3::new(rr(5) / rr(2), y.clone(), rr(0)),
                Point3::new(rr(5) / rr(2), -y, rr(0)),
            ],
            Crossing,
        ),
    );
    add(
        "29 circle cylinder four",
        k(),
        cyl(pi([0, 0, 0]), vi([1, 0, 0]), n(3)),
        rational_hits(&[[4, 3, 0], [-4, 3, 0], [-4, -3, 0], [4, -3, 0]], Crossing),
    );
    add(
        "29 circle contained plane",
        k(),
        axis_plane(2, n(0)),
        Hits::Contained,
    );
    add(
        "29 circle contained cylinder without fact",
        k(),
        zc(5),
        Hits::Contained,
    );
    // Eliminate y² from x²+y²=25 and y²+2(x-4)²=8:
    // x²-16x+49=0. Only 8-sqrt(15) is in [-5,5].
    let x = rr(8) - rr(15).sqrt();
    let y = (rr(25) - x.square()).sqrt();
    add(
        "30 tilted cylinder secant",
        k(),
        cyl(pi([4, 0, 0]), vi([0, 1, 1]), n(2)),
        points(
            vec![
                Point3::new(x.clone(), y.clone(), rr(0)),
                Point3::new(x, -y, rr(0)),
            ],
            Crossing,
        ),
    );
    add(
        "30 tilted cylinder tangent",
        k(),
        cyl(pi([3, 0, 0]), vi([0, 1, 1]), n(2)),
        rational_hits(&[[5, 0, 0]], Tangent),
    );
    add(
        "31 ellipse major vertices",
        e(),
        axis_plane(0, n(0)),
        points(
            vec![
                p(["0", "-5", "15/4"]).to_real(),
                p(["0", "5", "-15/4"]).to_real(),
            ],
            Crossing,
        ),
    );
    add(
        "31 ellipse minor vertices",
        e(),
        axis_plane(1, n(0)),
        rational_hits(&[[5, 0, 0], [-5, 0, 0]], Crossing),
    );
    let y = rr(3).sqrt() / rr(2);
    add(
        "32 upper branch secant",
        viviani(true),
        axis_plane(2, n(1)),
        points(
            vec![
                Point3::new(-rr(1) / rr(2), -y.clone(), rr(1)),
                Point3::new(-rr(1) / rr(2), y, rr(1)),
            ],
            Crossing,
        ),
    );
    add(
        "32 upper branch node",
        viviani(true),
        axis_plane(2, n(0)),
        rational_hits(&[[-1, 0, 0]], Tangent),
    );
    add(
        "32 upper branch disjoint",
        viviani(true),
        axis_plane(2, n(3)),
        Hits::Disjoint,
    );
    add(
        "32 wrong branch",
        viviani(false),
        axis_plane(2, n(1)),
        Hits::Disjoint,
    );
    out
}
