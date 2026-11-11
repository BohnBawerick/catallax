//! Rational patches and untrusted oracle candidates shared by section tests.
use catallax_numbers::{Decider, Rational, Real};
use catallax_shapes::{
    BSplineSurface, Cylinder, Point3, RationalPatch, SectionCandidate, SectionSurface, Sphere,
    Vector3,
};

pub fn r(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
pub fn p(x: &str, y: &str, z: &str) -> Point3 {
    Point3::new(r(x).into(), r(y).into(), r(z).into())
}

/// One span of bidegree (2,3): a quadratic profile carried through four stations.
pub fn loft_like(dec: &mut Decider) -> RationalPatch {
    let bulge = ["1", "2", "1/2", "3/2"];
    let mut points = Vec::new();
    let mut weights = Vec::new();
    for (i, x) in ["0", "1", "2"].iter().enumerate() {
        for (j, z) in ["0", "1", "2", "3"].iter().enumerate() {
            let y = if i == 1 { bulge[j] } else { "0" };
            points.push(p(x, y, z));
            weights.push(if i == 1 {
                Real::from_rational(r("3/4"))
            } else {
                Real::one()
            });
        }
    }
    let carrier = BSplineSurface::new(
        [2, 3],
        [
            ["0", "0", "0", "1", "1", "1"].map(r).to_vec(),
            ["0", "0", "0", "0", "1", "1", "1", "1"].map(r).to_vec(),
        ],
        [3, 4],
        points,
        weights,
        dec,
    )
    .unwrap();
    RationalPatch::new(carrier, dec).unwrap()
}

pub fn height_patch(hill: bool, dec: &mut Decider) -> RationalPatch {
    let mut controls = Vec::new();
    for i in 0..6 {
        for j in 0..6 {
            let z = if hill {
                Rational::from_i64(i64::from((2..=3).contains(&i) && (2..=3).contains(&j)))
            } else {
                Rational::from_f64_exact(
                    include_str!("../fixtures/patch-plane/wavy-heights.txt")
                        .split_whitespace()
                        .skip(1)
                        .nth((i * 6 + j) as usize)
                        .unwrap()
                        .parse()
                        .unwrap(),
                )
                .unwrap()
            };
            controls.push(Point3::new(Real::from_i64(i), Real::from_i64(j), z.into()));
        }
    }
    let knots = ["0", "0", "0", "0", "1", "2", "3", "3", "3", "3"]
        .map(r)
        .to_vec();
    RationalPatch::new(
        BSplineSurface::new(
            [3, 3],
            [knots.clone(), knots],
            [6, 6],
            controls,
            vec![Real::one(); 36],
            dec,
        )
        .unwrap(),
        dec,
    )
    .unwrap()
}

pub fn equal(a: &Point3, b: &Point3, dec: &mut Decider) {
    assert!(
        dec.sign(&a.minus(b).norm_squared()).is_coincident(),
        "{a:?} != {b:?}"
    );
}

pub fn candidates(text: &str) -> Vec<SectionCandidate> {
    text.split("CURVE ")
        .skip(1)
        .map(|curve| {
            assert!(curve.lines().next().unwrap().ends_with("DEGREE 3"));
            let rows = curve.lines().collect::<Vec<_>>();
            let values = |name: &str| {
                rows.iter()
                    .find(|row| row.starts_with(name))
                    .unwrap()
                    .split_whitespace()
                    .skip(1)
                    .map(|s| Rational::from_f64_exact(s.parse::<f64>().unwrap()).unwrap())
                    .collect::<Vec<_>>()
            };
            let knots = values("KNOTS");
            let points = values("POINTS");
            assert_eq!(points.len() % 3, 0);
            let controls = points
                .as_chunks::<3>()
                .0
                .iter()
                .map(|p| {
                    assert_eq!(p[2], Rational::one());
                    [p[0].clone(), p[1].clone()]
                })
                .collect::<Vec<_>>();
            assert_eq!(knots.len(), controls.len() + 4);
            assert_eq!((controls.len() - 1) % 3, 0);
            assert!(knots[..4].iter().all(|k| *k == knots[0]));
            assert!(
                knots[knots.len() - 4..]
                    .iter()
                    .all(|k| *k == knots[knots.len() - 1])
            );
            for triple in knots[4..knots.len() - 4].as_chunks::<3>().0.iter() {
                assert!(triple.iter().all(|k| *k == triple[0]));
            }
            SectionCandidate {
                pieces: (0..(controls.len() - 1) / 3)
                    .map(|i| std::array::from_fn(|j| controls[3 * i + j].clone()))
                    .collect(),
            }
        })
        .collect()
}

pub fn axis(x: i64, y: i64, z: i64) -> Vector3 {
    Vector3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
}

/// Loft cuts have radius 5/2 about the origin, with exact ends (0,0,5/2)
/// and (2,0,3/2). Hill cuts have radius 1 and their lowest point at 3/10,
/// one millionth below the hilltop, at it, or one millionth above it.
pub fn cutter(surface: &str, cut: &str) -> SectionSurface {
    let peak = r("225/256");
    let low = match cut {
        "loft" => Rational::zero(),
        "hill" => r("3/10"),
        "tiny" => &peak - &r("1/1000000"),
        "tangent" => peak,
        _ => &peak + &r("1/1000000"),
    };
    let center_z = Real::from_rational(&low + &Rational::one());
    match (surface, cut) {
        ("cylinder", "loft") => SectionSurface::Cylinder(Cylinder {
            origin: p("0", "0", "0"),
            axis: axis(0, 1, 0),
            radius_sq: Real::from_rational(r("25/4")),
            ref_dir: axis(0, 0, 1),
        }),
        ("cylinder", _) => SectionSurface::Cylinder(Cylinder {
            origin: catallax_shapes::Point3::new(Real::zero(), r("5/2").into(), center_z),
            axis: axis(1, 0, 0),
            radius_sq: Real::one(),
            ref_dir: axis(0, 1, 0),
        }),
        ("sphere", "loft") => SectionSurface::Sphere(Sphere {
            center: p("0", "0", "0"),
            radius: Real::from_rational(r("5/2")),
            axis: axis(0, 0, 1),
            ref_dir: axis(1, 0, 0),
        }),
        _ => SectionSurface::Sphere(Sphere {
            center: catallax_shapes::Point3::new(r("5/2").into(), r("5/2").into(), center_z),
            radius: Real::one(),
            axis: axis(0, 0, 1),
            ref_dir: axis(1, 0, 0),
        }),
    }
}
