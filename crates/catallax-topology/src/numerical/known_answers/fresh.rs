//! Seeded random edges, each with one hidden bump of known size and place.
//!
//! This is the first guard against tuning a checker to the fixed edges: every
//! run draws a new seed, so the bump widths, positions, knot vectors, degrees
//! and surfaces differ each time. The seed is printed in the receipt, and
//! `CATALLAX_GAP_SEED` replays it. None of the builders below is shared with
//! the fixed edges; only the rational tools are.

use super::exact::*;
use super::*;

/// SplitMix64 (Steele, Lea and Flood, "Fast splittable pseudorandom number
/// generators", OOPSLA 2014), with the constants of Vigna's public-domain
/// `splitmix64.c`. The standard library has no seedable generator.
pub(super) struct Random(u64);

impl Random {
    pub(super) fn new(seed: u64) -> Self {
        Self(seed)
    }
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// A whole number in `lo..=hi`.
    fn between(&mut self, lo: i64, hi: i64) -> i64 {
        lo + (self.next() % (hi - lo + 1) as u64) as i64
    }
    /// A dyadic fraction strictly inside `(0, 1)` with `bits` binary places.
    fn share(&mut self, bits: i64) -> Rational {
        qi(self.between(1, (1 << bits) - 1)).mul_pow2(-bits)
    }
    fn pick<T: Copy>(&mut self, options: &[T]) -> T {
        options[self.between(0, options.len() as i64 - 1) as usize]
    }
}

/// The seed of this run: `CATALLAX_GAP_SEED` when set, otherwise the clock
/// mixed with the process id.
pub(super) fn seed() -> u64 {
    if let Ok(text) = std::env::var("CATALLAX_GAP_SEED") {
        return text.parse().expect("CATALLAX_GAP_SEED is a whole number");
    }
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("clock after 1970");
    Random::new(now.as_nanos() as u64 ^ (u64::from(std::process::id()) << 32)).next()
}

/// A bump size between 1e-7 and 1e-2 mm.
fn bump(random: &mut Random) -> Rational {
    qi(random.between(10, 99)) * qi(10).pow(-(random.between(4, 8) as i32)).expect("power")
}

/// `least..=most` distinct sorted dyadic knots strictly inside `(0, 1)`.
fn random_knots(random: &mut Random, least: i64, most: i64) -> Vec<Rational> {
    let count = random.between(least, most) as usize;
    let mut knots = std::collections::BTreeSet::new();
    while knots.len() < count {
        knots.insert(random.share(10));
    }
    knots.into_iter().collect()
}

/// `degree + 2` distinct knots of random spacing in a window 2^-3 to 2^-20
/// wide that lies strictly inside one span of the interior `knots`.
fn cluster(random: &mut Random, degree: usize, knots: &[Rational]) -> Vec<Rational> {
    loop {
        let width = pow2(-random.between(3, 20));
        let start = random.share(30);
        let end = &start + &width;
        let clear = end < Rational::one() && knots.iter().all(|knot| knot < &start || knot > &end);
        if !clear {
            continue;
        }
        let mut shares = std::collections::BTreeSet::new();
        while shares.len() < degree {
            shares.insert(random.share(7));
        }
        let mut cluster = vec![start.clone()];
        cluster.extend(shares.iter().map(|share| &start + &(&width * share)));
        cluster.push(end);
        return cluster;
    }
}

/// An axis-aligned frame `[e1, e2, normal]` of a plane, as a cyclic
/// permutation of the coordinate axes, so the frame is exactly rational.
fn frame(random: &mut Random) -> [[i64; 3]; 3] {
    let turn = random.between(0, 2) as usize;
    std::array::from_fn(|k| std::array::from_fn(|axis| i64::from(axis == (k + turn) % 3)))
}

fn offset(random: &mut Random) -> [Rational; 3] {
    let reach = random.pick(&[10, 1_000, 1_000_000]);
    std::array::from_fn(|_| qi(random.between(-reach, reach)))
}

fn on_frame(origin: &[Rational; 3], frame: &[[i64; 3]; 3], local: [&Rational; 3]) -> [Rational; 3] {
    std::array::from_fn(|axis| {
        (0..3).fold(origin[axis].clone(), |sum, k| {
            sum + qi(frame[k][axis]) * local[k].clone()
        })
    })
}

fn fresh(
    name: String,
    source: &'static str,
    parts: (Curve, CaseSurface, FaceCopy, Truth),
) -> EdgeCase {
    let (curve, surface, face_copy, truth) = parts;
    EdgeCase {
        name,
        source,
        geometry: EdgeGeometry {
            curve,
            parameters: [Rational::zero(), Rational::one()],
            surface,
            face_copy,
        },
        truth,
        domain_exit: None,
        zero: None,
        zero_ceiling: Rational::zero(),
    }
}

/// Family A. A curved face copy of degree 2 or 3 on a random plane. The 3D
/// copy is its exact lift with extra knots inserted (Boehm) and one control of
/// the refined curve moved along the normal. The face copy has no knot near
/// the bump.
fn plane_lift(random: &mut Random) -> (Curve, CaseSurface, FaceCopy, Truth) {
    let degree = random.between(2, 3) as usize;
    let interior = random_knots(random, 2, 11);
    let flat_knots = clamped(degree, &interior, &Rational::zero(), &Rational::one());
    let flat: Vec<[Rational; 2]> = (0..flat_knots.len() - degree - 1)
        .map(|i| {
            [
                qi(8 * i as i64 + random.between(0, 5)),
                qi(random.between(-20, 20)),
            ]
        })
        .collect();
    let (origin, frame) = (offset(random), frame(random));
    let mut knots = flat_knots.clone();
    let mut controls: Vec<Vec<Rational>> = flat.iter().map(|uv| uv.to_vec()).collect();
    let cluster = cluster(random, degree, &interior);
    for knot in &cluster {
        insert_knot(degree, &mut knots, &mut controls, knot);
    }
    let index = knots
        .iter()
        .position(|knot| knot == &cluster[0])
        .expect("inserted");
    let mut lifts = vec![Rational::zero(); controls.len()];
    lifts[index] = bump(random);
    let lifted = controls
        .iter()
        .zip(&lifts)
        .map(|(uv, lift)| on_frame(&origin, &frame, [&uv[0], &uv[1], lift]))
        .collect();
    (
        Curve::BSpline(spline(degree, knots.clone(), lifted, None)),
        CaseSurface::Kernel(plane(origin, frame[2], frame[0])),
        FaceCopy::Native(pcurve(degree, flat_knots, flat, None)),
        spline_peak(&knots, degree, &lifts),
    )
}

/// Family B. A random spline patch and its own iso-curve as the 3D copy.
/// One control of the patch is then lifted, so the bump lives on the surface:
/// its size is the lift times the v basis value at the iso-parameter times
/// the u basis function along the edge.
fn surface_bump(random: &mut Random) -> (Curve, CaseSurface, FaceCopy, Truth) {
    let degree_v = random.between(1, 3) as usize;
    let mut interior = random_knots(random, 1, 6);
    let cluster = cluster(random, 3, &interior);
    interior.extend(cluster.iter().cloned());
    interior.sort();
    let knots_u = clamped(3, &interior, &Rational::zero(), &Rational::one());
    let knots_v = clamped(
        degree_v,
        &random_knots(random, 0, 3),
        &Rational::zero(),
        &Rational::one(),
    );
    let (us, vs) = (greville(&knots_u, 3), greville(&knots_v, degree_v));
    let origin = offset(random);
    let mut net: Vec<Vec<[Rational; 3]>> = us
        .iter()
        .map(|u| {
            vs.iter()
                .map(|v| {
                    [
                        &origin[0] + &(qi(100) * u.clone()),
                        &origin[1] + &(qi(20) * v.clone()),
                        &origin[2] + &ratio(random.between(-64, 64), 64),
                    ]
                })
                .collect()
        })
        .collect();
    let v = loop {
        let v = random.share(12);
        if !knots_v.contains(&v) {
            break v;
        }
    };
    let shares: Vec<_> = (0..vs.len())
        .map(|j| basis(&knots_v, degree_v, j, &v))
        .collect();
    let row = net
        .iter()
        .map(|column| {
            std::array::from_fn(|axis| {
                column
                    .iter()
                    .zip(&shares)
                    .fold(Rational::zero(), |sum, (p, share)| sum + share * &p[axis])
            })
        })
        .collect();
    let i = knots_u
        .iter()
        .position(|knot| knot == &cluster[0])
        .expect("cluster");
    let j = (0..vs.len())
        .max_by(|a, b| shares[*a].cmp(&shares[*b]))
        .expect("columns");
    let lift = bump(random);
    net[i][j][2] = &net[i][j][2] + &lift;
    let mut along = vec![Rational::zero(); us.len()];
    along[i] = &lift * &shares[j];
    (
        Curve::BSpline(spline(3, knots_u.clone(), row, None)),
        CaseSurface::Carrier(carrier(
            [3, degree_v],
            [knots_u.clone(), knots_v],
            net.into_iter().flatten().collect(),
        )),
        FaceCopy::Native(straight(
            [Rational::zero(), v.clone()],
            [Rational::one(), v],
        )),
        spline_peak(&knots_u, 3, &along),
    )
}

/// Family C. A single polynomial piece of degree 6 to 24 on a random plane,
/// straight in the plane, with one Bernstein function lifted. No knot marks
/// the bump.
fn bernstein_bump(random: &mut Random) -> (Curve, CaseSurface, FaceCopy, Truth) {
    let degree = random.between(6, 24) as usize;
    let index = random.between(1, degree as i64 - 1) as usize;
    let (origin, frame) = (offset(random), frame(random));
    let from = [qi(random.between(-50, 50)), qi(random.between(-50, 50))];
    let to = [
        &from[0] + &qi(random.between(20, 90)),
        qi(random.between(-50, 50)),
    ];
    let lift = bump(random);
    let controls = (0..=degree)
        .map(|i| {
            let share = ratio(i as i64, degree as i64);
            let uv: [Rational; 2] =
                std::array::from_fn(|axis| &from[axis] + &(&share * &(&to[axis] - &from[axis])));
            let z = if i == index {
                lift.clone()
            } else {
                Rational::zero()
            };
            on_frame(&origin, &frame, [&uv[0], &uv[1], &z])
        })
        .collect();
    (
        Curve::BSpline(spline(
            degree,
            clamped(degree, &[], &Rational::zero(), &Rational::one()),
            controls,
            None,
        )),
        CaseSurface::Kernel(plane(origin, frame[2], frame[0])),
        FaceCopy::Native(straight(from, to)),
        Truth::exact(
            &lift * &bernstein_peak(degree, index),
            ratio(index as i64, degree as i64),
        ),
    )
}

/// Family D. A circular arc on its cylinder. The face copy keeps the angle
/// and lifts the height by one basis function of a random knot vector.
fn cylinder_arc(random: &mut Random) -> (Curve, CaseSurface, FaceCopy, Truth) {
    let radius = qi(random.between(2, 50));
    let height = qi(random.between(-30, 30));
    let origin = offset(random);
    let frame = frame(random);
    let mut interior = random_knots(random, 1, 8);
    let cluster = cluster(random, 3, &interior);
    interior.extend(cluster.iter().cloned());
    interior.sort();
    let knots = clamped(3, &interior, &Rational::zero(), &Rational::one());
    let (start, sweep) = (
        ratio(random.between(-300, 0), 100),
        ratio(random.between(30, 300), 100),
    );
    let mut lifts = vec![Rational::zero(); knots.len() - 4];
    let index = knots
        .iter()
        .position(|knot| knot == &cluster[0])
        .expect("cluster");
    lifts[index] = bump(random);
    let controls = greville(&knots, 3)
        .iter()
        .zip(&lifts)
        .map(|(t, lift)| [&start + &(&sweep * t), &height + lift])
        .collect();
    let centre = on_frame(
        &origin,
        &frame,
        [&Rational::zero(), &Rational::zero(), &height],
    );
    // The circle's parameter is its angle, so the face copy's knots are
    // mapped from [0, 1] onto the angle range. The peak value does not change.
    let angle_knots: Vec<_> = knots.iter().map(|t| &start + &(&sweep * t)).collect();
    let mut truth = spline_peak(&knots, 3, &lifts);
    truth.at = &start + &(&sweep * &truth.at);
    (
        circle(centre, frame[2], &radius, frame[0]),
        CaseSurface::Kernel(cylinder(origin, frame[2], &radius, frame[0])),
        FaceCopy::Native(pcurve(3, angle_knots, controls, None)),
        truth,
    )
}

/// `per_family` new edges from each family, with the parameters of circular
/// edges set to their angle range.
pub(super) fn cases(seed: u64, per_family: usize) -> Vec<EdgeCase> {
    let mut random = Random::new(seed);
    let mut cases = Vec::new();
    for round in 0..per_family {
        cases.push(fresh(
            format!("fresh-plane-lift-{round}"),
            "seeded: Boehm-refined lift with a hidden normal bump",
            plane_lift(&mut random),
        ));
        cases.push(fresh(
            format!("fresh-surface-bump-{round}"),
            "seeded: hidden bump on one surface control",
            surface_bump(&mut random),
        ));
        cases.push(fresh(
            format!("fresh-bernstein-bump-{round}"),
            "seeded: hidden bump inside one polynomial piece",
            bernstein_bump(&mut random),
        ));
        let mut arc = fresh(
            format!("fresh-cylinder-arc-{round}"),
            "seeded: hidden height bump on a cylinder arc",
            cylinder_arc(&mut random),
        );
        if let FaceCopy::Native(curve) = &arc.geometry.face_copy {
            let [lo, hi] = curve.domain();
            arc.geometry.parameters = [lo.clone(), hi.clone()];
        }
        cases.push(arc);
    }
    cases
}

/// A twin of a case under an exact scale and shift. It has the same degrees,
/// knots, weights and surface kind, and its true gap is the scale times the
/// original. A checker whose point count depends only on shape plans the same
/// points for both; this is the handle for the second guard.
pub(super) fn twin(case: &EdgeCase, scale: &Rational, shift: &[Rational; 3]) -> EdgeCase {
    let moved = |p: &catallax_shapes::Point3| {
        let coordinates = [&p.x, &p.y, &p.z];
        point(std::array::from_fn(|axis| {
            scale * &coordinates[axis].as_rational().expect("rational fixture")
                + shift[axis].clone()
        }))
    };
    let scaled_sq = |r: &Real| Real::from(scale * scale * r.as_rational().expect("rational"));
    let respline =
        |curve: &BSplineCurve,
         map: &dyn Fn(&catallax_shapes::Point3) -> catallax_shapes::Point3| {
            BSplineCurve::new(
                curve.degree(),
                curve.knots().to_vec(),
                curve.controls().map(map).collect(),
                curve.weights().cloned().collect(),
                &mut Decider::default(),
            )
            .expect("twin curve")
        };
    let curve = match &case.geometry.curve {
        Curve::Line(line) => Curve::Line(catallax_shapes::Line {
            origin: moved(&line.origin),
            dir: line.dir.scaled(&Real::from(scale.clone())),
        }),
        Curve::Circle(circle) => Curve::Circle(catallax_shapes::Circle {
            center: moved(&circle.center),
            radius_sq: scaled_sq(&circle.radius_sq),
            ..circle.clone()
        }),
        Curve::BSpline(curve) => Curve::BSpline(respline(curve, &moved)),
        _ => unreachable!("fixtures use lines, circles and splines"),
    };
    // Plane parameters are lengths and scale with the body. A cylinder keeps
    // its angle and scales its height. Carrier parameters do not change.
    let factors = match &case.geometry.surface {
        CaseSurface::Kernel(Surface::Plane(_)) => [scale.clone(), scale.clone()],
        CaseSurface::Kernel(Surface::Cylinder(_)) => [Rational::one(), scale.clone()],
        _ => [Rational::one(), Rational::one()],
    };
    let surface = match &case.geometry.surface {
        CaseSurface::Kernel(Surface::Plane(plane)) => {
            CaseSurface::Kernel(Surface::Plane(catallax_shapes::Plane {
                origin: moved(&plane.origin),
                ..plane.clone()
            }))
        }
        CaseSurface::Kernel(Surface::Cylinder(cylinder)) => {
            CaseSurface::Kernel(Surface::Cylinder(catallax_shapes::Cylinder {
                origin: moved(&cylinder.origin),
                radius_sq: scaled_sq(&cylinder.radius_sq),
                ..cylinder.clone()
            }))
        }
        CaseSurface::Carrier(carrier) => CaseSurface::Carrier(
            BSplineSurface::new(
                carrier.degrees(),
                carrier.knots().map(|knots| knots.to_vec()),
                carrier.control_count(),
                carrier.controls().iter().map(moved).collect(),
                carrier.weights().to_vec(),
                &mut Decider::default(),
            )
            .expect("twin surface"),
        ),
        CaseSurface::Kernel(_) => unreachable!("fixtures use planes and cylinders"),
    };
    let face_copy = match &case.geometry.face_copy {
        FaceCopy::Native(curve) => FaceCopy::Native(respline(curve, &|p| {
            point([
                &factors[0] * &p.x.as_rational().expect("rational"),
                &factors[1] * &p.y.as_rational().expect("rational"),
                Rational::zero(),
            ])
        })),
        FaceCopy::Linear(ends) => FaceCopy::Linear(Box::new(
            ends.each_ref()
                .map(|uv| std::array::from_fn(|axis| &factors[axis] * &uv[axis])),
        )),
    };
    EdgeCase {
        name: format!("{}-twin", case.name),
        geometry: EdgeGeometry {
            curve,
            parameters: case.geometry.parameters.clone(),
            surface,
            face_copy,
        },
        truth: Truth {
            lo: &case.truth.lo * scale,
            hi: &case.truth.hi * scale,
            at: case.truth.at.clone(),
        },
        zero_ceiling: &case.zero_ceiling * scale,
        ..case.clone()
    }
}

#[test]
fn a_seed_replays_the_same_edges_and_their_truths_hold() {
    let seed = seed();
    println!("fresh edge seed={seed}");
    let first = cases(seed, 3);
    let again = cases(seed, 3);
    assert_eq!(first.len(), 12);
    for (a, b) in first.iter().zip(&again) {
        assert_eq!(a.name, b.name);
        assert_eq!(a.truth.lo, b.truth.lo, "seed {seed}");
        assert_eq!(a.truth.at, b.truth.at, "seed {seed}");
    }
    for case in &first {
        assert!(case.truth.lo.sign() == Sign::Positive, "seed {seed}");
        truth_is_attained_and_bounds(case, 96);
    }
    assert_ne!(
        cases(seed ^ 1, 1)[0].truth.at,
        first[0].truth.at,
        "a different seed moves the bump"
    );
}

#[test]
fn a_twin_keeps_the_shape_and_scales_the_truth() {
    let seed = seed();
    println!("twin seed={seed}");
    let scale = ratio(7, 3);
    let shift = [qi(11), qi(-5), qi(2)];
    // The domain-exit case is left out: its oracle check stops at the exit.
    let originals = cases(seed, 1).into_iter().chain(
        super::fixed::cases()
            .into_iter()
            .filter(|case| case.domain_exit.is_none()),
    );
    for case in originals {
        let twin = twin(&case, &scale, &shift);
        assert_eq!(twin.truth.lo, &case.truth.lo * &scale);
        truth_is_attained_and_bounds(&twin, 24);
    }
}
