//! A regular rational sphere window; its omitted pole is proved off the cone.

use super::*;

fn source(
    radius: &Rational,
    bound: &Rational,
    dec: &mut Decider,
) -> Result<BSplineSurface, Refusal> {
    let half = bound / &r(2);
    let positions = [-bound, -&half, r(0), half, bound.clone()];
    let square = bound * bound;
    let squares = [square.clone(), r(0), r(0), r(0), square];
    let twice_radius = &r(2) * radius;
    let mut controls = Vec::with_capacity(25);
    let mut weights = Vec::with_capacity(25);
    // Splitting u and v at zero makes every Bernstein weight positive.
    // X=2ru, Y=2rv, Z=r(1-u²-v²), W=1+u²+v². The natural
    // Du x Dv = (4r/W²)P points outwards. Only -r*z is omitted.
    for u in 0..5 {
        for v in 0..5 {
            let sum = &squares[u] + &squares[v];
            let weight = &r(1) + &sum;
            controls.push(Point3::new(
                Real::from(&(&twice_radius * &positions[u]) / &weight),
                Real::from(&(&twice_radius * &positions[v]) / &weight),
                Real::from(&(radius * &(&r(1) - &sum)) / &weight),
            ));
            weights.push(Real::from(weight));
        }
    }
    let knots = Vec::from([
        -bound,
        -bound,
        -bound,
        r(0),
        r(0),
        bound.clone(),
        bound.clone(),
        bound.clone(),
    ]);
    BSplineSurface::new(
        [2, 2],
        [knots.clone(), knots],
        [5, 5],
        controls,
        weights,
        dec,
    )
}

fn pole_window(
    delta: &RVector3,
    cone_axis: &RVector3,
    slope: &Rational,
    radius: &Rational,
    dec: &mut Decider,
) -> Result<(RVector3, Rational), Refusal> {
    const POLES: [[(i64, i64); 3]; 9] = [
        [(1, 1), (0, 1), (0, 1)],
        [(-1, 1), (0, 1), (0, 1)],
        [(0, 1), (1, 1), (0, 1)],
        [(0, 1), (-1, 1), (0, 1)],
        [(0, 1), (0, 1), (1, 1)],
        [(0, 1), (0, 1), (-1, 1)],
        [(3, 5), (4, 5), (0, 1)],
        [(3, 5), (0, 1), (4, 5)],
        [(0, 1), (3, 5), (4, 5)],
    ];
    let base = &delta.norm_squared() + &(radius * radius);
    let twice_radius = &r(2) * radius;
    let cone_factor = &r(1) + &(slope * slope);
    let height = delta.dot(cone_axis);
    let mut best: Option<(RVector3, Rational)> = None;
    for coordinates in POLES {
        dec.check_deadline()?;
        let [x, y, z] = coordinates
            .map(|(n, d)| Rational::from_ratio(n, d).expect("positive pole denominator"));
        let direction = RVector3 { x, y, z };
        if direction.norm_squared() != r(1) {
            return Err(refused(
                "the spherical exclusion pole needs an exact unit direction",
            ));
        }
        let axial_height = &height + &(radius * &direction.dot(cone_axis));
        let linear = &twice_radius * &delta.dot(&direction);
        let axial = &cone_factor * &(&axial_height * &axial_height);
        let value = (&(&base + &linear) - &axial).abs();
        if best.as_ref().is_none_or(|(_, previous)| &value > previous) {
            best = Some((direction, value));
        }
    }
    let (pole, value) = best.expect("the fixed pole set is nonempty");
    if value.is_zero() {
        return Err(Refusal::new(
            "intersect.trace-native-pole",
            "no fixed spherical pole has proved strict cone clearance",
        ));
    }
    // F=|P-A|²-(1+s²)((P-A).d)². On the sphere's convex ball,
    // |grad F| <= 2 (|C-A|+r) (2+s²). Thus F(P)=0 implies
    // |P-pole| >= |F(pole)|/L. For this stereographic chart,
    // |P-pole|=2r/sqrt(1+u²+v²), so |u|,|v|<2rL/|F(pole)|+1.
    let displacement = &(&delta.x.abs() + &delta.y.abs()) + &delta.z.abs();
    let distance = &displacement + radius;
    let gradient = &(&(&r(2) * &distance) * &(&r(2) + &(slope * slope))) + &r(1);
    let bound = &(&(&twice_radius * &gradient) / &value) + &r(1);
    Ok((pole, bound))
}

pub(in crate::pair) fn cone_sphere(
    cone: &Cone,
    sphere: &Sphere,
    sigma: Real,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let mut stage = "bounded native window";
    let result = (|| -> Result<Plan, Refusal> {
        let slope = rational(&cone.slope, dec)?;
        let radius = rational(&sphere.radius, dec)?;
        if slope.sign() != Sign::Positive || radius.sign() != Sign::Positive {
            return Err(refused(
                "the complete cone/sphere window needs positive slope and radius",
            ));
        }
        let direction = axis(&cone.axis, dec)?;
        let apex = point(&cone.apex, dec)?;
        let center = point(&sphere.center, dec)?;
        let delta = center.minus(&apex);
        let height = delta.dot(&direction);
        let upper = &height + &radius;
        if upper.sign() == Sign::Negative {
            return Ok(Plan::disjoint());
        }
        if upper.is_zero() {
            return Err(Refusal::new(
                "intersect.trace-native-apex",
                "the only possible retained cone height reaches its singular apex",
            ));
        }
        let lower = &height - &radius;
        let lower = if lower.sign() == Sign::Positive {
            &lower / &r(2)
        } else {
            let gap = (&delta.norm_squared() - &(&radius * &radius)).abs();
            if gap.is_zero() {
                return Err(Refusal::new(
                    "intersect.trace-native-apex",
                    "the cone apex lies on the sphere and needs a singular native chart proof",
                ));
            }
            let displacement = &(&delta.x.abs() + &delta.y.abs()) + &delta.z.abs();
            &gap / &(&r(2) * &(&(&(&displacement + &radius) + &r(1)) * &(&r(1) + &slope)))
        };
        if lower >= upper {
            return Ok(Plan::disjoint());
        }
        stage = "strict spherical pole exclusion and window bound";
        let (pole, bound) = pole_window(&delta, &direction, &slope, &radius, dec)?;
        let north = RVector3 {
            x: -&pole.x,
            y: -&pole.y,
            z: -&pole.z,
        };
        let sphere_map = frame(&north, &center)?;
        let cone_map = frame(&direction, &apex)?;
        let cone_carrier = SectionSurface::Cone(Cone {
            apex: zero_point(),
            axis: z_axis(),
            slope: Real::from(slope.clone()),
            ref_dir: Vector3::new(s(1), s(0), s(0)),
        });
        let sphere_carrier = SectionSurface::Sphere(Sphere {
            center: zero_point(),
            axis: z_axis(),
            radius: Real::from(radius.clone()),
            ref_dir: Vector3::new(s(1), s(0), s(0)),
        });
        stage = "exact oriented native-carrier equivalence";
        same_oriented_carrier(
            &cone_map.surface(&cone_carrier.surface()),
            &Surface::Cone(cone.clone()),
            dec,
        )?;
        same_oriented_carrier(
            &sphere_map.surface(&sphere_carrier.surface()),
            &Surface::Sphere(sphere.clone()),
            dec,
        )?;
        let cone_source = cone_source(&slope, [&lower, &(&upper + &r(1))], dec)?;
        let sphere_source = source(&radius, &bound, dec)?;
        complete_pair(
            [
                NativeWindow {
                    source: cone_source,
                    periodic: [true, false],
                    carrier: cone_carrier,
                    map: cone_map,
                    admission_stage: "complete canonical cone admission",
                },
                NativeWindow {
                    source: sphere_source,
                    periodic: [false, false],
                    carrier: sphere_carrier,
                    map: sphere_map,
                    admission_stage: "complete pole-excluded spherical window admission",
                },
            ],
            &sigma,
            &mut stage,
            dec,
        )
    })();
    result.map_err(|mut error| {
        error.detail = format!("traced cone/sphere {stage}: {}", error.detail);
        error
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_spherical_window_has_exact_incidence_and_outward_orientation() {
        let mut dec = Decider::default();
        let source = source(&r(2), &Rational::from_ratio(1, 4).unwrap(), &mut dec).unwrap();
        let atlases =
            RationalAtlas::partition_with_polynomial_bounds(source, [false, false], &mut dec)
                .unwrap();
        let carrier = SectionSurface::Sphere(Sphere {
            center: zero_point(),
            axis: z_axis(),
            radius: s(2),
            ref_dir: Vector3::new(s(1), s(0), s(0)),
        });
        assert_eq!(atlases.len(), 1);
        atlases[0]
            .certify_analytic_carrier(&carrier, &mut dec)
            .unwrap();
        let patch = atlases[0].charts()[0].patch();
        let at = std::array::from_fn(|axis| {
            let [lo, hi] = &patch.domain()[axis];
            Real::from(&(lo + hi) / &r(2))
        });
        let point = patch.point_at(&at, &mut dec).unwrap();
        let [du, dv] = patch.derivatives_at(&at, &mut dec).unwrap();
        assert_eq!(
            point.minus(&zero_point()).norm_squared().rational_ref(),
            Some(&r(4))
        );
        assert_eq!(
            sign(&du.cross(&dv).dot(&point.minus(&zero_point())), &mut dec).unwrap(),
            Sign::Positive
        );
    }

    #[test]
    fn a_sphere_through_the_apex_does_not_gain_a_regular_native_cover() {
        let cone = Cone {
            apex: zero_point(),
            axis: z_axis(),
            slope: Real::one(),
            ref_dir: Vector3::new(s(1), s(0), s(0)),
        };
        let sphere = Sphere {
            center: Point3::new(s(1), s(0), s(0)),
            axis: z_axis(),
            radius: Real::one(),
            ref_dir: Vector3::new(s(1), s(0), s(0)),
        };
        let refusal =
            cone_sphere(&cone, &sphere, Real::one(), &mut Decider::default()).unwrap_err();
        assert_eq!(refusal.code, "intersect.trace-native-apex");
    }
}
