//! Finite cone bands and a common cylinder frame for numerical face consumers.

use super::*;

/// Searches a finite positive cone-height band against a finite cylinder band.
/// Heights and cylinder coordinates are signed millimetres along their unit
/// axes from the original analytic origins. In either operand order, the same
/// cylinder and interval retain the same source and UV frame across cone faces.
///
/// All sixteen quarter pairs are searched. Every reported 3D/UV definition,
/// including duplicate boundary contacts, is retained without parameter changes.
/// Samples check requested rectangle ownership and physical incidence at 1e-7 mm.
/// This is neither finite Boolean topology nor a complete branch-count proof.
pub fn numerical_native_finite_intersection(
    surfaces: [&Surface; 2],
    cone_heights_mm: [Rational; 2],
    cylinder_parameters_mm: [Rational; 2],
    dec: &mut Decider,
) -> Result<Option<NumericalNativeIntersection>, Refusal> {
    let run =
        |dec: &mut Decider| build_finite(surfaces, cone_heights_mm, cylinder_parameters_mm, dec);
    if operation_budget_active() {
        return run(dec);
    }
    let (result, work) = with_operation_budget(20_000_000, || run(dec));
    if work.exhausted {
        return Err(Refusal::new(
            "numbers.operation-budget",
            "finite numerical native intersection exhausted its operation allowance",
        ));
    }
    result
}

fn build_finite(
    surfaces: [&Surface; 2],
    heights: [Rational; 2],
    cylinder_parameters: [Rational; 2],
    dec: &mut Decider,
) -> Result<Option<NumericalNativeIntersection>, Refusal> {
    let started = catallax_numbers::current_operation_work().used;
    let (cone, cylinder) = match surfaces {
        [Surface::Cone(cone), Surface::Cylinder(cylinder)]
        | [Surface::Cylinder(cylinder), Surface::Cone(cone)] => (cone, cylinder),
        _ => return Ok(None),
    };
    dec.charge_work(1)?;
    if dec.policy().max_precision_bits < BITS {
        return Err(bad("native numerical windows require 128-bit enclosures"));
    }
    if heights[0].sign() != Sign::Positive || heights[0] >= heights[1] {
        return Err(bad(
            "finite cone heights must be increasing and strictly positive",
        ));
    }
    if cylinder_parameters[0] >= cylinder_parameters[1] {
        return Err(bad(
            "finite cylinder coordinates must be strictly increasing",
        ));
    }
    if enclosure(&cone.slope, dec)?.lo().sign() != Sign::Positive
        || enclosure(&cylinder.radius(), dec)?.lo().sign() != Sign::Positive
    {
        return Err(bad(
            "finite native sources require positive slope and radius",
        ));
    }
    let cone_frame = Frame::new(cone.apex.clone(), &cone.axis, &cone.ref_dir, dec)?;
    let cylinder_frame = Frame::new(
        cylinder.origin.clone(),
        &cylinder.axis,
        &cylinder.ref_dir,
        dec,
    )?;
    let checks = vec![NumericalCheck::Passed {
        name: "native.retained-requested-finite-bands".into(),
        scope: "cone-cylinder-original-origins".into(),
    }];
    let windows = assemble(
        surfaces,
        heights,
        cylinder_parameters,
        cone_frame,
        cylinder_frame,
        checks,
        dec,
    )?
    .expect("cone/cylinder family checked above");
    let after_windows = catallax_numbers::current_operation_work().used;
    checked::search(windows, true, started, after_windows, dec).map(Some)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn axial(offset: Rational) -> [Surface; 2] {
        [
            Surface::Cone(Cone {
                apex: Point3::new(s(0), s(-11), s(0)),
                axis: Vector3::new(s(0), s(1), s(0)),
                slope: s(1),
                ref_dir: Vector3::new(s(1), s(0), s(0)),
            }),
            Surface::Cylinder(Cylinder {
                origin: Point3::new(offset.into(), s(12), s(0)),
                axis: Vector3::new(s(0), s(-1), s(0)),
                radius_sq: s(1),
                ref_dir: Vector3::new(s(1), s(0), s(0)),
            }),
        ]
    }

    #[test]
    fn finite_bands_distinguish_separation_touch_and_crossing_in_both_orders() {
        use catallax_trace::surface_curves::CurveKind;
        for (offset, expected) in [
            (Rational::parse("9.999").unwrap(), 0),
            (r(10), 4),
            (Rational::parse("10.01").unwrap(), 4),
        ] {
            for reverse in [false, true] {
                let surfaces = axial(offset.clone());
                let ordered = if reverse {
                    [&surfaces[1], &surfaces[0]]
                } else {
                    [&surfaces[0], &surfaces[1]]
                };
                let (result, work) = with_operation_budget(20_000_000, || {
                    numerical_native_finite_intersection(
                        ordered,
                        [r(11), r(12)],
                        [r(10), r(13)],
                        &mut Decider::default(),
                    )
                });
                let result = result.unwrap().unwrap();
                assert!(!work.exhausted, "{work:?}");
                assert_eq!(
                    result.pieces().len(),
                    expected,
                    "offset={offset}, reverse={reverse}"
                );
                let crossings = result
                    .pieces()
                    .iter()
                    .filter(|p| p.kind() == CurveKind::Crossing)
                    .count();
                assert_eq!(crossings, if offset > r(10) { 2 } else { 0 });
                assert_eq!(result.windows().cylinder_parameters(), &[r(10), r(13)]);
                assert!(result.windows().cylinder_origin().x.as_rational() == Some(offset.clone()));
                let stages = result.work_stages();
                assert_eq!(
                    stages.windows + stages.setup + stages.tracing + stages.fit_checks,
                    work.used
                );
                for piece in result.pieces() {
                    assert_eq!(piece.curve().domain(), piece.uv()[0].domain());
                    assert_eq!(piece.curve().domain(), piece.uv()[1].domain());
                }
            }
        }
    }

    #[test]
    fn finite_windows_do_not_require_global_parallel_generator_or_apex_separation() {
        for (offset, axis) in [
            (r(1), Vector3::new(s(0), s(-1), s(0))),
            (r(10), Vector3::new(s(1), s(-1), s(0))),
        ] {
            let mut surfaces = axial(offset);
            let Surface::Cylinder(cylinder) = &mut surfaces[1] else {
                unreachable!()
            };
            cylinder.axis = axis;
            let (result, work) = with_operation_budget(20_000_000, || {
                numerical_native_finite_intersection(
                    [&surfaces[0], &surfaces[1]],
                    [r(11), r(12)],
                    [r(10), r(13)],
                    &mut Decider::default(),
                )
            });
            assert!(result.unwrap().is_some());
            assert!(!work.exhausted, "{work:?}");
        }
    }

    #[test]
    fn invalid_bands_and_caller_exhaustion_do_not_fall_back_to_natural_domains() {
        let surfaces = axial(r(10));
        for (cone, cylinder) in [
            ([r(0), r(1)], [r(10), r(13)]),
            ([r(12), r(11)], [r(10), r(13)]),
            ([r(11), r(12)], [r(13), r(10)]),
            ([r(11), r(12)], [r(10), r(10)]),
        ] {
            assert!(
                numerical_native_finite_intersection(
                    [&surfaces[0], &surfaces[1]],
                    cone,
                    cylinder,
                    &mut Decider::default(),
                )
                .is_err()
            );
        }
        let (result, work) = with_operation_budget(0, || {
            numerical_native_finite_intersection(
                [&surfaces[0], &surfaces[1]],
                [r(11), r(12)],
                [r(10), r(13)],
                &mut Decider::default(),
            )
        });
        assert!(result.is_err());
        assert!(work.exhausted);
    }
}
