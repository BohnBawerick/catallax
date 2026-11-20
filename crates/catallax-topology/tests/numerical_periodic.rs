//! Complete circular rims need an explicit paired seam before a cylinder chart closes.
#[path = "common/periodic.rs"]
mod fixture;

use catallax_numbers::Decider;
use catallax_topology::numerical::{self as n, NumericalOptions};

#[test]
fn complete_rims_and_one_paired_seam_admit_a_closed_bore() {
    let body = n::build(
        fixture::bore(),
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let counts = body.counts();
    assert_eq!((counts.faces, counts.edges, counts.vertices), (7, 15, 10));
    assert_eq!((counts.shells, counts.lumps), (1, 1));
    assert_eq!(body.evidence().tolerance_mm(), &fixture::q("1/100"));
}

#[test]
fn closed_partial_arcs_extra_turns_and_implicit_closure_refuse() {
    for kind in 0..3 {
        let mut definition = fixture::bore();
        let rim = definition
            .edges
            .iter()
            .position(|e| e.key == "rim-0")
            .unwrap();
        match kind {
            0 => definition.edges[rim].parameters[1] = fixture::q("1"),
            1 => {
                definition.edges[rim].parameters[1] =
                    &definition.edges[rim].parameters[1] * &fixture::q("2")
            }
            _ => {
                let mut coincident = definition.vertices[8].clone();
                coincident.key = "unsewn-coincident-anchor".into();
                definition.edges[rim].vertices[1] = definition.vertices.len();
                definition.vertices.push(coincident);
            }
        }
        assert!(
            n::build(
                definition,
                NumericalOptions::default(),
                &mut Decider::default()
            )
            .is_err()
        );
    }
}

#[test]
fn separate_seam_edges_wrong_winding_and_double_period_cells_refuse() {
    for kind in 0..3 {
        let mut definition = fixture::bore();
        let wall = definition.faces.len() - 1;
        match kind {
            0 => {
                let mut duplicate = definition.edges.last().unwrap().clone();
                duplicate.key = "unshared-seam".into();
                definition.faces[wall].loops[0][3].edge = definition.edges.len();
                definition.edges.push(duplicate);
            }
            1 => definition.faces[wall].loops[0][3].sense = catallax_topology::Sense::Same,
            _ => {
                for trim in &mut definition.faces[wall].loops[0] {
                    for uv in &mut trim.uv {
                        uv[0] = &uv[0] * &fixture::q("2");
                    }
                }
            }
        }
        assert!(
            n::build(
                definition,
                NumericalOptions::default(),
                &mut Decider::default()
            )
            .is_err()
        );
    }
}

#[test]
fn translated_seam_branches_may_slope_within_the_numerical_tolerance() {
    let mut definition = fixture::bore();
    for trim in definition
        .faces
        .last_mut()
        .unwrap()
        .loops
        .iter_mut()
        .flatten()
    {
        for uv in &mut trim.uv {
            uv[0] = &uv[0] + &(&uv[1] * &fixture::q("1/100000000"));
        }
    }
    n::build(
        definition,
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
}

#[test]
fn a_shared_seam_cannot_fall_back_to_a_near_turn_chart() {
    let mut definition = fixture::bore();
    for trim in definition
        .faces
        .last_mut()
        .unwrap()
        .loops
        .iter_mut()
        .flatten()
    {
        for uv in &mut trim.uv {
            if uv[0] > fixture::q("1") {
                uv[0] = &uv[0] - &fixture::q("1/1000");
            }
        }
    }
    let error = n::build(
        definition,
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code, "numerical.invalid-trim");
    assert!(error.detail.contains("paired seam"));
}

#[test]
fn continuous_circle_clearance_respects_the_requested_tolerance() {
    let definition = fixture::bore_with_radius("399/200");
    assert!(
        n::build(
            definition.clone(),
            NumericalOptions::default(),
            &mut Decider::default()
        )
        .is_err()
    );
    let body = n::build(
        definition,
        NumericalOptions {
            tolerance_mm: fixture::q("1/1000"),
        },
        &mut Decider::default(),
    )
    .unwrap();
    assert!(body.evidence().checks().iter().any(|check| matches!(check,
        catallax_numbers::NumericalCheck::Passed { name,.. } if name=="trim.continuous-circle-hole-clearance")));
}

#[test]
fn multiple_bores_require_pairwise_clearance_under_the_requested_tolerance() {
    let body = n::build(
        fixture::bores(&["-0.5", "0.5"], "0.4"),
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let counts = body.counts();
    assert_eq!((counts.faces, counts.edges, counts.vertices), (8, 18, 12));
    assert_eq!(body.definition().faces[0].loops.len(), 3);
    for (tolerance, passes) in [("0.01", false), ("0.001", true)] {
        let result = n::build(
            fixture::bores(&["-0.4025", "0.4025"], "0.4"),
            NumericalOptions {
                tolerance_mm: fixture::q(tolerance),
            },
            &mut Decider::default(),
        );
        if passes {
            result.unwrap();
        } else {
            assert_eq!(result.unwrap_err().code, "numerical.invalid-trim");
        }
    }
    for centers in [["0", "0"], ["-0.2", "0.2"], ["0", "2"]] {
        let error = n::build(
            fixture::bores(&centers, "0.4"),
            NumericalOptions::default(),
            &mut Decider::default(),
        )
        .unwrap_err();
        assert_eq!(error.code, "numerical.invalid-trim");
    }
}

#[test]
fn a_thin_annulus_uses_its_actual_rims_not_crossing_inscribed_chords() {
    let mut definition = fixture::annulus("1.03");
    let reference =
        catallax_shapes::RVector3::new(fixture::q("3/5"), fixture::q("4/5"), fixture::q("0"))
            .to_real();
    for vertex in &mut definition.vertices[..2] {
        vertex.position[0] = fixture::q("3/5");
        vertex.position[1] = fixture::q("4/5");
    }
    for edge in &mut definition.edges[..2] {
        let catallax_shapes::Curve::Circle(circle) = &mut edge.curve else {
            unreachable!()
        };
        circle.ref_dir = reference.clone();
    }
    let catallax_shapes::Curve::Line(seam) = &mut definition.edges[2].curve else {
        unreachable!()
    };
    seam.origin =
        catallax_shapes::RPoint3::new(fixture::q("3/5"), fixture::q("4/5"), fixture::q("0"))
            .to_real();
    for cap in &mut definition.faces[..2] {
        cap.loops[1][0].uv = [
            [fixture::q("3/5"), fixture::q("4/5")],
            [fixture::q("3/5"), fixture::q("4/5")],
        ];
    }
    let catallax_shapes::Surface::Cylinder(wall) = &mut definition.faces[2].surface else {
        unreachable!()
    };
    wall.ref_dir = reference;
    let body = n::build(
        definition,
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let counts = body.counts();
    assert_eq!((counts.faces, counts.edges, counts.vertices), (4, 6, 4));
    for (tolerance, passes) in [("0.01", false), ("0.001", true)] {
        let result = n::build(
            fixture::annulus("1.015"),
            NumericalOptions {
                tolerance_mm: fixture::q(tolerance),
            },
            &mut Decider::default(),
        );
        if passes {
            result.unwrap();
        } else {
            assert_eq!(result.unwrap_err().code, "numerical.invalid-trim");
        }
    }
}

#[test]
fn rational_patch_native_circles_retain_only_sampled_boundary_evidence() {
    use catallax_numbers::{NumericalCheck, Real};
    use catallax_shapes::{
        BSplineCurve, BSplineSurface, Point3, RationalPatch, SplineSupport, Surface,
    };
    let mut definition = fixture::annulus("1.03");
    let mut dec = Decider::default();
    let patch = BSplineSurface::new(
        [1, 1],
        [
            vec![
                fixture::q("-2"),
                fixture::q("-2"),
                fixture::q("2"),
                fixture::q("2"),
            ],
            vec![
                fixture::q("-2"),
                fixture::q("-2"),
                fixture::q("2"),
                fixture::q("2"),
            ],
        ],
        [2, 2],
        [[-2, -2, 0], [-2, 2, 0], [2, -2, 0], [2, 2, 0]]
            .into_iter()
            .map(|[x, y, z]| Point3::new(x.into(), y.into(), z.into()))
            .collect(),
        vec![Real::one(); 4],
        &mut dec,
    )
    .unwrap();
    definition.faces[0].surface = Surface::BSpline(SplineSupport::RationalPatch(
        RationalPatch::new(patch, &mut dec).unwrap(),
    ));
    let period = catallax_numbers::Ball::pi(128).mid() * fixture::q("2");
    let quarter = &period / &fixture::q("4");
    let knots = std::iter::repeat_n(fixture::q("0"), 4)
        .chain((1..4).flat_map(|i| {
            std::iter::repeat_n(&quarter * &catallax_numbers::Rational::from_i64(i), 3)
        }))
        .chain(std::iter::repeat_n(period, 4))
        .collect::<Vec<_>>();
    const K: &str = "0.5522847498307936";
    const M: &str = "-0.5522847498307936";
    const CONTROLS: [[&str; 2]; 13] = [
        ["1", "0"],
        ["1", K],
        [K, "1"],
        ["0", "1"],
        [M, "1"],
        ["-1", K],
        ["-1", "0"],
        ["-1", M],
        [M, "-1"],
        ["0", "-1"],
        [K, "-1"],
        ["1", M],
        ["1", "0"],
    ];
    let mut bindings = Vec::new();
    for (boundary, radius) in [(0, "1.03"), (1, "1")] {
        let radius = fixture::q(radius);
        let controls = CONTROLS
            .into_iter()
            .map(|[u, v]| {
                Point3::new(
                    (&fixture::q(u) * &radius).into(),
                    (&fixture::q(v) * &radius).into(),
                    0.into(),
                )
            })
            .collect();
        bindings.push(n::NumericalNativeTrim {
            tag: catallax_topology::numerical::CopyTag::Independent,
            face: definition.faces[0].key.clone(),
            boundary,
            use_index: 0,
            edge: definition.edges[definition.faces[0].loops[boundary][0].edge]
                .key
                .clone(),
            curve: BSplineCurve::new(3, knots.clone(), controls, vec![Real::one(); 13], &mut dec)
                .unwrap(),
        });
    }
    let body =
        n::build_with_native_trims(definition, bindings, NumericalOptions::default(), &mut dec)
            .unwrap();
    let passed = |name: &str, scope: &str| {
        body.evidence().checks().iter().any(|check| {
            matches!(check, NumericalCheck::Passed { name: actual, scope: actual_scope }
            if actual == name && actual_scope == scope)
        })
    };
    assert!(passed("trim.sampled-hole-contained-disjoint", "f0"));
    for name in [
        "trim.sampled-boundaries-with-continuous-circle-clearance",
        "trim.continuous-circle-hole-clearance",
        "trim.continuous-circle-outer-clearance",
    ] {
        assert!(!passed(name, "f0"));
    }
    assert!(passed("trim.continuous-circle-hole-clearance", "f1"));
    assert!(passed("trim.continuous-circle-outer-clearance", "f1"));
}
