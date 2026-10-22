//! Numbered exact cases, refusals, empty identities and mark provenance.
mod common;
use common::*;

macro_rules! pair_test {
    ($(#[$attr:meta])* $name:ident,$number:literal) => {
        blind_test!($(#[$attr])* $name, {
            for &id in IDS
                .iter()
                .filter(|id| id[..2].parse::<usize>().unwrap() == $number)
            {
                let c = case(id);
                for (op, reverse, want) in expectations(id) {
                    run_case(&c, op, reverse, Some(&want), 8);
                }
            }
        });
    };
}
pair_test!(case_01_general_blocks, 1);
pair_test!(case_02_face_contact, 2);
pair_test!(case_03_partial_face_contact, 3);
pair_test!(case_04_edge_contact, 4);
pair_test!(case_05_corner_contact, 5);
pair_test!(case_06_through_blocks, 6);
pair_test!(case_07_equal_size_blocks, 7);
pair_test!(case_08_identical_blocks, 8);
pair_test!(case_09_near_miss_blocks, 9);
pair_test!(case_10_through_cylinder, 10);
pair_test!(case_11_blind_hole, 11);
pair_test!(case_12_tangent_line, 12);
pair_test!(case_13_tangent_near_misses, 13);
pair_test!(case_14_inscribed_cylinder, 14);
pair_test!(case_15_cylinder_seam, 15);
pair_test!(case_16_cylinder_void, 16);
pair_test!(case_17_corner_sphere, 17);
pair_test!(case_18_face_sphere, 18);
pair_test!(case_19_external_tangent_sphere, 19);
pair_test!(case_20_internal_tangent_sphere, 20);
pair_test!(case_21_inscribed_sphere, 21);
pair_test!(case_22_sphere_seam_and_poles, 22);
pair_test!(case_23_crossing_cylinders, 23);
pair_test!(case_24_steinmetz, 24);
pair_test!(case_25_parallel_overlap, 25);
pair_test!(case_26_parallel_tangent, 26);
pair_test!(case_27_coaxial_cylinders, 27);
pair_test!(case_28_cylinder_node, 28);
pair_test!(case_29_point_tangent_cylinders, 29);
pair_test!(case_30_coaxial_cylinder_sphere, 30);
pair_test!(case_31_equator_contact, 31);
pair_test!(case_32_viviani_node, 32);
pair_test!(case_33_sphere_bite, 33);
pair_test!(case_34_point_tangent_cylinder_sphere, 34);
pair_test!(case_35_cylinder_seam_sphere, 35);
pair_test!(case_36_crossing_spheres, 36);
pair_test!(case_37_external_tangent_spheres, 37);
pair_test!(case_38_internal_tangent_spheres, 38);
pair_test!(case_39_concentric_spheres, 39);
pair_test!(case_40_identical_spheres, 40);

blind_test!(case_45_empty_set_identities_preserve_inputs, {
    for &id in IDS {
        for (op, reverse, keeps_a) in [
            (BoolOp::Unite, false, true),
            (BoolOp::Unite, true, true),
            (BoolOp::Intersect, false, false),
            (BoolOp::Intersect, true, false),
            (BoolOp::Subtract, false, true),
            (BoolOp::Subtract, true, false),
        ] {
            let mut m = Model::new();
            let a = case(id).a.build(&mut m);
            let e = empty(&mut m);
            let before = topology_snapshot(&m, &[a, e]);
            let want = if keeps_a {
                counts(&m, a)
            } else {
                Counts::empty()
            };
            let mut dec = Decider::default();
            let (x, y) = if reverse { (e, a) } else { (a, e) };
            let answer = boolean(&mut m, x, y, op, &mut dec).unwrap();
            assert_ne!(answer.body, a);
            assert_ne!(answer.body, e);
            assert_eq!(before, topology_snapshot(&m, &[a, e]));
            check_counts(&m, answer.body, &want);
            assert!(answer.marks.is_empty());
            if !answer.report.stages.is_empty() {
                check_report(&answer, &DecisionStats::default(), dec.stats());
            }
            for point in [
                p("0", "0", "0"),
                p("1/2", "1/2", "1/2"),
                p("100", "-10", "7"),
            ] {
                let expected = if keeps_a {
                    case(id).a.locate(&point)
                } else {
                    Location::Outside
                };
                query(&mut m, answer.body, &point, expected, &mut dec);
            }
        }
    }
    let mut m = Model::new();
    let a = empty(&mut m);
    let b = empty(&mut m);
    for op in [BoolOp::Unite, BoolOp::Intersect, BoolOp::Subtract] {
        let answer = boolean(&mut m, a, b, op, &mut Decider::default()).unwrap();
        check_counts(&m, answer.body, &Counts::empty());
    }
});

blind_test!(case_47_connected_curve_contact_refuses_at_sew, {
    let mut m = Model::new();
    let a = block_at([10; 3], [0; 3]).build(&mut m);
    let b = cyl_at(5, 6, [5, 5, 2], false).build(&mut m);
    let before = topology_snapshot(&m, &[a, b]);
    let refusal = boolean(&mut m, a, b, BoolOp::Subtract, &mut Decider::default()).unwrap_err();
    assert_eq!(refusal.stage, Stage::Sew);
    assert_eq!(refusal.refusal.code, "boolean.self-touching-shell");
    assert_eq!(before, topology_snapshot(&m, &[a, b]));
});

blind_test!(
    marks_follow_operation_then_first_then_second_even_with_empty,
    {
        for op in [BoolOp::Unite, BoolOp::Intersect, BoolOp::Subtract] {
            for empty_first in [false, true] {
                let mut m = Model::new();
                let make = |m: &mut Model, body: BodyId, n: i64| {
                    let mut copy = m.body(body).clone();
                    for bound in [n, n + 1] {
                        copy.marks.push(Mark::new(
                            "numbers.budget-exhausted",
                            BoundKind::Distance,
                            Rational::from_ratio(bound, 100).unwrap(),
                            format!("input {n}"),
                        ));
                    }
                    m.add_body(copy)
                };
                let raw_a = if empty_first {
                    empty(&mut m)
                } else {
                    case("01").a.build(&mut m)
                };
                let a = make(&mut m, raw_a, 1);
                let raw_b = case("01").b.build(&mut m);
                let b = make(&mut m, raw_b, 3);
                let mut expected = MarkLog::new();
                expected.extend(&m.body(a).marks);
                expected.extend(&m.body(b).marks);
                let answer = boolean(&mut m, a, b, op, &mut Decider::default()).unwrap();
                assert_eq!(answer.marks, expected);
                assert_eq!(m.body(answer.body).marks, expected);
                assert_eq!(
                    answer.report.stages.iter().map(|s| s.marks).sum::<usize>(),
                    0
                );
            }
        }
    }
);

blind_test!(face_point_refuses_a_carrier_with_no_rational_points, {
    // On the cylinder about (1,1,1), project a rational point onto x+y+z=0.
    // Write that projection (a,b,-a-b). Radius 1 requires 2(a²+ab+b²)=1.
    // Clearing reduced denominators gives 2(A²+AB+B²)=D². The quadratic is
    // odd unless A and B are both even, so its 2-adic valuation is even.
    // The left side has odd valuation, the square even: no rational solution.
    let mut m = Model::new();
    let surface = m.geom_mut().add_surface(Surface::Cylinder(
        Cylinder::new(
            RPoint3::origin().into(),
            v(1, 1, 1).into(),
            r(1).into(),
            v(1, -1, 0).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let f = m.add_face(Face {
        surface,
        sense: Sense::Same,
        loops: vec![],
    });
    let refusal = face_point(&mut m, f, &mut Decider::default()).unwrap_err();
    assert_eq!(refusal.code, "boolean.no-rational-point");
});
