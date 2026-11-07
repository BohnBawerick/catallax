//! Independent section integrals, endpoint orientation and exact apex consumers.
use std::collections::BTreeMap;

use catallax_boolean::{Containment, classify_point};
use catallax_features::{Built, evaluate};
use catallax_measure::{Accuracy, Enclosure};
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::{Ball, BudgetPolicy, Decider, Outcome, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::{RPoint3, SplineSupport, Surface};
use catallax_topology::{LinearFilletBand, Sense};
use serde_json::json;

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn b(n: i64) -> Ball {
    Ball::from_rational(&Rational::from_i64(n), 192)
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 100_000_000,
        ..Default::default()
    })
}
fn recipe(start: i64, end: i64, backwards: bool) -> Recipe {
    let side = |name: &str| json!({"step":"b","role":{"side":{"segment":name}}});
    let top = json!({"step":"b","role":"end"});
    Recipe::parse(&json!({"format":"catallax-recipe/5","name":"linear control","result":"round","steps":[
        {"id":"b","op":"block","size":["20","20","10"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
        {"id":"round","op":"fillet_variable","body":"b","edges":[{
            "edge":{"between":[top,side("x_max")]},
            "from":{"faces":[top,side("x_max"),side(if backwards {"y_max"} else {"y_min"})]},
            "radius":{"kind":"linear","start":start.to_string(),"end":end.to_string()}}]}]
    }).to_string()).unwrap()
}
fn build(recipe: &Recipe) -> Built {
    evaluate(recipe, &mut decider()).unwrap()
}
fn encloses(actual: &Enclosure, expected: &Ball) {
    assert!(
        actual.lo <= expected.lo() && actual.hi >= expected.hi(),
        "{actual:?} does not enclose {expected:?}"
    );
}

fn sum<const N: usize>(terms: [&Ball; N]) -> Ball {
    terms.into_iter().fold(b(0), |total, value| &total + value)
}
fn product<const N: usize>(factors: [&Ball; N]) -> Ball {
    factors
        .into_iter()
        .fold(b(1), |total, value| &total * value)
}
fn quotient(numerator: &Ball, denominator: &Ball) -> Ball {
    numerator * &denominator.recip().unwrap()
}

// Integrate the unit square minus its tangent quarter disk. Its area, inward
// first moment, inward square moment and product moment are k, m, n and p.
// This section calculation is independent of the rational carrier and flux code.
fn oracle(a: i64, z: i64) -> (Ball, [Ball; 3], [[Ball; 3]; 3]) {
    let pi = Ball::pi(192);
    let k = &b(1) - &(quotient(&pi, &b(4)));
    let m = &(quotient(&b(5), &b(6))) - &(quotient(&pi, &b(4)));
    let n = &b(1) - &(quotient(&(&b(5) * &pi), &b(16)));
    let p = &(quotient(&b(19), &b(24))) - &(quotient(&pi, &b(4)));
    let (a, z, d) = (b(a), b(z), b(z - a));
    let (a2, z2, d2) = (&a * &a, &z * &z, &d * &d);
    let (a3, z3, d3) = (&a2 * &a, &z2 * &z, &d2 * &d);
    let i2 = quotient(&product([&b(20), &sum([&a2, &(&a * &z), &z2])]), &b(3));
    let i3 = quotient(
        &product([&b(20), &sum([&a3, &(&a2 * &z), &(&a * &z2), &z3])]),
        &b(4),
    );
    let i4 = &b(4)
        * &sum([
            &(&a3 * &a),
            &(&a3 * &z),
            &(&a2 * &z2),
            &(&a * &z3),
            &(&z3 * &z),
        ]);
    let iy2 = &b(400)
        * &sum([
            &(quotient(&a2, &b(2))),
            &(quotient(&product([&b(2), &a, &d]), &b(3))),
            &(quotient(&d2, &b(4))),
        ]);
    let iy3 = &b(400)
        * &sum([
            &(quotient(&a3, &b(2))),
            &(&a2 * &d),
            &(quotient(&product([&b(3), &a, &d2]), &b(4))),
            &(quotient(&d3, &b(5))),
        ]);
    let iyy2 = &b(8000)
        * &sum([
            &(quotient(&a2, &b(3))),
            &(quotient(&(&a * &d), &b(2))),
            &(quotient(&d2, &b(5))),
        ]);
    let volume = &b(4000) - &(&k * &i2);
    let centroid = [
        quotient(
            &sum([&b(40000), &(-product([&b(20), &k, &i2])), &(&m * &i3)]),
            &volume,
        ),
        quotient(&(&b(40000) - &(&k * &iy2)), &volume),
        quotient(
            &sum([&b(20000), &(-product([&b(10), &k, &i2])), &(&m * &i3)]),
            &volume,
        ),
    ];
    let xx = sum([
        &(quotient(&b(1600000), &b(3))),
        &(-product([&b(400), &k, &i2])),
        &product([&b(40), &m, &i3]),
        &(-(&n * &i4)),
    ]);
    let yy = &(quotient(&b(1600000), &b(3))) - &(&k * &iyy2);
    let zz = sum([
        &(quotient(&b(400000), &b(3))),
        &(-product([&b(100), &k, &i2])),
        &product([&b(20), &m, &i3]),
        &(-(&n * &i4)),
    ]);
    let xy = sum([&b(400000), &(-product([&b(20), &k, &iy2])), &(&m * &iy3)]);
    let xz = sum([
        &b(200000),
        &(-product([&b(200), &k, &i2])),
        &product([&b(30), &m, &i3]),
        &(-(&p * &i4)),
    ]);
    let yz = sum([&b(200000), &(-product([&b(10), &k, &iy2])), &(&m * &iy3)]);
    let raw = [
        [xx, xy.clone(), xz.clone()],
        [xy, yy, yz.clone()],
        [xz, yz, zz],
    ];
    let centered: [[Ball; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| &raw[i][j] - &product([&volume, &centroid[i], &centroid[j]]))
    });
    let inertia = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            if i == j {
                (0..3)
                    .filter(|k| *k != i)
                    .fold(b(0), |sum, k| &sum + &centered[k][k])
            } else {
                -&centered[i][j]
            }
        })
    });
    (volume, centroid, inertia)
}

#[test]
fn linear_sections_match_independent_volume_centroid_and_inertia() {
    for (start, end) in [(1, 2), (0, 2)] {
        let built = build(&recipe(start, end, false));
        assert!(built.model.check(built.result).is_empty());
        let answer = catallax_measure::mass_properties(
            &built.model,
            built.result,
            &Accuracy {
                relative: Rational::zero(),
                absolute: q("1e-7"),
            },
            &mut decider(),
        )
        .unwrap();
        assert!(answer.marks.is_empty(), "{:?}", answer.marks);
        let (volume, centroid, inertia) = oracle(start, end);
        encloses(&answer.value.volume, &volume);
        for (actual, expected) in answer.value.centroid.as_ref().unwrap().iter().zip(centroid) {
            encloses(actual, &expected);
        }
        for (actual, expected) in answer.value.inertia.as_ref().unwrap().iter().zip(inertia) {
            for (actual, expected) in actual.iter().zip(expected) {
                encloses(actual, &expected);
            }
        }
    }
}

#[test]
fn named_endpoint_reversal_changes_material_at_the_authored_station() {
    for (backwards, want) in [(false, Containment::Inside), (true, Containment::Outside)] {
        let mut built = build(&recipe(1, 2, backwards));
        let point = RPoint3::new(q("19.55"), q("5"), q("9.55"));
        let got = classify_point(&mut built.model, built.result, &point, &mut decider());
        assert!(
            matches!(got,Outcome::Proven(found) if found == want),
            "{got:?}, wanted {want:?}"
        );
    }
}

#[test]
fn either_zero_end_has_one_shared_apex_and_a_closed_bounded_display() {
    for (start, end, backwards) in [(0, 2, false), (2, 0, false), (0, 2, true), (2, 0, true)] {
        let mut built = build(&recipe(start, end, backwards));
        let face = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .find(|f| {
                matches!(
                    built.model.geom().surface(built.model.face(*f).surface),
                    Surface::BSpline(SplineSupport::LinearFillet(_))
                )
            })
            .unwrap();
        assert!(
            LinearFilletBand::of_face(&built.model, face, &mut decider())
                .unwrap()
                .is_some()
        );
        let mut wrong = built.model.face(face).clone();
        wrong.sense = if wrong.sense == Sense::Same {
            Sense::Reversed
        } else {
            Sense::Same
        };
        let wrong = built.model.add_face(wrong);
        assert!(LinearFilletBand::of_face(&built.model, wrong, &mut decider()).is_err());
        let options = MeshOptions {
            chord_tolerance: q("0.01"),
            ..Default::default()
        };
        let mesh = mesh_body(&built.model, built.result, &options, &mut decider()).unwrap();
        assert!(mesh.marks.is_empty(), "{:?}", mesh.marks);
        let mut incidence = BTreeMap::<_, (usize, i32)>::new();
        for triangle in &mesh.value.triangles {
            let [a, b, c] = triangle.vertices;
            assert!(a != b && a != c && b != c, "collapsed apex triangle");
            for (a, b) in [(a, b), (b, c), (c, a)] {
                let uses = incidence.entry((a.min(b), a.max(b))).or_default();
                uses.0 += 1;
                uses.1 += if a < b { 1 } else { -1 };
            }
        }
        assert!(
            incidence.values().all(|v| *v == (2, 0)),
            "unshared or misoriented boundary"
        );
        assert!(
            mesh.value
                .faces
                .iter()
                .all(|f| f.max_chord_error <= options.chord_tolerance)
        );
        assert_eq!(
            mesh.value
                .faces
                .iter()
                .map(|f| f.face)
                .collect::<std::collections::BTreeSet<_>>(),
            built
                .model
                .faces_of_body(built.result)
                .into_iter()
                .collect::<std::collections::BTreeSet<_>>(),
            "display omitted or mislabelled a physical face"
        );
        let apex = [
            20.0,
            if (start == 0) != backwards { 0.0 } else { 20.0 },
            10.0,
        ];
        assert_eq!(
            mesh.value.positions.iter().filter(|p| **p == apex).count(),
            1,
            "split zero-radius apex"
        );
        let mut bytes = Vec::new();
        let report = catallax_step::write_step(
            &built.model,
            &[built.result],
            &Default::default(),
            &mut decider(),
            &mut bytes,
        )
        .unwrap();
        assert!(report.marks.is_empty(), "{:?}", report.marks);
    }
}

#[test]
fn acute_planar_contact_has_the_exact_nonquadrant_section_volume() {
    let side = |segment: &str| json!({"step":"b","role":{"side":{"segment":segment}}});
    let doc = json!({"format":"catallax-recipe/5","name":"acute linear law","result":"round","steps":[
        {"id":"s","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[
            {"kind":"line","id":"base","to":["20","0"]},
            {"kind":"line","id":"diagonal","to":["0","20"]},
            {"kind":"line","id":"west","to":["0","0"]}]]},
        {"id":"b","op":"extrude","profile":"s","distance":"20"},
        {"id":"round","op":"fillet_variable","body":"b","edges":[{
            "edge":{"between":[side("base"),side("diagonal")]},
            "from":{"faces":[{"step":"b","role":"start"},side("base"),side("diagonal")]},
            "radius":{"kind":"linear","start":"1","end":"2"}}]}]});
    let built = build(&Recipe::parse(&doc.to_string()).unwrap());
    let expected = b(4000)
        - quotient(&b(140), &b(3))
            * (b(1) + b(2).sqrt().unwrap() - quotient(&(b(3) * Ball::pi(192)), &b(8)));
    let actual = catallax_measure::volume(
        &built.model,
        built.result,
        &Accuracy {
            relative: Rational::zero(),
            absolute: q("1e-8"),
        },
        &mut decider(),
    )
    .unwrap();
    assert!(actual.marks.is_empty());
    encloses(&actual.value, &expected);
}

#[test]
fn disjoint_linear_and_constant_laws_share_one_checked_construction() {
    for (start, end, other) in [(1, 2, 1), (1, 1, 2)] {
        let mut doc = serde_json::to_value(recipe(start, end, false)).unwrap();
        let side = |name: &str| json!({"step":"b","role":{"side":{"segment":name}}});
        let bottom = json!({"step":"b","role":"start"});
        doc["steps"][1]["edges"]
            .as_array_mut()
            .unwrap()
            .push(json!({
                "edge":{"between":[bottom,side("x_min")]},
                "from":{"faces":[bottom,side("x_min"),side("y_min")]},
                "radius":{"kind":"linear","start":other.to_string(),"end":other.to_string()}
            }));
        let mut built = build(&Recipe::parse(&doc.to_string()).unwrap());
        let coefficient = &b(1) - &(quotient(&Ball::pi(192), &b(4)));
        let average_square =
            &(quotient(&b(start * start + start * end + end * end), &b(3))) + &b(other * other);
        let expected = &b(4000) - &product([&b(20), &average_square, &coefficient]);
        let volume = catallax_measure::volume(
            &built.model,
            built.result,
            &Accuracy {
                relative: Rational::zero(),
                absolute: q("1e-8"),
            },
            &mut decider(),
        )
        .unwrap();
        assert!(volume.marks.is_empty());
        encloses(&volume.value, &expected);
        for (point, want) in [
            (["0.1", "5", "0.1"], Containment::Outside),
            (["1", "5", "1"], Containment::Inside),
        ] {
            let [x, y, z] = point.map(q);
            let got = classify_point(
                &mut built.model,
                built.result,
                &RPoint3::new(x, y, z),
                &mut decider(),
            );
            assert!(
                matches!(got,Outcome::Proven(found) if found == want),
                "{got:?}"
            );
        }
    }
}

#[test]
fn shared_endpoint_laws_reject_conflicting_values_before_corner_construction() {
    let mut doc = serde_json::to_value(recipe(1, 2, false)).unwrap();
    let side = |name: &str| json!({"step":"b","role":{"side":{"segment":name}}});
    let top = json!({"step":"b","role":"end"});
    doc["steps"][1]["edges"]
        .as_array_mut()
        .unwrap()
        .push(json!({
            "edge":{"between":[top,side("y_min")]},
            "from":{"faces":[top,side("x_max"),side("y_min")]},
            "radius":{"kind":"linear","start":"3","end":"4"}
        }));
    let recipe = Recipe::parse(&doc.to_string()).unwrap();
    let error = evaluate(&recipe, &mut decider()).unwrap_err();
    let catallax_features::BuildError::Geometry { refusal, .. } = error else {
        panic!("{error:?}");
    };
    assert_eq!(refusal.code, "fillet.incompatible-radius-laws");
}

#[test]
fn rotated_and_reflected_zero_end_copies_keep_volume_membership_and_band_orientation() {
    use catallax_shapes::{Isometry, Pose, RVector3, Rotation};
    let pose = Pose::from_rational(
        Rotation::from_quaternion(q("1"), q("2"), q("3"), q("4")).unwrap(),
        RVector3::new(q("1000"), q("-500"), q("10")),
    );
    let maps = [
        Isometry::from(&pose),
        Isometry::reflection(&RPoint3::origin(), &RVector3::new(q("1"), q("2"), q("3"))).unwrap(),
    ];
    let mut built = build(&recipe(0, 2, false));
    let (expected, _, _) = oracle(0, 2);
    for map in maps {
        let copy =
            catallax_topology::transformed_body(&mut built.model, built.result, &map, Clone::clone);
        assert!(built.model.check(copy).is_empty());
        assert!(LinearFilletBand::in_body(&built.model, copy, &mut decider()).unwrap());
        let measured = catallax_measure::volume(
            &built.model,
            copy,
            &Accuracy {
                relative: Rational::zero(),
                absolute: q("1e-8"),
            },
            &mut decider(),
        )
        .unwrap();
        assert!(measured.marks.is_empty());
        encloses(&measured.value, &expected);
        for (p, want) in [
            (["19.87", "5", "9.87"], Containment::Outside),
            (["19.8", "5", "9.8"], Containment::Inside),
        ] {
            let [x, y, z] = p.map(q);
            let point = map.apply_rational_point(&RPoint3::new(x, y, z)).unwrap();
            let found = classify_point(&mut built.model, copy, &point, &mut decider());
            assert!(
                matches!(found,Outcome::Proven(value) if value == want),
                "{found:?}"
            );
        }
    }
}

#[test]
fn constant_linear_law_keeps_consumed_face_history_at_the_face_limit() {
    let recipe = recipe(10, 10, false);
    let built = build(&recipe);
    let source: catallax_recipe::Origin = serde_json::from_value(json!({
        "step":"b","role":{"side":{"segment":"x_max"}}
    }))
    .unwrap();
    assert!(
        built.consumed.iter().any(|record| {
            record.step == recipe.result
                && record.source.contains(&source)
                && record.by.iter().any(|origin| {
                    origin.step == recipe.result
                        && matches!(&origin.role, catallax_recipe::Role::Blend { .. })
                })
        }),
        "consumed source face lost its replacement history"
    );
    assert!(
        built
            .model
            .faces_of_body(built.result)
            .iter()
            .all(|face| { !built.model.face_origins(*face).contains(&source) }),
        "removed source origin became selectable on a replacement"
    );
    let expected = &b(2000) + &(&b(500) * &Ball::pi(192));
    let actual = catallax_measure::volume(
        &built.model,
        built.result,
        &Accuracy {
            relative: Rational::zero(),
            absolute: q("1e-7"),
        },
        &mut decider(),
    )
    .unwrap();
    assert!(actual.marks.is_empty());
    encloses(&actual.value, &expected);
}
