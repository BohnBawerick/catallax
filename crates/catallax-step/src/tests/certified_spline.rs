use std::collections::{BTreeMap, BTreeSet};

use catallax_features::{CertifiedBuilt, evaluate_certified};
use catallax_numbers::{Ball, Decider, Outcome, Rational, Real, StrictSign, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{BSplineCurve, BSplineSurface, Point3, Vector3};

use super::{cartesian, entities, fields, q, references, step_real};
use crate::{CertifiedSplineStepOptions, write_spline_prism_step};

fn build(dec: &mut Decider) -> Result<CertifiedBuilt, Box<dyn std::error::Error>> {
    let recipe = Recipe::parse(
        r#"{
        "format":"catallax-recipe/3", "name":"Certified spline STEP prism",
        "steps":[
            {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[{"kind":"circle","id":"wall","centre":["0","0"],"radius":"2"}]]},
            {"id":"sweep","op":"extrude","profile":"profile","distance":"10","direction":["3","0","4"]},
            {"id":"window","op":"block","size":["8","8","4"],"pose":{"rotation":["3","0","1","0"],"translation":["-2","-4","4"]}},
            {"id":"source","op":"boolean","kind":"intersect","target":"sweep","tool":"window"},
            {"id":"offset","op":"offset","body":"source","distance":"-1/4"}
        ], "result":"offset"
    }"#,
    )?;
    Ok(evaluate_certified(&recipe, dec)?)
}

fn expanded(multiplicities: &str, knots: &str) -> Vec<Rational> {
    fields(multiplicities.trim_matches(['(', ')']))
        .into_iter()
        .zip(fields(knots.trim_matches(['(', ')'])))
        .flat_map(|(multiplicity, knot)| {
            std::iter::repeat_n(step_real(knot), multiplicity.parse().unwrap())
        })
        .collect()
}
fn real_point(values: [Rational; 3]) -> Point3 {
    let [x, y, z] = values.map(Real::from_rational);
    Point3::new(x, y, z)
}
fn error(a: &[Ball; 3], b: &[Ball; 3]) -> Rational {
    a.iter()
        .zip(b)
        .map(|(a, b)| {
            let d = a - b;
            d.lo().abs().max(d.hi().abs())
        })
        .fold(Rational::zero(), |a, b| a + b)
}
fn dot(a: &[Rational; 3], b: &[Rational; 3]) -> Rational {
    (0..3)
        .map(|i| &a[i] * &b[i])
        .fold(Rational::zero(), |a, b| a + b)
}

type Entities<'a> = [(&'a str, Vec<&'a str>)];
type Chart = (Vec<Vec<usize>>, BSplineSurface);

fn written_chart(
    es: &Entities<'_>,
    surface: &[&str],
    dec: &mut Decider,
) -> Result<Chart, Box<dyn std::error::Error>> {
    let rows: Vec<_> = fields(
        surface[3]
            .strip_prefix('(')
            .unwrap()
            .strip_suffix(')')
            .unwrap(),
    )
    .into_iter()
    .map(references)
    .collect();
    let points: Vec<_> = rows
        .iter()
        .flat_map(|row| row.iter().map(|id| real_point(cartesian(&es[*id - 1].1))))
        .collect();
    let chart = BSplineSurface::new(
        [surface[1].parse()?, surface[2].parse()?],
        [
            expanded(surface[8], surface[10]),
            expanded(surface[9], surface[11]),
        ],
        [rows.len(), rows[0].len()],
        points,
        vec![Real::one(); rows.len() * rows[0].len()],
        dec,
    )?;
    Ok((rows, chart))
}

#[test]
fn actual_step_geometry_retains_shared_caps_periodic_c1_and_oriented_solid_incidence() {
    let (result, work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let built = build(&mut dec)?;
        let mut bytes = Vec::new();
        let options = CertifiedSplineStepOptions::default();
        let report = write_spline_prism_step(
            &built.built().model,
            built.built().result,
            &options,
            &mut dec,
            &mut bytes,
        )?
        .value;
        assert!(report.max_deviation <= options.file_tolerance);
        let file = String::from_utf8(bytes)?;
        let es = entities(&file);
        let mut charts = es
            .iter()
            .filter(|(name, _)| *name == "B_SPLINE_SURFACE_WITH_KNOTS")
            .map(|(_, surface)| written_chart(&es, surface, &mut dec))
            .collect::<Result<Vec<_>, _>>()?;
        charts.sort_by(|a, b| a.1.knots()[0][0].cmp(&b.1.knots()[0][0]));
        let [(lower, first), (upper, second)] = charts.as_slice() else {
            panic!("expected two open side charts");
        };
        for (chart, u) in [(first, ["0", "2"]), (second, ["2", "4"])] {
            assert_eq!(chart.domain()[0], [&q(u[0]), &q(u[1])]);
            assert_eq!(chart.domain()[1], [&q("0"), &q("1")]);
        }
        let source = built.result_certificate().spline().unwrap().carrier();
        for (chart, parameters) in [(first, ["1/3", "7/6"]), (second, ["19/7", "11/3"])] {
            for parameter in parameters {
                let u = Ball::from_rational(&q(parameter), 128);
                let v = Ball::from_rational(&q("2/3"), 128);
                assert!(
                    error(
                        &chart.evaluate(&u, &v, &dec)?,
                        &source.evaluate(&u, &v, &dec)?
                    ) <= report.max_deviation
                );
            }
        }
        // Read the actual serialized knot cover, not the producer's table.
        let v = Ball::from_rational(&q("1/3"), 128);
        for chart in [first, second] {
            for knot in chart.knots()[0].iter().collect::<BTreeSet<_>>() {
                let derivatives = chart.derivatives(&Ball::from_rational(knot, 128), &v, &dec)?;
                assert!(
                    derivatives[0]
                        .iter()
                        .all(|value| value.hi() - value.lo() < q("1e-30"))
                );
            }
        }
        // Both joins keep the position and C1 tangent of the closed profile.
        for (end, start, a, b) in [(first, second, "2", "2"), (second, first, "4", "0")] {
            let a = Ball::from_rational(&q(a), 128);
            let b = Ball::from_rational(&q(b), 128);
            assert!(
                error(&end.evaluate(&a, &v, &dec)?, &start.evaluate(&b, &v, &dec)?) < q("1e-30")
            );
            assert!(
                error(
                    &end.derivatives(&a, &v, &dec)?[0],
                    &start.derivatives(&b, &v, &dec)?[0]
                ) < q("1e-30")
            );
        }
        let curves: Vec<_> = es
            .iter()
            .filter(|(name, _)| *name == "B_SPLINE_CURVE_WITH_KNOTS")
            .map(|(_, curve)| curve)
            .collect();
        assert_eq!(curves.len(), 6);
        assert_eq!(curves.iter().filter(|curve| curve[1] == "3").count(), 4);
        for curve in curves.iter().filter(|curve| curve[1] == "3") {
            let ids = references(curve[2]);
            let (rows, chart, parameters) = if lower.iter().any(|row| row.contains(&ids[1])) {
                (lower, first, ["0", "5/7", "13/9", "2"])
            } else {
                (upper, second, ["2", "13/6", "23/7", "4"])
            };
            let row = rows[0].iter().position(|id| *id == ids[0]).unwrap();
            assert!(
                ids.iter()
                    .zip(rows.iter())
                    .all(|(id, pair)| *id == pair[row])
            );
            let controls: Vec<_> = ids
                .iter()
                .map(|id| real_point(cartesian(&es[*id - 1].1)))
                .collect();
            let count = controls.len();
            let curve = BSplineCurve::new(
                curve[1].parse()?,
                expanded(curve[6], curve[7]),
                controls,
                vec![Real::one(); count],
                &mut dec,
            )?;
            for parameter in parameters {
                let u = Ball::from_rational(&q(parameter), 128);
                let v = Ball::from_rational(&Rational::from_i64(row as i64), 128);
                assert!(
                    error(&curve.evaluate(&u, &dec)?, &chart.evaluate(&u, &v, &dec)?) < q("1e-30")
                );
            }
        }
        // The two axial edges are both charts' shared end rows.
        let joins: BTreeSet<_> = curves
            .iter()
            .filter(|curve| curve[1] == "1")
            .map(|curve| references(curve[2]))
            .collect();
        assert_eq!(
            joins,
            BTreeSet::from([lower[0].to_vec(), upper[0].to_vec()])
        );
        assert_eq!(lower.last(), upper.first());
        assert_eq!(upper.last(), lower.first());
        let rows: Vec<_> = lower.iter().chain(&upper[1..]).collect();
        let faces: Vec<_> = es
            .iter()
            .filter(|(name, _)| *name == "ADVANCED_FACE")
            .map(|(_, face)| face)
            .collect();
        assert_eq!(faces.len(), 4);
        assert_eq!(
            es.iter()
                .filter(|(name, _)| *name == "VERTEX_POINT")
                .count(),
            4
        );
        assert_eq!(
            (report.faces, report.edges, report.vertices),
            (4, 6, 4),
            "exported topology, not the native model's"
        );
        assert_eq!((report.added_vertices, report.added_edges), (4, 4));
        let side = built.built().model.faces_of_body(built.built().result);
        assert_eq!(side.len(), 3);
        assert_eq!(report.face_map.len(), 4);
        for (written, face) in report.face_map.iter().zip(&faces) {
            assert_eq!(face[0], format!("'{}'", written.file_region));
            assert!(side.contains(&written.face));
        }
        let [_, _, a, b] = report.face_map.as_slice() else {
            panic!("expected two caps and two side patches");
        };
        assert_eq!(a.face, b.face);
        assert_eq!(a.file_domain, Some([[q("0"), q("2")], [q("0"), q("1")]]));
        assert_eq!(b.file_domain, Some([[q("2"), q("4")], [q("0"), q("1")]]));
        let mut cap_rows = BTreeSet::new();
        let mut incidences = BTreeMap::<usize, (usize, i32)>::new();
        for face in faces {
            let surface = references(face[2])[0];
            if es[surface - 1].0 == "PLANE" {
                let placement = references(es[surface - 1].1[1])[0];
                let placement = &es[placement - 1].1;
                let origin = cartesian(&es[references(placement[1])[0] - 1].1);
                let normal = cartesian(&es[references(placement[2])[0] - 1].1);
                let matching: Vec<_> = (0..2)
                    .filter(|row| {
                        rows.iter().all(|pair| {
                            let p = cartesian(&es[pair[*row] - 1].1);
                            dot(&std::array::from_fn(|i| &p[i] - &origin[i]), &normal).is_zero()
                        })
                    })
                    .collect();
                assert_eq!(matching.len(), 1);
                let row = matching[0];
                cap_rows.insert(row);
                let [x, y, z] = normal.map(Real::from_rational);
                let normal = Vector3::new(x, y, z);
                let material = normal.dot(built.result_certificate().target().axis())
                    * Real::from_i64(if face[3] == ".T." { 1 } else { -1 });
                let expected = if row == 0 {
                    StrictSign::Negative
                } else {
                    StrictSign::Positive
                };
                assert!(
                    matches!(dec.sign(&material), Outcome::Proven(actual) if actual == expected)
                );
            }
            for bound in references(face[1]) {
                let bound = &es[bound - 1].1;
                let loop_ = references(bound[1])[0];
                for oriented in references(es[loop_ - 1].1[1]) {
                    let oriented = &es[oriented - 1].1;
                    let edge = references(oriented[3])[0];
                    let incidence = incidences.entry(edge).or_default();
                    incidence.0 += 1;
                    incidence.1 += if (oriented[4] == ".T.") == (bound[2] == ".T.") {
                        1
                    } else {
                        -1
                    };
                    let edge = &es[edge - 1].1;
                    assert_ne!(
                        references(edge[1]),
                        references(edge[2]),
                        "every chart edge is open between distinct vertices"
                    );
                }
            }
        }
        assert_eq!(cap_rows, BTreeSet::from([0, 1]));
        assert_eq!(incidences.len(), 6);
        assert!(incidences.values().all(|uses| *uses == (2, 0)));
        Ok::<_, Box<dyn std::error::Error>>(())
    });
    eprintln!("certified actual STEP consumer request: {work:?}");
    result.unwrap();
    assert!(!work.exhausted);
}

#[test]
fn coarse_rounding_refuses_instead_of_widening_accuracy_or_publishing_partial_geometry() {
    let (result, work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let built = build(&mut dec)?;
        let options = CertifiedSplineStepOptions {
            decimal_places: 0,
            file_tolerance: q("1e-12"),
            ..CertifiedSplineStepOptions::default()
        };
        let mut bytes = b"existing sink contents".to_vec();
        let refusal = write_spline_prism_step(
            &built.built().model,
            built.built().result,
            &options,
            &mut dec,
            &mut bytes,
        )
        .unwrap_err();
        assert_eq!(refusal.code, "step.tolerance-not-reached");
        assert_eq!(bytes, b"existing sink contents");
        Ok::<_, Box<dyn std::error::Error>>(())
    });
    eprintln!("certified coarse STEP request: {work:?}");
    result.unwrap();
    assert!(!work.exhausted);
}
