//! Non-gating measurements of the contract's analytic solids.

#[path = "../src/fixtures.rs"]
mod fixtures;

use catallax_measure::{Accuracy, Enclosure, area, face_area, mass_properties, volume};
use catallax_numbers::{Answered, Decider, Rational};
use catallax_shapes::{Pose, RVector3, Rotation};
use catallax_topology::{BodyId, Model};
use std::time::Instant;

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}

fn primitive(kind: usize) -> (Model, BodyId) {
    let mut model = Model::new();
    let pose = Pose::identity();
    let body = match kind {
        0 => catallax_topology::block(&mut model, [r(2), r(3), r(5)], &pose).unwrap(),
        1 => catallax_topology::block(
            &mut model,
            [
                Rational::parse("1/3").unwrap(),
                r(7),
                Rational::parse("1e-3").unwrap(),
            ],
            &pose,
        )
        .unwrap(),
        2 => catallax_topology::cylinder(&mut model, r(2), r(5), &pose).unwrap(),
        3 => catallax_topology::sphere(&mut model, r(5), &pose).unwrap(),
        4 => catallax_topology::cone(&mut model, r(2), r(1), r(2), &pose).unwrap(),
        5 => catallax_topology::cone(&mut model, r(3), r(0), r(4), &pose).unwrap(),
        6 => catallax_topology::torus(&mut model, r(3), r(1), &pose).unwrap(),
        7 => catallax_topology::empty(&mut model),
        _ => unreachable!(),
    };
    (model, body)
}

fn timed<T>(
    case: &str,
    property: &str,
    accuracy: &Accuracy,
    f: impl FnOnce(&mut Decider) -> Answered<T>,
    widths: impl FnOnce(&T) -> Vec<Rational>,
) {
    let mut dec = Decider::default();
    let start = Instant::now();
    let answer = f(&mut dec);
    let elapsed = start.elapsed().as_nanos();
    match answer {
        Ok(answer) => println!(
            "case={case} property={property} relative={} absolute={} ns={elapsed} marks={} widths={:?} decider={:?}",
            accuracy.relative,
            accuracy.absolute,
            answer.marks.len(),
            widths(&answer.value),
            dec.stats()
        ),
        Err(refusal) => println!(
            "case={case} property={property} relative={} absolute={} ns={elapsed} refusal={} decider={:?}",
            accuracy.relative,
            accuracy.absolute,
            refusal.code,
            dec.stats()
        ),
    }
}

fn scalar_width(enclosure: &Enclosure) -> Vec<Rational> {
    vec![enclosure.width()]
}

fn main() {
    let names = [
        "block",
        "fractional-block",
        "cylinder",
        "sphere",
        "frustum",
        "cone-apex",
        "torus",
        "empty",
    ];
    let mut cases = names
        .into_iter()
        .enumerate()
        .map(|(kind, name)| {
            let (model, body) = primitive(kind);
            (name, model, body)
        })
        .collect::<Vec<_>>();
    let half_negative = Pose::new(
        Rotation::from_quaternion(r(0), r(1), r(0), r(0)).unwrap(),
        (RVector3::new(r(0), r(0), r(3))).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap();
    for (name, (model, body)) in [
        (
            "through-hole",
            fixtures::through_hole(&r(1), &Pose::identity()),
        ),
        ("void", fixtures::void_block()),
        ("hemisphere", fixtures::cap(2, 0)),
        ("sphere-cap", fixtures::cap(5, 3)),
        ("octant", fixtures::octant(&Pose::identity())),
        (
            "half-cylinder-positive",
            fixtures::half_cylinder(&Pose::identity()),
        ),
        (
            "half-cylinder-negative",
            fixtures::half_cylinder(&half_negative),
        ),
        ("steinmetz", fixtures::steinmetz(&Pose::identity())),
        ("half-torus", fixtures::half_torus()),
        ("coplanar-split", fixtures::coplanar_split_block()),
    ] {
        cases.push((name, model, body));
    }
    let filter = std::env::var("CATALLAX_MEASURE_BENCH_CASE").ok();
    let quick = std::env::var_os("CATALLAX_MEASURE_BENCH_QUICK").as_deref()
        == Some(std::ffi::OsStr::new("1"));
    let accuracies = if quick {
        vec![Accuracy::default()]
    } else {
        [6, 12, 18, 24, 30]
            .into_iter()
            .map(|digits| Accuracy {
                relative: Rational::parse(&format!("1e-{digits}")).unwrap(),
                absolute: Rational::parse("1e-40").unwrap(),
            })
            .collect()
    };
    for (name, model, body) in cases {
        if filter.as_ref().is_some_and(|filter| !name.contains(filter)) {
            continue;
        }
        for accuracy in &accuracies {
            timed(
                name,
                "volume",
                accuracy,
                |dec| volume(&model, body, accuracy, dec),
                scalar_width,
            );
            timed(
                name,
                "area",
                accuracy,
                |dec| area(&model, body, accuracy, dec),
                scalar_width,
            );
            timed(
                name,
                "mass_properties",
                accuracy,
                |dec| mass_properties(&model, body, accuracy, dec),
                |value| {
                    let mut widths = vec![value.volume.width(), value.area.width()];
                    if let Some(center) = &value.centroid {
                        widths.extend(center.iter().map(Enclosure::width));
                    }
                    widths
                },
            );
            for (index, face) in model.faces_of_body(body).into_iter().enumerate() {
                timed(
                    name,
                    &format!("face_area.{index}"),
                    accuracy,
                    |dec| face_area(&model, face, accuracy, dec),
                    scalar_width,
                );
            }
        }
    }
}
