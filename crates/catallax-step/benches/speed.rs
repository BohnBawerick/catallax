//! Timings only. No speed threshold determines success.

use std::hint::black_box;
use std::time::Instant;

use catallax_intersect::{
    Component, CurveForm, SurfacePair, curve_form, intersect_surfaces, trace,
};
use catallax_numbers::{Decider, Rational};
use catallax_shapes::{Cylinder, GeomStore, Pose, RPoint3, RVector3, Sphere, Surface};
use catallax_step::{StepOptions, write_step};
use catallax_topology::{Model, block, cone, cylinder, sphere, torus};

// Compile the actual private fitter here instead of maintaining a second
// implementation or adding a public benchmark-only API to catallax-step.
#[path = "../src/fit.rs"]
mod fit;
#[path = "../src/number.rs"]
#[expect(
    dead_code,
    reason = "The fitting benchmark uses only part of the shared numeric module; the library checks all of it"
)]
mod number;

fn vector(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(x.into(), y.into(), z.into())
}

fn point(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(x.into(), y.into(), z.into())
}

fn cylinder_surface(radius: i64, offset: i64, crossing: bool) -> Surface {
    let (origin, axis, reference) = if crossing {
        (point(0, offset, 0), vector(1, 0, 0), vector(0, 0, 1))
    } else {
        (RPoint3::origin(), vector(0, 0, 1), vector(1, 0, 0))
    };
    Surface::Cylinder(Cylinder::new(origin, axis, radius.into(), reference).unwrap())
}

fn sphere_surface(x: i64, radius: i64) -> Surface {
    Surface::Sphere(
        Sphere::new(
            point(x, 0, 0),
            radius.into(),
            vector(0, 0, 1),
            vector(1, 0, 0),
        )
        .unwrap(),
    )
}

fn primitives() {
    let mut model = Model::new();
    let pose = Pose::identity();
    let bodies = [
        (
            "block",
            block(&mut model, [2.into(), 3.into(), 5.into()], &pose).unwrap(),
        ),
        (
            "cylinder",
            cylinder(&mut model, 2.into(), 5.into(), &pose).unwrap(),
        ),
        ("sphere", sphere(&mut model, 5.into(), &pose).unwrap()),
        (
            "cone",
            cone(&mut model, 3.into(), 0.into(), 4.into(), &pose).unwrap(),
        ),
        (
            "torus",
            torus(&mut model, 3.into(), 1.into(), &pose).unwrap(),
        ),
    ];
    println!("primitive,copies,faces,bytes,microseconds");
    for (name, body) in bodies {
        for count in [1, 16, 128] {
            let bodies = vec![body; count];
            let mut out = Vec::new();
            let started = Instant::now();
            let answer = write_step(
                &model,
                &bodies,
                &StepOptions::default(),
                &mut Decider::default(),
                &mut out,
            )
            .unwrap();
            let elapsed = started.elapsed();
            println!(
                "{name},{count},{},{},{}",
                answer.value.faces,
                out.len(),
                elapsed.as_micros()
            );
            black_box(out);
        }
    }
}

fn branches() {
    let cases = [
        (
            17,
            cylinder_surface(2, 0, false),
            cylinder_surface(1, 0, true),
        ),
        (
            18,
            cylinder_surface(2, 0, false),
            cylinder_surface(1, 2, true),
        ),
        (
            19,
            cylinder_surface(2, 0, false),
            cylinder_surface(1, 1, true),
        ),
        (
            20,
            cylinder_surface(2, 0, false),
            // Case 20 is a point at offset 3; its specified near-tangent
            // variant supplies the branch that this fitting benchmark needs.
            Surface::Cylinder(
                Cylinder::new(
                    RPoint3::new(
                        0.into(),
                        Rational::from_i64(3) - Rational::parse("1e-30").unwrap(),
                        0.into(),
                    ),
                    vector(1, 0, 0),
                    1.into(),
                    vector(0, 0, 1),
                )
                .unwrap(),
            ),
        ),
        (
            21,
            cylinder_surface(1, 0, false),
            cylinder_surface(1, 0, true),
        ),
        (
            22,
            cylinder_surface(1, 0, false),
            Surface::Cylinder(
                Cylinder::new(
                    RPoint3::origin(),
                    vector(3, 0, 4),
                    1.into(),
                    vector(0, 1, 0),
                )
                .unwrap(),
            ),
        ),
        (
            23,
            cylinder_surface(5, 0, false),
            Surface::Sphere(
                Sphere::new(point(0, 0, 1), 6.into(), vector(0, 0, 1), vector(1, 0, 0)).unwrap(),
            ),
        ),
        (24, cylinder_surface(1, 0, false), sphere_surface(1, 2)),
        (25, cylinder_surface(1, 0, false), sphere_surface(2, 2)),
        (26, cylinder_surface(1, 0, false), sphere_surface(1, 3)),
    ];
    println!(
        "case,component,tolerance,control_points,spans,trace_microseconds,fit_microseconds,max_deviation,status"
    );
    let selected = std::env::args().nth(1).and_then(|s| s.parse::<i32>().ok());
    for (case, a, b) in cases {
        if selected.is_some_and(|selected| selected != case) {
            continue;
        }
        let mut geom = GeomStore::new();
        let a = geom.add_surface(a);
        let b = geom.add_surface(b);
        let SurfacePair::Meets(components) =
            intersect_surfaces(&mut geom, a, b, &mut Decider::default())
                .unwrap()
                .value
        else {
            continue;
        };
        for (index, component) in components.iter().enumerate() {
            let Component::Curve { curve, .. } = component else {
                println!("{case},{index},-,0,0,0,0,0,point");
                continue;
            };
            if !matches!(
                curve_form(&geom, *curve, &mut Decider::default()).unwrap(),
                CurveForm::Branch
            ) {
                println!("{case},{index},-,0,0,0,0,0,analytic");
                continue;
            }
            for exponent in 6..=12 {
                let options = StepOptions {
                    curve_tolerance: Rational::parse(&format!("1e-{exponent}")).unwrap(),
                    ..StepOptions::default()
                };
                let mut dec = Decider::default();
                let start = geom.curve(*curve).parameter_origin(&geom);
                let start = number::point(&start, options.significant_digits, &mut dec).unwrap();
                let began = Instant::now();
                let result = trace(&geom, *curve, None, None, &mut dec);
                let trace_time = began.elapsed();
                let trace = match result {
                    Ok(answer) => answer.value,
                    Err(refusal) => {
                        println!(
                            "{case},{index},1e-{exponent},0,0,{},0,0,{}",
                            trace_time.as_micros(),
                            refusal.code
                        );
                        continue;
                    }
                };
                let began = Instant::now();
                let fit = fit::fit(&trace, &start, &start, false, &options, &dec).unwrap();
                let elapsed = began.elapsed();
                println!(
                    "{case},{index},1e-{exponent},{},{},{},{},{},{}",
                    fit.points.len(),
                    fit.knots.len() - 1,
                    trace_time.as_micros(),
                    elapsed.as_micros(),
                    fit.error,
                    if fit.error <= options.curve_tolerance {
                        "within-tolerance"
                    } else {
                        "marked"
                    }
                );
                black_box((fit.degree, fit.multiplicities, fit.bounds, fit.points));
            }
        }
    }
}

fn main() {
    primitives();
    branches();
}
