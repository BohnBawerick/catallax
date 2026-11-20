//! Finite cylindrical pcurves with source vertices slightly off their carriers.
#[path = "../../catallax-topology/tests/common/analytic.rs"]
mod fixture;

use catallax_numbers::{Ball, Decider, Real};
use catallax_shapes::{Curve, Point3, RVector3, Surface, Vector3};
use catallax_step::numerical::{self as output, NumericalStepOptions};
use catallax_topology::numerical as n;
use fixture::q;

fn rotated(mut d: n::NumericalDefinition) -> n::NumericalDefinition {
    let c = q("20/29");
    let s = q("21/29");
    let rotate = |x: &Real, y: &Real| {
        [
            &(x * &Real::from(c.clone())) - &(y * &Real::from(s.clone())),
            &(x * &Real::from(s.clone())) + &(y * &Real::from(c.clone())),
        ]
    };
    let point = |p: &Point3| {
        let [x, y] = rotate(&p.x, &p.y);
        Point3::new(x, y, p.z.clone())
    };
    let vector = |v: &Vector3| {
        let [x, y] = rotate(&v.x, &v.y);
        Vector3::new(x, y, v.z.clone())
    };
    for vertex in &mut d.vertices {
        let [x, y, _] = vertex.position.clone();
        vertex.position[0] = &(&x * &c) - &(&y * &s);
        vertex.position[1] = &(&x * &s) + &(&y * &c);
    }
    for edge in &mut d.edges {
        match &mut edge.curve {
            Curve::Line(line) => {
                line.origin = point(&line.origin);
                line.dir = vector(&line.dir);
            }
            Curve::Circle(circle) => {
                circle.center = point(&circle.center);
                circle.normal = vector(&circle.normal);
                circle.ref_dir = vector(&circle.ref_dir);
            }
            _ => unreachable!(),
        }
    }
    for face in &mut d.faces {
        match &mut face.surface {
            Surface::Plane(plane) => {
                plane.origin = point(&plane.origin);
                plane.normal = vector(&plane.normal);
                plane.ref_dir = vector(&plane.ref_dir);
            }
            Surface::Cylinder(cylinder) => {
                cylinder.origin = point(&cylinder.origin);
                cylinder.axis = vector(&cylinder.axis);
                cylinder.ref_dir = vector(&cylinder.ref_dir);
            }
            _ => unreachable!(),
        }
    }
    d
}

fn chart(phase: usize) -> n::NumericalDefinition {
    let mut d = fixture::quarter();
    let surface = d.faces.last_mut().unwrap();
    let Surface::Cylinder(cylinder) = &mut surface.surface else {
        unreachable!()
    };
    let reference = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]][phase];
    cylinder.ref_dir = RVector3::new(
        reference[0].into(),
        reference[1].into(),
        reference[2].into(),
    )
    .to_real();
    let offset = &(Ball::pi(128).mid() / q("2"))
        * &catallax_numbers::Rational::from_i64(phase.try_into().unwrap());
    for uv in surface.loops.iter_mut().flatten().flat_map(|t| &mut t.uv) {
        uv[0] = &uv[0] - &offset;
    }
    d
}

#[test]
fn finite_cylinder_exchange_with_toleranced_source_vertex_incidence() {
    let displacements = ["-1e-10", "-1e-12", "-1e-14", "0", "1e-14", "1e-12", "1e-10"];
    for oblique in [false, true] {
        for phase in 0..4 {
            for vertex in [1, 2, 4, 5] {
                for (noise, displacement) in displacements.iter().enumerate() {
                    let mut definition = chart(phase);
                    if oblique {
                        definition = rotated(definition);
                    }
                    definition.vertices[vertex].position[2] =
                        &definition.vertices[vertex].position[2] + &q(displacement);
                    let body = n::build(
                        definition,
                        n::NumericalOptions::default(),
                        &mut Decider::default(),
                    )
                    .unwrap();
                    let file = output::build(
                        &body,
                        &NumericalStepOptions::default(),
                        &mut Decider::default(),
                    )
                    .unwrap();
                    file.check_written(file.bytes(), &mut Decider::default())
                        .unwrap();
                    assert_eq!(file.report().counts, [5, 9, 6]);
                    if let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE") {
                        let name = format!(
                            "incidence-oblique-{oblique}-phase-{phase}-vertex-{vertex}-noise-{noise}.step"
                        );
                        std::fs::write(std::path::Path::new(&path).join(name), file.bytes())
                            .unwrap();
                    }
                }
            }
        }
    }
}

#[test]
fn source_vertex_displacement_above_tolerance_is_refused() {
    for displacement in ["-0.02", "0.02"] {
        let mut definition = chart(1);
        definition.vertices[1].position[2] = q(displacement);
        assert!(
            n::build(
                definition,
                n::NumericalOptions::default(),
                &mut Decider::default(),
            )
            .is_err()
        );
    }
}

#[test]
fn tiny_tilted_carriers_keep_existing_flat_pcurves() {
    for oblique in [false, true] {
        for phase in 0..4 {
            for (noise, tilt) in ["-1e-14", "-1e-16", "1e-16", "1e-14"].iter().enumerate() {
                let mut definition = chart(phase);
                for edge in &mut definition.edges {
                    if let Curve::Circle(circle) = &mut edge.curve {
                        circle.normal = RVector3::new(q("0"), q(tilt), q("1")).to_real();
                    }
                }
                let face = definition.faces.last_mut().unwrap();
                let Surface::Cylinder(cylinder) = &mut face.surface else {
                    unreachable!()
                };
                let offset = q("37.3578654636477");
                cylinder.axis = RVector3::new(q("0"), q("0"), q("-1")).to_real();
                cylinder.origin.z = Real::from(-&offset);
                for uv in face
                    .loops
                    .iter_mut()
                    .flatten()
                    .flat_map(|trim| &mut trim.uv)
                {
                    uv[0] = -&uv[0];
                    uv[1] = -(&uv[1] + &offset);
                }
                if oblique {
                    definition = rotated(definition);
                }
                let body = n::build(
                    definition,
                    n::NumericalOptions::default(),
                    &mut Decider::default(),
                )
                .unwrap();
                let file = output::build(
                    &body,
                    &NumericalStepOptions::default(),
                    &mut Decider::default(),
                )
                .unwrap();
                file.check_written(file.bytes(), &mut Decider::default())
                    .unwrap();
                assert_eq!(file.report().counts, [5, 9, 6]);
                if let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE") {
                    let name = format!("tilted-oblique-{oblique}-phase-{phase}-noise-{noise}.step");
                    std::fs::write(std::path::Path::new(&path).join(name), file.bytes()).unwrap();
                }
            }
        }
    }
}

#[test]
fn circle_tilt_above_source_tolerance_is_refused() {
    for tilt in ["-0.02", "0.02"] {
        let mut definition = chart(1);
        for edge in &mut definition.edges {
            if let Curve::Circle(circle) = &mut edge.curve {
                circle.normal = RVector3::new(q("0"), q(tilt), q("1")).to_real();
            }
        }
        assert!(
            n::build(
                definition,
                n::NumericalOptions::default(),
                &mut Decider::default(),
            )
            .is_err()
        );
    }
}

// Rotate the cylinder axis away from every coordinate axis. The earlier
// XY-rotated controls keep their cylinder axis on Z and cannot pin this case.
fn oblique_axis(mut d: n::NumericalDefinition) -> n::NumericalDefinition {
    let c = q("20/29");
    let s = q("21/29");
    let rotate = |x: &Real, z: &Real| {
        [
            &(x * &Real::from(c.clone())) - &(z * &Real::from(s.clone())),
            &(x * &Real::from(s.clone())) + &(z * &Real::from(c.clone())),
        ]
    };
    let point = |p: &Point3| {
        let [x, z] = rotate(&p.x, &p.z);
        Point3::new(x, p.y.clone(), z)
    };
    let vector = |p: &Vector3| {
        let [x, z] = rotate(&p.x, &p.z);
        Vector3::new(x, p.y.clone(), z)
    };
    for vertex in &mut d.vertices {
        let [x, _, z] = vertex.position.clone();
        vertex.position[0] = &(&x * &c) - &(&z * &s);
        vertex.position[2] = &(&x * &s) + &(&z * &c);
    }
    for edge in &mut d.edges {
        match &mut edge.curve {
            Curve::Line(line) => {
                line.origin = point(&line.origin);
                line.dir = vector(&line.dir);
            }
            Curve::Circle(circle) => {
                circle.center = point(&circle.center);
                circle.normal = vector(&circle.normal);
                circle.ref_dir = vector(&circle.ref_dir);
            }
            _ => unreachable!(),
        }
    }
    for face in &mut d.faces {
        match &mut face.surface {
            Surface::Plane(plane) => {
                plane.origin = point(&plane.origin);
                plane.normal = vector(&plane.normal);
                plane.ref_dir = vector(&plane.ref_dir);
            }
            Surface::Cylinder(cylinder) => {
                cylinder.origin = point(&cylinder.origin);
                cylinder.axis = vector(&cylinder.axis);
                cylinder.ref_dir = vector(&cylinder.ref_dir);
            }
            _ => unreachable!(),
        }
    }
    d
}

#[test]
fn emitted_circle_pcurves_keep_constant_height_on_oblique_parallel_and_antiparallel_axes() {
    // Follow references in the actual emitted STEP, rather than inspecting the
    // writer's internal endpoints or inferring the claim from BRep validity.
    fn references(entity: &str) -> Vec<usize> {
        entity
            .split('#')
            .skip(1)
            .map(|s| {
                s[..s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len())]
                    .parse()
                    .unwrap()
            })
            .collect()
    }
    fn pair(entity: &str) -> [catallax_numbers::Rational; 2] {
        let (_, values) = entity.split_once(",(").unwrap();
        let (values, _) = values.split_once(')').unwrap();
        let values: Vec<_> = values.split(',').map(step_number).collect();
        values.try_into().unwrap()
    }
    fn step_number(value: &str) -> catallax_numbers::Rational {
        // STEP permits 0. and 1.E+2; the general rational parser requires a
        // digit after the decimal point. Normalize this emitted syntax only.
        let mut normalized = value.replace(".E", ".0E").replace(".e", ".0e");
        if normalized.ends_with('.') {
            normalized.push('0');
        }
        q(&normalized)
    }
    for antiparallel in [false, true] {
        let mut definition = chart(1);
        for index in [1, 4] {
            definition.vertices[index].position[2] =
                &definition.vertices[index].position[2] + &q("1e-10");
        }
        let offset = q("37.3578654636477");
        let face = definition.faces.last_mut().unwrap();
        let Surface::Cylinder(cylinder) = &mut face.surface else {
            unreachable!()
        };
        cylinder.origin.z = Real::from(-&offset);
        if antiparallel {
            cylinder.axis.z = Real::from(-1);
        }
        for uv in face
            .loops
            .iter_mut()
            .flatten()
            .flat_map(|trim| &mut trim.uv)
        {
            uv[1] = &uv[1] + &offset;
            if antiparallel {
                uv[0] = -&uv[0];
                uv[1] = -&uv[1];
            }
        }
        let body = n::build(
            oblique_axis(definition),
            n::NumericalOptions::default(),
            &mut Decider::default(),
        )
        .unwrap();
        let file = output::build(
            &body,
            &NumericalStepOptions::default(),
            &mut Decider::default(),
        )
        .unwrap();
        file.check_written(file.bytes(), &mut Decider::default())
            .unwrap();
        let text = std::str::from_utf8(file.bytes()).unwrap();
        let entities: std::collections::BTreeMap<usize, &str> = text
            .lines()
            .filter_map(|line| {
                let (id, value) = line.strip_prefix('#')?.split_once('=')?;
                Some((id.parse().unwrap(), value))
            })
            .collect();
        let mut heights = Vec::new();
        for entity in entities
            .values()
            .filter(|e| e.starts_with("SURFACE_CURVE("))
        {
            let refs = references(entity);
            if !entities[&refs[0]].starts_with("CIRCLE(") {
                continue;
            }
            for pcurve in refs.iter().skip(1).map(|id| entities[id]) {
                let refs = references(pcurve);
                if !entities[&refs[0]].starts_with("CYLINDRICAL_SURFACE(") {
                    continue;
                }
                let placement = references(entities[&refs[0]])[0];
                let axis = references(entities[&placement])[1];
                let (_, components) = entities[&axis].split_once(",(").unwrap();
                let components: Vec<_> = components
                    .split(')')
                    .next()
                    .unwrap()
                    .split(',')
                    .map(step_number)
                    .collect();
                assert!(!components[0].is_zero() && !components[2].is_zero());
                let line = references(entities[&refs[1]])[0];
                let refs = references(entities[&line]);
                let origin = pair(entities[&refs[0]]);
                let direction = references(entities[&refs[1]])[0];
                let direction = pair(entities[&direction]);
                assert!(!direction[0].is_zero());
                assert!(
                    direction[1].is_zero(),
                    "emitted pcurve V slope: {direction:?}"
                );
                heights.push(if antiparallel {
                    -&origin[1]
                } else {
                    origin[1].clone()
                });
            }
        }
        heights.sort();
        assert_eq!(heights.len(), 2);
        for (actual, expected) in heights.iter().zip([offset.clone(), &offset + &q("1")]) {
            assert!(
                (actual - &expected).abs() < q("1e-10"),
                "{actual} != {expected}"
            );
        }
        if let Ok(path) = std::env::var("CATALLAX_OBLIQUE_HEIGHT_EVIDENCE") {
            std::fs::write(
                std::path::Path::new(&path)
                    .join(format!("oblique-axis-antiparallel-{antiparallel}.step")),
                file.bytes(),
            )
            .unwrap();
        }
    }
}

#[test]
fn tilted_circle_skip_path_retains_its_checked_endpoint_height_slope() {
    use catallax_numbers::NumericalCheck;
    let mut definition = chart(1);
    for index in [1, 4] {
        definition.vertices[index].position[2] =
            &definition.vertices[index].position[2] + &q("1e-10");
    }
    for edge in &mut definition.edges {
        if let Curve::Circle(circle) = &mut edge.curve {
            circle.normal.y = Real::from(q("1e-12"));
        }
    }
    let circle_keys: Vec<_> = definition
        .edges
        .iter()
        .filter(|e| matches!(e.curve, Curve::Circle(_)))
        .map(|e| e.key.clone())
        .collect();
    // The two end heights of each circle along the cylinder axis, read off
    // its vertices before the rigid rotation, which keeps them. The emitted
    // admission record holds one measured sample, not the two ends.
    let end_heights: Vec<_> = definition
        .edges
        .iter()
        .filter(|e| matches!(e.curve, Curve::Circle(_)))
        .map(|e| {
            let mut heights = e
                .vertices
                .map(|v| definition.vertices[v].position[2].clone());
            heights.sort();
            heights
        })
        .collect();
    let face_key = definition.faces.last().unwrap().key.clone();
    let body = n::build(
        oblique_axis(definition),
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let file = output::build(
        &body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    file.check_written(file.bytes(), &mut Decider::default())
        .unwrap();
    for (key, original) in circle_keys.iter().zip(&end_heights) {
        let checks = file.report().emitted_evidence.checks();
        let actual = checks
            .iter()
            .find_map(|check| match check {
                NumericalCheck::Samples {
                    name,
                    scope,
                    samples,
                } if name == "file.actual-decimal-pcurve-edge-surface-distance"
                    && *scope == format!("{key} on {face_key}") =>
                {
                    Some(samples)
                }
                _ => None,
            })
            .unwrap();
        let heights = |samples: &[catallax_numbers::NumericalSample]| {
            let mut values = [
                samples.first().unwrap().parameters[2].clone(),
                samples.last().unwrap().parameters[2].clone(),
            ];
            values.sort();
            values
        };
        let written = heights(actual);
        assert_ne!(
            original[0], original[1],
            "control did not exercise a varying-V map"
        );
        assert_ne!(
            written[0], written[1],
            "tilted skip path flattened the pcurve"
        );
        for (before, after) in original.iter().zip(&written) {
            assert!((before - after).abs() < q("1e-12"));
        }
    }
}
