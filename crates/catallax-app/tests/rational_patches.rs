//! Public-kernel carrier consumers. Loft and sweep recipe entry are separate work.
mod rational_patch_cases;

use catallax_boolean::Containment;
use catallax_measure::{Accuracy, Enclosure};
use catallax_mesh::{Mesh, MeshOptions};
use catallax_numbers::{Ball, BoundKind, Decider};
use catallax_shapes::{Isometry, RPoint3, RVector3};
use catallax_topology::RationalPatchFace;
use rational_patch_cases::{Roof, p, r, solid};
use std::collections::{BTreeMap, BTreeSet};

fn accuracy() -> Accuracy {
    Accuracy {
        relative: r("0"),
        absolute: r("1e-6"),
    }
}
fn contains(actual: &Enclosure, expected: &str) {
    assert!(
        actual.contains(&r(expected)),
        "{actual:?} excludes {expected}"
    );
    assert!(
        actual.width() <= accuracy().absolute,
        "unmet requested accuracy: {actual:?}"
    );
}

#[test]
fn complete_named_faces_have_exact_witnesses_and_finite_membership() {
    for kind in [
        Roof::Bilinear,
        Roof::Rational,
        Roof::NonuniformBump,
        Roof::JoinedDegreeElevated,
        Roof::JoinedNonlinear,
    ] {
        let mut fixture = solid(kind);
        assert!(fixture.model.check(fixture.body).is_empty());
        let mut dec = Decider::default();
        for &face in &fixture.roofs {
            assert!(
                RationalPatchFace::of_face(&fixture.model, face, &mut dec)
                    .unwrap()
                    .is_some()
            );
            let witness = catallax_boolean::face_point(&mut fixture.model, face, &mut dec)
                .unwrap()
                .value;
            assert!(
                catallax_boolean::face_contains(&mut fixture.model, face, &witness, &mut dec)
                    .unwrap()
            );
            let mut outside = witness;
            outside.z = &outside.z + &r("1e-20").into();
            assert!(
                !catallax_boolean::face_contains(&mut fixture.model, face, &outside, &mut dec)
                    .unwrap()
            );
        }
        let inside = catallax_boolean::classify_point(
            &mut fixture.model,
            fixture.body,
            &p("1/2", "1/2", "1/2"),
            &mut dec,
        );
        assert_eq!(
            inside.proven(),
            Some(&Containment::Inside),
            "{kind:?}: {inside:?}"
        );
        let outside = catallax_boolean::classify_point(
            &mut fixture.model,
            fixture.body,
            &p("1/2", "1/2", "3"),
            &mut dec,
        );
        assert_eq!(
            outside.proven(),
            Some(&Containment::Outside),
            "{kind:?}: {outside:?}"
        );
    }
}

fn exact_area(bump: bool) -> Ball {
    let c = |value| Ball::from_rational(&catallax_numbers::Rational::from_i64(value), 160);
    let root_two = c(2).sqrt().unwrap();
    if bump {
        let logarithm = (&c(1) + &root_two).ln().unwrap();
        &c(16).checked_div(&c(3)).unwrap() + &(&root_two + &logarithm).checked_div(&c(2)).unwrap()
    } else {
        let root_three = c(3).sqrt().unwrap();
        let logarithm = (&c(1) + &root_three)
            .checked_div(&root_two)
            .unwrap()
            .ln()
            .unwrap();
        let numerator =
            &(&root_three + &(&c(4) * &logarithm)) - &Ball::pi(160).checked_div(&c(6)).unwrap();
        &c(6) + &numerator.checked_div(&c(3)).unwrap()
    }
}

fn contains_area(actual: &Enclosure, bump: bool) {
    let exact = exact_area(bump);
    assert!(
        actual.lo <= exact.lo() && actual.hi >= exact.hi(),
        "{actual:?} excludes the closed-form area"
    );
}

#[test]
fn tight_bilinear_area_keeps_truth_and_reports_the_unmet_width() {
    let fixture = solid(Roof::Bilinear);
    let requested = Accuracy {
        relative: r("0"),
        absolute: r("1e-7"),
    };
    let answer = catallax_measure::mass_properties(
        &fixture.model,
        fixture.body,
        &requested,
        &mut Decider::default(),
    )
    .unwrap();
    contains_area(&answer.value.area, false);
    let width = answer.value.area.width();
    assert!(width > requested.absolute);
    let mark = answer
        .marks
        .iter()
        .find(|mark| mark.code == "measure.accuracy-not-reached")
        .expect("an unmet area request must remain marked");
    assert_eq!(mark.kind, BoundKind::Value);
    assert_eq!(mark.bound, width);
}

#[test]
fn rational_parameterisation_preserves_mass_moments_and_geometric_bounds() {
    for kind in [Roof::Bilinear, Roof::Rational, Roof::NonuniformBump] {
        let fixture = solid(kind);
        let bump = matches!(kind, Roof::NonuniformBump);
        let mut dec = Decider::default();
        let mut requested = accuracy();
        if bump {
            requested.absolute = r("1e-2");
        }
        let answer =
            catallax_measure::mass_properties(&fixture.model, fixture.body, &requested, &mut dec)
                .unwrap();
        assert!(
            answer.marks.is_empty(),
            "{kind:?}: area width {:.6e}, {:?}",
            answer.value.area.width().to_f64_nearest(),
            answer.marks
        );
        let mass = answer.value;
        contains_area(&mass.area, bump);
        assert!(mass.area.width() <= requested.absolute);
        contains(&mass.volume, if bump { "7/6" } else { "5/4" });
        let center = mass.centroid.unwrap();
        for (actual, expected) in center.iter().zip(if bump {
            ["1/2", "1/2", "41/70"]
        } else {
            ["8/15", "8/15", "29/45"]
        }) {
            contains(actual, expected);
        }
        let expected = if bump {
            [
                ["733/3150", "0", "0"],
                ["0", "159/700", "0"],
                ["0", "0", "17/90"],
            ]
        } else {
            [
                ["1937/6480", "-1/180", "-31/1080"],
                ["-1/180", "1937/6480", "-31/1080"],
                ["-31/1080", "-31/1080", "37/180"],
            ]
        };
        let inertia = mass.inertia.unwrap();
        for i in 0..3 {
            for j in 0..3 {
                contains(&inertia[i][j], expected[i][j]);
            }
        }
        let bounds =
            catallax_measure::bounding_box(&fixture.model, fixture.body, &accuracy(), &mut dec)
                .unwrap();
        assert!(bounds.marks.is_empty());
        let bounds = bounds.value.unwrap();
        for axis in 0..3 {
            let max = if axis == 2 {
                if bump { r("5/4") } else { r("2") }
            } else {
                r("1")
            };
            assert!(bounds.min[axis] <= r("0") && bounds.min[axis] >= -accuracy().absolute);
            assert!(bounds.max[axis] >= max && &bounds.max[axis] - &max <= accuracy().absolute);
        }
    }
}

fn mesh_proof(mesh: &Mesh, options: &MeshOptions, expected_volume: f64) {
    let mut incidence = BTreeMap::<_, (usize, i32)>::new();
    let mut volume = 0.0;
    for triangle in &mesh.triangles {
        let [a, b, c] = triangle.vertices.map(|v| mesh.positions[v]);
        volume += (a[0] * (b[1] * c[2] - b[2] * c[1])
            + a[1] * (b[2] * c[0] - b[0] * c[2])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))
            / 6.0;
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let count = incidence.entry((a.min(b), a.max(b))).or_default();
            count.0 += 1;
            count.1 += if a < b { 1 } else { -1 };
        }
    }
    assert!(
        incidence
            .values()
            .all(|&(uses, winding)| uses == 2 && winding == 0)
    );
    assert!(
        (volume - expected_volume).abs() < 0.02,
        "oriented mesh volume {volume}"
    );
    for face in &mesh.faces {
        assert!(face.max_chord_error <= options.chord_tolerance);
        assert!(
            face.max_angular_deviation
                .as_ref()
                .is_some_and(|bound| bound <= options.angular_deviation.as_ref().unwrap())
        );
    }
    assert!(
        mesh.unavailable_controls.is_empty(),
        "{:?}",
        mesh.unavailable_controls
    );
}

#[test]
fn shared_rows_sew_even_with_degree_changes_and_nonlinear_parameters() {
    let options = MeshOptions {
        chord_tolerance: r("1/50"),
        angular_deviation: Some(r("1/3")),
        max_edge_length: None,
        max_aspect_ratio: None,
        max_vertices: 200_000,
    };
    for kind in [
        Roof::Bilinear,
        Roof::Rational,
        Roof::NonuniformBump,
        Roof::JoinedDegreeElevated,
        Roof::JoinedNonlinear,
    ] {
        let fixture = solid(kind);
        let mesh = catallax_mesh::mesh_body(
            &fixture.model,
            fixture.body,
            &options,
            &mut Decider::default(),
        )
        .unwrap_or_else(|error| panic!("{kind:?}: {error:?}"))
        .value;
        let expected: BTreeSet<_> = fixture
            .model
            .faces_of_body(fixture.body)
            .into_iter()
            .collect();
        assert_eq!(
            mesh.faces
                .iter()
                .map(|face| face.face)
                .collect::<BTreeSet<_>>(),
            expected
        );
        mesh_proof(
            &mesh,
            &options,
            if matches!(kind, Roof::NonuniformBump) {
                7.0 / 6.0
            } else {
                1.25
            },
        );
    }
}

#[test]
fn step_rounding_preserves_the_source_equations_and_shared_row_charts() {
    for kind in [
        Roof::Bilinear,
        Roof::Rational,
        Roof::NonuniformBump,
        Roof::JoinedDegreeElevated,
    ] {
        let fixture = solid(kind);
        let mut bytes = Vec::new();
        let mut dec = Decider::default();
        let options = catallax_step::StepOptions::default();
        let report = catallax_step::write_step(
            &fixture.model,
            &[fixture.body],
            &options,
            &mut dec,
            &mut bytes,
        )
        .unwrap()
        .value;
        assert!(report.max_deviation <= options.curve_tolerance);
        let parsed = catallax_import::read_step(&bytes).unwrap();
        assert_eq!(parsed.solids.len(), 1);
        assert_eq!(parsed.solids[0].faces.len(), report.faces);
        let mut patches = 0;
        let half = r("1/2");
        let one = Ball::from_rational(&r("1"), 160);
        let allowed = &report.max_deviation * &r("5") + r("1e-25");
        for face in &parsed.solids[0].faces {
            let catallax_import::Surface::BSpline(surface) = &face.surface else {
                continue;
            };
            patches += 1;
            let surface = surface.geometry.to_surface(&mut dec).unwrap();
            let point = surface
                .evaluate_box([&half, &half], [&half, &half], 160, &dec)
                .unwrap();
            let exact_z = if matches!(kind, Roof::NonuniformBump) {
                &one + &(&point[0] * &(&one - &point[0]))
            } else {
                &one + &(&point[0] * &point[1])
            };
            let residual = &point[2] - &exact_z;
            assert!(residual.lo() >= -allowed.clone() && residual.hi() <= allowed);
        }
        assert_eq!(patches, fixture.patches.len());
    }
    let fixture = solid(Roof::JoinedNonlinear);
    let mut bytes = Vec::new();
    let refusal = catallax_step::write_step(
        &fixture.model,
        &[fixture.body],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap_err();
    assert_eq!(refusal.code, "step.rational-patch-unproven");
    assert!(bytes.is_empty());
}

#[test]
fn warmed_definition_caches_do_not_leak_through_reflection() {
    let mut fixture = solid(Roof::Rational);
    let mut dec = Decider::default();
    let original =
        catallax_measure::mass_properties(&fixture.model, fixture.body, &accuracy(), &mut dec)
            .unwrap()
            .value;
    let names = |model: &catallax_topology::Model, body| {
        model
            .faces_of_body(body)
            .into_iter()
            .map(|face| {
                model
                    .face_origins(face)
                    .iter()
                    .map(|origin| origin.key())
                    .collect::<Vec<_>>()
            })
            .collect::<BTreeSet<_>>()
    };
    let before = names(&fixture.model, fixture.body);
    let map = Isometry::reflection(
        &RPoint3::new(r("2"), r("0"), r("0")),
        &RVector3::new(r("1"), r("0"), r("0")),
    )
    .unwrap();
    let copy =
        catallax_topology::transformed_body(&mut fixture.model, fixture.body, &map, Clone::clone);
    assert_eq!(before, names(&fixture.model, copy));
    assert!(fixture.model.check(copy).is_empty());
    let transformed =
        catallax_measure::mass_properties(&fixture.model, copy, &accuracy(), &mut dec).unwrap();
    assert!(transformed.marks.is_empty());
    contains(&transformed.value.volume, "5/4");
    for (actual, expected) in transformed
        .value
        .centroid
        .unwrap()
        .iter()
        .zip(["52/15", "8/15", "29/45"])
    {
        contains(actual, expected);
    }
    contains(&original.volume, "5/4");
}
