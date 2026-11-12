//! Imported definitions answer field queries and requested two-kernel checks.

#[path = "../../catallax-import/tests/common/mod.rs"]
mod common;

use catallax_compare::{Controller, Options, Verdict};
use catallax_features::ImportFiles;
use catallax_field::{Direction, FaceState, Location, Plan, Workspace};
use catallax_numbers::{Ball, Decider, Rational};
use catallax_recipe::{FileRef, Operation, Pose, Recipe, Step, StepId};
use catallax_shapes::{Pose as ShapePose, RPoint3};
use catallax_topology::{Model, cone, cylinder, sphere, torus};
use std::collections::BTreeSet;

fn q(n: i64) -> Rational {
    n.into()
}
fn point(v: [i64; 3]) -> RPoint3 {
    let [x, y, z] = v.map(q);
    RPoint3::new(x, y, z)
}
fn recipe(bytes: &[u8], pose: Pose) -> Recipe {
    let id = StepId::new("imported").unwrap();
    Recipe::new(
        "imported field",
        vec![Step {
            id: id.clone(),
            op: Operation::ImportStep {
                file: FileRef {
                    name: "part.step".into(),
                    sha256: catallax_import::sha256_hex(bytes),
                },
                pose,
                decisions: Default::default(),
            },
        }],
        id,
    )
}
fn sources(bytes: &[u8]) -> ImportFiles<'_> {
    ImportFiles::from([(catallax_import::sha256_hex(bytes), bytes)])
}
fn exported(kind: &str) -> Vec<u8> {
    let mut model = Model::new();
    let pose = ShapePose::identity();
    let body = match kind {
        "cylinder" => cylinder(&mut model, 3.into(), 7.into(), &pose),
        "cone" => cone(&mut model, 3.into(), 1.into(), 2.into(), &pose),
        "sphere" => sphere(&mut model, 3.into(), &pose),
        "torus" => torus(&mut model, 4.into(), 1.into(), &pose),
        _ => unreachable!(),
    }
    .unwrap();
    let mut bytes = Vec::new();
    catallax_step::write_step(
        &model,
        &[body],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    if kind == "cone" {
        // The current importer accepts exactly 45 degrees in a degree context.
        let text = String::from_utf8(bytes).unwrap();
        let radian = text
            .lines()
            .find(|line| line.contains("PLANE_ANGLE_UNIT()") && line.contains(".RADIAN."))
            .unwrap()
            .split_once('=')
            .unwrap()
            .0;
        let mut changed = String::new();
        for line in text.lines() {
            if line.contains("CONICAL_SURFACE") {
                let (prefix, _) = line.rsplit_once(',').unwrap();
                changed.push_str(&format!("{prefix},45.);\n"));
            } else {
                if line.contains("GLOBAL_UNIT_ASSIGNED_CONTEXT") {
                    changed.push_str(&line.replace(&format!("{radian},"), "#999999,"));
                } else {
                    changed.push_str(line);
                }
                changed.push('\n');
            }
        }
        return changed
            .replace(
                "ENDSEC;\nEND-ISO-10303-21;",
                &format!("#999997=DIMENSIONAL_EXPONENTS(0.,0.,0.,0.,0.,0.,0.);\n#999998=PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.017453292519943295),{radian});\n#999999=(CONVERSION_BASED_UNIT('DEGREE',#999998)NAMED_UNIT(#999997)PLANE_ANGLE_UNIT());\nENDSEC;\nEND-ISO-10303-21;"),
            )
            .into_bytes();
    }
    bytes
}
#[test]
fn imported_box_cylinder_and_cone_keep_names_and_agree_on_request() {
    for (kind, bytes, inside, wall, volume_factor) in [
        (
            "box",
            common::box_step(common::BoxFixture::new(["2", "3", "4"])).into_bytes(),
            [1, 1, 1],
            [0, 1, 1],
            q(24),
        ),
        (
            "cylinder",
            exported("cylinder"),
            [0, 0, 3],
            [3, 0, 3],
            q(63),
        ),
        (
            "cone",
            exported("cone"),
            [0, 0, 1],
            [2, 0, 1],
            Rational::parse("26/3").unwrap(),
        ),
    ] {
        let recipe = recipe(&bytes, Pose::identity());
        let files = sources(&bytes);
        let field =
            Plan::compile_with_files(&recipe, &files).unwrap_or_else(|e| panic!("{kind}: {e:?}"));
        let built =
            catallax_features::evaluate_with_files(&recipe, &files, &mut Decider::default())
                .unwrap();
        let brep_names: BTreeSet<_> = built
            .model
            .faces_of_body(built.result)
            .iter()
            .flat_map(|id| built.model.face_origins(*id))
            .cloned()
            .collect();
        let field_names: BTreeSet<_> = field
            .faces()
            .iter()
            .map(|face| face.origin.clone())
            .collect();
        assert_eq!(field_names, brep_names, "{kind}");
        for origin in &field_names {
            assert!(field.select(origin).is_some(), "{kind}: {origin:?}");
        }
        let inside_query = field.proven(&point(inside), &mut Workspace::default());
        assert!(inside_query.distance.millimetres.hi().is_finite());
        let signed = inside_query.value.enclosure();
        let expected_distance = if kind == "cylinder" { -3.0 } else { -1.0 };
        assert!(signed.lo() <= expected_distance && expected_distance <= signed.hi());
        assert_eq!(
            field
                .proven(&point(inside), &mut Workspace::default())
                .location,
            Location::Inside,
            "{kind}"
        );
        assert_eq!(
            field
                .proven(&point([20, 20, 20]), &mut Workspace::default())
                .location,
            Location::Outside,
            "{kind}"
        );
        let boundary = point(wall);
        assert_eq!(
            field.proven(&boundary, &mut Workspace::default()).location,
            Location::Boundary,
            "{kind}"
        );
        assert!(
            field
                .fast(inside.map(|x| x as f64), &mut Workspace::default())
                .unwrap()
                .value
                < 0.0
        );
        let mut covered = BTreeSet::new();
        field
            .visit_candidates(2, |p, id| {
                if field
                    .face_query(id, &p, &mut Workspace::default())
                    .unwrap()
                    .state
                    == FaceState::OnFace
                {
                    covered.insert(field.faces()[id.index()].origin.clone());
                }
                true
            })
            .unwrap();
        assert_eq!(covered, field_names, "{kind}");
        assert!(field.faces().iter().enumerate().any(|(index, _)| {
            let query = field
                .face_query(
                    field.result_face(index).unwrap(),
                    &boundary,
                    &mut Workspace::default(),
                )
                .unwrap();
            query.state == FaceState::OnFace && matches!(query.direction, Direction::Unique(_))
        }));
        let beyond_trim = point([wall[0], wall[1], wall[2] + 50]);
        for index in 0..field.faces().len() {
            let query = field
                .face_query(
                    field.result_face(index).unwrap(),
                    &beyond_trim,
                    &mut Workspace::default(),
                )
                .unwrap();
            assert_eq!(query.state, FaceState::NotOnFace, "{kind}");
        }
        let volume =
            catallax_field::mass::volume(&field, &Default::default(), &mut Decider::default())
                .unwrap();
        let expected = if kind == "box" {
            Ball::from_rational(&volume_factor, 128)
        } else {
            Ball::pi(128) * Ball::from_rational(&volume_factor, 128)
        };
        assert!(
            volume.lo <= expected.hi() && expected.lo() <= volume.hi,
            "{kind}: {volume:?}"
        );
        let controller = Controller::default();
        controller.begin().check_with_files(
            &recipe,
            &files,
            &Options {
                divisions: 2,
                ..Default::default()
            },
        );
        let report = controller.latest().unwrap();
        assert_eq!(report.verdict, Verdict::Agree, "{kind}: {report:?}");
        let imported = report.import_agreement.unwrap();
        assert_eq!(imported.points.len(), 125);
        assert!(imported.points.iter().all(|p| p.brep == p.field));
        assert!(imported.points.iter().any(|p| p.field == Location::Inside));
        assert!(imported.points.iter().any(|p| p.field == Location::Outside));
        assert!(imported.volume_difference_upper <= imported.volume_difference_allowed);
    }
}
#[test]
fn imported_plan_queries_repeat_with_one_two_and_four_workers() {
    for (bytes, inside) in [
        (
            common::box_step(common::BoxFixture::new(["2", "3", "4"])).into_bytes(),
            [11, 21, 31],
        ),
        (exported("torus"), [14, 20, 30]),
    ] {
        repeat_with_workers(&bytes, inside);
    }
}

fn repeat_with_workers(bytes: &[u8], inside: [i64; 3]) {
    let mut pose = Pose::identity();
    pose.translation.0 = [10.into(), 20.into(), 30.into()];
    let field = Plan::compile_with_files(&recipe(bytes, pose), &sources(bytes)).unwrap();
    let points: Vec<_> = (8..=14)
        .flat_map(|x| (18..=25).flat_map(move |y| (28..=36).map(move |z| point([x, y, z]))))
        .collect();
    let expected: Vec<_> = points
        .iter()
        .map(|p| field.proven(p, &mut Workspace::default()))
        .collect();
    assert_eq!(
        field
            .proven(&point(inside), &mut Workspace::default())
            .location,
        Location::Inside
    );
    for workers in [1, 2, 4] {
        let actual = std::thread::scope(|scope| {
            let handles: Vec<_> = points
                .chunks(points.len().div_ceil(workers))
                .map(|chunk| {
                    let field = &field;
                    scope.spawn(move || {
                        let mut storage = Workspace::default();
                        chunk
                            .iter()
                            .map(|p| field.proven(p, &mut storage))
                            .collect::<Vec<_>>()
                    })
                })
                .collect();
            handles
                .into_iter()
                .flat_map(|handle| handle.join().unwrap())
                .collect::<Vec<_>>()
        });
        assert_eq!(actual, expected, "{workers} workers");
    }
}
#[test]
fn import_refusals_preserve_missing_bytes_hash_and_family_limits() {
    let bytes = common::box_step(common::BoxFixture::new(["2", "3", "4"])).into_bytes();
    let part = recipe(&bytes, Pose::identity());
    assert_eq!(
        Plan::compile(&part).unwrap_err().code,
        "field.import-missing-source"
    );
    let mut files = sources(&bytes);
    files.insert(catallax_import::sha256_hex(&bytes), b"changed");
    assert_eq!(
        Plan::compile_with_files(&part, &files).unwrap_err().code,
        "import.hash-mismatch"
    );
    let spline = include_bytes!("../../catallax-import/tests/fixtures/spline-roof.step");
    let refusal =
        Plan::compile_with_files(&recipe(spline, Pose::identity()), &sources(spline)).unwrap_err();
    assert_eq!(refusal.code, "field.import-spline");
    assert!(refusal.detail.contains("M10"));
}

fn names(field: &Plan, built: &catallax_features::Built) -> [BTreeSet<catallax_recipe::Origin>; 2] {
    [
        field
            .faces()
            .iter()
            .map(|face| face.origin.clone())
            .collect(),
        built
            .model
            .faces_of_body(built.result)
            .iter()
            .flat_map(|id| built.model.face_origins(*id))
            .cloned()
            .collect(),
    ]
}

fn raw_exported(build: impl FnOnce(&mut Model) -> catallax_topology::BodyId) -> Vec<u8> {
    let mut model = Model::new();
    let body = build(&mut model);
    let mut bytes = Vec::new();
    catallax_step::write_step(
        &model,
        &[body],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    bytes
}

#[test]
fn imported_cones_at_any_angle_build_through_their_rims_and_agree_on_request() {
    let identity = ShapePose::identity();
    // Radian semi-angles as the exporter writes them, and a degree file.
    let frustum = raw_exported(|m| cone(m, 3.into(), 1.into(), 2.into(), &identity).unwrap());
    let steep = raw_exported(|m| cone(m, 2.into(), 1.into(), 5.into(), &identity).unwrap());
    let pointed = raw_exported(|m| cone(m, 3.into(), 0.into(), 2.into(), &identity).unwrap());
    let mut degrees = common::Soup::new();
    degrees.cone(["0", "0", "0"], "1", "30.", Ok(("3", "2.732050808")));
    let degrees = common::in_degrees(&degrees.finish()).into_bytes();
    for (kind, bytes, inside, outside, faces) in [
        ("radian frustum", frustum, [0, 0, 1], [3, 0, 1], 3),
        ("steep radian frustum", steep, [0, 0, 4], [2, 0, 4], 3),
        ("pointed radian cone", pointed, [0, 0, 1], [0, 0, 3], 2),
        ("thirty-degree frustum", degrees, [0, 0, 2], [3, 0, 1], 3),
    ] {
        let recipe = recipe(&bytes, Pose::identity());
        let files = sources(&bytes);
        let field =
            Plan::compile_with_files(&recipe, &files).unwrap_or_else(|e| panic!("{kind}: {e:?}"));
        let built =
            catallax_features::evaluate_with_files(&recipe, &files, &mut Decider::default())
                .unwrap_or_else(|e| panic!("{kind}: {e}"));
        let [field_names, brep_names] = names(&field, &built);
        assert_eq!(field_names, brep_names, "{kind}");
        assert_eq!(field_names.len(), faces, "{kind}");
        let record = built
            .import_record(&StepId::new("imported").unwrap())
            .unwrap();
        assert_eq!(record.face_tolerances.len(), 1, "{kind}");
        assert_eq!(record.face_tolerances[0].approved_tolerance_mm, "0.01");
        for (p, expected) in [(inside, Location::Inside), (outside, Location::Outside)] {
            assert_eq!(
                field.proven(&point(p), &mut Workspace::default()).location,
                expected,
                "{kind} at {p:?}"
            );
        }
        let controller = Controller::default();
        controller.begin().check_with_files(
            &recipe,
            &files,
            &Options {
                divisions: 2,
                ..Default::default()
            },
        );
        let report = controller.latest().unwrap();
        assert_eq!(
            report.verdict,
            Verdict::Agree,
            "{kind}: {:?} {:?}",
            report.refusal,
            report.findings.first()
        );
        let imported = report.import_agreement.unwrap();
        assert!(imported.points.iter().all(|p| p.brep == p.field), "{kind}");
        assert!(imported.volume_agrees, "{kind}");
    }
}

#[test]
fn imported_spheres_and_tori_keep_names_and_agree_on_request() {
    let mut vertex_loop = common::Soup::new();
    vertex_loop.vertex_loop_sphere(["1", "2", "3"], "2");
    // A quarter turn about x carries the torus axis onto -y.
    let turned = Pose {
        rotation: [1.into(), 1.into(), 0.into(), 0.into()],
        translation: catallax_recipe::Vector3([10.into(), 0.into(), 0.into()]),
    };
    let pi = Ball::pi(128);
    for (kind, bytes, pose, inside, boundary, outside, volume) in [
        (
            "sphere",
            exported("sphere"),
            Pose::identity(),
            [0, 0, 0],
            [0, 3, 0],
            [0, 0, 4],
            &pi * &Ball::from_rational(&q(36), 128),
        ),
        (
            "vertex-loop sphere",
            vertex_loop.finish().into_bytes(),
            Pose::identity(),
            [1, 2, 3],
            [1, 2, 5],
            [3, 4, 3],
            &pi * &Ball::from_rational(&Rational::parse("32/3").unwrap(), 128),
        ),
        (
            "torus",
            exported("torus"),
            turned,
            [14, 0, 0],
            [14, 1, 0],
            [10, 0, 0],
            &(&pi * &pi) * &Ball::from_rational(&q(8), 128),
        ),
    ] {
        let recipe = recipe(&bytes, pose);
        let files = sources(&bytes);
        let field =
            Plan::compile_with_files(&recipe, &files).unwrap_or_else(|e| panic!("{kind}: {e:?}"));
        let built =
            catallax_features::evaluate_with_files(&recipe, &files, &mut Decider::default())
                .unwrap();
        let [field_names, brep_names] = names(&field, &built);
        assert_eq!(field_names, brep_names, "{kind}");
        assert_eq!(field_names.len(), 1, "{kind}");
        for origin in &field_names {
            assert!(field.select(origin).is_some(), "{kind}: {origin:?}");
        }
        for (p, expected) in [
            (inside, Location::Inside),
            (boundary, Location::Boundary),
            (outside, Location::Outside),
        ] {
            assert_eq!(
                field.proven(&point(p), &mut Workspace::default()).location,
                expected,
                "{kind} at {p:?}"
            );
        }
        let id = field.result_face(0).unwrap();
        let query = field
            .face_query(id, &point(boundary), &mut Workspace::default())
            .unwrap();
        assert_eq!(query.state, FaceState::OnFace, "{kind}");
        let measured =
            catallax_field::mass::volume(&field, &Default::default(), &mut Decider::default())
                .unwrap();
        assert!(
            measured.lo <= volume.hi() && volume.lo() <= measured.hi,
            "{kind}: {measured:?}"
        );
        let controller = Controller::default();
        controller.begin().check_with_files(
            &recipe,
            &files,
            &Options {
                divisions: 2,
                ..Default::default()
            },
        );
        let report = controller.latest().unwrap();
        assert_eq!(
            report.verdict,
            Verdict::Agree,
            "{kind}: {:?} {:?}",
            report.refusal,
            report.findings.first()
        );
        let imported = report.import_agreement.unwrap();
        assert_eq!(imported.points.len(), 125);
        assert!(imported.points.iter().all(|p| p.brep == p.field));
        assert!(imported.points.iter().any(|p| p.field == Location::Inside));
        assert!(imported.points.iter().any(|p| p.field == Location::Outside));
    }
}

#[test]
fn nonconvex_planar_import_keeps_the_missing_corner_outside() {
    let polygon = [
        ["0", "0"],
        ["2", "0"],
        ["2", "1"],
        ["1", "1"],
        ["1", "2"],
        ["0", "2"],
    ];
    let mut bottom: common::Ring = polygon.iter().map(|p| [p[0], p[1], "0"]).collect();
    bottom.reverse();
    let top: common::Ring = polygon.iter().map(|p| [p[0], p[1], "2"]).collect();
    let mut faces = vec![bottom, top];
    for i in 0..polygon.len() {
        let a = polygon[i];
        let b = polygon[(i + 1) % polygon.len()];
        faces.push(vec![
            [a[0], a[1], "0"],
            [b[0], b[1], "0"],
            [b[0], b[1], "2"],
            [a[0], a[1], "2"],
        ]);
    }
    let bytes = common::polyhedra_step(&[faces]).into_bytes();
    let recipe = recipe(&bytes, Pose::identity());
    let field = Plan::compile_with_files(&recipe, &sources(&bytes)).unwrap();
    let rational_point = |x: &str, y: &str| {
        RPoint3::new(
            Rational::parse(x).unwrap(),
            Rational::parse(y).unwrap(),
            q(1),
        )
    };
    assert_eq!(
        field
            .proven(&rational_point("1/2", "3/2"), &mut Workspace::default())
            .location,
        Location::Inside
    );
    assert_eq!(
        field
            .proven(&rational_point("3/2", "3/2"), &mut Workspace::default())
            .location,
        Location::Outside
    );
    let volume =
        catallax_field::mass::volume(&field, &Default::default(), &mut Decider::default()).unwrap();
    assert_eq!(volume.lo, q(6));
    assert_eq!(volume.hi, q(6));
}

#[test]
fn inward_planar_shell_is_refused_instead_of_published_as_empty() {
    let text = common::box_step(common::BoxFixture::new(["2", "3", "4"]));
    let mut inward = String::new();
    for line in text.lines() {
        if line.contains("ADVANCED_FACE") {
            if let Some(prefix) = line.strip_suffix(".T.);") {
                inward.push_str(&format!("{prefix}.F.);"));
            } else {
                inward.push_str(&line.replace(".F.);", ".T.);"));
            }
        } else {
            inward.push_str(line);
        }
        inward.push('\n');
    }
    let bytes = inward.as_bytes();
    let refusal =
        Plan::compile_with_files(&recipe(bytes, Pose::identity()), &sources(bytes)).unwrap_err();
    assert_eq!(refusal.code, "field.import-domain");
    assert!(refusal.detail.contains("0 of 6"));
}

#[test]
fn imported_coplanar_face_names_keep_their_original_finite_trims() {
    for elbow in [false, true] {
        let mut left = common::box_rings(["0", "0", "0"], ["1", "1", "2"]);
        if elbow {
            left.remove(3);
        }
        left.remove(1);
        let mut right = common::box_rings(["1", "0", "0"], ["2", "1", "2"]);
        right.remove(0);
        left.extend(right);
        if elbow {
            let mut upper = common::box_rings(["0", "1", "0"], ["1", "2", "2"]);
            upper.remove(2);
            left.extend(upper);
        }
        let bytes = common::polyhedra_step(&[left]).into_bytes();
        let recipe = recipe(&bytes, Pose::identity());
        let files = sources(&bytes);
        let field = Plan::compile_with_files(&recipe, &files).unwrap();
        let mut built =
            catallax_features::evaluate_with_files(&recipe, &files, &mut Decider::default())
                .unwrap();
        let pose = ShapePose::from_rational(
            catallax_shapes::Rotation::from_quaternion(q(1), q(1), q(2), q(3)).unwrap(),
            catallax_shapes::RVector3::new(q(7), q(-3), q(5)),
        );
        let moved = field.transformed(&pose, &mut Decider::default()).unwrap();
        for x in ["1/2", "1", "3/2"] {
            let point = RPoint3::new(
                Rational::parse(x).unwrap(),
                Rational::parse("1/2").unwrap(),
                q(2),
            );
            let transformed = pose.apply_rational_point(&point).unwrap();
            for face in built.model.faces_of_body(built.result) {
                let expected = catallax_boolean::face_contains(
                    &mut built.model,
                    face,
                    &point,
                    &mut Decider::default(),
                )
                .unwrap();
                for origin in built.model.face_origins(face) {
                    for (plan, p) in [(&field, &point), (&moved, &transformed)] {
                        let selected = plan.select(origin).unwrap();
                        let actual = plan
                            .face_query(selected, p, &mut Workspace::default())
                            .unwrap()
                            .state
                            == FaceState::OnFace;
                        assert_eq!(actual, expected, "{origin:?} at x={x}, elbow={elbow}");
                    }
                }
            }
        }
    }
}

#[test]
fn previews_declared_cones_at_other_angles_in_degrees_and_radians() {
    use catallax_import::Surface;
    use catallax_mesh::{PreviewOptions, preview_step};
    let mut model = Model::new();
    let body = cone(
        &mut model,
        2.into(),
        1.into(),
        5.into(),
        &ShapePose::identity(),
    )
    .unwrap();
    let mut bytes = Vec::new();
    catallax_step::write_step(
        &model,
        &[body],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    let radians = String::from_utf8(bytes).unwrap();
    // The same cone with its semi-angle, atan(1/5), restated in degrees.
    let line = radians
        .lines()
        .find(|l| l.contains("CONICAL_SURFACE"))
        .unwrap();
    let (head, _) = line.rsplit_once(',').unwrap();
    let degrees =
        common::in_degrees(&radians.replace(line, &format!("{head},11.309932474020213);")));
    for text in [radians, degrees] {
        let preview = preview_step(
            text.as_bytes(),
            &PreviewOptions::default(),
            &mut Decider::default(),
        )
        .unwrap_or_else(|e| panic!("{e:?}"));
        let wall = preview.source.solids[0]
            .faces
            .iter()
            .find(|f| matches!(f.surface, Surface::Cone(..)))
            .unwrap();
        assert!(preview.mesh.faces.iter().any(|f| f.face == wall.name));
        for p in &preview.mesh.positions {
            let radius = (p[0] * p[0] + p[1] * p[1]).sqrt();
            assert!(
                radius <= 2.0 + 1e-6 && p[2] >= -1e-9 && p[2] <= 5.0 + 1e-9,
                "{p:?}"
            );
        }
    }
}
