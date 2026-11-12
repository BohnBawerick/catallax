//! Imported parts with spline faces: the flat bilinear patches the shared
//! importer admits build in the field under the file's face names, and every
//! other spline face stays a typed refusal.

#[path = "../../catallax-import/tests/common/mod.rs"]
mod common;

use catallax_compare::{Controller, Options, Verdict};
use catallax_features::ImportFiles;
use catallax_field::{FaceState, Location, Plan, Workspace};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{FileRef, Operation, Pose, Recipe, Step, StepId};
use catallax_shapes::{Pose as ShapePose, RPoint3, RVector3, Rotation};
use std::collections::BTreeSet;

fn q(text: &str) -> Rational {
    Rational::parse(text).unwrap()
}
fn point(v: [&str; 3]) -> RPoint3 {
    let [x, y, z] = v.map(q);
    RPoint3::new(x, y, z)
}
fn recipe(bytes: &[u8], pose: Pose) -> Recipe {
    let id = StepId::new("imported").unwrap();
    Recipe::new(
        "imported spline",
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

/// How the roof patch is written: its four corner heights over the box base
/// (`[x0y0, x1y0, x1y1, x0y1]`), whether its parameters are swapped, whether the
/// face carries the opposite sense to match, and its weights.
struct Roof {
    heights: [&'static str; 4],
    swapped: bool,
    weights: Option<[&'static str; 4]>,
}

/// A 2 by 3 base, a floor at 0 and a roof written as a bilinear spline patch.
/// The side faces are planar, so a roof that is not planar is not admitted.
fn roof_part(roof: &Roof) -> Vec<u8> {
    let [h00, h10, h11, h01] = roof.heights;
    let mut rings = common::box_rings(["0", "0", "0"], ["2", "3", "4"]);
    for ring in &mut rings {
        for corner in ring.iter_mut() {
            if corner[2] == "4" {
                corner[2] = match (corner[0], corner[1]) {
                    ("0", "0") => h00,
                    ("2", "0") => h10,
                    ("2", "3") => h11,
                    _ => h01,
                };
            }
        }
    }
    let text = common::polyhedra_step(&[rings]);
    let top = text
        .lines()
        .filter(|line| line.contains("ADVANCED_FACE"))
        .nth(5)
        .unwrap();
    let surface_id = top.split(',').nth(2).unwrap().trim_start_matches('#');
    let plane = text
        .lines()
        .find(|line| line.starts_with(&format!("#{surface_id}=PLANE(")))
        .unwrap();
    let corners = [
        ["0", "0", h00],
        ["2", "0", h10],
        ["2", "3", h11],
        ["0", "3", h01],
    ];
    // Controls run u-major: u along x and v along y, or the other way round.
    let order: [usize; 4] = if roof.swapped {
        [0, 1, 3, 2]
    } else {
        [0, 3, 1, 2]
    };
    let mut points = String::new();
    for (slot, corner) in order.iter().enumerate() {
        let c = corners[*corner];
        points.push_str(&format!(
            "#{}=CARTESIAN_POINT('',({},{},{}));\n",
            900001 + slot,
            common::dec(&common::r(c[0])),
            common::dec(&common::r(c[1])),
            common::dec(&common::r(c[2]))
        ));
    }
    let surface = match roof.weights {
        None => format!(
            "#{surface_id}=B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((#900001,#900002),(#900003,#900004)),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.);"
        ),
        Some(w) => format!(
            "#{surface_id}=(BOUNDED_SURFACE() B_SPLINE_SURFACE(1,1,((#900001,#900002),(#900003,#900004)),.UNSPECIFIED.,.F.,.F.,.F.) B_SPLINE_SURFACE_WITH_KNOTS((2,2),(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.) GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE((({},{}),({},{}))) REPRESENTATION_ITEM('') SURFACE());",
            w[0], w[1], w[2], w[3]
        ),
    };
    let mut text = text.replace(plane, &surface);
    if roof.swapped {
        // Swapping the parameters reverses the patch normal; the face says so.
        let flipped = top.replace(",.T.);", ",.F.);");
        text = text.replace(top, &flipped);
    }
    text.replace("DATA;", &format!("DATA;\n{points}"))
        .into_bytes()
}

/// A box whose first edge is written as a straight B-spline curve.
fn spline_edge_box() -> Vec<u8> {
    let text = common::box_step(common::BoxFixture::new(["2", "3", "4"]));
    let file = catallax_import::read_step(text.as_bytes()).unwrap();
    let edge = &file.solids[0].faces[0].loops[0].edges[0];
    let edge_line = text
        .lines()
        .find(|line| line.starts_with(&format!("{}=", edge.edge)))
        .unwrap();
    let curve_id = edge_line.split(',').nth(3).unwrap();
    let line = text
        .lines()
        .find(|line| line.starts_with(&format!("{curve_id}=")))
        .unwrap();
    let (a, b) = (&edge.start.point, &edge.end.point);
    let mid: [Rational; 3] = std::array::from_fn(|i| (&a[i] + &b[i]).mul_pow2(-1));
    let mut points = String::new();
    for (slot, p) in [a, &mid, b].into_iter().enumerate() {
        points.push_str(&format!(
            "#{}=CARTESIAN_POINT('',({},{},{}));\n",
            900001 + slot,
            common::dec(&p[0]),
            common::dec(&p[1]),
            common::dec(&p[2])
        ));
    }
    text.replace(line, &format!("{curve_id}=B_SPLINE_CURVE_WITH_KNOTS('',2,(#900001,#900002,#900003),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,1.),.UNSPECIFIED.);"))
        .replace("DATA;", &format!("DATA;\n{points}"))
        .into_bytes()
}

struct Case {
    kind: &'static str,
    bytes: Vec<u8>,
    inside: [&'static str; 3],
    boundary: [&'static str; 3],
    outside: [&'static str; 3],
    volume: i64,
}
fn cases() -> Vec<Case> {
    let level = Roof {
        heights: ["4", "4", "4", "4"],
        swapped: false,
        weights: None,
    };
    let sloped = |swapped, weights| Roof {
        heights: ["4", "4", "5", "5"],
        swapped,
        weights,
    };
    vec![
        Case {
            kind: "straight spline edge",
            bytes: spline_edge_box(),
            inside: ["1", "1", "1"],
            boundary: ["1", "1", "4"],
            outside: ["1", "1", "5"],
            volume: 24,
        },
        Case {
            kind: "level flat patch",
            bytes: roof_part(&level),
            inside: ["1", "1", "1"],
            boundary: ["1", "1", "4"],
            outside: ["1", "1", "5"],
            volume: 24,
        },
        Case {
            kind: "sloped flat patch",
            bytes: roof_part(&sloped(false, None)),
            inside: ["1", "1", "1"],
            boundary: ["1", "1", "13/3"],
            outside: ["1", "1", "6"],
            volume: 27,
        },
        Case {
            kind: "sloped flat patch, parameters swapped",
            bytes: roof_part(&sloped(true, None)),
            inside: ["1", "1", "1"],
            boundary: ["1", "1", "13/3"],
            outside: ["1", "1", "6"],
            volume: 27,
        },
        Case {
            kind: "sloped flat patch with equal weights",
            bytes: roof_part(&sloped(false, Some(["2", "2", "2", "2"]))),
            inside: ["1", "1", "1"],
            boundary: ["1", "1", "13/3"],
            outside: ["1", "1", "6"],
            volume: 27,
        },
    ]
}

fn names(field: &Plan, built: &catallax_features::Built) -> [BTreeSet<catallax_recipe::Origin>; 2] {
    [
        field.faces().iter().map(|f| f.origin.clone()).collect(),
        built
            .model
            .faces_of_body(built.result)
            .iter()
            .flat_map(|id| built.model.face_origins(*id))
            .cloned()
            .collect(),
    ]
}

fn turned() -> Pose {
    Pose {
        rotation: [1.into(), 1.into(), 0.into(), 0.into()],
        translation: catallax_recipe::Vector3([10.into(), 0.into(), 0.into()]),
    }
}
fn placement(pose: &Pose) -> ShapePose {
    let [w, x, y, z] = pose.rotation.each_ref().map(|n| n.value().clone());
    let [tx, ty, tz] = pose.translation.0.each_ref().map(|n| n.value().clone());
    ShapePose::from_rational(
        Rotation::from_quaternion(w, x, y, z).unwrap(),
        RVector3::new(tx, ty, tz),
    )
}

fn check(case: &Case, pose: Pose) {
    let kind = case.kind;
    let recipe = recipe(&case.bytes, pose.clone());
    let files = sources(&case.bytes);
    let field =
        Plan::compile_with_files(&recipe, &files).unwrap_or_else(|e| panic!("{kind}: {e:?}"));
    let built =
        catallax_features::evaluate_with_files(&recipe, &files, &mut Decider::default()).unwrap();
    let [named, brep] = names(&field, &built);
    assert_eq!(named, brep, "{kind}");
    assert_eq!(named.len(), 6, "{kind}");
    for origin in &named {
        assert!(field.select(origin).is_some(), "{kind}: {origin:?}");
    }
    let place = placement(&pose);
    let at = |v: [&str; 3]| place.apply_rational_point(&point(v)).unwrap();
    let (inside, boundary, outside) = (at(case.inside), at(case.boundary), at(case.outside));
    for (p, expected) in [
        (&inside, Location::Inside),
        (&boundary, Location::Boundary),
        (&outside, Location::Outside),
    ] {
        assert_eq!(
            field.proven(p, &mut Workspace::default()).location,
            expected,
            "{kind} at {p:?}"
        );
    }
    let holding = |p: &RPoint3| {
        (0..field.faces().len())
            .filter(|index| {
                field
                    .face_query(
                        field.result_face(*index).unwrap(),
                        p,
                        &mut Workspace::default(),
                    )
                    .unwrap()
                    .state
                    == FaceState::OnFace
            })
            .count()
    };
    assert_eq!(holding(&boundary), 1, "{kind}: one face holds the point");
    assert_eq!(holding(&inside), 0, "{kind}");
    let measured =
        catallax_field::mass::volume(&field, &Default::default(), &mut Decider::default()).unwrap();
    let volume = Rational::from_i64(case.volume);
    assert!(
        measured.lo <= volume && volume <= measured.hi,
        "{kind}: {measured:?}"
    );

    // On request, the two kernels agree: membership and volume, then every
    // sampled surface point and direction.
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

#[test]
fn flat_spline_patches_keep_names_and_agree_on_request() {
    for case in cases() {
        check(&case, Pose::identity());
    }
}

#[test]
fn placed_flat_spline_patches_keep_names_and_agree_on_request() {
    let case = cases().remove(2);
    check(&case, turned());
}

#[test]
fn spline_part_answers_repeat_with_one_two_four_and_all_workers() {
    let all = std::thread::available_parallelism().map_or(4, |n| n.get());
    for case in cases() {
        let mut pose = Pose::identity();
        pose.translation.0 = [10.into(), 20.into(), 30.into()];
        let field =
            Plan::compile_with_files(&recipe(&case.bytes, pose), &sources(&case.bytes)).unwrap();
        let points: Vec<_> = (-2..=2)
            .flat_map(|i| {
                (-2..=2).flat_map(move |j| {
                    (-2..=2).map(move |k| {
                        RPoint3::new(
                            Rational::from_ratio(10 * 4 + i * 5, 4).unwrap(),
                            Rational::from_ratio(20 * 4 + j * 5, 4).unwrap(),
                            Rational::from_ratio(30 * 4 + k * 5, 4).unwrap(),
                        )
                    })
                })
            })
            .collect();
        let expected: Vec<_> = points
            .iter()
            .map(|p| field.proven(p, &mut Workspace::default()))
            .collect();
        for workers in [1, 2, 4, all] {
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
            assert_eq!(actual, expected, "{}: {workers} workers", case.kind);
        }
    }
}

fn refused(bytes: &[u8]) -> catallax_numbers::Refusal {
    Plan::compile_with_files(&recipe(bytes, Pose::identity()), &sources(bytes)).unwrap_err()
}

#[test]
fn spline_faces_the_field_cannot_hold_stay_typed_refusals() {
    // A saddle roof: the shared importer admits it, the field does not yet.
    let saddle = roof_part(&Roof {
        heights: ["4", "4", "9/2", "4"],
        swapped: false,
        weights: None,
    });
    let recipe = recipe(&saddle, Pose::identity());
    catallax_features::evaluate_with_files(&recipe, &sources(&saddle), &mut Decider::default())
        .expect("the B-rep builds the saddle roof");
    let refusal = refused(&saddle);
    assert_eq!(refusal.code, "field.import-spline");
    assert!(
        refusal.detail.contains("curved spline face"),
        "{}",
        refusal.detail
    );
    // Unequal weights make a rational patch, which the field does not read.
    let weighted = roof_part(&Roof {
        heights: ["4", "4", "5", "5"],
        swapped: false,
        weights: Some(["1", "2", "1", "1"]),
    });
    let refusal = refused(&weighted);
    assert_eq!(refusal.code, "field.import-spline");
    // A spline file the shared importer does not certify keeps its M10 reason.
    let roof = include_bytes!("../../catallax-import/tests/fixtures/spline-roof.step");
    let refusal = refused(roof);
    assert_eq!(refusal.code, "field.import-spline");
    assert!(refusal.detail.contains("M10"), "{}", refusal.detail);
}
