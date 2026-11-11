//! Imported bodies of revolution with trimmed spheres and tori build in the field,
//! keep the file's face names and agree with the B-rep on request.

#[path = "../../catallax-import/tests/common/mod.rs"]
mod common;

use catallax_compare::{Controller, Options, Verdict};
use catallax_features::{ImportFiles, Profile, revolve};
use catallax_field::{FaceState, Location, Plan, Workspace};
use catallax_numbers::{Ball, Decider, Rational};
use catallax_recipe::{
    Axis, ElementId, FileRef, Loop, Operation, Point2, Point3, Pose, Recipe, Segment, Step, StepId,
    Turn, Vector3,
};
use catallax_shapes::{Pose as ShapePose, RPoint3, RVector3, Rotation};
use catallax_topology::Model;
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
        "imported revolution",
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

/// A quarter turn about x, then a shift along x.
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

fn pi_times(factor: &str) -> Ball {
    &Ball::pi(128) * &Ball::from_rational(&q(factor), 128)
}
fn pi_squared_times(factor: &str) -> Ball {
    &(&Ball::pi(128) * &Ball::pi(128)) * &Ball::from_rational(&q(factor), 128)
}

/// A cylinder of radius 4 and height 5 about y whose top edge is rounded by a
/// ball of radius 1, written by the exporter: a plane, a cylinder, a torus, a plane.
fn rounded_cylinder() -> Vec<u8> {
    let line = |id: &str, x: i64, y: i64| Segment::Line {
        id: ElementId::new(id).unwrap(),
        to: Point2([x.into(), y.into()]),
    };
    let loops = vec![Loop(vec![
        line("floor", 4, 0),
        line("wall", 4, 4),
        Segment::Arc {
            id: ElementId::new("round").unwrap(),
            centre: Point2([3.into(), 4.into()]),
            to: Point2([3.into(), 5.into()]),
            turn: Turn::Ccw,
        },
        line("roof", 0, 5),
        line("axis", 0, 0),
    ])];
    let mut dec = Decider::default();
    let profile = Profile::new(&Pose::identity(), &loops, &mut dec).unwrap();
    let mut model = Model::new();
    let axis = Axis {
        point: Point3([0.into(), 0.into(), 0.into()]),
        direction: Vector3([0.into(), 1.into(), 0.into()]),
    };
    let body = revolve(
        &mut model,
        &profile,
        &StepId::new("turn").unwrap(),
        &axis,
        &360.into(),
        &mut dec,
    )
    .unwrap();
    let mut bytes = Vec::new();
    catallax_step::write_step(&model, &[body], &Default::default(), &mut dec, &mut bytes).unwrap();
    bytes
}

struct Case {
    kind: &'static str,
    bytes: Vec<u8>,
    inside: [&'static str; 3],
    boundary: [&'static str; 3],
    outside: [&'static str; 3],
    volume: Volume,
    faces: usize,
}
enum Volume {
    Exact(Ball),
    /// A closed form with an inverse cosine: its value, compared to 1e-9.
    Approximate(f64),
}

fn cases() -> Vec<Case> {
    let mut hemisphere = common::Soup::new();
    hemisphere.hemisphere(["0", "0", "0"], 2, "3");
    let mut half_torus = common::Soup::new();
    half_torus.half_torus(["0", "0", "0"], 2, "4", "1");
    let mut domed = common::Soup::new();
    domed.domed_cylinder(["0", "0", "0"], 2, "3", "2");
    let mut lens = common::Soup::new();
    lens.lens(["0", "0", "0"], 2, "4", "3", "5");
    let mut zone = common::Soup::new();
    zone.sphere_zone(["1", "2", "3"], 0, "5", ("-3", "4"), ("4", "3"));
    let mut cap = common::Soup::new();
    cap.torus_cap(["0", "0", "0"], 2, "4", "1", "3/5", ("16/5", "24/5"));
    let mut tilted = common::Soup::new();
    tilted.hemisphere_about(["0", "0", "0"], ["0", "3/5", "4/5"], ["1", "0", "0"], "5");
    let (acos, pi) = ((3.0_f64 / 5.0).acos(), std::f64::consts::PI);
    vec![
        Case {
            kind: "hemisphere",
            bytes: hemisphere.finish().into_bytes(),
            inside: ["0", "0", "1"],
            boundary: ["0", "0", "3"],
            outside: ["0", "0", "-1"],
            volume: Volume::Exact(pi_times("18")),
            faces: 2,
        },
        Case {
            kind: "half torus",
            bytes: half_torus.finish().into_bytes(),
            inside: ["4", "0", "1/2"],
            boundary: ["4", "0", "1"],
            outside: ["4", "0", "-1/2"],
            volume: Volume::Exact(pi_squared_times("4")),
            faces: 2,
        },
        Case {
            kind: "domed cylinder",
            bytes: domed.finish().into_bytes(),
            inside: ["0", "0", "1"],
            boundary: ["2", "0", "1"],
            outside: ["0", "0", "6"],
            volume: Volume::Exact(pi_times("52/3")),
            faces: 3,
        },
        Case {
            kind: "lens",
            bytes: lens.finish().into_bytes(),
            inside: ["0", "0", "0"],
            boundary: ["0", "0", "1"],
            outside: ["0", "0", "2"],
            volume: Volume::Exact(pi_times("28/3")),
            faces: 2,
        },
        Case {
            kind: "sphere zone about x",
            bytes: zone.finish().into_bytes(),
            inside: ["1", "2", "3"],
            boundary: ["5", "2", "3"],
            outside: ["7", "2", "3"],
            volume: Volume::Exact(pi_times("434/3")),
            faces: 3,
        },
        Case {
            kind: "torus cap",
            bytes: cap.finish().into_bytes(),
            inside: ["4", "0", "4/5"],
            boundary: ["4", "0", "1"],
            outside: ["4", "0", "1/2"],
            volume: Volume::Approximate(8.0 * pi * (acos - 12.0 / 25.0)),
            faces: 2,
        },
        Case {
            kind: "hemisphere about a tilted axis",
            bytes: tilted.finish().into_bytes(),
            inside: ["0", "1", "1"],
            boundary: ["0", "3", "4"],
            outside: ["0", "-1", "-1"],
            volume: Volume::Exact(pi_times("250/3")),
            faces: 2,
        },
        Case {
            kind: "rounded cylinder from the exporter",
            bytes: rounded_cylinder(),
            inside: ["0", "2", "0"],
            boundary: ["4", "2", "0"],
            outside: ["0", "6", "0"],
            volume: Volume::Exact(&pi_times("221/3") + &pi_squared_times("3/2")),
            faces: 4,
        },
    ]
}

fn check(case: &Case, pose: Pose) {
    let Case { kind, bytes, .. } = case;
    let recipe = recipe(bytes, pose.clone());
    let files = sources(bytes);
    let field =
        Plan::compile_with_files(&recipe, &files).unwrap_or_else(|e| panic!("{kind}: {e:?}"));
    let built =
        catallax_features::evaluate_with_files(&recipe, &files, &mut Decider::default()).unwrap();

    // The file's face names survive, one for one.
    let brep: BTreeSet<_> = built
        .model
        .faces_of_body(built.result)
        .iter()
        .flat_map(|id| built.model.face_origins(*id))
        .cloned()
        .collect();
    let named: BTreeSet<_> = field.faces().iter().map(|f| f.origin.clone()).collect();
    assert_eq!(named, brep, "{kind}");
    assert_eq!(named.len(), case.faces, "{kind}");
    for origin in &named {
        assert!(field.select(origin).is_some(), "{kind}: {origin:?}");
    }

    // Exact membership, with the points carried by the same placement.
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
    let on_face = |p: &RPoint3| {
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
    assert!(
        on_face(&boundary) > 0,
        "{kind}: no face holds the boundary point"
    );
    assert_eq!(on_face(&inside), 0, "{kind}");
    assert_eq!(on_face(&outside), 0, "{kind}");

    // The field's volume encloses the closed form.
    let measured =
        catallax_field::mass::volume(&field, &Default::default(), &mut Decider::default()).unwrap();
    match &case.volume {
        Volume::Exact(volume) => assert!(
            measured.lo <= volume.hi() && volume.lo() <= measured.hi,
            "{kind}: {measured:?}"
        ),
        Volume::Approximate(value) => {
            let (lo, hi) = (measured.lo.to_f64_nearest(), measured.hi.to_f64_nearest());
            assert!(
                lo - 1e-9 <= *value && *value <= hi + 1e-9,
                "{kind}: {lo} {hi} {value}"
            );
        }
    }

    // On request, the two kernels agree on membership and volume.
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
    assert!(imported.points.iter().any(|p| p.field == Location::Inside));
    assert!(imported.points.iter().any(|p| p.field == Location::Outside));
    assert!(imported.volume_agrees, "{kind}");
}

fn named(kind: &str) -> Case {
    cases()
        .into_iter()
        .find(|case| case.kind == kind)
        .unwrap_or_else(|| panic!("no case {kind}"))
}

macro_rules! agree_on_request {
    ($($test:ident: $kind:expr, $pose:expr;)*) => {$(
        #[test]
        fn $test() {
            check(&named($kind), $pose);
        }
    )*};
}
agree_on_request! {
    hemisphere_keeps_names_and_agrees: "hemisphere", Pose::identity();
    half_torus_keeps_names_and_agrees: "half torus", Pose::identity();
    domed_cylinder_keeps_names_and_agrees: "domed cylinder", Pose::identity();
    lens_keeps_names_and_agrees: "lens", Pose::identity();
    sphere_zone_about_x_keeps_names_and_agrees: "sphere zone about x", Pose::identity();
    torus_cap_keeps_names_and_agrees: "torus cap", Pose::identity();
    tilted_hemisphere_keeps_names_and_agrees: "hemisphere about a tilted axis", Pose::identity();
    exported_rounded_cylinder_keeps_names_and_agrees: "rounded cylinder from the exporter", Pose::identity();
    // A quarter turn and a shift carry each axis off its original direction.
    placed_domed_cylinder_keeps_names_and_agrees: "domed cylinder", turned();
    placed_tilted_hemisphere_keeps_names_and_agrees: "hemisphere about a tilted axis", turned();
}

#[test]
fn imported_revolution_answers_repeat_with_one_two_four_and_all_workers() {
    let all = std::thread::available_parallelism().map_or(4, |n| n.get());
    for case in cases() {
        let mut pose = Pose::identity();
        pose.translation.0 = [10.into(), 20.into(), 30.into()];
        let recipe = recipe(&case.bytes, pose);
        let field = Plan::compile_with_files(&recipe, &sources(&case.bytes)).unwrap();
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
    let recipe = recipe(bytes, Pose::identity());
    let files = sources(bytes);
    // The B-rep builds each of these; only the field declines.
    catallax_features::evaluate_with_files(&recipe, &files, &mut Decider::default())
        .unwrap_or_else(|e| panic!("{e}"));
    Plan::compile_with_files(&recipe, &files).unwrap_err()
}

#[test]
fn a_sphere_cut_about_two_axes_is_not_a_body_of_revolution_and_stays_refused() {
    let mut soup = common::Soup::new();
    soup.sphere_with_flats(
        ["0", "0", "0"],
        "5",
        [
            (["0", "0", "1"], "4", "3", ["1", "0", "0"]),
            (["1", "0", "0"], "4", "3", ["0", "1", "0"]),
        ],
    );
    let refusal = refused(soup.finish().as_bytes());
    assert_eq!(refusal.code, "field.import-domain");
    assert!(
        refusal.detail.contains("one axis direction"),
        "{}",
        refusal.detail
    );
}

/// Flips every face's same-sense flag, or only the flat faces'.
fn reversed(text: &str, flats_only: bool) -> String {
    let mut out = String::new();
    for line in text.lines() {
        let flat = line.contains("ADVANCED_FACE") && {
            let surface = line.rsplit_once(",#").map_or("", |(_, rest)| rest);
            let id: String = surface.chars().take_while(char::is_ascii_digit).collect();
            text.lines()
                .any(|l| l.starts_with(&format!("#{id}=PLANE(")))
        };
        if line.contains("ADVANCED_FACE") && (!flats_only || flat) {
            if let Some(prefix) = line.strip_suffix(".T.);") {
                out.push_str(&format!("{prefix}.F.);"));
            } else {
                out.push_str(&line.replace(".F.);", ".T.);"));
            }
        } else {
            out.push_str(line);
        }
        out.push('\n');
    }
    out
}

#[test]
fn inward_and_mixed_orientations_are_refused_instead_of_published_as_solids() {
    let mut soup = common::Soup::new();
    soup.domed_cylinder(["0", "0", "0"], 2, "3", "2");
    let text = soup.finish();
    for (flats_only, what) in [(false, "every face"), (true, "the flat face only")] {
        let bytes = reversed(&text, flats_only).into_bytes();
        let recipe = recipe(&bytes, Pose::identity());
        let refusal = Plan::compile_with_files(&recipe, &sources(&bytes)).unwrap_err();
        assert_eq!(refusal.code, "field.import-domain", "{what}");
        assert!(
            refusal
                .detail
                .contains("faces with spherical or toroidal trims"),
            "{what}: {}",
            refusal.detail
        );
    }
}
