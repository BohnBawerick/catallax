//! A proper STEP trim over one of the independently exported M27 loft faces.

use catallax_boolean::{Containment, classify_point};
use catallax_features::{ImportFiles, evaluate_with_files};
use catallax_field::{Location, Plan, Value, Workspace};
use catallax_numbers::{Decider, Outcome, Rational, with_operation_budget};
use catallax_recipe::{FileRef, Operation, Pose, Recipe, Step, StepId};
use catallax_shapes::{RPoint3, SplineSupport, Surface};

mod import_support;
use import_support::trimmed_loft;
const WORK: u64 = 20_000_000;
// The STEP header's UNCERTAINTY_MEASURE_WITH_UNIT in millimetres.
const TOLERANCE_MM: f64 = 1e-7;

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn point(at: [&str; 3]) -> RPoint3 {
    let [x, y, z] = at.map(q);
    RPoint3::new(x, y, z)
}

#[test]
fn a_trimmed_occt_step_solid_has_the_brep_sign_and_distance() {
    let text = trimmed_loft();
    assert!(text.contains("LENGTH_MEASURE(1.E-07)"));
    let bytes = text.as_bytes();
    let hash = catallax_import::sha256_hex(bytes);
    let id = StepId::new("loft").unwrap();
    let recipe = Recipe::new(
        "Trimmed OpenCASCADE loft",
        vec![Step {
            id: id.clone(),
            op: Operation::ImportStep {
                file: FileRef {
                    name: "trimmed-loft.step".into(),
                    sha256: hash.clone(),
                },
                pose: Pose::identity(),
                source_solid: None,
                decisions: Default::default(),
            },
        }],
        id,
    );
    let sources = ImportFiles::from([(hash, bytes)]);
    let (brep, work) = with_operation_budget(WORK, || {
        evaluate_with_files(&recipe, &sources, &mut Decider::default())
    });
    let mut brep = brep.unwrap();
    assert!(!work.exhausted);
    assert!(brep.model.check(brep.result).is_empty());
    let brep_work = work.used;

    let source = catallax_import::read_step(bytes).unwrap();
    let source_face = &source.solids[0].faces[0];
    let catallax_import::Surface::BSpline(source_patch) = &source_face.surface else {
        panic!("the selected STEP face must carry a spline");
    };
    assert_eq!(
        source_patch.source_format,
        catallax_import::SourceFormat::Step
    );
    assert_eq!(*source_patch.geometry.domain()[1][0], q("-1"));

    let first = brep.model.faces_of_body(brep.result)[0];
    let Surface::BSpline(SplineSupport::RationalPatch(patch)) =
        brep.model.geom().surface(brep.model.face(first).surface)
    else {
        panic!("the imported B-rep must retain a rational patch");
    };
    assert_eq!(patch.domain()[1], [q("0"), q("1")]);

    let (field, work) = with_operation_budget(WORK, || {
        Plan::compile_with_files_checked(&recipe, &sources, &mut Decider::default())
    });
    let field = field.unwrap();
    assert!(!work.exhausted);
    assert!(field.result_patch_shell().is_some());
    println!("build work. B-rep {brep_work}. field {}", work.used);

    let near_patch = 69.0_f64.sqrt() / 10_000.0;
    // The chosen patch point is [3/4, 1/2, 1], with normal [8, 2, 1].
    // Write the patch's displacement there in centred unit coordinates s, t.
    // Its y/z displacement squared is >= 2s² + 3t², while the normal dot
    // displacement is -8st. The small normal offsets have that point as
    // their unique closest patch point. The other side faces are separated
    // by x + z/2 <= 1, y <= 0, or x <= 0, and each cap is one unit away.
    let witnesses = [
        (point(["0", "0", "1/100"]), Location::Inside, 0.01),
        (point(["0", "0", "-1/100"]), Location::Outside, 0.01),
        (
            point(["3746/5000", "2499/5000", "9999/10000"]),
            Location::Inside,
            near_patch,
        ),
        (
            point(["3754/5000", "2501/5000", "10001/10000"]),
            Location::Outside,
            near_patch,
        ),
        // This point lies on the discarded extension. The closest admitted
        // point is on the bottom edge at [1, -1/4, 0]. The entire B-rep is
        // on the x <= 1, z >= 0 side of those two planes.
        (
            point(["9/8", "-1/4", "-1/2"]),
            Location::Outside,
            17.0_f64.sqrt() / 8.0,
        ),
    ];
    let mut workspace = Workspace::default();
    for (point, expected, reference_distance) in witnesses {
        let (brep_sign, brep_work) = with_operation_budget(WORK, || {
            classify_point(
                &mut brep.model,
                brep.result,
                &point,
                &mut Decider::default(),
            )
        });
        assert!(!brep_work.exhausted);
        let brep_sign = match brep_sign {
            Outcome::Proven(Containment::Inside) => Location::Inside,
            Outcome::Proven(Containment::Outside) => Location::Outside,
            other => panic!("B-rep did not classify {point:?}: {other:?}"),
        };
        assert_eq!(brep_sign, expected, "{point:?}");
        let (query, work) = with_operation_budget(WORK, || {
            field.proven_checked(&point, &mut workspace, &mut Decider::default())
        });
        let query = query.unwrap();
        assert!(!work.exhausted);
        assert_eq!(query.location, brep_sign, "{point:?}");
        let Value::TrueDistance(enclosure) = query.value else {
            panic!("distance must be an enclosure, got {:?}", query.value);
        };
        let reference = if expected == Location::Inside {
            -reference_distance
        } else {
            reference_distance
        };
        assert!(enclosure.lo() <= reference && reference <= enclosure.hi());
        assert!((enclosure.hi() - enclosure.lo()) <= TOLERANCE_MM);
        println!("{point:?}: {expected:?}, {enclosure:?}, work {}", work.used);
    }
}
