use catallax_import::{ImportRefusal, Surface, read_step};
use catallax_numbers::{Decider, Rational};

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn source(axis: &str) -> String {
    let text = std::str::from_utf8(include_bytes!("fixtures/spline-roof.step")).unwrap();
    text.lines()
        .map(|line| {
            if line.starts_with("#116=") {
                "#116=SURFACE_OF_LINEAR_EXTRUSION('',#45,#23);".into()
            } else if line.starts_with("#18=") {
                format!("#18=DIRECTION('',{axis});")
            } else {
                line.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn spline_extrusion_retains_profile_weights_and_the_whole_bounded_strip() {
    let file = read_step(source("(0.,2.,0.)").as_bytes()).unwrap();
    let Surface::BSpline(surface) = &file.solids[0].faces.last().unwrap().surface else {
        panic!("extrusion must map to a tensor spline");
    };
    assert_eq!(surface.geometry.degrees(), [2, 1]);
    assert_eq!(surface.geometry.control_count(), [3, 2]);
    assert_eq!(
        surface.geometry.weights(),
        &[q("1"), q("1"), q("2"), q("2"), q("1"), q("1")]
    );
    let mut dec = Decider::default();
    let carrier = surface.geometry.to_surface(&mut dec).unwrap();
    for (u, v, expected) in [
        ("0.1", "0", ["0", "0", "1"]),
        ("0.9", "1", ["1", "1", "1"]),
        ("0.5", "0.5", ["0.5", "0.5", "1"]),
    ] {
        let (u, v) = (q(u), q(v));
        let actual = carrier.evaluate_box([&u, &u], [&v, &v], 128, &dec).unwrap();
        for (actual, expected) in actual.iter().zip(expected.map(q)) {
            assert!(actual.lo() <= expected && actual.hi() >= expected);
        }
    }
}

#[test]
fn extrusion_zero_vector_is_malformed() {
    let source = source("(0.,1.,0.)").replace("#23=VECTOR('',#18,1.);", "#23=VECTOR('',#18,0.);");
    assert!(matches!(
        read_step(source.as_bytes()),
        Err(ImportRefusal::Malformed { .. })
    ));
}

#[test]
fn extrusion_without_a_finite_boundary_refuses_without_inventing_a_domain() {
    let source = source("(0.,1.,0.)").replace(
        "#117=ADVANCED_FACE('',(#106),#116,.T.);",
        "#117=ADVANCED_FACE('',(),#116,.T.);",
    );
    assert!(matches!(
        read_step(source.as_bytes()),
        Err(ImportRefusal::NotSolid { .. })
    ));
}
