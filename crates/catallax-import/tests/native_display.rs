//! Display carriers sample in binary64 and say when binary64 is not enough.
use catallax_import::{
    geometry_source::{AngleUnit, ParsedStep},
    tolerant::native::NativeSurface,
};
use catallax_numbers::{Decider, Rational, with_operation_budget};

fn sphere(centre: &str) -> NativeSurface {
    let text = format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('display'),'2;1');\nFILE_NAME('display.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n#1=CARTESIAN_POINT('centre',({centre},0.,0.));\n#2=AXIS2_PLACEMENT_3D('frame',#1,$,$);\n#3=SPHERICAL_SURFACE('sphere',#2,2.5);\nENDSEC;\nEND-ISO-10303-21;\n"
    );
    let mut dec = Decider::default();
    let parsed = ParsedStep::read(text.as_bytes(), &mut dec).unwrap();
    NativeSurface::from_source(
        &parsed,
        3,
        &Rational::parse("1").unwrap(),
        AngleUnit::Radians,
        &mut dec,
    )
    .unwrap()
}

const STATIONS: [[f64; 2]; 4] = [[0.0, 0.0], [0.7, -0.4], [3.1, 1.2], [5.9, -1.5]];

#[test]
fn a_display_sample_is_the_enclosed_sample_at_a_fraction_of_its_work() {
    let native = sphere("0.");
    let display = native.for_display(&mut Decider::default()).unwrap();
    let (enclosed, enclosing) = with_operation_budget(1_000_000, || {
        STATIONS.map(|uv| {
            let mut dec = Decider::default();
            native
                .sample(native.location(uv).unwrap(), &mut dec)
                .unwrap()
        })
    });
    let (rounded, rounding) = with_operation_budget(1_000_000, || {
        STATIONS.map(|uv| {
            let mut dec = Decider::default();
            display
                .sample(display.location(uv).unwrap(), &mut dec)
                .unwrap()
        })
    });
    for (enclosed, rounded) in enclosed.iter().zip(&rounded) {
        for axis in 0..3 {
            let miss = (enclosed.point_mm[axis] - rounded.point_mm[axis]).abs();
            assert!(
                miss <= rounded.uncertainty_mm,
                "{miss} over {}",
                rounded.uncertainty_mm
            );
        }
        assert!(rounded.uncertainty_mm < 1e-13);
    }
    assert!(
        rounding.used * 10 < enclosing.used,
        "{} against {}",
        rounding.used,
        enclosing.used
    );
}

// At a thousand kilometres from the origin binary64 places a point of this
// sphere to about 2e-7 mm, coarser than the width projection stops at. The
// display carrier then answers with the enclosing evaluation.
#[test]
fn a_display_sample_binary64_cannot_place_is_enclosed() {
    let native = sphere("1000000000.");
    let display = native.for_display(&mut Decider::default()).unwrap();
    for uv in STATIONS {
        let mut dec = Decider::default();
        let enclosed = native
            .sample(native.location(uv).unwrap(), &mut dec)
            .unwrap();
        let shown = display
            .sample(display.location(uv).unwrap(), &mut dec)
            .unwrap();
        assert_eq!(
            enclosed.point_mm.map(f64::to_bits),
            shown.point_mm.map(f64::to_bits)
        );
        assert_eq!(
            enclosed.uncertainty_mm.to_bits(),
            shown.uncertainty_mm.to_bits()
        );
    }
}
