//! The square tube with a spline top, written to STEP: a file whose u/v
//! curves the gap rule for imports reads. Shared by `check_once.rs` and by
//! the import crate's own test of an edge's work cap.
#[path = "../../../catallax-step/tests/common/n03_tube.rs"]
pub mod tube;

use catallax_numbers::{Decider, Rational};
use catallax_step::numerical::{self, NumericalStepOptions};
use catallax_topology::numerical::{self as n, NumericalBody, NumericalOptions};
use tube::q;

fn options() -> NumericalOptions {
    NumericalOptions {
        tolerance_mm: q("0.01"),
        vertex_tolerance_mm: None,
    }
}

/// The body's STEP file and the uncertainty it records.
pub fn written(body: &NumericalBody) -> (Vec<u8>, Rational) {
    let file = numerical::build(
        body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    (file.bytes().to_vec(), file.report().uncertainty_mm.clone())
}

/// The same file as another system's: the stamp no longer names a Catallax
/// numerical export, so it records nothing.
pub fn unstamped(bytes: &[u8]) -> Vec<u8> {
    let text = std::str::from_utf8(bytes).unwrap();
    assert!(text.contains("catallax numerical approximate;"));
    text.replacen(
        "catallax numerical approximate;",
        "another system approximate;",
        1,
    )
    .into_bytes()
}

/// The square tube with its top face a flat spline patch a micrometre above
/// its edges. A spline face with a real gap is taken by neither the certified
/// import nor the analytic one, so the file opens as a native tolerant body,
/// the way the R2 rim does.
pub fn spline_topped_tube() -> NumericalBody {
    use catallax_shapes::{BSplineSurface, Point3, RationalPatch, Surface};
    let mut definition = tube::tube();
    let mut dec = Decider::default();
    let lifted = catallax_numbers::Real::from(q("1.000001"));
    let controls = [(0, 0), (0, 4), (4, 0), (4, 4)]
        .map(|(x, y)| Point3::new(x.into(), y.into(), lifted.clone()))
        .to_vec();
    let knots = vec![q("0"), q("0"), q("1"), q("1")];
    let surface = BSplineSurface::new(
        [1, 1],
        [knots.clone(), knots],
        [2, 2],
        controls,
        vec![catallax_numbers::Real::one(); 4],
        &mut dec,
    )
    .unwrap();
    definition.faces[0].surface =
        Surface::BSpline(RationalPatch::new(surface, &mut dec).unwrap().into());
    for trim in definition.faces[0].loops.iter_mut().flatten() {
        for uv in trim.uv.iter_mut().flatten() {
            *uv = uv.clone() / q("4");
        }
    }
    n::build(definition, options(), &mut dec).unwrap()
}

/// One straight u/v curve of the spline face, written again as a quadratic
/// that traces the same segment at another speed: its middle control sits
/// four tenths of the way along, so the copy lags its edge by up to a
/// twentieth of the segment at the same parameter.
pub fn with_one_uv_curve_at_another_speed(text: &str) -> String {
    let curve = text
        .lines()
        .find(|line| {
            line.contains("=B_SPLINE_CURVE_WITH_KNOTS('',1,") && line.contains("(2,2),(0.,1.),")
        })
        .expect("a straight u/v curve over the edge's interval");
    let controls = curve
        .split("'',1,(")
        .nth(1)
        .and_then(|rest| rest.split(')').next())
        .expect("two controls");
    let ends: Vec<&str> = controls.split(',').collect();
    assert_eq!(ends.len(), 2, "{curve}");
    let uv = |name: &str| -> [f64; 2] {
        let line = text
            .lines()
            .find(|line| line.starts_with(&format!("{name}=CARTESIAN_POINT(")))
            .expect("a control point");
        let inner = line
            .rsplit("(")
            .next()
            .and_then(|rest| rest.split(')').next())
            .expect("coordinates");
        let values: Vec<f64> = inner.split(',').map(|x| x.parse().unwrap()).collect();
        [values[0], values[1]]
    };
    let (first, last) = (uv(ends[0]), uv(ends[1]));
    let middle = [0, 1].map(|axis| first[axis] + 0.4 * (last[axis] - first[axis]));
    let next = text
        .lines()
        .filter_map(|line| {
            line.strip_prefix('#')?
                .split('=')
                .next()?
                .parse::<u64>()
                .ok()
        })
        .max()
        .expect("an entity")
        + 1;
    let quadratic = curve
        .replace(
            &format!("'',1,({controls})"),
            &format!("'',2,({},#{next},{})", ends[0], ends[1]),
        )
        .replace("(2,2),(0.,1.),", "(3,3),(0.,1.),");
    assert_ne!(quadratic, curve);
    text.replace(
        curve,
        &format!(
            "#{next}=CARTESIAN_POINT('',({:.17},{:.17}));\n{quadratic}",
            middle[0], middle[1]
        ),
    )
}

/// The last straight u/v curve of the spline face, written over half its
/// edge's interval. The gap rule for imports does not read a u/v curve whose
/// interval does not hold its edge's parameter range, so as written the body
/// is refused for it.
#[allow(dead_code)]
pub fn with_one_uv_curve_over_half_its_interval(text: &str) -> String {
    let curve = text
        .lines()
        .rfind(|line| {
            line.contains("=B_SPLINE_CURVE_WITH_KNOTS('',1,") && line.contains("(2,2),(0.,1.),")
        })
        .expect("a straight u/v curve over the edge's interval");
    let halved = curve.replace("(2,2),(0.,1.),", "(2,2),(0.,0.5),");
    assert_ne!(halved, curve);
    text.replace(curve, &halved)
}
