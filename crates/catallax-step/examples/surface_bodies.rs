//! Trims and extends a curved spline sheet, then thickens it and writes actual STEP.
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_shapes::spline_tools::RationalSplineSurface;
use catallax_step::numerical::{self, NumericalStepOptions};
use catallax_topology::{
    BodyKind,
    numerical::{NumericalBody, NumericalOptions, surface::SplineSurfaceBody},
};
use serde_json::json;
use std::{fs, path::Path};
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn export(out: &Path, name: &str, body: &NumericalBody, dec: &mut Decider) -> serde_json::Value {
    let artifact = numerical::build(body, &NumericalStepOptions::default(), dec).unwrap();
    let path = out.join(format!("{name}.step"));
    fs::write(&path, artifact.bytes()).unwrap();
    artifact
        .check_written(&fs::read(&path).unwrap(), dec)
        .unwrap();
    assert_eq!(artifact.report().kind, body.kind());
    let count = body.counts();
    json!({"name":name,"kind":format!("{:?}",body.kind()),"faces":count.faces,"edges":count.edges,"vertices":count.vertices,
        "tolerance_mm":body.options().tolerance_mm.to_string(),"history":body.history().len(),"checks":body.evidence().checks().len(),"bytes":artifact.bytes().len()})
}
fn run(out: &Path) {
    fs::create_dir_all(out).unwrap();
    let mut dec = Decider::default();
    let mut controls = Vec::new();
    for i in 0..3 {
        for j in 0..3 {
            let u = Rational::from_ratio(i, 2).unwrap();
            let v = Rational::from_ratio(j, 2).unwrap();
            let bow = if j == 1 { q("1/20") } else { q("0") };
            controls.push([
                &u * &q("10"),
                &v * &q("8"),
                &u * &q("1/50") + bow + &u * &v * q("1/100"),
            ]);
        }
    }
    let knots = vec![q("0"), q("0"), q("0"), q("1"), q("1"), q("1")];
    let surface = RationalSplineSurface::new(
        [2, 2],
        [knots.clone(), knots],
        [3, 3],
        controls,
        vec![q("1"); 9],
        &dec,
    )
    .unwrap();
    let sheet = SplineSurfaceBody::new(
        "curved-panel",
        surface,
        NumericalOptions {
            tolerance_mm: q("1/100000"),
        },
        &mut dec,
    )
    .unwrap();
    assert_eq!(sheet.body().kind(), BodyKind::Surface);
    let trimmed = sheet
        .trim([[q("1/5"), q("4/5")], [q("1/4"), q("3/4")]], &mut dec)
        .unwrap();
    assert_eq!(
        trimmed.surface().domain(),
        [[&q("1/5"), &q("4/5")], [&q("1/4"), &q("3/4")]]
    );
    let extended = trimmed.extend(0, q("0"), &mut dec).unwrap();
    assert_eq!(extended.surface().domain()[0], [&q("0"), &q("4/5")]);
    let solid = extended.thicken(q("2/5"), &mut dec).unwrap();
    assert_eq!(solid.kind(), BodyKind::Solid);
    assert_eq!(solid.counts().faces, 6);
    assert_eq!(solid.counts().edges, 12);
    assert!(
        solid
            .definition()
            .faces
            .iter()
            .any(|face| face.key == "curved-panel")
    );
    let exports = vec![
        export(out, "source-sheet", sheet.body(), &mut dec),
        export(out, "trimmed-sheet", trimmed.body(), &mut dec),
        export(out, "extended-sheet", extended.body(), &mut dec),
        export(out, "thickened-solid", &solid, &mut dec),
    ];
    let report = json!({"format":"catallax-surface-body-proof/1","checks":"surface kind retained; exact parameter trim; tangential extension source samples; normal offset samples; closed thickened topology; retained source face key","exports":exports});
    fs::write(
        out.join("surface-body-report.json"),
        serde_json::to_vec_pretty(&report).unwrap(),
    )
    .unwrap();
    println!("{report}");
}
fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: surface_bodies OUTPUT_DIRECTORY");
    let (_, work) = with_operation_budget(100_000_000, || run(Path::new(&path)));
    assert!(!work.exhausted);
    println!("work={}", work.used);
}
