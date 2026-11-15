//! Missing cylinder seams are explicit topology repairs, retaining source rims.
#[path = "../../catallax-topology/tests/common/periodic.rs"]
mod fixture;
use catallax_import::{
    numerical::{self, StepImport},
    sha256_hex,
};
use catallax_numbers::Decider;
use catallax_recipe::StepId;
use catallax_step::numerical::{self as output, NumericalStepOptions};
use catallax_topology::numerical as n;

fn refs(line: &str) -> Vec<u64> {
    line.split_once('=')
        .unwrap()
        .1
        .split('#')
        .skip(1)
        .map(|s| {
            s.chars()
                .take_while(char::is_ascii_digit)
                .collect::<String>()
                .parse()
                .unwrap()
        })
        .collect()
}
fn source(delta: &str) -> Vec<u8> {
    let body = n::build(
        fixture::bore(),
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
    let mut text = std::str::from_utf8(file.bytes()).unwrap().to_owned();
    let line = |id| {
        text.lines()
            .find(|s| s.starts_with(&format!("#{id}=")))
            .unwrap()
            .to_owned()
    };
    let face = text
        .lines()
        .find(|s| s.contains("=ADVANCED_FACE('bore-wall',"))
        .unwrap()
        .to_owned();
    let ids = refs(&face);
    let bound = line(ids[0]);
    let edge_loop = line(refs(&bound)[0]);
    let coedges = refs(&edge_loop);
    let max = text
        .lines()
        .filter_map(|s| s.strip_prefix('#'))
        .map(|s| s.split_once('=').unwrap().0.parse::<u64>().unwrap())
        .max()
        .unwrap();
    let replaced_loop = format!(
        "{}=EDGE_LOOP('',(#{}));",
        edge_loop.split_once('=').unwrap().0,
        coedges[0]
    );
    let replaced_face = face.replace(
        &format!("(#{})", ids[0]),
        &format!("(#{},#{})", ids[0], max + 2),
    );
    let extra = format!(
        "#{}=EDGE_LOOP('',(#{}));\n#{}=FACE_OUTER_BOUND('',#{},.T.);\n",
        max + 1,
        coedges[2],
        max + 2,
        max + 1
    );
    let plane = text.lines().find(|s| s.contains("=PLANE('f2',")).unwrap();
    let placement = line(refs(plane)[0]);
    let origin = line(refs(&placement)[0]);
    text = text
        .replace(&edge_loop, &replaced_loop)
        .replace(&face, &replaced_face)
        .replace(
            &origin,
            &format!(
                "{}=CARTESIAN_POINT('',(0.,{delta},0.));",
                origin.split_once('=').unwrap().0
            ),
        );
    let end = text.rfind("ENDSEC;").unwrap();
    text.insert_str(end, &extra);
    text.into_bytes()
}
fn import(bytes: &[u8]) -> Result<StepImport, numerical::ImportError> {
    numerical::import_step(
        bytes,
        &sha256_hex(bytes),
        &StepId::new("periodic-source").unwrap(),
        &n::NumericalOptions::default(),
        &mut Decider::default(),
    )
}
#[test]
fn missing_seam_joins_retained_rims_and_records_the_derived_edge() {
    let bytes = source("-1.99999999");
    let before = bytes.clone();
    let recognized = catallax_import::read_step(&bytes).unwrap();
    let source_vertices = recognized.solids[0]
        .faces
        .iter()
        .flat_map(|f| &f.loops)
        .flat_map(|b| &b.edges)
        .flat_map(|e| [&e.start, &e.end])
        .map(|v| (v.name.clone(), v.point.clone()))
        .collect::<std::collections::BTreeMap<_, _>>();
    let StepImport::Numerical(result) = import(&bytes).unwrap() else {
        panic!("expected numerical admission")
    };
    let record = result.source_admission().unwrap();
    assert_eq!(record.topology_repairs.len(), 1);
    assert_eq!(record.intersection_checks.len(), 21);
    let repair = &record.topology_repairs[0];
    assert_eq!(repair.method, "paired_cylinder_generatrix/1");
    let body = &result.bodies()[0].1;
    assert_eq!(
        (
            body.counts().faces,
            body.counts().edges,
            body.counts().vertices
        ),
        (7, 15, 10)
    );
    assert_eq!(
        body.definition()
            .vertices
            .iter()
            .map(|v| (v.key.clone(), v.position.clone()))
            .collect::<std::collections::BTreeMap<_, _>>(),
        source_vertices
    );
    assert!(
        body.definition()
            .edges
            .iter()
            .any(|e| e.key == repair.derived_seam_edge)
    );
    assert_eq!(bytes, before);
    output::build(
        body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
}
#[test]
fn over_cap_source_incidence_and_exhausted_requests_return_no_body() {
    assert!(import(&source("-1.98")).is_err());
    let bytes = source("-1.99999999");
    let (result, work) = catallax_numbers::with_operation_budget(0, || import(&bytes));
    assert!(result.is_err() && work.exhausted);
}

#[test]
fn separated_rim_anchors_do_not_trigger_nearest_seam_guessing() {
    let bytes = source("-1.99999999");
    let text = std::str::from_utf8(&bytes).unwrap();
    let vertex = text
        .lines()
        .find(|s| s.contains("=VERTEX_POINT('v9',"))
        .unwrap();
    let point = refs(vertex)[0];
    let old = text
        .lines()
        .find(|s| s.starts_with(&format!("#{point}=")))
        .unwrap();
    let changed = text.replace(old, &format!("#{point}=CARTESIAN_POINT('',(0.,1.,1.));"));
    assert!(import(changed.as_bytes()).is_err());
}
