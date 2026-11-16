//! Material endpoints, source units and transactional refusal at a crossed face.
mod common;
use catallax_import::{
    Surface,
    direct::{DirectSolid, Edit},
    numerical::ImportError,
    sha256_hex,
    tolerant::opened,
};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::StepId;
use common::{BoxFixture, Unit, box_step, r};

fn source() -> DirectSolid {
    let mut fixture = BoxFixture::new(["0.002", "0.003", "0.004"]);
    fixture.unit = Unit::Metre;
    let text = box_step(fixture);
    let opened = opened::open_step(
        text.as_bytes(),
        &sha256_hex(text.as_bytes()),
        &StepId::new("direct-source").unwrap(),
        Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    DirectSolid::from_opened(std::sync::Arc::new(opened), 0, &mut Decider::default()).unwrap()
}
fn positive_x(solid: &DirectSolid) -> String {
    solid
        .faces()
        .iter()
        .find(|f| {
            matches!(&f.surface,
        Surface::Plane(p) if p.origin[0] == r("2") && !p.axis[0].is_zero())
        })
        .unwrap()
        .name
        .clone()
}
fn max_x(solid: &DirectSolid) -> String {
    let mut dec = Decider::default();
    let bytes = match solid.body() {
        opened::BodyView::Certified { model, body } => {
            let mut bytes = Vec::new();
            catallax_step::write_step(model, &[body], &Default::default(), &mut dec, &mut bytes)
                .unwrap();
            bytes
        }
        opened::BodyView::Numerical(body) => {
            catallax_step::numerical::build(body, &Default::default(), &mut dec)
                .unwrap()
                .bytes()
                .to_vec()
        }
        _ => panic!("material required"),
    };
    let file = catallax_import::read_step_with_decider(&bytes, &mut dec).unwrap();
    file.solids[0]
        .faces
        .iter()
        .flat_map(|f| &f.loops)
        .flat_map(|b| &b.edges)
        .flat_map(|e| [&e.start.point[0], &e.end.point[0]])
        .max()
        .unwrap()
        .to_string()
}

#[test]
fn imported_metre_units_do_not_scale_a_millimetre_edit() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let before = source();
        let after = before
            .apply(
                Edit::OffsetFace {
                    face: positive_x(&before),
                    distance_mm: r("1"),
                },
                Default::default(),
                Default::default(),
                &mut Decider::default(),
            )
            .unwrap();
        assert_eq!(max_x(&before), "2");
        assert_eq!(max_x(&after), "3");
        assert_eq!(after.source_faces(), before.faces());
        assert_eq!(after.source_sha256(), before.source_sha256());
        assert!(before.history().is_empty());
    });
    assert!(!work.exhausted);
}

#[test]
fn crossing_the_opposite_support_refuses_without_corrupting_the_input() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let before = source();
        let key = positive_x(&before);
        let refusal = before
            .apply(
                Edit::OffsetFace {
                    face: key.clone(),
                    distance_mm: r("-3"),
                },
                Default::default(),
                Default::default(),
                &mut Decider::default(),
            )
            .unwrap_err();
        assert!(matches!(refusal, ImportError::Numerical(r) if r.code == "direct.topology-change"));
        let after = before
            .apply(
                Edit::MoveFace {
                    face: key,
                    translation_mm: [r("1"), r("0"), r("0")],
                },
                Default::default(),
                Default::default(),
                &mut Decider::default(),
            )
            .unwrap();
        assert_eq!(max_x(&before), "2");
        assert_eq!(max_x(&after), "3");
        assert_eq!(before.faces(), before.source_faces());
    });
    assert!(!work.exhausted);
}

#[test]
fn inherited_line_gap_is_bounded_and_reconstructed_before_step_export() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut text = box_step(BoxFixture::new(["2", "3", "4"]));
        let line = text
            .lines()
            .find(|line| line.contains("=LINE("))
            .unwrap()
            .to_owned();
        let (key, args) = line.split_once("=LINE(").unwrap();
        let args: Vec<_> = args.strip_suffix(");").unwrap().split(',').collect();
        let origin_key = args[1].trim();
        let origin = text
            .lines()
            .find(|line| line.starts_with(&format!("{origin_key}=")))
            .unwrap();
        let coordinates = origin
            .split_once(",(")
            .unwrap()
            .1
            .strip_suffix("));")
            .unwrap();
        let coordinates: Vec<_> = coordinates.split(',').map(r).collect();
        let next = 1 + text
            .lines()
            .filter_map(|line| {
                line.strip_prefix('#')?
                    .split_once('=')?
                    .0
                    .parse::<u64>()
                    .ok()
            })
            .max()
            .unwrap();
        let noise = [r("0.000001"), r("0.000002"), r("0.000003")];
        let point = format!(
            "#{next}=CARTESIAN_POINT('',({},{},{}));",
            common::dec(&(&coordinates[0] + &noise[0])),
            common::dec(&(&coordinates[1] + &noise[1])),
            common::dec(&(&coordinates[2] + &noise[2]))
        );
        text = text.replace(
            &line,
            &format!("{key}=LINE({},#{next},{});", args[0], args[2]),
        );
        text = text.replace(
            "ENDSEC;\nEND-ISO-10303-21;",
            &format!("{point}\nENDSEC;\nEND-ISO-10303-21;"),
        );
        let opened = opened::open_step(
            text.as_bytes(),
            &sha256_hex(text.as_bytes()),
            &StepId::new("line-gap").unwrap(),
            Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        let before =
            DirectSolid::from_opened(std::sync::Arc::new(opened), 0, &mut Decider::default())
                .unwrap();
        let cap = before
            .faces()
            .iter()
            .find(|face| {
                matches!(&face.surface,
            Surface::Plane(p) if p.origin[2]==r("4") && !p.axis[2].is_zero())
            })
            .unwrap()
            .name
            .clone();
        let after = before
            .apply(
                Edit::OffsetFace {
                    face: cap,
                    distance_mm: r("1"),
                },
                Default::default(),
                Default::default(),
                &mut Decider::default(),
            )
            .unwrap();
        let repairs = &after.history()[0].line_repairs;
        assert!(
            repairs
                .iter()
                .any(|repair| !repair.error_upper_mm.is_zero())
        );
        assert!(
            repairs
                .iter()
                .all(|repair| repair.error_upper_mm <= r("0.01"))
        );
        let opened::BodyView::Numerical(body) = after.body() else {
            panic!("numerical edit required")
        };
        let file =
            catallax_step::numerical::build(body, &Default::default(), &mut Decider::default())
                .unwrap();
        file.check_written(file.bytes(), &mut Decider::default())
            .unwrap();
        let emitted =
            catallax_import::read_step_with_decider(file.bytes(), &mut Decider::default()).unwrap();
        let highest = emitted.solids[0]
            .faces
            .iter()
            .flat_map(|f| &f.loops)
            .flat_map(|b| &b.edges)
            .flat_map(|e| [&e.start.point[2], &e.end.point[2]])
            .max()
            .unwrap();
        assert_eq!(highest, &r("5"));
        for edge in emitted.solids[0]
            .faces
            .iter()
            .flat_map(|f| &f.loops)
            .flat_map(|b| &b.edges)
        {
            let catallax_import::Curve::Line(origin, direction) = &edge.curve else {
                panic!("line box")
            };
            for vertex in [&edge.start, &edge.end] {
                let delta: [_; 3] = std::array::from_fn(|i| &vertex.point[i] - &origin[i]);
                for i in 0..3 {
                    let j = (i + 1) % 3;
                    assert_eq!(&delta[i] * &direction[j], &delta[j] * &direction[i]);
                }
            }
        }
        assert_eq!(before.source_faces(), after.source_faces());
    });
    assert!(!work.exhausted);
}
