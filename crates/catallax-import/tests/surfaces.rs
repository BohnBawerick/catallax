mod common;

use catallax_boolean::{BoolOp, boolean};
use catallax_import::{ImportRefusal, build_file, canonical_bytes, read_step};
use catallax_measure::{Accuracy, area, bounding_box, volume};
use catallax_numbers::{Decider, Rational};
use catallax_topology::{Body, BodyKind, Model, Shell};
use common::{BoxFixture, box_step};

const SQUARE: &str = include_str!("fixtures/surface-square.step");

fn build(source: &str) -> (Model, Vec<catallax_topology::BodyId>) {
    let source = read_step(source.as_bytes()).unwrap();
    let mut model = Model::new();
    let bodies = build_file(&mut model, &source, None).unwrap();
    for body in &bodies {
        assert!(model.check(*body).is_empty());
    }
    (model, bodies)
}

fn round_trip(model: &Model, bodies: &[catallax_topology::BodyId]) -> catallax_import::StepFile {
    let mut bytes = Vec::new();
    let report = catallax_step::write_step(
        model,
        bodies,
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    assert_eq!(
        report.value.solids,
        bodies
            .iter()
            .filter(|body| model.body(**body).kind == BodyKind::Solid)
            .count()
    );
    catallax_step::check_seams(&bytes).unwrap();
    let source = read_step(&bytes).unwrap();
    let mut reimported = Model::new();
    let result = build_file(&mut reimported, &source, None).unwrap();
    assert_eq!(result.len(), bodies.len());
    for (a, b) in bodies.iter().zip(result) {
        assert_eq!(model.body(*a).kind, reimported.body(b).kind);
        assert_eq!(model.counts(*a), reimported.counts(b));
        assert!(reimported.check(b).is_empty());
    }
    source
}

#[test]
fn open_surface_has_area_and_bounds_but_refuses_volume_and_every_boolean_shortcut() {
    let file = read_step(SQUARE.as_bytes()).unwrap();
    assert!(file.solids.is_empty());
    assert_eq!(file.surface_models.len(), 1);
    assert!(!file.surface_models[0].shells[0].closed);
    let (mut model, bodies) = build(SQUARE);
    let body = bodies[0];
    assert_eq!(model.body(body).kind, BodyKind::Surface);
    let measured = area(&model, body, &Accuracy::default(), &mut Decider::default()).unwrap();
    assert_eq!(measured.value.lo, Rational::from_i64(6));
    assert_eq!(measured.value.hi, Rational::from_i64(6));
    let bounds = bounding_box(&model, body, &Accuracy::default(), &mut Decider::default())
        .unwrap()
        .value
        .unwrap();
    assert_eq!(bounds.min, std::array::from_fn(|_| Rational::zero()));
    assert_eq!(
        bounds.max,
        [
            Rational::from_i64(2),
            Rational::from_i64(3),
            Rational::zero()
        ]
    );
    assert_eq!(
        volume(&model, body, &Accuracy::default(), &mut Decider::default())
            .unwrap_err()
            .code,
        "measure.not-solid"
    );
    let empty = model.add_body(Body {
        kind: BodyKind::Empty,
        lumps: Vec::new(),
        marks: Default::default(),
    });
    for op in [BoolOp::Unite, BoolOp::Subtract, BoolOp::Intersect] {
        for (a, b) in [(body, body), (body, empty), (empty, body)] {
            assert_eq!(
                boolean(&mut model, a, b, op, &mut Decider::default())
                    .unwrap_err()
                    .refusal
                    .code,
                "boolean.not-solid"
            );
        }
    }
    let exported = round_trip(&model, &bodies);
    assert!(exported.solids.is_empty());
    assert!(!exported.surface_models[0].shells[0].closed);
}

#[test]
fn closed_shell_declared_as_a_surface_stays_a_surface() {
    let source = box_step(BoxFixture::new(["1", "2", "3"]));
    let source = source
        .lines()
        .map(|line| {
            if line.contains("MANIFOLD_SOLID_BREP") {
                let (id, _) = line.split_once('=').unwrap();
                let shell = line.split(",#").nth(1).unwrap().trim_end_matches(");");
                format!("{id}=SHELL_BASED_SURFACE_MODEL('',(#{shell}));")
            } else {
                line.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    let (model, bodies) = build(&source);
    assert_eq!(model.body(bodies[0]).kind, BodyKind::Surface);
    assert_eq!(
        area(
            &model,
            bodies[0],
            &Accuracy::default(),
            &mut Decider::default()
        )
        .unwrap()
        .value
        .lo,
        Rational::from_i64(22)
    );
    assert_eq!(
        volume(
            &model,
            bodies[0],
            &Accuracy::default(),
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "measure.not-solid"
    );
    let source = round_trip(&model, &bodies);
    assert!(source.surface_models[0].shells[0].closed);
}

#[test]
fn two_connected_faces_form_an_open_surface_shell() {
    let solid = box_step(BoxFixture::new(["1", "2", "3"]));
    let source = solid
        .lines()
        .map(|line| {
            if let Some((id, faces)) = line.split_once("CLOSED_SHELL('',(") {
                let faces: Vec<_> = faces.strip_suffix("));").unwrap().split(',').collect();
                format!("{id}OPEN_SHELL('',({},{}));", faces[0], faces[2])
            } else if line.contains("MANIFOLD_SOLID_BREP") {
                let (id, _) = line.split_once('=').unwrap();
                let shell = line.split(",#").nth(1).unwrap().trim_end_matches(");");
                format!("{id}=SHELL_BASED_SURFACE_MODEL('',(#{shell}));")
            } else {
                line.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    let (model, bodies) = build(&source);
    assert_eq!(model.counts(bodies[0]).faces, 2);
    assert_eq!(model.body(bodies[0]).kind, BodyKind::Surface);
    round_trip(&model, &bodies);
}

#[test]
fn multiple_surface_shells_are_kept_and_affect_the_canonical_definition() {
    let source = SQUARE.replace(
        "#32=SHELL_BASED_SURFACE_MODEL('',(#31));",
        "#36=OPEN_SHELL('',(#30));\n#32=SHELL_BASED_SURFACE_MODEL('',(#31,#36));",
    );
    let file = read_step(source.as_bytes()).unwrap();
    assert_ne!(
        canonical_bytes(&file),
        canonical_bytes(&read_step(SQUARE.as_bytes()).unwrap())
    );
    let (model, bodies) = build(&source);
    assert_eq!(model.body(bodies[0]).lumps.len(), 2);
    assert_eq!(model.faces_of_body(bodies[0]).len(), 2);
    assert_eq!(
        round_trip(&model, &bodies).surface_models[0].shells.len(),
        2
    );
}

#[test]
fn a_mixed_file_keeps_separate_solid_and_surface_bodies() {
    let surface_data = SQUARE
        .split("DATA;\n")
        .nth(1)
        .unwrap()
        .split("ENDSEC;")
        .next()
        .unwrap()
        .replace('#', "#900");
    let solid = box_step(BoxFixture::new(["1", "2", "3"]));
    let source = solid.replace(
        "ENDSEC;\nEND-ISO-10303-21;",
        &format!("{surface_data}ENDSEC;\nEND-ISO-10303-21;"),
    );
    assert_ne!(source, solid);
    let (model, bodies) = build(&source);
    assert_eq!(bodies.len(), 2);
    assert_eq!(model.body(bodies[0]).kind, BodyKind::Solid);
    assert_eq!(model.body(bodies[1]).kind, BodyKind::Surface);
    let exported = round_trip(&model, &bodies);
    assert_eq!(exported.solids.len(), 1);
    assert_eq!(exported.surface_models.len(), 1);
}

#[test]
fn an_annulus_retains_both_boundary_cycles() {
    let extra = "#40=CIRCLE('',#8,2.);\n#41=CIRCLE('',#8,1.);\n#42=VERTEX_POINT('',#2);\n#43=CARTESIAN_POINT('',(1.,0.,0.));\n#44=VERTEX_POINT('',#43);\n#45=EDGE_CURVE('',#42,#42,#40,.T.);\n#46=EDGE_CURVE('',#44,#44,#41,.T.);\n#47=ORIENTED_EDGE('',*,*,#45,.T.);\n#48=ORIENTED_EDGE('',*,*,#46,.F.);\n#49=EDGE_LOOP('',(#47));\n#50=EDGE_LOOP('',(#48));\n#51=FACE_OUTER_BOUND('',#49,.T.);\n#52=FACE_BOUND('',#50,.T.);\n#53=ADVANCED_FACE('',(#51,#52),#9,.T.);\n";
    let source = SQUARE
        .replace("#31=OPEN_SHELL('',(#30));", "#31=OPEN_SHELL('',(#53));")
        .replace(
            "ENDSEC;\nEND-ISO-10303-21;",
            &format!("{extra}ENDSEC;\nEND-ISO-10303-21;"),
        );
    let (model, bodies) = build(&source);
    let measured = area(
        &model,
        bodies[0],
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert!(measured.lo >= Rational::parse("9.42").unwrap());
    assert!(measured.hi <= Rational::parse("9.43").unwrap());
    let source = round_trip(&model, &bodies);
    assert_eq!(source.surface_models[0].shells[0].faces[0].loops.len(), 2);
}

#[test]
fn surface_topology_rejects_nonmanifold_edge_use_and_empty_models() {
    let (mut model, bodies) = build(SQUARE);
    let face = model.faces_of_body(bodies[0])[0];
    let shell = model.add_shell(Shell {
        faces: vec![face; 3],
    });
    let lump = model.add_lump(catallax_topology::Lump {
        shells: vec![shell],
    });
    let invalid = model.add_body(Body {
        kind: BodyKind::Surface,
        lumps: vec![lump],
        marks: Default::default(),
    });
    assert!(
        model
            .check(invalid)
            .iter()
            .any(|fault| fault.code == "topology.edge-use")
    );
    for source in [
        SQUARE.replace("(#31));", "());"),
        SQUARE.replace("#31=OPEN_SHELL('',(#30));", "#31=OPEN_SHELL('',());"),
    ] {
        let file = read_step(source.as_bytes()).unwrap();
        assert!(matches!(
            build_file(&mut Model::new(), &file, None),
            Err(ImportRefusal::NotSolid { .. })
        ));
    }
}
