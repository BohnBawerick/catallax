//! Free-wire geometry survives journal relocation and orientation-reversing maps.
use catallax_numbers::{Decider, MarkLog, Rational, Real};
use catallax_shapes::{
    Curve, Isometry, Pose, RPoint3, RVector3, spline_tools::RationalSplineCurve,
};
use catallax_topology::{Body, BodyId, BodyKind, Coedge, Model, Sense, Wire};

fn points(model: &Model, body: BodyId) -> Vec<[Rational; 3]> {
    let mut dec = Decider::default();
    model
        .vertices_of_body(body)
        .into_iter()
        .map(|vertex| {
            let point = model.geom().position(model.vertex(vertex).point);
            [&point.x, &point.y, &point.z].map(|value| dec.rational_value(value).unwrap().unwrap())
        })
        .collect()
}

#[test]
fn wire_replay_and_reflection_retain_the_directed_geometry() {
    let mut base = Model::new();
    let block = catallax_topology::block(
        &mut base,
        [2, 3, 5].map(Rational::from_i64),
        &Pose::identity(),
    )
    .unwrap();
    let edges = base.edges_of_body(block);
    let mut branch = base.clone();
    branch.journal_begin();
    let coedge = branch.add_coedge(Coedge {
        edge: edges[0],
        sense: Sense::Same,
    });
    let wire = branch.add_body(Body {
        kind: BodyKind::Wire,
        lumps: Vec::new(),
        wires: vec![Wire {
            coedges: vec![coedge],
        }],
        marks: MarkLog::new(),
    });
    let journal = branch.journal_end();
    base.add_coedge(Coedge {
        edge: edges[1],
        sense: Sense::Same,
    });
    let remap = base.replay(&journal);
    let wire = remap.body(wire);
    assert!(base.check(wire).is_empty());
    assert_eq!(
        points(&base, wire),
        vec![[0, 0, 0], [0, 0, 5]]
            .into_iter()
            .map(|p| p.map(Rational::from_i64))
            .collect::<Vec<_>>()
    );
    let reflection = Isometry::reflection(
        &RPoint3::origin(),
        &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
    )
    .unwrap();
    let mirrored = catallax_topology::transformed_body(&mut base, wire, &reflection, Clone::clone);
    assert!(base.check(mirrored).is_empty());
    assert_eq!(
        points(&base, mirrored),
        vec![[0, 0, 0], [0, 0, -5]]
            .into_iter()
            .map(|p| p.map(Rational::from_i64))
            .collect::<Vec<_>>()
    );
    assert_eq!(
        points(&base, wire),
        vec![[0, 0, 0], [0, 0, 5]]
            .into_iter()
            .map(|p| p.map(Rational::from_i64))
            .collect::<Vec<_>>()
    );
}

#[test]
fn nonclamped_wire_uses_exact_active_domain_endpoints() {
    let mut dec = Decider::default();
    let curve = RationalSplineCurve::new(
        2,
        (0..6).map(Rational::from_i64).collect(),
        [[0, 0, 0], [2, 0, 0], [4, 2, 0]]
            .map(|p| p.map(Rational::from_i64))
            .to_vec(),
        vec![Rational::one(); 3],
        &dec,
    )
    .unwrap();
    let (model, body) = catallax_topology::spline_wire(&[curve], &mut dec).unwrap();
    assert_eq!(
        points(&model, body),
        [[1, 0, 0], [3, 1, 0]].map(|p| p.map(Rational::from_i64))
    );
    let edge = model.edge(model.edges_of_body(body)[0]);
    let Curve::BSpline(carrier) = model.geom().curve(edge.curve) else {
        panic!("expected an exact spline carrier");
    };
    let midpoint = carrier
        .point_at(&Real::from(Rational::parse("5/2").unwrap()), &mut dec)
        .unwrap();
    let actual =
        [&midpoint.x, &midpoint.y, &midpoint.z].map(|p| dec.rational_value(p).unwrap().unwrap());
    assert_eq!(
        actual,
        ["2", "1/4", "0"].map(|p| Rational::parse(p).unwrap())
    );
}
