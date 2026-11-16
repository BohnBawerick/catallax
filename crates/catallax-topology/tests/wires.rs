use catallax_numbers::{Decider, MarkLog, Rational};
use catallax_shapes::{Isometry, Pose, Position, RPoint3, RVector3};
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
