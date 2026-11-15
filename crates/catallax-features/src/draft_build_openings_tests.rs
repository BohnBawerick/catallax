use super::*;
use catallax_numbers::with_operation_budget;
use catallax_recipe::Recipe;
use std::cell::Cell;
use std::time::{Duration, Instant};

fn source() -> (crate::Built, Vec<FaceId>, LoopId) {
    let mut recipe: serde_json::Value =
        serde_json::from_str(include_str!("../tests/fixtures/draft-mixed-openings.json")).unwrap();
    recipe["steps"].as_array_mut().unwrap().pop();
    recipe["result"] = serde_json::json!("pierced");
    let built = crate::evaluate(
        &Recipe::parse(&recipe.to_string()).unwrap(),
        &mut Decider::default(),
    )
    .unwrap();
    let mut sides = Vec::new();
    let mut wire = None;
    for face in built.model.faces_of_body(built.result) {
        let f = built.model.face(face);
        if let Surface::Plane(plane) = built.model.geom().surface(f.surface) {
            if plane.normal.z.as_rational().is_some_and(|z| z.is_zero()) {
                sides.push(face);
            } else {
                wire = f.loops.iter().copied().find(|wire| {
                    matches!(built.model.loop_(*wire), Loop::Coedges(coedges) if coedges.len() == 6)
                });
            }
        }
    }
    (built, sides, wire.unwrap())
}

#[test]
fn indexed_walls_preserve_missing_and_ambiguous_owner_refusals() {
    let (built, sides, wire) = source();
    let mut dec = Decider::default();
    let mut index = SideIndex::new(&built.model, sides.iter().copied(), &mut dec).unwrap();
    let cap = index.resolve(&built.model, wire, &mut dec).unwrap();
    assert_eq!(cap.walls.len(), 6);
    let first = cap.edges[0];
    let owner = index.owners.remove(&first).unwrap();
    let missing = index.resolve(&built.model, wire, &mut dec).err().unwrap();
    assert_eq!(missing.detail, "a polygonal opening edge has no side wall");
    index.owners.insert(first, owner);
    // Adjacent side walls really share vertical edges. They must remain
    // ambiguous, even though a cap's horizontal edges have unique owners.
    let side_wire = built.model.face(sides[0]).loops[0];
    let ambiguous = index
        .resolve(&built.model, side_wire, &mut dec)
        .err()
        .unwrap();
    assert_eq!(
        ambiguous.detail,
        "a polygonal opening edge has several side walls"
    );
}

#[test]
fn wall_inventory_and_wire_resolution_obey_work_and_deadline_limits() {
    let (built, sides, wire) = source();
    let visited = Cell::new(0);
    let (answer, work) = with_operation_budget(10, || {
        SideIndex::new(
            &built.model,
            sides
                .iter()
                .copied()
                .inspect(|_| visited.set(visited.get() + 1)),
            &mut Decider::default(),
        )
    });
    assert!(answer.is_err() && work.exhausted);
    assert!(visited.get() < sides.len());
    let index =
        SideIndex::new(&built.model, sides.iter().copied(), &mut Decider::default()).unwrap();
    let (answer, work) = with_operation_budget(2, || {
        index.resolve(&built.model, wire, &mut Decider::default())
    });
    assert!(answer.is_err() && work.exhausted);
    let mut expired = Decider::default();
    expired.set_deadline(Some(Instant::now() - Duration::from_secs(1)));
    assert!(SideIndex::new(&built.model, sides.iter().copied(), &mut expired).is_err());
    assert!(index.resolve(&built.model, wire, &mut expired).is_err());
}
