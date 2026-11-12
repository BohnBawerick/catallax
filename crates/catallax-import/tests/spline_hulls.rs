//! Spline faces kept apart by their control hulls: a ruled loft written by
//! OpenCASCADE, whose four side faces are bilinear B-spline patches, imports
//! as a solid, with every pair of side faces proved to meet only along their
//! shared edge or not at all.

use catallax_import::{Surface, build_file, read_step};
use catallax_numbers::with_operation_budget;
use catallax_topology::Model;

/// The M27 loft: the square of half-width 1 at height 0, ruled to the
/// quadrilateral `(1,0), (0,2), (-2,0), (0,-1)` at height 2. Written by
/// OpenCASCADE 8.0 from `BRepOffsetAPI_ThruSections` (ruled) with its default
/// AP214 STEP writer.
const LOFT: &[u8] = include_bytes!("fixtures/m27-loft-occt.step");

#[test]
fn an_occt_ruled_loft_imports_with_its_spline_faces() {
    let file = read_step(LOFT).unwrap();
    let faces = &file.solids[0].faces;
    let splines = faces
        .iter()
        .filter(|face| matches!(face.surface, Surface::BSpline(_)))
        .count();
    assert_eq!((file.solids.len(), faces.len(), splines), (1, 6, 4));
    let mut model = Model::new();
    let (bodies, work) = with_operation_budget(20_000_000, || build_file(&mut model, &file, None));
    let bodies = bodies.unwrap();
    println!("OCCT loft import: {} units", work.used);
    let [body] = bodies[..] else {
        panic!("one body, not {}", bodies.len());
    };
    assert_eq!(model.check(body), vec![]);
    assert_eq!(
        (
            model.faces_of_body(body).len(),
            model.edges_of_body(body).len()
        ),
        (6, 12)
    );
}
