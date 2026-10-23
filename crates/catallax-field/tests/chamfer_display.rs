//! Real default native-display limits apply to every elementary bevel variant.
use catallax_field::{
    Plan,
    display::{Options, mesh_plan},
};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::{Recipe, Role};
use serde_json::{Value, json};
use std::collections::BTreeSet;

#[test]
fn default_p2_and_fractional_edit_publish_all_six_native_faces() {
    for distance in ["1", "3/4"] {
        let mut value: Value =
            serde_json::from_str(include_str!("../../../parts/P2-collar.json")).unwrap();
        value["steps"][2]["distance"] = json!(distance);
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let options = Options::default();
        assert_eq!(options.tolerance, Rational::from_ratio(1, 100).unwrap());
        assert_eq!(options.max_cells, 2_000_000);
        assert_eq!(options.max_vertices, 1_000_000);
        let (result, work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::default();
            let plan = Plan::compile_checked(&recipe, &mut dec)?;
            let mesh = mesh_plan(&plan, &options, &mut dec)?;
            Ok::<_, catallax_numbers::Refusal>((plan, mesh))
        });
        assert!(!work.exhausted, "{distance}: {work:?}");
        let (plan, mesh) = result.unwrap();
        assert!(!mesh.positions.is_empty() && !mesh.triangles.is_empty());
        assert!(mesh.max_deviation <= options.tolerance);
        assert!(mesh.cells <= options.max_cells && mesh.positions.len() <= options.max_vertices);
        let faces: BTreeSet<_> = mesh
            .triangles
            .iter()
            .map(|t| plan.faces()[t.face.index()].origin.clone())
            .collect();
        assert_eq!(faces.len(), 6, "{distance}: {faces:?}");
        assert_eq!(
            faces
                .iter()
                .filter(|f| matches!(f.role, Role::Bevel { .. }))
                .count(),
            2
        );
        let mut obj = Vec::new();
        mesh.write_obj(&mut obj).unwrap();
        let obj = String::from_utf8(obj).unwrap();
        assert_eq!(
            obj.lines().filter(|l| l.starts_with("v ")).count(),
            mesh.positions.len()
        );
        assert_eq!(
            obj.lines().filter(|l| l.starts_with("f ")).count(),
            mesh.triangles.len()
        );
        eprintln!(
            "distance {distance}: {} cells, {} vertices, {} triangles, {} work",
            mesh.cells,
            mesh.positions.len(),
            mesh.triangles.len(),
            work.used
        );
    }
}
