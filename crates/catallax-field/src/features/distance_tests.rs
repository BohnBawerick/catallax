use super::*;
use crate::{Instruction, Plan};
use catallax_recipe::{Operation, Recipe};
use catallax_shapes::{EllipticCylinder, distance::Options};
use serde_json::json;

#[test]
fn native_elliptic_cylinder_distances_use_the_metric_and_closed_sweep_trim() {
    // Recipe emission of oblique circular sweeps is a separate capability.
    // Assemble its exact native definition so this query contract does not
    // depend on the geometry builder accepting that operation.
    let recipe = Recipe::parse(&json!({
        "format":"catallax-recipe/3", "name":"ellipse-distance",
        "steps":[
            {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[
                {"kind":"circle","id":"wall","centre":["0","0"],"radius":"1"}
            ]]},
            {"id":"stock","op":"extrude","profile":"profile","distance":"5"}
        ], "result":"stock"
    }).to_string()).unwrap();
    let Operation::Sketch { plane, loops } = &recipe.steps[0].op else {
        unreachable!()
    };
    let mut dec = Decider::default();
    let profile = Profile::new(plane, loops, &mut dec).unwrap();
    let step = &recipe.steps[1].id;
    let q = Rational::from_i64;
    let point = |x, y, z| RPoint3::new(q(x), q(y), q(z)).to_real();
    let vector = |x, y, z| RVector3::new(q(x), q(y), q(z));
    let mut faces = vec![
        face(
            Surface::EllipticCylinder(EllipticCylinder {
                origin: point(0, 0, 0),
                axis: vector(3, 0, 4).into(),
                major_dir: vector(0, 1, 0).into(),
                major_sq: q(1).into(),
                minor_sq: Rational::from_ratio(16, 25).unwrap().into(),
            }),
            step,
            side(profile.loops[0][0].ids.first().unwrap()),
            &mut dec,
        )
        .unwrap(),
        face(
            Surface::Plane(Plane {
                origin: point(0, 0, 0),
                normal: vector(0, 0, 1).into(),
                ref_dir: vector(1, 0, 0).into(),
            }),
            step,
            Role::Start,
            &mut dec,
        )
        .unwrap(),
        face(
            Surface::Plane(Plane {
                origin: point(3, 0, 4),
                normal: vector(0, 0, 1).into(),
                ref_dir: vector(1, 0, 0).into(),
            }),
            step,
            Role::End,
            &mut dec,
        )
        .unwrap(),
    ];
    let source = Feature::new(
        &profile,
        Sweep::Extrude {
            vector: vector(3, 0, 4),
            through: false,
        },
        &faces,
        q(0),
        Vec::new(),
    );
    let mut support = Support::new(&faces[0].surface);
    let reflection =
        Isometry::reflection(&RPoint3::new(q(3), q(0), q(1)), &vector(2, 3, 0)).unwrap();
    support.transform(&reflection);
    let target = reflection.apply_point(
        &RPoint3::new(
            Rational::parse("29/10").unwrap(),
            q(0),
            Rational::parse("17/10").unwrap(),
        )
        .to_real(),
    );
    let enclosure =
        [&target.x, &target.y, &target.z].map(|x| geometry::enclose(x, 128, &mut dec).unwrap());
    let residual = support.proven(&enclosure, &mut dec).unwrap();
    let half = Rational::from_ratio(1, 2).unwrap();
    assert!(residual.lo() <= half && residual.hi() >= half);
    assert!(residual.lo() > Rational::parse("0.499").unwrap());
    assert!(residual.hi() < Rational::parse("0.501").unwrap());
    for face in &mut faces {
        face.leaf = 0;
        face.in_result = true;
    }
    let plan = Plan {
        instructions: vec![Instruction::Leaf(0)],
        leaves: vec![geometry::Leaf::Feature(source)],
        faces,
        root: 0,
        leaf_registers: vec![Some(0)],
        definition_roles: Default::default(),
    };
    let id = plan.select(&plan.faces[0].origin).unwrap();
    let prepared = plan.prepare_distances(&mut dec).unwrap();
    let options = Options {
        width: Rational::from_ratio(1, 1000).unwrap(),
        max_cells: 4096,
    };
    for (xyz, expected) in [
        (["29/10", "0", "17/10"], "1/2"),
        (["23/5", "0", "24/5"], "1"),
    ] {
        let point = RPoint3::new(
            Rational::parse(xyz[0]).unwrap(),
            Rational::parse(xyz[1]).unwrap(),
            Rational::parse(xyz[2]).unwrap(),
        )
        .to_real();
        let answer = prepared
            .point_to_face(id, &point, &options, &mut dec)
            .unwrap();
        let expected = Rational::parse(expected).unwrap();
        let upper = answer.upper.as_ref().unwrap();
        assert!(answer.lower <= expected && upper >= &expected, "{answer:?}");
        assert!(
            answer.converged && upper - &answer.lower <= options.width,
            "{answer:?}"
        );
    }
}
