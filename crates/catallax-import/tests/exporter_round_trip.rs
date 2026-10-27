//! Files written by Catallax's own exporter read back as the same primitives.

use catallax_import::{Curve, Surface, read_step};
use catallax_numbers::{Decider, Rational};
use catallax_shapes::Pose;
use catallax_step::{StepOptions, write_step};
use catallax_topology::{Model, block, cylinder};

fn export(model: &Model, body: catallax_topology::BodyId) -> Vec<u8> {
    let mut out = Vec::new();
    let mut dec = Decider::default();
    write_step(model, &[body], &StepOptions::default(), &mut dec, &mut out).expect("export");
    out
}

#[test]
fn exported_block_reads_back_with_its_dimensions() {
    let mut m = Model::new();
    let b = block(&mut m, [2, 3, 5].map(Rational::from_i64), &Pose::identity()).unwrap();
    let file = read_step(&export(&m, b)).unwrap();
    let solid = &file.solids[0];
    assert_eq!(solid.faces.len(), 6);
    let coord = |i: usize| {
        solid
            .faces
            .iter()
            .flat_map(|f| f.loops.iter().flat_map(|l| l.edges.iter()))
            .map(|e| e.start.point[i].clone())
            .max()
            .unwrap()
    };
    assert_eq!(
        [coord(0), coord(1), coord(2)].map(|v| v.max(Rational::zero())),
        [2, 3, 5].map(Rational::from_i64)
    );
}

#[test]
fn exported_cylinder_reads_back_as_cylinder_and_circles() {
    let mut m = Model::new();
    let b = cylinder(
        &mut m,
        Rational::from_i64(3),
        Rational::from_i64(7),
        &Pose::identity(),
    )
    .unwrap();
    let file = read_step(&export(&m, b)).unwrap();
    let faces = &file.solids[0].faces;
    assert_eq!(faces.len(), 3);
    assert!(
        faces
            .iter()
            .any(|f| matches!(&f.surface, Surface::Cylinder(_, r) if *r == Rational::from_i64(3)))
    );
    assert!(
        faces
            .iter()
            .flat_map(|f| f.loops.iter().flat_map(|l| l.edges.iter()))
            .any(|e| matches!(e.curve, Curve::Circle(..)))
    );
}
