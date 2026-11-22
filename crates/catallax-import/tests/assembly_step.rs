//! Nested STEP export, rebuild and re-export with shared definitions and PMI.

use std::{num::NonZeroUsize, path::PathBuf};

use catallax_import::assembly::read_step_assembly;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Role, StepId};
use catallax_shapes::{Pose, RVector3, Rotation, Surface};
use catallax_step::{
    StepOptions,
    assembly::{Assembly, Occurrence, Part, Product, ProductContent, write_step_assembly},
    pmi::{
        Characteristic, Datum, GeometricTolerance, Item, ItemKind, PartNotes, Pmi, Size, SizeKind,
    },
};
use catallax_topology::{BodyId, FaceId, Model, block, cylinder};

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn annotation(model: &Model, body: BodyId, note: &str) -> Pmi<FaceId> {
    let faces = model.faces_of_body(body);
    Pmi {
        items: vec![
            Item {
                id: "Mounting datum".into(),
                kind: ItemKind::Datum(Datum {
                    letter: "A".into(),
                    feature: vec![faces[0]],
                }),
            },
            Item {
                id: "Mount parallelism".into(),
                kind: ItemKind::Tolerance(GeometricTolerance {
                    characteristic: Characteristic::Parallelism,
                    value: q("0.01"),
                    diameter_zone: false,
                    material: None,
                    datums: vec!["A".into()],
                    features: vec![vec![faces[1]]],
                }),
            },
        ],
        part: PartNotes {
            notes: vec![note.into()],
            ..PartNotes::default()
        },
    }
}

fn instance(name: &str, product: usize, translation: [&str; 3], rotation: Rotation) -> Occurrence {
    let [x, y, z] = translation.map(q);
    Occurrence {
        name: name.into(),
        product,
        pose: Pose::from_rational(rotation, RVector3::new(x, y, z)),
    }
}

fn export_and_read(threads: usize) -> (Vec<u8>, Vec<u8>, u64) {
    catallax_numbers::set_threads(NonZeroUsize::new(threads));
    let mut model = Model::new();
    let base = block(&mut model, [q("6"), q("8"), q("3")], &Pose::identity()).unwrap();
    let shaft = cylinder(&mut model, q("1"), q("5"), &Pose::identity()).unwrap();
    let base_pmi = annotation(&model, base, "Base note #17");
    let mut shaft_pmi = annotation(&model, shaft, "Shaft note #8");
    let barrel = model
        .faces_of_body(shaft)
        .into_iter()
        .find(|&face| {
            matches!(
                model.geom().surface(model.face(face).surface),
                Surface::Cylinder(_)
            )
        })
        .unwrap();
    shaft_pmi.items.push(Item {
        id: "Shaft diameter".into(),
        kind: ItemKind::Size(Size {
            kind: SizeKind::Diameter,
            features: vec![vec![barrel]],
            limits: None,
            through: false,
            thread: None,
            nominal: None,
        }),
    });
    let quarter_turn = Rotation::from_quaternion(q("1"), q("0"), q("0"), q("1")).unwrap();
    let identity = Rotation::identity();
    let products = vec![
        Product {
            name: "Drive assembly".into(),
            content: ProductContent::Assembly(vec![
                instance("Rotated carrier", 2, ["10", "20", "30"], quarter_turn),
                instance("Fixed base", 1, ["4", "6", "8"], identity.clone()),
            ]),
        },
        Product {
            name: "Base's plate".into(),
            content: ProductContent::Part(Part {
                model: &model,
                body: base,
                colour: Some([q("0.125"), q("0.25"), q("0.5")]),
                pmi: Some(&base_pmi),
            }),
        },
        Product {
            name: "Carrier".into(),
            content: ProductContent::Assembly(vec![
                instance("Left shaft", 3, ["2", "0", "0"], identity.clone()),
                instance("Right shaft", 3, ["0", "3", "0"], identity),
            ]),
        },
        Product {
            name: "Shaft #17".into(),
            content: ProductContent::Part(Part {
                model: &model,
                body: shaft,
                colour: Some([q("0.75"), q("0.25"), q("0")]),
                pmi: Some(&shaft_pmi),
            }),
        },
    ];
    let mut dec = Decider::default();
    let mut first = Vec::new();
    let options = StepOptions::default();
    let report = write_step_assembly(
        &Assembly { root: 0, products },
        &options,
        &mut dec,
        &mut first,
    )
    .unwrap();
    assert_eq!(report.parts.len(), 2);
    assert_eq!(report.placements.len(), 4);
    assert!(
        report.placements.iter().all(|placement| {
            placement.origin_error.is_zero() && placement.basis_error.is_zero()
        })
    );
    assert!(report.marks.is_empty(), "{:?}", report.marks);

    let imported = read_step_assembly(&first, &mut dec).unwrap();
    assert_eq!(imported.products.len(), 4);
    assert_eq!(imported.products[imported.root].name, "Drive assembly");
    let ProductContent::Assembly(root) = &imported.products[imported.root].content else {
        panic!("missing root assembly");
    };
    assert_eq!(
        root.iter()
            .map(|occurrence| occurrence.name.as_str())
            .collect::<Vec<_>>(),
        ["Rotated carrier", "Fixed base"]
    );
    let carrier = &imported.products[root[0].product];
    assert_eq!(carrier.name, "Carrier");
    let ProductContent::Assembly(nested) = &carrier.content else {
        panic!("missing nested assembly");
    };
    assert_eq!(nested.len(), 2);
    assert_eq!(nested[0].product, nested[1].product);
    for (occurrence, (name, translation)) in nested.iter().zip([
        ("Left shaft", ["10", "22", "30"]),
        ("Right shaft", ["7", "20", "30"]),
    ]) {
        assert_eq!(occurrence.name, name);
        let pose = occurrence.pose.then(&root[0].pose);
        let [x, y, z] = translation.map(q);
        assert_eq!(
            pose.translation.as_rational().unwrap(),
            RVector3::new(x, y, z)
        );
    }

    let mut recovered_model = Model::new();
    let mut sources = Vec::new();
    for (index, product) in imported.products.iter().enumerate() {
        let ProductContent::Part(part) = &product.content else {
            continue;
        };
        assert_eq!(part.geometry.solids.len(), 1);
        let (name, colour, pmi) = if index == root[1].product {
            ("Base's plate", ["0.125", "0.25", "0.5"], &base_pmi)
        } else {
            assert_eq!(index, nested[0].product);
            ("Shaft #17", ["0.75", "0.25", "0"], &shaft_pmi)
        };
        assert_eq!(product.name, name);
        assert_eq!(part.colour, Some(colour.map(q)));
        assert!(part.pmi.skipped.is_empty(), "{:?}", part.pmi.skipped);
        assert_eq!(part.pmi.pmi.part, pmi.part);
        let step = StepId::new("import").unwrap();
        let body = catallax_import::build_solid(
            &mut recovered_model,
            &part.geometry.solids[0],
            Some(&step),
        )
        .unwrap();
        let pmi = part
            .pmi
            .pmi
            .try_map(|_, entity| {
                recovered_model
                    .faces_of_body(body)
                    .into_iter()
                    .find(|face| {
                        recovered_model.face_origins(*face).iter().any(|origin| {
                    matches!(&origin.role, Role::Imported { key } if key == entity)
                })
                    })
                    .ok_or("PMI face was not rebuilt")
            })
            .unwrap();
        assert_eq!(
            pmi.items.len(),
            if index == root[1].product { 2 } else { 3 }
        );
        let item = &pmi.items[0];
        assert_eq!(item.id, "Mounting datum");
        assert!(
            matches!(&item.kind, ItemKind::Datum(datum) if datum.letter == "A" && datum.feature.len() == 1)
        );
        let tolerance = pmi
            .items
            .iter()
            .find(|item| item.id == "Mount parallelism")
            .unwrap();
        let ItemKind::Tolerance(tolerance) = &tolerance.kind else {
            panic!("missing parallelism");
        };
        assert_eq!(tolerance.characteristic, Characteristic::Parallelism);
        assert_eq!(tolerance.value, q("0.01"));
        assert_eq!(tolerance.datums, ["A"]);
        assert_ne!(tolerance.features[0], *item.kind.features()[0]);
        if index == nested[0].product {
            let diameter = pmi
                .items
                .iter()
                .find(|item| item.id == "Shaft diameter")
                .unwrap();
            let ItemKind::Size(size) = &diameter.kind else {
                panic!("missing diameter");
            };
            assert_eq!(size.nominal, Some(q("2")));
        }
        sources.push((index, body, pmi));
    }

    let products = imported
        .products
        .into_iter()
        .enumerate()
        .map(|(index, product)| {
            let content = match product.content {
                ProductContent::Assembly(children) => ProductContent::Assembly(children),
                ProductContent::Part(part) => {
                    let (_, body, pmi) = sources.iter().find(|source| source.0 == index).unwrap();
                    ProductContent::Part(Part {
                        model: &recovered_model,
                        body: *body,
                        colour: part.colour,
                        pmi: Some(pmi),
                    })
                }
            };
            Product {
                name: product.name,
                content,
            }
        })
        .collect();
    let mut second = Vec::new();
    write_step_assembly(
        &Assembly {
            root: imported.root,
            products,
        },
        &options,
        &mut dec,
        &mut second,
    )
    .unwrap();
    assert_eq!(
        read_step_assembly(&second, &mut dec)
            .unwrap()
            .products
            .len(),
        4
    );
    let work = dec.stats().work;
    assert!(work < 20_000_000, "{work}");
    (first, second, work)
}

#[test]
fn nested_step_assembly_round_trips_part_placement_and_pmi() {
    let first = export_and_read(1);
    for threads in [2, 4] {
        assert_eq!(first, export_and_read(threads));
    }
    println!("Assembly exchange used {} work units", first.2);
    if let Some(directory) = std::env::var_os("CATALLAX_ASSEMBLY_FIXTURES") {
        let directory = PathBuf::from(directory);
        std::fs::create_dir_all(&directory).unwrap();
        std::fs::write(directory.join("first.step"), first.0).unwrap();
        std::fs::write(directory.join("roundtrip.step"), first.1).unwrap();
    }
}

#[test]
fn a_cyclic_assembly_leaves_the_sink_untouched() {
    let assembly = Assembly {
        root: 0,
        products: vec![Product {
            name: "Cycle".into(),
            content: ProductContent::Assembly(vec![instance(
                "Self",
                0,
                ["0", "0", "0"],
                Rotation::identity(),
            )]),
        }],
    };
    let mut out = b"prior output".to_vec();
    let refusal = write_step_assembly(
        &assembly,
        &StepOptions::default(),
        &mut Decider::default(),
        &mut out,
    )
    .err()
    .unwrap();
    assert_eq!(refusal.code, "step.assembly-invalid");
    assert_eq!(out, b"prior output");
}
