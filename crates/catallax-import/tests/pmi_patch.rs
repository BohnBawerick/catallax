//! PMI-only edits of a STEP file's own records: no edit gives the same bytes,
//! an edit changes only the item's records and the counts that include it, and
//! an edit that cannot be proven safe is refused before anything is written.

use std::collections::BTreeSet;

use catallax_import::assembly::read_step_assembly;
use catallax_import::{PmiEdit, patch_step_pmi, read_step_pmi};
use catallax_numbers::{Decider, Rational};
use catallax_shapes::{Pose, RVector3, Rotation, Surface};
use catallax_step::{
    StepOptions,
    assembly::{Assembly, Occurrence, Part, Product, ProductContent, write_step_assembly},
    pmi::{
        Camera, Characteristic, Datum, GeometricTolerance, Item, ItemKind, PartNotes, Pmi, Process,
        SavedView, Size, SizeKind, SurfaceTexture,
    },
};
use catallax_topology::{Model, block, cylinder};

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn tolerance(
    characteristic: Characteristic,
    value: &str,
    datums: &[&str],
    face: catallax_topology::FaceId,
) -> ItemKind<catallax_topology::FaceId> {
    ItemKind::Tolerance(GeometricTolerance {
        characteristic,
        value: q(value),
        diameter_zone: false,
        material: None,
        datums: datums.iter().map(|letter| (*letter).to_owned()).collect(),
        features: vec![vec![face]],
    })
}

/// A two-part assembly with colours and PMI: a base block with datum A, a
/// flatness, a parallelism shown in a saved view and a note, and a shaft with
/// a diameter. Then records no PMI reader interprets are added: a document
/// reference on the base product, a property, and a group that names the
/// parallelism tolerance.
fn source() -> String {
    let mut model = Model::new();
    let base = block(&mut model, [q("6"), q("8"), q("3")], &Pose::identity()).unwrap();
    let shaft = cylinder(&mut model, q("1"), q("5"), &Pose::identity()).unwrap();
    let faces = model.faces_of_body(base);
    let base_pmi = Pmi {
        items: vec![
            Item {
                id: "Base datum".into(),
                kind: ItemKind::Datum(Datum {
                    letter: "A".into(),
                    feature: vec![faces[0]],
                }),
            },
            Item {
                id: "Base flatness".into(),
                kind: tolerance(Characteristic::Flatness, "0.05", &[], faces[2]),
            },
            Item {
                id: "Base parallelism".into(),
                kind: tolerance(Characteristic::Parallelism, "0.01", &["A"], faces[1]),
            },
            Item {
                id: "Base profile".into(),
                kind: tolerance(Characteristic::SurfaceProfile, "0.2", &[], faces[3]),
            },
            Item {
                id: "Base finish".into(),
                kind: ItemKind::Texture(SurfaceTexture {
                    parameter: "Ra".into(),
                    value: q("3.2"),
                    process: Process::Any,
                    feature: vec![faces[5]],
                }),
            },
            Item {
                id: "Base note".into(),
                kind: ItemKind::FaceNote {
                    text: "Mask before painting".into(),
                    feature: vec![faces[4]],
                },
            },
        ],
        part: PartNotes::default(),
        views: vec![SavedView {
            name: "Profile".into(),
            camera: Camera {
                eye: [q("3"), q("-50"), q("1.5")],
                direction: [q("0"), q("1"), q("0")],
                up: [q("0"), q("0"), q("1")],
                distance: q("54"),
                window: [q("20"), q("12")],
                perspective: false,
            },
            items: vec!["Base profile".into()],
        }],
    };
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
    let shaft_pmi = Pmi {
        items: vec![Item {
            id: "Shaft diameter".into(),
            kind: ItemKind::Size(Size {
                kind: SizeKind::Diameter,
                features: vec![vec![barrel]],
                limits: None,
                through: false,
                thread: None,
                nominal: None,
            }),
        }],
        part: PartNotes::default(),
        views: Vec::new(),
    };
    let at =
        |x: &str| Pose::from_rational(Rotation::identity(), RVector3::new(q(x), q("0"), q("0")));
    let products = vec![
        Product {
            name: "Fixture".into(),
            content: ProductContent::Assembly(vec![
                Occurrence {
                    name: "Base".into(),
                    product: 1,
                    pose: at("0"),
                },
                Occurrence {
                    name: "Pin".into(),
                    product: 2,
                    pose: at("10"),
                },
            ]),
        },
        Product {
            name: "Base".into(),
            content: ProductContent::Part(Part {
                model: &model,
                body: base,
                colour: Some([q("0.125"), q("0.25"), q("0.5")]),
                pmi: Some(&base_pmi),
            }),
        },
        Product {
            name: "Pin".into(),
            content: ProductContent::Part(Part {
                model: &model,
                body: shaft,
                colour: Some([q("0.75"), q("0.25"), q("0")]),
                pmi: Some(&shaft_pmi),
            }),
        },
    ];
    let mut bytes = Vec::new();
    write_step_assembly(
        &Assembly { root: 0, products },
        &StepOptions::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    let file = String::from_utf8(bytes).unwrap();
    let highest: u64 = file
        .lines()
        .filter_map(|line| line.strip_prefix('#')?.split('=').next()?.parse().ok())
        .max()
        .unwrap();
    let parallelism = file
        .lines()
        .find(|line| line.contains("PARALLELISM_TOLERANCE('Base parallelism'"))
        .and_then(|line| line.split('=').next())
        .unwrap()
        .to_owned();
    let product = file
        .lines()
        .find(|line| line.contains("=PRODUCT_DEFINITION('design'"))
        .and_then(|line| line.split('=').next())
        .unwrap()
        .to_owned();
    let (a, b, c, d) = (highest + 1, highest + 2, highest + 3, highest + 4);
    let opaque = format!(
        "#{a}=DOCUMENT('DRW-7','Inspection drawing','',#{b});\n\
         #{b}=DOCUMENT_TYPE('configuration controlled document version');\n\
         #{c}=APPLIED_DOCUMENT_REFERENCE(#{a},'review',({product}));\n\
         #{d}=APPLIED_GROUP_ASSIGNMENT(#{}, ({parallelism}));\n\
         #{}=GROUP('Supplier review','kept as written');\n",
        highest + 5,
        highest + 5,
    );
    file.replacen("ENDSEC;\nEND-ISO", &format!("{opaque}ENDSEC;\nEND-ISO"), 1)
}

fn lines(file: &str) -> BTreeSet<&str> {
    file.lines().collect()
}

#[test]
fn no_edit_gives_the_same_bytes() {
    let file = source();
    assert_eq!(
        patch_step_pmi(file.as_bytes(), &[]).unwrap(),
        file.as_bytes()
    );
}

#[test]
fn a_tolerance_value_changes_its_own_records_only() {
    let file = source();
    let before = read_step_pmi(file.as_bytes()).unwrap();
    assert_eq!(before.skipped, Vec::new());
    let out = patch_step_pmi(
        file.as_bytes(),
        &[PmiEdit::ToleranceValue {
            item: "Base flatness".into(),
            value: q("0.02"),
        }],
    )
    .unwrap();
    let out = String::from_utf8(out).unwrap();
    if let Ok(directory) = std::env::var("CATALLAX_PMI_OUT") {
        let directory = std::path::Path::new(&directory);
        std::fs::write(directory.join("c4d-source.step"), &file).unwrap();
        std::fs::write(directory.join("c4d-flatness-0.02.step"), &out).unwrap();
    }

    // Exactly the tolerance record changed, its old value went and a new value
    // came under a fresh number.
    let (old, new) = (lines(&file), lines(&out));
    let gone: Vec<&&str> = old.difference(&new).collect();
    let came: Vec<&&str> = new.difference(&old).collect();
    assert_eq!(gone.len(), 2, "{gone:?}");
    assert_eq!(came.len(), 2, "{came:?}");
    assert!(
        gone.iter()
            .any(|line| line.contains("FLATNESS_TOLERANCE('Base flatness'"))
    );
    assert!(
        gone.iter()
            .any(|line| line.contains("LENGTH_MEASURE(0.05)"))
    );
    assert!(
        came.iter()
            .any(|line| line.contains("LENGTH_MEASURE(0.02)"))
    );
    let highest_before: u64 = file
        .lines()
        .filter_map(|line| line.strip_prefix('#')?.split('=').next()?.parse().ok())
        .max()
        .unwrap();
    assert!(
        came.iter()
            .any(|line| line.starts_with(&format!("#{}=", highest_before + 1)))
    );

    // The PMI reads back edited, the rest as before.
    let after = read_step_pmi(out.as_bytes()).unwrap();
    let mut expected = before.pmi.clone();
    for item in &mut expected.items {
        if let (true, ItemKind::Tolerance(tolerance)) = (item.id == "Base flatness", &mut item.kind)
        {
            tolerance.value = q("0.02");
        }
    }
    assert_eq!(after.pmi, expected);
    assert_eq!(after.skipped, before.skipped);

    // The assembly, its names, colours and per-part PMI read as before.
    let mut dec = Decider::default();
    let original = read_step_assembly(file.as_bytes(), &mut dec).unwrap();
    let patched = read_step_assembly(out.as_bytes(), &mut dec).unwrap();
    assert_eq!(original.products.len(), patched.products.len());
    for (a, b) in original.products.iter().zip(&patched.products) {
        assert_eq!(a.name, b.name);
        match (&a.content, &b.content) {
            (ProductContent::Part(a), ProductContent::Part(b)) => {
                assert_eq!(a.colour, b.colour);
                assert_eq!(a.pmi.pmi.items.len(), b.pmi.pmi.items.len());
            }
            (ProductContent::Assembly(a), ProductContent::Assembly(b)) => {
                assert_eq!(a.len(), b.len());
            }
            _ => panic!("{} changed kind", a.name),
        }
    }
    for line in file.lines().filter(|line| {
        [
            "DOCUMENT(",
            "DOCUMENT_TYPE(",
            "APPLIED_DOCUMENT_REFERENCE(",
            "GROUP(",
            "COLOUR_RGB(",
        ]
        .iter()
        .any(|kind| line.contains(kind))
    }) {
        assert!(out.contains(line), "{line}");
    }
}

#[test]
fn an_unsafe_edit_is_refused_before_writing() {
    let file = source();
    for (edit, reason) in [
        (
            PmiEdit::ToleranceValue {
                item: "Base profile".into(),
                value: q("0.1"),
            },
            "(DRAUGHTING_MODEL_ITEM_ASSOCIATION) may still depend on the old value",
        ),
        (
            PmiEdit::ToleranceValue {
                item: "Shaft diameter".into(),
                value: q("2"),
            },
            "only a geometric tolerance's value is patched",
        ),
        (
            PmiEdit::ToleranceValue {
                item: "Base flatness".into(),
                value: q("1/3"),
            },
            "not a decimal",
        ),
        (
            PmiEdit::ToleranceValue {
                item: "No such item".into(),
                value: q("2"),
            },
            "no such item was read",
        ),
    ] {
        let refusal = patch_step_pmi(file.as_bytes(), &[edit]).unwrap_err();
        assert!(
            refusal.detail.contains(reason),
            "{} lacks {reason}",
            refusal.detail
        );
    }
    // A tolerance named by the opaque group keeps its number, so its value
    // can change while the group still names it.
    let out = patch_step_pmi(
        file.as_bytes(),
        &[PmiEdit::ToleranceValue {
            item: "Base parallelism".into(),
            value: q("0.03"),
        }],
    )
    .unwrap();
    let out = String::from_utf8(out).unwrap();
    let group = file
        .lines()
        .find(|line| line.contains("APPLIED_GROUP_ASSIGNMENT("))
        .unwrap();
    assert!(out.contains(group));
}

#[test]
fn tolerance_patch_preserves_literals() {
    let file = source();
    let before = read_step_pmi(file.as_bytes()).unwrap();
    let flatness = &before.sources["Base flatness"];
    let flatness_record = file
        .lines()
        .find(|line| line.starts_with(&format!("{flatness}=")))
        .unwrap();
    let value = flatness_record.split(',').nth(2).unwrap();
    let description = format!("'Audit refers to {value}'");
    let annotated = flatness_record.replacen("''", &description, 1);
    let annotated =
        annotated.replacen(&format!(",{value},"), &format!(",/* {value} */{value},"), 1);
    let annotated_file = file.replacen(flatness_record, &annotated, 1);
    let edit = PmiEdit::ToleranceValue {
        item: "Base flatness".into(),
        value: q("0.02"),
    };
    let patched =
        String::from_utf8(patch_step_pmi(annotated_file.as_bytes(), &[edit]).unwrap()).unwrap();
    assert!(patched.contains(&description));
    assert!(patched.contains(&format!("/* {value} */#")));
}

#[test]
fn tolerance_patch_refuses_zone_properties() {
    let file = source();
    let before = read_step_pmi(file.as_bytes()).unwrap();
    let flatness = &before.sources["Base flatness"];
    let part = file
        .lines()
        .find_map(|line| {
            line.split_once("=PRODUCT_DEFINITION_SHAPE(")
                .map(|(id, _)| id)
        })
        .unwrap();
    let dependent = format!(
        "#900000=PROPERTY_DEFINITION('zone validation','',#900001);\n\
         #900001=TOLERANCE_ZONE('','',{part},.F.,({flatness}),#900002);\n\
         #900002=TOLERANCE_ZONE_FORM('between two parallel planes');\n"
    );
    let with_dependency = file.replacen(
        "ENDSEC;\nEND-ISO",
        &format!("{dependent}ENDSEC;\nEND-ISO"),
        1,
    );
    let Err(refusal) = patch_step_pmi(
        with_dependency.as_bytes(),
        &[PmiEdit::ToleranceValue {
            item: "Base flatness".into(),
            value: q("0.02"),
        }],
    ) else {
        panic!("a property using the tolerance zone must prevent the edit");
    };
    assert!(refusal.detail.contains("#900000"), "{}", refusal.detail);
}
