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

fn count(file: &str, name: &str) -> Vec<String> {
    file.lines()
        .filter(|line| line.contains(&format!("INTEGER_REPRESENTATION_ITEM('{name}'")))
        .map(|line| line.rsplit(',').next().unwrap().to_owned())
        .collect()
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
fn removing_an_item_takes_its_records_and_lowers_its_counts() {
    let file = source();
    let out = patch_step_pmi(
        file.as_bytes(),
        &[PmiEdit::Remove {
            item: "Shaft diameter".into(),
        }],
    )
    .unwrap();
    let out = String::from_utf8(out).unwrap();
    let after = read_step_pmi(out.as_bytes()).unwrap();
    assert!(
        after
            .pmi
            .items
            .iter()
            .all(|item| item.id != "Shaft diameter")
    );
    assert_eq!(
        after.pmi.items.len(),
        read_step_pmi(file.as_bytes()).unwrap().pmi.items.len() - 1
    );
    assert!(!out.contains("DIMENSIONAL_SIZE("));
    // The pin's own counts drop to zero; the base keeps its own.
    let mut before = count(&file, "number of dimensional sizes");
    before.sort();
    let mut now = count(&out, "number of dimensional sizes");
    now.sort();
    assert_eq!(before, ["0.);", "1.);"]);
    assert_eq!(now, ["0.);", "0.);"]);
    let mut elements = count(&out, "number of semantic pmi elements");
    elements.sort();
    assert_eq!(elements, ["0.);", "5.);"]);
    // No record of the file refers to a number that is gone.
    let numbers: BTreeSet<String> = out
        .lines()
        .filter_map(|line| Some(line.strip_prefix('#')?.split('=').next()?.to_owned()))
        .collect();
    for line in out.lines().filter(|line| line.starts_with('#')) {
        let body = line.split_once('=').unwrap().1;
        for reference in body.split('#').skip(1) {
            let digits: String = reference.chars().take_while(char::is_ascii_digit).collect();
            assert!(numbers.contains(&digits), "{line} names #{digits}");
        }
    }
}

#[test]
fn an_unsafe_edit_is_refused_before_writing() {
    let file = source();
    for (edit, reason) in [
        (
            PmiEdit::Remove {
                item: "Base parallelism".into(),
            },
            "(APPLIED_GROUP_ASSIGNMENT) still refers to its records",
        ),
        (
            PmiEdit::Remove {
                item: "Base datum".into(),
            },
            "Base parallelism refers to datum A",
        ),
        (
            PmiEdit::Remove {
                item: "Base profile".into(),
            },
            "presents it in a callout that would stay",
        ),
        (
            PmiEdit::ToleranceValue {
                item: "Base profile".into(),
                value: q("0.1"),
            },
            "(DRAUGHTING_MODEL_ITEM_ASSOCIATION) would still show the old value",
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
            PmiEdit::Remove {
                item: "Base finish".into(),
            },
            "presented an edited item and would stay as graphics",
        ),
        (
            PmiEdit::Remove {
                item: "No such item".into(),
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
