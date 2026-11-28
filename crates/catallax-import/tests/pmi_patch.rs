//! PMI-only edits of a STEP file's own records: no edit gives the same bytes,
//! an edit changes only the item's magnitude records, and
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
    assembly_source(false)
}

fn assembly_source(repeat_items: bool) -> String {
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
                id: "Base position".into(),
                kind: tolerance(Characteristic::Position, "0.2", &[], faces[3]),
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
            items: vec!["Base position".into()],
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
    let mut shaft_pmi = Pmi {
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
    if repeat_items {
        let caps: Vec<_> = model
            .faces_of_body(shaft)
            .into_iter()
            .filter(|face| *face != barrel)
            .collect();
        shaft_pmi.items.extend([
            Item {
                id: "Base datum".into(),
                kind: ItemKind::Datum(Datum {
                    letter: "A".into(),
                    feature: vec![caps[0]],
                }),
            },
            Item {
                id: "Base parallelism".into(),
                kind: tolerance(Characteristic::Parallelism, "0.04", &["A"], caps[1]),
            },
        ]);
    }
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
         #{d}=APPLIED_GROUP_ASSIGNMENT(#{}, ({product}));\n\
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
        patch_step_pmi(file.as_bytes(), None).unwrap(),
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
        Some(&PmiEdit::ToleranceValue {
            source: read_step_pmi(file.as_bytes()).unwrap().sources["Base flatness"].clone(),
            value: q("0.02"),
        }),
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
            "APPLIED_GROUP_ASSIGNMENT(",
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
                source: read_step_pmi(file.as_bytes()).unwrap().sources["Base position"].clone(),
                value: q("0.1"),
            },
            "(DRAUGHTING_MODEL_ITEM_ASSOCIATION) may still depend on the old value",
        ),
        (
            PmiEdit::ToleranceValue {
                source: read_step_pmi(file.as_bytes()).unwrap().sources["Shaft diameter"].clone(),
                value: q("2"),
            },
            "only a geometric tolerance's value is patched",
        ),
        (
            PmiEdit::ToleranceValue {
                source: read_step_pmi(file.as_bytes()).unwrap().sources["Base flatness"].clone(),
                value: q("1/3"),
            },
            "not a decimal",
        ),
        (
            PmiEdit::ToleranceValue {
                source: "#900001".into(),
                value: q("2"),
            },
            "no item on one part was read",
        ),
    ] {
        let refusal = patch_step_pmi(file.as_bytes(), Some(&edit)).unwrap_err();
        assert!(
            refusal.detail.contains(reason),
            "{} lacks {reason}",
            refusal.detail
        );
    }
    let parallelism = &read_step_pmi(file.as_bytes()).unwrap().sources["Base parallelism"];
    let group = file
        .lines()
        .find(|line| line.contains("=APPLIED_GROUP_ASSIGNMENT("))
        .unwrap();
    let (assignment, _) = group.rsplit_once(", (").unwrap();
    let unsafe_group = file.replacen(group, &format!("{assignment}, ({parallelism}));"), 1);
    let refusal = patch_step_pmi(
        unsafe_group.as_bytes(),
        Some(&PmiEdit::ToleranceValue {
            source: parallelism.clone(),
            value: q("0.03"),
        }),
    )
    .unwrap_err();
    assert!(refusal.detail.contains("APPLIED_GROUP_ASSIGNMENT"));
}

#[test]
fn an_edit_refuses_a_callout_linked_to_the_controlled_feature() {
    let file = source();
    let read = read_step_pmi(file.as_bytes()).unwrap();
    let root = &read.sources["Base position"];
    let tolerance = file
        .lines()
        .find(|line| line.starts_with(&format!("{root}=")))
        .unwrap();
    let target = tolerance.split(',').nth(3).unwrap().trim_end_matches(");");
    let association = file
        .lines()
        .find(|line| {
            line.contains("=DRAUGHTING_MODEL_ITEM_ASSOCIATION(")
                && line.contains(&format!(",{root},"))
        })
        .unwrap();
    let attached = file.replacen(
        association,
        &association.replace(&format!(",{root},"), &format!(",{target},")),
        1,
    );
    let Err(refusal) = patch_step_pmi(
        attached.as_bytes(),
        Some(&PmiEdit::ToleranceValue {
            source: root.clone(),
            value: q("0.1"),
        }),
    ) else {
        panic!("patched a tolerance that still has a feature-linked callout");
    };
    assert!(
        refusal.detail.contains("DRAUGHTING_MODEL_ITEM_ASSOCIATION"),
        "{refusal:?}"
    );
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
        source: read_step_pmi(file.as_bytes()).unwrap().sources["Base flatness"].clone(),
        value: q("0.02"),
    };
    let patched =
        String::from_utf8(patch_step_pmi(annotated_file.as_bytes(), Some(&edit)).unwrap()).unwrap();
    assert!(patched.contains(&description));
    assert!(patched.contains(&format!("/* {value} */#")));
}

#[test]
fn tolerance_patch_refuses_retained_magnitude_text() {
    let file = source();
    let original = read_step_pmi(file.as_bytes()).unwrap();
    for (name, carrier, magnitude) in [
        ("Base flatness", "FLATNESS_TOLERANCE", "0.05"),
        ("Base parallelism", "GEOMETRIC_TOLERANCE", "0.01"),
    ] {
        let source = &original.sources[name];
        let plain = file
            .lines()
            .find(|line| line.starts_with(&format!("{source}=")))
            .unwrap();
        let line = if carrier == "GEOMETRIC_TOLERANCE" {
            let (attrs, datums) = plain
                .strip_prefix(&format!("{source}=PARALLELISM_TOLERANCE("))
                .and_then(|rest| rest.strip_suffix(");"))
                .and_then(|rest| rest.rsplit_once(",("))
                .unwrap();
            format!(
                "{source}=(GEOMETRIC_TOLERANCE({attrs})GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE(({datums})PARALLELISM_TOLERANCE());"
            )
        } else {
            plain.to_owned()
        };
        let edit = PmiEdit::ToleranceValue {
            source: source.clone(),
            value: q("0.02"),
        };
        let baseline = file.replacen(plain, &line, 1);
        patch_step_pmi(baseline.as_bytes(), Some(&edit)).unwrap();
        let base = format!("{carrier}('{name}','',");
        for text in [
            format!(
                "{carrier}('{name}','MAX {} mm',",
                magnitude.replace('.', ",")
            ),
            format!(
                "{carrier}('{name}','MAX {} mm',",
                magnitude.replace('.', " , ")
            ),
            format!("{carrier}('{name}','MAX \"{magnitude}\"',"),
            format!("{carrier}('{name}','MAX ''{magnitude}''',"),
            format!("{carrier}('{name} {magnitude} mm','',"),
            format!("{carrier}('{name}','{magnitude} mm maximum',"),
            format!("{carrier}('{name}','MAX{magnitude}mm',"),
            format!(
                "{carrier}('{name}','{}0mm',",
                magnitude.trim_start_matches('0')
            ),
            format!(
                "{carrier}('{name}','{} um',",
                (&q(magnitude) * &q("1000")).to_decimal_exact().unwrap()
            ),
        ] {
            let changed = line.replace(&base, &text);
            assert_ne!(changed, line);
            let changed = baseline.replacen(&line, &changed, 1);
            let refusal = patch_step_pmi(changed.as_bytes(), Some(&edit)).unwrap_err();
            assert!(
                refusal.detail.contains(&format!("{source} (")),
                "{refusal:?}"
            );
        }
    }
    let source = &original.sources["Base flatness"];
    let named = file.replacen(
        "FLATNESS_TOLERANCE('Base flatness','',",
        "FLATNESS_TOLERANCE('','',",
        1,
    );
    let named = named.replacen(
        "ENDSEC;\nEND-ISO",
        &format!("#900000=ID_ATTRIBUTE('MAX0.05mm',{source});\nENDSEC;\nEND-ISO"),
        1,
    );
    let read = read_step_pmi(named.as_bytes()).unwrap();
    assert_eq!(&read.sources["MAX0.05mm"], source);
    let error = patch_step_pmi(
        named.as_bytes(),
        Some(&PmiEdit::ToleranceValue {
            source: source.clone(),
            value: q("0.02"),
        }),
    )
    .unwrap_err();
    assert!(error.detail.contains("#900000 (ID_ATTRIBUTE)"), "{error:?}");
}

#[test]
fn tolerance_patch_checks_attached_length_measures_in_their_units() {
    let file = source();
    let original = read_step_pmi(file.as_bytes()).unwrap();
    let source = &original.sources["Base flatness"];
    let record = |id: &str| {
        file.lines()
            .find(|line| line.starts_with(&format!("{id}=")))
            .unwrap()
    };
    let tolerance = record(source);
    let aspect = tolerance.rsplit(',').next().unwrap().trim_end_matches(");");
    let old_measure = tolerance.split(',').nth(2).unwrap();
    let unit = record(old_measure)
        .rsplit(',')
        .next()
        .unwrap()
        .trim_end_matches(");");
    let representation = file
        .lines()
        .find(|line| line.contains("=ADVANCED_BREP_SHAPE_REPRESENTATION("))
        .unwrap();
    let context = representation
        .rsplit(',')
        .next()
        .unwrap()
        .trim_end_matches(");");
    let edit = PmiEdit::ToleranceValue {
        source: source.clone(),
        value: q("0.02"),
    };
    let context_record = record(context);
    let mut precise_context = context_record.to_owned();
    if let Some(start) = precise_context.find("GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT(") {
        let end = start + precise_context[start..].find("))").unwrap() + 2;
        precise_context.replace_range(start..end, "");
    }
    precise_context = precise_context.replacen(
        "GLOBAL_UNIT_ASSIGNED_CONTEXT",
        "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#900005))GLOBAL_UNIT_ASSIGNED_CONTEXT",
        1,
    );
    let file = file.replacen(context_record, &precise_context, 1);
    for (value, measure_unit) in [("0.05", unit), ("0.00005", "#900004"), ("0.050001", unit)] {
        let changes = format!(
            "#900000=PROPERTY_DEFINITION('pmi validation property','',{aspect});\n             #900001=REPRESENTATION('',(#900002),{context});\n             #900002=MEASURE_REPRESENTATION_ITEM('limit',POSITIVE_LENGTH_MEASURE({value}),{measure_unit});\n             #900003=PROPERTY_DEFINITION_REPRESENTATION(#900000,#900001);\n             #900004=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT($,.METRE.));\n\
             #900005=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(0.00001),{unit},'distance_accuracy_value','');\n"
        );
        let changed = file.replacen("ENDSEC;\nEND-ISO", &format!("{changes}ENDSEC;\nEND-ISO"), 1);
        let refusal = patch_step_pmi(changed.as_bytes(), Some(&edit)).unwrap_err();
        assert!(
            refusal
                .detail
                .contains("#900002 (MEASURE_REPRESENTATION_ITEM)"),
            "{refusal:?}"
        );
        let independent = changed.replace(
            &format!("POSITIVE_LENGTH_MEASURE({value}),{measure_unit}"),
            &format!("POSITIVE_LENGTH_MEASURE(0.02),{measure_unit}"),
        );
        patch_step_pmi(independent.as_bytes(), Some(&edit)).unwrap();
    }
}

#[test]
fn tolerance_patch_refuses_zone_properties() {
    let file = source();
    let before = read_step_pmi(file.as_bytes()).unwrap();
    let flatness = &before.sources["Base flatness"];
    let aspect = file
        .lines()
        .find_map(|line| line.strip_prefix(&format!("{flatness}=FLATNESS_TOLERANCE(")))
        .unwrap()
        .trim_end_matches(");")
        .rsplit(',')
        .next()
        .unwrap();
    let part = file
        .lines()
        .find_map(|line| line.strip_prefix(&format!("{aspect}=SHAPE_ASPECT(")))
        .unwrap()
        .split(',')
        .nth(2)
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
        Some(&PmiEdit::ToleranceValue {
            source: read_step_pmi(file.as_bytes()).unwrap().sources["Base flatness"].clone(),
            value: q("0.02"),
        }),
    ) else {
        panic!("a property using the tolerance zone must prevent the edit");
    };
    assert!(refusal.detail.contains("#900000"), "{}", refusal.detail);
    let zone_without_property = with_dependency.replace(
        "#900000=PROPERTY_DEFINITION('zone validation','',#900001);\n",
        "",
    );
    let stale_form = zone_without_property.replace(
        "TOLERANCE_ZONE_FORM('between two parallel planes')",
        "TOLERANCE_ZONE_FORM('0.05 mm maximum')",
    );
    let refused = patch_step_pmi(
        stale_form.as_bytes(),
        Some(&PmiEdit::ToleranceValue {
            source: flatness.clone(),
            value: q("0.02"),
        }),
    )
    .unwrap_err();
    assert!(refused.detail.contains("#900002"), "{refused:?}");
    let described =
        zone_without_property.replace("TOLERANCE_ZONE('','',", "TOLERANCE_ZONE('','0.05 mm',");
    let Err(refusal) = patch_step_pmi(
        described.as_bytes(),
        Some(&PmiEdit::ToleranceValue {
            source: flatness.clone(),
            value: q("0.02"),
        }),
    ) else {
        panic!("the patch retained the old value in a zone description");
    };
    assert!(refusal.detail.contains("#900001"), "{refusal:?}");
}

#[test]
fn nonempty_patches_refuse_unindexed_exchange_sections() {
    let file = source();
    let edit = PmiEdit::ToleranceValue {
        source: read_step_pmi(file.as_bytes()).unwrap().sources["Base flatness"].clone(),
        value: q("0.02"),
    };
    for (marker, section) in [
        (
            "ENDSEC;\nDATA;",
            "ENDSEC;\nANCHOR;\n<old-value>=#1;\nENDSEC;\nDATA;",
        ),
        (
            "ENDSEC;\nDATA;",
            "ENDSEC;\nREFERENCE;\n#900000=<ref.stp#target>;\nENDSEC;\nDATA;",
        ),
        (
            "END-ISO-10303-21;",
            "END-ISO-10303-21;\nSIGNATURE;\nYWJj\nENDSEC;",
        ),
    ] {
        let outstanding = file.replacen(marker, section, 1);
        assert_eq!(
            patch_step_pmi(outstanding.as_bytes(), None).unwrap(),
            outstanding.as_bytes()
        );
        assert!(
            patch_step_pmi(outstanding.as_bytes(), Some(&edit)).is_err(),
            "{section}"
        );
    }
}

#[test]
fn a_source_id_selects_the_second_parts_same_named_tolerance() {
    let file = assembly_source(true);
    let mut dec = Decider::default();
    let read = |bytes: &[u8], dec: &mut Decider| {
        let assembly = read_step_assembly(bytes, dec).unwrap();
        assembly
            .products
            .into_iter()
            .filter_map(|product| match product.content {
                ProductContent::Part(part) => Some((product.name, part.pmi)),
                ProductContent::Assembly(_) => None,
            })
            .collect::<std::collections::BTreeMap<_, _>>()
    };
    let original = read(file.as_bytes(), &mut dec);
    let base = &original["Base"];
    let pin = &original["Pin"];
    let name = "Base parallelism";
    assert_ne!(pin.sources[name], base.sources[name]);
    assert!(
        base.skipped
            .iter()
            .any(|item| item.reason.contains("assembly saved-view"))
    );
    assert!(
        base.skipped
            .iter()
            .any(|item| item.reason.contains("assembly graphical"))
    );
    assert!(pin.skipped.is_empty(), "{:?}", pin.skipped);
    let changed = patch_step_pmi(
        file.as_bytes(),
        Some(&PmiEdit::ToleranceValue {
            source: pin.sources[name].clone(),
            value: q("0.06"),
        }),
    )
    .unwrap();
    let after = read(&changed, &mut dec);
    assert_eq!(after["Base"], *base);
    let mut expected = pin.clone();
    let item = expected
        .pmi
        .items
        .iter_mut()
        .find(|item| item.id == name)
        .unwrap();
    let ItemKind::Tolerance(tol) = &mut item.kind else {
        panic!("expected a tolerance")
    };
    assert_eq!(tol.value, q("0.04"));
    tol.value = q("0.06");
    assert_eq!(after["Pin"], expected);
}

#[test]
fn a_qualified_magnitude_cannot_be_replaced() {
    let file = source();
    let original = read_step_pmi(file.as_bytes()).unwrap();
    let id = &original.sources["Base flatness"];
    let magnitude = file
        .lines()
        .find(|line| line.starts_with(&format!("{id}=")))
        .unwrap()
        .split(',')
        .nth(2)
        .unwrap();
    let qualified = file.replacen("ENDSEC;\nEND-ISO", &format!(
        "#900000=MEASURE_QUALIFICATION('', '', {magnitude}, (#900001));\n#900001=TYPE_QUALIFIER('maximum');\nENDSEC;\nEND-ISO"
    ), 1);
    let result = patch_step_pmi(
        qualified.as_bytes(),
        Some(&PmiEdit::ToleranceValue {
            source: id.clone(),
            value: q("0.02"),
        }),
    );
    assert!(result.is_err(), "a replacement detached the qualification");
    assert!(result.unwrap_err().detail.contains("MEASURE_QUALIFICATION"));
}

#[test]
fn an_unrelated_assembly_dimension_and_id_comments_survive_a_leaf_edit() {
    let file = assembly_source(true);
    let assembly = read_step_assembly(file.as_bytes(), &mut Decider::default()).unwrap();
    let ProductContent::Part(base) = &assembly
        .products
        .iter()
        .find(|p| p.name == "Base")
        .unwrap()
        .content
    else {
        panic!("expected Base part");
    };
    let source = base.pmi.sources["Base flatness"].clone();
    let refs = |text: &str| -> Vec<String> {
        text.split('#')
            .skip(1)
            .map(|tail| {
                format!(
                    "#{}",
                    tail.chars()
                        .take_while(char::is_ascii_digit)
                        .collect::<String>()
                )
            })
            .collect()
    };
    let named = |kind: &str| {
        file.lines()
            .filter(|line| line.contains(kind))
            .collect::<Vec<_>>()
    };
    let occurrences = named("=NEXT_ASSEMBLY_USAGE_OCCURRENCE(");
    let occurrences: Vec<_> = occurrences.iter().map(|record| refs(record)).collect();
    let root_product = &occurrences[0][1];
    let shape = named("=PRODUCT_DEFINITION_SHAPE(")
        .into_iter()
        .map(refs)
        .find(|ids| &ids[1] == root_product)
        .unwrap()[0]
        .clone();
    let datum_features = named("=DATUM_FEATURE(");
    let base_feature = refs(datum_features[0])[0].clone();
    let shaft_feature = refs(datum_features[1])[0].clone();
    let dimension = format!(
        "#900000=COMPONENT_PATH_SHAPE_ASPECT('','',{shape},.T.,{},{base_feature});
         #900001=COMPONENT_PATH_SHAPE_ASPECT('','',{shape},.T.,{},{shaft_feature});
         #900002=DIMENSIONAL_LOCATION('linear distance','',#900000,#900001);
",
        occurrences[0][0], occurrences[1][0]
    );
    let extended = file.replacen(
        "ENDSEC;\nEND-ISO",
        &format!("{dimension}ENDSEC;\nEND-ISO"),
        1,
    );
    let colour = file
        .lines()
        .find(|line| line.contains("=COLOUR_RGB("))
        .unwrap();
    let commented = colour.replacen('=', "/* supplier colour */=", 1);
    let extended = extended.replacen(colour, &commented, 1);
    let out = patch_step_pmi(
        extended.as_bytes(),
        Some(&PmiEdit::ToleranceValue {
            source,
            value: q("0.02"),
        }),
    )
    .unwrap();
    let out = String::from_utf8(out).unwrap();
    assert!(out.contains(&dimension));
    assert!(out.contains(&commented));
    let read = read_step_pmi(out.as_bytes()).unwrap();
    assert!(read.pmi.items.iter().any(|item| item.id == "Base flatness"
        && matches!(&item.kind, ItemKind::Tolerance(tolerance) if tolerance.value == q("0.02"))));
}

#[test]
fn new_patch_ids_leave_existing_dangling_references_unbound() {
    let file = source();
    let id = read_step_pmi(file.as_bytes()).unwrap().sources["Base flatness"].clone();
    let qualification = "#900000=MEASURE_QUALIFICATION('','',#900001,(#900002));";
    let extended = file.replacen(
        "ENDSEC;\nEND-ISO",
        &format!("{qualification}\nENDSEC;\nEND-ISO"),
        1,
    );
    let out = patch_step_pmi(
        extended.as_bytes(),
        Some(&PmiEdit::ToleranceValue {
            source: id,
            value: q("0.02"),
        }),
    )
    .unwrap();
    let out = String::from_utf8(out).unwrap();
    assert!(out.lines().any(|record| record == qualification));
    assert!(!out.lines().any(|record| record.starts_with("#900001=")));
    assert!(!out.lines().any(|record| record.starts_with("#900002=")));
    let read = read_step_pmi(out.as_bytes()).unwrap();
    assert!(read.pmi.items.iter().any(|item| item.id == "Base flatness"
        && matches!(&item.kind, ItemKind::Tolerance(tolerance) if tolerance.value == q("0.02"))));
}

#[test]
fn a_patch_refuses_a_length_with_a_derived_si_prefix() {
    let file = source();
    let id = read_step_pmi(file.as_bytes()).unwrap().sources["Base flatness"].clone();
    let malformed = file.replace("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT(*,.METRE.)");
    assert_ne!(malformed, file);
    let report = read_step_pmi(malformed.as_bytes()).unwrap();
    assert!(!report.sources.contains_key("Base flatness"));
    assert!(
        patch_step_pmi(
            malformed.as_bytes(),
            Some(&PmiEdit::ToleranceValue {
                source: id,
                value: q("0.02"),
            })
        )
        .is_err()
    );
}

/// The flatness edit spelled out on the source bytes alone: the tolerance's
/// value reference moves to a fresh record numbered one above the file's
/// highest, the old measure goes when nothing else names it, and the new
/// measure goes in just before the data section ends. Every other byte,
/// including record order and line ends, stays.
fn flatness_splice(file: &str, newline: &str) -> String {
    let record_start = |id: &str| {
        file.find(&format!("{newline}{id}="))
            .map(|at| at + newline.len())
            .unwrap()
    };
    let line_end = |start: usize| start + file[start..].find(newline).unwrap() + newline.len();
    let tolerance_at = file.find("FLATNESS_TOLERANCE('Base flatness'").unwrap();
    let tolerance_start = file[..tolerance_at].rfind(newline).unwrap() + newline.len();
    let tolerance_end = line_end(tolerance_start);
    let tolerance = &file[tolerance_start..tolerance_end];
    // FLATNESS_TOLERANCE(name, description, magnitude, toleranced shape)
    let magnitude = tolerance
        .split("FLATNESS_TOLERANCE(")
        .nth(1)
        .unwrap()
        .split(',')
        .nth(2)
        .unwrap()
        .to_owned();
    let measure_start = record_start(&magnitude);
    let measure_end = line_end(measure_start);
    let measure = &file[measure_start..measure_end];
    let unit = measure
        .rsplit_once(',')
        .unwrap()
        .1
        .trim_end_matches(&[')', ';', '\r', '\n'][..]);
    let highest: u64 = file
        .split(newline)
        .filter_map(|line| line.strip_prefix('#')?.split('=').next()?.parse().ok())
        .max()
        .unwrap();
    let fresh = format!("#{}", highest + 1);
    let named_elsewhere = file
        .split(newline)
        .filter(|line| !line.starts_with(&format!("{magnitude}=")))
        .filter(|line| {
            line.match_indices(magnitude.as_str()).any(|(at, _)| {
                !line[at + magnitude.len()..].starts_with(|c: char| c.is_ascii_digit())
            })
        })
        .count();
    assert_eq!(named_elsewhere, 1, "only the tolerance names {magnitude}");
    let rewritten = tolerance.replacen(&format!(",{magnitude},"), &format!(",{fresh},"), 1);
    let end = file.rfind(&format!("ENDSEC;{newline}END-ISO")).unwrap();
    let mut out = String::new();
    let (first, second) = if tolerance_start < measure_start {
        (
            (tolerance_start, tolerance_end),
            (measure_start, measure_end),
        )
    } else {
        (
            (measure_start, measure_end),
            (tolerance_start, tolerance_end),
        )
    };
    let piece = |range: (usize, usize)| {
        if range.0 == tolerance_start {
            rewritten.as_str()
        } else {
            ""
        }
    };
    out.push_str(&file[..first.0]);
    out.push_str(piece(first));
    out.push_str(&file[first.1..second.0]);
    out.push_str(piece(second));
    out.push_str(&file[second.1..end]);
    out.push_str(&format!(
        "{fresh}=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(0.02),{unit});{newline}"
    ));
    out.push_str(&file[end..]);
    out
}

#[test]
fn the_flatness_edit_is_exactly_its_byte_splice() {
    for newline in ["\n", "\r\n"] {
        let file = source().replace('\n', newline);
        let edit = PmiEdit::ToleranceValue {
            source: read_step_pmi(file.as_bytes()).unwrap().sources["Base flatness"].clone(),
            value: q("0.02"),
        };
        let out = String::from_utf8(patch_step_pmi(file.as_bytes(), Some(&edit)).unwrap()).unwrap();
        let expected = flatness_splice(&file, newline);
        assert!(
            out == expected,
            "line end {newline:?}: the patch differs from the splice at byte {}",
            out.bytes()
                .zip(expected.bytes())
                .position(|(a, b)| a != b)
                .unwrap_or(out.len().min(expected.len()))
        );
    }
}

#[test]
fn a_face_usage_names_only_its_own_parts_faces() {
    let file = source();
    let field = |record: &str, index: usize| -> String {
        let body = record.split_once('(').unwrap().1;
        body.split(',')
            .nth(index)
            .unwrap()
            .trim_end_matches([')', ';'])
            .to_owned()
    };
    let line = |id: &str| -> String {
        file.lines()
            .find(|line| line.starts_with(&format!("{id}=")))
            .unwrap()
            .to_owned()
    };
    let usage_of = |aspect: &str| -> String {
        file.lines()
            .find(|line| {
                line.contains("=GEOMETRIC_ITEM_SPECIFIC_USAGE(")
                    && line.contains(&format!(",{aspect},"))
            })
            .unwrap()
            .to_owned()
    };
    let flatness = file
        .lines()
        .find(|line| line.contains("FLATNESS_TOLERANCE('Base flatness'"))
        .unwrap();
    let base_aspect = field(flatness, 3);
    let base_usage = usage_of(&base_aspect);
    let size = file
        .lines()
        .find(|line| line.contains("=DIMENSIONAL_SIZE("))
        .unwrap();
    let pin_aspect = field(size, 0);
    let pin_usage = usage_of(&pin_aspect);
    // GEOMETRIC_ITEM_SPECIFIC_USAGE(name, description, aspect, representation, face)
    let (base_rep, base_face) = (field(&base_usage, 3), field(&base_usage, 4));
    let (pin_rep, pin_face) = (field(&pin_usage, 3), field(&pin_usage, 4));
    assert_ne!(base_rep, pin_rep);
    let _ = line(&pin_face);

    let mut dec = Decider::default();
    let original = read_step_assembly(file.as_bytes(), &mut dec).unwrap();
    let base_items = |assembly: &Assembly<catallax_import::assembly::Part>| {
        assembly
            .products
            .iter()
            .find(|product| product.name == "Base")
            .and_then(|product| match &product.content {
                ProductContent::Part(part) => Some(part.pmi.clone()),
                _ => None,
            })
            .unwrap()
    };
    let before = base_items(&original);
    for (changed, reason, linked) in [
        // The pin's face through the pin's representation.
        (
            base_usage.replace(
                &format!(",{base_rep},{base_face})"),
                &format!(",{pin_rep},{pin_face})"),
            ),
            "which is not of the aspect's part",
            false,
        ),
        (
            base_usage.replace(
                &format!(",{base_rep},{base_face})"),
                &format!(",{pin_rep},{pin_face})"),
            ),
            "which is not of the aspect's part",
            true,
        ),
        // The pin's face through the base's own representation.
        (
            base_usage.replace(&format!(",{base_face})"), &format!(",{pin_face})")),
            "does not hold",
            false,
        ),
    ] {
        let mut broken = file.replacen(&base_usage, &changed, 1);
        if linked {
            broken = broken.replacen("ENDSEC;\nEND-ISO", &format!(
                "#900999=SHAPE_REPRESENTATION_RELATIONSHIP('','',{base_rep},{pin_rep});\nENDSEC;\nEND-ISO"
            ), 1);
        }
        let assembly = read_step_assembly(broken.as_bytes(), &mut dec)
            .unwrap_or_else(|error| panic!("the assembly was refused: {error:?}"));
        let after = base_items(&assembly);
        let names = |pmi: &catallax_import::StepPmi| {
            pmi.pmi
                .items
                .iter()
                .map(|item| item.id.clone())
                .collect::<Vec<_>>()
        };
        let mut expected = names(&before);
        expected.retain(|id| id != "Base flatness");
        assert_eq!(names(&after), expected);
        assert!(
            after
                .skipped
                .iter()
                .any(|item| item.reason.contains(reason)),
            "{:?}",
            after.skipped
        );
    }
}
