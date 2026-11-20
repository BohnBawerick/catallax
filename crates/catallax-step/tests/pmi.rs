//! AP242 semantic PMI on the P1 bracket: what is written reads back the
//! same, on the same faces, with sizes measured from the geometry.
//!
//! Set `CATALLAX_PMI_OUT` to a directory to keep the written file, for a
//! check in another STEP reader.

use std::collections::BTreeMap;

use catallax_features::evaluate;
use catallax_import::read_step_pmi;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::Surface;
use catallax_step::pmi::{
    Anchor, At, Characteristic, Datum, Feature, GeometricTolerance, Item, ItemKind, Limits,
    Location, PartNotes, Pmi, Process, Size, SizeKind, SurfaceTexture, nominals,
};
use catallax_step::{StepOptions, write_step_with_pmi};
use catallax_topology::{BodyId, FaceId, Loop, Model};

fn q(text: &str) -> Rational {
    Rational::parse(text).unwrap()
}

struct Bracket {
    model: Model,
    body: BodyId,
}

fn bracket() -> Bracket {
    let recipe = Recipe::parse(include_str!("../../../parts/P1-bracket.json")).unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    Bracket {
        model: built.model,
        body: built.result,
    }
}

impl Bracket {
    fn faces(&self) -> Vec<FaceId> {
        self.model.faces_of_body(self.body)
    }

    /// The planar faces normal to coordinate `axis` at `at`.
    fn planes(&self, axis: usize, at: &str) -> Feature<FaceId> {
        let at = q(at);
        let mut found: Vec<FaceId> = self
            .faces()
            .into_iter()
            .filter(|face| {
                let Surface::Plane(plane) =
                    self.model.geom().surface(self.model.face(*face).surface)
                else {
                    return false;
                };
                let normal = plane.normal.rational_components().unwrap();
                let origin = plane.origin.rational_components().unwrap();
                (0..3).all(|i| (i == axis) != normal[i].is_zero()) && *origin[axis] == at
            })
            .collect();
        found.sort();
        assert!(!found.is_empty(), "no plane on axis {axis} at {at}");
        found
    }

    /// The cylindrical faces along coordinate `axis` with this radius, by
    /// the two other coordinates of their axis.
    fn cylinders(&self, axis: usize, radius: &str) -> Vec<(FaceId, [Rational; 2])> {
        let radius_sq = &q(radius) * &q(radius);
        let mut found = Vec::new();
        for face in self.faces() {
            let Surface::Cylinder(cylinder) =
                self.model.geom().surface(self.model.face(face).surface)
            else {
                continue;
            };
            let direction = cylinder.axis.rational_components().unwrap();
            let origin = cylinder.origin.rational_components().unwrap();
            if (0..3).all(|i| (i == axis) != direction[i].is_zero())
                && cylinder.radius_sq.as_rational() == Some(radius_sq.clone())
            {
                let others: Vec<Rational> = (0..3)
                    .filter(|i| *i != axis)
                    .map(|i| origin[i].clone())
                    .collect();
                found.push((face, [others[0].clone(), others[1].clone()]));
            }
        }
        found.sort_by(|a, b| a.1.cmp(&b.1));
        found
    }

    /// The round feature along `axis` with this radius whose axis passes
    /// through `(u, v)` in the two other coordinates.
    fn round(&self, axis: usize, radius: &str, at: [&str; 2]) -> Feature<FaceId> {
        let at = [q(at[0]), q(at[1])];
        let mut found: Vec<FaceId> = self
            .cylinders(axis, radius)
            .into_iter()
            .filter(|(_, centre)| *centre == at)
            .map(|(face, _)| face)
            .collect();
        found.sort();
        assert!(
            !found.is_empty(),
            "no cylinder of radius {radius} at {at:?}"
        );
        found
    }

    fn corners(&self, face: FaceId) -> usize {
        let mut vertices = std::collections::BTreeSet::new();
        for id in &self.model.face(face).loops {
            if let Loop::Coedges(coedges) = self.model.loop_(*id) {
                for coedge in coedges {
                    let edge = self.model.edge(self.model.coedge(*coedge).edge);
                    vertices.extend(edge.start);
                    vertices.extend(edge.end);
                }
            }
        }
        vertices.len()
    }
}

fn item(id: &str, kind: ItemKind<FaceId>) -> Item<FaceId> {
    Item {
        id: id.to_owned(),
        kind,
    }
}

fn deviations(upper: &str, lower: &str) -> Option<Limits> {
    Some(Limits::Deviations {
        upper: q(upper),
        lower: q(lower),
    })
}

fn at_feature(feature: &Feature<FaceId>) -> Anchor<FaceId> {
    Anchor {
        feature: feature.clone(),
        at: At::Feature,
    }
}

/// The fourteen items of the P1 bracket test as the entities a file holds, and
/// its three part notes, with the nominal each dimension must measure. The slot
/// item is three dimensions, the four holes go through, and the locating
/// dimensions from the datums follow the items.
fn bracket_pmi(part: &Bracket) -> (Pmi<FaceId>, BTreeMap<&'static str, &'static str>) {
    let a = part.planes(2, "0");
    let b = part.planes(1, "0");
    let c = part.planes(0, "0");
    let bore = part.round(2, "4", ["30", "25"]);
    let boss = part.round(2, "8", ["30", "25"]);
    let holes: Vec<Feature<FaceId>> = [["10", "12"], ["10", "32"], ["50", "12"], ["50", "32"]]
        .into_iter()
        .map(|at| part.round(2, "2", at))
        .collect();
    // Two planes lie at z = 15: the boss top, a ring between two circles, and
    // the slot's lower flat, a rectangle with four corners.
    let at_fifteen = part.planes(2, "15");
    let slot_low: Feature<FaceId> = at_fifteen
        .iter()
        .copied()
        .filter(|face| part.corners(*face) == 4)
        .collect();
    let boss_top: Feature<FaceId> = at_fifteen
        .iter()
        .copied()
        .filter(|face| !slot_low.contains(face))
        .collect();
    assert_eq!((slot_low.len(), boss_top.len()), (1, 1));
    let slot_high = part.planes(2, "21");
    let slot_ends = [
        part.round(1, "3", ["23", "18"]),
        part.round(1, "3", ["37", "18"]),
    ];
    let mut slot_flats = vec![slot_low[0], slot_high[0]];
    slot_flats.sort();

    let datum = |letter: &str, feature: &Feature<FaceId>| {
        item(
            &format!("Datum {letter}"),
            ItemKind::Datum(Datum {
                letter: letter.to_owned(),
                feature: feature.clone(),
            }),
        )
    };
    // The locating dimensions of the features under a position tolerance: a basic
    // distance from each datum plane their axes run along.
    let located = |name: &str, datum: &Feature<FaceId>, to: Anchor<FaceId>, limits| {
        item(
            name,
            ItemKind::Location(Location {
                from: at_feature(datum),
                to,
                limits,
                nominal: None,
            }),
        )
    };
    let mut locating = vec![
        located("Bore from B", &b, at_feature(&bore), Some(Limits::Basic)),
        located("Bore from C", &c, at_feature(&bore), Some(Limits::Basic)),
    ];
    for (index, hole) in holes.iter().enumerate() {
        for (letter, datum) in [("B", &b), ("C", &c)] {
            locating.push(located(
                &format!("Hole {} from {letter}", index + 1),
                datum,
                at_feature(hole),
                Some(Limits::Basic),
            ));
        }
    }
    let all = ["A", "B", "C"].map(str::to_owned).to_vec();
    let mut items = vec![
        datum("A", &a),
        datum("B", &b),
        datum("C", &c),
        item(
            "Bore diameter",
            ItemKind::Size(Size {
                kind: SizeKind::Diameter,
                features: vec![bore.clone()],
                limits: Some(Limits::Fit {
                    deviation: "H".into(),
                    grade: "7".into(),
                }),
                nominal: None,
                through: false,
                thread: None,
            }),
        ),
        item(
            "Bore finish",
            ItemKind::Texture(SurfaceTexture {
                parameter: "Ra".into(),
                value: q("1.6"),
                process: Process::Removal,
                feature: bore.clone(),
            }),
        ),
        item(
            "Bore position",
            ItemKind::Tolerance(GeometricTolerance {
                characteristic: Characteristic::Position,
                value: q("0.05"),
                diameter_zone: true,
                material: None,
                datums: all.clone(),
                features: vec![bore.clone()],
            }),
        ),
        item(
            "Boss diameter",
            ItemKind::Size(Size {
                kind: SizeKind::Diameter,
                features: vec![boss.clone()],
                limits: deviations("0.1", "-0.1"),
                nominal: None,
                through: false,
                thread: None,
            }),
        ),
        item(
            "Boss height",
            ItemKind::Location(Location {
                from: at_feature(&a),
                to: at_feature(&boss_top),
                limits: deviations("0.1", "-0.1"),
                nominal: None,
            }),
        ),
        item(
            "Hole diameter",
            ItemKind::Size(Size {
                kind: SizeKind::Diameter,
                features: holes.clone(),
                limits: deviations("0.1", "0"),
                nominal: None,
                through: true,
                thread: None,
            }),
        ),
        item(
            "Hole position",
            ItemKind::Tolerance(GeometricTolerance {
                characteristic: Characteristic::Position,
                value: q("0.2"),
                diameter_zone: true,
                material: None,
                datums: all,
                features: holes,
            }),
        ),
        item(
            "Slot width",
            ItemKind::Location(Location {
                from: at_feature(&slot_low),
                to: at_feature(&slot_high),
                limits: deviations("0.1", "0"),
                nominal: None,
            }),
        ),
        item(
            "Slot length",
            ItemKind::Location(Location {
                from: Anchor {
                    feature: slot_ends[0].clone(),
                    at: At::Outer,
                },
                to: Anchor {
                    feature: slot_ends[1].clone(),
                    at: At::Outer,
                },
                limits: deviations("0.1", "0"),
                nominal: None,
            }),
        ),
        item(
            "Slot end radius",
            ItemKind::Size(Size {
                kind: SizeKind::Radius,
                features: slot_ends.to_vec(),
                limits: None,
                nominal: None,
                through: false,
                thread: None,
            }),
        ),
        item(
            "Slot height",
            ItemKind::Location(Location {
                from: at_feature(&a),
                to: Anchor {
                    feature: slot_flats,
                    at: At::Middle,
                },
                limits: deviations("0.1", "-0.1"),
                nominal: None,
            }),
        ),
        item(
            "Flatness of A",
            ItemKind::Tolerance(GeometricTolerance {
                characteristic: Characteristic::Flatness,
                value: q("0.05"),
                diameter_zone: false,
                material: None,
                datums: Vec::new(),
                features: vec![a],
            }),
        ),
        item(
            "Perpendicularity of B",
            ItemKind::Tolerance(GeometricTolerance {
                characteristic: Characteristic::Perpendicularity,
                value: q("0.1"),
                diameter_zone: false,
                material: None,
                datums: vec!["A".into()],
                features: vec![b],
            }),
        ),
    ];
    items.extend(locating);
    let nominal = BTreeMap::from([
        ("Bore diameter", "8"),
        ("Boss diameter", "16"),
        ("Boss height", "15"),
        ("Hole diameter", "4"),
        ("Slot width", "6"),
        ("Slot length", "20"),
        ("Slot end radius", "3"),
        ("Slot height", "18"),
        ("Bore from B", "25"),
        ("Bore from C", "30"),
        ("Hole 1 from B", "12"),
        ("Hole 1 from C", "10"),
        ("Hole 2 from B", "32"),
        ("Hole 2 from C", "10"),
        ("Hole 3 from B", "12"),
        ("Hole 3 from C", "50"),
        ("Hole 4 from B", "32"),
        ("Hole 4 from C", "50"),
    ]);
    let pmi = Pmi {
        items,
        part: PartNotes {
            general_tolerance: Some("ISO 2768-m".into()),
            material: Some("6061-T6".into()),
            notes: vec!["Break sharp edges".into()],
        },
    };
    (pmi, nominal)
}

/// The model face of each `ADVANCED_FACE` in a written file, by entity name.
fn written_faces(file: &str, report: &catallax_step::StepReport) -> BTreeMap<String, FaceId> {
    let regions: BTreeMap<&str, FaceId> = report
        .face_map
        .iter()
        .map(|region| (region.file_region.as_str(), region.face))
        .collect();
    let mut faces = BTreeMap::new();
    for line in file.lines() {
        let Some((id, rest)) = line.split_once("=ADVANCED_FACE('") else {
            continue;
        };
        let name = rest.split('\'').next().unwrap();
        faces.insert(id.to_owned(), regions[name]);
    }
    faces
}

/// The same information with the faces of every feature in one order.
fn sorted(pmi: Pmi<FaceId>) -> Pmi<FaceId> {
    let order = |feature: &mut Feature<FaceId>| feature.sort();
    let mut pmi = pmi;
    for item in &mut pmi.items {
        match &mut item.kind {
            ItemKind::Datum(datum) => order(&mut datum.feature),
            ItemKind::Size(size) => size.features.iter_mut().for_each(order),
            ItemKind::Location(location) => {
                order(&mut location.from.feature);
                order(&mut location.to.feature);
            }
            ItemKind::Tolerance(tolerance) => tolerance.features.iter_mut().for_each(order),
            ItemKind::Texture(texture) => order(&mut texture.feature),
        }
    }
    pmi
}

fn write(part: &Bracket, pmi: &Pmi<FaceId>) -> (String, catallax_step::StepReport) {
    let mut bytes = Vec::new();
    let answer = write_step_with_pmi(
        &part.model,
        part.body,
        &StepOptions::default(),
        pmi,
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    (String::from_utf8(bytes).unwrap(), answer.value)
}

#[test]
fn the_bracket_pmi_reads_back_on_the_same_faces() {
    let part = bracket();
    let (pmi, nominal) = bracket_pmi(&part);

    let measured = nominals(&part.model, &pmi, &mut Decider::default()).unwrap();
    for (item, value) in pmi.items.iter().zip(&measured) {
        let expected = nominal.get(item.id.as_str()).map(|text| q(text));
        assert_eq!(
            value.as_ref().map(|value| value.as_rational().unwrap()),
            expected,
            "{}",
            item.id
        );
    }

    let (file, report) = write(&part, &pmi);
    if let Ok(directory) = std::env::var("CATALLAX_PMI_OUT") {
        std::fs::write(
            std::path::Path::new(&directory).join("p1-bracket-pmi.step"),
            &file,
        )
        .unwrap();
    }

    let read = read_step_pmi(file.as_bytes()).unwrap();
    assert_eq!(read.skipped, Vec::new());
    let faces = written_faces(&file, &report);
    let back = read
        .pmi
        .try_map(|_, key| faces.get(key).copied().ok_or_else(|| key.clone()))
        .unwrap();

    // What comes back is what went out, with each dimension's measured value.
    let mut expected = pmi.clone();
    for item in &mut expected.items {
        let value = nominal.get(item.id.as_str()).map(|text| q(text));
        match &mut item.kind {
            ItemKind::Size(size) => size.nominal = value,
            ItemKind::Location(location) => location.nominal = value,
            _ => {}
        }
    }
    assert_eq!(sorted(back), sorted(expected));

    // The geometry is the file the plain writer makes, with the faces named.
    let mut plain = Vec::new();
    catallax_step::write_step(
        &part.model,
        &[part.body],
        &StepOptions::default(),
        &mut Decider::default(),
        &mut plain,
    )
    .unwrap();
    let plain = String::from_utf8(plain).unwrap();
    let count = |text: &str, kind: &str| text.matches(kind).count();
    for kind in ["=ADVANCED_FACE(", "=EDGE_CURVE(", "=VERTEX_POINT("] {
        assert_eq!(count(&file, kind), count(&plain, kind), "{kind}");
    }
    assert_eq!(count(&plain, "SHAPE_ASPECT"), 0);
}

#[test]
fn a_changed_bore_keeps_its_items_and_reads_its_new_size() {
    // The bore edit of the bracket test: radius 4 becomes radius 5.
    let source = include_str!("../../../parts/P1-bracket.json");
    let mut recipe: serde_json::Value = serde_json::from_str(source).unwrap();
    for step in recipe["steps"].as_array_mut().unwrap() {
        if step["id"] == "bore" {
            step["radius"] = "5".into();
        }
    }
    let built = evaluate(
        &Recipe::parse(&recipe.to_string()).unwrap(),
        &mut Decider::default(),
    )
    .unwrap();
    let part = Bracket {
        model: built.model,
        body: built.result,
    };
    let bore = part.round(2, "5", ["30", "25"]);
    let datum = |letter: &str, feature: Feature<FaceId>| {
        item(
            &format!("Datum {letter}"),
            ItemKind::Datum(Datum {
                letter: letter.to_owned(),
                feature,
            }),
        )
    };
    // The bore's fit and its position tolerance, with the datums the tolerance needs.
    let pmi = Pmi {
        items: vec![
            datum("A", part.planes(2, "0")),
            datum("B", part.planes(1, "0")),
            datum("C", part.planes(0, "0")),
            item(
                "Bore diameter",
                ItemKind::Size(Size {
                    kind: SizeKind::Diameter,
                    features: vec![bore.clone()],
                    limits: Some(Limits::Fit {
                        deviation: "H".into(),
                        grade: "7".into(),
                    }),
                    nominal: None,
                    through: false,
                    thread: None,
                }),
            ),
            item(
                "Bore position",
                ItemKind::Tolerance(GeometricTolerance {
                    characteristic: Characteristic::Position,
                    value: q("0.05"),
                    diameter_zone: true,
                    material: None,
                    datums: ["A", "B", "C"].map(str::to_owned).to_vec(),
                    features: vec![bore],
                }),
            ),
        ],
        part: PartNotes::default(),
    };
    let (file, _) = write(&part, &pmi);
    let read = read_step_pmi(file.as_bytes()).unwrap();
    assert!(read.skipped.is_empty(), "{:?}", read.skipped);
    let find = |name: &str| {
        read.pmi
            .items
            .iter()
            .find(|item| item.id == name)
            .unwrap_or_else(|| panic!("{name} is in the file"))
    };
    // The size is the new one and still has its fit.
    let ItemKind::Size(size) = &find("Bore diameter").kind else {
        panic!("a size");
    };
    assert_eq!(size.nominal, Some(q("10")));
    assert_eq!(
        size.limits,
        Some(Limits::Fit {
            deviation: "H".into(),
            grade: "7".into(),
        })
    );
    // The position tolerance is still on the bore: the size's own face, with its zone
    // and its three datums.
    let ItemKind::Tolerance(position) = &find("Bore position").kind else {
        panic!("a geometric tolerance");
    };
    assert_eq!(position.features, size.features);
    assert_eq!(position.characteristic, Characteristic::Position);
    assert_eq!(
        (position.value.clone(), position.diameter_zone),
        (q("0.05"), true)
    );
    assert_eq!(position.datums, ["A", "B", "C"]);
}

#[test]
fn information_that_does_not_fit_the_part_is_refused() {
    let part = bracket();
    let a = part.planes(2, "0");
    let refused = |pmi: &Pmi<FaceId>| {
        write_step_with_pmi(
            &part.model,
            part.body,
            &StepOptions::default(),
            pmi,
            &mut Decider::default(),
            &mut Vec::new(),
        )
        .unwrap_err()
        .code
    };
    let one = |kind| Pmi {
        items: vec![item("X", kind)],
        part: PartNotes::default(),
    };

    // A diameter on a plane.
    assert_eq!(
        refused(&one(ItemKind::Size(Size {
            kind: SizeKind::Diameter,
            features: vec![a.clone()],
            limits: None,
            nominal: None,
            through: false,
            thread: None,
        }))),
        "step.pmi-not-measurable"
    );
    // A distance between faces that are not parallel.
    assert_eq!(
        refused(&one(ItemKind::Location(Location {
            from: at_feature(&a),
            to: at_feature(&part.planes(1, "0")),
            limits: None,
            nominal: None,
        }))),
        "step.pmi-not-measurable"
    );
    // A tolerance to a datum nobody defined.
    assert_eq!(
        refused(&one(ItemKind::Tolerance(GeometricTolerance {
            characteristic: Characteristic::Perpendicularity,
            value: q("0.1"),
            diameter_zone: false,
            material: None,
            datums: vec!["A".into()],
            features: vec![a.clone()],
        }))),
        "step.pmi-invalid"
    );
    // Limits the wrong way round.
    assert_eq!(
        refused(&one(ItemKind::Size(Size {
            kind: SizeKind::Diameter,
            features: vec![part.round(2, "4", ["30", "25"])],
            limits: deviations("0", "0.1"),
            nominal: None,
            through: false,
            thread: None,
        }))),
        "step.pmi-invalid"
    );
}

#[test]
fn a_threaded_hole_reads_back_with_its_designation() {
    let part = bracket();
    let holes: Vec<Feature<FaceId>> = [["10", "12"], ["10", "32"], ["50", "12"], ["50", "32"]]
        .into_iter()
        .map(|at| part.round(2, "2", at))
        .collect();
    let pmi = Pmi {
        items: vec![item(
            "Tapped holes",
            ItemKind::Size(Size {
                kind: SizeKind::Diameter,
                features: holes.clone(),
                limits: None,
                nominal: None,
                through: true,
                thread: Some("M5x0.8 - 6H".to_owned()),
            }),
        )],
        part: PartNotes::default(),
    };
    let (file, _) = write(&part, &pmi);
    assert!(
        file.contains("'thread'") && file.contains("'M5x0.8 - 6H'"),
        "the designation is in the file"
    );
    let read = read_step_pmi(file.as_bytes()).unwrap();
    assert!(read.skipped.is_empty(), "{:?}", read.skipped);
    let [
        Item {
            kind: ItemKind::Size(size),
            ..
        },
    ] = read.pmi.items.as_slice()
    else {
        panic!("one size expected, found {:?}", read.pmi.items);
    };
    assert_eq!(size.thread.as_deref(), Some("M5x0.8 - 6H"));
    assert!(size.through);
    assert_eq!(size.features.len(), 4, "on the four holes");

    // A radius is not threaded.
    let radius = Pmi {
        items: vec![item(
            "Threaded radius",
            ItemKind::Size(Size {
                kind: SizeKind::Radius,
                features: vec![holes[0].clone()],
                limits: None,
                nominal: None,
                through: false,
                thread: Some("M5".to_owned()),
            }),
        )],
        part: PartNotes::default(),
    };
    let refusal = write_step_with_pmi(
        &part.model,
        part.body,
        &StepOptions::default(),
        &radius,
        &mut Decider::default(),
        &mut Vec::new(),
    )
    .unwrap_err();
    assert_eq!(refusal.code, "step.pmi-invalid");
}
