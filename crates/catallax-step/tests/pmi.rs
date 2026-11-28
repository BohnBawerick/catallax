//! AP242 semantic PMI on the P1 bracket: what is written reads back the
//! same, on the same faces, with sizes measured from the geometry.
//!
//! Set `CATALLAX_PMI_OUT` to a directory to keep the written file, for a
//! check in another STEP reader.

mod pmi_strokes;

use std::collections::BTreeMap;

use catallax_features::evaluate;
use catallax_import::read_step_pmi;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{HolePart, Recipe, Role};
use catallax_shapes::{Pose, Surface};
use catallax_step::pmi::{
    Anchor, At, Camera, ChamferNominal, Characteristic, Datum, DatumRef, DatumTarget, Feature,
    GeometricTolerance, HoleHead, HoleSize, Item, ItemKind, Limits, Location, MaterialCondition,
    PartNotes, Pmi, PositionControl, PositionFrame, PositionSegment, Process, ProjectedZone,
    SavedView, Size, SizeKind, SurfaceTexture, TargetArea, TargetKind, nominals,
};
use catallax_step::{StepOptions, write_step_with_pmi};
use catallax_topology::{BodyId, FaceId, Loop, Model, cylinder};

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
            standard: None,
            general_tolerance: Some("ISO 2768-m".into()),
            material: Some("6061-T6".into()),
            notes: vec!["Break sharp edges".into()],
        },
        views: Vec::new(),
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
            ItemKind::Width { flats, .. } | ItemKind::BasicAngle { planes: flats, .. } => {
                flats.iter_mut().for_each(order);
            }
            ItemKind::Chamfer { bevel, sides, .. } | ItemKind::EdgeBreak { bevel, sides, .. } => {
                order(bevel);
                sides.iter_mut().for_each(order);
            }
            ItemKind::Tolerance(tolerance) => tolerance.features.iter_mut().for_each(order),
            ItemKind::Texture(texture) => order(&mut texture.feature),
            ItemKind::Position(control) => control.features.iter_mut().for_each(order),
            ItemKind::FaceNote { feature, .. } | ItemKind::Knurl { feature, .. } => order(feature),
            ItemKind::DatumTarget(_) => {}
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

/// The entities of a written file by number: type and parameter text.
fn entities(file: &str) -> BTreeMap<String, (String, String)> {
    file.lines()
        .filter_map(|line| {
            let (id, rest) = line.split_once('=')?;
            let rest = rest.strip_suffix(");")?;
            let (kind, fields) = rest.split_once('(')?;
            Some((id.to_owned(), (kind.to_owned(), fields.to_owned())))
        })
        .collect()
}

/// The entity references in a parameter text, in order.
fn references(fields: &str) -> Vec<String> {
    let mut found = Vec::new();
    let mut rest = fields;
    while let Some(at) = rest.find('#') {
        let digits: String = rest[at + 1..]
            .chars()
            .take_while(char::is_ascii_digit)
            .collect();
        found.push(format!("#{digits}"));
        rest = &rest[at + 1 + digits.len()..];
    }
    found
}

/// The numbers of a parameter text such as `'',(30.,25.,7.5)`.
fn numbers(fields: &str) -> Vec<f64> {
    fields
        .split(|c: char| !(c.is_ascii_digit() || c == '.' || c == '-' || c == 'E'))
        .filter(|text| text.chars().any(|c| c.is_ascii_digit()))
        .map(|text| text.parse().unwrap())
        .collect()
}

/// The `pmi validation property` items of each entity they are on.
fn validation_items(file: &str) -> BTreeMap<String, Vec<(String, String)>> {
    let all = entities(file);
    let mut found = BTreeMap::<String, Vec<(String, String)>>::new();
    for (kind, fields) in all.values() {
        if kind != "PROPERTY_DEFINITION_REPRESENTATION" {
            continue;
        }
        let linked = references(fields);
        let [property, representation] = &linked[..] else {
            continue;
        };
        let (_, property) = &all[property];
        if !property.starts_with("'pmi validation property',''") {
            continue;
        }
        let of = references(property).pop().unwrap();
        // The representation's items, then its context.
        let mut items = references(&all[representation].1);
        items.pop();
        for item in items {
            found
                .entry(of.clone())
                .or_default()
                .push(all[&item].clone());
        }
    }
    found
}

#[test]
fn the_bracket_validation_counts_are_step_reals() {
    let part = bracket();
    let (pmi, _) = bracket_pmi(&part);
    let (file, _) = write(&part, &pmi);
    let all = entities(&file);
    let properties = validation_items(&file);

    // One property on the part counts every type of item, and the callout.
    let shape = all
        .iter()
        .find(|(_, (kind, _))| kind == "PRODUCT_DEFINITION_SHAPE")
        .map(|(id, _)| id.clone())
        .unwrap();
    let expected = [
        ("number of semantic pmi elements", "26."),
        ("number of dimensional locations", "14."),
        ("number of dimensional sizes", "4."),
        ("number of geometric tolerances", "4."),
        ("number of composite tolerances", "0."),
        ("number of datums", "3."),
        ("number of datum targets", "0."),
        ("number of annotations", "1."),
    ]
    .map(|(name, value)| {
        (
            "INTEGER_REPRESENTATION_ITEM".to_owned(),
            format!("'{name}',{value}"),
        )
    });
    assert_eq!(properties[&shape], expected);
}

#[test]
fn the_bracket_file_states_its_equivalent_unicode_strings() {
    let part = bracket();
    let (pmi, _) = bracket_pmi(&part);
    let (file, _) = write(&part, &pmi);
    let all = entities(&file);
    let properties = validation_items(&file);

    // Each datum and distance states its text; edition 1 has no property on a
    // size or a tolerance.
    let strings: Vec<(String, String)> = properties
        .iter()
        .filter(|(of, _)| all[*of].0 != "PRODUCT_DEFINITION_SHAPE")
        .filter_map(|(of, items)| {
            let [(kind, fields)] = &items[..] else {
                return None;
            };
            assert_eq!(kind, "DESCRIPTIVE_REPRESENTATION_ITEM");
            let text = fields.strip_prefix("'equivalent unicode string',")?;
            Some((all[of].0.clone(), text.to_owned()))
        })
        .collect();
    assert_eq!(strings.len(), 3 + 14);
    for (kind, text) in [
        ("DATUM", r"'DTM\\wA'"),
        ("DATUM", r"'DTM\\wC'"),
        ("DIMENSIONAL_LOCATION", r"'DIM\\w15\X2\00B1\X0\0.1'"),
        ("DIMENSIONAL_LOCATION", r"'DIM\\w6\\u+0.1\\n0\\u'"),
        ("DIMENSIONAL_LOCATION", r"'DIM\\w20\\u+0.1\\n0\\u'"),
        ("DIMENSIONAL_LOCATION", r"'DIM\\w[25]'"),
        ("DIMENSIONAL_LOCATION", r"'DIM\\w[50]'"),
    ] {
        assert!(
            strings.contains(&(kind.to_owned(), text.to_owned())),
            "{kind} {text} in {strings:?}"
        );
    }
    assert!(
        !all.values()
            .any(|(kind, fields)| kind == "PROPERTY_DEFINITION"
                && fields.starts_with("'pmi validation property'")
                && references(fields)
                    .iter()
                    .any(|of| all[of].0 == "DIMENSIONAL_SIZE")),
        "no property on a size"
    );
}

#[test]
fn the_bore_finish_is_drawn_on_the_bore() {
    let part = bracket();
    let (pmi, _) = bracket_pmi(&part);
    let (file, report) = write(&part, &pmi);
    let all = entities(&file);
    let faces = written_faces(&file, &report);
    let bore = part.round(2, "4", ["30", "25"]);

    // The callout hangs on an aspect of its own whose only use names the bore.
    let associations: Vec<&String> = all
        .values()
        .filter(|(kind, _)| kind == "DRAUGHTING_MODEL_ITEM_ASSOCIATION")
        .map(|(_, fields)| fields)
        .collect();
    assert_eq!(associations.len(), 2, "{associations:?}");
    let association = associations
        .iter()
        .find(|fields| fields.starts_with("'PMI presentation to geometry link',''"))
        .unwrap();
    let semantic = associations
        .iter()
        .find(|fields| fields.starts_with("'PMI representation to presentation link',''"))
        .unwrap();
    let semantic_source = read_step_pmi(file.as_bytes()).unwrap().sources["Bore finish"].clone();
    assert_eq!(references(semantic)[0], semantic_source);
    let [aspect, model, callout] = &references(association)[..] else {
        panic!("{association}");
    };
    assert_eq!(all[aspect].0, "SHAPE_ASPECT");
    assert!(all[aspect].1.starts_with("'Bore finish',''"));
    assert_eq!(all[model].0, "DRAUGHTING_MODEL");
    assert_eq!(all[callout].0, "DRAUGHTING_CALLOUT");
    let uses: Vec<&String> = all
        .values()
        .filter(|(kind, fields)| {
            kind.ends_with("_USAGE") && references(fields).first() == Some(aspect)
        })
        .map(|(_, fields)| fields)
        .collect();
    let [usage] = uses[..] else {
        panic!("{uses:?}");
    };
    let mut on: Vec<FaceId> = references(usage)[2..]
        .iter()
        .map(|face| faces[face])
        .collect();
    on.sort();
    assert_eq!(on, bore);

    let link = line_with(
        &file,
        &[
            "DRAUGHTING_MODEL_ITEM_ASSOCIATION",
            "'PMI representation to presentation link'",
        ],
    );
    let legacy = mutate(&file, &format!("{link}={}\n", records(&file)[&link]), "");
    let legacy = mutate(
        &legacy,
        "'PMI presentation to geometry link'",
        "'Bore finish'",
    );
    let read = read_step_pmi(legacy.as_bytes()).unwrap();
    assert_eq!(
        read.pmi
            .items
            .iter()
            .filter(|item| item.id == "Bore finish")
            .count(),
        1
    );
    assert!(read.graphics.is_empty(), "{:?}", read.graphics);
    assert!(read.skipped.is_empty(), "{:?}", read.skipped);

    // The callout is in the model's annotation plane, which faces the front.
    let (plane, (_, fields)) = all
        .iter()
        .find(|(_, (kind, _))| kind == "ANNOTATION_PLANE")
        .unwrap();
    assert_eq!(references(fields).last(), Some(callout));
    assert_eq!(&references(&all[model].1)[0], plane);

    // Its lines: the leader starts on the bore, every point is in one plane
    // parallel to the annotation plane, and the stated length and centre are
    // those of the lines.
    let occurrence = &references(&all[callout].1)[0];
    let set = references(&all[occurrence].1).pop().unwrap();
    assert!(all[&set].1.starts_with("'surface roughness'"));
    let lines: Vec<Vec<[f64; 3]>> = references(&all[&set].1)
        .iter()
        .map(|line| {
            references(&all[line].1)
                .iter()
                .map(|point| {
                    let xyz = numbers(&all[point].1);
                    [xyz[0], xyz[1], xyz[2]]
                })
                .collect()
        })
        .collect();
    let tip = lines[0][0];
    let off_axis = ((tip[0] - 30.0).powi(2) + (tip[1] - 25.0).powi(2)).sqrt();
    assert!((off_axis - 4.0).abs() < 1e-6, "{tip:?}");
    assert!((0.0..=20.0).contains(&tip[2]), "{tip:?}");
    assert!(
        lines
            .iter()
            .flatten()
            .all(|point| (point[1] - tip[1]).abs() < 1e-9)
    );
    let mut length = 0.0;
    let mut moment = [0.0; 3];
    for line in &lines {
        for pair in line.windows(2) {
            let piece = (0..3)
                .map(|i| (pair[1][i] - pair[0][i]).powi(2))
                .sum::<f64>()
                .sqrt();
            length += piece;
            for i in 0..3 {
                moment[i] += piece * (pair[0][i] + pair[1][i]) / 2.0;
            }
        }
    }
    let properties = validation_items(&file);
    let within = all
        .iter()
        .find(|(_, (kind, fields))| {
            kind == "CHARACTERIZED_ITEM_WITHIN_REPRESENTATION"
                && references(fields) == [callout.clone(), model.clone()]
        })
        .map(|(id, _)| id)
        .unwrap();
    let items = &properties[within];
    assert_eq!(
        items[0],
        (
            "DESCRIPTIVE_REPRESENTATION_ITEM".to_owned(),
            r"'equivalent unicode string','SRFMRRN\\w\\w\\w\\wRa 1.6'".to_owned()
        )
    );
    assert_eq!(items[1].0, "MEASURE_REPRESENTATION_ITEM");
    assert!(
        items[1]
            .1
            .starts_with("'polyline curve length',POSITIVE_LENGTH_MEASURE(")
    );
    let stated = numbers(&items[1].1)[0];
    assert!((stated - length).abs() < 1e-5, "{stated} {length}");
    assert_eq!(items[2].0, "CARTESIAN_POINT");
    assert!(items[2].1.starts_with("'polyline centre point'"));
    let centre = numbers(&items[2].1);
    for i in 0..3 {
        assert!((centre[i] - moment[i] / length).abs() < 1e-5, "{centre:?}");
    }

    // Catallax still reads the bracket's fourteen items and nothing else.
    let read = read_step_pmi(file.as_bytes()).unwrap();
    assert_eq!(read.skipped, Vec::new());
    assert_eq!(read.pmi.items.len(), pmi.items.len());
}

#[test]
fn each_process_condition_draws_its_symbol() {
    let part = bracket();
    let a = part.planes(2, "0");
    let strings = [Process::Any, Process::Removal, Process::NoRemoval].map(|process| {
        let pmi = Pmi {
            items: vec![item(
                "Base finish",
                ItemKind::Texture(SurfaceTexture {
                    parameter: "Rz".into(),
                    value: q("12.5"),
                    process,
                    feature: a.clone(),
                }),
            )],
            part: PartNotes::default(),
            views: Vec::new(),
        };
        let (file, _) = write(&part, &pmi);
        let all = entities(&file);
        let polylines = all.values().filter(|(kind, _)| kind == "POLYLINE").count();
        // The leader ends on a corner of the base.
        let tip = all
            .values()
            .filter(|(kind, _)| kind == "POLYLINE")
            .map(|(_, fields)| numbers(&all[&references(fields)[0]].1))
            .find(|xyz| xyz[2].abs() < 1e-9)
            .expect("a leader on the base");
        assert!(file.lines().any(|line| line.contains("=VERTEX_POINT(")
            && numbers(&all[&references(line.split_once('=').unwrap().1)[0]].1) == tip));
        let string = file
            .lines()
            .find_map(|line| line.split_once("'equivalent unicode string','SRF"))
            .map(|(_, rest)| rest.split('\'').next().unwrap().to_owned())
            .unwrap();
        (string, polylines)
    });
    assert_eq!(strings[0].0, r"APAN\\w\\w\\w\\wRz 12.5");
    assert_eq!(strings[1].0, r"MRRN\\w\\w\\w\\wRz 12.5");
    assert_eq!(strings[2].0, r"NMRN\\w\\w\\w\\wRz 12.5");
    // Removal adds a bar, no removal a circle.
    assert_eq!(strings[1].1, strings[0].1 + 1);
    assert_eq!(strings[2].1, strings[0].1 + 1);
}

#[test]
fn a_conflicting_datum_does_not_reserve_its_letter() {
    let mut model = Model::new();
    let body = cylinder(&mut model, q("4"), q("5"), &Pose::identity()).unwrap();
    let part = Bracket { model, body };
    let bottom = part.planes(2, "0");
    let top = part.planes(2, "5");
    let barrel = part.round(2, "4", ["0", "0"]);
    let datum = |letter: &str, feature| {
        item(
            &format!("Datum {letter}"),
            ItemKind::Datum(Datum {
                letter: letter.into(),
                feature,
            }),
        )
    };
    let pmi = Pmi {
        items: vec![
            datum("A", bottom.clone()),
            datum("B", top),
            item(
                "Diameter",
                ItemKind::Size(Size {
                    kind: SizeKind::Diameter,
                    features: vec![barrel.clone()],
                    limits: None,
                    nominal: Some(q("8")),
                    through: false,
                    thread: None,
                }),
            ),
            item(
                "Parallelism to B",
                ItemKind::Tolerance(GeometricTolerance {
                    characteristic: Characteristic::Parallelism,
                    value: q("0.1"),
                    diameter_zone: false,
                    material: None,
                    datums: vec!["B".into()],
                    features: vec![bottom],
                }),
            ),
        ],
        part: PartNotes::default(),
        views: Vec::new(),
    };
    let (clean, _) = write(&part, &pmi);
    if let Some(directory) = std::env::var_os("CATALLAX_PMI_OUT") {
        std::fs::write(
            std::path::Path::new(&directory).join("datum-conflict-cylinder.step"),
            &clean,
        )
        .unwrap();
    }

    // Put a datum before the real B, then make it claim B's letter and A's faces.
    let mut mixed = pmi.clone();
    mixed.items.insert(1, datum("X", barrel));
    let (file, report) = write(&part, &mixed);
    let datum = |letter| {
        let suffix = format!(",'{letter}');");
        file.lines()
            .filter_map(|line| line.split_once("=DATUM("))
            .find(|(_, fields)| fields.ends_with(&suffix))
            .unwrap()
    };
    let (a, _) = datum("A");
    let (conflict, fields) = datum("X");
    let (b, _) = datum("B");
    let number = |id: &str| id[1..].parse::<u64>().unwrap();
    assert!(number(a) < number(conflict) && number(conflict) < number(b));
    let relationship = |datum| {
        let suffix = format!(",{datum});");
        file.lines()
            .find(|line| line.contains("=SHAPE_ASPECT_RELATIONSHIP(") && line.ends_with(&suffix))
            .unwrap()
    };
    let (_, a_relationship) = relationship(a).split_once('=').unwrap();
    let (a_feature, _) = a_relationship.rsplit_once(',').unwrap();
    let (link, _) = relationship(conflict).split_once('=').unwrap();
    let file = file
        .replacen(
            &format!("{conflict}=DATUM({fields}"),
            &format!("{conflict}=DATUM({}", fields.replace("'X'", "'B'")),
            1,
        )
        .replacen(
            relationship(conflict),
            &format!("{link}={a_feature},{conflict});"),
            1,
        );
    let read = read_step_pmi(file.as_bytes()).unwrap();
    assert_eq!(
        read.skipped,
        [catallax_import::SkippedPmi {
            entity: conflict.into(),
            kind: "DATUM".into(),
            reason: "its feature already carries a datum".into(),
        }]
    );
    let faces = written_faces(&file, &report);
    let retained = read
        .pmi
        .try_map(|_, face| faces.get(face).copied().ok_or_else(|| face.clone()))
        .unwrap();
    assert_eq!(sorted(retained.clone()), sorted(pmi));
    assert_eq!(write(&part, &retained).0, clean);

    let wrong_source = file
        .lines()
        .find(|line| {
            line.contains("=DATUM_REFERENCE_COMPARTMENT(")
                && references(line).contains(&b.to_owned())
        })
        .unwrap();
    let broken = mutate(
        &file,
        wrong_source,
        &wrong_source.replace(&format!(",{b},"), &format!(",{conflict},")),
    );
    let read = read_step_pmi(broken.as_bytes()).unwrap();
    assert!(read.pmi.items.iter().any(|item| item.id == "Datum B"));
    assert!(
        read.pmi
            .items
            .iter()
            .all(|item| item.id != "Parallelism to B")
    );
    assert!(
        read.skipped
            .iter()
            .any(|item| item.reason.contains("datum B, which was not read"))
    );

    // With no valid B, the tolerance still has the same diagnostic as before.
    let without_b = file.replacen(relationship(b), "", 1);
    let read = read_step_pmi(without_b.as_bytes()).unwrap();
    let reasons: Vec<_> = read
        .skipped
        .iter()
        .filter(|skip| skip.kind == "PARALLELISM_TOLERANCE")
        .map(|skip| skip.reason.as_str())
        .collect();
    assert_eq!(reasons, ["datum B is not defined"]);

    let malformed = file.replace(",'B');", ",'b');");
    let read = read_step_pmi(malformed.as_bytes()).unwrap();
    let reasons: Vec<_> = read
        .skipped
        .iter()
        .map(|skip| skip.reason.as_str())
        .collect();
    assert_eq!(
        reasons,
        [
            "a datum letter is one or two capital letters",
            "a datum letter is one or two capital letters",
            "it refers to datum b, which was not read",
        ]
    );

    // A newly recovered datum must not recreate the old collision with the
    // temporary name used to validate each later item.
    let named_item = file.replacen("=ID_ATTRIBUTE('Datum B',", "=ID_ATTRIBUTE('item',", 1);
    let read = read_step_pmi(named_item.as_bytes()).unwrap();
    assert_eq!(read.skipped.len(), 1);
    let names: Vec<_> = read.pmi.items.iter().map(|item| item.id.as_str()).collect();
    assert_eq!(names, ["Datum A", "item", "Diameter", "Parallelism to B"]);
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
        views: Vec::new(),
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
        views: Vec::new(),
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
        views: Vec::new(),
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
        views: Vec::new(),
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

fn hole_fixture() -> (Bracket, Pmi<FaceId>) {
    use serde_json::json;

    let mut steps = vec![json!({
        "id":"stock", "op":"block", "size":["40","20","12"],
        "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}
    })];
    let cases = [
        (
            "counterbore",
            "5",
            json!("through"),
            Some((
                json!({"counterbore":{"diameter":"7.5","depth":"2.25"}}),
                HolePart::CounterboreWall,
                HoleHead::Counterbore {
                    diameter: q("7.5"),
                    depth: q("2.25"),
                },
            )),
            None,
            None,
        ),
        (
            "countersink",
            "15",
            json!("through"),
            Some((
                json!({"countersink":{"diameter":"6.25","angle":"90"}}),
                HolePart::Countersink,
                HoleHead::Countersink {
                    diameter: q("6.25"),
                    included_angle: q("90"),
                },
            )),
            None,
            None,
        ),
        (
            "blind",
            "25",
            json!({"blind":"6.5"}),
            None,
            Some(q("6.5")),
            None,
        ),
        (
            "tapped",
            "35",
            json!({"blind":"8.25"}),
            None,
            Some(q("8.25")),
            Some(q("5.5")),
        ),
    ];
    let mut target = "stock";
    for (id, x, depth, head, ..) in &cases {
        let mut step = json!({
            "id":id, "op":"hole", "target":target, "diameter":"2", "depth":depth, "tip":"flat",
            "at":[{"id":"centre","point":[x,"10","0"],"direction":["0","0","1"]}]
        });
        if let Some((recipe, ..)) = head {
            step["head"] = recipe.clone();
        }
        steps.push(step);
        target = id;
    }
    let recipe = Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"hole-pmi","steps":steps,"result":target})
            .to_string(),
    )
    .unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let part = Bracket {
        model: built.model,
        body: built.result,
    };
    let faces = part.faces();
    let named = |id: &str, part_kind: HolePart| {
        let found: Feature<FaceId> = faces
            .iter()
            .copied()
            .filter(|face| {
                part.model.face_origins(*face).iter().any(|origin| {
                    origin.step.as_str() == id
                        && matches!(origin.role, Role::Hole { part, .. } if part == part_kind)
                })
            })
            .collect();
        assert!(!found.is_empty(), "{id}: missing {part_kind:?}");
        found
    };
    let pmi = Pmi {
        items: cases
            .into_iter()
            .map(|(id, _, _, head, drill_depth, thread_depth)| {
                let mut features = vec![named(id, HolePart::Wall)];
                let head = head.map(|(_, role, head)| {
                    features.push(named(id, role));
                    head
                });
                item(
                    id,
                    ItemKind::Size(Size {
                        through: drill_depth.is_none(),
                        thread: thread_depth.as_ref().map(|_| "M2.5x0.45 - 6H".to_owned()),
                        kind: SizeKind::Hole(Box::new(HoleSize {
                            head,
                            drill_depth,
                            thread_depth,
                        })),
                        features,
                        limits: None,
                        nominal: None,
                    }),
                )
            })
            .collect(),
        part: PartNotes::default(),
        views: Vec::new(),
    };
    (part, pmi)
}

#[test]
fn compound_and_blind_holes_keep_separate_faces_and_depths() {
    let (part, pmi) = hole_fixture();
    let mut with_view = pmi.clone();
    with_view.views.push(SavedView {
        name: "Hole detail".into(),
        camera: Camera {
            eye: v("10", "-50", "5"),
            direction: v("0", "1", "0"),
            up: v("0", "0", "1"),
            distance: q("50"),
            window: [q("80"), q("40")],
            perspective: false,
        },
        items: vec![pmi.items[0].id.clone()],
    });
    assert!(
        try_write(&part, &with_view)
            .unwrap_err()
            .contains("has no callout")
    );
    let (file, report) = write(&part, &pmi);
    if let Ok(directory) = std::env::var("CATALLAX_PMI_OUT") {
        std::fs::write(
            std::path::Path::new(&directory).join("hole-pmi.step"),
            &file,
        )
        .unwrap();
    }
    let read = |file: &str, report: &catallax_step::StepReport| {
        let imported = read_step_pmi(file.as_bytes()).unwrap();
        assert_eq!(imported.skipped, Vec::new());
        let faces = written_faces(file, report);
        sorted(
            imported
                .pmi
                .try_map(|_, key| Ok::<_, ()>(faces[key]))
                .unwrap(),
        )
    };
    let mut expected = pmi;
    for item in &mut expected.items {
        let ItemKind::Size(size) = &mut item.kind else {
            unreachable!()
        };
        size.nominal = Some(q("2"));
    }
    let expected = sorted(expected);
    let imported = read(&file, &report);
    assert_eq!(imported, expected);
    let (exported, report) = write(&part, &imported);
    assert_eq!(read(&exported, &report), expected);

    let broken = file.replace("'counterbore','',#", "'unknown head','',#");
    let imported = read_step_pmi(broken.as_bytes()).unwrap();
    assert_eq!(imported.skipped.len(), 1, "{:?}", imported.skipped);
    assert!(
        imported.skipped[0].reason.contains("head relationship"),
        "{:?}",
        imported.skipped
    );
    assert_eq!(
        imported
            .pmi
            .items
            .iter()
            .map(|item| item.id.as_str())
            .collect::<Vec<_>>(),
        ["countersink", "blind", "tapped"]
    );
}

#[test]
fn a_non_equivalent_group_is_not_a_simple_diameter() {
    let (part, mut pmi) = hole_fixture();
    let pilot = |index: usize| {
        let ItemKind::Size(size) = &pmi.items[index].kind else {
            unreachable!()
        };
        size.features[0].clone()
    };
    let features = vec![pilot(0), pilot(1)];
    pmi.items.push(item(
        "pattern",
        ItemKind::Size(Size {
            kind: SizeKind::Diameter,
            features,
            limits: None,
            through: false,
            thread: None,
            nominal: None,
        }),
    ));
    let (file, _) = write(&part, &pmi);
    let all = entities(&file);
    let group = all
        .iter()
        .find(|(_, (kind, _))| kind == "COMPOSITE_GROUP_SHAPE_ASPECT")
        .unwrap()
        .0;
    let head = all
        .values()
        .find(|(kind, fields)| {
            kind == "SHAPE_ASPECT_RELATIONSHIP" && fields.starts_with("'counterbore',")
        })
        .and_then(|(_, fields)| references(fields).pop())
        .unwrap();
    let (id, (_, fields)) = all
        .iter()
        .find(|(_, (kind, fields))| {
            kind == "SHAPE_ASPECT_RELATIONSHIP" && references(fields).first() == Some(group)
        })
        .unwrap();
    let bad = file.replace(
        &format!("{id}=SHAPE_ASPECT_RELATIONSHIP({fields});"),
        &format!("{id}=SHAPE_ASPECT_RELATIONSHIP('','',{group},{head});"),
    );
    let imported = read_step_pmi(bad.as_bytes()).unwrap();
    assert_eq!(imported.skipped.len(), 1, "{:?}", imported.skipped);
    assert!(
        imported.skipped[0].reason.contains("different diameters"),
        "{:?}",
        imported.skipped
    );
    assert_eq!(imported.pmi.items.len(), 4);
    assert!(!imported.pmi.items.iter().any(|item| item.id == "pattern"));
}

fn planar_fixture(width: &str, first: &str, second: &str) -> (Bracket, [Feature<FaceId>; 4]) {
    use serde_json::json;

    let recipe = Recipe::parse(
        &json!({
            "format":"catallax-recipe/3", "name":"planar-pmi", "steps":[
                {"id":"base","op":"block","size":[width,"12","8"],
                 "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
                {"id":"bevel","op":"chamfer","body":"base","distance":first,"other_distance":second,
                 "edges":[{"between":[{"step":"base","role":{"side":{"segment":"x_min"}}},
                                      {"step":"base","role":{"side":{"segment":"y_min"}}}]}]}
            ], "result":"bevel"
        })
        .to_string(),
    )
    .unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let part = Bracket {
        model: built.model,
        body: built.result,
    };
    let bevel: Feature<FaceId> = part
        .faces()
        .into_iter()
        .filter(|id| {
            let Surface::Plane(plane) = part.model.geom().surface(part.model.face(*id).surface)
            else {
                return false;
            };
            let [x, y, z] = plane.normal.rational_components().unwrap();
            !x.is_zero() && !y.is_zero() && z.is_zero()
        })
        .collect();
    assert_eq!(bevel.len(), 1);
    let sides = [
        part.planes(0, "0"),
        part.planes(1, "0"),
        part.planes(0, width),
        bevel,
    ];
    (part, sides)
}

fn planar_pmi([x_min, y_min, x_max, bevel]: &[Feature<FaceId>; 4]) -> Pmi<FaceId> {
    Pmi {
        items: vec![
            item(
                "opposed flats",
                ItemKind::Width {
                    flats: [x_min.clone(), x_max.clone()],
                    limits: Some(Limits::Deviations {
                        lower: q("-0.1"),
                        upper: q("0.2"),
                    }),
                    nominal: None,
                },
            ),
            item(
                "chamfer",
                ItemKind::Chamfer {
                    bevel: bevel.clone(),
                    sides: [x_min.clone(), y_min.clone()],
                    nominal: None,
                },
            ),
            item(
                "edge break",
                ItemKind::EdgeBreak {
                    bevel: bevel.clone(),
                    sides: [x_min.clone(), y_min.clone()],
                    min_mm: q("0.5"),
                    max_mm: q("2.5"),
                    nominal: None,
                },
            ),
            item(
                "basic angle",
                ItemKind::BasicAngle {
                    planes: [x_min.clone(), bevel.clone()],
                    nominal: None,
                },
            ),
        ],
        part: PartNotes::default(),
        views: Vec::new(),
    }
}

fn measured_planar(
    mut pmi: Pmi<FaceId>,
    [width, setback, angle, largest, included]: [&str; 5],
) -> Pmi<FaceId> {
    for item in &mut pmi.items {
        match &mut item.kind {
            ItemKind::Width { nominal, .. } => *nominal = Some(q(width)),
            ItemKind::Chamfer { nominal, .. } => {
                *nominal = Some(ChamferNominal {
                    setback_mm: q(setback),
                    angle_deg: q(angle),
                })
            }
            ItemKind::EdgeBreak { nominal, .. } => *nominal = Some(q(largest)),
            ItemKind::BasicAngle { nominal, .. } => *nominal = Some(q(included)),
            _ => unreachable!(),
        }
    }
    pmi
}

fn read_planar(file: &str, report: &catallax_step::StepReport) -> Pmi<FaceId> {
    let read = read_step_pmi(file.as_bytes()).unwrap();
    assert!(read.skipped.is_empty(), "{:?}", read.skipped);
    let faces = written_faces(file, report);
    read.pmi
        .try_map(|_, face| faces.get(face).copied().ok_or_else(|| face.clone()))
        .unwrap()
}

#[test]
fn planar_dimensions_keep_their_roles_and_recompute_after_edit() {
    let (original, roles) = planar_fixture("10", "2", "2");
    let pmi = planar_pmi(&roles);
    let (file, report) = write(&original, &pmi);
    let read = read_planar(&file, &report);
    assert_eq!(
        sorted(read.clone()),
        sorted(measured_planar(pmi, ["10", "2", "45", "2", "135"]))
    );
    assert!(
        file.lines()
            .any(|line| line.contains("=ANGULAR_LOCATION(") && line.ends_with(",.LARGE.);"))
    );
    assert!(file.contains("PLANE_ANGLE_MEASURE(135.)"));
    if let Ok(directory) = std::env::var("CATALLAX_PMI_OUT") {
        std::fs::write(
            std::path::Path::new(&directory).join("planar-pmi.step"),
            &file,
        )
        .unwrap();
    }

    let (edited, edited_roles) = planar_fixture("14", "3", "4");
    // Carry the original PMI, including its now stale values, to the faces
    // of the edited recipe. The edge break is now outside its explicit bounds.
    let faces: BTreeMap<_, _> = roles
        .iter()
        .zip(&edited_roles)
        .flat_map(|(old, new)| {
            assert_eq!(old.len(), new.len());
            old.iter().copied().zip(new.iter().copied())
        })
        .collect();
    let moved = read
        .try_map(|_, face| faces.get(face).copied().ok_or(*face))
        .unwrap();
    let values = nominals(&edited.model, &moved, &mut Decider::default()).unwrap();
    let values: Vec<_> = values
        .into_iter()
        .map(|value| value.unwrap().as_rational().unwrap())
        .collect();
    assert_eq!(values, ["14", "3", "4", "126.869897645844"].map(q));
    let (file, report) = write(&edited, &moved);
    assert_eq!(
        sorted(read_planar(&file, &report)),
        sorted(measured_planar(
            moved,
            ["14", "3", "53.130102354156", "4", "126.869897645844"]
        )),
    );
}

#[test]
fn malformed_planar_pmi_does_not_hide_other_items() {
    let (part, roles) = planar_fixture("10", "2", "2");
    let pmi = planar_pmi(&roles);
    let (file, _) = write(&part, &pmi);
    let broken = file
        .replace(
            "PROPERTY_DEFINITION('maximum setback',",
            "PROPERTY_DEFINITION('unknown setback',",
        )
        .replace(",.LARGE.);", ",.SMALL.);");
    let imported = read_step_pmi(broken.as_bytes()).unwrap();
    assert_eq!(
        imported
            .pmi
            .items
            .iter()
            .map(|item| item.id.as_str())
            .collect::<Vec<_>>(),
        ["opposed flats", "chamfer"]
    );
    assert_eq!(imported.skipped.len(), 2);
    assert!(
        imported
            .skipped
            .iter()
            .any(|skip| skip.reason.contains("maximum setback"))
    );
    assert!(
        imported
            .skipped
            .iter()
            .any(|skip| skip.reason.contains("angle selection"))
    );

    let [x_min, y_min, x_max, bevel] = roles;
    let break_sides = [x_min.clone(), y_min.clone()];
    for kind in [
        ItemKind::Width {
            flats: [x_min.clone(), y_min.clone()],
            limits: None,
            nominal: None,
        },
        ItemKind::Chamfer {
            bevel: bevel.clone(),
            sides: [x_max.clone(), y_min],
            nominal: None,
        },
        ItemKind::BasicAngle {
            planes: [x_min, x_max],
            nominal: None,
        },
    ] {
        let invalid = Pmi {
            items: vec![item("wrong faces", kind)],
            part: PartNotes::default(),
            views: Vec::new(),
        };
        let refusal = nominals(&part.model, &invalid, &mut Decider::default()).unwrap_err();
        assert_eq!(refusal.code, "step.pmi-not-measurable");
    }
    let invalid = Pmi {
        items: vec![item(
            "wrong bounds",
            ItemKind::EdgeBreak {
                bevel,
                sides: break_sides,
                min_mm: q("2"),
                max_mm: q("1"),
                nominal: None,
            },
        )],
        part: PartNotes::default(),
        views: Vec::new(),
    };
    assert_eq!(invalid.validate().unwrap_err().code, "step.pmi-invalid");
}

/// The faces the datum target and position fixture uses on the P1 bracket.
struct GdtFaces {
    underside: FaceId,
    front: Feature<FaceId>,
    left: Feature<FaceId>,
    bore: Feature<FaceId>,
    boss_top: FaceId,
    flange_top: FaceId,
    holes: Vec<Feature<FaceId>>,
}

fn gdt_faces(part: &Bracket) -> GdtFaces {
    let underside = part.planes(2, "0");
    assert_eq!(underside.len(), 1);
    let boss_top: Vec<FaceId> = part
        .planes(2, "15")
        .into_iter()
        .filter(|face| part.corners(*face) != 4)
        .collect();
    assert_eq!(boss_top.len(), 1);
    let flange_top = part.planes(2, "30");
    assert_eq!(flange_top.len(), 1);
    GdtFaces {
        underside: underside[0],
        front: part.planes(1, "0"),
        left: part.planes(0, "0"),
        bore: part.round(2, "4", ["30", "25"]),
        boss_top: boss_top[0],
        flange_top: flange_top[0],
        holes: [["10", "12"], ["10", "32"], ["50", "12"], ["50", "32"]]
            .into_iter()
            .map(|at| part.round(2, "2", at))
            .collect(),
    }
}

fn v(x: &str, y: &str, z: &str) -> [Rational; 3] {
    [q(x), q(y), q(z)]
}

fn bounded(letter: &str, material: MaterialCondition) -> DatumRef {
    DatumRef {
        letter: letter.to_owned(),
        material: Some(material),
    }
}

fn segment(
    value: &str,
    material: Option<MaterialCondition>,
    datums: Vec<DatumRef>,
) -> PositionSegment {
    PositionSegment {
        value_mm: q(value),
        diameter_zone: true,
        material,
        datums,
    }
}

fn target(number: u32, face: FaceId, at: [Rational; 3], kind: TargetKind) -> ItemKind<FaceId> {
    ItemKind::DatumTarget(DatumTarget {
        datum: "A".into(),
        number,
        face,
        point_mm: at,
        kind,
    })
}

/// Datum A on the plate's underside, set up by a point, a line, a circular and a
/// rectangular target clear of its five holes; datums B, C and the bore D; a
/// composite frame on the bolt holes whose segments carry their own modifiers
/// and datum boundaries; a zone projected 10 up out of the boss top on the bore,
/// and one projected 8 down out of the underside on the bolt holes.
fn gdt_pmi(f: &GdtFaces) -> Pmi<FaceId> {
    use MaterialCondition::{Least, Maximum};
    let datum = |letter: &str, feature: &Feature<FaceId>| {
        item(
            &format!("Datum {letter}"),
            ItemKind::Datum(Datum {
                letter: letter.to_owned(),
                feature: feature.clone(),
            }),
        )
    };
    let x = || v("1", "0", "0");
    Pmi {
        items: vec![
            datum("A", &vec![f.underside]),
            item(
                "A1",
                target(1, f.underside, v("5", "5", "0"), TargetKind::Point),
            ),
            item(
                "A2",
                target(
                    2,
                    f.underside,
                    v("30", "5", "0"),
                    TargetKind::Line {
                        length_mm: q("20"),
                        direction: x(),
                    },
                ),
            ),
            item(
                "A3",
                target(
                    3,
                    f.underside,
                    v("55", "36", "0"),
                    TargetKind::Area(TargetArea::Circular {
                        diameter_mm: q("6"),
                    }),
                ),
            ),
            item(
                "A4",
                target(
                    4,
                    f.underside,
                    v("30", "36", "0"),
                    TargetKind::Area(TargetArea::Rectangular {
                        length_mm: q("10"),
                        width_mm: q("4"),
                        direction: x(),
                    }),
                ),
            ),
            datum("B", &f.front),
            datum("C", &f.left),
            datum("D", &f.bore),
            item(
                "Bore",
                ItemKind::Size(Size {
                    kind: SizeKind::Diameter,
                    features: vec![f.bore.clone()],
                    limits: Some(Limits::Fit {
                        deviation: "H".into(),
                        grade: "7".into(),
                    }),
                    through: true,
                    thread: None,
                    nominal: None,
                }),
            ),
            item(
                "Bolt pattern",
                ItemKind::Position(PositionControl {
                    features: f.holes.clone(),
                    frame: PositionFrame::Composite {
                        upper: segment(
                            "0.4",
                            Some(Maximum),
                            vec![
                                DatumRef::new("A"),
                                bounded("D", Maximum),
                                DatumRef::new("B"),
                            ],
                        ),
                        lower: segment(
                            "0.1",
                            Some(Least),
                            vec![DatumRef::new("A"), DatumRef::new("D")],
                        ),
                    },
                }),
            ),
            item(
                "Bore projection",
                ItemKind::Position(PositionControl {
                    features: vec![f.bore.clone()],
                    frame: PositionFrame::Single {
                        segment: segment(
                            "0.05",
                            None,
                            vec![DatumRef::new("B"), DatumRef::new("C")],
                        ),
                        projected_zone: Some(ProjectedZone {
                            height_mm: q("10"),
                            direction: v("0", "0", "1"),
                            entry_faces: vec![f.boss_top],
                        }),
                    },
                }),
            ),
            item(
                "Bolt projection",
                ItemKind::Position(PositionControl {
                    features: f.holes.clone(),
                    frame: PositionFrame::Single {
                        segment: segment(
                            "0.2",
                            Some(Maximum),
                            vec![DatumRef::new("A"), DatumRef::new("B"), DatumRef::new("C")],
                        ),
                        projected_zone: Some(ProjectedZone {
                            height_mm: q("8"),
                            direction: v("0", "0", "-1"),
                            entry_faces: vec![f.underside; 4],
                        }),
                    },
                }),
            ),
        ],
        part: PartNotes {
            standard: Some(catallax_step::pmi::DimensioningStandard::AsmeY14_5_2018),
            ..PartNotes::default()
        },
        views: Vec::new(),
    }
}

/// The items in name order, each feature's faces in one order.
fn by_name(pmi: Pmi<FaceId>) -> Pmi<FaceId> {
    let mut pmi = sorted(pmi);
    pmi.items.sort_by(|a, b| a.id.cmp(&b.id));
    pmi
}

fn assert_same(back: Pmi<FaceId>, expected: Pmi<FaceId>) {
    let (back, expected) = (by_name(back), by_name(expected));
    for (got, want) in back.items.iter().zip(&expected.items) {
        assert_eq!(got, want);
    }
    assert_eq!(back, expected);
}

/// What a file reads back as, with faces named by the model faces written.
fn read_back(
    file: &str,
    report: &catallax_step::StepReport,
) -> (Pmi<FaceId>, Vec<catallax_import::SkippedPmi>) {
    let read = read_step_pmi(file.as_bytes()).unwrap();
    let faces = written_faces(file, report);
    let pmi = read
        .pmi
        .try_map(|_, key| faces.get(key).copied().ok_or_else(|| key.clone()))
        .unwrap();
    (pmi, read.skipped)
}

fn try_write(part: &Bracket, pmi: &Pmi<FaceId>) -> Result<String, String> {
    let mut bytes = Vec::new();
    write_step_with_pmi(
        &part.model,
        part.body,
        &StepOptions::default(),
        pmi,
        &mut Decider::default(),
        &mut bytes,
    )
    .map(|_| String::from_utf8(bytes).unwrap())
    .map_err(|refusal| refusal.detail)
}

fn kind_of<'a, F>(pmi: &'a Pmi<F>, id: &str) -> &'a ItemKind<F> {
    &pmi.items
        .iter()
        .find(|item| item.id == id)
        .unwrap_or_else(|| panic!("{id} was not read"))
        .kind
}

#[test]
fn datum_targets_and_position_frames_read_back_on_their_faces() {
    let part = bracket();
    let f = gdt_faces(&part);
    let pmi = gdt_pmi(&f);
    pmi.validate().unwrap();
    let (file, report) = write(&part, &pmi);
    if let Ok(directory) = std::env::var("CATALLAX_PMI_OUT") {
        std::fs::write(
            std::path::Path::new(&directory).join("c2b-targets-positions.step"),
            &file,
        )
        .unwrap();
    }
    let (back, skipped) = read_back(&file, &report);
    assert_eq!(skipped, Vec::new());

    // Each target on the underside, with its own number, point and sizes.
    for (id, number) in [("A1", 1), ("A2", 2), ("A3", 3), ("A4", 4)] {
        let ItemKind::DatumTarget(target) = kind_of(&back, id) else {
            panic!("{id} is not a datum target");
        };
        assert_eq!(
            (target.datum.as_str(), target.number, target.face),
            ("A", number, f.underside)
        );
    }
    assert!(matches!(
        kind_of(&back, "A4"),
        ItemKind::DatumTarget(DatumTarget {
            kind: TargetKind::Area(TargetArea::Rectangular { length_mm, width_mm, .. }),
            ..
        }) if *length_mm == q("10") && *width_mm == q("4")
    ));

    // The composite keeps its order, each segment's own modifier and its own
    // datum boundaries: D(M) above, plain D below.
    let ItemKind::Position(PositionControl {
        features,
        frame: PositionFrame::Composite { upper, lower },
    }) = kind_of(&back, "Bolt pattern")
    else {
        panic!("the bolt pattern is not a composite frame");
    };
    assert_eq!(features.len(), 4);
    assert_eq!(
        (upper.value_mm.clone(), upper.material),
        (q("0.4"), Some(MaterialCondition::Maximum))
    );
    assert_eq!(
        (lower.value_mm.clone(), lower.material),
        (q("0.1"), Some(MaterialCondition::Least))
    );
    assert_eq!(upper.datums[1], bounded("D", MaterialCondition::Maximum));
    assert_eq!(lower.datums, vec![DatumRef::new("A"), DatumRef::new("D")]);

    // Each projected zone starts on its own entry face and points out of it.
    for (id, entry, height, up) in [
        ("Bore projection", f.boss_top, "10", "1"),
        ("Bolt projection", f.underside, "8", "-1"),
    ] {
        let ItemKind::Position(PositionControl {
            frame:
                PositionFrame::Single {
                    projected_zone: Some(zone),
                    ..
                },
            ..
        }) = kind_of(&back, id)
        else {
            panic!("{id} has no projected zone");
        };
        assert!(zone.entry_faces.iter().all(|face| *face == entry), "{id}");
        assert_eq!(zone.height_mm, q(height), "{id}");
        assert!(
            zone.direction[0].is_zero() && zone.direction[1].is_zero(),
            "{id}"
        );
        assert_eq!(
            zone.direction[2] > Rational::zero(),
            q(up) > Rational::zero(),
            "{id}"
        );
    }

    // Everything else is what went out, the bore with its measured size.
    let mut expected = pmi.clone();
    for item in &mut expected.items {
        if let ItemKind::Size(size) = &mut item.kind {
            size.nominal = Some(q("8"));
        }
    }
    assert_same(back.clone(), expected);

    // The file counts one composite of two tolerances among four, and four targets.
    let counts: Vec<String> = validation_items(&file)
        .into_values()
        .flatten()
        .map(|(_, fields)| fields)
        .collect();
    for count in [
        "'number of semantic pmi elements',12.",
        "'number of geometric tolerances',4.",
        "'number of composite tolerances',1.",
        "'number of datums',4.",
        "'number of datum targets',4.",
    ] {
        assert!(counts.iter().any(|fields| fields == count), "{count}");
    }

    // Writing what was read writes the same information again.
    let mut again = back;
    for item in &mut again.items {
        if let ItemKind::Size(size) = &mut item.kind {
            size.nominal = None;
        }
    }
    let (file, report) = write(&part, &again);
    let (twice, skipped) = read_back(&file, &report);
    assert_eq!(skipped, Vec::new());
    let mut expected = pmi;
    for item in &mut expected.items {
        if let ItemKind::Size(size) = &mut item.kind {
            size.nominal = Some(q("8"));
        }
    }
    assert_same(twice, expected);
}

#[test]
fn malformed_targets_and_frames_are_refused_before_writing() {
    let part = bracket();
    let f = gdt_faces(&part);
    let wrong = std::cell::RefCell::new(Vec::new());
    let refused = |change: &dyn Fn(&mut Pmi<FaceId>), reason: &str| {
        let mut pmi = gdt_pmi(&f);
        change(&mut pmi);
        match try_write(&part, &pmi) {
            Err(detail) if detail.contains(reason) => {}
            Err(detail) => wrong
                .borrow_mut()
                .push(format!("{reason}: refused with {detail}")),
            Ok(_) => wrong.borrow_mut().push(format!("{reason}: written")),
        }
    };
    let set = |pmi: &mut Pmi<FaceId>, id: &str, kind: ItemKind<FaceId>| {
        pmi.items
            .iter_mut()
            .find(|item| item.id == id)
            .unwrap()
            .kind = kind;
    };
    let with_frame =
        |pmi: &mut Pmi<FaceId>, id: &str, change: &dyn Fn(&mut PositionControl<FaceId>)| {
            let item = pmi.items.iter_mut().find(|item| item.id == id).unwrap();
            let ItemKind::Position(control) = &mut item.kind else {
                unreachable!()
            };
            change(control);
        };
    let composite = |change: fn(&mut PositionSegment, &mut PositionSegment)| {
        move |control: &mut PositionControl<FaceId>| {
            if let PositionFrame::Composite { upper, lower } = &mut control.frame {
                change(upper, lower);
            }
        }
    };
    let projection = |change: fn(&mut ProjectedZone<FaceId>)| {
        move |control: &mut PositionControl<FaceId>| {
            if let PositionFrame::Single {
                projected_zone: Some(zone),
                ..
            } = &mut control.frame
            {
                change(zone);
            }
        }
    };

    // Composite frames.
    refused(
        &|pmi| {
            with_frame(
                pmi,
                "Bolt pattern",
                &composite(|_, lower| lower.value_mm = q("0.5")),
            )
        },
        "cannot be larger than the upper one",
    );
    refused(
        &|pmi| {
            with_frame(
                pmi,
                "Bolt pattern",
                &composite(|_, lower| lower.datums = vec![DatumRef::new("B")]),
            )
        },
        "repeat the upper segment's first datums in order",
    );
    refused(
        &|pmi| {
            with_frame(
                pmi,
                "Bolt pattern",
                &composite(|upper, _| upper.datums.clear()),
            )
        },
        "the upper segment of a composite needs datums",
    );
    refused(
        &|pmi| {
            with_frame(
                pmi,
                "Bolt pattern",
                &composite(|upper, _| upper.diameter_zone = false),
            )
        },
        "zone is cylindrical",
    );
    // Datum boundaries: only on a feature of size that no target sets up.
    refused(
        &|pmi| {
            with_frame(
                pmi,
                "Bolt pattern",
                &composite(|upper, _| {
                    upper.datums[2] = bounded("B", MaterialCondition::Maximum);
                }),
            )
        },
        "datum B is not a feature of size",
    );
    refused(
        &|pmi| {
            with_frame(
                pmi,
                "Bolt pattern",
                &composite(|upper, _| {
                    upper.datums[0] = bounded("A", MaterialCondition::Least);
                }),
            )
        },
        "datum A is set up by datum targets",
    );
    // Projected zones.
    refused(
        &|pmi| {
            with_frame(
                pmi,
                "Bore projection",
                &projection(|zone| zone.height_mm = q("0")),
            )
        },
        "a projected zone needs a positive height",
    );
    refused(
        &|pmi| {
            with_frame(
                pmi,
                "Bore projection",
                &projection(|zone| zone.direction = v("0", "0", "-1")),
            )
        },
        "does not point out of the part",
    );
    refused(
        &|pmi| {
            with_frame(
                pmi,
                "Bore projection",
                &projection(|zone| zone.direction = v("1", "0", "1")),
            )
        },
        "does not point out of the part",
    );
    refused(
        &|pmi| {
            let flange = f.flange_top;
            with_frame(
                pmi,
                "Bore projection",
                &move |control: &mut PositionControl<FaceId>| {
                    if let PositionFrame::Single {
                        projected_zone: Some(zone),
                        ..
                    } = &mut control.frame
                    {
                        zone.entry_faces = vec![flange];
                    }
                },
            )
        },
        "needs exactly one projected zone entry face it meets",
    );
    refused(
        &|pmi| {
            with_frame(
                pmi,
                "Bolt projection",
                &projection(|zone| {
                    zone.entry_faces.pop();
                }),
            )
        },
        "one entry face per controlled feature",
    );
    refused(
        &|pmi| {
            with_frame(pmi, "Bore projection", &|control: &mut PositionControl<
                FaceId,
            >| {
                if let PositionFrame::Single { projected_zone, .. } = &mut control.frame {
                    *projected_zone = None;
                }
            })
        },
        "is a geometric tolerance",
    );
    // Datum targets.
    let line = |direction: [Rational; 3]| TargetKind::Line {
        length_mm: q("10"),
        direction,
    };
    let cases: Vec<(&str, ItemKind<FaceId>, &str)> = vec![
        (
            "A1",
            target(0, f.underside, v("5", "5", "0"), TargetKind::Point),
            "number is positive",
        ),
        (
            "A1",
            target(2, f.underside, v("5", "5", "0"), TargetKind::Point),
            "datum target A2 is defined twice",
        ),
        (
            "A1",
            target(1, f.boss_top, v("30", "25", "15"), TargetKind::Point),
            "not on its datum's feature",
        ),
        (
            "A1",
            target(1, f.underside, v("5", "5", "1"), TargetKind::Point),
            "not on its face's plane",
        ),
        (
            "A1",
            target(1, f.underside, v("10", "12", "0"), TargetKind::Point),
            "the target point is outside its face",
        ),
        (
            "A2",
            target(2, f.underside, v("10", "12", "0"), line(v("1", "0", "0"))),
            "outside its face",
        ),
        (
            "A2",
            target(2, f.underside, v("10", "8", "0"), line(v("0", "0", "1"))),
            "does not run along its face",
        ),
        (
            "A2",
            target(2, f.underside, v("10", "6", "0"), line(v("0", "1", "0"))),
            "the target line leaves its face",
        ),
        (
            "A2",
            target(
                2,
                f.underside,
                v("10", "16", "0"),
                TargetKind::Line {
                    length_mm: q("12"),
                    direction: v("0", "1", "0"),
                },
            ),
            "the target line leaves its face",
        ),
        (
            "A3",
            target(
                3,
                f.underside,
                v("14", "12", "0"),
                TargetKind::Area(TargetArea::Circular {
                    diameter_mm: q("6"),
                }),
            ),
            "the target area crosses its face's boundary",
        ),
        (
            "A3",
            target(
                3,
                f.underside,
                v("58", "36", "0"),
                TargetKind::Area(TargetArea::Circular {
                    diameter_mm: q("6"),
                }),
            ),
            "the target area crosses its face's boundary",
        ),
        (
            "A4",
            target(
                4,
                f.underside,
                v("30", "21", "0"),
                TargetKind::Area(TargetArea::Rectangular {
                    length_mm: q("14"),
                    width_mm: q("4"),
                    direction: v("1", "0", "0"),
                }),
            ),
            "the target area crosses its face's boundary",
        ),
        (
            "A4",
            target(
                4,
                f.underside,
                v("30", "30", "0"),
                TargetKind::Area(TargetArea::Rectangular {
                    length_mm: q("20"),
                    width_mm: q("20"),
                    direction: v("1", "0", "0"),
                }),
            ),
            "the target area crosses its face's boundary",
        ),
        (
            "A3",
            target(
                3,
                f.underside,
                v("55", "36", "0"),
                TargetKind::Area(TargetArea::Circular {
                    diameter_mm: q("0"),
                }),
            ),
            "sizes must be more than zero",
        ),
    ];
    for (id, kind, reason) in cases {
        refused(&|pmi| set(pmi, id, kind.clone()), reason);
    }
    refused(
        &|pmi| {
            set(
                pmi,
                "A1",
                ItemKind::DatumTarget(DatumTarget {
                    datum: "E".into(),
                    number: 1,
                    face: f.underside,
                    point_mm: v("5", "5", "0"),
                    kind: TargetKind::Point,
                }),
            )
        },
        "datum E is not defined",
    );
    assert_eq!(wrong.into_inner(), Vec::<String>::new());
}

/// The file with `from` replaced once by `to`.
fn mutate(file: &str, from: &str, to: &str) -> String {
    assert_eq!(file.matches(from).count(), 1, "{from}");
    file.replacen(from, to, 1)
}

/// The entity number of the one line holding all of `parts`.
fn line_with(file: &str, parts: &[&str]) -> String {
    let lines: Vec<&str> = file
        .lines()
        .filter(|line| parts.iter().all(|part| line.contains(part)))
        .collect();
    assert_eq!(lines.len(), 1, "{parts:?}");
    lines[0].split('=').next().unwrap().to_owned()
}

#[test]
fn a_bad_target_or_frame_on_import_keeps_the_independent_items() {
    let part = bracket();
    let f = gdt_faces(&part);
    let (file, report) = write(&part, &gdt_pmi(&f));
    let ids = |pmi: &Pmi<FaceId>| {
        let mut ids: Vec<String> = pmi.items.iter().map(|item| item.id.clone()).collect();
        ids.sort();
        ids
    };
    let reasons = |skipped: &[catallax_import::SkippedPmi]| {
        skipped
            .iter()
            .map(|skip| skip.reason.clone())
            .collect::<Vec<_>>()
            .join("; ")
    };

    // A1 loses its number: datum A and everything that needs it are left out,
    // with their reasons; B, C, D, the bore size and the bore projection stay.
    let a1 = line_with(&file, &["PLACED_DATUM_TARGET_FEATURE", "'A1'"]);
    let record = file
        .lines()
        .find(|line| line.starts_with(&format!("{a1}=")))
        .unwrap();
    let broken = mutate(&file, record, &record.replace(",'1');", ",'0');"));
    let (back, skipped) = read_back(&broken, &report);
    assert_eq!(
        ids(&back),
        ["Bore", "Bore projection", "Datum B", "Datum C", "Datum D"],
        "{}",
        reasons(&skipped)
    );
    let text = reasons(&skipped);
    for reason in [
        "its target number is not a positive whole number",
        "was not read",
        "it refers to datum A, which was not read",
        "lower segment of the composite frame",
    ] {
        assert!(text.contains(reason), "{text} lacks {reason}");
    }
    assert!(matches!(
        kind_of(&back, "Bore projection"),
        ItemKind::Position(_)
    ));

    let dangling = mutate(&file, &format!("{record}\n"), "");
    let (back, skipped) = read_back(&dangling, &report);
    assert_eq!(
        ids(&back),
        ["Bore", "Bore projection", "Datum B", "Datum C", "Datum D"]
    );
    assert!(reasons(&skipped).contains(&format!("its datum target {a1} was not read")));

    let relationship = file
        .lines()
        .find(|line| {
            line.contains("SHAPE_ASPECT_RELATIONSHIP(")
                && line.contains("'datum target'")
                && line.contains(&format!(",{a1},"))
        })
        .unwrap();
    for orphaned in [
        mutate(&file, &format!("{relationship}\n"), ""),
        file.lines()
            .filter(|line| !line.contains("=SHAPE_ASPECT_RELATIONSHIP('','datum target'"))
            .collect::<Vec<_>>()
            .join("\n"),
    ] {
        let (back, skipped) = read_back(&orphaned, &report);
        assert_eq!(
            ids(&back),
            ["Bore", "Bore projection", "Datum B", "Datum C", "Datum D"],
            "{}",
            reasons(&skipped)
        );
        assert!(reasons(&skipped).contains("does not belong to exactly one datum"));
    }
    let endpoint = references(relationship.split_once('=').unwrap().1);
    let reversed = mutate(
        &file,
        relationship,
        &relationship.replace(
            &format!("{},{}", endpoint[0], endpoint[1]),
            &format!("{},{}", endpoint[1], endpoint[0]),
        ),
    );
    let (back, skipped) = read_back(&reversed, &report);
    assert_eq!(
        ids(&back),
        ["Bore", "Bore projection", "Datum B", "Datum C", "Datum D"]
    );
    assert!(reasons(&skipped).contains("datum target relationship is unreadable"));

    let written = records(&file);
    let rim = written
        .iter()
        .find_map(|(id, record)| {
            if !record.starts_with("CIRCLE(") || !record.ends_with(",2.);") {
                return None;
            }
            let frame = &references(record)[0];
            let location = &references(&written[frame])[0];
            (written[location] == "CARTESIAN_POINT('',(10.,12.,0.));").then_some(id)
        })
        .unwrap();
    let bad_rim = mutate(
        &file,
        &format!("{rim}={}", written[rim]),
        &format!("{rim}={}", written[rim].replace(",2.);", ",1.);")),
    );
    let bad_rim = mutate(
        &bad_rim,
        "CARTESIAN_POINT('',(5.,5.,0.))",
        "CARTESIAN_POINT('',(11.5,12.,0.))",
    );
    let (back, skipped) = read_back(&bad_rim, &report);
    assert!(!ids(&back).contains(&"A1".to_owned()));
    assert!(
        skipped
            .iter()
            .any(|item| item.entity == a1 && item.reason.contains("boundary"))
    );
    let file_faces = written_faces(&file, &report);
    for (face, affected) in [
        (f.holes[0][0], "Bolt pattern"),
        (f.bore[0], "Bore projection"),
    ] {
        let id = file_faces
            .iter()
            .find(|(_, native)| **native == face)
            .unwrap()
            .0;
        let surface = references(&written[id]).pop().unwrap();
        let cylinder = &written[&surface];
        let (attributes, radius) = cylinder.rsplit_once(',').unwrap();
        let broken = mutate(
            &file,
            &format!("{surface}={cylinder}"),
            &format!("{surface}={attributes},0.);"),
        );
        assert!(radius.ends_with(");"));
        let (back, skipped) = read_back(&broken, &report);
        assert!(!ids(&back).contains(&affected.to_owned()));
        assert!(reasons(&skipped).contains("each one cylinder"));

        let placement = references(cylinder)[0].clone();
        let broken = mutate(
            &file,
            &format!("{surface}={cylinder}"),
            &format!(
                "{surface}={}",
                cylinder.replace(&format!(",{placement},"), &format!(",{id},"))
            ),
        );
        let (back, _) = read_back(&broken, &report);
        assert!(!ids(&back).contains(&affected.to_owned()));
        let placement_record = &written[&placement];
        let axes = references(placement_record);
        let broken = mutate(
            &file,
            &format!("{placement}={placement_record}"),
            &format!(
                "{placement}={}",
                placement_record.replace(
                    &format!("{},{}", axes[1], axes[2]),
                    &format!("{0},{0}", axes[1])
                )
            ),
        );
        let (back, _) = read_back(&broken, &report);
        assert!(!ids(&back).contains(&affected.to_owned()));
    }

    // A projected height of zero leaves out that control alone.
    let zones: Vec<&str> = file
        .lines()
        .filter(|line| line.contains("PROJECTED_ZONE_DEFINITION"))
        .collect();
    assert_eq!(zones.len(), 2);
    let bore_zone = zones
        .iter()
        .copied()
        .find(|zone| {
            let height = zone.trim_end_matches(");").rsplit('#').next().unwrap();
            file.lines()
                .any(|line| line.starts_with(&format!("#{height}=")) && line.contains("(10.)"))
        })
        .unwrap();
    let zone_refs = references(bore_zone.split_once('=').unwrap().1);
    let swapped = mutate(
        &file,
        bore_zone,
        &bore_zone.replace(
            &format!("{},(),{}", zone_refs[0], zone_refs[1]),
            &format!("{},(),{}", zone_refs[1], zone_refs[0]),
        ),
    );
    let (back, skipped) = read_back(&swapped, &report);
    assert!(!ids(&back).contains(&"Bore projection".into()));
    assert!(reasons(&skipped).contains("must name one tolerance zone"));

    let zone = &zone_refs[0];
    let zone_record = &written[zone];
    let links = references(zone_record);
    let [_, tol, form] = links.as_slice() else {
        panic!("{zone_record}");
    };
    let broken = mutate(
        &file,
        &format!("{zone}={zone_record}"),
        &format!(
            "{zone}={}",
            zone_record.replace(&format!("({tol}),{form}"), &format!("({form}),{tol}"))
        ),
    );
    let (back, skipped) = read_back(&broken, &report);
    assert!(!ids(&back).contains(&"Bore projection".into()));
    assert!(reasons(&skipped).contains("invalid defining tolerance list"));

    let composite = file
        .lines()
        .find(|line| line.contains("GEOMETRIC_TOLERANCE_RELATIONSHIP('composite',"))
        .unwrap();
    let endpoints = references(composite.split_once('=').unwrap().1);
    let broken = mutate(
        &file,
        composite,
        &composite.replace(
            &format!("'',{},{}", endpoints[0], endpoints[1]),
            &format!("{},'',{}", endpoints[0], endpoints[1]),
        ),
    );
    let read = read_step_pmi(broken.as_bytes()).unwrap();
    assert!(
        endpoints
            .iter()
            .all(|id| read.skipped.iter().any(|item| &item.entity == id))
    );
    assert!(!read.sources.contains_key("Bolt pattern"));

    let height = bore_zone.trim_end_matches(");").rsplit('#').next().unwrap();
    let record = file
        .lines()
        .find(|line| line.starts_with(&format!("#{height}=")))
        .unwrap();
    let broken = mutate(&file, record, &record.replace("(10.)", "(0.)"));
    let (back, skipped) = read_back(&broken, &report);
    assert!(
        !ids(&back).contains(&"Bore projection".to_owned()),
        "{}",
        reasons(&skipped)
    );
    assert_eq!(back.items.len(), 11, "{}", reasons(&skipped));
    assert!(reasons(&skipped).contains("a projected zone needs a positive height"));

    for select in ["PLANE_ANGLE_MEASURE(10.)", "LENGTH_MEASURE(10.,11.)"] {
        let broken = mutate(
            &file,
            record,
            &record.replace("LENGTH_MEASURE(10.)", select),
        );
        let (back, skipped) = read_back(&broken, &report);
        assert!(!ids(&back).contains(&"Bore projection".to_owned()));
        assert!(reasons(&skipped).contains("height is not a length"));
    }

    let mut unspecified = gdt_pmi(&f);
    unspecified.part.standard = None;
    assert!(
        try_write(&part, &unspecified)
            .unwrap_err()
            .contains("explicit selection")
    );
    let scope = &f.front;
    let incompatible = item(
        "ISO symmetry",
        ItemKind::Tolerance(GeometricTolerance {
            characteristic: Characteristic::Symmetry,
            value: q("0.1"),
            diameter_zone: false,
            material: None,
            datums: vec!["A".into()],
            features: vec![scope.clone()],
        }),
    );
    unspecified.items.retain(|item| item.id != "Bolt pattern");
    unspecified.items.push(incompatible.clone());
    assert!(
        try_write(&part, &unspecified)
            .unwrap()
            .contains("DOCUMENT('ISO 1101:2017'")
    );
    let mut mixed = gdt_pmi(&f);
    mixed.items.push(incompatible);
    assert!(
        try_write(&part, &mixed)
            .unwrap_err()
            .contains("ISO symmetry")
    );

    let iso = mutate(&file, "ASME Y14.5-2018", "ISO 1101:2017");
    let (back, skipped) = read_back(&iso, &report);
    assert!(!ids(&back).contains(&"Bolt pattern".to_owned()));
    assert!(ids(&back).contains(&"Bore projection".to_owned()));
    assert!(reasons(&skipped).contains("requires a part declaration of ASME Y14.5-2018"));
    let all = records(&file);
    let (target_id, target) = all
        .iter()
        .find(|(_, fields)| fields.starts_with("PLACED_DATUM_TARGET_FEATURE('A1'"))
        .unwrap();
    let owner = &references(target)[0];
    let attachment = all
        .values()
        .find(|fields| {
            fields.starts_with("FEATURE_FOR_DATUM_TARGET_RELATIONSHIP(")
                && references(fields)[1] == *target_id
        })
        .unwrap();
    let aspect_id = &references(attachment)[0];
    let aspect = &all[aspect_id];
    let bad_owner = mutate(&file, target, &target.replace(owner, "$"));
    let bad_owner = mutate(
        &bad_owner,
        &format!("{aspect_id}={aspect}"),
        &format!("{aspect_id}={}", aspect.replace(owner, "$")),
    );
    let (back, skipped) = read_back(&bad_owner, &report);
    assert!(!ids(&back).contains(&"Datum A".to_owned()));
    assert!(ids(&back).contains(&"Bore projection".to_owned()));
    assert!(reasons(&skipped).contains("target has no owning part"));

    let orientation = all
        .values()
        .find(|fields| {
            fields.starts_with("AXIS2_PLACEMENT_3D('orientation'")
                && all[&references(fields)[0]].ends_with("(30.,5.,0.));")
        })
        .unwrap();
    let ortho = &references(orientation)[2];
    let tilted = mutate(&file, orientation, &orientation.replace(ortho, "#900010"));
    let tilted = mutate(
        &tilted,
        "ENDSEC;\nEND-ISO",
        "#900010=DIRECTION('',(1.,0.,1.));\nENDSEC;\nEND-ISO",
    );
    let (back, skipped) = read_back(&tilted, &report);
    assert!(skipped.is_empty(), "{skipped:?}");
    let ItemKind::DatumTarget(target) = kind_of(&back, "A2") else {
        panic!("expected line target")
    };
    assert!(
        matches!(&target.kind, TargetKind::Line { direction, .. } if direction == &v("1", "0", "0"))
    );

    let projection = all
        .values()
        .find(|fields| {
            fields.starts_with("PROJECTED_ZONE_DEFINITION(")
                && all[&references(fields)[2]].contains("(10.)")
        })
        .unwrap();
    let zone = &references(projection)[0];
    let wrapped = mutate(
        &file,
        projection,
        &projection.replacen(zone, &format!("({zone})"), 1),
    );
    let (back, skipped) = read_back(&wrapped, &report);
    assert!(!ids(&back).contains(&"Bore projection".to_owned()));
    assert!(ids(&back).contains(&"Bolt pattern".to_owned()));
    assert!(reasons(&skipped).contains("must name one tolerance zone"));
    let specialized = mutate(
        &file,
        &all[zone],
        &all[zone].replace("TOLERANCE_ZONE(", "DIRECTED_TOLERANCE_ZONE("),
    );
    let (back, skipped) = read_back(&specialized, &report);
    assert!(!ids(&back).contains(&"Bore projection".to_owned()));
    assert!(reasons(&skipped).contains("unsupported shape attributes"));
    let relationship = all
        .values()
        .find(|body| body.starts_with("GEOMETRIC_TOLERANCE_RELATIONSHIP('composite'"))
        .unwrap();
    let upper = &references(relationship)[0];
    let wrapped = mutate(
        &file,
        relationship,
        &relationship.replacen(upper, &format!("({upper})"), 1),
    );
    let (back, skipped) = read_back(&wrapped, &report);
    assert!(!ids(&back).contains(&"Bolt pattern".to_owned()));
    assert!(ids(&back).contains(&"Bore projection".to_owned()));
    assert_eq!(
        reasons(&skipped)
            .matches("no upper position segment")
            .count(),
        2
    );
    let lower = &all[&references(relationship)[1]];
    let marker = "GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE((";
    let datum = lower
        .split_once(marker)
        .unwrap()
        .1
        .split(')')
        .next()
        .unwrap();
    let poison = lower.replace(&format!("{marker}{datum}"), &format!("{marker}{datum},$"));
    let malformed = mutate(&file, lower, &poison);
    let (back, skipped) = read_back(&malformed, &report);
    assert!(!ids(&back).contains(&"Bolt pattern".to_owned()));
    assert!(ids(&back).contains(&"Bore projection".to_owned()));
    assert!(reasons(&skipped).contains("invalid datum reference aggregate"));

    // An unknown relationship between the segments leaves out both of them.
    let broken = mutate(
        &file,
        "GEOMETRIC_TOLERANCE_RELATIONSHIP('composite'",
        "GEOMETRIC_TOLERANCE_RELATIONSHIP('precedence'",
    );
    let (back, skipped) = read_back(&broken, &report);
    assert!(!ids(&back).contains(&"Bolt pattern".to_owned()));
    assert_eq!(back.items.len(), 11, "{}", reasons(&skipped));
    let text = reasons(&skipped);
    assert!(
        text.contains("`precedence` tolerance relationship"),
        "{text}"
    );
    assert_eq!(
        text.matches("`precedence` tolerance relationship").count(),
        2,
        "{text}"
    );

    for extra in [
        ".LEAST_MATERIAL_REQUIREMENT.",
        ".REGARDLESS_OF_FEATURE_SIZE.",
    ] {
        let upper = file
            .lines()
            .find(|line| {
                line.contains("GEOMETRIC_TOLERANCE('Bolt pattern'")
                    && line.contains(".MAXIMUM_MATERIAL_REQUIREMENT.")
            })
            .unwrap();
        let bad = mutate(
            &file,
            upper,
            &upper.replace(
                ".MAXIMUM_MATERIAL_REQUIREMENT.",
                &format!(".MAXIMUM_MATERIAL_REQUIREMENT.,{extra}"),
            ),
        );
        let (back, skipped) = read_back(&bad, &report);
        assert!(!ids(&back).contains(&"Bolt pattern".to_owned()));
        assert!(reasons(&skipped).contains("material modifier is stated more than once"));
    }

    let composite = file
        .lines()
        .find(|line| line.contains("GEOMETRIC_TOLERANCE_RELATIONSHIP('composite'"))
        .unwrap();
    let [upper, lower] =
        <[String; 2]>::try_from(references(composite.split_once('=').unwrap().1)).unwrap();
    let second_upper = &read_step_pmi(file.as_bytes()).unwrap().sources["Bore projection"];
    let owner = file
        .lines()
        .find_map(|line| {
            line.split_once("=PRODUCT_DEFINITION_SHAPE(")
                .map(|(id, _)| id)
        })
        .unwrap();
    for (additional, refused, reason) in [
        (
            format!(
                "#900000=GEOMETRIC_TOLERANCE_RELATIONSHIP('composite','',{second_upper},{lower});\n"
            ),
            "Bolt pattern",
            "two uniquely paired segments",
        ),
        (
            format!(
                "#900000=TOLERANCE_ZONE('','',{owner},.F.,({upper}),#900001);\n#900001=TOLERANCE_ZONE_FORM('between two parallel planes');\n"
            ),
            "Bolt pattern",
            "more than one tolerance zone",
        ),
    ] {
        let broken = mutate(
            &file,
            "ENDSEC;\nEND-ISO",
            &format!("{additional}ENDSEC;\nEND-ISO"),
        );
        let (back, skipped) = read_back(&broken, &report);
        assert!(!ids(&back).contains(&refused.to_owned()));
        assert!(reasons(&skipped).contains(reason));
    }

    let projection = references(bore_zone.split_once('=').unwrap().1);
    let start = &projection[1];
    let values = records(&file);
    let zone = &projection[0];
    let definition = bore_zone.split_once('=').unwrap().0;
    let [owner, tolerance, form] = <[String; 3]>::try_from(references(&values[zone])).unwrap();
    let complex_zone =
        format!("{zone}=(SHAPE_ASPECT('','',{owner},.F.)TOLERANCE_ZONE(({tolerance}),{form}));");
    let complex_definition = format!(
        "{definition}=(PROJECTED_ZONE_DEFINITION({start},{})TOLERANCE_ZONE_DEFINITION({zone},()));",
        projection[2],
    );
    let complex = mutate(&file, &format!("{zone}={}", values[zone]), &complex_zone);
    let complex = mutate(&complex, bore_zone, &complex_definition);
    let (back, skipped) = read_back(&complex, &report);
    assert!(skipped.is_empty(), "{skipped:?}");
    assert!(matches!(
        kind_of(&back, "Bore projection"),
        ItemKind::Position(_)
    ));

    let additional = format!("#900000=TOLERANCE_ZONE_DEFINITION({zone},());\n");
    let mixed = mutate(
        &complex,
        "ENDSEC;\nEND-ISO",
        &format!("{additional}ENDSEC;\nEND-ISO"),
    );
    let (back, skipped) = read_back(&mixed, &report);
    assert!(!ids(&back).contains(&"Bore projection".to_owned()));
    assert!(reasons(&skipped).contains("zone definition uses TOLERANCE_ZONE_DEFINITION"));

    let circle = mutate(&file, bore_zone, "");
    let circle = mutate(
        &circle,
        &format!("{form}={}", values[&form]),
        &format!("{form}=TOLERANCE_ZONE_FORM('within a circle');"),
    );
    let (back, skipped) = read_back(&circle, &report);
    assert!(skipped.is_empty(), "{skipped:?}");
    assert!(
        matches!(kind_of(&back, "Bore projection"), ItemKind::Tolerance(tol) if tol.diameter_zone)
    );

    let usage = file
        .lines()
        .find(|line| {
            line.contains("=GEOMETRIC_ITEM_SPECIFIC_USAGE(")
                && references(line.split_once('=').unwrap().1).first() == Some(start)
        })
        .unwrap();
    let [aspect, shape, entry] =
        <[String; 3]>::try_from(references(usage.split_once('=').unwrap().1)).unwrap();
    let wrong = written_faces(&file, &report)
        .into_iter()
        .find_map(|(entity, face)| (face == f.flange_top).then_some(entity))
        .unwrap();
    assert_ne!(wrong, entry);
    let usage_id = usage.split_once('=').unwrap().0;
    let additional = format!(
        "{usage_id}=ITEM_IDENTIFIED_REPRESENTATION_USAGE('','',{aspect},{shape},SET_REPRESENTATION_ITEM(({entry},{wrong})));"
    );
    let broken = mutate(&file, usage, &additional);
    let (_, skipped) = read_back(&broken, &report);
    assert!(reasons(&skipped).contains("unused start face"));
    let malformed = format!(
        "{usage_id}=ITEM_IDENTIFIED_REPRESENTATION_USAGE('','',{aspect},{shape},SET_REPRESENTATION_ITEM(({entry},42)));"
    );
    let read = read_step_pmi(mutate(&file, usage, &malformed).as_bytes()).unwrap();
    assert!(read.sources.contains_key("Bore"));
    assert!(!read.sources.contains_key("Bore projection"));
    assert!(
        read.skipped
            .iter()
            .any(|skip| skip.reason.contains("identified face")),
        "{:?}",
        read.skipped
    );

    let target = &a1;
    let (property, _) = values
        .iter()
        .find(|(_, record)| {
            record.starts_with("PROPERTY_DEFINITION(") && references(record).contains(target)
        })
        .unwrap();
    let representation = references(
        values
            .values()
            .find(|record| {
                record.starts_with("SHAPE_DEFINITION_REPRESENTATION(")
                    && references(record).first() == Some(property)
            })
            .unwrap(),
    )[1]
    .clone();
    let target_rep = &values[&representation];
    let old_context = references(target_rep).pop().unwrap();
    let changed_rep = target_rep.replace(&format!("{old_context});"), "#900004);");
    let mixed_units = mutate(
        &file,
        &format!("{representation}={target_rep}"),
        &format!("{representation}={changed_rep}"),
    );
    let units = "#900001=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\n\
                 #900002=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#900001);\n\
                 #900003=(CONVERSION_BASED_UNIT('inch',#900002)LENGTH_UNIT()NAMED_UNIT(*));\n\
                 #900004=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#900003))REPRESENTATION_CONTEXT('',''));\n";
    let mixed_units = mutate(
        &mixed_units,
        "ENDSEC;\nEND-ISO",
        &format!("{units}ENDSEC;\nEND-ISO"),
    );
    let (_, skipped) = read_back(&mixed_units, &report);
    assert!(reasons(&skipped).contains("the target point is outside its face"));
    let valid = mutate(
        &mixed_units,
        "CARTESIAN_POINT('',(5.,5.,0.))",
        "CARTESIAN_POINT('',(0.2,0.2,0.))",
    );
    let (back, skipped) = read_back(&valid, &report);
    assert!(skipped.is_empty(), "{skipped:?}");
    let ItemKind::DatumTarget(target) = kind_of(&back, "A1") else {
        panic!("expected datum target");
    };
    assert_eq!(target.point_mm, v("5.08", "5.08", "0"));

    // A target moved into a hole is refused on the file's own face boundary.
    let broken = mutate(
        &file,
        "CARTESIAN_POINT('',(55.,36.,0.))",
        "CARTESIAN_POINT('',(50.,32.,0.))",
    );
    let (_, skipped) = read_back(&broken, &report);
    assert!(
        reasons(&skipped).contains("the target point is outside its face"),
        "{}",
        reasons(&skipped)
    );
}

#[test]
fn face_usages_accept_a_representation_wrapper() {
    let part = bracket();
    let faces = gdt_faces(&part);
    let (file, report) = write(&part, &gdt_pmi(&faces));
    let all = records(&file);
    let shape = line_with(&file, &["=ADVANCED_BREP_SHAPE_REPRESENTATION("]);
    let link = line_with(
        &file,
        &["=SHAPE_DEFINITION_REPRESENTATION(", &format!(",{shape});")],
    );
    let context = references(&all[&shape]).pop().unwrap();
    let item = &references(&all[&shape])[0];
    let wrapper = format!(
        "#900001=SHAPE_REPRESENTATION('',({item}),{context});\n\
         #900002=SHAPE_REPRESENTATION_RELATIONSHIP('','',{shape},#900001);\n"
    );
    let wrapped = mutate(
        &file,
        &format!("{link}={}", all[&link]),
        &format!("{link}={}", all[&link].replace(&shape, "#900001")),
    );
    let wrapped = mutate(
        &wrapped,
        "ENDSEC;\nEND-ISO",
        &format!("{wrapper}ENDSEC;\nEND-ISO"),
    );
    let (expected, _) = read_back(&file, &report);
    let (actual, skipped) = read_back(&wrapped, &report);
    assert!(skipped.is_empty(), "{skipped:?}");
    assert_eq!(actual, expected);
    let two_hops = mutate(
        &wrapped,
        &format!("SHAPE_REPRESENTATION_RELATIONSHIP('','',{shape},#900001);"),
        &format!(
            "SHAPE_REPRESENTATION_RELATIONSHIP('','',{shape},#900003);\n\
                  #900003=SHAPE_REPRESENTATION('',({item}),{context});\n\
                  #900004=SHAPE_REPRESENTATION_RELATIONSHIP('','',#900003,#900001);"
        ),
    );
    let read = read_step_pmi(two_hops.as_bytes()).unwrap();
    assert!(!read.sources.contains_key("A1"));
    assert!(
        read.skipped
            .iter()
            .any(|item| item.reason.contains("not of the aspect's part"))
    );
}

#[test]
fn unavailable_or_closed_projection_entries_are_refused() {
    let mut model = Model::new();
    let body = cylinder(&mut model, q("4"), q("5"), &Pose::identity()).unwrap();
    let part = Bracket { model, body };
    let barrel = part.round(2, "4", ["0", "0"]);
    let top = part.planes(2, "5")[0];
    let mut pmi = Pmi {
        items: vec![item(
            "Projection",
            ItemKind::Position(PositionControl {
                features: vec![barrel],
                frame: PositionFrame::Single {
                    segment: segment("0.05", None, Vec::new()),
                    projected_zone: Some(ProjectedZone {
                        height_mm: q("8"),
                        direction: v("0", "0", "1"),
                        entry_faces: vec![top],
                    }),
                },
            }),
        )],
        part: PartNotes::default(),
        views: Vec::new(),
    };
    let err = try_write(&part, &pmi).unwrap_err();
    assert!(err.contains("axis meets a closed end"), "{err}");
    let ItemKind::Position(control) = &mut pmi.items[0].kind else {
        panic!("expected a position control");
    };
    let PositionFrame::Single {
        projected_zone: Some(zone),
        ..
    } = &mut control.frame
    else {
        panic!("expected a projected zone");
    };
    zone.entry_faces[0] = top.moved(100_000);
    let err =
        catallax_step::pmi::check_on_model(&part.model, part.body, &pmi, &mut Decider::default())
            .unwrap_err();
    assert!(err.detail.contains("Projection"), "{}", err.detail);
    assert!(
        try_write(&part, &pmi)
            .unwrap_err()
            .contains("not in the file")
    );
}

#[test]
fn duplicate_target_number_invalidates_only_its_datum_and_dependents() {
    let part = bracket();
    let f = gdt_faces(&part);
    let (file, report) = write(&part, &gdt_pmi(&f));
    let a2 = file
        .lines()
        .find(|line| line.contains("PLACED_DATUM_TARGET_FEATURE('A2'"))
        .unwrap();
    let broken = mutate(&file, a2, &a2.replace(",'2');", ",'1');"));
    let (back, skipped) = read_back(&broken, &report);
    let mut names: Vec<_> = back.items.iter().map(|item| item.id.as_str()).collect();
    names.sort_unstable();
    assert_eq!(
        names,
        ["Bore", "Bore projection", "Datum B", "Datum C", "Datum D"]
    );
    assert!(
        skipped
            .iter()
            .any(|item| item.reason.contains("defined twice"))
    );
}

#[test]
fn straight_datum_target_can_lie_on_the_face_boundary() {
    let part = bracket();
    let f = gdt_faces(&part);
    let mut pmi = gdt_pmi(&f);
    for (centre, length) in [(v("30", "0", "0"), q("20")), (v("10", "10", "0"), q("6"))] {
        let ItemKind::DatumTarget(line) = &mut pmi.items[2].kind else {
            panic!("A2 should be the line target");
        };
        line.point_mm = centre;
        line.kind = TargetKind::Line {
            length_mm: length,
            direction: v("1", "0", "0"),
        };
        let (file, report) = write(&part, &pmi);
        let (back, skipped) = read_back(&file, &report);
        assert!(skipped.is_empty(), "{skipped:?}");
        assert_eq!(kind_of(&back, "A2"), kind_of(&pmi, "A2"));
    }
}

/// The C-2B fixture plus a knurl on the boss, a two-line note on the flange top
/// and a second note on the knurled boss itself.
fn text_pmi(part: &Bracket, f: &GdtFaces) -> Pmi<FaceId> {
    let boss = part.round(2, "8", ["30", "25"]);
    let mut pmi = gdt_pmi(f);
    pmi.items.extend([
        item(
            "Knurl",
            ItemKind::Knurl {
                designation: "RGE 0.8 DIN 82".into(),
                feature: boss.clone(),
            },
        ),
        item(
            "Flange note",
            ItemKind::FaceNote {
                text: "Deburr after knurling\nNo marking on this face".into(),
                feature: vec![f.flange_top],
            },
        ),
        item(
            "Boss note",
            ItemKind::FaceNote {
                text: "Protect the knurl during plating".into(),
                feature: boss,
            },
        ),
    ]);
    pmi
}

#[test]
fn knurl_and_face_notes_stay_apart_on_their_faces() {
    let part = bracket();
    let f = gdt_faces(&part);
    let boss = part.round(2, "8", ["30", "25"]);
    let pmi = text_pmi(&part, &f);
    pmi.validate().unwrap();
    let (file, report) = write(&part, &pmi);
    if let Ok(directory) = std::env::var("CATALLAX_PMI_OUT") {
        std::fs::write(
            std::path::Path::new(&directory).join("c2n-knurl-notes.step"),
            &file,
        )
        .unwrap();
    }
    let (back, skipped) = read_back(&file, &report);
    assert_eq!(skipped, Vec::new());

    // The knurl and the note on the same boss stay two items of their own kinds.
    let mut sorted_boss = boss.clone();
    sorted_boss.sort();
    match kind_of(&back, "Knurl") {
        ItemKind::Knurl {
            designation,
            feature,
        } => {
            let mut feature = feature.clone();
            feature.sort();
            assert_eq!(
                (designation.as_str(), feature),
                ("RGE 0.8 DIN 82", sorted_boss.clone())
            );
        }
        other => panic!("the knurl reads back as {other:?}"),
    }
    match kind_of(&back, "Boss note") {
        ItemKind::FaceNote { text, feature } => {
            let mut feature = feature.clone();
            feature.sort();
            assert_eq!(
                (text.as_str(), feature),
                ("Protect the knurl during plating", sorted_boss)
            );
        }
        other => panic!("the boss note reads back as {other:?}"),
    }
    assert_eq!(
        kind_of(&back, "Flange note"),
        &ItemKind::FaceNote {
            text: "Deburr after knurling\nNo marking on this face".into(),
            feature: vec![f.flange_top],
        }
    );
    let mut expected = pmi.clone();
    for item in &mut expected.items {
        if let ItemKind::Size(size) = &mut item.kind {
            size.nominal = Some(q("8"));
        }
    }
    assert_same(back, expected);

    // Text is counted apart from semantic PMI, and none of it is drawn: the
    // carrier is editable text, not a graphic callout or a tolerance.
    let counts: Vec<String> = validation_items(&file)
        .into_values()
        .flatten()
        .map(|(_, fields)| fields)
        .collect();
    for count in [
        "'number of semantic pmi elements',12.",
        "'number of annotations',0.",
        "'number of semantic text notes',3.",
        "'number of semantic text notes on part',0.",
        "'number of semantic text notes on geometry',3.",
        "'number of semantic text notes on PMI',0.",
    ] {
        assert!(counts.iter().any(|fields| fields == count), "{count}");
    }
    assert_eq!(
        file.matches("PROPERTY_DEFINITION('semantic text'").count(),
        3
    );
    assert_eq!(
        file.matches("DESCRIPTIVE_REPRESENTATION_ITEM('knurl designation','RGE 0.8 DIN 82')")
            .count(),
        1
    );
    assert!(!file.contains("DRAUGHTING_CALLOUT"));

    // The P1 bracket, which has no text, writes no text counts.
    let (p1, _) = bracket_pmi(&part);
    let (p1, _) = write(&part, &p1);
    assert!(!p1.contains("semantic text"));
}

#[test]
fn bad_text_is_refused_and_unknown_text_is_reported() {
    let part = bracket();
    let f = gdt_faces(&part);
    for (kind, reason) in [
        (
            ItemKind::FaceNote {
                text: " ".into(),
                feature: vec![f.flange_top],
            },
            "a face note needs text",
        ),
        (
            ItemKind::Knurl {
                designation: " ".into(),
                feature: vec![f.flange_top],
            },
            "a knurl designation needs text",
        ),
        (
            ItemKind::Knurl {
                designation: "RGE 0.8 DIN 82".into(),
                feature: Vec::new(),
            },
            "a feature has no face",
        ),
    ] {
        let mut pmi = gdt_pmi(&f);
        pmi.items.push(item("Text", kind));
        let detail = try_write(&part, &pmi).expect_err(reason);
        assert!(detail.contains(reason), "{detail} lacks {reason}");
    }

    let mut text = text_pmi(&part, &f);
    let note = "One\n\nThree, café 😀";
    let designation = format!("RGE 'Ø' \\ DIN 82 {}", "x".repeat(60));
    for item in &mut text.items {
        match &mut item.kind {
            ItemKind::FaceNote { text, .. } if item.id == "Flange note" => *text = note.into(),
            ItemKind::Knurl {
                designation: value, ..
            } => *value = designation.clone(),
            _ => {}
        }
    }
    let (file, report) = write(&part, &text);
    let legacy = mutate(&file, r"\X4\0001F600\X0\", r"\X2\D83DDE00\X0\");
    let (legacy, skipped) = read_back(&legacy, &report);
    assert!(skipped.is_empty(), "{skipped:?}");
    assert_eq!(
        kind_of(&legacy, "Flange note"),
        kind_of(&text, "Flange note")
    );
    let unsupported = mutate(&file, r"\X4\0001F600\X0\", r"\PE\\S\*\S\U\S\b");
    let unsupported = read_step_pmi(unsupported.as_bytes()).unwrap();
    assert!(!unsupported.sources.contains_key("Flange note"));
    assert!(
        unsupported
            .skipped
            .iter()
            .any(|item| item.reason.contains("unsupported STEP encoding"))
    );

    let (back, skipped) = read_back(&file, &report);
    assert!(skipped.is_empty(), "{skipped:?}");
    let (back, text) = (sorted(back), sorted(text));
    for name in ["Flange note", "Knurl"] {
        assert_eq!(kind_of(&back, name), kind_of(&text, name));
    }

    let records = records(&file);
    let source = read_step_pmi(file.as_bytes()).unwrap().sources["Flange note"].clone();
    let link = records
        .values()
        .find(|body| {
            body.starts_with("PROPERTY_DEFINITION_REPRESENTATION(") && references(body)[0] == source
        })
        .unwrap();
    let rep = &references(link)[1];
    let body = &records[rep];
    let value = &references(body)[0];
    let unordered = mutate(
        &file,
        &format!("{rep}={body}"),
        &format!(
            "{rep}={}",
            body.replace(&format!("({value})"), &format!("({value},#900000)"))
        ),
    );
    let unordered = mutate(
        &unordered,
        "ENDSEC;\nEND-ISO",
        "#900000=DESCRIPTIVE_REPRESENTATION_ITEM('','Append after plating');\nENDSEC;\nEND-ISO",
    );
    let read = read_step_pmi(unordered.as_bytes()).unwrap();
    assert!(!read.sources.contains_key("Flange note"));
    assert!(
        read.skipped
            .iter()
            .any(|item| item.reason.contains("text has no defined order"))
    );

    let original = read_step_pmi(file.as_bytes()).unwrap();
    let sources = &original.sources;
    for (damaged, independent) in [("Flange note", "Knurl"), ("Knurl", "Flange note")] {
        let link = records
            .values()
            .find(|body| {
                body.starts_with("PROPERTY_DEFINITION_REPRESENTATION(")
                    && references(body)[0] == sources[damaged]
            })
            .unwrap();
        let rep = &references(link)[1];
        let body = &records[rep];
        let value = &references(body)[0];
        for bad in [
            body.replace(&format!("({value})"), &format!("({value},$)")),
            format!("GEOMETRIC_CURVE_SET('',({value}));"),
        ] {
            let changed = mutate(&file, &format!("{rep}={body}"), &format!("{rep}={bad}"));
            let read = read_step_pmi(changed.as_bytes()).unwrap();
            assert!(!read.sources.contains_key(damaged));
            assert_eq!(
                kind_of(&read.pmi, independent),
                kind_of(&original.pmi, independent)
            );
            assert!(
                read.skipped
                    .iter()
                    .any(|skip| skip.entity == sources[damaged])
            );
        }
    }

    // A text variable this reader does not know is reported, and the rest stays.
    let (file, report) = write(&part, &text_pmi(&part, &f));
    let broken = mutate(&file, "'knurl designation'", "'surface treatment'");
    let (back, skipped) = read_back(&broken, &report);
    assert_eq!(back.items.len(), 14);
    assert!(
        skipped
            .iter()
            .any(|skip| skip.reason == "a text variable `surface treatment` is not read"),
        "{skipped:?}"
    );
}

fn camera(
    eye: [&str; 3],
    direction: [&str; 3],
    up: [&str; 3],
    distance: &str,
    window: [&str; 2],
    perspective: bool,
) -> Camera {
    Camera {
        eye: eye.map(q),
        direction: direction.map(q),
        up: up.map(q),
        distance: q(distance),
        window: window.map(q),
        perspective,
    }
}

/// The P1 bracket with two saved views: the bore from the front, and the
/// datum frame and hole pattern from above in perspective. The bore position
/// is in both.
fn viewed_pmi(part: &Bracket) -> (Pmi<FaceId>, BTreeMap<&'static str, &'static str>) {
    let (mut pmi, nominal) = bracket_pmi(part);
    let names = |names: &[&str]| names.iter().map(|name| (*name).to_owned()).collect();
    pmi.views = vec![
        SavedView {
            name: "Bore".into(),
            camera: camera(
                ["30", "-120", "10"],
                ["0", "1", "0"],
                ["0", "0", "1"],
                "145",
                ["100", "60"],
                false,
            ),
            items: names(&["Datum A", "Bore diameter", "Bore finish", "Bore position"]),
        },
        SavedView {
            name: "Top".into(),
            camera: camera(
                ["30", "20", "150"],
                ["0", "0", "-1"],
                ["0", "1", "0"],
                "150",
                ["90", "60"],
                true,
            ),
            items: names(&[
                "Datum B",
                "Datum C",
                "Hole position",
                "Bore position",
                "Bore finish",
            ]),
        },
    ];
    (pmi, nominal)
}

/// The file's entities by number: the full record text.
fn records(file: &str) -> BTreeMap<String, String> {
    file.lines()
        .filter_map(|line| {
            let (id, rest) = line.split_once('=')?;
            id.starts_with('#')
                .then(|| (id.to_owned(), rest.to_owned()))
        })
        .collect()
}

fn callout_associations(file: &BTreeMap<String, String>, item: &str) -> Vec<(String, String)> {
    file.iter()
        .filter(|(_, record)| record.starts_with("DRAUGHTING_MODEL_ITEM_ASSOCIATION("))
        .filter(|(_, record)| {
            let callout = references(record).pop().unwrap();
            file[&callout].starts_with(&format!("DRAUGHTING_CALLOUT('{item}'"))
        })
        .map(|(id, record)| (id.clone(), record.clone()))
        .collect()
}

fn check_edition_one_draughting_models(file: &str) {
    let entities = entities(file);
    for (id, (kind, fields)) in &entities {
        if kind == "CAMERA_MODEL_D3" {
            assert!(
                entities.values().any(|(kind, fields)| {
                    kind == "ITEM_DEFINED_TRANSFORMATION" && references(fields).first() == Some(id)
                }),
                "camera_model.WR1: {id}"
            );
        }
        if kind == "REPRESENTATION_MAP" {
            let [origin, representation] = <[String; 2]>::try_from(references(fields)).unwrap();
            assert!(
                references(&entities[&representation].1).contains(&origin),
                "representation_map.WR1: {id}"
            );
        }
    }
    let mut names = std::collections::BTreeSet::new();
    for (kind, fields) in entities.values() {
        if kind != "DRAUGHTING_MODEL" {
            continue;
        }
        let (name, _) = fields.split_once(",(").unwrap();
        assert!(names.insert(name), "DRAUGHTING_MODEL.UR1: {name}");
        for member in references(fields).into_iter().rev().skip(1) {
            let (kind, attrs) = &entities[&member];
            if kind != "STYLED_ITEM" {
                continue;
            }
            let links = references(attrs);
            let (target, styles) = links.split_last().unwrap();
            if entities[target].0 != "MAPPED_ITEM" {
                continue;
            }
            let mapping = &references(&entities[target].1)[0];
            let representation = &references(&entities[mapping].1)[1];
            assert!(entities[representation].0.ends_with("SHAPE_REPRESENTATION"));
            let rejected_assignments = styles
                .iter()
                .filter(|style| {
                    let (_, choices) = &entities[*style];
                    let non_curve = references(choices)
                        .iter()
                        .filter(|id| entities[*id].0 != "CURVE_STYLE")
                        .count()
                        + usize::from(choices.contains("NULL_STYLE"));
                    non_curve != 1
                })
                .count();
            assert_eq!(rejected_assignments, 1, "DRAUGHTING_MODEL.WR2: {member}");
        }
    }
}

#[test]
fn p1_position_patch_refuses_a_dimension_linked_callout() {
    use catallax_import::{PmiEdit, patch_step_pmi};
    let part = bracket();
    let (mut pmi, _) = viewed_pmi(&part);
    pmi.items
        .retain(|item| !matches!(&item.kind, ItemKind::Texture(_)));
    for view in &mut pmi.views {
        view.items.retain(|item| item == "Bore position");
    }
    let (mut file, _) = write(&part, &pmi);
    let all = records(&file);
    let associations = callout_associations(&all, "Bore position");
    assert_eq!(associations.len(), 2);
    let redirect_size = |source: &str, root: &str, size: &str| {
        let all = records(source);
        let record = &all[root];
        let attributes = record
            .split_once("GEOMETRIC_TOLERANCE(")
            .map_or(record.as_str(), |(_, fields)| {
                fields.split_once(')').unwrap().0
            });
        let target = &references(attributes)[1];
        mutate(
            source,
            &format!("{root}={record}"),
            &format!("{root}={}", record.replace(target, size)),
        )
    };
    let root = &references(&associations[0].1)[0];
    let original = read_step_pmi(file.as_bytes()).unwrap();
    let size = &original.sources["Bore diameter"];
    file = redirect_size(&file, root, size);
    for (id, association) in &associations {
        let redirected = association.replace(&format!(",{root},"), &format!(",{size},"));
        file = mutate(
            &file,
            &format!("{id}={association}"),
            &format!("{id}={redirected}"),
        );
    }
    let read = read_step_pmi(file.as_bytes()).unwrap();
    assert_eq!(&read.sources["Bore position"], root);
    assert!(
        matches!(kind_of(&read.pmi, "Bore position"), ItemKind::Tolerance(item) if item.value == q("0.05"))
    );
    let edit = PmiEdit::ToleranceValue {
        source: root.clone(),
        value: q("0.1"),
    };
    let error = patch_step_pmi(file.as_bytes(), std::slice::from_ref(&edit)).unwrap_err();
    assert!(
        associations.iter().any(|(id, _)| {
            error
                .detail
                .contains(&format!("{id} (DRAUGHTING_MODEL_ITEM_ASSOCIATION)"))
        }),
        "{error:?}"
    );
    pmi.views.clear();
    let (unpresented, _) = write(&part, &pmi);
    let plain = read_step_pmi(unpresented.as_bytes()).unwrap();
    let source = &plain.sources["Bore position"];
    let unpresented = redirect_size(&unpresented, source, &plain.sources["Bore diameter"]);
    let patched = patch_step_pmi(
        unpresented.as_bytes(),
        &[PmiEdit::ToleranceValue {
            source: source.clone(),
            value: q("0.1"),
        }],
    )
    .unwrap();
    let read = read_step_pmi(&patched).unwrap();
    assert!(
        matches!(kind_of(&read.pmi, "Bore position"), ItemKind::Tolerance(item) if item.value == q("0.1"))
    );
    assert_eq!(
        kind_of(&read.pmi, "Bore diameter"),
        kind_of(&original.pmi, "Bore diameter")
    );
}

#[test]
fn a_wrong_plane_cannot_pair_a_shared_callout() {
    let part = bracket();
    let (pmi, _) = viewed_pmi(&part);
    let (file, _) = write(&part, &pmi);
    let all = records(&file);
    let top = line_with(&file, &["=DRAUGHTING_MODEL('Top',"]);
    let members = references(&all[&top]);
    let camera = members
        .iter()
        .find(|id| all[*id].starts_with("CAMERA_MODEL_D3("))
        .unwrap();
    let context = members.last().unwrap();
    let (plane, callout) = members
        .iter()
        .filter(|id| all[*id].starts_with("ANNOTATION_PLANE("))
        .find_map(|id| {
            let callout = references(&all[id]).pop()?;
            all[&callout]
                .starts_with("DRAUGHTING_CALLOUT('Bore position',")
                .then_some((id, callout))
        })
        .unwrap();
    let style = &references(&all[plane])[0];
    let datum = &references(&callout_associations(&all, "Datum A")[0].1)[0];
    let extra = format!(
        "#899901=CARTESIAN_POINT('',(0.,0.,0.));\n\
         #899902=DIRECTION('',(1.,0.,0.));\n\
         #899903=AXIS2_PLACEMENT_3D('',#899901,#899902,$);\n\
         #899904=PLANE('',#899903);\n\
         #899905=ANNOTATION_PLANE('',({style}),#899904,({callout}));\n\
         #899906=DRAUGHTING_MODEL('Side',(#899905,#899908),{context});\n\
         #899907=DRAUGHTING_MODEL_ITEM_ASSOCIATION('','',{datum},#899906,{callout});\n\
         #899908={camera};\n",
        camera = all[camera].replace("'Top'", "'Side'"),
    );
    let changed = mutate(
        &file,
        "ENDSEC;\nEND-ISO",
        &format!("{extra}ENDSEC;\nEND-ISO"),
    );
    let read = read_step_pmi(changed.as_bytes()).unwrap();
    assert_eq!(read.pmi.views[1], pmi.views[1]);
    assert!(
        read.skipped
            .iter()
            .any(|item| item.entity == callout && item.reason.contains("vary in depth"))
    );
    let valid_side = mutate(
        &changed,
        "#899902=DIRECTION('',(1.,0.,0.));",
        "#899902=DIRECTION('',(0.,0.,1.));",
    );
    let read = read_step_pmi(valid_side.as_bytes()).unwrap();
    let top = read
        .pmi
        .views
        .iter()
        .find(|view| view.name == "Top")
        .unwrap();
    assert_eq!(top, &pmi.views[1]);
    let side = read
        .pmi
        .views
        .iter()
        .find(|view| view.name == "Side")
        .unwrap();
    assert_eq!(side.items, ["Datum A"]);
}

#[test]
fn saved_views_write_their_callouts_and_read_back() {
    let part = bracket();
    let (pmi, nominal) = viewed_pmi(&part);
    pmi.validate().unwrap();
    let (file, report) = write(&part, &pmi);
    check_edition_one_draughting_models(&file);
    let directory = std::env::var_os("CATALLAX_PMI_OUT").map(std::path::PathBuf::from);
    if let Some(directory) = &directory {
        std::fs::write(directory.join("c4c-saved-views.step"), &file).unwrap();
    }
    for (view, filename) in pmi
        .views
        .iter()
        .zip(["c4c-bore-view.step", "c4c-top-view.step"])
    {
        let mut one_view = pmi.clone();
        one_view.views = vec![view.clone()];
        let (fixture, _) = write(&part, &one_view);
        check_edition_one_draughting_models(&fixture);
        let read = read_step_pmi(fixture.as_bytes()).unwrap();
        assert_eq!(read.pmi.views.as_slice(), std::slice::from_ref(view));
        if let Some(directory) = &directory {
            std::fs::write(directory.join(filename), fixture).unwrap();
        }
    }
    let read = read_step_pmi(file.as_bytes()).unwrap();
    assert_eq!(read.skipped, Vec::new());
    // Every callout presents an item that was read, so none is graphics.
    assert_eq!(read.graphics, Vec::new());
    assert_eq!(read.pmi.views, pmi.views);
    let faces = written_faces(&file, &report);
    let back = read
        .pmi
        .try_map(|_, key| faces.get(key).copied().ok_or_else(|| key.clone()))
        .unwrap();
    let mut expected = pmi.clone();
    for item in &mut expected.items {
        let value = nominal.get(item.id.as_str()).map(|text| q(text));
        match &mut item.kind {
            ItemKind::Size(size) => size.nominal = value,
            ItemKind::Location(location) => location.nominal = value,
            _ => {}
        }
    }
    assert_eq!(back.items.len(), pmi.items.len());
    assert_same(back, expected);

    let all = records(&file);
    assert_eq!(file.matches("=DRAUGHTING_CALLOUT(").count(), 9);
    let counts: Vec<String> = validation_items(&file)
        .into_values()
        .flatten()
        .map(|(_, fields)| fields)
        .collect();
    assert!(
        counts
            .iter()
            .any(|fields| fields == "'number of annotations',9.")
    );
    for (item, entity, count) in [
        ("Datum A", "DATUM_FEATURE(", 1),
        ("Bore diameter", "DIMENSIONAL_SIZE(", 1),
        ("Bore position", "POSITION_TOLERANCE(", 2),
        ("Hole position", "POSITION_TOLERANCE(", 1),
    ] {
        let associations = callout_associations(&all, item);
        assert_eq!(associations.len(), count, "{item}");
        for (_, association) in associations {
            assert!(association.starts_with(
                "DRAUGHTING_MODEL_ITEM_ASSOCIATION('PMI representation to presentation link'"
            ));
            let definition = references(&association)[0].clone();
            assert!(
                all[&definition].contains(entity),
                "{item}: {}",
                all[&definition]
            );
        }
    }
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
    let cross = |a: &[f64], b: &[f64]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let unit = |axis: Vec<f64>| {
        let norm = dot(&axis, &axis).sqrt();
        axis.into_iter().map(|x| x / norm).collect::<Vec<_>>()
    };
    let style = |entity: &str, expected: &str| {
        let assignment = references(entity)[0].clone();
        let style = references(&all[&assignment])[0].clone();
        assert!(all[&style].starts_with(expected), "{}", all[&style]);
    };
    assert_eq!(
        all.values()
            .filter(|body| body.starts_with("ANNOTATION_PLANE("))
            .count(),
        9
    );
    for view in ["Bore", "Top"] {
        let model = all
            .values()
            .find(|body| body.starts_with(&format!("DRAUGHTING_MODEL('{view}'")))
            .unwrap();
        let members = references(model);
        let camera = members
            .iter()
            .find(|id| all[*id].starts_with("CAMERA_MODEL_D3("))
            .unwrap();
        let camera = references(&all[camera]);
        let axes = references(&all[&camera[0]]);
        let eye = numbers(&all[&axes[0]]);
        let depth = unit(numbers(&all[&axes[1]]));
        let across = [1.0, 0.0, 0.0];
        let up = if view == "Bore" {
            [0.0, 0.0, 1.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let camera_x = unit(numbers(&all[&axes[2]]));
        assert!(dot(&cross(&depth, &camera_x), &up) > 1.0 - 1e-9);
        let volume = &all[&camera[1]];
        let perspective = volume.contains(".CENTRAL.");
        let distance = volume.split(',').nth(2).unwrap().parse::<f64>().unwrap();
        let window = references(volume).pop().unwrap();
        let fields: Vec<_> = all[&window].split(',').collect();
        let half = [
            fields[1].parse::<f64>().unwrap() / 2.0,
            fields[2].parse::<f64>().unwrap() / 2.0,
        ];
        let window_placement = &all[&references(&all[&window])[0]];
        assert!(window_placement.starts_with("AXIS2_PLACEMENT_3D("));
        assert_eq!(
            numbers(&all[&references(window_placement)[0]]),
            [-half[0], -half[1], distance]
        );
        let mapped: Vec<_> = members
            .iter()
            .filter(|id| all[*id].starts_with("MAPPED_ITEM("))
            .collect();
        assert_eq!(mapped.len(), 1);
        let mapped_style = members
            .iter()
            .find(|id| {
                all[*id].starts_with("STYLED_ITEM(")
                    && references(&all[*id]).last() == Some(mapped[0])
            })
            .unwrap();
        let assignment = &references(&all[mapped_style])[0];
        let curve = &references(&all[assignment])[0];
        let colour = &references(&all[curve])[0];
        assert_eq!(all[curve], format!("CURVE_STYLE('',$,$,{colour});"));
        assert_eq!(all[colour], "COLOUR();");
        let planes: Vec<_> = members
            .into_iter()
            .filter(|id| all[id].starts_with("ANNOTATION_PLANE("))
            .collect();
        assert_eq!(
            planes.len(),
            pmi.views
                .iter()
                .find(|entry| entry.name == view)
                .unwrap()
                .items
                .len()
        );
        let mut preceding_bottom = half[1];
        for id in planes {
            style(&all[&id], "FILL_AREA_STYLE(");
            let keys = references(&all[&id]);
            let surface = &all[&keys[1]];
            let frame = references(&all[&references(surface)[0]]);
            let origin = numbers(&all[&frame[0]]);
            let normal = unit(numbers(&all[&frame[1]]));
            assert!(dot(&normal, &depth) < -1.0 + 1e-9);
            let plane_right = unit(numbers(&all[&frame[2]]));
            assert!(dot(&plane_right, &across) > 1.0 - 1e-9);
            assert!(dot(&cross(&normal, &plane_right), &up) > 1.0 - 1e-9);
            let callout = &keys[2];
            let name = all[callout].split('\'').nth(1).unwrap();
            let occurrence = references(&all[callout]).pop().unwrap();
            let set = references(&all[&occurrence]).pop().unwrap();
            let lines = references(&all[&set]);
            let leader = references(&all[&lines[0]]);
            assert_eq!(numbers(&all[&leader[0]]), origin);
            let squared_radius = |x: f64, y: f64| (origin[0] - x).powi(2) + (origin[1] - y).powi(2);
            match name {
                "Datum A" => assert_eq!(origin[2], 0.0),
                "Datum B" => assert_eq!(origin[1], 0.0),
                "Datum C" => assert_eq!(origin[0], 0.0),
                "Bore diameter" | "Bore finish" | "Bore position" => {
                    assert!((squared_radius(30.0, 25.0) - 16.0).abs() < 1e-9)
                }
                "Hole position" => assert!(
                    [(10.0, 12.0), (10.0, 32.0), (50.0, 12.0), (50.0, 32.0)]
                        .into_iter()
                        .any(|(x, y)| (squared_radius(x, y) - 4.0).abs() < 1e-9)
                ),
                _ => panic!("unexpected saved-view label {name}"),
            }
            let wire: Vec<_> = lines[1..].iter().map(|id| references(&all[id])).collect();
            let symbol = &wire[0];
            let junction = usize::from(name == "Bore finish");
            assert_eq!(numbers(&all[&leader[1]]), numbers(&all[&symbol[junction]]));
            let strokes_in_plane: Vec<Vec<[f64; 2]>> = wire
                .iter()
                .map(|path| {
                    path.iter()
                        .map(|id| {
                            let xyz = numbers(&all[id]);
                            [dot(&xyz, &across), dot(&xyz, &up)]
                        })
                        .collect()
                })
                .collect();
            let (decoded, step) = pmi_strokes::label(&strokes_in_plane, name == "Bore finish");
            let expected = match name {
                "Datum A" => "A",
                "Datum B" => "B",
                "Datum C" => "C",
                "Bore diameter" => "⌀8 H7",
                "Bore finish" => "Ra 1.6",
                "Bore position" => "POS|⌀0.05|A|B|C",
                "Hole position" => "POS|⌀0.2|A|B|C",
                _ => panic!("unexpected saved-view label {name}"),
            };
            assert_eq!(decoded, expected, "{view} {name}");
            if name == "Bore diameter" {
                let mut altered = strokes_in_plane;
                let [x, y] = altered[4][0];
                altered[4] = vec![[x + step, y], [x + 3.0 * step, y]];
                altered[5] = vec![[x + 2.0 * step, y], [x + 2.0 * step, y + 6.0 * step]];
                altered[6] = vec![[x + step, y + 6.0 * step], [x + 3.0 * step, y + 6.0 * step]];
                assert_eq!(pmi_strokes::label(&altered, false).0, "⌀8 I7");
            }
            let expected_faces: std::collections::BTreeSet<_> = pmi
                .items
                .iter()
                .find(|item| item.id == name)
                .unwrap()
                .kind
                .faces()
                .into_iter()
                .copied()
                .collect();
            let definitions: Vec<_> = callout_associations(&all, name)
                .into_iter()
                .filter_map(|(_, link)| {
                    let ids = references(&link);
                    (ids.last() == Some(callout)).then(|| ids[0].clone())
                })
                .collect();
            let mut pending = definitions;
            let mut actual_faces = std::collections::BTreeSet::new();
            while let Some(id) = pending.pop() {
                let record = &all[&id];
                if record.contains("POSITION_TOLERANCE(") || record.starts_with("DIMENSIONAL_SIZE(")
                {
                    pending.push(
                        references(record)[if record.contains("POSITION_TOLERANCE(") {
                            1
                        } else {
                            0
                        }]
                        .clone(),
                    );
                }
                for usage in all.values() {
                    let refs = references(usage);
                    if usage.contains("_USAGE(") && refs.first() == Some(&id) {
                        actual_faces
                            .extend(refs[2..].iter().filter_map(|id| faces.get(id).copied()));
                    } else if usage.starts_with("SHAPE_ASPECT_RELATIONSHIP(")
                        && refs.first() == Some(&id)
                        && all[&refs[1]].starts_with("SHAPE_ASPECT(")
                    {
                        pending.push(refs[1].clone());
                    }
                }
            }
            assert_eq!(actual_faces, expected_faces, "{view} {name}");
            let project = |id: String| {
                let world = numbers(&all[&id]);
                let from_eye: Vec<_> = world.iter().zip(&eye).map(|(p, eye)| p - eye).collect();
                let from_plane: Vec<_> = world
                    .iter()
                    .zip(&origin)
                    .map(|(p, plane)| p - plane)
                    .collect();
                assert!(dot(&from_plane, &normal).abs() < 1e-9);
                let scale = if perspective {
                    distance / dot(&from_eye, &depth)
                } else {
                    1.0
                };
                [dot(&from_eye, &across) * scale, dot(&from_eye, &up) * scale]
            };
            let strokes: Vec<_> = lines
                .iter()
                .skip(1)
                .flat_map(|line| references(&all[line]))
                .map(project)
                .collect();
            assert!(
                strokes[0][0] < strokes[1][0],
                "{view}: the first stroke is mirrored"
            );
            for point in &strokes {
                assert!(
                    point
                        .iter()
                        .zip(half)
                        .all(|(value, bound)| value.abs() < bound),
                    "{view}: {point:?}"
                );
            }
            for axis in 0..2 {
                let (lo, hi) = strokes
                    .iter()
                    .map(|p| p[axis])
                    .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
                        (lo.min(x), hi.max(x))
                    });
                assert!(hi - lo > 2.0, "{view}: callout collapsed along {axis}");
                if axis == 1 {
                    assert!(hi < preceding_bottom);
                    preceding_bottom = lo;
                }
            }
        }
    }
    for (id, expected) in [
        ("Datum A", r"DTM\wA"),
        ("Bore diameter", r"DIM\w⌀8 H7"),
        ("Bore position", r"FCF\w⌖\w⌀0.05\wA\wB\wC"),
    ] {
        let association = callout_associations(&all, id).pop().unwrap().1;
        let callout = references(&association).pop().unwrap();
        let within = all
            .iter()
            .find(|(_, body)| {
                body.starts_with("CHARACTERIZED_ITEM_WITHIN_REPRESENTATION(")
                    && references(body)[0] == callout
            })
            .unwrap()
            .0;
        let properties = validation_items(&file);
        let actual = &properties[within][0].1;
        let expected = expected.replace('\\', r"\\");
        let expected = expected
            .replace('⌀', r"\X2\2300\X0\")
            .replace('⌖', r"\X2\2316\X0\");
        assert_eq!(*actual, format!("'equivalent unicode string','{expected}'"));
    }
    // Two views: a model named for each, its camera and its relation to the
    // global model.
    assert_eq!(file.matches("=CAMERA_MODEL_D3(").count(), 2);
    assert_eq!(
        file.matches("=MECHANICAL_DESIGN_AND_DRAUGHTING_RELATIONSHIP(")
            .count(),
        2
    );
    assert!(file.contains("=DRAUGHTING_MODEL('Top',"));
    assert!(file.contains("=VIEW_VOLUME(.CENTRAL.,"));

    // A view that names an item with no callout, or a camera looking along
    // its own up direction, is refused.
    for unsupported in ["Slot width", "Slot end radius", "Flatness of A"] {
        let mut bad = pmi.clone();
        bad.views[0].items.push(unsupported.into());
        let reason = try_write(&part, &bad).unwrap_err();
        assert!(
            reason.contains(&format!("item {unsupported} has no callout")),
            "{reason}"
        );
    }
    let mut bad = pmi.clone();
    bad.views[1].camera.up = v("0", "0", "1");
    assert!(
        try_write(&part, &bad)
            .unwrap_err()
            .contains("up direction square to it")
    );
}

#[test]
fn graphical_only_callouts_are_read_apart_from_semantic_items() {
    let part = bracket();
    let (pmi, _) = viewed_pmi(&part);
    let (file, report) = write(&part, &pmi);
    let all = records(&file);
    let association = |item| callout_associations(&all, item).pop().unwrap();

    let (_, size_link) = association("Bore diameter");
    let (datum_link_id, datum_link) = association("Datum A");
    let size_callout = references(&size_link).pop().unwrap();
    let size_occurrence = references(&all[&size_callout]).pop().unwrap();
    let curve_set = references(&all[&size_occurrence]).pop().unwrap();
    let stroke = references(&all[&curve_set])[0].clone();
    let malformed = mutate(
        &file,
        &format!("{stroke}=POLYLINE("),
        &format!("{stroke}=CIRCLE("),
    );
    let read = read_step_pmi(malformed.as_bytes()).unwrap();
    assert!(read.sources.contains_key("Bore diameter"));
    assert!(
        !read.pmi.views[0]
            .items
            .iter()
            .any(|id| id == "Bore diameter")
    );
    assert!(
        read.skipped
            .iter()
            .any(|item| item.entity == size_callout && item.reason.contains("CIRCLE"))
    );

    let view_model = all
        .iter()
        .find(|(_, record)| record.starts_with("DRAUGHTING_MODEL('Bore',"))
        .unwrap()
        .0
        .clone();
    let mapped = references(&all[&view_model])
        .into_iter()
        .find(|member| all[member].starts_with("MAPPED_ITEM("))
        .unwrap();
    let map_target = references(&all[&mapped]).pop().unwrap();
    let target_point = references(&all[&map_target])[0].clone();
    let shifted = mutate(
        &file,
        &format!("{target_point}={}", all[&target_point]),
        &format!("{target_point}=CARTESIAN_POINT('',(100.,0.,0.));"),
    );
    let read = read_step_pmi(shifted.as_bytes()).unwrap();
    assert_eq!(
        read.pmi
            .views
            .iter()
            .map(|view| view.name.as_str())
            .collect::<Vec<_>>(),
        ["Top"]
    );
    assert!(
        read.skipped
            .iter()
            .any(|item| item.entity == view_model && item.reason.contains("nonidentity"))
    );
    // The refused view's mapping is checked at that view alone: its callout
    // is not refused with it, and no other view shows the size.
    assert!(
        !read.skipped.iter().any(|item| item.entity == size_callout),
        "{:?}",
        read.skipped
    );
    assert!(
        read.pmi
            .views
            .iter()
            .all(|view| !view.items.iter().any(|id| id == "Bore diameter"))
    );

    let datum_callout = references(&datum_link).pop().unwrap();
    let new_link = datum_link.replace(&datum_callout, &size_occurrence);
    let shared = mutate(
        &file,
        &format!("{datum_link_id}={datum_link}"),
        &format!("{datum_link_id}={new_link}"),
    );
    let read = read_step_pmi(shared.as_bytes()).unwrap();
    for name in ["Datum A", "Bore diameter"] {
        assert_eq!(
            read.pmi.views[0]
                .items
                .iter()
                .filter(|id| *id == name)
                .count(),
            1
        );
    }

    let finish = &read_step_pmi(file.as_bytes()).unwrap().sources["Bore finish"];
    let link = callout_associations(&all, "Bore finish")
        .into_iter()
        .find(|(_, link)| references(link).first() == Some(finish))
        .unwrap();
    let preserved = mutate(&file, &format!("{}={}\n", link.0, link.1), "");
    let read = read_step_pmi(preserved.as_bytes()).unwrap();
    assert!(
        read.graphics
            .iter()
            .any(|graphic| graphic.name == "Bore finish" && graphic.faces.is_some())
    );
    assert!(
        !read.pmi.views[0]
            .items
            .iter()
            .any(|item| item == "Bore finish")
    );

    // Without its association the datum callout presents nothing read: it is
    // graphics in the Bore view, and the datum is still one semantic item.
    let (id, record) = association("Datum A");
    let broken = mutate(&file, &format!("{id}={record}\n"), "");
    let (back, skipped) = read_back(&broken, &report);
    assert_eq!(skipped, Vec::new());
    assert_eq!(back.items.len(), pmi.items.len());
    let read = read_step_pmi(broken.as_bytes()).unwrap();
    let [graphic] = read.graphics.as_slice() else {
        panic!("{:?}", read.graphics);
    };
    assert_eq!(
        (
            graphic.name.as_str(),
            graphic.kind.as_str(),
            graphic.text.as_deref(),
            &graphic.presents[..]
        ),
        ("Datum A", "datum", Some(r"DTM\wA"), &[][..])
    );
    assert_eq!(graphic.views, vec!["Bore".to_owned()]);
    assert!(!read.pmi.views[0].items.contains(&"Datum A".to_owned()));

    // A callout whose size is not read stays graphics that names the size.
    let (_, record) = association("Bore diameter");
    let size = references(&record)[0].clone();
    let size_record = &all[&size];
    let broken = mutate(
        &file,
        &format!("{size}={size_record}"),
        &format!("{size}={}", size_record.replace("'diameter'", "'girth'")),
    );
    let read = read_step_pmi(broken.as_bytes()).unwrap();
    assert!(read.pmi.items.iter().all(|item| item.id != "Bore diameter"));
    let [graphic] = read.graphics.as_slice() else {
        panic!("{:?}", read.graphics);
    };
    let bore = part.round(2, "4", ["30", "25"]);
    let expected: Vec<_> = written_faces(&file, &report)
        .into_iter()
        .filter_map(|(id, face)| bore.contains(&face).then_some(id))
        .collect();
    assert_eq!(graphic.faces.as_ref(), Some(&expected));
    let mixed = mutate(
        &shared,
        &format!("{size}={size_record}"),
        &format!("{size}={}", size_record.replace("'diameter'", "'girth'")),
    );
    let mixed = read_step_pmi(mixed.as_bytes()).unwrap();
    assert!(mixed.pmi.views[0].items.contains(&"Datum A".to_owned()));
    let retained = mixed
        .graphics
        .iter()
        .find(|graphic| graphic.name == "Bore diameter")
        .unwrap();
    assert_eq!(retained.presents.as_slice(), std::slice::from_ref(&size));
    assert_eq!(retained.faces.as_ref(), Some(&expected));
    assert_eq!(retained.views, ["Bore"]);
    assert_eq!(graphic.presents, [size]);
    assert_eq!(graphic.text.as_deref(), Some(r"DIM\w⌀8 H7"));
    assert_eq!(graphic.kind, "diameter dimension");
}

#[test]
fn graphical_only_fixture_preserves_strokes_plane_face_and_view() {
    let source = include_str!("fixtures/graphical-only-pmi.step");
    let read = read_step_pmi(source.as_bytes()).unwrap();
    assert!(read.skipped.is_empty(), "{:?}", read.skipped);
    assert!(read.pmi.items.is_empty());
    assert_eq!(read.pmi.views[0].name, "Top");
    let [callout] = read.graphics.as_slice() else {
        panic!("expected one callout, got {:?}", read.graphics);
    };
    assert_eq!(callout.name, "Plate mark");
    assert_eq!(callout.kind, "surface roughness");
    assert_eq!(callout.text, None);
    assert_eq!(callout.faces, Some(vec!["#19".to_owned()]));
    assert_eq!(callout.views, ["Top"]);
    assert_eq!(
        callout.polylines_mm,
        [vec![v("2", "2", "5"), v("4", "4", "5"), v("8", "4", "5")]]
    );
    assert_eq!(
        callout.planes,
        [catallax_import::GraphicPlane {
            origin_mm: v("0", "0", "5"),
            normal: v("0", "0", "1"),
            reference: v("1", "0", "0"),
        }]
    );

    for (from, to) in [
        ("(#30,#31,#32)", "(#30,$,#32)"),
        ("(4.,4.,5.)", "(4.,4.,6.)"),
        ("'surface roughness',(#33)", "'surface roughness',(#33,$)"),
        ("'Plate mark',(#35)", "'Plate mark',(#35,$)"),
        ("#16,(#36)", "#16,(#36,$)"),
    ] {
        let broken = mutate(source, from, to);
        let broken = read_step_pmi(broken.as_bytes()).unwrap();
        assert!(broken.graphics.is_empty(), "{from}");
        assert_eq!(broken.pmi.views, read.pmi.views);
        assert!(
            broken.skipped.iter().any(|skip| skip.entity == "#36"),
            "{from}"
        );
    }
    let repeated = mutate(source, "(#30,#31,#32)", "(#30,#31,#31,#32)");
    let repeated = read_step_pmi(repeated.as_bytes()).unwrap();
    assert!(repeated.skipped.is_empty(), "{:?}", repeated.skipped);
    assert_eq!(
        repeated.graphics[0].polylines_mm[0],
        [
            v("2", "2", "5"),
            v("4", "4", "5"),
            v("4", "4", "5"),
            v("8", "4", "5")
        ]
    );

    let unsupported_callout = mutate(source, "DRAUGHTING_CALLOUT(", "UNSUPPORTED_CALLOUT(");
    let unsupported_callout = read_step_pmi(unsupported_callout.as_bytes()).unwrap();
    assert!(unsupported_callout.graphics.is_empty());
    assert_eq!(unsupported_callout.pmi.views[0].name, "Top");

    let subtype = mutate(source, "DRAUGHTING_CALLOUT(", "SURFACE_CONDITION_CALLOUT(");
    let subtype = read_step_pmi(subtype.as_bytes()).unwrap();
    assert!(subtype.graphics.is_empty());
    assert_eq!(subtype.pmi.views, read.pmi.views);
    assert!(
        subtype
            .skipped
            .iter()
            .any(|item| item.kind == "SURFACE_CONDITION_CALLOUT")
    );

    let complex = mutate(
        source,
        "#50=DRAUGHTING_MODEL('Top',(#38,#49),#9);",
        "#50=(CHARACTERIZED_OBJECT(*,*)CHARACTERIZED_REPRESENTATION()DRAUGHTING_MODEL()REPRESENTATION('Top',(#38,#49),#9));",
    );
    let complex = mutate(
        &complex,
        "#39=DRAUGHTING_MODEL('',(#38),#9);",
        "#39=(CHARACTERIZED_OBJECT(*,*)CHARACTERIZED_REPRESENTATION()DRAUGHTING_MODEL()REPRESENTATION('',(#38),#9));",
    );
    let complex = read_step_pmi(complex.as_bytes()).unwrap();
    assert!(complex.skipped.is_empty(), "{:?}", complex.skipped);
    assert_eq!(complex.pmi.views, read.pmi.views);
    assert_eq!(complex.graphics, read.graphics);

    let two_units = mutate(
        source,
        "ENDSEC;\nEND-ISO",
        "#90=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT($,.METRE.));\nENDSEC;\nEND-ISO",
    );
    for (units, usable) in [("(#8,#90)", false), ("(#90,#8)", false), ("(#90)", true)] {
        let changed = mutate(
            &two_units,
            "GLOBAL_UNIT_ASSIGNED_CONTEXT((#8))",
            &format!("GLOBAL_UNIT_ASSIGNED_CONTEXT({units})"),
        );
        let read = read_step_pmi(changed.as_bytes()).unwrap();
        if usable {
            assert!(read.skipped.is_empty(), "{:?}", read.skipped);
            assert_eq!(read.pmi.views[0].camera.distance, q("15000"));
            assert_eq!(
                read.graphics[0].polylines_mm[0][0],
                v("2000", "2000", "5000")
            );
        } else {
            assert!(read.graphics.is_empty());
            assert!(read.pmi.views.is_empty());
            assert_eq!(
                read.skipped
                    .iter()
                    .filter(|item| item.reason.contains("conflicting length units"))
                    .count(),
                2
            );
        }
    }

    let boxed = mutate(
        source,
        "#38=ANNOTATION_PLANE('',(#37),#16,(#36));",
        "#38=ANNOTATION_PLANE('',(#54),#55,(#36));",
    );
    let boxed = mutate(
        &boxed,
        "ENDSEC;\nEND-ISO",
        "#54=PRESENTATION_STYLE_ASSIGNMENT((#28));\n#55=PLANAR_BOX('',10.,10.,#15);\nENDSEC;\nEND-ISO",
    );
    let boxed = read_step_pmi(boxed.as_bytes()).unwrap();
    assert!(boxed.skipped.is_empty());
    assert_eq!(boxed.graphics, read.graphics);

    let repeated = mutate(
        source,
        "#50=DRAUGHTING_MODEL('Top',(#38,#49),#9);",
        "#50=DRAUGHTING_MODEL('',(#38,#49,#54,#55,#56),#9);",
    );
    let repeated = mutate(
        &repeated,
        "ENDSEC;\nEND-ISO",
        "#54=CAMERA_MODEL_D3('Top',#43,#48);\n#55=CAMERA_MODEL_D3('Top',#43,#48);\n#56=CAMERA_MODEL_D3(' ',#43,#48);\nENDSEC;\nEND-ISO",
    );
    let repeated = read_step_pmi(repeated.as_bytes()).unwrap();
    assert!(repeated.skipped.is_empty());
    assert_eq!(
        repeated.graphics[0].views,
        ["Top", "Top (2)", "Top (3)", "View #56"]
    );

    let unsupported = mutate(
        source,
        "ENDSEC;\nEND-ISO",
        r"#900000=CHARACTERIZED_ITEM_WITHIN_REPRESENTATION('','',#36,#39);
#900001=PROPERTY_DEFINITION('pmi validation property','',#900000);
#900002=DESCRIPTIVE_REPRESENTATION_ITEM('equivalent unicode string','\PE\\S\*\S\U\S\b');
#900003=REPRESENTATION('',(#900002),#9);
#900004=PROPERTY_DEFINITION_REPRESENTATION(#900001,#900003);
ENDSEC;
END-ISO",
    );
    let unsupported = read_step_pmi(unsupported.as_bytes()).unwrap();
    assert_eq!(unsupported.graphics, read.graphics);
    assert!(
        unsupported
            .skipped
            .iter()
            .any(|item| item.reason.contains("unsupported STEP encoding"))
    );

    let scoped = mutate(
        source,
        "ENDSEC;\nEND-ISO",
        "#54=CAMERA_MODEL_D3('Side',#43,#48);\n\
         #55=DRAUGHTING_MODEL('Side',(#38,#54),#9);\n\
         #56=CHARACTERIZED_ITEM_WITHIN_REPRESENTATION('','',#36,#55);\n\
         #57=PROPERTY_DEFINITION('pmi validation property','',#56);\n\
         #58=DESCRIPTIVE_REPRESENTATION_ITEM('equivalent unicode string','Side mark');\n\
         #59=REPRESENTATION('',(#58),#9);\n\
         #60=PROPERTY_DEFINITION_REPRESENTATION(#57,#59);\n\
         ENDSEC;\nEND-ISO",
    );
    let scoped = read_step_pmi(scoped.as_bytes()).unwrap();
    assert!(scoped.skipped.is_empty(), "{:?}", scoped.skipped);
    let [top, side] = scoped.graphics.as_slice() else {
        panic!("{:?}", scoped.graphics);
    };
    assert_eq!(top, callout);
    assert_eq!(side.views, ["Side"]);
    assert_eq!(side.text.as_deref(), Some("Side mark"));
    assert_eq!(side.faces, None);

    let defaults = mutate(source, "#41,#42,#14", "#41,#42,$");
    let defaults = read_step_pmi(defaults.as_bytes()).unwrap();
    assert!(defaults.skipped.is_empty(), "{:?}", defaults.skipped);
    assert_eq!(defaults.pmi.views[0].camera.up, v("0", "-1", "0"));

    for (from, to, reason) in [
        ("(-10.,-8.,15.)", "(-9.,-8.,15.)", "off centre"),
        ("(-10.,-8.,15.)", "(-10.,-8.,14.)", "window depth"),
        (
            "AXIS2_PLACEMENT_3D('',#45,$,$)",
            "AXIS2_PLACEMENT_2D('',#45,$)",
            "3D placement",
        ),
        (
            "AXIS2_PLACEMENT_3D('',#45,$,$)",
            "AXIS2_PLACEMENT_3D('',#45,#42,#14)",
            "rotated",
        ),
        ("15.,.F.,.F.,#47", "15.,.F.,.T.,#47", "clipping"),
        ("15.,15.,.F.", "0.,15.,.F.", "positive view distance"),
        ("#41,#42,#14", "#41,#42,#13", "placement axes are parallel"),
        ("(0.,0.,0.)", "(1.,0.,0.)", "projection point"),
    ] {
        let broken = mutate(source, from, to);
        let read = read_step_pmi(broken.as_bytes()).unwrap();
        assert!(read.pmi.views.is_empty());
        assert!(
            read.skipped.iter().any(|item| item.reason.contains(reason)),
            "{:?}",
            read.skipped
        );
    }
}

#[test]
fn datum_target_survives_oblique_straight_trim() {
    use catallax_shapes::{RPoint3, RVector3};
    use catallax_topology::block_in_frame;

    let mut model = Model::new();
    let vector = |x, y, z| RVector3::new(q(x), q(y), q(z));
    let body = block_in_frame(
        &mut model,
        &RPoint3::origin(),
        &[
            vector("1", "3", "0"),
            vector("-3", "1", "0"),
            vector("0", "0", "5"),
        ],
    )
    .unwrap();
    let part = Bracket { model, body };
    let [face] = part.planes(2, "5").try_into().unwrap();
    let normal_target = target(1, face, v("-1", "2", "5"), TargetKind::Point);
    let pmi = Pmi {
        items: vec![
            item(
                "Datum A",
                ItemKind::Datum(Datum {
                    letter: "A".into(),
                    feature: vec![face],
                }),
            ),
            item("A1", normal_target.clone()),
        ],
        ..Pmi::default()
    };
    let (file, report) = write(&part, &pmi);
    let (back, skipped) = read_back(&file, &report);
    assert!(skipped.is_empty(), "{skipped:?}");
    assert_eq!(kind_of(&back, "A1"), &normal_target);

    let all = records(&file);
    let top = written_faces(&file, &report)
        .into_iter()
        .find_map(|(id, native)| (native == face).then_some(id))
        .unwrap();
    let bound = &references(&all[&top])[0];
    let rim = &references(&all[bound])[0];
    let oriented_edge = &references(&all[rim])[0];
    let edge = &references(&all[oriented_edge])[0];
    let line = &references(&all[edge])[2];
    let vector = &references(&all[line])[1];
    let direction = &references(&all[vector])[0];
    let broken = mutate(
        &file,
        &format!("{direction}={}", all[direction]),
        &format!("{direction}=DIRECTION('',(1.,0.,0.));"),
    );
    let read = read_step_pmi(broken.as_bytes()).unwrap();
    assert!(!read.sources.contains_key("A1"));
    assert!(!read.sources.contains_key("Datum A"));
}

#[test]
fn target_refuses_a_face_moved_by_step_rounding() {
    use catallax_shapes::{RVector3, Rotation};
    use catallax_topology::block;

    let z = "0.1234567890123456789";
    let pose = Pose::from_rational(Rotation::identity(), RVector3::new(q("0"), q("0"), q(z)));
    let mut model = Model::new();
    let body = block(&mut model, [q("10"), q("10"), q("10")], &pose).unwrap();
    let part = Bracket { model, body };
    let face = part.planes(2, z)[0];
    let pmi = Pmi {
        items: vec![
            item(
                "Datum A",
                ItemKind::Datum(Datum {
                    letter: "A".into(),
                    feature: vec![face],
                }),
            ),
            item("A1", target(1, face, v("5", "5", z), TargetKind::Point)),
        ],
        part: PartNotes::default(),
        views: Vec::new(),
    };
    catallax_step::pmi::check_on_model(&part.model, part.body, &pmi, &mut Decider::default())
        .unwrap();
    let mut out = Vec::new();
    let refusal = write_step_with_pmi(
        &part.model,
        part.body,
        &StepOptions::default(),
        &pmi,
        &mut Decider::default(),
        &mut out,
    )
    .unwrap_err();
    assert!(
        refusal.detail.contains("A1") && refusal.detail.contains("written"),
        "{refusal:?}"
    );
    assert!(out.is_empty());
}

#[test]
fn written_position_needs_coaxial_trimmed_cylinders() {
    use catallax_numbers::Real;
    use catallax_shapes::{RVector3, Rotation, Vector3};

    let mut model = Model::new();
    let mut dec = Decider::default();
    let vector = |x, y, z| Vector3::new(Real::from(q(x)), Real::from(q(y)), Real::from(q(z)));
    let rotation =
        Rotation::from_frame(&vector("3", "4", "0"), &vector("0", "0", "1"), &mut dec).unwrap();
    let a = cylinder(
        &mut model,
        q("2"),
        q("0.5"),
        &Pose::from_rational(rotation.clone(), RVector3::zero()),
    )
    .unwrap();
    let b = cylinder(
        &mut model,
        q("2"),
        q("0.5"),
        &Pose::from_rational(rotation, RVector3::new(q("1"), q("4/3"), q("0"))),
    )
    .unwrap();
    let axes = [a, b].map(|body| {
        model
            .faces_of_body(body)
            .into_iter()
            .find(|face| {
                matches!(
                    model.geom().surface(model.face(*face).surface),
                    Surface::Cylinder(_)
                )
            })
            .unwrap()
    });
    let mut body = model.body(a).clone();
    body.lumps.extend(&model.body(b).lumps);
    let body = model.add_body(body);
    let part = Bracket { model, body };
    let support = part
        .faces()
        .into_iter()
        .find(|face| !axes.contains(face))
        .unwrap();
    for (feature, datum, frame) in [
        (
            axes.to_vec(),
            vec![axes[0]],
            PositionFrame::Single {
                segment: segment("0.05", None, vec![bounded("A", MaterialCondition::Maximum)]),
                projected_zone: None,
            },
        ),
        (
            axes.to_vec(),
            vec![support],
            PositionFrame::Composite {
                upper: segment("0.05", None, vec![DatumRef::new("A")]),
                lower: segment("0.025", None, Vec::new()),
            },
        ),
        (
            vec![axes[0]],
            axes.to_vec(),
            PositionFrame::Single {
                segment: segment("0.05", None, vec![bounded("A", MaterialCondition::Maximum)]),
                projected_zone: None,
            },
        ),
    ] {
        let pmi = Pmi {
            items: vec![
                item(
                    "Datum A",
                    ItemKind::Datum(Datum {
                        letter: "A".into(),
                        feature: datum,
                    }),
                ),
                item(
                    "Bore position",
                    ItemKind::Position(PositionControl {
                        features: vec![feature],
                        frame,
                    }),
                ),
            ],
            part: PartNotes {
                standard: Some(catallax_step::pmi::DimensioningStandard::AsmeY14_5_2018),
                ..PartNotes::default()
            },
            ..Pmi::default()
        };
        catallax_step::pmi::check_on_model(&part.model, part.body, &pmi, &mut dec).unwrap();
        let error = try_write(&part, &pmi).unwrap_err();
        assert!(
            error.contains("written faces do not prove one cylinder"),
            "{error}"
        );
    }
}

#[test]
fn rounded_step_geometry_cannot_reverse_a_projection() {
    use catallax_shapes::{Isometry, RVector3, Rotation};
    let (mut part, _) = hole_fixture();
    let rotation = Rotation::from_quaternion(q("3"), q("-1"), q("0"), q("0")).unwrap();
    let pose = Pose::from_rational(
        rotation,
        RVector3::new(q("0"), Rational::from_ratio(-2, 3).unwrap(), q("1")),
    );
    part.body = catallax_topology::transformed_body(
        &mut part.model,
        part.body,
        &Isometry::from(&pose),
        Clone::clone,
    );
    let head: Vec<_> = part
        .faces()
        .into_iter()
        .filter(|&face| {
            part.model.face_origins(face).iter().any(|origin| {
                origin.step.as_str() == "counterbore"
                    && matches!(
                        origin.role,
                        Role::Hole {
                            part: HolePart::CounterboreWall,
                            ..
                        }
                    )
            })
        })
        .collect();
    let entry = part
        .faces()
        .into_iter()
        .find(|&face| {
            let Surface::Plane(plane) = part.model.geom().surface(part.model.face(face).surface)
            else {
                return false;
            };
            let [_, y, z] = plane.origin.rational_components().unwrap();
            let [nx, ny, nz] = plane.normal.rational_components().unwrap();
            nx.is_zero()
                && ((&q("0.8") * ny) - (&q("0.6") * nz)).is_zero()
                && &(&q("0.6") * y) + &(&q("0.8") * z) == q("0.4")
        })
        .unwrap();
    let pmi = Pmi {
        items: vec![item(
            "Projection",
            ItemKind::Position(PositionControl {
                features: vec![head],
                frame: PositionFrame::Single {
                    segment: segment("0.05", None, Vec::new()),
                    projected_zone: Some(ProjectedZone {
                        height_mm: q("8"),
                        direction: v("0", "-0.6", "-0.8"),
                        entry_faces: vec![entry],
                    }),
                },
            }),
        )],
        ..Pmi::default()
    };
    catallax_step::pmi::check_on_model(&part.model, part.body, &pmi, &mut Decider::default())
        .unwrap();
    let error = try_write(&part, &pmi).unwrap_err();
    assert!(
        error.contains("projected zone cannot be proved on the written geometry"),
        "{error}"
    );
}

#[test]
fn rectangle_checks_the_rounded_corner_arc_instead_of_its_circle() {
    use catallax_numbers::Real;
    use catallax_shapes::{
        Circle, Plane, Point3, Vector3,
        distance::{PlanarBoundary, PlanarEdge},
    };
    let point = |x, y| Point3::new(Real::from(q(x)), Real::from(q(y)), Real::zero());
    let normal = Vector3::new(Real::zero(), Real::zero(), Real::one());
    let along = Vector3::new(Real::one(), Real::zero(), Real::zero());
    let boundary = PlanarBoundary {
        plane: Plane {
            origin: point("0", "0"),
            normal: normal.clone(),
            ref_dir: along.clone(),
        },
        vertices: vec![
            point("1", "0"),
            point("4", "0"),
            point("4", "4"),
            point("0", "4"),
            point("0", "1"),
        ],
        edges: vec![
            PlanarEdge::Segment([0, 1]),
            PlanarEdge::Segment([1, 2]),
            PlanarEdge::Segment([2, 3]),
            PlanarEdge::Segment([3, 4]),
            PlanarEdge::Arc {
                circle: Circle {
                    center: point("1", "1"),
                    normal,
                    ref_dir: along,
                    radius_sq: Real::one(),
                },
                ends: Some([4, 0]),
            },
        ],
    };
    let rectangle = |x, y, width| DatumTarget {
        datum: "A".into(),
        number: 1,
        face: (),
        point_mm: v(x, y, "0"),
        kind: TargetKind::Area(TargetArea::Rectangular {
            length_mm: q("0.6"),
            width_mm: q(width),
            direction: v("1", "0", "0"),
        }),
    };
    for (target, accepted) in [
        (rectangle("1.9", "1", "0.2"), true),
        (rectangle("0.5", "0.4", "0.3"), false),
    ] {
        let reason =
            catallax_step::pmi::off_face(&boundary, &target, &mut Decider::default()).unwrap();
        assert_eq!(reason.is_none(), accepted, "{reason:?}");
    }
}

#[test]
fn a_counterbore_shoulder_is_not_an_outward_entry() {
    let (part, pmi) = hole_fixture();
    let ItemKind::Size(size) = kind_of(&pmi, "counterbore") else {
        panic!("expected counterbore");
    };
    let bottom = part.planes(2, "0")[0];
    let shoulder = part.planes(2, "2.25")[0];
    let projection = |entry| {
        ItemKind::Position(PositionControl {
            features: vec![size.features[1].clone()],
            frame: PositionFrame::Single {
                segment: segment("0.05", None, Vec::new()),
                projected_zone: Some(ProjectedZone {
                    height_mm: q("8"),
                    direction: v("0", "0", "-1"),
                    entry_faces: vec![entry],
                }),
            },
        })
    };
    let mut pmi = Pmi {
        items: vec![
            item(
                "Datum A",
                ItemKind::Datum(Datum {
                    letter: "A".into(),
                    feature: vec![shoulder],
                }),
            ),
            item("Projection", projection(bottom)),
        ],
        part: PartNotes::default(),
        views: Vec::new(),
    };
    let (file, report) = write(&part, &pmi);
    let (back, skipped) = read_back(&file, &report);
    assert!(skipped.is_empty(), "{skipped:?}");
    assert_eq!(back, pmi);

    pmi.items[1].kind = projection(shoulder);
    let err = try_write(&part, &pmi).unwrap_err();
    assert!(
        err.contains("points back into its controlled cylinder"),
        "{err}"
    );

    let all = records(&file);
    let (datum, _) = all
        .iter()
        .find(|(_, body)| body.starts_with("DATUM_FEATURE("))
        .unwrap();
    let projected = all
        .values()
        .find(|body| body.starts_with("PROJECTED_ZONE_DEFINITION("))
        .unwrap();
    let start = &references(projected)[1];
    let modified = mutate(&file, projected, &projected.replace(start, datum));
    let (back, skipped) = read_back(&modified, &report);
    assert_eq!(back.items, vec![pmi.items[0].clone()]);
    assert!(
        skipped.iter().any(|item| item
            .reason
            .contains("points back into its controlled cylinder")),
        "{skipped:?}"
    );
}

#[test]
fn a_coplanar_disk_closes_a_projected_entry() {
    use catallax_numbers::Real;
    use catallax_shapes::{Circle, Curve};
    use catallax_topology::{Coedge, Edge, Lump, Sense, Shell};

    let mut model = Model::new();
    let body = cylinder(&mut model, q("4"), q("5"), &Pose::identity()).unwrap();
    let mut part = Bracket { model, body };
    let barrel = part.round(2, "4", ["0", "0"]);
    let top = part.planes(2, "5")[0];
    let mut annulus = part.model.face(top).clone();
    let model = &mut part.model;
    let Surface::Plane(plane) = model.geom().surface(annulus.surface).clone() else {
        panic!("expected top plane")
    };
    let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
        center: plane.origin,
        normal: plane.normal,
        ref_dir: plane.ref_dir,
        radius_sq: Real::one(),
    }));
    let edge = model.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: None,
        end: None,
        tolerance: None,
    });
    let inner = model.add_coedge(Coedge {
        edge,
        sense: Sense::Reversed,
    });
    let outer = model.add_coedge(Coedge {
        edge,
        sense: Sense::Same,
    });
    let inner = model.add_loop(Loop::Coedges(vec![inner]));
    let outer = model.add_loop(Loop::Coedges(vec![outer]));
    annulus.loops = annulus
        .loops
        .iter()
        .map(|id| {
            let Loop::Coedges(coedges) = model.loop_(*id).clone() else {
                panic!("expected top ring")
            };
            let coedges = coedges
                .into_iter()
                .map(|coedge| model.add_coedge(model.coedge(coedge).clone()))
                .collect();
            model.add_loop(Loop::Coedges(coedges))
        })
        .collect();
    annulus.loops.push(inner);
    let mut disk = annulus.clone();
    disk.loops = vec![outer];
    let annulus = model.add_face(annulus);
    let disk = model.add_face(disk);
    let mut faces = model.faces_of_body(body);
    faces.retain(|face| *face != top);
    faces.extend([annulus, disk]);
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let mut split = model.body(body).clone();
    split.lumps = vec![lump];
    part.body = model.add_body(split);

    let mut pmi = Pmi {
        items: vec![
            item(
                "Datum A",
                ItemKind::Datum(Datum {
                    letter: "A".into(),
                    feature: vec![annulus],
                }),
            ),
            item(
                "Projection",
                ItemKind::Position(PositionControl {
                    features: vec![barrel.clone()],
                    frame: PositionFrame::Single {
                        segment: segment("0.05", None, Vec::new()),
                        projected_zone: Some(ProjectedZone {
                            height_mm: q("8"),
                            direction: v("0", "0", "1"),
                            entry_faces: vec![annulus],
                        }),
                    },
                }),
            ),
        ],
        part: PartNotes::default(),
        views: Vec::new(),
    };
    let err = try_write(&part, &pmi).unwrap_err();
    assert!(err.contains("axis meets a closed end"), "{err}");

    pmi.items[1].kind = ItemKind::Tolerance(GeometricTolerance {
        characteristic: Characteristic::Position,
        value: q("0.05"),
        diameter_zone: true,
        material: None,
        datums: Vec::new(),
        features: vec![barrel],
    });
    let (file, _) = write(&part, &pmi);
    let all = records(&file);
    let (zone, _) = all
        .iter()
        .find(|(_, body)| body.starts_with("TOLERANCE_ZONE("))
        .unwrap();
    let (datum, _) = all
        .iter()
        .find(|(_, body)| body.starts_with("DATUM_FEATURE("))
        .unwrap();
    let (_, tol) = all
        .iter()
        .find(|(_, body)| body.starts_with("POSITION_TOLERANCE('Projection'"))
        .unwrap();
    let measure = references(tol)[0].clone();
    let addition = format!("#900000=PROJECTED_ZONE_DEFINITION({zone},(),{datum},{measure});\n");
    let source = mutate(
        &file,
        "ENDSEC;\nEND-ISO",
        &format!("{addition}ENDSEC;\nEND-ISO"),
    );
    let read = read_step_pmi(source.as_bytes()).unwrap();
    assert!(read.pmi.items.iter().all(|item| item.id != "Projection"));
    assert!(
        read.skipped
            .iter()
            .any(|item| item.reason.contains("axis meets a closed end")),
        "{:?}",
        read.skipped
    );
}

#[test]
fn a_boss_blocks_projection_across_a_recess() {
    for raised in [false, true] {
        use catallax_boolean::{BoolOp, boolean};
        use catallax_shapes::{RVector3, Rotation};
        let mut model = Model::new();
        let mut dec = Decider::default();
        let at = |z| Pose::from_rational(Rotation::identity(), RVector3::new(q("0"), q("0"), q(z)));
        let stock = cylinder(&mut model, q("5"), q("6"), &Pose::identity()).unwrap();
        let outer_cut = cylinder(&mut model, q("4"), q("5"), &at("2")).unwrap();
        let boss = cylinder(&mut model, q("1"), q("7"), &at("1")).unwrap();
        let stock = if raised {
            boolean(&mut model, stock, boss, BoolOp::Unite, &mut dec)
                .unwrap()
                .body
        } else {
            stock
        };
        let annular_cut = boolean(&mut model, outer_cut, boss, BoolOp::Subtract, &mut dec)
            .unwrap()
            .body;
        let body = boolean(&mut model, stock, annular_cut, BoolOp::Subtract, &mut dec)
            .unwrap()
            .body;
        let body = if raised {
            let other = cylinder(&mut model, q("0.5"), q("1"), &at("7.5")).unwrap();
            let mut mixed = model.body(body).clone();
            mixed.lumps.extend(&model.body(other).lumps);
            model.add_body(mixed)
        } else {
            body
        };
        let part = Bracket { model, body };
        let rim = part
            .planes(2, "6")
            .into_iter()
            .find(|face| part.model.face(*face).loops.len() == 2)
            .unwrap();
        let wall = part.round(2, "4", ["0", "0"]);
        let make_pmi = |kind| Pmi {
            items: vec![
                item(
                    "Rim",
                    ItemKind::Datum(Datum {
                        letter: "R".into(),
                        feature: vec![rim],
                    }),
                ),
                item("Recess wall", kind),
            ],
            ..Pmi::default()
        };
        let pmi = make_pmi(ItemKind::Position(PositionControl {
            features: vec![wall.clone()],
            frame: PositionFrame::Single {
                segment: segment("0.05", None, Vec::new()),
                projected_zone: Some(ProjectedZone {
                    height_mm: q("8"),
                    direction: v("0", "0", "1"),
                    entry_faces: vec![rim],
                }),
            },
        }));
        let error =
            catallax_step::pmi::check_on_model(&part.model, part.body, &pmi, &mut dec).unwrap_err();
        assert!(
            error.detail.contains("axis meets a closed end"),
            "{error:?}"
        );
        assert!(
            try_write(&part, &pmi)
                .unwrap_err()
                .contains("axis meets a closed end")
        );
        if raised {
            let mut short = pmi.clone();
            let ItemKind::Position(control) = &mut short.items[1].kind else {
                unreachable!()
            };
            let PositionFrame::Single {
                projected_zone: Some(zone),
                ..
            } = &mut control.frame
            else {
                unreachable!()
            };
            zone.height_mm = q("1");
            let error =
                catallax_step::pmi::check_on_model(&part.model, part.body, &short, &mut dec)
                    .unwrap_err();
            assert!(
                error.detail.contains("axis starts inside the part"),
                "{error:?}"
            );
        }
        let (file, _) = write(
            &part,
            &make_pmi(ItemKind::Tolerance(GeometricTolerance {
                characteristic: Characteristic::Position,
                value: q("0.05"),
                diameter_zone: true,
                material: None,
                datums: Vec::new(),
                features: vec![wall.clone()],
            })),
        );
        let all = records(&file);
        let start = line_with(&file, &["=DATUM_FEATURE("]);
        let zone = line_with(&file, &["=TOLERANCE_ZONE("]);
        let control = line_with(&file, &["=POSITION_TOLERANCE("]);
        let height = &references(&all[&control])[0];
        let file = mutate(
            &file,
            "ENDSEC;\nEND-ISO",
            &format!(
                "#900000=PROJECTED_ZONE_DEFINITION({zone},(),{start},{height});\nENDSEC;\nEND-ISO"
            ),
        );
        let read = read_step_pmi(file.as_bytes()).unwrap();
        assert_eq!(read.sources.keys().collect::<Vec<_>>(), [&"Rim".to_owned()]);
        // The written height is the tolerance's own 0.05 mm: too short to
        // reach the raised boss's cap, so that zone is refused for starting in
        // the boss, and the flush one for the cap at its start.
        let reason = if raised {
            "axis starts inside the part"
        } else {
            "axis meets a closed end"
        };
        assert!(
            read.skipped.iter().any(|item| item.reason.contains(reason)),
            "{:?}",
            read.skipped
        );
    }
}

#[test]
fn a_void_target_round_trips_on_the_inner_wall() {
    use catallax_boolean::{BoolOp, boolean};
    use catallax_shapes::{RVector3, Rotation};
    let mut model = Model::new();
    let outer =
        catallax_topology::block(&mut model, [q("10"), q("10"), q("10")], &Pose::identity())
            .unwrap();
    let inner = catallax_topology::block(
        &mut model,
        [q("2"), q("2"), q("2")],
        &Pose::from_rational(Rotation::identity(), RVector3::new(q("4"), q("4"), q("4"))),
    )
    .unwrap();
    let body = boolean(
        &mut model,
        outer,
        inner,
        BoolOp::Subtract,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    let part = Bracket { model, body };
    let face = part.planes(0, "4")[0];
    let pmi = Pmi {
        items: vec![
            item(
                "Void datum",
                ItemKind::Datum(Datum {
                    letter: "A".into(),
                    feature: vec![face],
                }),
            ),
            item(
                "A1",
                ItemKind::DatumTarget(DatumTarget {
                    datum: "A".into(),
                    number: 1,
                    face,
                    point_mm: v("4", "5", "5"),
                    kind: TargetKind::Point,
                }),
            ),
        ],
        ..Pmi::default()
    };
    let (file, _) = write(&part, &pmi);
    let all = entities(&file);
    let (_, oriented) = all
        .values()
        .find(|(kind, _)| kind == "ORIENTED_CLOSED_SHELL")
        .unwrap();
    let shell = &references(oriented)[0];
    let wall = references(&all[shell].1)
        .into_iter()
        .find(|face| {
            let plane = &references(&all[face].1).pop().unwrap();
            let placement = &references(&all[plane].1)[0];
            let frame = references(&all[placement].1);
            let origin = numbers(&all[&frame[0]].1);
            let axis = numbers(&all[&frame[1]].1);
            (origin[0] - 4.0).abs() < 1e-9 && (axis[0].abs() - 1.0).abs() < 1e-9
        })
        .unwrap();
    let read = read_step_pmi(file.as_bytes()).unwrap();
    assert!(read.skipped.is_empty(), "{:?}", read.skipped);
    assert_eq!(kind_of(&read.pmi, "Void datum").faces(), [&wall]);
    let ItemKind::DatumTarget(target) = kind_of(&read.pmi, "A1") else {
        panic!("expected a point target");
    };
    assert_eq!(target.point_mm, v("4", "5", "5"));
    assert_eq!(target.face, wall);
}

#[test]
fn an_oblique_camera_and_target_write_exact_direction_ratios() {
    let source = include_str!("fixtures/graphical-only-pmi.step");
    let source = mutate(source, "#41,#42,#14", "#41,#42,$");
    let source = mutate(
        &source,
        "#42=DIRECTION('',(0.,0.,-1.));",
        "#42=DIRECTION('',(1.,1.,-1.));",
    );
    let imported = read_step_pmi(source.as_bytes())
        .unwrap()
        .pmi
        .views
        .pop()
        .unwrap();
    let third = &q("1") / &q("3");
    let expected_up = [q("0"), -third.clone(), -third.clone()];
    assert_eq!(imported.camera.up, expected_up);
    let mut model = Model::new();
    let body = catallax_topology::block(&mut model, [q("10"), q("10"), q("5")], &Pose::identity())
        .unwrap();
    let part = Bracket { model, body };
    let top = part.planes(2, "5")[0];
    let pmi = Pmi {
        views: vec![SavedView {
            items: vec!["Datum A".into()],
            ..imported
        }],
        items: vec![
            item(
                "Datum A",
                ItemKind::Datum(Datum {
                    letter: "A".into(),
                    feature: vec![top],
                }),
            ),
            item(
                "A1",
                ItemKind::DatumTarget(DatumTarget {
                    datum: "A".into(),
                    number: 1,
                    face: top,
                    point_mm: v("5", "5", "5"),
                    kind: TargetKind::Line {
                        length_mm: q("2"),
                        direction: [third.clone(), &third * &q("2"), q("0")],
                    },
                }),
            ),
        ],
        ..Pmi::default()
    };
    let (file, report) = write(&part, &pmi);
    let (back, skipped) = read_back(&file, &report);
    assert!(skipped.is_empty(), "{skipped:?}");
    let result = &back.views[0];
    assert_eq!(result.items, ["Datum A"]);
    assert_eq!(result.camera.direction, v("1", "1", "-1"));
    assert!(result.camera.up[0].is_zero());
    assert!(result.camera.up[1] < q("0"));
    assert_eq!(result.camera.up[1], result.camera.up[2]);
    let ItemKind::DatumTarget(target) = kind_of(&back, "A1") else {
        panic!("expected a datum target");
    };
    let TargetKind::Line {
        direction,
        length_mm,
    } = &target.kind
    else {
        panic!("expected a line");
    };
    assert_eq!(target.face, top);
    assert_eq!(*length_mm, q("2"));
    assert!(direction[0] > q("0"));
    assert_eq!(direction[1], &direction[0] * &q("2"));
    assert!(direction[2].is_zero());
}

#[test]
fn authored_view_names_are_reserved_before_suffixes() {
    let source = include_str!("fixtures/graphical-only-pmi.step");
    // Two more views of the same plane, named Top and Top (2). The repeat
    // must not take Top (2), which a later camera authored.
    let more = "#60=CAMERA_MODEL_D3('Top',#43,#48);\n\
                #61=DRAUGHTING_MODEL('Top',(#38,#60),#9);\n\
                #62=MECHANICAL_DESIGN_AND_DRAUGHTING_RELATIONSHIP('','',#61,#39);\n\
                #63=CAMERA_MODEL_D3('Top (2)',#43,#48);\n\
                #64=DRAUGHTING_MODEL('Top (2)',(#38,#63),#9);\n\
                #65=MECHANICAL_DESIGN_AND_DRAUGHTING_RELATIONSHIP('','',#64,#39);\n";
    let named = mutate(
        source,
        "ENDSEC;\nEND-ISO",
        &format!("{more}ENDSEC;\nEND-ISO"),
    );
    let read = read_step_pmi(named.as_bytes()).unwrap();
    assert!(read.skipped.is_empty(), "{:?}", read.skipped);
    let names: Vec<&str> = read
        .pmi
        .views
        .iter()
        .map(|view| view.name.as_str())
        .collect();
    assert_eq!(names, ["Top", "Top (3)", "Top (2)"]);
    let [callout] = read.graphics.as_slice() else {
        panic!("expected one callout, got {:?}", read.graphics);
    };
    assert_eq!(callout.views, ["Top", "Top (3)", "Top (2)"]);
}

#[test]
fn a_refused_view_keeps_shared_callouts_in_other_views() {
    let source = include_str!("fixtures/graphical-only-pmi.step");
    // A second view of the same annotation plane maps the part geometry
    // through a 5 mm translation, which the reader refuses for that view.
    let side = "#70=CARTESIAN_POINT('',(0.,0.,0.));\n\
                #71=AXIS2_PLACEMENT_3D('',#70,#13,#14);\n\
                #72=CARTESIAN_POINT('',(5.,0.,0.));\n\
                #73=AXIS2_PLACEMENT_3D('',#72,#13,#14);\n\
                #74=REPRESENTATION_MAP(#71,#22);\n\
                #75=MAPPED_ITEM('',#74,#73);\n\
                #76=CAMERA_MODEL_D3('Side',#43,#48);\n\
                #77=DRAUGHTING_MODEL('Side',(#38,#75,#76),#9);\n\
                #78=MECHANICAL_DESIGN_AND_DRAUGHTING_RELATIONSHIP('','',#77,#39);\n";
    let changed = mutate(
        source,
        "ENDSEC;\nEND-ISO",
        &format!("{side}ENDSEC;\nEND-ISO"),
    );
    let read = read_step_pmi(changed.as_bytes()).unwrap();
    let names: Vec<&str> = read
        .pmi
        .views
        .iter()
        .map(|view| view.name.as_str())
        .collect();
    assert_eq!(names, ["Top"]);
    assert!(
        read.skipped.iter().any(
            |item| item.entity == "#77" && item.reason.contains("nonidentity geometry mapping")
        ),
        "{:?}",
        read.skipped
    );
    let [callout] = read.graphics.as_slice() else {
        panic!("the shared callout was lost: {:?}", read.skipped);
    };
    assert_eq!(callout.name, "Plate mark");
    assert_eq!(callout.views, ["Top"]);
    assert_eq!(callout.faces, Some(vec!["#19".to_owned()]));

    let unresolved = mutate(&changed, "(#38,#75,#76),#9", "(#38,#75,#76),#900001");
    let unresolved = read_step_pmi(unresolved.as_bytes()).unwrap();
    assert_eq!(unresolved.pmi.views, read.pmi.views);
    assert_eq!(unresolved.graphics, read.graphics);
    assert!(
        unresolved
            .skipped
            .iter()
            .any(|item| item.entity == "#77" && item.reason.contains("context")),
        "{:?}",
        unresolved.skipped
    );
    let extra = "#79=AXIS2_PLACEMENT_3D('',#70,#14,#13);\n\
                 #80=PLANE('',#79);\n\
                 #81=ANNOTATION_PLANE('',(#37),#80,(#36));\n\
                 #82=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT($,.METRE.));\n\
                 #83=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#82))REPRESENTATION_CONTEXT('',''));\n\
                 #84=DRAUGHTING_MODEL_ITEM_ASSOCIATION('Plate mark','',#24,#77,#36);\n";
    for (members, units) in [("(#81,#75,#76)", "#9"), ("(#38,#76)", "#83")] {
        let changed = mutate(&changed, "(#38,#75,#76),#9", &format!("{members},{units}"));
        let changed = mutate(
            &changed,
            "ENDSEC;\nEND-ISO",
            &format!("{extra}ENDSEC;\nEND-ISO"),
        );
        let retained = read_step_pmi(changed.as_bytes()).unwrap();
        assert_eq!(retained.pmi.views, read.pmi.views);
        assert_eq!(retained.graphics, read.graphics);
        assert!(retained.skipped.iter().any(|item| item.entity == "#77"));
    }
}

#[test]
fn only_fitted_diameters_get_saved_view_callouts() {
    let part = bracket();
    let (pmi, _) = viewed_pmi(&part);
    for limits in [
        None,
        Some(Limits::Basic),
        deviations("0.1", "-0.1"),
        Some(Limits::Fit {
            deviation: "g".into(),
            grade: "6".into(),
        }),
    ] {
        let mut pmi = pmi.clone();
        let bore = pmi
            .items
            .iter_mut()
            .find(|item| item.id == "Bore diameter")
            .unwrap();
        let ItemKind::Size(size) = &mut bore.kind else {
            panic!("expected a diameter");
        };
        if matches!(&limits, Some(Limits::Fit { .. })) {
            let (mut shaft, _) = bracket_pmi(&part);
            let ItemKind::Size(shaft) = shaft
                .items
                .remove(
                    shaft
                        .items
                        .iter()
                        .position(|item| item.id == "Boss diameter")
                        .unwrap(),
                )
                .kind
            else {
                panic!("expected boss diameter");
            };
            size.features = shaft.features;
            size.through = false;
        }
        size.limits = limits.clone();
        let Err(reason) = try_write(&part, &pmi) else {
            panic!("saved view accepted {limits:?}");
        };
        assert!(reason.contains("Bore diameter has no callout"), "{reason}");
        pmi.views[0].items.retain(|item| item != "Bore diameter");
        let (file, _) = write(&part, &pmi);
        let read = read_step_pmi(file.as_bytes()).unwrap();
        assert!(read.skipped.is_empty(), "{:?}", read.skipped);
        assert!(
            matches!(kind_of(&read.pmi, "Bore diameter"), ItemKind::Size(size) if size.limits == limits)
        );
    }
}

#[test]
fn datum_carriers_follow_the_ap242_rules() {
    let part = bracket();
    let f = gdt_faces(&part);
    let (file, report) = write(&part, &gdt_pmi(&f));
    let ids = |pmi: &Pmi<FaceId>| {
        let mut ids: Vec<String> = pmi.items.iter().map(|item| item.id.clone()).collect();
        ids.sort();
        ids
    };
    let reasons = |skipped: &[catallax_import::SkippedPmi]| {
        skipped
            .iter()
            .map(|skip| skip.reason.clone())
            .collect::<Vec<_>>()
            .join("; ")
    };
    let (all_items, skipped) = read_back(&file, &report);
    assert!(skipped.is_empty(), "{}", reasons(&skipped));
    let record = |patterns: &[&str]| {
        file.lines()
            .find(|line| patterns.iter().all(|pattern| line.contains(pattern)))
            .unwrap_or_else(|| panic!("no record with {patterns:?}"))
            .to_owned()
    };
    // Without datum A: the same items the A1 target case keeps.
    let without_a = ["Bore", "Bore projection", "Datum B", "Datum C", "Datum D"];

    let a1 = record(&["=PLACED_DATUM_TARGET_FEATURE(", "'A1'"]);
    let datum_a = record(&["=DATUM('", "'A');"]);
    let modified = record(&[
        "=DATUM_REFERENCE_COMPARTMENT(",
        "SIMPLE_DATUM_REFERENCE_MODIFIER",
    ]);
    let modified_id = modified.split('=').next().unwrap().to_owned();
    let cases: Vec<(String, &str)> = vec![
        // A placed target is product definitional.
        (
            mutate(&file, &a1, &a1.replace(",.T.,", ",.F.,")),
            "its placed target needs .T.",
        ),
        // A datum has an empty name and is not product definitional.
        (
            mutate(&file, &datum_a, &datum_a.replace(",.F.,", ",.T.,")),
            "needs an empty name and .F.",
        ),
        (
            mutate(
                &file,
                &datum_a,
                &datum_a.replacen("DATUM(''", "DATUM('A'", 1),
            ),
            "needs an empty name and .F.",
        ),
    ];
    for (broken, reason) in &cases {
        let (back, skipped) = read_back(broken, &report);
        assert_eq!(ids(&back), without_a, "{}", reasons(&skipped));
        assert!(
            reasons(&skipped).contains(reason),
            "{} lacks {reason}",
            reasons(&skipped)
        );
    }

    // A datum feature is product definitional: datum B goes, with what needs it.
    let feature_b = {
        let b = record(&["=DATUM('", "'B');"]);
        let b = b.split('=').next().unwrap().to_owned();
        let link = record(&["SHAPE_ASPECT_RELATIONSHIP(", &format!(",{b});")]);
        references(link.split_once('=').unwrap().1)[0].clone()
    };
    let feature_record = record(&[&format!("{feature_b}=DATUM_FEATURE(")]);
    let broken = mutate(
        &file,
        &feature_record,
        &feature_record.replace(",.T.)", ",.F.)"),
    );
    let (back, skipped) = read_back(&broken, &report);
    assert!(!ids(&back).contains(&"Datum B".to_owned()));
    assert!(ids(&back).contains(&"Datum C".to_owned()));
    assert!(reasons(&skipped).contains(&format!("its datum feature {feature_b} needs .T.")));

    // A datum reference compartment has an empty name, no description, no
    // identifier and is not product definitional. Only the composite whose
    // upper segment holds it goes.
    let id_attribute = format!("#999999=ID_ATTRIBUTE('D',{modified_id});\n");
    for broken in [
        mutate(&file, &modified, &modified.replacen(",.F.,", ",.T.,", 1)),
        mutate(
            &file,
            &modified,
            &modified.replacen("COMPARTMENT('',$,", "COMPARTMENT('D',$,", 1),
        ),
        mutate(
            &file,
            &modified,
            &modified.replacen("COMPARTMENT('',$,", "COMPARTMENT('','D',", 1),
        ),
        mutate(
            &file,
            "ENDSEC;\nEND-ISO",
            &format!("{id_attribute}ENDSEC;\nEND-ISO"),
        ),
    ] {
        let (back, skipped) = read_back(&broken, &report);
        let mut expected = ids(&all_items);
        expected.retain(|id| id != "Bolt pattern");
        assert_eq!(ids(&back), expected, "{}", reasons(&skipped));
        assert!(
            reasons(&skipped).contains(&format!("compartment {modified_id} needs an empty name")),
            "{}",
            reasons(&skipped)
        );
    }
}

#[test]
fn a_projected_zone_is_checked_over_its_own_height() {
    use catallax_boolean::{BoolOp, boolean};
    use catallax_shapes::{RVector3, Rotation};
    let mut model = Model::new();
    let mut dec = Decider::default();
    let at = |x, y, z| Pose::from_rational(Rotation::identity(), RVector3::new(q(x), q(y), q(z)));
    // A C-shaped block: a 10 mm gap between two plates joined at the back, and
    // a bore through the lower plate into the gap.
    use catallax_topology::block;
    let stock = block(&mut model, [q("20"), q("20"), q("20")], &Pose::identity()).unwrap();
    let gap = block(&mut model, [q("20"), q("15"), q("10")], &at("0", "0", "5")).unwrap();
    let bore = cylinder(&mut model, q("2"), q("5"), &at("10", "7", "0")).unwrap();
    let body = boolean(&mut model, stock, gap, BoolOp::Subtract, &mut dec)
        .unwrap()
        .body;
    let body = boolean(&mut model, body, bore, BoolOp::Subtract, &mut dec)
        .unwrap()
        .body;
    let off_axis = block(&mut model, [q("2"), q("2"), q("2")], &at("10", "30", "30")).unwrap();
    let mut separated = model.body(body).clone();
    separated.lumps.extend(&model.body(off_axis).lumps);
    let body = model.add_body(separated);
    let part = Bracket { model, body };
    let floor = part
        .planes(2, "5")
        .into_iter()
        .find(|face| part.model.face(*face).loops.len() == 2)
        .unwrap();
    let wall = part.round(2, "2", ["10", "7"]);
    let pmi = |height: &str| Pmi {
        items: vec![item(
            "Bore projection",
            ItemKind::Position(PositionControl {
                features: vec![wall.clone()],
                frame: PositionFrame::Single {
                    segment: segment("0.05", None, Vec::new()),
                    projected_zone: Some(ProjectedZone {
                        height_mm: q(height),
                        direction: v("0", "0", "1"),
                        entry_faces: vec![floor],
                    }),
                },
            }),
        )],
        ..Pmi::default()
    };

    // 2 mm into the gap, and exactly up to the far plate: the zone is free.
    // An endless ray would meet the far plate and refuse both.
    for height in ["2", "10"] {
        catallax_step::pmi::check_on_model(&part.model, part.body, &pmi(height), &mut dec)
            .unwrap_or_else(|error| panic!("{height} mm: {error:?}"));
    }
    let (file, report) = write(&part, &pmi("2"));
    let (back, skipped) = read_back(&file, &report);
    assert!(skipped.is_empty(), "{skipped:?}");
    let ItemKind::Position(control) = kind_of(&back, "Bore projection") else {
        panic!("expected the position control");
    };
    let PositionFrame::Single {
        projected_zone: Some(zone),
        ..
    } = &control.frame
    else {
        panic!("expected a projected zone");
    };
    assert_eq!(zone.height_mm, q("2"));

    // 12 mm reaches into the far plate: refused on the model, on export and
    // on import.
    let error = catallax_step::pmi::check_on_model(&part.model, part.body, &pmi("12"), &mut dec)
        .unwrap_err();
    assert!(
        error.detail.contains("axis meets a closed end"),
        "{error:?}"
    );
    assert!(
        try_write(&part, &pmi("12"))
            .unwrap_err()
            .contains("axis meets a closed end")
    );
    let all = records(&file);
    let definition = line_with(&file, &["=PROJECTED_ZONE_DEFINITION("]);
    let height = references(&all[&definition]).pop().unwrap();
    let measure = &all[&height];
    let long = mutate(
        &file,
        &format!("{height}={measure}"),
        &format!(
            "{height}={}",
            measure.replace("LENGTH_MEASURE(2.)", "LENGTH_MEASURE(12.)")
        ),
    );
    let read = read_step_pmi(long.as_bytes()).unwrap();
    assert!(read.pmi.items.is_empty());
    assert!(
        read.skipped
            .iter()
            .any(|item| item.reason.contains("axis meets a closed end")),
        "{:?}",
        read.skipped
    );

    // A boss standing up through the entry plane puts the zone's start in
    // material: refused however short the zone, though nothing lies within
    // its height.
    let mut model = Model::new();
    let at_z = |z| Pose::from_rational(Rotation::identity(), RVector3::new(q("0"), q("0"), q(z)));
    let stock = cylinder(&mut model, q("5"), q("6"), &Pose::identity()).unwrap();
    let outer_cut = cylinder(&mut model, q("4"), q("5"), &at_z("2")).unwrap();
    let boss = cylinder(&mut model, q("1"), q("7"), &at_z("1")).unwrap();
    let stock = boolean(&mut model, stock, boss, BoolOp::Unite, &mut dec)
        .unwrap()
        .body;
    let annular_cut = boolean(&mut model, outer_cut, boss, BoolOp::Subtract, &mut dec)
        .unwrap()
        .body;
    let body = boolean(&mut model, stock, annular_cut, BoolOp::Subtract, &mut dec)
        .unwrap()
        .body;
    let part = Bracket { model, body };
    let rim = part
        .planes(2, "6")
        .into_iter()
        .find(|face| part.model.face(*face).loops.len() == 2)
        .unwrap();
    let recess = part.round(2, "4", ["0", "0"]);
    let short = Pmi {
        items: vec![item(
            "Recess wall",
            ItemKind::Position(PositionControl {
                features: vec![recess],
                frame: PositionFrame::Single {
                    segment: segment("0.05", None, Vec::new()),
                    projected_zone: Some(ProjectedZone {
                        height_mm: q("1"),
                        direction: v("0", "0", "1"),
                        entry_faces: vec![rim],
                    }),
                },
            }),
        )],
        ..Pmi::default()
    };
    let error =
        catallax_step::pmi::check_on_model(&part.model, part.body, &short, &mut dec).unwrap_err();
    assert!(error.detail.contains("starts inside the part"), "{error:?}");
    assert!(
        try_write(&part, &short)
            .unwrap_err()
            .contains("starts inside the part")
    );
}

#[test]
fn projected_zone_checks_the_other_solids_in_its_part() {
    use catallax_boolean::{BoolOp, boolean};
    use catallax_shapes::{RVector3, Rotation};
    use catallax_topology::block;

    let mut model = Model::new();
    let mut dec = Decider::default();
    let at = |x, y, z| Pose::from_rational(Rotation::identity(), RVector3::new(q(x), q(y), q(z)));
    let stock = block(&mut model, [q("10"), q("10"), q("5")], &Pose::identity()).unwrap();
    let bore = cylinder(&mut model, q("2"), q("5"), &at("5", "5", "0")).unwrap();
    let plate = boolean(&mut model, stock, bore, BoolOp::Subtract, &mut dec)
        .unwrap()
        .body;
    let blocker = block(&mut model, [q("4"), q("4"), q("2")], &at("3", "3", "6")).unwrap();
    let mut body = model.body(plate).clone();
    body.lumps.extend(&model.body(blocker).lumps);
    let body = model.add_body(body);
    let part = Bracket { model, body };
    let [entry] = part.planes(2, "5").try_into().unwrap();
    let pmi = |height: &str| Pmi {
        items: vec![item(
            "Bore projection",
            ItemKind::Position(PositionControl {
                features: vec![part.round(2, "2", ["5", "5"])],
                frame: PositionFrame::Single {
                    segment: segment("0.05", None, Vec::new()),
                    projected_zone: Some(ProjectedZone {
                        height_mm: q(height),
                        direction: v("0", "0", "1"),
                        entry_faces: vec![entry],
                    }),
                },
            }),
        )],
        ..Pmi::default()
    };
    let (file, _) = write(&part, &pmi("0.5"));
    assert!(read_step_pmi(file.as_bytes()).unwrap().skipped.is_empty());
    let all = records(&file);
    let definition = line_with(&file, &["=PROJECTED_ZONE_DEFINITION("]);
    let measure = references(&all[&definition]).pop().unwrap();
    let old = &all[&measure];
    let blocked = mutate(
        &file,
        &format!("{measure}={old}"),
        &format!(
            "{measure}={}",
            old.replace("LENGTH_MEASURE(0.5)", "LENGTH_MEASURE(2.)")
        ),
    );
    assert_ne!(blocked, file);
    let read = read_step_pmi(blocked.as_bytes()).unwrap();
    assert!(read.pmi.items.is_empty(), "{:?}", read.pmi.items);
    assert!(
        read.skipped
            .iter()
            .any(|skip| skip.reason.contains("axis meets a closed end"))
    );
    let representation = line_with(&file, &["=ADVANCED_BREP_SHAPE_REPRESENTATION("]);
    let record = &all[&representation];
    let [_, _, blocker, context] = <[String; 4]>::try_from(references(record)).unwrap();
    let separate = mutate(
        &blocked,
        &format!("{representation}={record}"),
        &format!(
            "{representation}={}",
            record.replace(&format!(",{blocker})"), ")")
        ),
    );
    let separate = mutate(
        &separate,
        "ENDSEC;\nEND-ISO",
        &format!(
            "#900001=ADVANCED_BREP_SHAPE_REPRESENTATION('',({blocker}),{context});\n\
                  #900002=SHAPE_REPRESENTATION_RELATIONSHIP('','',#900001,{representation});\nENDSEC;\nEND-ISO"
        ),
    );
    let read = read_step_pmi(separate.as_bytes()).unwrap();
    assert!(read.pmi.items.is_empty(), "{:?}", read.pmi.items);
    assert!(
        read.skipped
            .iter()
            .any(|skip| skip.reason.contains("axis meets a closed end"))
    );
    let frame = &references(record)[0];
    let mapped = mutate(
        &separate,
        &format!("#900002=SHAPE_REPRESENTATION_RELATIONSHIP('','',#900001,{representation});"),
        &format!(
            "#900002=SHAPE_REPRESENTATION_RELATIONSHIP('','',#900003,{representation});\n\
             #900003=SHAPE_REPRESENTATION('',(#900004),{context});\n\
             #900004=MAPPED_ITEM('',#900005,{frame});\n\
             #900005=REPRESENTATION_MAP({frame},#900001);"
        ),
    );
    let read = read_step_pmi(mapped.as_bytes()).unwrap();
    assert!(read.pmi.items.is_empty());
    assert!(
        read.skipped
            .iter()
            .any(|item| item.reason.contains("unchecked geometry"))
    );

    let nested = mutate(
        &separate,
        &format!("#900002=SHAPE_REPRESENTATION_RELATIONSHIP('','',#900001,{representation});"),
        &format!(
            "#900002=SHAPE_REPRESENTATION_RELATIONSHIP('','',#900001,#900003);\n\
                 #900003=SHAPE_REPRESENTATION('',({}),{context});\n\
                 #900004=SHAPE_REPRESENTATION_RELATIONSHIP('','',#900003,{representation});",
            references(record)[0]
        ),
    );
    let read = read_step_pmi(nested.as_bytes()).unwrap();
    assert!(read.pmi.items.is_empty());
    assert!(
        read.skipped
            .iter()
            .any(|item| item.reason.contains("unchecked linked geometry"))
    );
    let start = &references(&all[&definition])[1];
    let usage = line_with(
        &file,
        &["=GEOMETRIC_ITEM_SPECIFIC_USAGE(", &format!(",{start},")],
    );
    let face = &references(&all[&usage])[2];
    let product = line_with(&file, &["=PRODUCT_DEFINITION('"]);
    let [frame, plate_solid, _, context] = <[String; 4]>::try_from(references(record)).unwrap();
    let foreign = mutate(
        &separate,
        &format!("{definition}={}", all[&definition]),
        &format!(
            "{definition}={}",
            all[&definition].replace(&format!(",{start},"), ",#900012,")
        ),
    );
    let foreign = mutate(
        &foreign,
        "ENDSEC;\nEND-ISO",
        &format!(
            "#900010=PRODUCT_DEFINITION_SHAPE('','',{product});\n\
                  #900011=ADVANCED_BREP_SHAPE_REPRESENTATION('',({frame},{plate_solid}),{context});\n\
                  #900012=SHAPE_ASPECT('','',#900010,.T.);\n\
                  #900013=SHAPE_DEFINITION_REPRESENTATION(#900010,#900011);\n\
                  #900014=GEOMETRIC_ITEM_SPECIFIC_USAGE('','',#900012,#900011,{face});\n\
                  ENDSEC;\nEND-ISO"
        ),
    );
    let read = read_step_pmi(foreign.as_bytes()).unwrap();
    assert!(read.pmi.items.is_empty());
    assert!(
        read.skipped
            .iter()
            .any(|item| item.reason.contains("start belongs to another part"))
    );
    catallax_step::pmi::check_on_model(&part.model, plate, &pmi("2"), &mut dec).unwrap();
    assert!(
        try_write(&part, &pmi("2"))
            .unwrap_err()
            .contains("axis meets a closed end")
    );
    assert!(
        catallax_step::pmi::check_on_model(&part.model, part.body, &pmi("2"), &mut dec)
            .unwrap_err()
            .detail
            .contains("axis meets a closed end")
    );
}

#[test]
fn projected_zone_checks_a_spherical_face_above_its_trim() {
    use catallax_boolean::{BoolOp, boolean};
    use catallax_shapes::{RVector3, Rotation};
    use catallax_topology::{block, sphere};

    let mut model = Model::new();
    let mut dec = Decider::default();
    let at = |x, y, z| Pose::from_rational(Rotation::identity(), RVector3::new(q(x), q(y), q(z)));
    let stock = block(&mut model, [q("20"), q("10"), q("5")], &Pose::identity()).unwrap();
    let dome = sphere(&mut model, q("5"), &at("14", "5", "9")).unwrap();
    let stock = boolean(&mut model, stock, dome, BoolOp::Unite, &mut dec)
        .unwrap()
        .body;
    let bore = cylinder(&mut model, q("0.5"), q("5"), &at("10", "5", "0")).unwrap();
    let body = boolean(&mut model, stock, bore, BoolOp::Subtract, &mut dec)
        .unwrap()
        .body;
    let part = Bracket { model, body };
    let [entry] = part.planes(2, "5").try_into().unwrap();
    let pmi = Pmi {
        items: vec![item(
            "Bore projection",
            ItemKind::Position(PositionControl {
                features: vec![part.round(2, "0.5", ["10", "5"])],
                frame: PositionFrame::Single {
                    segment: segment("0.05", None, Vec::new()),
                    projected_zone: Some(ProjectedZone {
                        height_mm: q("2"),
                        direction: v("0", "0", "1"),
                        entry_faces: vec![entry],
                    }),
                },
            }),
        )],
        ..Pmi::default()
    };
    let error =
        catallax_step::pmi::check_on_model(&part.model, part.body, &pmi, &mut dec).unwrap_err();
    assert!(error.detail.contains("projected zone"), "{error:?}");
    assert!(
        try_write(&part, &pmi)
            .unwrap_err()
            .contains("projected zone")
    );
}

#[test]
fn projected_zone_ignores_elliptic_trim_on_an_unreachable_plane() {
    use catallax_boolean::{BoolOp, boolean};
    use catallax_shapes::{RVector3, Rotation};
    use catallax_topology::block;

    let mut model = Model::new();
    let mut dec = Decider::default();
    let at = |rotation, x, y, z| Pose::from_rational(rotation, RVector3::new(q(x), q(y), q(z)));
    let stock = block(&mut model, [q("10"), q("10"), q("5")], &Pose::identity()).unwrap();
    let vertical = cylinder(
        &mut model,
        q("1"),
        q("5"),
        &at(Rotation::identity(), "5", "5", "0"),
    )
    .unwrap();
    let rotation = Rotation::from_quaternion(q("4"), q("0"), q("3"), q("0")).unwrap();
    let slanted = cylinder(
        &mut model,
        q("0.5"),
        q("14"),
        &at(rotation, "-1", "2", "0.7"),
    )
    .unwrap();
    let plate = boolean(&mut model, stock, vertical, BoolOp::Subtract, &mut dec)
        .unwrap()
        .body;
    let body = boolean(&mut model, plate, slanted, BoolOp::Subtract, &mut dec)
        .unwrap()
        .body;
    let part = Bracket { model, body };
    let [entry] = part.planes(2, "5").try_into().unwrap();
    let pmi = Pmi {
        items: vec![item(
            "Bore projection",
            ItemKind::Position(PositionControl {
                features: vec![part.round(2, "1", ["5", "5"])],
                frame: PositionFrame::Single {
                    segment: segment("0.05", None, Vec::new()),
                    projected_zone: Some(ProjectedZone {
                        height_mm: q("2"),
                        direction: v("0", "0", "1"),
                        entry_faces: vec![entry],
                    }),
                },
            }),
        )],
        ..Pmi::default()
    };
    catallax_step::pmi::check_on_model(&part.model, part.body, &pmi, &mut dec).unwrap();
    let (file, report) = write(&part, &pmi);
    assert!(
        records(&file)
            .values()
            .any(|record| record.starts_with("ELLIPSE("))
    );
    let (back, skipped) = read_back(&file, &report);
    assert!(skipped.is_empty(), "{skipped:?}");
    assert_eq!(back, pmi);
}
