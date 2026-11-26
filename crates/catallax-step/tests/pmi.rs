//! AP242 semantic PMI on the P1 bracket: what is written reads back the
//! same, on the same faces, with sizes measured from the geometry.
//!
//! Set `CATALLAX_PMI_OUT` to a directory to keep the written file, for a
//! check in another STEP reader.

use std::collections::BTreeMap;

use catallax_features::evaluate;
use catallax_import::read_step_pmi;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{HolePart, Recipe, Role};
use catallax_shapes::{Pose, Surface};
use catallax_step::pmi::{
    Anchor, At, Characteristic, Datum, Feature, GeometricTolerance, HoleHead, HoleSize, Item,
    ItemKind, Limits, Location, PartNotes, Pmi, Process, Size, SizeKind, SurfaceTexture, nominals,
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
fn the_bracket_file_counts_its_items_and_states_their_strings() {
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
    let counts: Vec<(String, String)> = properties[&shape].clone();
    let expected = [
        ("number of semantic pmi elements", "26"),
        ("number of dimensional locations", "14"),
        ("number of dimensional sizes", "4"),
        ("number of geometric tolerances", "4"),
        ("number of composite tolerances", "0"),
        ("number of datums", "3"),
        ("number of datum targets", "0"),
        ("number of annotations", "1"),
    ]
    .map(|(name, value)| {
        (
            "INTEGER_REPRESENTATION_ITEM".to_owned(),
            format!("'{name}',{value}"),
        )
    });
    assert_eq!(counts, expected);

    // Each datum and distance states its text; edition 1 has no property on a
    // size or a tolerance.
    let strings: Vec<(String, String)> = properties
        .iter()
        .filter(|(of, _)| **of != shape)
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
    let [association] = associations[..] else {
        panic!("one callout: {associations:?}");
    };
    assert!(association.starts_with("'Bore finish',''"));
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
            "its feature already carries a datum",
            "datum b is defined twice",
            "datum b is not defined",
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
    };
    (part, pmi)
}

#[test]
fn compound_and_blind_holes_keep_separate_faces_and_depths() {
    let (part, pmi) = hole_fixture();
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
