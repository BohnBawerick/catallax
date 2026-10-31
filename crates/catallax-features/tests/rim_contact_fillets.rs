//! P5 flanged housing: the junction fillet continues over the four countersink
//! mouths as torus pieces between the holes and one checked rolling-envelope
//! face per hole, for the authored radius 2, its edit to 5/2, and the translated,
//! scaled and rotated definitions, all inside one application allowance. Below
//! radius 1 the plane rail clears every mouth and one torus joins the housing.
//! At radius 1 the rail touches every mouth at one shared vertex.

use std::collections::{BTreeMap, BTreeSet};

use catallax_boolean::{Containment, classify_point};
use catallax_features::evaluate;
use catallax_measure::{
    Accuracy, blend_definition_properties, face_area, mass_properties, volume_area_centroid,
};
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, Real, with_operation_budget};
use catallax_recipe::{ElementId, Origin, Recipe, Role, StepId};
use catallax_shapes::{Curve, Point3, ROLLING_ENVELOPE_TARGET_NM, SplineSupport, Surface};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, VertexId};
use serde_json::{Value, json};

const APPLICATION_ALLOWANCE: u64 = 20_000_000;

fn canonical() -> Value {
    serde_json::from_str(include_str!("../../../parts/P5-flanged-housing.json")).unwrap()
}
fn step<'a>(recipe: &'a mut Value, id: &str) -> &'a mut Value {
    recipe["steps"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|step| step["id"] == id)
        .unwrap()
}
fn number(value: &Value) -> Rational {
    Rational::parse(value.as_str().unwrap()).unwrap()
}
fn rational(value: &Rational) -> Value {
    json!(value.to_string())
}
fn with_radius(radius: &str) -> Value {
    let mut recipe = canonical();
    step(&mut recipe, "filleted")["radius"] = json!(radius);
    recipe
}
// Definition points of the sketch plane, both tool poses and both axes.
fn placed_points(recipe: &mut Value, mut edit: impl FnMut(&mut Value)) {
    edit(&mut step(recipe, "section")["plane"]["translation"]);
    edit(&mut step(recipe, "housing")["axis"]["point"]);
    edit(&mut step(recipe, "drill")["pose"]["translation"]);
    edit(&mut step(recipe, "sink")["pose"]["translation"]);
    edit(&mut step(recipe, "hole-pattern")["layout"]["circular"]["axis"]["point"]);
}
fn translated(by: [i64; 3]) -> Value {
    let mut recipe = canonical();
    placed_points(&mut recipe, |point| {
        for (i, by) in by.iter().enumerate() {
            point[i] = rational(&(number(&point[i]) + Rational::from_i64(*by)));
        }
    });
    recipe
}
// A definition edit: every length of the part, including the fillet radius.
fn scaled(factor: &Rational) -> Value {
    let times = |value: &mut Value| *value = rational(&(number(value) * factor.clone()));
    let mut recipe = canonical();
    placed_points(&mut recipe, |point| {
        point.as_array_mut().unwrap().iter_mut().for_each(times);
    });
    let section = step(&mut recipe, "section");
    for segment in section["loops"][0].as_array_mut().unwrap() {
        segment["to"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .for_each(times);
    }
    for (id, lengths) in [
        ("drill", &["radius", "height"][..]),
        ("sink", &["base_radius", "top_radius", "height"]),
        ("filleted", &["radius"]),
    ] {
        for length in lengths {
            times(&mut step(&mut recipe, id)[*length]);
        }
    }
    recipe
}
// An exact rigid placement: rotation matrix over its denominator and the same
// rotation as a quaternion.
struct Placement {
    matrix: [[i64; 3]; 3],
    over: i64,
    quaternion: [i64; 4],
}
const RZ: Placement = Placement {
    matrix: [[3, -4, 0], [4, 3, 0], [0, 0, 5]],
    over: 5,
    quaternion: [2, 0, 0, 1],
};
const RQ: Placement = Placement {
    matrix: [[-20, 4, 22], [20, -10, 20], [10, 28, 4]],
    over: 30,
    quaternion: [1, 2, 3, 4],
};
const T: [i64; 3] = [7, 11, -13];
const NO_MOVE: [i64; 3] = [0, 0, 0];
impl Placement {
    fn point(&self, point: [Rational; 3], by: [i64; 3]) -> [Rational; 3] {
        std::array::from_fn(|i| {
            let row = (0..3).fold(Rational::zero(), |sum, j| {
                sum + point[j].clone() * ratio(self.matrix[i][j], self.over)
            });
            row + Rational::from_i64(by[i])
        })
    }
    // The placed definition: every frame, point and axis direction moves before
    // evaluation, so nothing rotates the finished body.
    fn recipe(&self, by: [i64; 3], radius: &str) -> Value {
        let point = |vector: &mut Value, by: [i64; 3]| {
            let old = std::array::from_fn(|i| number(&vector[i]));
            let new = self.point(old, by);
            *vector = json!(new.each_ref().map(rational));
        };
        let compose = |rotation: &mut Value| {
            let [a, b, c, d] = self.quaternion;
            let [w, x, y, z] =
                [0, 1, 2, 3].map(|i| rotation[i].as_str().unwrap().parse::<i64>().unwrap());
            *rotation = json!(
                [
                    a * w - b * x - c * y - d * z,
                    a * x + b * w + c * z - d * y,
                    a * y - b * z + c * w + d * x,
                    a * z + b * y - c * x + d * w,
                ]
                .map(|part| part.to_string())
            );
        };
        let mut recipe = with_radius(radius);
        placed_points(&mut recipe, |vector| point(vector, by));
        compose(&mut step(&mut recipe, "section")["plane"]["rotation"]);
        compose(&mut step(&mut recipe, "drill")["pose"]["rotation"]);
        compose(&mut step(&mut recipe, "sink")["pose"]["rotation"]);
        point(
            &mut step(&mut recipe, "housing")["axis"]["direction"],
            NO_MOVE,
        );
        point(
            &mut step(&mut recipe, "hole-pattern")["layout"]["circular"]["axis"]["direction"],
            NO_MOVE,
        );
        recipe
    }
}
fn parsed(recipe: &Value) -> Recipe {
    Recipe::parse(&recipe.to_string()).unwrap()
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: APPLICATION_ALLOWANCE,
        ..Default::default()
    })
}

// The exact contact law one row must retain, in millimetres.
struct Contact {
    origin: [Rational; 3],
    axis: [Rational; 3],
    radius: Rational,
    scale: Rational,
}

fn ratio(numerator: i64, denominator: i64) -> Rational {
    Rational::from_ratio(numerator, denominator).unwrap()
}

#[test]
fn p5_junction_continues_over_countersink_mouths() {
    let target = Rational::from_ratio(ROLLING_ENVELOPE_TARGET_NM, 1_000_000).unwrap();
    let whole = Rational::one();
    let up = [0, 1, 0].map(Rational::from_i64);
    let plain = |origin: [i64; 3], radius: Rational, scale: &Rational| Contact {
        origin: origin.map(Rational::from_i64),
        axis: up.clone(),
        radius,
        scale: scale.clone(),
    };
    let moved = |placement: &Placement, by: [i64; 3], radius: Rational| Contact {
        origin: placement.point([0, 5, 0].map(Rational::from_i64), by),
        axis: placement.point(up.clone(), NO_MOVE),
        radius,
        scale: whole.clone(),
    };
    let half = ratio(1, 2);
    let double = Rational::from_i64(2);
    let rows = [
        (canonical(), plain([0, 5, 0], ratio(2, 1), &whole)),
        (with_radius("5/2"), plain([0, 5, 0], ratio(5, 2), &whole)),
        (
            translated([7, 11, -13]),
            plain([7, 16, -13], ratio(2, 1), &whole),
        ),
        (
            scaled(&half),
            Contact {
                origin: [Rational::zero(), ratio(5, 2), Rational::zero()],
                axis: up.clone(),
                radius: Rational::one(),
                scale: half.clone(),
            },
        ),
        (scaled(&double), plain([0, 10, 0], ratio(4, 1), &double)),
        (RZ.recipe(NO_MOVE, "2"), moved(&RZ, NO_MOVE, ratio(2, 1))),
        (RZ.recipe(T, "5/2"), moved(&RZ, T, ratio(5, 2))),
        (RQ.recipe(NO_MOVE, "2"), moved(&RQ, NO_MOVE, ratio(2, 1))),
    ];
    let mut names: Option<BTreeSet<Origin>> = None;
    for (row, (recipe, contact)) in rows.iter().enumerate() {
        let recipe = parsed(recipe);
        let (all, work) = with_operation_budget(APPLICATION_ALLOWANCE, || {
            let mut dec = decider();
            let built = evaluate(&recipe, &mut dec).unwrap();
            let model = &built.model;
            assert_eq!(model.check(built.result), vec![]);
            let faces = model.faces_of_body(built.result);
            assert_eq!(faces.len(), 22);
            let mut kinds = BTreeMap::new();
            let mut all = BTreeSet::new();
            let mut blends = BTreeSet::new();
            for face in &faces {
                let surface = model.geom().surface(model.face(*face).surface);
                let blend = model
                    .face_origins(*face)
                    .iter()
                    .any(|o| matches!(o.role, Role::Blend { .. }));
                let kind = match surface {
                    Surface::Plane(_) => "plane",
                    Surface::Cylinder(_) => "cylinder",
                    Surface::Cone(_) => "cone",
                    Surface::Torus(_) => "torus",
                    Surface::BSpline(_) => "envelope",
                    other => panic!("P5 has no {other:?} face"),
                };
                *kinds.entry(kind).or_insert(0) += 1;
                match surface {
                    Surface::Torus(_) => assert!(blend),
                    Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) => {
                        assert!(blend);
                        assert!(envelope.bound().upper <= target);
                        let definition = envelope.definition();
                        let exact = |value: &catallax_numbers::Real| value.as_rational().unwrap();
                        let point = definition.origin();
                        assert_eq!(
                            [exact(&point.x), exact(&point.y), exact(&point.z)],
                            contact.origin
                        );
                        let axis = definition.axis();
                        assert_eq!(
                            [exact(&axis.x), exact(&axis.y), exact(&axis.z)],
                            contact.axis
                        );
                        assert_eq!(exact(definition.radius()), contact.radius);
                        let length = |mm: i64| Rational::from_i64(mm) * contact.scale.clone();
                        assert_eq!(exact(definition.body()), length(12));
                        assert_eq!(exact(definition.mouth()), length(4));
                        assert_eq!(exact(definition.offset()), length(17));
                    }
                    _ => assert!(!blend, "only the junction strip carries the blend origin"),
                }
                all.extend(model.face_origins(*face).iter().cloned());
                if blend {
                    blends.extend(model.face_origins(*face).iter().cloned());
                }
            }
            let expected = [
                ("cone", 4),
                ("cylinder", 7),
                ("envelope", 4),
                ("plane", 3),
                ("torus", 4),
            ];
            assert_eq!(kinds, BTreeMap::from(expected));
            // Every piece keeps the one blend name, unchanged by the radius edit.
            assert_eq!(blends.len(), 1);

            // Every boundary vertex lies exactly on every face it bounds, and the
            // shell is the genus-five surface of a bored housing with four holes.
            let mut vertices = BTreeSet::new();
            let mut edges = BTreeSet::new();
            let mut rings = BTreeSet::new();
            let mut loops = 0i64;
            let geom = model.geom();
            for face in &faces {
                let surface = model.face(*face).surface;
                for lp in &model.face(*face).loops {
                    loops += 1;
                    let Loop::Coedges(cs) = model.loop_(*lp) else {
                        continue;
                    };
                    for c in cs {
                        let edge = model.coedge(*c).edge;
                        edges.insert(edge);
                        let e = model.edge(edge);
                        if e.start.is_none() {
                            rings.insert(edge);
                        }
                        for v in [e.start, e.end].into_iter().flatten() {
                            vertices.insert(v);
                            let point = geom.position(model.vertex(v).point);
                            assert!(matches!(
                                geom.side_of(surface, &point, &mut dec),
                                Outcome::ProvenCoincident(_)
                            ));
                        }
                    }
                }
            }
            let euler = vertices.len() as i64 + rings.len() as i64 - edges.len() as i64
                + 2 * faces.len() as i64
                - loops;
            assert_eq!(euler, 2 * (1 - 5));
            all
        });
        assert!(!work.exhausted);
        // Face origins, the source names included, survive every row.
        if let Some(previous) = &names {
            assert_eq!(previous, &all, "row {row}");
        }
        names = Some(all);
    }
}

#[test]
fn a_ball_that_would_roll_onto_the_countersink_cone_refuses() {
    // Beyond 2 + sqrt(2) the contact leaves the mouth edge for the cone face,
    // a different continuation; it must refuse rather than build this one.
    let refusal = evaluate(&parsed(&with_radius("7/2")), &mut decider()).unwrap_err();
    let text = format!("{refusal:?}");
    assert!(text.contains("shapes.rolling-envelope"), "{text}");
}

#[test]
fn uncovered_contact_strata_keep_their_typed_refusal() {
    // Each case is a separate coverage gap, not impossible geometry: the moved
    // contact-change radius, and tubes too short to carry the cylinder rail, with
    // a touched mouth, a mouth on the rail, and none.
    let shortened = |mut recipe: Value, top: &str| {
        for segment in step(&mut recipe, "section")["loops"][0]
            .as_array_mut()
            .unwrap()
        {
            if segment["id"] == "body-wall" || segment["id"] == "body-top" {
                segment["to"][1] = json!(top);
            }
        }
        recipe
    };
    for (recipe, code) in [
        (with_radius("7/2"), "shapes.rolling-envelope"),
        (RZ.recipe(T, "7/2"), "shapes.rolling-envelope"),
        (
            shortened(RZ.recipe(T, "2"), "6"),
            "fillet.unsupported-overflow",
        ),
        (
            shortened(with_radius("1"), "11/2"),
            "fillet.unsupported-overflow",
        ),
        (
            shortened(with_radius("999/1000"), "11/2"),
            "fillet.unsupported-overflow",
        ),
    ] {
        let (refusal, work) = with_operation_budget(APPLICATION_ALLOWANCE, || {
            evaluate(&parsed(&recipe), &mut decider()).unwrap_err()
        });
        assert!(!work.exhausted);
        let text = format!("{refusal:?}");
        assert!(text.contains(code), "{text}");
    }
}

// One separated row: how its definition sits relative to the canonical part.
struct Clear {
    recipe: Value,
    // Fillet radius before scaling.
    radius: Rational,
    scale: Rational,
    pose: Option<(&'static Placement, [i64; 3])>,
}

impl Clear {
    fn place(&self, point: [Rational; 3]) -> [Rational; 3] {
        let point = point.map(|part| part * self.scale.clone());
        match self.pose {
            Some((placement, by)) => placement.point(point, by),
            None => point,
        }
    }
    fn point(&self, point: [Rational; 3]) -> Point3 {
        let [x, y, z] = self.place(point).map(Real::from_rational);
        Point3::new(x, y, z)
    }
    // Local rail radius.
    fn rail(&self) -> Rational {
        Rational::from_i64(12) + self.radius.clone()
    }
    // The enclosure of A*pi - B*pi^2 from the stated pi bounds, in the row's scale.
    fn volume(&self) -> [Rational; 2] {
        let [low, high] = pi();
        let r = self.radius.clone();
        let cube = r.clone() * r.clone() * r.clone();
        let a =
            ratio(11932, 3) + Rational::from_i64(24) * r.clone() * r.clone() + ratio(5, 3) * cube;
        let b = self.rail() * r.clone() * r / Rational::from_i64(2);
        let scaled = self.scale.clone() * self.scale.clone() * self.scale.clone();
        [
            (a.clone() * low.clone() - b.clone() * high.clone() * high.clone()) * scaled.clone(),
            (a * high - b * low.clone() * low) * scaled,
        ]
    }
    // The enclosure of (A + 48 sqrt(2)) pi + B pi^2, all coefficients positive.
    fn area(&self) -> [Rational; 2] {
        let [low, high] = pi();
        let root = sqrt_two();
        let r = self.radius.clone();
        let a = Rational::from_i64(2108)
            - Rational::from_i64(48) * r.clone()
            - Rational::from_i64(3) * r.clone() * r.clone();
        let b = self.rail() * r;
        let square = self.scale.clone() * self.scale.clone();
        let forty_eight = Rational::from_i64(48);
        [
            (low.clone() * (a.clone() + forty_eight.clone() * root[0].clone())
                + b.clone() * low.clone() * low)
                * square.clone(),
            (high.clone() * (a + forty_eight * root[1].clone()) + b * high.clone() * high) * square,
        ]
    }
}

fn pi() -> [Rational; 2] {
    let digits = "314159265358979323846264338327950288419716939937510";
    let denominator = format!("1{}", "0".repeat(50));
    let low = Rational::parse(&format!("{digits}/{denominator}")).unwrap();
    let high = low.clone() + Rational::parse(&format!("1/{denominator}")).unwrap();
    [low, high]
}

fn sqrt_two() -> [Rational; 2] {
    let digits = "141421356237309504880168872420969807856967187537694";
    let denominator = format!("1{}", "0".repeat(50));
    let low = Rational::parse(&format!("{digits}/{denominator}")).unwrap();
    let high = low.clone() + Rational::parse(&format!("1/{denominator}")).unwrap();
    assert!(low.clone() * low.clone() < Rational::from_i64(2));
    assert!(Rational::from_i64(2) < high.clone() * high.clone());
    [low, high]
}

fn clear_rows() -> [Clear; 4] {
    let near = ratio(999, 1000);
    let mut half = scaled(&ratio(1, 2));
    step(&mut half, "filleted")["radius"] = json!("999/2000");
    [
        Clear {
            recipe: with_radius("999/1000"),
            radius: near.clone(),
            scale: Rational::one(),
            pose: None,
        },
        Clear {
            recipe: with_radius("1/2"),
            radius: ratio(1, 2),
            scale: Rational::one(),
            pose: None,
        },
        Clear {
            recipe: RQ.recipe(T, "999/1000"),
            radius: near.clone(),
            scale: Rational::one(),
            pose: Some((&RQ, T)),
        },
        Clear {
            recipe: half,
            radius: near,
            scale: ratio(1, 2),
            pose: None,
        },
    ]
}

impl Clear {
    // Both native enclosures overlap the independent references.
    fn assert_properties(&self, model: &Model, body: BodyId) {
        let accuracy = Accuracy {
            relative: Rational::parse("1e-9").unwrap(),
            absolute: Rational::parse("1e-9").unwrap(),
        };
        let measured = volume_area_centroid(model, body, &accuracy, &mut decider())
            .unwrap()
            .value;
        let [low, high] = self.volume();
        let volume = measured.volume;
        assert!(volume.lo <= high, "{volume:?}");
        assert!(low <= volume.hi, "{volume:?}");
        assert!(volume.width() <= Rational::parse("1e-4").unwrap());
        let [low, high] = self.area();
        let area = measured.area;
        assert!(area.lo <= high, "{area:?}");
        assert!(low <= area.hi, "{area:?}");
        assert!(area.width() <= Rational::parse("1e-5").unwrap());
    }
    fn assert_sides(
        &self,
        model: &mut Model,
        body: BodyId,
        probes: Vec<([Rational; 3], Containment)>,
    ) {
        for (point, expected) in probes {
            match classify_point(model, body, &self.point(point.clone()), &mut decider()) {
                Outcome::Proven(found) => assert_eq!(found, expected, "{point:?}"),
                other => panic!("{point:?}: {other:?}"),
            }
        }
    }
}

// Strict hand-reasoned points of the definition, with the strip point's side.
fn strict_points(strip: Containment) -> [([Rational; 3], Containment); 6] {
    [
        (
            [ratio(10, 1), ratio(15, 1), ratio(0, 1)],
            Containment::Inside,
        ),
        (
            [ratio(0, 1), ratio(15, 1), ratio(0, 1)],
            Containment::Outside,
        ),
        (
            [ratio(17, 1), ratio(2, 1), ratio(0, 1)],
            Containment::Outside,
        ),
        (
            [ratio(0, 1), ratio(2, 1), ratio(20, 1)],
            Containment::Inside,
        ),
        ([ratio(15, 2), ratio(51, 10), ratio(10, 1)], strip),
        (
            [ratio(15, 2), ratio(26, 5), ratio(10, 1)],
            Containment::Outside,
        ),
    ]
}

fn kind(surface: &Surface) -> &'static str {
    match surface {
        Surface::Plane(_) => "plane",
        Surface::Cylinder(_) => "cylinder",
        Surface::Cone(_) => "cone",
        Surface::Torus(_) => "torus",
        other => panic!("a separated P5 has no {other:?} face"),
    }
}

fn blended(model: &Model, face: FaceId) -> bool {
    model
        .face_origins(face)
        .iter()
        .any(|o| matches!(o.role, Role::Blend { .. }))
}

fn edges_of(model: &Model, face: FaceId) -> Vec<EdgeId> {
    let mut edges = Vec::new();
    for lp in &model.face(face).loops {
        if let Loop::Coedges(cs) = model.loop_(*lp) {
            edges.extend(cs.iter().map(|c| model.coedge(*c).edge));
        }
    }
    edges
}

fn exact_point(point: &Point3) -> [Rational; 3] {
    [&point.x, &point.y, &point.z].map(|value| value.as_rational().unwrap())
}

#[test]
fn a_rail_clear_of_every_mouth_keeps_one_analytic_torus() {
    let mut fitted = decider();
    let fitted = evaluate(&parsed(&canonical()), &mut fitted).unwrap();
    let fitted_blends: BTreeSet<_> = fitted
        .model
        .faces_of_body(fitted.result)
        .into_iter()
        .filter(|face| blended(&fitted.model, *face))
        .flat_map(|face| fitted.model.face_origins(face).to_vec())
        .collect();
    for (row, clear) in clear_rows().iter().enumerate() {
        let recipe = parsed(&clear.recipe);
        let (built, work) = with_operation_budget(APPLICATION_ALLOWANCE, || {
            evaluate(&recipe, &mut decider()).unwrap()
        });
        assert!(!work.exhausted);
        let mut built = built;
        let model = &built.model;
        assert_eq!(model.check(built.result), vec![], "row {row}");
        let body = model.body(built.result);
        assert!(body.marks.is_empty(), "row {row}: {:?}", body.marks);
        assert_eq!(body.lumps.len(), 1);
        assert_eq!(model.lump(body.lumps[0]).shells.len(), 1);

        let faces = model.faces_of_body(built.result);
        let mut kinds = BTreeMap::new();
        for face in &faces {
            let surface = model.geom().surface(model.face(*face).surface);
            *kinds.entry(kind(surface)).or_insert(0) += 1;
            assert_eq!(matches!(surface, Surface::Torus(_)), blended(model, *face));
        }
        assert_eq!(
            kinds,
            BTreeMap::from([("cone", 4), ("cylinder", 7), ("plane", 3), ("torus", 1)]),
            "row {row}"
        );

        // A bored housing with four open holes: every edge is a ring, used
        // twice, so only faces and loops count towards genus five.
        let loops: i64 = faces
            .iter()
            .map(|face| model.face(*face).loops.len() as i64)
            .sum();
        assert_eq!(2 * faces.len() as i64 - loops, 2 * (1 - 5), "row {row}");

        // The blend joins its two selected supports along its exact rails.
        let torus = *faces.iter().find(|face| blended(model, **face)).unwrap();
        let rails = edges_of(model, torus);
        assert_eq!(rails.len(), 2);
        let mut found = BTreeMap::new();
        for rail in rails {
            let Curve::Circle(circle) = model.geom().curve(model.edge(rail).curve) else {
                panic!("a rail is a circle");
            };
            let neighbour: Vec<_> = faces
                .iter()
                .copied()
                .filter(|face| *face != torus && edges_of(model, *face).contains(&rail))
                .collect();
            assert_eq!(neighbour.len(), 1);
            let support = kind(model.geom().surface(model.face(neighbour[0]).surface));
            found.insert(
                support,
                (
                    exact_point(&circle.center),
                    circle.radius_sq.as_rational().unwrap(),
                    model.face_origins(neighbour[0]).to_vec(),
                ),
            );
        }
        let side = |segment: &str| Origin {
            step: StepId::new("housing").unwrap(),
            role: Role::Side {
                segment: ElementId::new(segment).unwrap(),
                path: None,
            },
        };
        let square = |value: Rational| value.clone() * value;
        let y = |height: Rational| [Rational::zero(), height, Rational::zero()];
        assert_eq!(
            found["plane"],
            (
                clear.place(y(Rational::from_i64(5))),
                square(clear.rail() * clear.scale.clone()),
                vec![side("flange-top")],
            ),
            "row {row}"
        );
        assert_eq!(
            found["cylinder"],
            (
                clear.place(y(Rational::from_i64(5) + clear.radius.clone())),
                square(Rational::from_i64(12) * clear.scale.clone()),
                vec![side("body-wall")],
            ),
            "row {row}"
        );

        // The one blend name matches the fitted pieces. Every drilled source
        // keeps its support, and the fillet adds nothing else.
        let blends: BTreeSet<_> = model.face_origins(torus).iter().cloned().collect();
        assert_eq!(blends.len(), 1);
        assert_eq!(blends, fitted_blends, "row {row}");
        let support = |model: &Model, face: FaceId| {
            (
                model.face_origins(face).to_vec(),
                format!("{:?}", model.geom().surface(model.face(face).surface)),
            )
        };
        let drilled = built.bodies[&StepId::new("drilled").unwrap()];
        let source: BTreeSet<_> = model
            .faces_of_body(drilled)
            .into_iter()
            .map(|face| support(model, face))
            .collect();
        let kept: BTreeSet<_> = faces
            .iter()
            .filter(|face| !blended(model, **face))
            .map(|face| support(model, *face))
            .collect();
        assert_eq!(source.len(), 14);
        assert_eq!(kept, source, "row {row}");

        clear.assert_properties(model, built.result);

        // Strict points of the definition, then one point on the plane rail.
        let near = clear.radius == ratio(999, 1000);
        let inside_strip = if near {
            Containment::Inside
        } else {
            Containment::Outside
        };
        let mut probes = strict_points(inside_strip).to_vec();
        probes.extend([
            (
                [ratio(25999, 2000), ratio(50001, 10000), ratio(0, 1)],
                Containment::Outside,
            ),
            (
                [ratio(25999, 2000), ratio(49999, 10000), ratio(0, 1)],
                Containment::Inside,
            ),
        ]);
        let result = built.result;
        clear.assert_sides(&mut built.model, result, probes);
        let rail = clear.rail();
        let on_rail = [
            rail.clone() * ratio(3, 5),
            Rational::from_i64(5),
            rail * ratio(4, 5),
        ];
        assert!(
            matches!(
                classify_point(
                    &mut built.model,
                    result,
                    &clear.point(on_rail),
                    &mut decider()
                ),
                Outcome::ProvenCoincident(_)
            ),
            "row {row}"
        );
    }
}

#[test]
fn volume_encloses_the_definition_and_the_fitted_carrier() {
    let half = Rational::from_ratio(1, 2).unwrap();
    for (recipe, reference, margin) in [
        (canonical(), 12558.681699123754, 0.02),
        (with_radius("5/2"), 12591.90560450861, 0.03),
        (translated([7, 11, -13]), 12558.681699123754, 0.02),
        (scaled(&half), 1569.8352123904692, 0.005),
        (RZ.recipe(NO_MOVE, "2"), 12558.681699123754, 0.02),
        (RZ.recipe(T, "5/2"), 12591.90560450861, 0.03),
        (RQ.recipe(NO_MOVE, "2"), 12558.681699123754, 0.02),
    ] {
        let mut dec = decider();
        let built = evaluate(&parsed(&recipe), &mut dec).unwrap();
        let accuracy = Accuracy {
            relative: Rational::parse("1e-9").unwrap(),
            absolute: Rational::parse("1e-9").unwrap(),
        };
        let definition =
            blend_definition_properties(&built.model, built.result, &accuracy, &mut dec)
                .unwrap()
                .value;
        let represented = volume_area_centroid(&built.model, built.result, &accuracy, &mut dec)
            .unwrap()
            .value;
        let middle = (&definition.volume.lo + &definition.volume.hi) / Rational::from_i64(2);
        assert!((middle.to_f64_nearest() - reference).abs() < 1e-5);
        let rounding = Rational::parse("1e-8").unwrap();
        let exact = Rational::parse(&reference.to_string()).unwrap();
        assert!(definition.volume.lo <= exact.clone() + rounding.clone());
        assert!(exact - rounding <= definition.volume.hi);
        assert!(definition.volume.width() <= Rational::parse("1e-4").unwrap());
        assert!(represented.volume.lo <= definition.volume.lo);
        assert!(definition.volume.hi <= represented.volume.hi);
        assert!(represented.volume.width() <= Rational::parse(&margin.to_string()).unwrap());
        let refusal = mass_properties(&built.model, built.result, &accuracy, &mut dec).unwrap_err();
        assert_eq!(refusal.code, "measure.unsupported-surface");
        // Exactly the five faces a rail trims refuse, before any work.
        let refused = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .filter(|face| face_area(&built.model, *face, &accuracy, &mut dec).is_err())
            .count();
        assert_eq!(refused, 5);
    }
}

fn tangent_rows() -> [Clear; 4] {
    let whole = Rational::one();
    let mut half = scaled(&ratio(1, 2));
    step(&mut half, "filleted")["radius"] = json!("1/2");
    [
        Clear {
            recipe: with_radius("1"),
            radius: whole.clone(),
            scale: whole.clone(),
            pose: None,
        },
        Clear {
            recipe: RQ.recipe(T, "1"),
            radius: whole.clone(),
            scale: whole.clone(),
            pose: Some((&RQ, T)),
        },
        Clear {
            recipe: RZ.recipe(T, "1"),
            radius: whole.clone(),
            scale: whole.clone(),
            pose: Some((&RZ, T)),
        },
        Clear {
            recipe: half,
            radius: whole,
            scale: ratio(1, 2),
            pose: None,
        },
    ]
}

// The faces around a vertex, one sector each: joined through the edges at the
// vertex they must form a single cycle, or the surface pinches there.
fn assert_one_cycle(model: &Model, faces: &[FaceId], vertex: VertexId) -> BTreeSet<FaceId> {
    let mut sectors = BTreeMap::<FaceId, Vec<FaceId>>::new();
    let mut joins = 0;
    let edges: BTreeSet<_> = faces
        .iter()
        .flat_map(|face| edges_of(model, *face))
        .collect();
    for edge in edges {
        let e = model.edge(edge);
        if e.start != Some(vertex) && e.end != Some(vertex) {
            continue;
        }
        let [a, b] = faces
            .iter()
            .copied()
            .filter(|face| edges_of(model, *face).contains(&edge))
            .collect::<Vec<_>>()[..]
        else {
            panic!("an edge at a contact joins two faces");
        };
        sectors.entry(a).or_default().push(b);
        sectors.entry(b).or_default().push(a);
        joins += 1;
    }
    assert_eq!(joins, sectors.len());
    let first = *sectors.keys().next().unwrap();
    let (mut previous, mut at, mut seen) = (first, sectors[&first][0], 1);
    while at != first {
        let around = &sectors[&at];
        assert_eq!(around.len(), 2);
        let next = if around[0] == previous {
            around[1]
        } else {
            around[0]
        };
        (previous, at, seen) = (at, next, seen + 1);
    }
    assert_eq!(seen, sectors.len(), "one vertex link");
    sectors.into_keys().collect()
}

fn circle_of(model: &Model, edge: EdgeId) -> ([Rational; 3], Rational) {
    let Curve::Circle(circle) = model.geom().curve(model.edge(edge).curve) else {
        panic!("a tangent arc is a circle");
    };
    (
        exact_point(&circle.center),
        circle.radius_sq.as_rational().unwrap(),
    )
}

#[test]
fn a_rail_touching_every_mouth_keeps_one_exact_torus() {
    let mut fitted = decider();
    let fitted = evaluate(&parsed(&canonical()), &mut fitted).unwrap();
    let fitted_blends: BTreeSet<_> = fitted
        .model
        .faces_of_body(fitted.result)
        .into_iter()
        .filter(|face| blended(&fitted.model, *face))
        .flat_map(|face| fitted.model.face_origins(face).to_vec())
        .collect();
    let accuracy = Accuracy {
        relative: Rational::parse("1e-9").unwrap(),
        absolute: Rational::parse("1e-9").unwrap(),
    };
    for (row, tangent) in tangent_rows().iter().enumerate() {
        let recipe = parsed(&tangent.recipe);
        let (mut built, work) = with_operation_budget(APPLICATION_ALLOWANCE, || {
            evaluate(&recipe, &mut decider()).unwrap()
        });
        assert!(!work.exhausted);
        let result = built.result;
        let model = &built.model;
        assert_eq!(model.check(result), vec![], "row {row}");
        let body = model.body(result);
        assert!(body.marks.is_empty(), "row {row}: {:?}", body.marks);
        assert_eq!(body.lumps.len(), 1);
        assert_eq!(model.lump(body.lumps[0]).shells.len(), 1);
        let faces = model.faces_of_body(result);
        let scale = |value: i64| Rational::from_i64(value) * tangent.scale.clone();
        let square = |value: Rational| value.clone() * value;
        let y = |height: i64| {
            [
                Rational::zero(),
                Rational::from_i64(height),
                Rational::zero(),
            ]
        };

        // Fifteen support groups, with the plane's seams explained: the flange
        // top alone is cut, into one main patch and a cap per contact. Each
        // patch keeps the group's source names, and the rest is unchanged.
        let support = |model: &Model, face: FaceId| {
            (
                model.face_origins(face).to_vec(),
                format!("{:?}", model.geom().surface(model.face(face).surface)),
            )
        };
        let drilled = built.bodies[&StepId::new("drilled").unwrap()];
        let source: BTreeSet<_> = model
            .faces_of_body(drilled)
            .into_iter()
            .map(|face| support(model, face))
            .collect();
        assert_eq!(source.len(), 14);
        let mut groups = BTreeMap::<_, Vec<FaceId>>::new();
        for face in &faces {
            groups.entry(support(model, *face)).or_default().push(*face);
        }
        assert_eq!(groups.len(), 15, "row {row}");
        let side = |segment: &str| Origin {
            step: StepId::new("housing").unwrap(),
            role: Role::Side {
                segment: ElementId::new(segment).unwrap(),
                path: None,
            },
        };
        let mut flange = Vec::new();
        let mut torus = None;
        for (key, patches) in &groups {
            if key.0 == vec![side("flange-top")] {
                flange = patches.clone();
            } else {
                assert_eq!(patches.len(), 1, "row {row}");
            }
            if blended(model, patches[0]) {
                assert_eq!(
                    key.0.iter().cloned().collect::<BTreeSet<_>>(),
                    fitted_blends,
                    "row {row}"
                );
                assert!(torus.replace(patches[0]).is_none());
            } else {
                assert!(source.contains(key), "row {row}");
            }
        }
        assert_eq!(flange.len(), 5, "row {row}");
        let torus = torus.unwrap();

        // Positive patches that cover the flange top with nothing missing or
        // doubled: the annulus 13 to 22 without four radius-4 discs is 251 pi.
        let [low, high] = pi();
        let area = |face: FaceId| {
            face_area(model, face, &accuracy, &mut decider())
                .unwrap()
                .value
        };
        let mut covered = [Rational::zero(), Rational::zero()];
        for patch in &flange {
            let area = area(*patch);
            assert!(area.lo > Rational::zero(), "row {row}");
            covered = [covered[0].clone() + area.lo, covered[1].clone() + area.hi];
        }
        let target = square(tangent.scale.clone()) * Rational::from_i64(251);
        assert!(covered[0] <= target.clone() * high, "row {row}");
        assert!(target * low <= covered[1], "row {row}");

        // The exact torus and both exact rails.
        let Surface::Torus(carrier) = model.geom().surface(model.face(torus).surface) else {
            panic!("the blend is a torus");
        };
        assert_eq!(
            exact_point(&carrier.center),
            tangent.place(y(6)),
            "row {row}"
        );
        assert_eq!(carrier.major.as_rational().unwrap(), scale(13));
        assert_eq!(carrier.minor.as_rational().unwrap(), scale(1));
        assert_eq!(model.face(torus).loops.len(), 2);
        let on_both = |edge: EdgeId| {
            faces
                .iter()
                .copied()
                .filter(|face| edges_of(model, *face).contains(&edge))
                .collect::<Vec<_>>()
        };
        let rails = edges_of(model, torus);
        let mut arcs = 0;
        for rail in &rails {
            let neighbour = on_both(*rail)
                .into_iter()
                .find(|face| *face != torus)
                .unwrap();
            let (center, radius_sq) = circle_of(model, *rail);
            if model.edge(*rail).start.is_none() {
                assert_eq!(center, tangent.place(y(6)), "row {row}");
                assert_eq!(radius_sq, square(scale(12)));
                assert_eq!(model.face_origins(neighbour), [side("body-wall")]);
            } else {
                assert_eq!(center, tangent.place(y(5)), "row {row}");
                assert_eq!(radius_sq, square(scale(13)));
                assert_eq!(model.face_origins(neighbour), [side("flange-top")]);
                arcs += 1;
            }
        }
        assert_eq!(arcs, rails.len() - 1);

        // Every vertex lies exactly on every face it bounds. Genus five comes
        // from the real face domains, with annuli and ring edges counted.
        let mut dec = decider();
        let mut vertices = BTreeSet::new();
        let mut edges = BTreeSet::new();
        let mut rings = BTreeSet::new();
        let mut loops = 0i64;
        let geom = model.geom();
        for face in &faces {
            let surface = model.face(*face).surface;
            loops += model.face(*face).loops.len() as i64;
            for edge in edges_of(model, *face) {
                edges.insert(edge);
                let e = model.edge(edge);
                if e.start.is_none() {
                    rings.insert(edge);
                }
                for v in [e.start, e.end].into_iter().flatten() {
                    vertices.insert(v);
                    let point = geom.position(model.vertex(v).point);
                    assert!(matches!(
                        geom.side_of(surface, &point, &mut dec),
                        Outcome::ProvenCoincident(_)
                    ));
                }
            }
        }
        let euler = vertices.len() as i64 + rings.len() as i64 - edges.len() as i64
            + 2 * faces.len() as i64
            - loops;
        assert_eq!(euler, 2 * (1 - 5), "row {row}");

        // Each contact is one shared vertex of the rail and its mouth, and the
        // surface around it is one cycle: torus, flange top, cone, flange top.
        let contacts = [[13, 0], [0, 13], [-13, 0], [0, -13]];
        for [x, z] in contacts {
            let local = [x, 5, z].map(Rational::from_i64);
            let exact = tangent.place(local);
            let matching: Vec<_> = vertices
                .iter()
                .copied()
                .filter(|v| exact_point(&model.geom().position(model.vertex(*v).point)) == exact)
                .collect();
            let [contact] = matching[..] else {
                panic!("row {row}: one vertex at {exact:?}");
            };
            let around = assert_one_cycle(model, &faces, contact);
            assert_eq!(around.len(), 4, "row {row}");
            assert!(around.contains(&torus));
            let planes = around.iter().filter(|face| flange.contains(face)).count();
            assert_eq!(planes, 2, "row {row}");
            assert_eq!(
                around
                    .iter()
                    .filter(|face| kind(model.geom().surface(model.face(**face).surface)) == "cone")
                    .count(),
                1,
                "row {row}"
            );
            let [mx, mz] = [x, z].map(|value| value * 17 / 13);
            let mouth = tangent.place([mx, 5, mz].map(Rational::from_i64));
            let mut carried = BTreeSet::new();
            for edge in edges.iter().copied() {
                let e = model.edge(edge);
                if e.start == Some(contact) || e.end == Some(contact) {
                    carried.insert(circle_of(model, edge));
                }
            }
            assert_eq!(
                carried,
                BTreeSet::from([
                    (tangent.place(y(5)), square(scale(13))),
                    (mouth, square(scale(4))),
                ]),
                "row {row}"
            );
        }

        tangent.assert_properties(model, result);

        // Strict points, then the four around the first contact: a margin of
        // over 50 nm on each side of the plane, the rail and the mouth.
        let mut probes = strict_points(Containment::Inside).to_vec();
        probes.extend([
            (
                [ratio(1299, 100), ratio(5, 1), ratio(0, 1)],
                Containment::Inside,
            ),
            (
                [ratio(1301, 100), ratio(5, 1), ratio(0, 1)],
                Containment::Outside,
            ),
            (
                [ratio(13, 1), ratio(499, 100), ratio(0, 1)],
                Containment::Inside,
            ),
            (
                [ratio(13, 1), ratio(501, 100), ratio(0, 1)],
                Containment::Outside,
            ),
        ]);
        tangent.assert_sides(&mut built.model, result, probes);

        // Exact boundary: all four contacts, and a cusp point the plane keeps
        // beyond the rail and the mouth, so neither was widened.
        let mut boundary: Vec<_> = contacts
            .iter()
            .map(|[x, z]| [x, &5, z].map(|value| Rational::from_i64(*value)))
            .collect();
        boundary.push([ratio(13, 1), ratio(5, 1), ratio(1, 100)]);
        for point in boundary {
            assert!(
                matches!(
                    classify_point(
                        &mut built.model,
                        result,
                        &tangent.point(point.clone()),
                        &mut decider()
                    ),
                    Outcome::ProvenCoincident(_)
                ),
                "row {row}: {point:?}"
            );
        }
    }
}

// The tangent housing with a far sphere and a far separated housing beside it.
fn with_separate_bodies() -> Value {
    let mut recipe = with_radius("1");
    let mut far = translated([100, 0, 0]);
    step(&mut far, "filleted")["radius"] = json!("1/2");
    let mut steps = far["steps"].to_string();
    for id in [
        "section",
        "housing",
        "drill",
        "sink",
        "countersink",
        "hole-pattern",
        "drilled",
        "filleted",
    ] {
        steps = steps.replace(&format!("\"{id}\""), &format!("\"far-{id}\""));
    }
    let list = recipe["steps"].as_array_mut().unwrap();
    list.extend(serde_json::from_str::<Vec<Value>>(&steps).unwrap());
    list.push(json!({
        "id": "ball",
        "op": "sphere",
        "radius": "3",
        "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "100"]},
    }));
    recipe
}

#[test]
fn every_body_written_beside_a_tangent_ring_keeps_its_seams() {
    let built = evaluate(&parsed(&with_separate_bodies()), &mut decider()).unwrap();
    // Seven cylinders, four cones and the torus of each housing carry one seam.
    for (other, seams) in [("ball", 12), ("far-filleted", 24)] {
        let body = built.bodies[&StepId::new(other).unwrap()];
        let mut file = Vec::new();
        catallax_step::write_step(
            &built.model,
            &[built.result, body],
            &Default::default(),
            &mut decider(),
            &mut file,
        )
        .unwrap();
        let file = String::from_utf8(file).unwrap();
        assert_eq!(file.matches("SEAM_CURVE").count(), seams, "{other}");
        assert_eq!(file.matches("MANIFOLD_SOLID_BREP").count(), 2, "{other}");
    }
}
