//! P5 flanged housing: the junction fillet continues over the four countersink
//! mouths as torus pieces between the holes and one checked rolling-envelope
//! face per hole, for the authored radius 2, its edit to 5/2, and the translated,
//! scaled and rotated definitions, all inside one application allowance. Below
//! radius 1 the plane rail clears every mouth and one torus joins the housing.
//! At radius 1 the rail touches every mouth at one shared vertex.

use std::collections::{BTreeMap, BTreeSet};

use catallax_boolean::{Containment, classify_point, face_contains};
use catallax_features::evaluate;
use catallax_intersect::{Contact as Touch, CurveSurface, intersect_curve_surface};
use catallax_measure::{
    Accuracy, blend_definition_properties, face_area, mass_properties, volume_area_centroid,
};
use catallax_numbers::{
    BudgetPolicy, Decider, Outcome, Rational, Real, Sign, with_operation_budget,
};
use catallax_recipe::{ElementId, Origin, Recipe, Role, StepId};
use catallax_shapes::{
    BSplineCurve, Cone, Curve, CurveId, IntersectionCurve, Line, Meet, MeetOf, Plane, Point3,
    PointDef, Position, ROLLING_ENVELOPE_TARGET_NM, RollingEnvelope, SplineSupport, Surface,
    Vector3,
};
use catallax_topology::{
    BodyId, Coedge, Edge, EdgeId, Face, FaceId, Loop, Lump, Model, Sense, Shell, Vertex, VertexId,
    envelope_retains_chart, envelope_rows,
};
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
        (RQ.recipe(T, "5/2"), moved(&RQ, T, ratio(5, 2))),
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

// Hand-derived points of the fitted body. Only the last depends on the radius:
// b(2) < 4/5 < b(5/2) for the torus height b(r) above the flange at radius 25/2.
fn fitted_witnesses(edit: Containment) -> [([Rational; 3], Containment); 6] {
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
        (
            [ratio(15, 2), ratio(11, 2), ratio(10, 1)],
            Containment::Inside,
        ),
        ([ratio(15, 2), ratio(29, 5), ratio(10, 1)], edit),
    ]
}
// V8 J2 strict points past the fit bound: E2, O and Q-. The mouth point
// N = (853/65, 5, 64/65) lies strictly between that mouth's transition ends.
fn reaching_witnesses(edit: Containment) -> [([Rational; 3], Containment); 3] {
    [
        ([ratio(25, 2), ratio(109, 20), ratio(0, 1)], edit),
        (
            [ratio(25, 2), ratio(23, 4), ratio(0, 1)],
            Containment::Outside,
        ),
        (
            [ratio(3347, 260), ratio(19, 4), ratio(64, 65)],
            Containment::Inside,
        ),
    ]
}
fn fitted_supports(model: &Model, body: BodyId) -> Vec<(FaceId, RollingEnvelope)> {
    let support = |face: &FaceId| match model.geom().surface(model.face(*face).surface) {
        Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) => {
            Some((*face, envelope.clone()))
        }
        _ => None,
    };
    model
        .faces_of_body(body)
        .iter()
        .filter_map(support)
        .collect()
}
// What a query must leave alone. Auxiliary points and lines may join the store.
fn query_snapshot(model: &Model, body: BodyId) -> String {
    let mut text = format!("{:?}\n", model.body(body));
    for face in model.faces_of_body(body) {
        let record = model.face(face);
        text += &format!("{face:?} {record:?} {:?}\n", model.face_origins(face));
        for lp in &record.loops {
            text += &format!("{lp:?} {:?}\n", model.loop_(*lp));
        }
    }
    for edge in model.edges_of_body(body) {
        let record = model.edge(edge);
        text += &format!("{edge:?} {record:?} {:?}", model.edge_constructions(edge));
        for coedge in model.coedges_of_edge(edge) {
            text += &format!(" {coedge:?} {:?}", model.coedge(coedge));
        }
        for vertex in [record.start, record.end].into_iter().flatten() {
            text += &format!(" {vertex:?} {:?}", model.vertex(vertex));
        }
        text += "\n";
    }
    text
}

#[test]
fn fitted_housing_answers_native_containment() {
    let halved = ratio(1, 2);
    let rows = [
        (canonical(), ratio(2, 1), None, Containment::Outside),
        (with_radius("5/2"), ratio(5, 2), None, Containment::Inside),
        (
            RQ.recipe(T, "5/2"),
            ratio(5, 2),
            Some((&RQ, T)),
            Containment::Inside,
        ),
        (scaled(&halved), ratio(2, 1), None, Containment::Outside),
    ];
    for (row, (recipe, radius, pose, edit)) in rows.into_iter().enumerate() {
        let scale = if row == 3 {
            halved.clone()
        } else {
            Rational::one()
        };
        let parsed_recipe = parsed(&recipe);
        let place = Clear {
            recipe,
            radius,
            scale,
            pose,
        };
        let mut built = evaluate(&parsed_recipe, &mut decider()).unwrap();
        let body = built.result;
        assert_eq!(built.model.check(body), vec![], "row {row}");
        assert_eq!(built.model.faces_of_body(body).len(), 22, "row {row}");
        assert_eq!(fitted_supports(&built.model, body).len(), 4, "row {row}");
        let before = query_snapshot(&built.model, body);
        for (point, expected) in fitted_witnesses(edit) {
            let point = place.point(point);
            let (answer, work) = with_operation_budget(APPLICATION_ALLOWANCE, || {
                classify_point(&mut built.model, body, &point, &mut decider())
            });
            assert!(!work.exhausted, "row {row}");
            assert!(
                matches!(answer, Outcome::Proven(found) if found == expected),
                "row {row}: {point:?}: {answer:?}"
            );
        }
        // Strict points whose rays really reach the +x patch. E2 is the edit.
        for (point, expected) in reaching_witnesses(edit) {
            let point = place.point(point);
            let (answer, work) = with_operation_budget(APPLICATION_ALLOWANCE, || {
                classify_point(&mut built.model, body, &point, &mut decider())
            });
            assert!(!work.exhausted, "row {row}");
            assert!(
                matches!(answer, Outcome::Proven(found) if found == expected),
                "row {row}: {point:?}: {answer:?}"
            );
        }
        // The closest point of the first countersink mouth stays exact incidence.
        let mouth = place.point([ratio(13, 1), ratio(5, 1), ratio(0, 1)]);
        let (answer, work) = with_operation_budget(APPLICATION_ALLOWANCE, || {
            classify_point(&mut built.model, body, &mouth, &mut decider())
        });
        assert!(!work.exhausted, "row {row}");
        assert!(
            matches!(answer, Outcome::ProvenCoincident(_)),
            "row {row}: {answer:?}"
        );
        assert_eq!(before, query_snapshot(&built.model, body), "row {row}");
    }
}

#[test]
fn finite_queries_prove_a_miss_and_keep_a_real_contact() {
    let mut built = evaluate(&parsed(&canonical()), &mut decider()).unwrap();
    let body = built.result;
    let supports = fitted_supports(&built.model, body);
    assert_eq!(supports.len(), 4);
    let point = |coordinates: [Rational; 3]| {
        let [x, y, z] = coordinates.map(Real::from_rational);
        Point3::new(x, y, z)
    };
    let tube = point([ratio(10, 1), ratio(15, 1), ratio(0, 1)]);
    let torus = point([ratio(15, 2), ratio(11, 2), ratio(10, 1)]);
    let mouth = point([ratio(13, 1), ratio(5, 1), ratio(0, 1)]);
    let up = || Vector3::new(Real::from_i64(0), Real::from_i64(1), Real::from_i64(0));
    let mut dec = decider();
    let mut line_through = |model: &mut Model, through: &Point3| {
        let line = Line::new(through.clone(), up(), &mut dec).unwrap();
        model.geom_mut().add_curve(Curve::Line(line))
    };
    let tube_line = line_through(&mut built.model, &tube);
    let torus_line = line_through(&mut built.model, &torus);
    let mouth_line = line_through(&mut built.model, &mouth);
    for (face, _) in &supports {
        assert!(!face_contains(&mut built.model, *face, &tube, &mut decider()).unwrap());
        let surface = built.model.face(*face).surface;
        for line in [tube_line, torus_line] {
            let answer =
                intersect_curve_surface(built.model.geom_mut(), line, surface, &mut decider())
                    .unwrap();
            assert!(matches!(answer.value, CurveSurface::Disjoint));
            assert!(answer.marks.is_empty());
        }
    }
    // The vertical line through the +x mouth point is a tested J2 contact now.
    let surface = supports
        .iter()
        .find(|(_, envelope)| toward_x(envelope))
        .map(|(face, _)| built.model.face(*face).surface)
        .unwrap();
    let answer =
        intersect_curve_surface(built.model.geom_mut(), mouth_line, surface, &mut decider())
            .unwrap();
    assert!(answer.marks.is_empty());
    let CurveSurface::Meets(hits) = answer.value else {
        panic!("{:?}", answer.value);
    };
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].contact, Touch::Crossing);
    let gap = built.model.geom().position(hits[0].point).minus(&mouth);
    assert!(decider().sign(&gap.norm_squared()).is_coincident());
}

fn toward_x(envelope: &RollingEnvelope) -> bool {
    let toward = envelope.definition().toward();
    [&toward.x, &toward.y, &toward.z].map(|part| part.as_rational().unwrap())
        == [ratio(1, 1), ratio(0, 1), ratio(0, 1)]
}

// One V8 J2 row: a definition, and how its local points and directions land.
struct LineRow {
    recipe: Value,
    scale: Rational,
    pose: Option<(&'static Placement, [i64; 3])>,
    // A's independent height window about the flange mouth, canonical frame.
    window: [Rational; 2],
}
impl LineRow {
    fn point(&self, local: [Rational; 3]) -> [Rational; 3] {
        let scaled = local.map(|part| part * self.scale.clone());
        match self.pose {
            Some((placement, by)) => placement.point(scaled, by),
            None => scaled,
        }
    }
    fn direction(&self, local: [Rational; 3]) -> [Rational; 3] {
        match self.pose {
            Some((placement, _)) => placement.point(local, NO_MOVE),
            None => local,
        }
    }
    // Back to the canonical frame, for a hit's oracle window.
    fn height(&self, world: &[catallax_numbers::Ball; 3]) -> catallax_numbers::Ball {
        match self.pose {
            Some((placement, by)) => {
                let shifted: [_; 3] = std::array::from_fn(|i| {
                    &world[i]
                        - &catallax_numbers::Ball::from_rational(&Rational::from_i64(by[i]), 128)
                });
                (0..3).fold(
                    catallax_numbers::Ball::from_rational(&Rational::zero(), 128),
                    |sum, j| {
                        let cell = ratio(placement.matrix[j][1], placement.over);
                        sum + &catallax_numbers::Ball::from_rational(&cell, 128) * &shifted[j]
                    },
                )
            }
            None => {
                world[1].clone()
                    * catallax_numbers::Ball::from_rational(
                        &(Rational::one() / self.scale.clone()),
                        128,
                    )
            }
        }
    }
}
fn line_rows() -> [LineRow; 4] {
    let narrow = [ratio(27, 5), ratio(271, 50)];
    let wide = [ratio(137, 25), ratio(138, 25)];
    [
        LineRow {
            recipe: canonical(),
            scale: Rational::one(),
            pose: None,
            window: narrow.clone(),
        },
        LineRow {
            recipe: with_radius("5/2"),
            scale: Rational::one(),
            pose: None,
            window: wide.clone(),
        },
        LineRow {
            recipe: RQ.recipe(T, "5/2"),
            scale: Rational::one(),
            pose: Some((&RQ, T)),
            window: wide,
        },
        LineRow {
            recipe: scaled(&ratio(1, 2)),
            scale: ratio(1, 2),
            pose: None,
            window: narrow,
        },
    ]
}
fn rational_point(coordinates: [Rational; 3]) -> Point3 {
    let [x, y, z] = coordinates.map(Real::from_rational);
    Point3::new(x, y, z)
}
fn rational_vector(coordinates: [Rational; 3]) -> Vector3 {
    let [x, y, z] = coordinates.map(Real::from_rational);
    Vector3::new(x, y, z)
}
// The support whose retained frame points at this row's +x mouth.
fn mouth_support(model: &Model, body: BodyId, row: &LineRow) -> (FaceId, RollingEnvelope) {
    let expected = row.direction([ratio(1, 1), ratio(0, 1), ratio(0, 1)]);
    let mut matching = fitted_supports(model, body)
        .into_iter()
        .filter(|(_, envelope)| {
            let toward = envelope.definition().toward();
            let found = [&toward.x, &toward.y, &toward.z].map(|part| part.as_rational().unwrap());
            found == expected
        });
    let support = matching.next().unwrap();
    assert!(matching.next().is_none());
    support
}

// A, BM and BN meet the +x support once. A meets inside it, BM and BN on its
// mouth at exactly M and N, and the other three supports stay clear. X misses
// every support, and T, tangent along the rail, stays a typed refusal.
#[test]
fn finite_lines_cross_the_fitted_patch_exactly_once() {
    let along_y = || [ratio(0, 1), ratio(1, 1), ratio(0, 1)];
    let m = || [ratio(13, 1), ratio(5, 1), ratio(0, 1)];
    let n = || [ratio(853, 65), ratio(5, 1), ratio(64, 65)];
    for (index, row) in line_rows().iter().enumerate() {
        let mut built = evaluate(&parsed(&row.recipe), &mut decider()).unwrap();
        let body = built.result;
        let (mouth_face, mouth) = mouth_support(&built.model, body, row);
        let lines = [
            ([ratio(25, 2), ratio(0, 1), ratio(0, 1)], along_y(), None),
            (
                [ratio(13, 1), ratio(0, 1), ratio(0, 1)],
                along_y(),
                Some(m()),
            ),
            (n(), [ratio(1, 1), ratio(1, 1), ratio(0, 1)], Some(n())),
        ];
        for (line, (through, direction, closure)) in lines.into_iter().enumerate() {
            let label = format!("row {index} line {line}");
            let through = rational_point(row.point(through));
            let direction = rational_vector(row.direction(direction));
            let carrier = Line::new(through, direction, &mut decider()).unwrap();
            let curve = built.model.geom_mut().add_curve(Curve::Line(carrier));
            for (face, _) in fitted_supports(&built.model, body) {
                let surface = built.model.face(face).surface;
                let (answer, work) = with_operation_budget(APPLICATION_ALLOWANCE, || {
                    intersect_curve_surface(built.model.geom_mut(), curve, surface, &mut decider())
                });
                assert!(!work.exhausted, "{label}");
                let answer = answer.unwrap();
                assert!(answer.marks.is_empty(), "{label}");
                if face != mouth_face {
                    assert!(matches!(answer.value, CurveSurface::Disjoint), "{label}");
                    continue;
                }
                let CurveSurface::Meets(hits) = answer.value else {
                    panic!("{label}: {:?}", answer.value);
                };
                assert_eq!(hits.len(), 1, "{label}");
                assert_eq!(hits[0].contact, Touch::Crossing, "{label}");
                let position = built.model.geom().position(hits[0].point);
                // Every hit is on the represented carrier, not only on its target.
                let width = ratio(1, 1_000_000);
                let mut enclose = decider();
                let at = built
                    .model
                    .geom()
                    .enclose_point(hits[0].point, &ratio(1, 4_000_000), &mut enclose)
                    .unwrap();
                let support = SplineSupport::RollingEnvelope(mouth.clone());
                let [u, v] = support.parameters_of(&at, &enclose).unwrap().unwrap();
                let image = mouth.carrier().evaluate(&u, &v, &enclose).unwrap();
                assert!(
                    (0..3).all(|i| image[i].lo() <= at[i].hi() && at[i].lo() <= image[i].hi()),
                    "{label}"
                );
                let on_row = mouth.on_boundary(&position, &mut decider()).unwrap();
                assert_eq!(on_row, closure.is_some(), "{label}");
                if let Some(expected) = &closure {
                    let expected = rational_point(row.point(expected.clone()));
                    let gap = position.minus(&expected).norm_squared();
                    assert!(
                        matches!(decider().sign(&gap), Outcome::ProvenCoincident(_)),
                        "{label}"
                    );
                } else {
                    let height = row.height(&at);
                    assert!(
                        row.window[0] < height.lo() && height.hi() < row.window[1],
                        "{label}: {height:?}"
                    );
                    assert!(height.hi() - height.lo() <= width, "{label}");
                }
            }
        }
        let vertical = |model: &mut Model, x: Rational| {
            let through = rational_point(row.point([x, ratio(0, 1), ratio(0, 1)]));
            let direction = rational_vector(row.direction(along_y()));
            let line = Line::new(through, direction, &mut decider()).unwrap();
            model.geom_mut().add_curve(Curve::Line(line))
        };
        let beyond = vertical(&mut built.model, ratio(27, 2));
        let tangent = vertical(&mut built.model, ratio(12, 1));
        for (face, _) in fitted_supports(&built.model, body) {
            let surface = built.model.face(face).surface;
            let answer =
                intersect_curve_surface(built.model.geom_mut(), beyond, surface, &mut decider())
                    .unwrap();
            assert!(
                matches!(answer.value, CurveSurface::Disjoint),
                "row {index}"
            );
            assert!(answer.marks.is_empty(), "row {index}");
        }
        let surface = built.model.face(mouth_face).surface;
        let refusal =
            intersect_curve_surface(built.model.geom_mut(), tangent, surface, &mut decider())
                .unwrap_err();
        assert_eq!(refusal.code, "intersect.parametric-surface", "row {index}");
    }
}

// The mouth circle makes the section at x=13 lowest at z=0, so a line across it
// at a height above M crosses the patch twice, once on each side of the mirror.
#[test]
fn a_line_across_the_patch_gets_both_crossings_in_order() {
    let row = &line_rows()[0];
    let mut built = evaluate(&parsed(&row.recipe), &mut decider()).unwrap();
    let body = built.result;
    let (face, _) = mouth_support(&built.model, body, row);
    let surface = built.model.face(face).surface;
    let through = rational_point([ratio(13, 1), ratio(51, 10), ratio(-2, 1)]);
    let along = rational_vector([ratio(0, 1), ratio(0, 1), ratio(1, 1)]);
    let line = Line::new(through, along, &mut decider()).unwrap();
    let curve = built.model.geom_mut().add_curve(Curve::Line(line));
    let answer =
        intersect_curve_surface(built.model.geom_mut(), curve, surface, &mut decider()).unwrap();
    let CurveSurface::Meets(hits) = answer.value else {
        panic!("{:?}", answer.value);
    };
    assert_eq!(hits.len(), 2);
    assert!(hits.iter().all(|hit| hit.contact == Touch::Crossing));
    let z: Vec<_> = hits
        .iter()
        .map(|hit| {
            built
                .model
                .geom()
                .enclose_point(hit.point, &ratio(1, 1_000_000), &mut decider())
                .unwrap()[2]
                .clone()
        })
        .collect();
    let zero = Rational::zero();
    assert!(z[0].hi() < zero && zero < z[1].lo());
    let balance = (&z[0] + &z[1]).hi().abs().max((&z[0] + &z[1]).lo().abs());
    assert!(balance < ratio(1, 1000));
}

// The plane x=13 cuts the cone with this apex, axis (0,1,1) and slope 1 along
// two rays. The generator through `witness` runs along y from the apex at
// (13,apex,0), on the supporting line of BM, and stops there.
fn generator(model: &mut Model, apex: i64, witness: i64) -> CurveId {
    let geom = model.geom_mut();
    let plane = geom.add_surface(Surface::Plane(Plane {
        origin: rational_point([ratio(13, 1), ratio(0, 1), ratio(0, 1)]),
        normal: rational_vector([ratio(1, 1), ratio(0, 1), ratio(0, 1)]),
        ref_dir: rational_vector([ratio(0, 1), ratio(1, 1), ratio(0, 1)]),
    }));
    let cone = geom.add_surface(Surface::Cone(Cone {
        apex: rational_point([ratio(13, 1), ratio(apex, 1), ratio(0, 1)]),
        axis: rational_vector([ratio(0, 1), ratio(1, 1), ratio(1, 1)]),
        slope: Real::one(),
        ref_dir: rational_vector([ratio(1, 1), ratio(0, 1), ratio(0, 1)]),
    }));
    let witness = geom.add_point(PointDef::Meet(Meet {
        of: MeetOf::Surfaces(vec![plane, cone]),
        position: Position::Closed(rational_point([
            ratio(13, 1),
            ratio(witness, 1),
            ratio(0, 1),
        ])),
    }));
    geom.add_curve(Curve::Intersection(IntersectionCurve {
        surfaces: [plane, cone],
        witness,
    }))
}

// BM's only root is M at y=5. A cone generator on BM's line keeps M only when
// its ray reaches it, so the line's root is never published for a ray that stops
// short of M.
#[test]
fn a_cone_generator_meets_the_patch_only_on_its_retained_ray() {
    let row = &line_rows()[0];
    let mut built = evaluate(&parsed(&row.recipe), &mut decider()).unwrap();
    let body = built.result;
    let (face, _) = mouth_support(&built.model, body, row);
    let surface = built.model.face(face).surface;
    let stops_short = generator(&mut built.model, 10, 11);
    let answer =
        intersect_curve_surface(built.model.geom_mut(), stops_short, surface, &mut decider());
    let answer = answer.unwrap();
    assert!(
        matches!(answer.value, CurveSurface::Disjoint),
        "{:?}",
        answer.value
    );
    assert!(answer.marks.is_empty());

    let reaches = generator(&mut built.model, 0, 1);
    let answer =
        intersect_curve_surface(built.model.geom_mut(), reaches, surface, &mut decider()).unwrap();
    assert!(answer.marks.is_empty());
    let CurveSurface::Meets(hits) = answer.value else {
        panic!("{:?}", answer.value);
    };
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].contact, Touch::Crossing);
    let m = rational_point([ratio(13, 1), ratio(5, 1), ratio(0, 1)]);
    let gap = built.model.geom().position(hits[0].point).minus(&m);
    assert!(matches!(
        decider().sign(&gap.norm_squared()),
        Outcome::ProvenCoincident(_)
    ));
}

// A root in an end span also carries the irrational transition weight.
#[test]
fn a_crossing_in_an_end_span_stays_on_the_carrier() {
    let row = &line_rows()[0];
    let mut built = evaluate(&parsed(&row.recipe), &mut decider()).unwrap();
    let body = built.result;
    let (face, envelope) = mouth_support(&built.model, body, row);
    let surface = built.model.face(face).surface;
    let (u, v) = (ratio(509, 512), ratio(1, 2));
    let bits = 96;
    let (u, v) = (
        catallax_numbers::Ball::from_rational(&u, bits),
        catallax_numbers::Ball::from_rational(&v, bits),
    );
    let dec = decider();
    let target = envelope.carrier().evaluate(&u, &v, &dec).unwrap();
    let through = rational_point([target[0].mid(), ratio(0, 1), target[2].mid()]);
    let up = rational_vector([ratio(0, 1), ratio(1, 1), ratio(0, 1)]);
    let line = Line::new(through, up, &mut decider()).unwrap();
    let curve = built.model.geom_mut().add_curve(Curve::Line(line));
    let answer =
        intersect_curve_surface(built.model.geom_mut(), curve, surface, &mut decider()).unwrap();
    let CurveSurface::Meets(hits) = answer.value else {
        panic!("{:?}", answer.value);
    };
    assert_eq!(hits.len(), 1);
    let at = built
        .model
        .geom()
        .enclose_point(hits[0].point, &ratio(1, 1_000_000), &mut decider())
        .unwrap();
    let gap = (&at[1] - &target[1]).hi().abs();
    assert!(gap < ratio(1, 100_000), "{gap:?}");
}

// V8 J2 caller policy. The point and line lie a million millimetres above the
// whole carrier, so each answer is a hand-proved miss. A fresh enclosure above
// the ceiling would refuse at the budgeted boundary, and every fresh coordinate
// must cross that boundary: `max_bits` alone also passes without it.
#[test]
fn fresh_query_enclosures_stay_under_the_caller_ceiling() {
    let mut built = evaluate(&parsed(&canonical()), &mut decider()).unwrap();
    let body = built.result;
    let supports = fitted_supports(&built.model, body);
    let root = || Real::from_i64(2).sqrt();
    let far = || Point3::new(root(), Real::from_i64(1_000_000), Real::zero());
    let limited = |bits| {
        Decider::new(BudgetPolicy {
            max_work: APPLICATION_ALLOWANCE,
            max_precision_bits: bits,
        })
    };
    // Three query coordinates and twelve frame coordinates.
    const FRESH: u64 = 15;
    let mut dec = limited(64);
    let (answer, work) = with_operation_budget(APPLICATION_ALLOWANCE, || {
        face_contains(&mut built.model, supports[0].0, &far(), &mut dec)
    });
    assert!(!work.exhausted);
    assert!(matches!(answer, Ok(false)), "{answer:?}");
    assert!(dec.stats().refined >= FRESH, "{:?}", dec.stats());
    assert!(dec.stats().max_bits <= 64, "{:?}", dec.stats());
    let line = Line::new(
        far(),
        Vector3::new(Real::one(), Real::zero(), root()),
        &mut limited(64),
    )
    .unwrap();
    let line = built.model.geom_mut().add_curve(Curve::Line(line));
    for (face, _) in &supports {
        let surface = built.model.face(*face).surface;
        let mut dec = limited(64);
        let answer =
            intersect_curve_surface(built.model.geom_mut(), line, surface, &mut dec).unwrap();
        assert!(matches!(answer.value, CurveSurface::Disjoint));
        assert!(answer.marks.is_empty());
        assert!(dec.stats().refined >= FRESH + 3, "{:?}", dec.stats());
        assert!(dec.stats().max_bits <= 64, "{:?}", dec.stats());
    }
    let answer = face_contains(&mut built.model, supports[0].0, &far(), &mut limited(1));
    assert!(matches!(answer, Ok(false)), "{answer:?}");
    let answer = face_contains(&mut built.model, supports[0].0, &far(), &mut limited(0));
    assert_eq!(answer.unwrap_err().code, "numbers.precision-limit");
}

// A valid face keeping only the local z<=0 half of the +x patch. Support incidence
// must not make the mouth point N, at z=64/65, a member of it.
#[test]
fn a_truncated_envelope_face_keeps_its_trimmed_domain() {
    let row = &line_rows()[0];
    let mut built = evaluate(&parsed(&row.recipe), &mut decider()).unwrap();
    let body = built.result;
    let (face, envelope) = mouth_support(&built.model, body, row);
    let mut dec = decider();
    assert!(envelope_retains_chart(&built.model, face, &mut dec).unwrap());
    let model = &mut built.model;
    let surface = model.face(face).surface;
    let carrier = envelope.carrier();
    // Column of the exact chart section at the symmetry knot s=0.
    let [columns, rows] = carrier.control_count();
    assert_eq!(rows, 3);
    let column = (columns - 1) / 2;
    let section: Vec<_> = carrier.controls()[column * 3..column * 3 + 3].to_vec();
    let weights: Vec<_> = carrier.weights()[column * 3..column * 3 + 3].to_vec();
    let [rail_end, _, mouth_end] = [&section[0], &section[1], &section[2]];
    assert!(envelope.on_boundary(rail_end, &mut dec).unwrap());
    assert!(envelope.on_boundary(mouth_end, &mut dec).unwrap());
    let rail_end = model
        .geom_mut()
        .add_point(PointDef::Explicit(rail_end.clone()));
    let mouth_end = model
        .geom_mut()
        .add_point(PointDef::Explicit(mouth_end.clone()));
    let rail_end = model.add_vertex(Vertex { point: rail_end });
    let mouth_end = model.add_vertex(Vertex { point: mouth_end });
    let [foot, plane] = envelope.definition().corners(Sign::Negative);
    let at_corner = |model: &Model, vertex: VertexId, corner: &Point3| {
        let gap = model
            .geom()
            .position(model.vertex(vertex).point)
            .minus(corner)
            .norm_squared();
        matches!(decider().sign(&gap), Outcome::ProvenCoincident(_))
    };
    // Each kept row: its curve and sense, with the kept end first along the edge.
    let mut kept = BTreeMap::new();
    let rows_of_face = envelope_rows(model, face, &mut dec).unwrap();
    for (edge, row, _) in rows_of_face {
        kept.insert(format!("{row:?}"), model.edge(edge).clone());
    }
    let mut piece = |model: &mut Model, row: &str, from: &Point3, to: Option<VertexId>| {
        let edge = kept.remove(row).unwrap();
        let [start, end] = [edge.start.unwrap(), edge.end.unwrap()];
        let forward = at_corner(model, start, from);
        assert!(forward || at_corner(model, end, from));
        let [start, end] = match (forward, to) {
            (_, None) => [start, end],
            (true, Some(to)) => [start, to],
            (false, Some(to)) => [to, end],
        };
        (
            model.add_edge(Edge {
                curve: edge.curve,
                sense: edge.sense,
                start: Some(start),
                end: Some(end),
            }),
            if forward {
                Sense::Same
            } else {
                Sense::Reversed
            },
        )
    };
    let meridian = piece(model, "Negative", &foot, None);
    let mouth = piece(model, "Mouth", &plane, Some(mouth_end));
    let rail = piece(model, "Rail", &foot, Some(rail_end));
    let knots = [0, 0, 0, 1, 1, 1].map(Rational::from_i64).to_vec();
    let symmetry = BSplineCurve::new(2, knots, section, weights, &mut dec).unwrap();
    let symmetry = model.geom_mut().add_curve(Curve::BSpline(symmetry));
    let symmetry = model.add_edge(Edge {
        curve: symmetry,
        sense: Sense::Same,
        start: Some(rail_end),
        end: Some(mouth_end),
    });
    // Foot to mouth along the meridian, on to the section, down it and back.
    let flip = |sense| match sense {
        Sense::Same => Sense::Reversed,
        Sense::Reversed => Sense::Same,
    };
    let walk = [
        meridian,
        mouth,
        (symmetry, Sense::Reversed),
        (rail.0, flip(rail.1)),
    ];
    let ends = |model: &Model, (edge, sense): (EdgeId, Sense)| {
        let record = model.edge(edge);
        match sense {
            Sense::Same => (record.start.unwrap(), record.end.unwrap()),
            Sense::Reversed => (record.end.unwrap(), record.start.unwrap()),
        }
    };
    for i in 0..4 {
        assert_eq!(ends(model, walk[i]).1, ends(model, walk[(i + 1) % 4]).0);
    }
    let coedges = walk.map(|(edge, sense)| model.add_coedge(Coedge { edge, sense }));
    let boundary = model.add_loop(Loop::Coedges(coedges.to_vec()));
    let half = model.add_face(Face {
        surface,
        sense: model.face(face).sense,
        loops: vec![boundary],
    });
    assert!(!envelope_retains_chart(model, half, &mut dec).unwrap());
    let n = rational_point(row.point([ratio(853, 65), ratio(5, 1), ratio(64, 65)]));
    assert!(face_contains(model, face, &n, &mut decider()).unwrap());
    match face_contains(model, half, &n, &mut decider()) {
        Ok(found) => assert!(!found),
        Err(refusal) => assert_eq!(refusal.code, "shapes.bspline-domain"),
    }
    // A ray's recorded hit at z>0 is on the support, yet not on the half face.
    let shell = model.lump(model.body(body).lumps[0]).shells[0];
    let mut faces = model.shell(shell).faces.clone();
    let at = faces.iter().position(|kept| *kept == face).unwrap();
    faces[at] = half;
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let mut cut = model.body(body).clone();
    cut.lumps = vec![lump];
    let cut = model.add_body(cut);
    let point = rational_point([ratio(3347, 260), ratio(19, 4), ratio(64, 65)]);
    let answer = classify_point(model, cut, &point, &mut decider());
    match answer {
        Outcome::Refused(refusal) => assert_eq!(refusal.code, "shapes.bspline-domain"),
        answer => panic!("{answer:?}"),
    }
}
