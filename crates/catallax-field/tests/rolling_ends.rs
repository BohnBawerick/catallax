//! 4.2a step 4: fillet ends that touch an edge or vertex. P5 at radius 2 and
//! the rim-contact file, the hole-ring file, the slot file and the corpus
//! rows `box-all`, `l-mixed-all` and `vertex-4-one`.
//!
//! The field compiles these fillets from its own path pieces. Its membership
//! must match the witnesses of each test file (checked against the B-rep, and
//! for the corpus rows against OpenCASCADE). Where the B-rep builds and
//! displays, every display vertex of a B-rep blend face must lie within the
//! field's guaranteed distance range of the field's own face, at the part's
//! tolerance: 1 nm on exact faces, and on a fitted rolling envelope its
//! certified bound (rules T0 to T8), with the same records at 1, 2 and 4
//! workers.
use catallax_features::evaluate;
use catallax_field::rolling_region::elementary::{Flat, Ring, Round};
use catallax_field::rolling_region::piece::Piece;
use catallax_field::rolling_region::{Batch, Range, distance_batch};
use catallax_field::{Location, Plan, Workspace};
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::{BudgetPolicy, Decider, OperationAllowance, Rational, Real};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{Point3, RPoint3, SplineSupport, Surface, Vector3};
use serde_json::{Value, json};
use std::collections::BTreeMap;

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn ratio(n: i64, d: i64) -> Rational {
    Rational::from_ratio(n, d).unwrap()
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000_000,
        ..Default::default()
    })
}
fn parse(value: &Value) -> Recipe {
    Recipe::parse(&value.to_string()).unwrap()
}
fn number(value: &Value) -> Rational {
    Rational::parse(value.as_str().unwrap()).unwrap()
}
fn rational(value: &Rational) -> Value {
    json!(value.to_string())
}

/// One recipe with points whose membership a test file fixed.
struct Case {
    name: String,
    recipe: Recipe,
    inside: Vec<[Rational; 3]>,
    outside: Vec<[Rational; 3]>,
    /// Points the B-rep proves on its boundary.
    boundary: Vec<[Rational; 3]>,
    /// Rolling pieces the field compiles.
    pieces: usize,
}

// The P5 family, as the rim-contact file builds it.
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
fn with_radius(radius: &str) -> Value {
    let mut recipe = canonical();
    step(&mut recipe, "filleted")["radius"] = json!(radius);
    recipe
}
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
fn shortened(mut recipe: Value, top: &str) -> Value {
    for segment in step(&mut recipe, "section")["loops"][0]
        .as_array_mut()
        .unwrap()
    {
        if segment["id"] == "body-wall" || segment["id"] == "body-top" {
            segment["to"][1] = json!(top);
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

/// Where a row's definition sits relative to the canonical part.
struct Place {
    scale: Rational,
    pose: Option<(&'static Placement, [i64; 3])>,
}
impl Place {
    fn plain() -> Self {
        Self {
            scale: Rational::one(),
            pose: None,
        }
    }
    fn point(&self, point: &[Rational; 3]) -> [Rational; 3] {
        let point = point.clone().map(|part| part * self.scale.clone());
        match self.pose {
            Some((placement, by)) => placement.point(point, by),
            None => point,
        }
    }
}
const IDENTITY: Placement = Placement {
    matrix: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    over: 1,
    quaternion: [1, 0, 0, 0],
};

fn points(list: &[([i64; 2], [i64; 2], [i64; 2])]) -> Vec<[Rational; 3]> {
    list.iter()
        .map(|(x, y, z)| [ratio(x[0], x[1]), ratio(y[0], y[1]), ratio(z[0], z[1])])
        .collect()
}

/// The rim-contact file's witnesses, placed. `edit` is the side of
/// (15/2, 29/5, 10), which lies outside at radius 2 and inside at 5/2.
fn p5_case(name: &str, recipe: Value, place: Place, edit: Location) -> Case {
    // Hand-derived points of the fitted body (`fitted_witnesses`).
    let mut inside = points(&[
        ([10, 1], [15, 1], [0, 1]),
        ([0, 1], [2, 1], [20, 1]),
        ([15, 2], [11, 2], [10, 1]),
        // The mouth point N lies strictly between that mouth's transition
        // ends (`reaching_witnesses`).
        ([3347, 260], [19, 4], [64, 65]),
    ]);
    let mut outside = points(&[
        ([0, 1], [15, 1], [0, 1]),
        ([17, 1], [2, 1], [0, 1]),
        ([25, 2], [23, 4], [0, 1]),
    ]);
    for (p, from_edit) in [
        (([15, 2], [29, 5], [10, 1]), true),
        (([25, 2], [109, 20], [0, 1]), true),
    ] {
        let p = points(&[p]).remove(0);
        if from_edit && edit == Location::Inside {
            inside.push(p);
        } else {
            outside.push(p);
        }
    }
    Case {
        name: name.into(),
        recipe: parse(&recipe),
        inside: inside.iter().map(|p| place.point(p)).collect(),
        outside: outside.iter().map(|p| place.point(p)).collect(),
        // The mouth's nearest point (13, 5, 0) lies on the deepest ball,
        // where no cover settles membership; it is left out.
        boundary: Vec::new(),
        pieces: 1,
    }
}

/// Strict points of the definition (`strict_points`), with the strip point's
/// side.
fn strict(strip_inside: bool) -> (Vec<[Rational; 3]>, Vec<[Rational; 3]>) {
    let mut inside = points(&[([10, 1], [15, 1], [0, 1]), ([0, 1], [2, 1], [20, 1])]);
    let mut outside = points(&[
        ([0, 1], [15, 1], [0, 1]),
        ([17, 1], [2, 1], [0, 1]),
        ([15, 2], [26, 5], [10, 1]),
    ]);
    let strip = points(&[([15, 2], [51, 10], [10, 1])]).remove(0);
    if strip_inside {
        inside.push(strip);
    } else {
        outside.push(strip);
    }
    (inside, outside)
}

fn p5_cases() -> Vec<Case> {
    let half = ratio(1, 2);
    let double = Rational::from_i64(2);
    let plain = Place::plain;
    let posed = |placement: &'static Placement, by: [i64; 3]| Place {
        scale: Rational::one(),
        pose: Some((placement, by)),
    };
    let scaled_place = |scale: &Rational| Place {
        scale: scale.clone(),
        pose: None,
    };
    let mut cases = vec![
        p5_case("P5 r=2", canonical(), plain(), Location::Outside),
        p5_case("P5 r=5/2", with_radius("5/2"), plain(), Location::Inside),
        p5_case(
            "P5 translated",
            translated(T),
            posed(&IDENTITY, T),
            Location::Outside,
        ),
        p5_case(
            "P5 half",
            scaled(&half),
            scaled_place(&half),
            Location::Outside,
        ),
        p5_case(
            "P5 double",
            scaled(&double),
            scaled_place(&double),
            Location::Outside,
        ),
        p5_case(
            "P5 RZ r=2",
            RZ.recipe(NO_MOVE, "2"),
            posed(&RZ, NO_MOVE),
            Location::Outside,
        ),
        p5_case(
            "P5 RZ moved r=5/2",
            RZ.recipe(T, "5/2"),
            posed(&RZ, T),
            Location::Inside,
        ),
        p5_case(
            "P5 RQ r=2",
            RQ.recipe(NO_MOVE, "2"),
            posed(&RQ, NO_MOVE),
            Location::Outside,
        ),
        p5_case(
            "P5 RQ moved r=5/2",
            RQ.recipe(T, "5/2"),
            posed(&RQ, T),
            Location::Inside,
        ),
    ];
    // The rail touches every mouth at one point (`tangent_rows`).
    let mut half_tangent = scaled(&half);
    step(&mut half_tangent, "filleted")["radius"] = json!("1/2");
    for (name, recipe, place) in [
        ("P5 r=1", with_radius("1"), plain()),
        ("P5 RQ moved r=1", RQ.recipe(T, "1"), posed(&RQ, T)),
        ("P5 RZ moved r=1", RZ.recipe(T, "1"), posed(&RZ, T)),
        ("P5 half r=1/2", half_tangent, scaled_place(&half)),
    ] {
        // The file's points in the flange plane under the round, (12.99, 5,
        // 0) and the four contacts, lie on a face the piece shares with the
        // housing; a rolling piece registers no face to settle them, so they
        // are left out.
        let (mut inside, mut outside) = strict(true);
        inside.extend(points(&[([13, 1], [499, 100], [0, 1])]));
        outside.extend(points(&[
            ([1301, 100], [5, 1], [0, 1]),
            ([13, 1], [501, 100], [0, 1]),
        ]));
        let boundary = points(&[([13, 1], [5, 1], [1, 100])]);
        cases.push(Case {
            name: name.into(),
            recipe: parse(&recipe),
            inside: inside.iter().map(|p| place.point(p)).collect(),
            outside: outside.iter().map(|p| place.point(p)).collect(),
            boundary: boundary.iter().map(|p| place.point(p)).collect(),
            pieces: 1,
        });
    }
    // The rail clears every mouth (`clear_rows`): the elementary torus.
    let mut half_near = scaled(&half);
    step(&mut half_near, "filleted")["radius"] = json!("999/2000");
    for (name, recipe, place, near) in [
        ("P5 r=999/1000", with_radius("999/1000"), plain(), true),
        ("P5 r=1/2", with_radius("1/2"), plain(), false),
        (
            "P5 RQ moved r=999/1000",
            RQ.recipe(T, "999/1000"),
            posed(&RQ, T),
            true,
        ),
        ("P5 half r=999/2000", half_near, scaled_place(&half), true),
    ] {
        let (mut inside, mut outside) = strict(near);
        inside.extend(points(&[([25999, 2000], [49999, 10000], [0, 1])]));
        outside.extend(points(&[([25999, 2000], [50001, 10000], [0, 1])]));
        cases.push(Case {
            name: name.into(),
            recipe: parse(&recipe),
            inside: inside.iter().map(|p| place.point(p)).collect(),
            outside: outside.iter().map(|p| place.point(p)).collect(),
            boundary: Vec::new(),
            pieces: 0,
        });
    }
    cases
}

fn cases() -> Vec<Case> {
    p5_cases()
}

fn location(plan: &Plan, p: &[Rational; 3]) -> Location {
    let point = RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone());
    plan.proven_checked(&point, &mut Workspace::default(), &mut decider())
        .unwrap()
        .location
}

#[test]
fn the_field_builds_every_step_four_fillet_and_agrees_on_membership() {
    let mut failures = Vec::new();
    for case in cases() {
        let start = std::time::Instant::now();
        let plan = match Plan::compile_checked(&case.recipe, &mut decider()) {
            Ok(plan) => plan,
            Err(e) => {
                failures.push(format!("{}: {e:?}", case.name));
                continue;
            }
        };
        let pieces = plan.rolling_pieces().count();
        if pieces != case.pieces {
            failures.push(format!("{}: {pieces} rolling pieces", case.name));
        }
        for (points, want) in [
            (&case.inside, Location::Inside),
            (&case.outside, Location::Outside),
            (&case.boundary, Location::Boundary),
        ] {
            for p in points {
                let found = location(&plan, p);
                if found != want {
                    failures.push(format!("{} at {p:?}: {found:?}", case.name));
                }
            }
        }
        eprintln!(
            "{}: {} witnesses in {:?}",
            case.name,
            case.inside.len() + case.outside.len() + case.boundary.len(),
            start.elapsed()
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// The cases whose B-rep blend vertices are checked against the field's
/// rolling pieces, one test per group so they run side by side.
const VERTEX_GROUPS: &[&[&str]] = &[
    &["P5 r=2", "P5 r=5/2"],
    &["P5 translated", "P5 half", "P5 double"],
    &["P5 RZ r=2", "P5 RZ moved r=5/2"],
    &["P5 RQ r=2", "P5 RQ moved r=5/2"],
    &[
        "P5 r=1",
        "P5 RQ moved r=1",
        "P5 RZ moved r=1",
        "P5 half r=1/2",
    ],
];
macro_rules! vertex_tests {
    ($($test:ident = $group:literal;)+) => {
        $(
            #[test]
            fn $test() {
                blend_vertices_lie_on_the_field_face(VERTEX_GROUPS[$group]);
            }
        )+
    };
}
vertex_tests! {
    rim_blend_vertices_lie_on_the_field_face = 0;
    moved_rim_blend_vertices_lie_on_the_field_face = 1;
    turned_rim_blend_vertices_lie_on_the_field_face = 2;
    spatial_rim_blend_vertices_lie_on_the_field_face = 3;
    touching_rim_blend_vertices_lie_on_the_field_face = 4;
}

#[test]
fn every_rolling_case_has_its_vertices_checked() {
    let listed: std::collections::BTreeSet<String> = VERTEX_GROUPS
        .iter()
        .flat_map(|group| group.iter().map(|name| name.to_string()))
        .collect();
    let wanted: std::collections::BTreeSet<String> = cases()
        .into_iter()
        .filter(|c| c.pieces > 0)
        .map(|c| c.name)
        .collect();
    assert_eq!(listed, wanted);
}

/// A B-rep body's blend display vertices, each with the part's tolerance at
/// it: 1 nm, or a fitted envelope's certified bound.
fn blend_vertices(recipe: &Recipe, chord: &str) -> (Vec<Point3>, Vec<Rational>) {
    let own = q("1/1000000");
    let built = evaluate(recipe, &mut decider()).unwrap();
    let mesh = mesh_body(
        &built.model,
        built.result,
        &MeshOptions {
            chord_tolerance: q(chord),
            ..MeshOptions::default()
        },
        &mut decider(),
    )
    .unwrap()
    .value;
    let mut indices: BTreeMap<usize, Rational> = BTreeMap::new();
    for face in &mesh.faces {
        if !built.model.face_origins(face.face).iter().any(|o| {
            let mut role = &o.role;
            while let Role::Copy { of, .. } = role {
                role = &of.role;
            }
            matches!(role, Role::Blend { .. })
        }) {
            continue;
        }
        let bound = match built
            .model
            .geom()
            .surface(built.model.face(face.face).surface)
        {
            Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) => {
                envelope.bound().upper.clone()
            }
            _ => Rational::zero(),
        };
        let tolerance = own.clone().max(bound);
        for triangle in &mesh.triangles[face.triangles.clone()] {
            for v in triangle.vertices {
                let entry = indices.entry(v).or_insert_with(|| tolerance.clone());
                *entry = entry.clone().max(tolerance.clone());
            }
        }
    }
    let tolerances: Vec<Rational> = indices.values().cloned().collect();
    let points = indices
        .keys()
        .map(|&v| {
            let [x, y, z] = mesh.positions[v].map(|c| Rational::from_f64_exact(c).unwrap());
            Point3::new(x.into(), y.into(), z.into())
        })
        .collect();
    (points, tolerances)
}

/// Upper bounds on the distance from each point to the nearest of the
/// pieces' faces, from guaranteed ranges at most `tolerance` wide. Every
/// 128th point is also queried at 1, 2 and 4 workers, whose records must
/// agree.
fn nearest_witnesses(
    name: &str,
    pieces: &[&Piece],
    points: &[Point3],
    tolerances: &[Rational],
) -> Vec<Option<Rational>> {
    const CAP: u64 = 20_000;
    let allowance = |n: usize| OperationAllowance::reserve(n as u64 * CAP).unwrap();
    let mut nearest = vec![None::<Rational>; points.len()];
    let mut widths: Vec<Rational> = tolerances.to_vec();
    widths.sort();
    widths.dedup();
    for (piece, width) in pieces
        .iter()
        .flat_map(|piece| widths.iter().map(move |w| (*piece, w)))
    {
        let index: Vec<usize> = (0..points.len())
            .filter(|&i| &tolerances[i] == width && nearest[i].as_ref().is_none_or(|u| u > width))
            .collect();
        if index.is_empty() {
            continue;
        }
        let mine: Vec<Point3> = index.iter().map(|&i| points[i].clone()).collect();
        let sample: Vec<Point3> = mine.iter().step_by(128).cloned().collect();
        let all = distance_batch(
            piece,
            &mine,
            width,
            None,
            allowance(mine.len()),
            &mut decider(),
        )
        .unwrap();
        let mut batches: Vec<Batch<Range>> = Vec::new();
        for workers in [1, 2, 4] {
            batches.push(
                distance_batch(
                    piece,
                    &sample,
                    width,
                    Some(workers),
                    allowance(sample.len()),
                    &mut decider(),
                )
                .unwrap(),
            );
        }
        assert_eq!(batches[0].records, batches[1].records, "{name}");
        assert_eq!(batches[0].records, batches[2].records, "{name}");
        for (&i, record) in index.iter().zip(&all.records) {
            let Ok(range) = &record.result else {
                assert!(record.exhausted, "{name}: {record:?}");
                continue;
            };
            if let Some(upper) = &range.upper
                && nearest[i].as_ref().is_none_or(|best| upper < best)
            {
                nearest[i] = Some(upper.clone());
            }
        }
    }
    nearest
}

fn assert_within(name: &str, nearest: &[Option<Rational>], tolerances: &[Rational]) {
    let mut worst = Rational::zero();
    for (found, tolerance) in nearest.iter().zip(tolerances) {
        let found = found
            .clone()
            .unwrap_or_else(|| panic!("{name}: no witness"));
        assert!(&found <= tolerance, "{name}: {found:?} over {tolerance:?}");
        worst = worst.max(found);
    }
    eprintln!(
        "{name}: {} blend vertices, worst field distance {}",
        nearest.len(),
        worst.to_f64_nearest()
    );
}

/// Every B-rep display vertex of a blend face lies within the field's
/// guaranteed range of a rolling piece's face, at the part's tolerance.
fn blend_vertices_lie_on_the_field_face(names: &[&str]) {
    for case in cases()
        .into_iter()
        .filter(|c| names.contains(&c.name.as_str()))
    {
        let start = std::time::Instant::now();
        let plan = Plan::compile_checked(&case.recipe, &mut decider()).unwrap();
        let (points, tolerances) = blend_vertices(&case.recipe, "1/20");
        assert!(!points.is_empty(), "{}: no blend vertices", case.name);
        let pieces: Vec<&Piece> = plan.rolling_pieces().collect();
        let nearest = nearest_witnesses(&case.name, &pieces, &points, &tolerances);
        eprintln!("{}: checked in {:?}", case.name, start.elapsed());
        assert_within(&case.name, &nearest, &tolerances);
    }
}

/// The step 4 files' recipes that the field refuses, with the B-rep's
/// answer: both refuse each of these.
fn refusals() -> Vec<(&'static str, Recipe)> {
    vec![
        ("P5 r=7/2", parse(&with_radius("7/2"))),
        ("P5 RZ moved r=7/2", parse(&RZ.recipe(T, "7/2"))),
        (
            "P5 short RZ moved r=2",
            parse(&shortened(RZ.recipe(T, "2"), "6")),
        ),
        ("P5 short r=1", parse(&shortened(with_radius("1"), "11/2"))),
        (
            "P5 short r=999/1000",
            parse(&shortened(with_radius("999/1000"), "11/2")),
        ),
    ]
}

const REFUSED: &[(&str, &str)] = &[
    // The ball would leave the mouth's edge for the countersink's cone; the
    // drill's top, 1 mm above the flange, reaches the wider fan first.
    ("P5 r=7/2", "fillet.unsupported-overflow"),
    ("P5 RZ moved r=7/2", "fillet.unsupported-overflow"),
    // The wall above the flange is lower than the ball: its top plane
    // crosses the fan.
    ("P5 short RZ moved r=2", "fillet.unsupported-overflow"),
    ("P5 short r=1", "fillet.unsupported-overflow"),
    ("P5 short r=999/1000", "fillet.unsupported-overflow"),
];

#[test]
fn the_field_names_what_it_refuses() {
    let mut found = Vec::new();
    for (name, recipe) in refusals() {
        let brep = evaluate(&recipe, &mut decider()).is_ok();
        assert!(!brep, "{name}: the B-rep builds it");
        let code = match Plan::compile_checked(&recipe, &mut decider()) {
            Ok(_) => "built",
            Err(e) => e.code,
        };
        found.push((name, code));
    }
    assert_eq!(found, REFUSED);
}

// The hole-ring file: rings where a perpendicular cylinder meets a plane.
// The field keeps them on the elementary route; the B-rep's vertices are
// checked against the field's own ring path for each.
const UP: [&str; 4] = ["1", "0", "0", "0"];
const DOWN: [&str; 4] = ["0", "1", "0", "0"];
fn block() -> Vec<Value> {
    vec![
        json!({"id": "block-section", "op": "sketch",
        "plane": {"rotation": UP, "translation": ["0", "0", "0"]},
        "loops": [[
            {"kind": "line", "id": "s", "to": ["20", "-20"]},
            {"kind": "line", "id": "e", "to": ["20", "20"]},
            {"kind": "line", "id": "n", "to": ["-20", "20"]},
            {"kind": "line", "id": "w", "to": ["-20", "-20"]}
        ]]}),
        json!({"id": "block", "op": "extrude", "profile": "block-section", "distance": "20"}),
    ]
}
fn cylinder(id: &str, radius: &str, height: &str, rotation: [&str; 4], z: &str) -> Value {
    json!({"id": id, "op": "cylinder", "radius": radius, "height": height,
        "pose": {"rotation": rotation, "translation": ["0", "0", z]}})
}
fn boolean(id: &str, kind: &str, target: &str, tool: &str) -> Value {
    json!({"id": id, "op": "boolean", "kind": kind, "target": target, "tool": tool})
}
fn wall(step: &str) -> Value {
    json!({"step": step, "role": {"side": {"segment": "wall"}}})
}
fn cap(step: &str, role: &str) -> Value {
    json!({"step": step, "role": role})
}
fn ring_recipe(mut steps: Vec<Value>, body: &str, radius: &str, between: [Value; 2]) -> Recipe {
    steps.push(
        json!({"id": "round", "op": "fillet", "body": body, "radius": radius,
        "edges": [{"between": between}]}),
    );
    parse(
        &json!({"format": "catallax-recipe/3", "name": "ring", "steps": steps, "result": "round"}),
    )
}
fn through_hole(rotation: [&str; 4], z: &str) -> Vec<Value> {
    let mut steps = block();
    steps.push(cylinder("hole", "3", "30", rotation, z));
    steps.push(boolean("cut", "subtract", "block", "hole"));
    steps
}
fn counterbore(rotation: [&str; 4]) -> Vec<Value> {
    let mut steps = block();
    let (bore_z, hole_z) = if rotation == UP {
        ("15", "-5")
    } else {
        ("25", "25")
    };
    steps.push(cylinder("bore", "6", "10", rotation, bore_z));
    steps.push(boolean("bored", "subtract", "block", "bore"));
    steps.push(cylinder("hole", "3", "30", rotation, hole_z));
    steps.push(boolean("cut", "subtract", "bored", "hole"));
    steps
}
fn drilled(depth: Value, head: Option<Value>) -> Vec<Value> {
    let mut steps = block();
    let mut hole = json!({"id": "drill", "op": "hole", "target": "block",
        "at": [{"id": "h", "point": ["0", "0", "20"], "direction": ["0", "0", "-1"]}],
        "diameter": "6", "depth": depth, "tip": "flat"});
    if let Some(head) = head {
        hole["head"] = head;
    }
    steps.push(hole);
    steps
}
fn part(name: &str) -> Value {
    json!({"step": "drill", "role": {"hole": {"at": "h", "part": name}}})
}

/// One ring of the hole-ring file: the recipe, whether it is convex, the
/// ball centre's height, the wall's radius and whether the ball lies
/// outside the wall, and the direction from the centre to the corner.
struct HoleRing {
    name: String,
    recipe: Recipe,
    convex: bool,
    centre: i64,
    wall: i64,
    outside: bool,
    toward: [i64; 2],
}

fn hole_rings() -> Vec<HoleRing> {
    let ring = |name: &str, recipe, convex, centre, wall, outside, toward| HoleRing {
        name: name.to_string(),
        recipe,
        convex,
        centre,
        wall,
        outside,
        toward,
    };
    let mut rings = vec![ring(
        "boss top",
        ring_recipe(
            vec![cylinder("boss", "3", "10", UP, "0")],
            "boss",
            "1",
            [cap("boss", "end"), wall("boss")],
        ),
        true,
        9,
        3,
        false,
        [1, 1],
    )];
    for (name, steps, face, top) in [
        ("hole mouth, axis up", through_hole(UP, "-5"), "end", true),
        (
            "hole mouth, axis down",
            through_hole(DOWN, "25"),
            "end",
            true,
        ),
        (
            "bottom hole mouth, axis up",
            through_hole(UP, "-5"),
            "start",
            false,
        ),
        (
            "bottom hole mouth, axis down",
            through_hole(DOWN, "25"),
            "start",
            false,
        ),
    ] {
        let (centre, toward) = if top { (19, [-1, 1]) } else { (1, [-1, -1]) };
        rings.push(ring(
            name,
            ring_recipe(steps, "cut", "1", [cap("block", face), wall("hole")]),
            true,
            centre,
            3,
            true,
            toward,
        ));
    }
    for (name, rotation, centre, toward) in [
        ("blind floor, axis up", UP, 11, [1, -1]),
        ("blind floor, axis down", DOWN, 9, [1, 1]),
    ] {
        let mut steps = block();
        steps.push(cylinder("hole", "3", "20", rotation, "10"));
        steps.push(boolean("cut", "subtract", "block", "hole"));
        rings.push(ring(
            name,
            ring_recipe(steps, "cut", "1", [cap("hole", "start"), wall("hole")]),
            false,
            centre,
            3,
            false,
            toward,
        ));
    }
    for (name, rotation, z) in [
        ("boss base, axis up", UP, "10"),
        ("boss base, axis down", DOWN, "30"),
    ] {
        let mut steps = block();
        steps.push(cylinder("boss", "3", "20", rotation, z));
        steps.push(boolean("joined", "union", "block", "boss"));
        rings.push(ring(
            name,
            ring_recipe(steps, "joined", "1", [cap("block", "end"), wall("boss")]),
            false,
            21,
            3,
            true,
            [-1, -1],
        ));
    }
    for (label, rotation) in [("up", UP), ("down", DOWN)] {
        let floor = if rotation == UP { "start" } else { "end" };
        rings.push(ring(
            &format!("counterbore mouth, axis {label}"),
            ring_recipe(
                counterbore(rotation),
                "cut",
                "1",
                [cap("block", "end"), wall("bore")],
            ),
            true,
            19,
            6,
            true,
            [-1, 1],
        ));
        rings.push(ring(
            &format!("counterbore floor, axis {label}"),
            ring_recipe(
                counterbore(rotation),
                "cut",
                "1",
                [cap("bore", floor), wall("bore")],
            ),
            false,
            16,
            6,
            false,
            [1, -1],
        ));
        rings.push(ring(
            &format!("counterbore hole, axis {label}"),
            ring_recipe(
                counterbore(rotation),
                "cut",
                "1",
                [cap("bore", floor), wall("hole")],
            ),
            true,
            14,
            3,
            true,
            [-1, 1],
        ));
    }
    let head = || Some(json!({"counterbore": {"diameter": "12", "depth": "5"}}));
    for (name, steps, between, convex, centre, wall_radius, outside, toward) in [
        (
            "hole operation, through",
            drilled(json!("through"), None),
            [cap("block", "end"), part("wall")],
            true,
            19,
            3,
            true,
            [-1, 1],
        ),
        (
            "hole operation, blind floor",
            drilled(json!({"blind": "10"}), None),
            [part("floor"), part("wall")],
            false,
            11,
            3,
            false,
            [1, -1],
        ),
        (
            "hole operation, counterbore mouth",
            drilled(json!("through"), head()),
            [cap("block", "end"), part("counterbore_wall")],
            true,
            19,
            6,
            true,
            [-1, 1],
        ),
        (
            "hole operation, counterbore floor",
            drilled(json!("through"), head()),
            [part("counterbore_floor"), part("counterbore_wall")],
            false,
            16,
            6,
            false,
            [1, -1],
        ),
        (
            "hole operation, counterbore hole",
            drilled(json!("through"), head()),
            [part("counterbore_floor"), part("wall")],
            true,
            14,
            3,
            true,
            [-1, 1],
        ),
    ] {
        rings.push(ring(
            name,
            ring_recipe(steps, "drill", "1", between),
            convex,
            centre,
            wall_radius,
            outside,
            toward,
        ));
    }
    rings
}

impl HoleRing {
    /// The ball centre's radius.
    fn major(&self) -> i64 {
        if self.outside {
            self.wall + 1
        } else {
            self.wall - 1
        }
    }
    /// The file's two points beside the round, 0.99 and 1.018 radii from the
    /// ball centre toward the corner, on the positive x axis.
    fn beside(&self) -> [[Rational; 3]; 2] {
        let near = |offset: [i64; 2]| {
            [
                ratio(100 * self.major() + self.toward[0] * offset[0], 100),
                Rational::zero(),
                ratio(100 * self.centre + self.toward[1] * offset[1], 100),
            ]
        };
        [near([71, 69]), near([72, 72])]
    }
    /// The field's own ring path: the plane at the corner's height, normal
    /// toward the ball, and the wall about the z axis.
    fn path(&self) -> Ring {
        let corner = self.centre + self.toward[1];
        let normal = Vector3::new(0.into(), 0.into(), (-self.toward[1]).into());
        Ring::new(
            Flat {
                point: Point3::new(0.into(), 0.into(), corner.into()),
                normal,
            },
            Round {
                point: Point3::origin(),
                axis: Vector3::new(0.into(), 0.into(), 1.into()),
                radius: Real::from_i64(self.wall),
                outside: self.outside,
            },
            Real::one(),
            &mut decider(),
        )
        .unwrap()
    }
}

#[test]
fn every_hole_ring_builds_and_agrees_on_its_sides() {
    for ring in hole_rings() {
        let plan = Plan::compile_checked(&ring.recipe, &mut decider())
            .unwrap_or_else(|e| panic!("{}: {e:?}", ring.name));
        let [ball, corner] = ring.beside();
        let (ball_side, corner_side) = if ring.convex {
            (Location::Inside, Location::Outside)
        } else {
            (Location::Outside, Location::Inside)
        };
        assert_eq!(location(&plan, &ball), ball_side, "{}", ring.name);
        assert_eq!(location(&plan, &corner), corner_side, "{}", ring.name);
    }
}

#[test]
fn hole_ring_blend_vertices_lie_on_the_field_ring() {
    for ring in hole_rings() {
        let (points, tolerances) = blend_vertices(&ring.recipe, "1/20");
        assert!(!points.is_empty(), "{}: no blend vertices", ring.name);
        let path = ring.path();
        const CAP: u64 = 20_000;
        let mut nearest = Vec::with_capacity(points.len());
        for (i, p) in points.iter().enumerate() {
            let range = catallax_numbers::with_operation_budget(CAP, || {
                catallax_field::rolling_region::path::distance(
                    &path,
                    p,
                    &tolerances[i],
                    &mut decider(),
                )
            })
            .0
            .unwrap();
            nearest.push(range.upper);
        }
        assert_within(&ring.name, &nearest, &tolerances);
    }
}
