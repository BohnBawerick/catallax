//! Fillets rolling over earlier fillets. A later fillet's tangent chain
//! continues over an earlier cylinder strip as a torus strip, or as a sphere
//! strip where the two radii match: on a boss top over rounded uprights, a
//! pocket floor over rounded corners, a box corner over a rounded edge and an
//! L's end over its inner round. Each is checked against an exact volume,
//! classified points, a closed display and STEP, with names kept across
//! radius edits; the cases not built yet refuse with typed codes.

use std::collections::{BTreeMap, BTreeSet};

use catallax_boolean::Containment;
use catallax_features::{BuildError, Built, evaluate};
use catallax_mesh::MeshOptions;
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational};
use catallax_recipe::{Origin, Recipe, Role};
use catallax_shapes::RPoint3;
use catallax_topology::FaceId;
use serde_json::{Value, json};

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 20_000_000,
        ..Default::default()
    })
}
/// Exact measurements refine further than a build needs: a body with eight
/// trimmed torus strips takes about two billion units to measure its volume.
fn measuring() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 2_000_000_000,
        ..Default::default()
    })
}
fn r(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn pose(at: [&str; 3]) -> Value {
    json!({"rotation":["1","0","0","0"],"translation":at})
}
fn side(step: &str, segment: &str) -> Value {
    json!({"step":step,"role":{"side":{"segment":segment}}})
}
fn cap(step: &str, role: &str) -> Value {
    json!({"step":step,"role":role})
}
/// Every face a step made, such as an earlier fillet's strips.
fn made(step: &str) -> Value {
    json!({"step":step})
}
fn edge(a: Value, b: Value) -> Value {
    json!({"between":[a,b]})
}
fn fillet(id: &str, body: &str, radius: &str, edges: Vec<Value>) -> Value {
    json!({"id":id,"op":"fillet","body":body,"radius":radius,"edges":edges})
}
fn recipe(format: u8, steps: Vec<Value>) -> Recipe {
    let result = steps.last().unwrap()["id"].clone();
    Recipe::parse(
        &json!({"format":format!("catallax-recipe/{format}"),"name":"fillet over fillet",
            "steps":steps,"result":result})
        .to_string(),
    )
    .unwrap()
}
fn build(recipe: &Recipe) -> Built {
    let built = evaluate(recipe, &mut decider()).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(built.model.check(built.result), vec![]);
    built
}
fn kinds(built: &Built) -> BTreeMap<String, usize> {
    let mut kinds = BTreeMap::new();
    for face in built.model.faces_of_body(built.result) {
        let surface = built.model.geom().surface(built.model.face(face).surface);
        *kinds.entry(format!("{:?}", surface.kind())).or_insert(0) += 1;
    }
    kinds
}
fn expect_kinds(built: &Built, want: &[(&str, usize)]) {
    let want: BTreeMap<_, _> = want.iter().map(|(k, n)| (k.to_string(), *n)).collect();
    assert_eq!(kinds(built), want);
}
/// The blend origins of one step, which a radius edit must keep.
fn names(built: &Built, step: &str) -> BTreeSet<Origin> {
    built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|face| built.model.face_origins(face).to_vec())
        .filter(|o| o.step.as_str() == step && matches!(o.role, Role::Blend { .. }))
        .collect()
}
fn faces_of(built: &Built, step: &str) -> Vec<FaceId> {
    built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|f| {
            built
                .model
                .face_origins(*f)
                .iter()
                .any(|o| o.step.as_str() == step)
        })
        .collect()
}
fn pi() -> (Rational, Rational) {
    let lo = r("3.1415926535897932384626433832795028841971693993751");
    let hi = &lo + &r("1e-49");
    (lo, hi)
}
/// The measured volume overlaps `a + b pi + c pi^2`, within `width`.
fn exact_volume(built: &Built, [a, b, c]: [Rational; 3], width: &str) {
    let volume = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut measuring(),
    )
    .unwrap()
    .value;
    let (lo, hi) = pi();
    let at = |p: &Rational| &(&a + &(&b * p)) + &(&(&c * p) * p);
    let (x, y) = (at(&lo), at(&hi));
    let (below, above) = if x <= y { (x, y) } else { (y, x) };
    assert!(
        volume.lo <= above && volume.hi >= below,
        "[{}, {}] misses [{below}, {above}]",
        volume.lo.to_f64_nearest(),
        volume.hi.to_f64_nearest()
    );
    assert!(&volume.hi - &volume.lo <= r(width), "{volume:?}");
}
fn point(p: [&str; 3]) -> RPoint3 {
    let [x, y, z] = p.map(r);
    RPoint3::new(x, y, z)
}
fn classify(built: &mut Built, points: &[([&str; 3], Containment)]) {
    for (p, want) in points {
        let actual = catallax_boolean::classify_point(
            &mut built.model,
            built.result,
            &point(*p),
            &mut decider(),
        );
        assert!(
            matches!(actual, Outcome::Proven(actual) if actual == *want),
            "B-rep at {p:?}: {actual:?}"
        );
    }
}
/// The kernel's default display and SolidDraft's 15 degree normal limit:
/// every face meshed, every triangle edge used twice in opposite directions.
fn closed_display(built: &Built) {
    let options = MeshOptions {
        angular_deviation: Some(r("0.2617993877991494")),
        ..MeshOptions::default()
    };
    for options in [MeshOptions::default(), options] {
        let mesh = catallax_mesh::mesh_body(&built.model, built.result, &options, &mut decider())
            .unwrap()
            .value;
        let meshed: BTreeSet<_> = mesh.triangles.iter().map(|t| t.face).collect();
        assert_eq!(
            meshed,
            built
                .model
                .faces_of_body(built.result)
                .into_iter()
                .collect()
        );
        let mut incidence = BTreeMap::<_, (usize, i32)>::new();
        for triangle in &mesh.triangles {
            for i in 0..3 {
                let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                let uses = incidence.entry((a.min(b), a.max(b))).or_default();
                uses.0 += 1;
                uses.1 += if a < b { 1 } else { -1 };
            }
        }
        assert!(incidence.values().all(|uses| *uses == (2, 0)));
        for face in &mesh.faces {
            assert!(face.max_chord_error <= options.chord_tolerance);
        }
    }
}
/// The body writes to STEP without a mark, every face at least once.
fn step(built: &Built, surface: &str) {
    let mut bytes = Vec::new();
    let answer = catallax_step::write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    assert!(answer.marks.is_empty(), "{:?}", answer.marks);
    let text = String::from_utf8(bytes).unwrap();
    // A periodic face may be written as several charts.
    assert!(
        text.matches("ADVANCED_FACE(").count() >= built.model.faces_of_body(built.result).len()
    );
    assert!(text.contains(surface), "{surface}");
}
/// The section `r^2 (1 - pi/4)` between two perpendicular supports and the
/// ball, as `[1, pi]` coefficients.
fn section(rho: &Rational) -> [Rational; 2] {
    let r2 = rho * rho;
    [r2.clone(), -&(&r2 / &r("4"))]
}
/// A quarter turn of that section about an axis at `axis` from the support
/// it leaves, its centroid `rho (5/6 - pi/4) / (1 - pi/4)` from both faces,
/// inside (`-1`) or outside (`1`) that support: `pi/2 (axis A -+ rho^3
/// (5/6 - pi/4))`, as `[pi, pi^2]` coefficients.
fn quarter_turn(rho: &Rational, axis: &Rational, outside: i64) -> [Rational; 2] {
    let [a0, a1] = section(rho);
    let r3 = &(rho * rho) * rho;
    let k = Rational::from_i64(outside);
    let half = r("1/2");
    [
        &half * &(&(axis * &a0) + &(&k * &(&r3 * &r("5/6")))),
        &half * &(&(axis * &a1) - &(&k * &(&r3 * &r("1/4")))),
    ]
}

fn block(id: &str, size: [&str; 3], at: [&str; 3]) -> Value {
    json!({"id":id,"op":"block","size":size,"pose":pose(at)})
}
const SIDES: [&str; 4] = ["x_min", "y_min", "x_max", "y_max"];
fn uprights(step: &str) -> Vec<Value> {
    (0..4)
        .map(|i| edge(side(step, SIDES[i]), side(step, SIDES[(i + 1) % 4])))
        .collect()
}

// A 20 x 20 x 10 block with its four uprights rounded by `upright`, then its
// top edges by `top`: one selected edge rolls over all four corner rounds.
fn boss(format: u8, upright: &str, top: &str, every: bool) -> Recipe {
    let edges = if every {
        let mut edges: Vec<_> = SIDES
            .iter()
            .map(|s| edge(cap("block", "end"), side("block", s)))
            .collect();
        edges.push(edge(cap("block", "end"), made("upright")));
        edges
    } else {
        vec![edge(cap("block", "end"), side("block", "x_max"))]
    };
    recipe(
        format,
        vec![
            block("block", ["20", "20", "10"], ["0", "0", "0"]),
            fillet("upright", "block", upright, uprights("block")),
            fillet("top", "upright", top, edges),
        ],
    )
}

#[test]
fn a_boss_top_rolls_over_its_rounded_uprights() {
    let big = r("3");
    for top in ["1", "2", "3"] {
        let rho = r(top);
        let built = build(&boss(4, "3", top, false));
        // A ball past half the round's radius sweeps the outer sheet of a
        // spindle torus, and one as wide as the round a sphere.
        let corner = match top {
            "1" => "Torus",
            "2" => "TorusSheet",
            _ => "Sphere",
        };
        expect_kinds(&built, &[("Plane", 6), ("Cylinder", 8), (corner, 4)]);
        // 4000 less the four upright sections over the height, less the top
        // section along runs 20 - 2R, less four quarter turns inside R.
        let runs = &r("80") - &(&r("8") * &big);
        let [s0, s1] = section(&rho);
        let [u0, u1] = section(&big);
        let [t1, t2] = quarter_turn(&rho, &big, -1);
        let four = r("4");
        exact_volume(
            &built,
            [
                &(&r("4000") - &(&(&four * &u0) * &r("10"))) - &(&s0 * &runs),
                &(&(-&(&(&four * &u1) * &r("10"))) - &(&s1 * &runs)) - &(&four * &t1),
                -&(&four * &t2),
            ],
            "1e-6",
        );
        closed_display(&built);
        step(
            &built,
            if top == "3" {
                "SPHERICAL_SURFACE("
            } else {
                "TOROIDAL_SURFACE("
            },
        );
        // The continued chain keeps the one selected edge's name; listing
        // every edge of the loop builds the same body with eight names.
        let [name] = &names(&built, "top").into_iter().collect::<Vec<_>>()[..] else {
            panic!("one continued name");
        };
        assert!(matches!(name.role, Role::Blend { .. }));
        assert_eq!(faces_of(&built, "top").len(), 8);
        let explicit = build(&boss(3, "3", top, true));
        assert_eq!(kinds(&explicit), kinds(&built));
        assert_eq!(names(&explicit, "top").len(), 8);
        assert!(names(&explicit, "top").contains(name));
    }
    // Radius edits of either fillet keep every name, through the sphere.
    let reference = names(&build(&boss(4, "3", "1", false)), "top");
    for (upright, top) in [("3", "2"), ("3", "3"), ("2", "2"), ("5/2", "1")] {
        assert_eq!(
            names(&build(&boss(4, upright, top, false)), "top"),
            reference
        );
        assert_eq!(
            names(&build(&boss(3, upright, top, true)), "top"),
            names(&build(&boss(3, "3", "1", true)), "top")
        );
    }
    use Containment::{Inside, Outside};
    let mut built = build(&boss(4, "3", "1", false));
    classify(
        &mut built,
        &[
            // Cut by the straight strip over the right face's top edge.
            (["199/10", "10", "99/10"], Outside),
            (["39/2", "10", "19/2"], Inside),
            // Cut by the torus over the corner round about (17, 17): 59/20
            // from its axis and 1/20 below the top, outside the tube.
            (["1877/100", "1936/100", "199/20"], Outside),
            (["37/2", "19", "19/2"], Inside),
            // The corner round itself, and uncut material.
            (["199/10", "199/10", "5"], Outside),
            (["10", "10", "5"], Inside),
        ],
    );
    let mut built = build(&boss(4, "3", "3", false));
    classify(
        &mut built,
        &[
            // The sphere about (17, 17, 7): 1/2 short of its radius, and
            // just outside it beyond the corner round's top.
            (["1847/100", "1896/100", "17/2"], Inside),
            (["1877/100", "1936/100", "199/20"], Outside),
            (["10", "10", "5"], Inside),
        ],
    );
}

// The 20 x 20 x 10 block with a 10 x 10 pocket 6 deep from z = 4, its four
// uprights rounded by `upright`, then its floor by `floor`.
fn pocket(upright: &str, floor: &str) -> Recipe {
    recipe(
        4,
        vec![
            block("block", ["20", "20", "10"], ["0", "0", "0"]),
            block("tool", ["10", "10", "10"], ["5", "5", "4"]),
            json!({"id":"pocket","op":"boolean","kind":"subtract","target":"block","tool":"tool"}),
            fillet("upright", "pocket", upright, uprights("tool")),
            fillet(
                "floor",
                "upright",
                floor,
                vec![edge(cap("tool", "start"), side("tool", "x_max"))],
            ),
        ],
    )
}

#[test]
fn a_pocket_floor_rolls_over_its_rounded_corners() {
    let big = r("3");
    for floor in ["1", "3"] {
        let rho = r(floor);
        let built = build(&pocket("3", floor));
        let corner = if floor == "1" { "Torus" } else { "Sphere" };
        expect_kinds(&built, &[("Plane", 11), ("Cylinder", 8), (corner, 4)]);
        // 3400 plus the four upright sections over the depth 6, the floor
        // section along runs 10 - 2R, and four quarter turns inside R.
        let runs = &r("40") - &(&r("8") * &big);
        let [s0, s1] = section(&rho);
        let [u0, u1] = section(&big);
        let [t1, t2] = quarter_turn(&rho, &big, -1);
        let four = r("4");
        exact_volume(
            &built,
            [
                &(&r("3400") + &(&(&four * &u0) * &r("6"))) + &(&s0 * &runs),
                &(&(&(&four * &u1) * &r("6")) + &(&s1 * &runs)) + &(&four * &t1),
                &four * &t2,
            ],
            "1e-6",
        );
        closed_display(&built);
        step(
            &built,
            if floor == "1" {
                "TOROIDAL_SURFACE("
            } else {
                "SPHERICAL_SURFACE("
            },
        );
        assert_eq!(names(&built, "floor").len(), 1);
        assert_eq!(faces_of(&built, "floor").len(), 8);
    }
    assert_eq!(
        names(&build(&pocket("3", "1")), "floor"),
        names(&build(&pocket("2", "2")), "floor")
    );
    use Containment::{Inside, Outside};
    let mut built = build(&pocket("3", "1"));
    classify(
        &mut built,
        &[
            // Added along the floor's right edge, and the ball's own space.
            (["299/20", "10", "81/20"], Inside),
            (["29/2", "10", "9/2"], Outside),
            // Added round the corner about (12, 12): 59/20 from its axis and
            // 1/20 above the floor, outside the tube; then inside the tube.
            (["1377/100", "1436/100", "81/20"], Inside),
            (["27/2", "14", "9/2"], Outside),
            (["10", "10", "7"], Outside),
            (["2", "2", "2"], Inside),
        ],
    );
    let mut built = build(&pocket("3", "3"));
    classify(
        &mut built,
        &[
            // Outside the sphere about (12, 12, 7) in the corner, and inside it.
            (["1377/100", "1436/100", "81/20"], Inside),
            (["27/2", "14", "6"], Outside),
        ],
    );
}

// The 20 x 20 x 10 block with its top front edge rounded by `front`, then a
// right-hand edge by `right`: the top one, or the round's own end arc.
fn corner(front: &str, right: &str, arc: bool) -> Recipe {
    let selected = if arc {
        edge(made("front"), side("block", "x_max"))
    } else {
        edge(cap("block", "end"), side("block", "x_max"))
    };
    recipe(
        4,
        vec![
            block("block", ["20", "20", "10"], ["0", "0", "0"]),
            fillet(
                "front",
                "block",
                front,
                vec![edge(cap("block", "end"), side("block", "y_min"))],
            ),
            fillet("right", "front", right, vec![selected]),
        ],
    )
}

#[test]
fn a_box_corner_rolls_over_an_earlier_edge_round() {
    let big = r("3");
    for (right, arc) in [("1", false), ("1", true), ("3", false), ("2", false)] {
        let rho = r(right);
        let built = build(&corner("3", right, arc));
        let carrier = match right {
            "1" => "Torus",
            "2" => "TorusSheet",
            _ => "Sphere",
        };
        expect_kinds(&built, &[("Plane", 6), ("Cylinder", 3), (carrier, 1)]);
        // 4000 less the front round over 20, the right chain's straight runs
        // 20 - R along the top and 10 - R down the front, and a quarter turn
        // inside R about the front round's axis.
        let runs = &r("30") - &(&r("2") * &big);
        let [s0, s1] = section(&rho);
        let [f0, f1] = section(&big);
        let [t1, t2] = quarter_turn(&rho, &big, -1);
        let twenty = r("20");
        exact_volume(
            &built,
            [
                &(&r("4000") - &(&twenty * &f0)) - &(&s0 * &runs),
                &(&(-&(&twenty * &f1)) - &(&s1 * &runs)) - &t1,
                -&t2,
            ],
            "1e-6",
        );
        closed_display(&built);
        step(&built, "CYLINDRICAL_SURFACE(");
        assert_eq!(names(&built, "right").len(), 1);
        assert_eq!(faces_of(&built, "right").len(), 3);
    }
    // The chain keeps the name of whichever edge was selected.
    assert_ne!(
        names(&build(&corner("3", "1", false)), "right"),
        names(&build(&corner("3", "1", true)), "right")
    );
    for arc in [false, true] {
        assert_eq!(
            names(&build(&corner("3", "1", arc)), "right"),
            names(&build(&corner("2", "2", arc)), "right")
        );
    }
    use Containment::{Inside, Outside};
    let mut built = build(&corner("3", "1", false));
    classify(
        &mut built,
        &[
            (["199/10", "15", "99/10"], Outside),
            (["39/2", "15", "19/2"], Inside),
            // About the front round's axis through (y, z) = (3, 7): 59/20 out
            // toward the top front, 1/20 inside the right face, outside the
            // tube; then inside it.
            (["399/20", "123/100", "936/100"], Outside),
            (["39/2", "3/2", "9"], Inside),
            (["199/10", "1/10", "3"], Outside),
            (["39/2", "1/2", "3"], Inside),
            (["10", "10", "5"], Inside),
        ],
    );
    let mut built = build(&corner("2", "2", false));
    classify(
        &mut built,
        &[
            // The sphere about (18, 2, 8).
            (["199/10", "43/50", "238/25"], Outside),
            (["19", "11/10", "46/5"], Inside),
        ],
    );
    // A third round rolls over the second: the bottom right edge continues
    // over the foot of the right chain's upright, a spindle sheet whose
    // major radius equals its tube's, and along the bottom front edge.
    let mut steps = vec![
        block("block", ["20", "20", "10"], ["0", "0", "0"]),
        fillet(
            "front",
            "block",
            "3",
            vec![edge(cap("block", "end"), side("block", "y_min"))],
        ),
        fillet(
            "right",
            "front",
            "2",
            vec![edge(cap("block", "end"), side("block", "x_max"))],
        ),
    ];
    steps.push(fillet(
        "bottom",
        "right",
        "1",
        vec![edge(cap("block", "start"), side("block", "x_max"))],
    ));
    let built = build(&recipe(4, steps));
    expect_kinds(&built, &[("Plane", 6), ("Cylinder", 5), ("TorusSheet", 2)]);
    // 4000 less the front round over 20; the right chain at radius 2 along
    // 17 and 7 with a quarter turn inside 3; the bottom chain at radius 1
    // along 18 and 18 with a quarter turn inside 2.
    let (one, two) = (r("1"), r("2"));
    let [f0, f1] = section(&big);
    let [m0, m1] = section(&two);
    let [b0, b1] = section(&one);
    let [p1, p2] = quarter_turn(&two, &big, -1);
    let [q1, q2] = quarter_turn(&one, &two, -1);
    let sum = |terms: &[Rational]| terms.iter().fold(Rational::zero(), |a, b| &a + b);
    exact_volume(
        &built,
        [
            sum(&[
                r("4000"),
                -&(&r("20") * &f0),
                -&(&r("24") * &m0),
                -&(&r("36") * &b0),
            ]),
            sum(&[
                -&(&r("20") * &f1),
                -&(&r("24") * &m1),
                -&(&r("36") * &b1),
                -&p1,
                -&q1,
            ]),
            -&(&p2 + &q2),
        ],
        "1e-6",
    );
    closed_display(&built);
    step(&built, "TOROIDAL_SURFACE(");
}

// An L: a 20 x 20 x 5 base and a 20 x 5 x 10 wall along its front, the
// inside corner rounded by `inner`, then the right edge of the base's top by
// `right`, which rolls over the inner round's end and up the wall.
fn bracket(inner: &str, right: &str) -> Recipe {
    recipe(
        4,
        vec![
            block("base", ["20", "20", "5"], ["0", "0", "0"]),
            block("wall", ["20", "5", "10"], ["0", "0", "0"]),
            json!({"id":"l","op":"boolean","kind":"union","target":"base","tool":"wall"}),
            fillet(
                "inner",
                "l",
                inner,
                vec![edge(cap("base", "end"), side("wall", "y_max"))],
            ),
            fillet(
                "right",
                "inner",
                right,
                vec![edge(cap("base", "end"), side("base", "x_max"))],
            ),
        ],
    )
}

#[test]
fn a_convex_round_rolls_over_a_concave_one() {
    let big = r("3");
    for right in ["1", "4"] {
        let rho = r(right);
        let built = build(&bracket("3", right));
        expect_kinds(&built, &[("Plane", 8), ("Cylinder", 3), ("Torus", 1)]);
        // 2500 plus the inner round over 20, less the right chain's runs
        // 15 - R along the base and 5 - R up the wall, and a quarter turn
        // outside R about the inner round's axis.
        let runs = &r("20") - &(&r("2") * &big);
        let [s0, s1] = section(&rho);
        let [i0, i1] = section(&big);
        let [t1, t2] = quarter_turn(&rho, &big, 1);
        let twenty = r("20");
        exact_volume(
            &built,
            [
                &(&r("2500") + &(&twenty * &i0)) - &(&s0 * &runs),
                &(&(&twenty * &i1) - &(&s1 * &runs)) - &t1,
                -&t2,
            ],
            "1e-6",
        );
        closed_display(&built);
        step(&built, "TOROIDAL_SURFACE(");
        assert_eq!(faces_of(&built, "right").len(), 3);
    }
    assert_eq!(
        names(&build(&bracket("3", "1")), "right"),
        names(&build(&bracket("2", "3/2")), "right")
    );
    use Containment::{Inside, Outside};
    let mut built = build(&bracket("3", "1"));
    classify(
        &mut built,
        &[
            (["199/10", "12", "49/10"], Outside),
            (["39/2", "12", "9/2"], Inside),
            // About the inner round's axis through (y, z) = (8, 8): 61/20 out
            // toward the inside corner, 1/20 inside the right face, outside
            // the tube of radius 1 about the circle of radius 4; then inside.
            (["399/20", "617/100", "139/25"], Outside),
            (["39/2", "59/10", "26/5"], Inside),
            (["199/10", "49/10", "9"], Outside),
            (["10", "10", "2"], Inside),
        ],
    );
}

fn refused(recipe: &Recipe, step: &str) -> catallax_numbers::Refusal {
    match evaluate(recipe, &mut decider()) {
        Err(BuildError::Geometry { step: at, refusal }) => {
            assert_eq!(at.as_str(), step);
            refusal
        }
        Ok(_) => panic!("step {step} built where a refusal was expected"),
        Err(other) => panic!("{other}"),
    }
}

#[test]
fn fillets_over_fillets_not_built_yet_refuse_with_typed_codes() {
    use catallax_features::fillet::{RADIUS_TOO_LARGE, radius_limit};
    // A ball wider than a round cannot roll along it: an edge of that round
    // chosen by name refuses, with that round's radius as the limit.
    for recipe in [boss(4, "1", "3", true), corner("1", "2", true)] {
        let refusal = refused(&recipe, recipe.steps.last().unwrap().id.as_str());
        assert_eq!(refusal.code, RADIUS_TOO_LARGE, "{}", refusal.detail);
        assert_eq!(radius_limit(&refusal), Some(r("1")));
    }
    // Without roll-over, the strip's carrier passes the earlier round into
    // the far faces, and an upright's carrier passes an earlier mitre.
    for (recipe, step) in [
        (boss(3, "3", "1", false), "top"),
        (mitred(3, "3", "3", false), "upright"),
    ] {
        let refusal = refused(&recipe, step);
        assert_eq!(
            refusal.code, "fillet.unsupported-overflow",
            "{}",
            refusal.detail
        );
    }
    // An earlier round's tangent boundary is not a sharp edge.
    let seam = recipe(
        4,
        vec![
            block("block", ["20", "20", "10"], ["0", "0", "0"]),
            fillet(
                "front",
                "block",
                "3",
                vec![edge(cap("block", "end"), side("block", "y_min"))],
            ),
            fillet(
                "seam",
                "front",
                "1",
                vec![edge(cap("block", "end"), made("front"))],
            ),
        ],
    );
    let refusal = refused(&seam, "seam");
    assert_eq!(
        refusal.code, "fillet.unsupported-topology",
        "{}",
        refusal.detail
    );
    assert!(
        refusal.detail.contains("tangent faces"),
        "{}",
        refusal.detail
    );
    // An earlier round meeting a leaning end face ends on an ellipse, and a
    // chain over that elliptic arc has no exact join yet.
    let leaning = recipe(
        4,
        vec![
            json!({"id":"profile","op":"sketch",
                "plane":{"rotation":["1","1","0","0"],"translation":["0","0","0"]},
                "loops":[[
                    {"kind":"line","id":"floor","to":["20","0"]},
                    {"kind":"line","id":"lean","to":["15","10"]},
                    {"kind":"line","id":"roof","to":["0","10"]},
                    {"kind":"line","id":"back","to":["0","0"]}]]}),
            json!({"id":"body","op":"extrude","profile":"profile","distance":"20"}),
            fillet(
                "front",
                "body",
                "3",
                vec![edge(side("body", "roof"), cap("body", "start"))],
            ),
            fillet(
                "right",
                "front",
                "1",
                vec![edge(side("body", "roof"), side("body", "lean"))],
            ),
        ],
    );
    let refusal = refused(&leaning, "right");
    assert_eq!(
        refusal.code,
        catallax_features::fillet::UNSUPPORTED_CORNER,
        "{}",
        refusal.detail
    );
    // The field compiler has no roll-over over earlier strips.
    for recipe in [
        boss(4, "3", "1", false),
        pocket("3", "3"),
        mitred(4, "3", "3", false),
    ] {
        let Err(refusal) = catallax_field::Plan::compile_checked(&recipe, &mut decider()) else {
            panic!("the field compiled a fillet over a fillet");
        };
        assert_eq!(
            refusal.code, "fillet.unsupported-overflow",
            "{}",
            refusal.detail
        );
    }
}

/// A box `h` long with an `x` by `y` section whose four long edges are
/// rounded by `radii`, and every edge of both ends by `r`, at most each of
/// them: the middle section over `h - 2r`, and at each end a section offset
/// inward by `r - sqrt(r^2 - s^2)` at depth `s`. As `[1, pi, pi^2]`.
fn rounded_box(h: &str, [x, y]: [&str; 2], radii: [&str; 4], r_end: &str) -> [Rational; 3] {
    let (h, x, y, rho) = (r(h), r(x), r(y), r(r_end));
    let two = r("2");
    let (x1, y1) = (&x - &(&two * &rho), &y - &(&two * &rho));
    let long = &h - &(&two * &rho);
    let r2 = &rho * &rho;
    let r3 = &r2 * &rho;
    let mut squares = Rational::zero();
    let mut s0 = Rational::zero();
    let mut s1 = Rational::zero();
    for radius in radii.map(r) {
        squares = &squares + &(&radius * &radius);
        let gap = &radius - &rho;
        s0 = &(&s0 + &(&(&gap * &gap) * &rho)) + &(&(&two * &r3) / &r("3"));
        s1 = &s1 + &(&(&gap * &r2) / &two);
    }
    let sum = |terms: &[Rational]| terms.iter().fold(Rational::zero(), |a, b| &a + b);
    [
        sum(&[
            &long * &(&x * &y),
            -&(&long * &squares),
            &two * &(&(&x1 * &y1) * &rho),
            &(&r("16") * &r3) / &r("3"),
            -&(&two * &s0),
        ]),
        sum(&[
            &(&long * &squares) / &r("4"),
            &(&x1 + &y1) * &r2,
            -&(&two * &s1),
            &s0 / &two,
        ]),
        &s1 / &two,
    ]
}

#[test]
fn every_edge_rounds_over_an_earlier_round() {
    // The smooth rails of the earlier round are not edges to round.
    let every = |front: &str, rest: &str| {
        recipe(
            4,
            vec![
                block("block", ["20", "20", "10"], ["0", "0", "0"]),
                fillet(
                    "front",
                    "block",
                    front,
                    vec![edge(cap("block", "end"), side("block", "y_min"))],
                ),
                json!({"id":"rest","op":"fillet","body":"front","radius":rest,"edges":"all"}),
            ],
        )
    };
    for (front, rest) in [("3", "1"), ("2", "2")] {
        let built = build(&every(front, rest));
        let ends = if front == rest { "Sphere" } else { "Torus" };
        let mut want = BTreeMap::from([
            ("Plane".to_string(), 6),
            ("Cylinder".to_string(), 12),
            ("Sphere".to_string(), 6),
        ]);
        *want.entry(ends.to_string()).or_default() += 2;
        assert_eq!(kinds(&built), want);
        // Slabs along x: the yz section has the front round and three
        // rounds of the rest.
        exact_volume(
            &built,
            rounded_box("20", ["20", "10"], [front, rest, rest, rest], rest),
            "1e-6",
        );
        closed_display(&built);
        step(&built, "SPHERICAL_SURFACE(");
    }
    // Rounded uprights, then every edge: both loops roll over all four.
    let boss = recipe(
        4,
        vec![
            block("block", ["20", "20", "10"], ["0", "0", "0"]),
            fillet("upright", "block", "3", uprights("block")),
            json!({"id":"rest","op":"fillet","body":"upright","radius":"1","edges":"all"}),
        ],
    );
    let built = build(&boss);
    expect_kinds(&built, &[("Plane", 6), ("Cylinder", 12), ("Torus", 8)]);
    exact_volume(
        &built,
        rounded_box("10", ["20", "20"], ["3", "3", "3", "3"], "1"),
        "1e-6",
    );
    closed_display(&built);
    // A top loop with mitred corners, then every edge: each upright closes
    // with its corner and the mitre above it goes without a round of its own.
    let mitred = recipe(
        4,
        vec![
            block("block", ["20", "20", "10"], ["0", "0", "0"]),
            fillet(
                "top",
                "block",
                "2",
                SIDES
                    .iter()
                    .map(|s| edge(cap("block", "end"), side("block", s)))
                    .collect(),
            ),
            json!({"id":"rest","op":"fillet","body":"top","radius":"2","edges":"all"}),
        ],
    );
    let built = build(&mitred);
    expect_kinds(&built, &[("Plane", 6), ("Cylinder", 12), ("Sphere", 8)]);
    exact_volume(
        &built,
        rounded_box("10", ["20", "20"], ["2", "2", "2", "2"], "2"),
        "1e-6",
    );
    closed_display(&built);
}

// The block's top edges along x_max and y_max, or its whole top loop, rounded
// together with mitred corners, then the upright under each mitre.
fn mitred(format: u8, top: &str, upright: &str, whole: bool) -> Recipe {
    let (tops, ups): (Vec<_>, Vec<_>) = if whole {
        (
            SIDES
                .iter()
                .map(|s| edge(cap("block", "end"), side("block", s)))
                .collect(),
            uprights("block"),
        )
    } else {
        (
            vec![
                edge(cap("block", "end"), side("block", "x_max")),
                edge(cap("block", "end"), side("block", "y_max")),
            ],
            vec![edge(side("block", "x_max"), side("block", "y_max"))],
        )
    };
    recipe(
        format,
        vec![
            block("block", ["20", "20", "10"], ["0", "0", "0"]),
            fillet("top", "block", top, tops),
            fillet("upright", "top", upright, ups),
        ],
    )
}

#[test]
fn an_upright_round_swallows_an_earlier_mitre() {
    // The upright's ball reaches both top strips at once and closes the
    // corner they would have closed together; the mitre between them goes.
    for rho in ["3", "2"] {
        let built = build(&mitred(4, rho, rho, false));
        expect_kinds(&built, &[("Plane", 6), ("Cylinder", 3), ("Sphere", 1)]);
        // 4000 less each section over its run up to the corner, 20 - r, 20 - r
        // and 10 - r, and the corner cube less its ball octant.
        let rho = r(rho);
        let [s0, s1] = section(&rho);
        let runs = &r("50") - &(&r("3") * &rho);
        let cube = &(&rho * &rho) * &rho;
        exact_volume(
            &built,
            [
                &(&r("4000") - &(&s0 * &runs)) - &cube,
                &(-&(&s1 * &runs)) + &(&cube / &r("6")),
                r("0"),
            ],
            "1e-6",
        );
        closed_display(&built);
        step(&built, "SPHERICAL_SURFACE(");
        let corners: Vec<_> = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .filter(|f| {
                built
                    .model
                    .face_origins(*f)
                    .iter()
                    .any(|o| matches!(o.role, Role::Corner { .. }))
            })
            .collect();
        assert_eq!(corners.len(), 1);
    }
    let corner_names = |recipe: &Recipe| -> BTreeSet<Origin> {
        let built = build(recipe);
        built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .flat_map(|f| built.model.face_origins(f).to_vec())
            .filter(|o| matches!(o.role, Role::Corner { .. }))
            .collect()
    };
    assert_eq!(
        corner_names(&mitred(4, "3", "3", false)),
        corner_names(&mitred(4, "2", "2", false))
    );
    // The whole top loop, then all four uprights: four corners.
    let built = build(&mitred(4, "3", "3", true));
    expect_kinds(&built, &[("Plane", 6), ("Cylinder", 8), ("Sphere", 4)]);
    exact_volume(&built, [r("3136"), r("207"), r("0")], "1e-6");
    closed_display(&built);
    use Containment::{Inside, Outside};
    let mut built = build(&mitred(4, "3", "3", false));
    classify(
        &mut built,
        &[
            // The ball about (17, 17, 7): inside it, and outside it where the
            // mitre stood.
            (["37/2", "37/2", "17/2"], Inside),
            (["199/10", "199/10", "9"], Outside),
            (["19", "19", "99/10"], Outside),
            (["10", "10", "5"], Inside),
        ],
    );
}

// A pocket's floor loop rounded with mitred corners, then its uprights.
fn mitred_pocket(floor: &str, upright: &str) -> Recipe {
    recipe(
        4,
        vec![
            block("block", ["20", "20", "10"], ["0", "0", "0"]),
            block("tool", ["10", "10", "10"], ["5", "5", "4"]),
            json!({"id":"pocket","op":"boolean","kind":"subtract","target":"block","tool":"tool"}),
            fillet(
                "floor",
                "pocket",
                floor,
                SIDES
                    .iter()
                    .map(|s| edge(cap("tool", "start"), side("tool", s)))
                    .collect(),
            ),
            fillet("upright", "floor", upright, uprights("tool")),
        ],
    )
}

#[test]
fn pocket_uprights_swallow_the_floor_loop_mitres() {
    let built = build(&mitred_pocket("1", "1"));
    expect_kinds(&built, &[("Plane", 11), ("Cylinder", 8), ("Sphere", 4)]);
    // The pocket with every inner edge rounded at once.
    exact_volume(&built, [r("3456"), -&r("41/3"), r("0")], "1e-6");
    closed_display(&built);
    step(&built, "SPHERICAL_SURFACE(");
}

// A 10 x 10 x 10 boss on a 30 x 30 x 5 plate, its uprights rounded by
// `upright`, then its base by `base`: a concave loop over convex rounds.
fn based(upright: &str, base: &str) -> Recipe {
    recipe(
        4,
        vec![
            block("plate", ["30", "30", "5"], ["0", "0", "0"]),
            block("boss", ["10", "10", "10"], ["10", "10", "5"]),
            json!({"id":"joined","op":"boolean","kind":"union","target":"plate","tool":"boss"}),
            fillet("upright", "joined", upright, uprights("boss")),
            fillet(
                "base",
                "upright",
                base,
                vec![edge(cap("plate", "end"), side("boss", "x_max"))],
            ),
        ],
    )
}

#[test]
fn a_boss_base_rolls_over_its_rounded_uprights() {
    let big = r("2");
    for base in ["1", "3"] {
        let rho = r(base);
        let built = build(&based("2", base));
        expect_kinds(&built, &[("Plane", 11), ("Cylinder", 8), ("Torus", 4)]);
        // 5500 less the four upright sections over the boss's height 10, plus
        // the base section along runs 10 - 2R and four quarter turns outside R.
        let runs = &r("40") - &(&r("8") * &big);
        let [s0, s1] = section(&rho);
        let [u0, u1] = section(&big);
        let [t1, t2] = quarter_turn(&rho, &big, 1);
        let four = r("4");
        exact_volume(
            &built,
            [
                &(&r("5500") - &(&(&four * &u0) * &r("10"))) + &(&s0 * &runs),
                &(&(-&(&(&four * &u1) * &r("10"))) + &(&s1 * &runs)) + &(&four * &t1),
                &four * &t2,
            ],
            "1e-6",
        );
        closed_display(&built);
        step(&built, "TOROIDAL_SURFACE(");
        assert_eq!(faces_of(&built, "base").len(), 8);
    }
    assert_eq!(
        names(&build(&based("2", "1")), "base"),
        names(&build(&based("3", "1")), "base")
    );
    use Containment::{Inside, Outside};
    let mut built = build(&based("2", "1"));
    classify(
        &mut built,
        &[
            // Added along the boss's right wall, and outside the ball there.
            (["201/10", "15", "51/10"], Inside),
            (["41/2", "15", "11/2"], Outside),
            // Added round the corner about (18, 18), 41/20 from its axis and
            // 1/20 above the plate, outside the tube about the circle of
            // radius 3; then inside the tube.
            (["1923/100", "1964/100", "101/20"], Inside),
            (["39/2", "20", "11/2"], Outside),
            (["15", "15", "10"], Inside),
        ],
    );
}

/// The block in a rational pose that turns it a third of a turn about the
/// diagonal (1, 1, 1) and moves it off the origin.
fn placed(recipe: &Recipe) -> Recipe {
    let mut value = serde_json::to_value(recipe).unwrap();
    for step in value["steps"].as_array_mut().unwrap() {
        if step["op"] == "block" {
            step["pose"] = json!({"rotation":["3","1","1","1"],"translation":["7","-2","5"]});
        }
    }
    Recipe::parse(&value.to_string()).unwrap()
}

#[test]
fn fillets_over_fillets_build_in_any_rational_placement() {
    for (recipe, kind) in [
        (boss(4, "3", "1", false), "Torus"),
        (boss(4, "3", "3", false), "Sphere"),
        (corner("3", "1", false), "Torus"),
        (mitred(4, "3", "3", false), "Sphere"),
    ] {
        let level = build(&recipe);
        let turned = build(&placed(&recipe));
        assert_eq!(kinds(&turned), kinds(&level), "{kind}");
        let volume = |built: &Built| {
            catallax_measure::volume(
                &built.model,
                built.result,
                &Default::default(),
                &mut measuring(),
            )
            .unwrap()
            .value
        };
        let (a, b) = (volume(&level), volume(&turned));
        assert!(
            a.lo <= b.hi && b.lo <= a.hi && &b.hi - &b.lo <= r("1e-6"),
            "{kind}: {a:?} {b:?}"
        );
        closed_display(&turned);
        assert_eq!(names(&turned, "top"), names(&level, "top"));
    }
}

#[test]
fn every_edge_leaves_out_an_oblique_bands_rails_by_its_definition() {
    // A cylinder of radius 5 cut off by a tilted plane, its oblique rim
    // rounded by an exact rational band, then every edge: the band's rails
    // are smooth by its retained definition, so only the bottom rim rounds.
    let tilt = json!(["24/25", "0", "7/25", "0"]);
    let steps = |every: bool| {
        let mut steps = vec![
            json!({"id":"wall","op":"cylinder","radius":"5","height":"16","pose":pose(["0","0","0"])}),
            json!({"id":"cut-section","op":"sketch",
                "plane":{"rotation":tilt,"translation":["0","0","8"]},
                "loops":[[{"kind":"line","id":"s","to":["30","-30"]},
                    {"kind":"line","id":"e","to":["30","30"]},
                    {"kind":"line","id":"n","to":["-30","30"]},
                    {"kind":"line","id":"w","to":["-30","-30"]}]]}),
            json!({"id":"cut","op":"extrude","profile":"cut-section","distance":"30"}),
            json!({"id":"joined","op":"boolean","kind":"subtract","target":"wall","tool":"cut"}),
            fillet(
                "round",
                "joined",
                "1",
                vec![edge(side("wall", "wall"), cap("cut", "start"))],
            ),
        ];
        if every {
            steps.push(
                json!({"id":"rest","op":"fillet","body":"round","radius":"1/2","edges":"all"}),
            );
        }
        recipe(4, steps)
    };
    let (band, rounded) = (build(&steps(false)), build(&steps(true)));
    let kinds_of = |built: &Built| kinds(built);
    let mut want = kinds_of(&band);
    *want.entry("Torus".to_string()).or_default() += 1;
    assert_eq!(kinds_of(&rounded), want);
    assert_eq!(faces_of(&rounded, "rest").len(), 1);
    // The bottom ring removes a full turn of its section inside radius 5.
    let volume = |built: &Built| {
        catallax_measure::volume(
            &built.model,
            built.result,
            &Default::default(),
            &mut measuring(),
        )
        .unwrap()
        .value
    };
    let (before, after) = (volume(&band), volume(&rounded));
    let [t1, t2] = quarter_turn(&r("1/2"), &r("5"), -1);
    let (lo, hi) = pi();
    let ring = |p: &Rational| &r("4") * &(&(&t1 * p) + &(&(&t2 * p) * p));
    let [a, b] = [ring(&lo), ring(&hi)];
    let removed_lo = &before.lo - &after.hi;
    let removed_hi = &before.hi - &after.lo;
    assert!(
        removed_lo <= a.clone().max(b.clone()) && removed_hi >= a.min(b),
        "removed [{}, {}]",
        removed_lo.to_f64_nearest(),
        removed_hi.to_f64_nearest()
    );
    closed_display(&rounded);
}

#[test]
fn a_wider_upright_takes_the_cliff_at_an_earlier_mitre() {
    // Radius-2 uprights after a radius-1 top loop: the ball cannot roll onto
    // the narrower rounds, so each upright keeps its planes, runs on to the
    // top and the rounds cut it (recipe section 4.5.1, the cliff rule).
    let mut built = build(&mitred(4, "1", "2", true));
    expect_kinds(&built, &[("Plane", 6), ("Cylinder", 8)]);
    assert_eq!(faces_of(&built, "upright").len(), 4);
    // The block with its top loop, 4000 - 220/3 + 18 pi, less at each corner
    // the section 4 - pi over nine units and, above, the section of the
    // inset square's corner outside the upright's disk. That integral has no
    // closed form here; 50-digit quadrature of it gives this volume.
    near_volume(&built, "3950.54249556093264192891406600233053869536782548");
    use Containment::{Inside, Outside};
    classify(
        &mut built,
        &[
            (["10", "10", "5"], Inside),
            // Removed below the cliff, as on any upright.
            (["3/10", "3/10", "5"], Outside),
            // Under the x_min round, outside the upright's disk: cut away.
            (["1/5", "1", "46/5"], Outside),
            // Inside both the round and the disk: kept.
            (["1/2", "3/2", "19/2"], Inside),
            // The mitre's outer part survives: a sphere corner would remove it.
            (["9/10", "9/10", "199/20"], Inside),
            (["3/2", "3/2", "999/100"], Inside),
        ],
    );
    closed_display(&built);
    step(&built, "CYLINDRICAL_SURFACE(");
    // One upright after two radius-3 top edges, at radius 4: the block less
    // those edges, 3685 + 153 pi / 2, less the upright's section 16 - 4 pi
    // over seven units and its quadrature under the rounds.
    let built = build(&mitred(4, "3", "4", false));
    expect_kinds(&built, &[("Plane", 6), ("Cylinder", 3)]);
    near_volume(&built, "3896.86178057071696975579481666599842658805044128");
    closed_display(&built);
    // Concave: a pocket's floor loop at radius 1, then its uprights at 3/2.
    // The floor adds 100/3 - 8 pi to the 3400 left by the pocket; each
    // upright adds its section over the five units above the floor round
    // and its quadrature inside the round's inset square.
    let built = build(&mitred_pocket("1", "3/2"));
    expect_kinds(&built, &[("Plane", 11), ("Cylinder", 8)]);
    near_volume(&built, "3418.73476474139205170497470828914411541290912172");
    closed_display(&built);
    step(&built, "CYLINDRICAL_SURFACE(");
}

/// The measured volume holds a value computed independently to 50 digits.
fn near_volume(built: &Built, expected: &str) {
    let volume = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut measuring(),
    )
    .unwrap()
    .value;
    let expected = r(expected);
    let slack = r("1e-9");
    assert!(
        volume.lo <= &expected + &slack && volume.hi >= &expected - &slack,
        "[{}, {}]",
        volume.lo.to_f64_nearest(),
        volume.hi.to_f64_nearest()
    );
    assert!(&volume.hi - &volume.lo <= r("1e-6"), "{volume:?}");
}

#[test]
fn a_chain_takes_the_cliff_at_a_narrower_round() {
    // The chain does not roll onto a round narrower than its ball: the strip
    // keeps its planes, runs on and the faces past the round cut it, as in
    // recipe section 4.5.1. Each volume is the input less (or, concave, plus)
    // the strip's section over its full run and a quadrature of the section
    // the round shrinks, computed independently to 50 digits.
    let built = build(&boss(4, "1", "3", false));
    expect_kinds(&built, &[("Plane", 6), ("Cylinder", 5)]);
    near_volume(&built, "3953.63873405415624765806822306421993832160480262");
    closed_display(&built);
    step(&built, "CYLINDRICAL_SURFACE(");
    // Concave: a pocket's floor edge at radius 3 after its uprights at 1.
    let built = build(&pocket("1", "3"));
    expect_kinds(&built, &[("Plane", 11), ("Cylinder", 5)]);
    near_volume(&built, "3423.61347126597327316119182659268688810962018647");
    closed_display(&built);
    // One end on a narrower round along the edge, the other on a plane.
    let mut built = build(&corner("1", "2", false));
    expect_kinds(&built, &[("Plane", 6), ("Cylinder", 2)]);
    near_volume(&built, "3978.79331451223336740924776465736081221666311832");
    closed_display(&built);
    step(&built, "CYLINDRICAL_SURFACE(");
    use Containment::{Inside, Outside};
    classify(
        &mut built,
        &[
            (["10", "10", "5"], Inside),
            // Cut by the strip's extension below the front round.
            (["199/10", "1/2", "93/10"], Outside),
            // Inside the strip's disk under the front round: kept.
            (["19", "1/2", "9"], Inside),
            (["19", "10", "99/10"], Outside),
        ],
    );
}

#[test]
fn a_narrower_upright_rolls_over_an_earlier_mitre() {
    // Radius-1/2 uprights after a radius-1 top loop: the ball is narrower
    // than the rounds, so it rolls over each mitre to its apex (recipe
    // section 4.5.1) along an exact pipe around the rounds' offset ellipse.
    let mut built = build(&mitred(4, "1", "1/2", true));
    expect_kinds(&built, &[("Plane", 6), ("Cylinder", 8), ("BSpline", 4)]);
    // The block with its top loop, 4000 - 220/3 + 18 pi, less the uprights'
    // section (1 - pi/4)/4 over nine units and, at each corner, the normal
    // sections between the ball and the rounds integrated along the pipe's
    // spine: 45-digit quadrature, checked by an opening computation.
    near_volume(&built, "3981.14478206068032165465886714390960314813664");
    use Containment::{Inside, Outside};
    classify(
        &mut built,
        &[
            (["10", "10", "5"], Inside),
            (["1/10", "1/10", "5"], Outside),
            (["3/10", "3/10", "5"], Inside),
            // On the pipe's symmetry plane, inside and outside its section.
            (["92547/200000", "92547/200000", "9720977/1000000"], Inside),
            (
                ["425993/1000000", "425993/1000000", "4897231/500000"],
                Outside,
            ),
            // Under the x round, short of the contact rail: kept.
            (
                ["646447/1000000", "211261/1000000", "2395407/250000"],
                Inside,
            ),
            // 0.002 inside the pipe, where OCCT's corner leaves no material.
            (["46843/200000", "46843/200000", "590191/62500"], Inside),
            (
                ["197913/1000000", "197913/1000000", "9341479/1000000"],
                Inside,
            ),
        ],
    );
    closed_display(&built);
    step(&built, "B_SPLINE_SURFACE");
    // One upright at radius 2 after two radius-3 top edges: the block less
    // those edges, 3685 + 153 pi / 2, less the upright's section 4 - pi over
    // seven units and the corner's band integral.
    let built = build(&mitred(4, "3", "2", false));
    expect_kinds(&built, &[("Plane", 6), ("Cylinder", 3), ("BSpline", 1)]);
    near_volume(&built, "3917.63762828103699951026186048602218403174678");
    closed_display(&built);
    // Concave: a pocket's uprights at 1/2 after its floor loop at radius 1.
    // The floor adds 100/3 - 8 pi to the 3400 left by the pocket, each
    // upright its section over the five units above the floor round, and each
    // corner the region between the ball and the rounds that the convex
    // corner removes.
    let built = build(&mitred_pocket("1", "1/2"));
    expect_kinds(&built, &[("Plane", 11), ("Cylinder", 8), ("BSpline", 4)]);
    near_volume(&built, "3409.41273712880740396843021007216492857803222");
    closed_display(&built);
    step(&built, "B_SPLINE_SURFACE");
}

#[test]
fn mirrored_narrower_uprights_keep_volume_material_and_step() {
    // A reflection turns each mitre band's chart normal toward its ball. The
    // material side, the volume, the display and STEP must not turn with it.
    use Containment::{Inside, Outside};
    use catallax_shapes::{Isometry, RVector3};
    let mirror =
        Isometry::reflection(&RPoint3::origin(), &RVector3::new(r("1"), r("2"), r("3"))).unwrap();
    let convex: Vec<([&str; 3], Containment)> = vec![
        (["10", "10", "5"], Inside),
        (["1/10", "1/10", "5"], Outside),
        (["92547/200000", "92547/200000", "9720977/1000000"], Inside),
        (
            ["425993/1000000", "425993/1000000", "4897231/500000"],
            Outside,
        ),
        (
            ["646447/1000000", "211261/1000000", "2395407/250000"],
            Inside,
        ),
        (["46843/200000", "46843/200000", "590191/62500"], Inside),
    ];
    // In the pocket's corner under both rounds, and on the floor round's
    // axis, clear of the upright round.
    let concave: Vec<([&str; 3], Containment)> = vec![
        (["2", "2", "5"], Inside),
        (["10", "10", "7"], Outside),
        (["51/10", "51/10", "41/10"], Inside),
        (["6", "6", "5"], Outside),
    ];
    for (source, volume, witnesses) in [
        (
            mitred(4, "1", "1/2", true),
            "3981.14478206068032165465886714390960314813664",
            convex,
        ),
        (
            mitred_pocket("1", "1/2"),
            "3409.41273712880740396843021007216492857803222",
            concave,
        ),
    ] {
        let mut doc = serde_json::to_value(&source).unwrap();
        let body = doc["result"].clone();
        doc["steps"].as_array_mut().unwrap().push(json!({
            "id":"mirrored", "op":"mirror", "body":body,
            "plane":{"point":["0","0","0"],"normal":["1","2","3"]},
            "keep_original":false
        }));
        doc["result"] = json!("mirrored");
        let mut built = build(&Recipe::parse(&doc.to_string()).unwrap());
        near_volume(&built, volume);
        for (p, want) in &witnesses {
            let at = mirror.apply_rational_point(&point(*p)).unwrap();
            let actual = catallax_boolean::classify_point(
                &mut built.model,
                built.result,
                &at,
                &mut decider(),
            );
            assert!(
                matches!(actual, Outcome::Proven(actual) if actual == *want),
                "mirrored B-rep at {p:?}: {actual:?}"
            );
        }
        closed_display(&built);
        step(&built, "B_SPLINE_SURFACE");
    }
}
