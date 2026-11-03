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
/// Exact measurements refine further than a build needs.
fn measuring() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 100_000_000,
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
        other => panic!("{other:?}"),
    }
}

#[test]
fn fillets_over_fillets_not_built_yet_refuse_with_typed_codes() {
    use catallax_features::fillet::{RADIUS_TOO_LARGE, radius_limit};
    // A ball wider than the round it rolls over cannot rest on it: the limit
    // is that round's radius, convex or concave.
    for recipe in [
        boss(4, "1", "3", false),
        pocket("1", "3"),
        corner("1", "2", false),
    ] {
        let refusal = refused(&recipe, recipe.steps.last().unwrap().id.as_str());
        assert_eq!(refusal.code, RADIUS_TOO_LARGE, "{}", refusal.detail);
        assert_eq!(radius_limit(&refusal), Some(r("1")));
    }
    // Without roll-over, the strip's carrier passes the earlier round into
    // the far faces.
    let refusal = refused(&boss(3, "3", "1", false), "top");
    assert_eq!(
        refusal.code, "fillet.unsupported-overflow",
        "{}",
        refusal.detail
    );
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
    // A ball as wide as a cylinder consumes the face inside its rim.
    let dome = recipe(
        4,
        vec![
            json!({"id":"post","op":"cylinder","radius":"3","height":"10","pose":pose(["0","0","0"])}),
            fillet(
                "dome",
                "post",
                "3",
                vec![edge(cap("post", "end"), side("post", "wall"))],
            ),
        ],
    );
    let refusal = refused(&dome, "dome");
    assert_eq!(
        refusal.code, "fillet.unsupported-contact",
        "{}",
        refusal.detail
    );
    assert!(refusal.detail.contains("consumes"), "{}", refusal.detail);
}
