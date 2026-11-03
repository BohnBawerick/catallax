//! Straight and arc strip ends and corners: strips trimmed by their own end
//! faces in pockets and on bosses, concave sphere corners, two strips of
//! unequal reach, arc strips ending at a meridian and straight strips ending
//! on a sphere, each against an exact or independent volume, the field where
//! it applies, a closed display and STEP; and typed refusals for the rest.

use std::collections::{BTreeMap, BTreeSet};

use catallax_boolean::Containment;
use catallax_features::{BuildError, Built, evaluate, evaluate_incremental};
use catallax_field::{Location, Plan, Workspace};
use catallax_mesh::MeshOptions;
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational, Real};
use catallax_recipe::{Origin, Recipe, Role};
use catallax_shapes::{Curve, Point3, RPoint3, Surface};
use catallax_topology::{FaceId, Sense};
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
fn edge(a: Value, b: Value) -> Value {
    json!({"between":[a,b]})
}
fn recipe(
    format: u8,
    mut steps: Vec<Value>,
    body: &str,
    radius: &str,
    edges: Vec<Value>,
) -> Recipe {
    steps.push(json!({"id":"round","op":"fillet","body":body,"radius":radius,"edges":edges}));
    Recipe::parse(
        &json!({"format":format!("catallax-recipe/{format}"),"name":"strip ends",
            "steps":steps,"result":"round"})
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
/// The blend and corner origins, which a radius edit must keep.
fn names(built: &Built) -> BTreeSet<Origin> {
    built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .flat_map(|face| built.model.face_origins(face).to_vec())
        .filter(|o| matches!(o.role, Role::Blend { .. } | Role::Corner { .. }))
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
/// The measured volume lies within `tolerance` of an independent oracle.
fn oracle_volume(built: &Built, want: &str, tolerance: &str) {
    let volume = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut measuring(),
    )
    .unwrap()
    .value;
    let (want, tolerance) = (r(want), r(tolerance));
    assert!(
        &want - &tolerance <= volume.lo && volume.hi <= &want + &tolerance,
        "[{}, {}] not within {tolerance} of {want}",
        volume.lo.to_f64_nearest(),
        volume.hi.to_f64_nearest()
    );
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
/// Both evaluators agree on every point of a lattice that both decide, and
/// both decide most of them. The lattice starts off any support plane.
fn agree(built: &mut Built, recipe: &Recipe, low: [&str; 3], step: &str, count: usize) {
    let plan = Plan::compile_checked(recipe, &mut decider()).unwrap();
    let blend: BTreeSet<_> = plan
        .faces()
        .iter()
        .filter(|f| matches!(f.origin.role, Role::Blend { .. } | Role::Corner { .. }))
        .map(|f| f.origin.clone())
        .collect();
    assert_eq!(blend, names(built));
    let mut workspace = Workspace::default();
    let (low, step) = (low.map(r), r(step));
    let mut decided = 0;
    for i in 0..count {
        for j in 0..count {
            for k in 0..count {
                let at =
                    |base: &Rational, n: usize| base + &(&step * &Rational::from_i64(n as i64));
                let p = RPoint3::new(at(&low[0], i), at(&low[1], j), at(&low[2], k));
                let field = plan
                    .proven_checked(&p, &mut workspace, &mut decider())
                    .unwrap()
                    .location;
                let brep = catallax_boolean::classify_point(
                    &mut built.model,
                    built.result,
                    &p,
                    &mut decider(),
                );
                match (field, &brep) {
                    (Location::Inside, Outcome::Proven(Containment::Inside))
                    | (Location::Outside, Outcome::Proven(Containment::Outside)) => decided += 1,
                    (Location::Inside | Location::Outside, Outcome::Proven(_)) => {
                        panic!("at {p:?}: field {field:?}, B-rep {brep:?}")
                    }
                    _ => {}
                }
            }
        }
    }
    assert!(decided * 5 >= count * count * count * 4, "{decided}");
}
/// SolidDraft's default display, a 15 degree normal limit: every face
/// meshed, every triangle edge used twice in opposite directions.
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
fn faces_with(built: &Built, role: fn(&Role) -> bool) -> Vec<FaceId> {
    built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|f| built.model.face_origins(*f).iter().any(|o| role(&o.role)))
        .collect()
}
fn exact(p: &Point3) -> [Rational; 3] {
    let p = p.as_rational().expect("rational point");
    [p.x, p.y, p.z]
}
/// The body's vertex exactly at a point.
fn vertex_at(built: &Built, at: &[Rational; 3]) -> Option<catallax_topology::VertexId> {
    let target = RPoint3::new(at[0].clone(), at[1].clone(), at[2].clone()).to_real();
    built
        .model
        .vertices_of_body(built.result)
        .into_iter()
        .find(|v| {
            let p = built.model.geom().position(built.model.vertex(*v).point);
            matches!(
                decider().sign(&p.minus(&target).norm_squared()),
                Outcome::ProvenCoincident(_)
            )
        })
}

// A 20 x 20 x 10 block with a 10 x 10 pocket 6 deep, from z = 4.
fn pocket() -> Vec<Value> {
    vec![
        json!({"id":"block","op":"block","size":["20","20","10"],"pose":pose(["0","0","0"])}),
        json!({"id":"tool","op":"block","size":["10","10","10"],"pose":pose(["5","5","4"])}),
        json!({"id":"pocket","op":"boolean","kind":"subtract","target":"block","tool":"tool"}),
    ]
}
fn floor(wall: &str) -> Value {
    edge(cap("tool", "start"), side("tool", wall))
}
fn upright(a: &str, b: &str) -> Value {
    edge(side("tool", a), side("tool", b))
}
fn rim(wall: &str) -> Value {
    edge(cap("block", "end"), side("tool", wall))
}
const WALLS: [&str; 4] = ["x_min", "y_min", "x_max", "y_max"];

#[test]
fn a_pocket_floor_corner_closes_with_a_concave_sphere_sector() {
    let corner = |radius: &str| {
        recipe(
            3,
            pocket(),
            "pocket",
            radius,
            vec![floor("x_min"), floor("y_min"), upright("x_min", "y_min")],
        )
    };
    for radius in ["1", "3/2"] {
        let built = build(&corner(radius));
        expect_kinds(&built, &[("Plane", 11), ("Cylinder", 3), ("Sphere", 1)]);
        let rho = r(radius);
        let [sphere] = faces_with(&built, |role| matches!(role, Role::Corner { .. }))[..] else {
            panic!("one corner face");
        };
        let face = built.model.face(sphere);
        // The material lies outside the ball, so the face turns toward its centre.
        assert_eq!(face.sense, Sense::Reversed);
        let Surface::Sphere(s) = built.model.geom().surface(face.surface) else {
            panic!("the corner is a sphere");
        };
        let at = |base: &str| &r(base) + &rho;
        assert_eq!(exact(&s.center), [at("5"), at("5"), at("4")]);
        assert_eq!(s.radius.as_rational(), Some(rho.clone()));
        // 3400 plus each strip's section r^2 (1 - pi/4) over its run beyond the
        // corner cube, and that cube less its ball octant.
        let (r2, r3) = (&rho * &rho, &(&rho * &rho) * &rho);
        let runs = &r("26") - &(&r("3") * &rho);
        exact_volume(
            &built,
            [
                &(&r("3400") + &(&r2 * &runs)) + &r3,
                -&(&(&(&r2 * &runs) / &r("4")) + &(&r3 / &r("6"))),
                r("0"),
            ],
            "1e-6",
        );
        closed_display(&built);
        step(&built, "SPHERICAL_SURFACE(");
    }
    assert_eq!(names(&build(&corner("1"))), names(&build(&corner("3/2"))));
    let recipe = corner("1");
    let mut built = build(&recipe);
    use Containment::{Inside, Outside};
    classify(
        &mut built,
        &[
            // Added: the corner cube outside the ball, and a floor strip.
            (["26/5", "26/5", "21/5"], Inside),
            (["26/5", "10", "21/5"], Inside),
            // The ball and the pocket stay empty.
            (["29/5", "29/5", "24/5"], Outside),
            (["11/2", "10", "9/2"], Outside),
            (["10", "10", "7"], Outside),
            // Uncut material.
            (["2", "2", "2"], Inside),
        ],
    );
    agree(&mut built, &recipe, ["36/7", "36/7", "29/7"], "1/4", 9);
}

#[test]
fn a_whole_pocket_rounds_every_inner_edge() {
    let edges = || {
        WALLS
            .iter()
            .map(|w| floor(w))
            .chain((0..4).map(|i| upright(WALLS[i], WALLS[(i + 1) % 4])))
            .collect::<Vec<_>>()
    };
    let rounded = |radius: &str| recipe(3, pocket(), "pocket", radius, edges());
    let built = build(&rounded("1"));
    expect_kinds(&built, &[("Plane", 11), ("Cylinder", 8), ("Sphere", 4)]);
    // Floor runs 10 - 2r, uprights 6 - r, and four corner cubes less octants.
    exact_volume(&built, [r("3456"), -&r("41/3"), r("0")], "1e-6");
    assert_eq!(names(&built), names(&build(&rounded("3/2"))));
    closed_display(&built);
    step(&built, "SPHERICAL_SURFACE(");
    let recipe = rounded("1");
    let mut built = built;
    agree(&mut built, &recipe, ["29/7", "29/7", "22/7"], "1/2", 6);
}

#[test]
fn pocket_rim_and_floor_loops_mitre_at_their_corners() {
    // Each loop's section r^2 (1 - pi/4), its centroid r (5/6 - pi/4) / (1 - pi/4)
    // from both faces, swept round a square of side 10 + 2u (rim) or 10 - 2u.
    for radius in ["1", "2"] {
        let rho = r(radius);
        let (r2, r3) = (&rho * &rho, &(&rho * &rho) * &rho);
        let section = [&r("40") * &r2, -&(&r("10") * &r2)];
        let moment = [&r("40/6") * &r3, -&(&r("2") * &r3)];
        let rim_loop = build(&recipe(
            3,
            pocket(),
            "pocket",
            radius,
            WALLS.map(rim).to_vec(),
        ));
        expect_kinds(&rim_loop, &[("Plane", 11), ("Cylinder", 4)]);
        exact_volume(
            &rim_loop,
            [
                &(&r("3400") - &section[0]) - &moment[0],
                &(-&section[1]) - &moment[1],
                r("0"),
            ],
            "1e-6",
        );
        closed_display(&rim_loop);
        let floor_loop = build(&recipe(
            3,
            pocket(),
            "pocket",
            radius,
            WALLS.map(floor).to_vec(),
        ));
        exact_volume(
            &floor_loop,
            [
                &(&r("3400") + &section[0]) - &moment[0],
                &section[1] - &moment[1],
                r("0"),
            ],
            "1e-6",
        );
        closed_display(&floor_loop);
    }
}

// A prism on (0,0) (12,0) (12,5) (0,10), 10 high: the slope meets the right
// wall at 112.6 degrees and the top at 90, so their strips reach the right
// wall's top edge at different points.
fn slope_prism() -> Vec<Value> {
    vec![
        json!({"id":"profile","op":"sketch","plane":pose(["0","0","0"]),
        "loops":[[
            {"kind":"line","id":"bottom","to":["12","0"]},
            {"kind":"line","id":"right","to":["12","5"]},
            {"kind":"line","id":"slope","to":["0","10"]},
            {"kind":"line","id":"left","to":["0","0"]}
        ]]}),
        json!({"id":"body","op":"extrude","profile":"profile","distance":"10"}),
    ]
}

#[test]
fn strips_of_unequal_reach_meet_on_an_exact_step() {
    // Independent OpenCASCADE booleans of the section 4.5 definition: the
    // prism less the two kite prisms, each cut by the other's outer face.
    let oracle = [("1", "896.4677608418"), ("2", "886.0479447309")];
    let mut previous = None;
    for (radius, volume) in oracle {
        for order in [0, 1] {
            let mut edges = vec![
                edge(side("body", "slope"), side("body", "right")),
                edge(side("body", "slope"), cap("body", "end")),
            ];
            if order == 1 {
                edges.reverse();
            }
            let recipe = recipe(3, slope_prism(), "body", radius, edges);
            let mut built = build(&recipe);
            expect_kinds(&built, &[("Plane", 6), ("Cylinder", 2)]);
            oracle_volume(&built, volume, "1e-7");
            // The shared rail point, the right rail's entry into the top
            // strip and the top strip's end on the right wall's top edge.
            let rho = r(radius);
            let p = [
                &r("12") - &(&r("8/13") * &rho),
                &r("5") + &(&r("10/39") * &rho),
                &r("10") - &rho,
            ];
            let m = [
                r("12"),
                &r("5") - &(&r("2/3") * &rho),
                &r("10") - &(&rho / &r("13")),
            ];
            let q = [r("12"), &r("5") - &(&r("13/12") * &rho), r("10")];
            assert!(vertex_at(&built, &p).is_some(), "{p:?}");
            let m = vertex_at(&built, &m).expect("the entry point");
            let q = vertex_at(&built, &q).expect("the step's foot");
            // The top strip is cut by the right wall from m to q: an ellipse
            // shared by those two faces.
            let model = &built.model;
            let step_edge = model
                .edges_of_body(built.result)
                .into_iter()
                .find(|e| {
                    let e = model.edge(*e);
                    [e.start, e.end] == [Some(m), Some(q)] || [e.start, e.end] == [Some(q), Some(m)]
                })
                .expect("an edge from m to q");
            assert!(matches!(
                model.geom().curve(model.edge(step_edge).curve),
                Curve::Ellipse(_)
            ));
            let owners: BTreeSet<_> = model
                .coedges_of_edge(step_edge)
                .iter()
                .map(|c| {
                    let face = model.owner_of_coedge(*c).unwrap().1;
                    format!(
                        "{:?}",
                        model.geom().surface(model.face(face).surface).kind()
                    )
                })
                .collect();
            assert_eq!(owners, ["Cylinder".to_string(), "Plane".to_string()].into());
            closed_display(&built);
            step(&built, "CYLINDRICAL_SURFACE(");
            let named = names(&built);
            assert_eq!(previous.get_or_insert_with(|| named.clone()), &named);
            if radius == "1" && order == 0 {
                use Containment::{Inside, Outside};
                classify(
                    &mut built,
                    &[
                        // Past the vertical strip's rail the top strip still
                        // cuts the right wall's top: the step.
                        (["1199/100", "21/5", "1999/200"], Outside),
                        (["599/50", "4", "999/100"], Inside),
                        (["11", "4", "9"], Inside),
                        (["119/10", "49/10", "99/10"], Outside),
                    ],
                );
                agree(&mut built, &recipe, ["71/7", "22/7", "57/7"], "1/4", 9);
            }
        }
    }
}

// The slope prism cut 6 deep into a block from z = 4: its floor meets the
// vertical right wall at 90 degrees and the sloped wall elsewhere.
fn drafted_pocket() -> Vec<Value> {
    vec![
        json!({"id":"block","op":"block","size":["30","30","10"],"pose":pose(["-5","-5","0"])}),
        json!({"id":"profile","op":"sketch","plane":pose(["0","0","4"]),
        "loops":[[
            {"kind":"line","id":"bottom","to":["12","0"]},
            {"kind":"line","id":"right","to":["12","5"]},
            {"kind":"line","id":"slope","to":["0","10"]},
            {"kind":"line","id":"left","to":["0","0"]}
        ]]}),
        json!({"id":"tool","op":"extrude","profile":"profile","distance":"10"}),
        json!({"id":"pocket","op":"boolean","kind":"subtract","target":"block","tool":"tool"}),
    ]
}

#[test]
fn concave_strips_of_unequal_reach_meet_on_an_exact_step() {
    // Independent OpenCASCADE booleans: the block plus both concave kite
    // prisms, each kept inside the pocket.
    for (radius, volume) in [("1", "8463.7789465127"), ("2", "8474.7802398620")] {
        let edges = vec![floor("right"), floor("slope")];
        let recipe = recipe(3, drafted_pocket(), "pocket", radius, edges);
        let mut built = build(&recipe);
        expect_kinds(&built, &[("Plane", 11), ("Cylinder", 2)]);
        oracle_volume(&built, volume, "1e-7");
        closed_display(&built);
        step(&built, "CYLINDRICAL_SURFACE(");
        if radius == "1" {
            agree(&mut built, &recipe, ["71/7", "22/7", "29/7"], "1/4", 9);
        }
    }
}

#[test]
fn every_edge_of_a_sloped_prism_rounds_with_disjoint_strip_hulls() {
    // The parallel body of the inward offset: (1,1) (11,1) (11,13/3) (1,17/2),
    // 8 high, area 325/6 and perimeter 95/3, so 795 + 151 pi / 6 by Steiner.
    let recipe = Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"slope","steps":[
            slope_prism()[0].clone(), slope_prism()[1].clone(),
            {"id":"round","op":"fillet","body":"body","radius":"1","edges":"all"}
        ],"result":"round"})
        .to_string(),
    )
    .unwrap();
    let mut built = build(&recipe);
    expect_kinds(&built, &[("Plane", 6), ("Cylinder", 12), ("Sphere", 8)]);
    exact_volume(&built, [r("795"), r("151/6"), r("0")], "1e-6");
    closed_display(&built);
    step(&built, "SPHERICAL_SURFACE(");
    agree(&mut built, &recipe, ["1/7", "1/7", "1/7"], "3/2", 8);
}

// A half disc of radius 5 about the origin, flat side on x = 0.
fn half_disc(at: &str) -> Value {
    json!({"id":"profile","op":"sketch","plane":pose(["0","0",at]),
    "loops":[[
        {"kind":"arc","id":"round","centre":["0","0"],"to":["0","5"],"turn":"ccw"},
        {"kind":"line","id":"flat","to":["0","-5"]}
    ]]})
}
fn d_boss() -> Vec<Value> {
    vec![
        half_disc("0"),
        json!({"id":"body","op":"extrude","profile":"profile","distance":"10"}),
    ]
}
fn d_pocket() -> Vec<Value> {
    vec![
        json!({"id":"block","op":"block","size":["20","20","10"],"pose":pose(["-10","-10","0"])}),
        half_disc("4"),
        json!({"id":"tool","op":"extrude","profile":"profile","distance":"10"}),
        json!({"id":"pocket","op":"boolean","kind":"subtract","target":"block","tool":"tool"}),
    ]
}
/// Half a turn of the ring section's first moment about the axis, the
/// section r^2 (1 - pi/4) at R -+ its centroid offset r (5/6 - pi/4) / (1 - pi/4):
/// pi (R r^2 (1 - pi/4) -+ r^3 (5/6 - pi/4)), as [pi, pi^2] coefficients.
fn half_ring(rho: &Rational, outside: bool) -> [Rational; 2] {
    let (r2, r3) = (rho * rho, &(rho * rho) * rho);
    let s = r(if outside { "1" } else { "-1" });
    [
        &(&r("5") * &r2) + &(&s * &(&r("5/6") * &r3)),
        -&(&(&r("5/4") * &r2) + &(&s * &(&r3 / &r("4")))),
    ]
}

#[test]
fn arc_strips_end_at_meridians_on_planes_through_their_axis() {
    let mut previous = None;
    for format in [3, 4] {
        for radius in ["1", "3/2"] {
            let rho = r(radius);
            let built = build(&recipe(
                format,
                d_boss(),
                "body",
                radius,
                vec![edge(cap("body", "end"), side("body", "round"))],
            ));
            expect_kinds(&built, &[("Plane", 3), ("Cylinder", 1), ("Torus", 1)]);
            // 125 pi less half a turn of the inside section.
            let [a, b] = half_ring(&rho, false);
            exact_volume(&built, [r("0"), &r("125") - &a, -&b], "1e-6");
            // Both ends are meridian quarter circles in x = 0 about the ball
            // centres (0, +-(5 - r), 10 - r).
            let model = &built.model;
            let [strip] = faces_with(&built, |role| matches!(role, Role::Blend { .. }))[..] else {
                panic!("one strip");
            };
            let mut meridians = 0;
            for lp in &model.face(strip).loops {
                let catallax_topology::Loop::Coedges(cs) = model.loop_(*lp) else {
                    continue;
                };
                for c in cs {
                    let Curve::Circle(circle) =
                        model.geom().curve(model.edge(model.coedge(*c).edge).curve)
                    else {
                        continue;
                    };
                    let zero = |v: Real| matches!(decider().sign(&v), Outcome::ProvenCoincident(_));
                    if zero(circle.normal.y.clone()) && zero(circle.normal.z.clone()) {
                        let (height, offset) = (&r("10") - &rho, &r("5") - &rho);
                        assert!(zero(circle.center.x.clone()));
                        assert!(zero(&circle.center.z - &Real::from(height)));
                        assert!(zero(
                            circle.center.y.square() - Real::from(&offset * &offset)
                        ));
                        assert!(zero(&circle.radius_sq - &Real::from(&rho * &rho)));
                        meridians += 1;
                    }
                }
            }
            assert_eq!(meridians, 2);
            closed_display(&built);
            step(&built, "TOROIDAL_SURFACE(");
            let named = names(&built);
            assert_eq!(previous.get_or_insert_with(|| named.clone()), &named);
        }
    }
    let mut built = build(&recipe(
        3,
        d_boss(),
        "body",
        "1",
        vec![edge(cap("body", "end"), side("body", "round"))],
    ));
    use Containment::{Inside, Outside};
    classify(
        &mut built,
        &[
            (["1/10", "49/10", "99/10"], Outside),
            (["49/10", "0", "99/10"], Outside),
            (["-1/10", "0", "5"], Outside),
            (["1/10", "47/10", "93/10"], Inside),
            (["5/2", "0", "99/10"], Inside),
        ],
    );
}

#[test]
fn arc_strips_in_a_half_disc_pocket_end_at_meridians() {
    // 4000 less the 75 pi pocket, plus the floor ring's inside section or
    // less the rim ring's outside one, each over half a turn.
    let [a, b] = half_ring(&r("1"), false);
    let floor = build(&recipe(3, d_pocket(), "pocket", "1", vec![floor("round")]));
    expect_kinds(&floor, &[("Plane", 8), ("Cylinder", 1), ("Torus", 1)]);
    exact_volume(&floor, [r("4000"), &a - &r("75"), b], "1e-6");
    closed_display(&floor);
    step(&floor, "TOROIDAL_SURFACE(");
    let [a, b] = half_ring(&r("1"), true);
    let rim = build(&recipe(3, d_pocket(), "pocket", "1", vec![rim("round")]));
    exact_volume(&rim, [r("4000"), -&(&a + &r("75")), -&b], "1e-6");
    closed_display(&rim);
    step(&rim, "TOROIDAL_SURFACE(");
}

#[test]
fn ring_rounds_point_their_torus_from_the_ball_to_the_plane_contact() {
    // A torus is the same either way round, but STEP charts a face that
    // meets its outer equator only on the axis side: a cylinder's bottom
    // edge and a blind pocket's floor used to write unreadable faces.
    let boss = vec![
        json!({"id":"boss","op":"cylinder","radius":"5","height":"10","pose":pose(["0","0","0"])}),
    ];
    let blind = vec![
        json!({"id":"block","op":"block","size":["20","20","10"],"pose":pose(["-10","-10","0"])}),
        json!({"id":"tool","op":"cylinder","radius":"5","height":"10","pose":pose(["0","0","4"])}),
        json!({"id":"pocket","op":"boolean","kind":"subtract","target":"block","tool":"tool"}),
    ];
    let [a, b] = half_ring(&r("1"), false);
    let full = |c: &Rational| c * &r("2");
    for (steps, body, edges, toward, volume) in [
        (
            boss.clone(),
            "boss",
            vec![edge(cap("boss", "start"), side("boss", "wall"))],
            "-1",
            [r("0"), &r("250") - &full(&a), -&full(&b)],
        ),
        (
            boss,
            "boss",
            vec![edge(cap("boss", "end"), side("boss", "wall"))],
            "1",
            [r("0"), &r("250") - &full(&a), -&full(&b)],
        ),
        (
            blind,
            "pocket",
            vec![edge(cap("tool", "start"), side("tool", "wall"))],
            "-1",
            [r("4000"), &full(&a) - &r("150"), full(&b)],
        ),
    ] {
        let built = build(&recipe(3, steps, body, "1", edges));
        let [strip] = faces_with(&built, |role| matches!(role, Role::Blend { .. }))[..] else {
            panic!("one ring");
        };
        let Surface::Torus(torus) = built.model.geom().surface(built.model.face(strip).surface)
        else {
            panic!("the ring is a torus");
        };
        assert_eq!(
            torus.axis.z.as_rational().map(|z| z.sign()),
            Some(r(toward).sign())
        );
        exact_volume(&built, volume, "1e-6");
        closed_display(&built);
        step(&built, "TOROIDAL_SURFACE(");
    }
}

#[test]
fn a_straight_strip_ends_on_a_spherical_wall() {
    // A 10 x 10 block under a ball of radius 12 about (0, 0, -2). The
    // vertical edge at (5, 5) ends on the sphere; its strip runs up to it.
    let steps = vec![
        json!({"id":"block","op":"block","size":["10","10","10"],"pose":pose(["-5","-5","0"])}),
        json!({"id":"ball","op":"sphere","radius":"12","pose":pose(["0","0","-2"])}),
        json!({"id":"cap","op":"boolean","kind":"intersect","target":"block","tool":"ball"}),
    ];
    let rounded = |radius| {
        recipe(
            3,
            steps.clone(),
            "cap",
            radius,
            vec![edge(side("block", "x_max"), side("block", "y_max"))],
        )
    };
    let mut built = build(&rounded("1"));
    expect_kinds(&built, &[("Plane", 5), ("Cylinder", 1), ("Sphere", 1)]);
    // Quadrature of the sphere's height over the square less the corner kite.
    oracle_volume(&built, "925.731827569374120214536685369", "1e-9");
    let model = &built.model;
    let [strip] = faces_with(&built, |role| matches!(role, Role::Blend { .. }))[..] else {
        panic!("one strip");
    };
    let curved_ends = model
        .face(strip)
        .loops
        .iter()
        .flat_map(|lp| match model.loop_(*lp) {
            catallax_topology::Loop::Coedges(cs) => cs.clone(),
            _ => Vec::new(),
        })
        .filter(|c| {
            matches!(
                model.geom().curve(model.edge(model.coedge(*c).edge).curve),
                Curve::Intersection(_)
            )
        })
        .count();
    assert_eq!(curved_ends, 1);
    closed_display(&built);
    step(&built, "SPHERICAL_SURFACE(");
    assert_eq!(names(&built), names(&build(&rounded("3/2"))));
    use Containment::{Inside, Outside};
    classify(
        &mut built,
        &[
            (["49/10", "49/10", "7"], Outside),
            (["9/2", "9/2", "7"], Inside),
            (["0", "0", "9"], Inside),
            (["49/10", "49/10", "1/10"], Outside),
        ],
    );
}

fn refused(recipe: &Recipe) -> (String, String) {
    match evaluate(recipe, &mut decider()) {
        Err(BuildError::Geometry { step, refusal }) => {
            assert_eq!(step.as_str(), "round");
            (refusal.code.to_string(), refusal.detail)
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn unfinished_ends_and_corners_refuse_with_their_own_codes() {
    const CORNER: &str = catallax_features::fillet::UNSUPPORTED_CORNER;
    const END: &str = catallax_features::fillet::UNSUPPORTED_END;
    // A step on the y = 0 face: low top, riser and front meet at (10, 0, 5),
    // where the low top's front edge is convex and its riser edge concave.
    let stepped = vec![
        json!({"id":"profile","op":"sketch","plane":{"rotation":["1","1","0","0"],"translation":["0","0","0"]},
        "loops":[[
            {"kind":"line","id":"floor","to":["20","0"]},
            {"kind":"line","id":"back","to":["20","10"]},
            {"kind":"line","id":"high","to":["10","10"]},
            {"kind":"line","id":"riser","to":["10","5"]},
            {"kind":"line","id":"low","to":["0","5"]},
            {"kind":"line","id":"front","to":["0","0"]}
        ]]}),
        json!({"id":"body","op":"extrude","profile":"profile","distance":"10"}),
    ];
    // Two arcs crossing at a sharp point.
    let lens = vec![
        json!({"id":"profile","op":"sketch","plane":pose(["0","0","0"]),
        "loops":[[
            {"kind":"arc","id":"east","centre":["-3","0"],"to":["0","4"],"turn":"ccw"},
            {"kind":"arc","id":"west","centre":["3","0"],"to":["0","-4"],"turn":"ccw"}
        ]]}),
        json!({"id":"body","op":"extrude","profile":"profile","distance":"10"}),
    ];
    // A quarter disc closed by a chord that misses the arc's axis.
    let chord = vec![
        json!({"id":"profile","op":"sketch","plane":pose(["0","0","0"]),
        "loops":[[
            {"kind":"arc","id":"round","centre":["0","0"],"to":["0","5"],"turn":"ccw"},
            {"kind":"line","id":"flat","to":["5","0"]}
        ]]}),
        json!({"id":"body","op":"extrude","profile":"profile","distance":"10"}),
    ];
    // A pocket whose right wall leans out by 5 in 12, cut across y.
    let leaning = vec![
        json!({"id":"block","op":"block","size":["30","30","10"],"pose":pose(["-5","-15","0"])}),
        json!({"id":"profile","op":"sketch","plane":{"rotation":["1","1","0","0"],"translation":["0","0","0"]},
        "loops":[[
            {"kind":"line","id":"bottom","to":["12","4"]},
            {"kind":"line","id":"lean","to":["17","16"]},
            {"kind":"line","id":"top","to":["0","16"]},
            {"kind":"line","id":"left","to":["0","4"]}
        ]]}),
        json!({"id":"tool","op":"extrude","profile":"profile","distance":"10"}),
        json!({"id":"pocket","op":"boolean","kind":"subtract","target":"block","tool":"tool"}),
    ];
    let base = || {
        evaluate_incremental(
            &Recipe::parse(
                &json!({"format":"catallax-recipe/3","name":"pocket",
        "steps":pocket(),"result":"pocket"})
                .to_string(),
            )
            .unwrap(),
            None,
            &mut decider(),
        )
        .unwrap()
    };
    let cases = [
        // Convex rim strips over a concave pocket edge.
        (
            recipe(
                3,
                pocket(),
                "pocket",
                "1",
                vec![rim("x_min"), rim("y_min"), upright("x_min", "y_min")],
            ),
            CORNER,
            "convex and concave",
        ),
        (
            recipe(
                3,
                stepped.clone(),
                "body",
                "1",
                vec![
                    edge(side("body", "low"), cap("body", "start")),
                    edge(side("body", "low"), side("body", "riser")),
                ],
            ),
            CORNER,
            "convex and a concave",
        ),
        (
            recipe(
                3,
                d_boss(),
                "body",
                "1",
                vec![
                    edge(cap("body", "end"), side("body", "round")),
                    edge(cap("body", "end"), side("body", "flat")),
                ],
            ),
            CORNER,
            "straight rails",
        ),
        (
            recipe(
                3,
                d_boss(),
                "body",
                "1",
                vec![
                    edge(cap("body", "end"), side("body", "round")),
                    edge(cap("body", "end"), side("body", "flat")),
                    edge(side("body", "round"), side("body", "flat")),
                ],
            ),
            CORNER,
            "three planar supports",
        ),
        (
            recipe(
                3,
                lens.clone(),
                "body",
                "1",
                vec![
                    edge(cap("body", "end"), side("body", "east")),
                    edge(cap("body", "end"), side("body", "west")),
                ],
            ),
            CORNER,
            "straight rails",
        ),
        (
            recipe(
                3,
                leaning,
                "pocket",
                "1",
                vec![
                    edge(cap("block", "end"), side("tool", "lean")),
                    edge(cap("block", "end"), cap("tool", "start")),
                ],
            ),
            CORNER,
            "other convexity",
        ),
        (
            recipe(
                3,
                lens,
                "body",
                "1",
                vec![edge(cap("body", "end"), side("body", "east"))],
            ),
            END,
            "planar end face",
        ),
        (
            recipe(
                3,
                chord,
                "body",
                "1",
                vec![edge(cap("body", "end"), side("body", "round"))],
            ),
            END,
            "through its ring axis",
        ),
    ];
    for (recipe, code, detail) in cases {
        let (found, text) = refused(&recipe);
        assert_eq!(
            (found.as_str(), text.contains(detail)),
            (code, true),
            "{text}"
        );
    }
    // A refusal leaves the retained input as it was.
    let retained = base();
    let bytes = {
        let mut bytes = Vec::new();
        catallax_step::write_step(
            &retained.built().model,
            &[retained.built().result],
            &Default::default(),
            &mut Decider::default(),
            &mut bytes,
        )
        .unwrap();
        bytes
    };
    let mixed = recipe(
        3,
        pocket(),
        "pocket",
        "1",
        vec![rim("x_min"), rim("y_min"), upright("x_min", "y_min")],
    );
    let Err(BuildError::Geometry { refusal, .. }) =
        evaluate_incremental(&mixed, Some(retained.clone()), &mut decider())
    else {
        panic!("the mixed corner built");
    };
    assert_eq!(refusal.code, CORNER);
    let mut again = Vec::new();
    catallax_step::write_step(
        &retained.built().model,
        &[retained.built().result],
        &Default::default(),
        &mut Decider::default(),
        &mut again,
    )
    .unwrap();
    assert_eq!(again, bytes);
    assert_eq!(
        retained.built().model.check(retained.built().result),
        vec![]
    );
}
