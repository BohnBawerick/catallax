//! 4.2a step 5: fillets that roll over an earlier straight round. The
//! corpus rows `tangent-boss`, `tangent-pocket` and `unequal-ordered`, and
//! the recipes of `fillet_over_fillet.rs` whose chain rolls over an earlier
//! strip with its ball on the strip's axis side (part 1), or with a ball as
//! wide as the strip, a sphere corner, and the recipes whose straight run
//! stops at a strip narrower than its ball and takes the cliff (part 2), and
//! the recipes whose straight run swallows an earlier mitre of two rounds as
//! wide as its ball, the corpus row `roll-over-mitres` among them, or that
//! takes the cliff at both feet under a mitre of narrower rounds, the corpus
//! row `unequal-mitres` (part 3).
//! The corpus row `box-mitre`, which the elementary route already builds, is
//! checked for membership.
//!
//! The field compiles each chain from its own path pieces: straight pieces
//! along the sharp runs and, over each earlier round, a ring sector about the
//! round's centre line (a ring, horn or spindle by the radii, or one ball on
//! the line when the radii are equal), and, where a run swallows an earlier
//! mitre, the trihedral corner of its three planes less the ball, whose
//! sphere face is registered under the corner's name. Its membership
//! must match each witness: the corpus row's (confirmed by OpenCASCADE), the
//! test file's (checked against the B-rep there) and points derived here,
//! each checked against the B-rep first. Every display vertex of the B-rep's
//! blend faces of the rolling step must lie within the field's guaranteed
//! distance range of a piece's face at 1 nm, with the same records at 1, 2
//! and 4 workers; every display vertex of a mitre corner's sphere lies within
//! 1 nm of the field's corner face. The step 5 recipes left for later parts
//! keep their refusals, each named.
use catallax_features::evaluate;
use catallax_field::rolling_region::piece::Piece;
use catallax_field::rolling_region::{Batch, Range, Relation, distance_batch};
use catallax_field::{Direction, FaceState, Location, Plan, Workspace};
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::{BudgetPolicy, Decider, OperationAllowance, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{Point3, RPoint3, Surface};
use serde_json::{Value, json};
use std::collections::BTreeMap;

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

// The recipes of `fillet_over_fillet.rs`.
fn side(step: &str, segment: &str) -> Value {
    json!({"step":step,"role":{"side":{"segment":segment}}})
}
fn cap(step: &str, role: &str) -> Value {
    json!({"step":step,"role":role})
}
fn made(step: &str) -> Value {
    json!({"step":step})
}
fn edge(a: Value, b: Value) -> Value {
    json!({"between":[a,b]})
}
fn fillet(id: &str, body: &str, radius: &str, edges: Vec<Value>) -> Value {
    json!({"id":id,"op":"fillet","body":body,"radius":radius,"edges":edges})
}
fn block(id: &str, size: [&str; 3], at: [&str; 3]) -> Value {
    json!({"id":id,"op":"block","size":size,
        "pose":{"rotation":["1","0","0","0"],"translation":at}})
}
fn recipe(format: u8, steps: Vec<Value>) -> Recipe {
    let result = steps.last().unwrap()["id"].clone();
    parse(
        &json!({"format":format!("catallax-recipe/{format}"),"name":"fillet over fillet",
        "steps":steps,"result":result}),
    )
}
const SIDES: [&str; 4] = ["x_min", "y_min", "x_max", "y_max"];
fn uprights(step: &str) -> Vec<Value> {
    (0..4)
        .map(|i| edge(side(step, SIDES[i]), side(step, SIDES[(i + 1) % 4])))
        .collect()
}
/// A 20 x 20 x 10 block with its uprights rounded, then one top edge, which
/// version 4 continues round the loop; or, in version 3, the whole loop
/// named.
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
/// The block with a 10 x 10 pocket 6 deep from z = 4, its uprights
/// rounded, then its floor edge along the pocket's x_max wall.
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
/// The block with its top front edge rounded, then a right-hand edge: the
/// top one, or the round's own end arc.
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
/// The blocks in a rational pose: a third of a turn about (1, 1, 1), then
/// moved by (7, -2, 5).
fn placed(recipe: &Recipe) -> Recipe {
    let mut value = serde_json::to_value(recipe).unwrap();
    for step in value["steps"].as_array_mut().unwrap() {
        if step["op"] == "block" {
            step["pose"] = json!({"rotation":["3","1","1","1"],"translation":["7","-2","5"]});
        }
    }
    parse(&value)
}
/// A point moved as `placed` moves the blocks: the quaternion (3, 1, 1, 1)
/// turns by (2 -1 2; 2 2 -1; -1 2 2) / 3.
fn place(p: &[Rational; 3]) -> [Rational; 3] {
    const M: [[i64; 3]; 3] = [[2, -1, 2], [2, 2, -1], [-1, 2, 2]];
    const T: [i64; 3] = [7, -2, 5];
    std::array::from_fn(|i| {
        let turned = (0..3).fold(Rational::zero(), |sum, j| {
            sum + Rational::from_i64(M[i][j]) * p[j].clone()
        });
        turned / Rational::from_i64(3) + Rational::from_i64(T[i])
    })
}

fn corpus_row(id: &str) -> Value {
    include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|row| row["id"] == id)
        .unwrap()
}

type Points = [([i64; 2], [i64; 2], [i64; 2])];
fn points(list: &Points) -> Vec<[Rational; 3]> {
    list.iter()
        .map(|(x, y, z)| [ratio(x[0], x[1]), ratio(y[0], y[1]), ratio(z[0], z[1])])
        .collect()
}

/// One recipe whose last step rolls over an earlier round.
struct Case {
    name: String,
    recipe: Recipe,
    /// The rolling step, whose B-rep blend faces are checked.
    step: &'static str,
    /// Rolling pieces the field compiles: one per straight run and ring
    /// sector of the chain.
    pieces: usize,
    inside: Vec<[Rational; 3]>,
    outside: Vec<[Rational; 3]>,
    /// Whether the points are derived here, so the test first checks each
    /// against the B-rep.
    derived: bool,
}
impl Case {
    fn new(
        name: &str,
        recipe: Recipe,
        step: &'static str,
        pieces: usize,
        inside: &Points,
        outside: &Points,
    ) -> Self {
        Self {
            name: name.to_string(),
            recipe,
            step,
            pieces,
            inside: points(inside),
            outside: points(outside),
            derived: true,
        }
    }
    /// The case with its recipe and points placed.
    fn placed(&self, name: &str) -> Self {
        Self {
            name: name.to_string(),
            recipe: placed(&self.recipe),
            step: self.step,
            pieces: self.pieces,
            inside: self.inside.iter().map(place).collect(),
            outside: self.outside.iter().map(place).collect(),
            derived: true,
        }
    }
}

// Points of the boss with uprights of radius 3 and a top of radius 1 (the
// test file's): cut by the straight strip over the right face's top edge
// and kept inside its ball; cut by the torus over the corner round about
// (17, 17), 59/20 from its axis and 1/20 below the top, and kept inside the
// tube; the corner round itself, and uncut material.
const BOSS_IN: &Points = &[
    ([39, 2], [10, 1], [19, 2]),
    ([37, 2], [19, 1], [19, 2]),
    ([10, 1], [10, 1], [5, 1]),
];
const BOSS_OUT: &Points = &[
    ([199, 10], [10, 1], [99, 10]),
    ([1877, 100], [1936, 100], [199, 20]),
    ([199, 10], [199, 10], [5, 1]),
];
// The pocket with uprights of radius 3 and a floor of radius 1 (the test
// file's): added along the floor's right edge and the ball's own space;
// added round the corner about (12, 12), 59/20 from its axis and 1/20 above
// the floor, then inside the tube; the pocket's air and uncut material.
const POCKET_IN: &Points = &[
    ([299, 20], [10, 1], [81, 20]),
    ([1377, 100], [1436, 100], [81, 20]),
    ([2, 1], [2, 1], [2, 1]),
];
const POCKET_OUT: &Points = &[
    ([29, 2], [10, 1], [9, 2]),
    ([27, 2], [14, 1], [9, 2]),
    ([10, 1], [10, 1], [7, 1]),
];
// The corner with a front round of radius 3 and a right chain of radius 1
// (the test file's): along the top, over the front round's axis through
// (y, z) = (3, 7), 59/20 out and 1/20 inside the right face, outside the
// tube and then inside it, and down the front.
const CORNER_IN: &Points = &[
    ([39, 2], [15, 1], [19, 2]),
    ([39, 2], [3, 2], [9, 1]),
    ([39, 2], [1, 2], [3, 1]),
    ([10, 1], [10, 1], [5, 1]),
];
const CORNER_OUT: &Points = &[
    ([199, 10], [15, 1], [99, 10]),
    ([399, 20], [123, 100], [936, 100]),
    ([199, 10], [1, 10], [3, 1]),
];

fn cases() -> Vec<Case> {
    let corpus =
        |id: &str, step: &'static str, pieces: usize, inside: &Points, outside: &Points| {
            let row = corpus_row(id);
            let mut case = Case::new(id, parse(&row["recipe"]), step, pieces, inside, outside);
            if let Some(witnesses) = row["witnesses"].as_array() {
                for witness in witnesses {
                    let p: [Rational; 3] = std::array::from_fn(|i| number(&witness["point"][i]));
                    match witness["state"].as_str().unwrap() {
                        "TopAbs_IN" => case.inside.push(p),
                        "TopAbs_OUT" => case.outside.push(p),
                        other => panic!("{id}: witness state {other}"),
                    }
                }
            }
            case
        };
    let boss_31 = Case::new(
        "boss 3/1",
        boss(4, "3", "1", false),
        "top",
        8,
        BOSS_IN,
        BOSS_OUT,
    );
    let corner_31 = Case::new(
        "corner 3/1",
        corner("3", "1", false),
        "right",
        3,
        CORNER_IN,
        CORNER_OUT,
    );
    let mut cases = vec![
        // The top loop rounded at once, mitred at each corner: the
        // elementary route's two strips meet in the mitre plane, as before
        // this step. Its strips register their faces, so it has no rolling
        // piece. Derived: kept inside both strips' balls at the corner over
        // (0, 0, 10), and cut by one strip or both.
        corpus(
            "box-mitre",
            "round",
            0,
            &[
                ([1, 2], [1, 2], [19, 2]),
                ([1, 5], [9, 10], [19, 2]),
                ([10, 1], [10, 1], [5, 1]),
            ],
            &[
                ([1, 10], [1, 10], [99, 10]),
                ([1, 5], [9, 10], [49, 5]),
                ([10, 1], [1, 10], [99, 10]),
            ],
        ),
        // Uprights of radius 2 and a top of radius 1: four horn sectors
        // (centre circle radius 1 = r) about (2, 2), (18, 2), (18, 18) and
        // (2, 18) at z = 9. Derived: cut and kept beside the right strip's
        // centre line (19, y, 9), and at 39/20 and 3/2 from the axis through
        // (18, 18) on the meridian (3/5, 4/5), 1/20 and 1/2 below the top.
        corpus(
            "tangent-boss",
            "top",
            8,
            &[
                ([39, 2], [10, 1], [19, 2]),
                ([189, 10], [96, 5], [19, 2]),
                ([10, 1], [10, 1], [5, 1]),
            ],
            &[
                ([199, 10], [10, 1], [99, 10]),
                ([1917, 100], [1956, 100], [199, 20]),
                ([199, 10], [199, 10], [5, 1]),
            ],
        ),
        // Pocket uprights of radius 2 and a floor of radius 1: four concave
        // horn sectors about (7, 7), (13, 7), (13, 13) and (7, 13) at z = 5.
        // Derived: added and air beside the floor strip's centre line
        // (14, y, 5), and at 39/20 and 3/2 from the axis through (13, 13) on
        // the meridian (3/5, 4/5), 1/20 and 1/2 above the floor.
        corpus(
            "tangent-pocket",
            "floor",
            8,
            &[
                ([299, 20], [10, 1], [81, 20]),
                ([1417, 100], [1456, 100], [81, 20]),
                ([2, 1], [2, 1], [2, 1]),
            ],
            &[
                ([29, 2], [10, 1], [9, 2]),
                ([139, 10], [71, 5], [9, 2]),
                ([10, 1], [10, 1], [7, 1]),
            ],
        ),
        // A front round of radius 2, then the right chain at radius 1: a
        // horn sector about the front round's axis (y, z) = (2, 8). The
        // corpus row's witnesses, checked by OpenCASCADE.
        corpus("unequal-ordered", "right", 3, &[], &[]),
        // A spindle sector (centre circle radius 1 < r = 2) about (17, 17)
        // at z = 8. Derived: cut and kept beside the right strip's centre
        // line (18, y, 8), and at 59/20 and 2 from the axis on the meridian
        // (3/5, 4/5).
        Case::new(
            "boss 3/2",
            boss(4, "3", "2", false),
            "top",
            8,
            &[
                ([19, 1], [10, 1], [9, 1]),
                ([91, 5], [93, 5], [9, 1]),
                ([10, 1], [10, 1], [5, 1]),
            ],
            &[
                ([199, 10], [10, 1], [99, 10]),
                ([1877, 100], [1936, 100], [199, 20]),
            ],
        ),
        Case::new(
            "boss 3/1 named loop",
            boss(3, "3", "1", true),
            "top",
            8,
            BOSS_IN,
            BOSS_OUT,
        ),
        Case::new(
            "pocket 3/1",
            pocket("3", "1"),
            "floor",
            8,
            POCKET_IN,
            POCKET_OUT,
        ),
        Case::new(
            "corner 3/1 from the arc",
            corner("3", "1", true),
            "right",
            3,
            CORNER_IN,
            CORNER_OUT,
        ),
        // A spindle sector about the front round's axis (y, z) = (3, 7) in the
        // plane x = 18. Derived: cut and kept beside the top strip's centre
        // line (18, y, 8), at 2 from the axis on the meridian (-3/5, 4/5) and
        // 1 inside the right face, and beside the front strip's centre line
        // (18, 2, z).
        Case::new(
            "corner 3/2",
            corner("3", "2", false),
            "right",
            3,
            &[
                ([19, 1], [15, 1], [9, 1]),
                ([19, 1], [9, 5], [43, 5]),
                ([19, 1], [1, 1], [3, 1]),
                ([10, 1], [10, 1], [5, 1]),
            ],
            &[
                ([199, 10], [15, 1], [99, 10]),
                ([399, 20], [123, 100], [936, 100]),
                ([199, 10], [1, 10], [3, 1]),
            ],
        ),
    ];
    cases.extend(sphere_corners());
    cases.extend(cliffs());
    // The top loop and then the uprights, all at radius 1: each upright's
    // ball swallows the mitre above it. Derived: kept inside an upright's
    // ball, and cut where the mitre stood at the corner over (0, 0, 10),
    // the corner ball about (1, 1, 9).
    cases.push(corpus(
        "roll-over-mitres",
        "upright",
        4,
        &[([1, 2], [3, 4], [5, 1]), ([3, 2], [1, 2], [19, 2])],
        &[([1, 4], [1, 4], [91, 10]), ([1, 10], [1, 10], [3, 1])],
    ));
    cases.extend(mitre_corners());
    // The top loop at radius 1, then the uprights at radius 2, as the test
    // file's mitred loop 1/2: each upright takes the cliff at both feet,
    // runs on to the top and both rounds cut it. The B-rep refuses the row
    // inside the corpus allowance and builds it without one. Derived, at the
    // corner over (20, 20, 10): cut outside the upright's disk under the
    // x_max round, kept inside both.
    cases.push(corpus(
        "unequal-mitres",
        "upright",
        4,
        &[([39, 2], [37, 2], [19, 2])],
        &[([99, 5], [19, 1], [46, 5])],
    ));
    cases.extend(mitre_cliffs());
    let boss_33 = cases.iter().find(|c| c.name == "boss 3/3").unwrap();
    cases.push(boss_33.placed("placed boss 3/3"));
    let mitred_33 = cases.iter().find(|c| c.name == "mitred 3/3").unwrap();
    cases.push(mitred_33.placed("placed mitred 3/3"));
    cases.push(boss_31.placed("placed boss 3/1"));
    cases.push(corner_31.placed("placed corner 3/1"));
    cases.push(boss_31);
    cases.push(corner_31);
    cases
}

/// A ball as wide as the round it rolls over (`C = 0`): the centre stays on
/// the round's axis, so each sector is one ball, a sphere corner. Derived
/// points: kept inside each strip's ball and the corner ball (on its axis
/// too, just short of the plane it touches), cut or added outside them in
/// the strip and in the corner cell; with the test file's points.
fn sphere_corners() -> Vec<Case> {
    // Ball centres: the boss at (17, 17, 7) and (18, 18, 8), the pocket at
    // (12, 12, 7) and (13, 13, 6), the corner at (17, 3, 7) and (18, 2, 8).
    // Corner points lie on the meridian (3/5, 4/5) of each round's axis.
    let right = |name: &str, recipe: Recipe, inside: &Points, outside: &Points| {
        Case::new(name, recipe, "right", 3, inside, outside)
    };
    vec![
        Case::new(
            "boss 3/3",
            boss(4, "3", "3", false),
            "top",
            8,
            &[
                ([18, 1], [10, 1], [8, 1]),
                ([91, 5], [93, 5], [9, 1]),
                ([1847, 100], [1896, 100], [17, 2]),
                ([17, 1], [17, 1], [999, 100]),
                ([10, 1], [10, 1], [5, 1]),
            ],
            &[
                ([199, 10], [10, 1], [99, 10]),
                ([1877, 100], [1936, 100], [199, 20]),
                ([199, 10], [199, 10], [5, 1]),
            ],
        ),
        Case::new(
            "boss 2/2",
            boss(4, "2", "2", false),
            "top",
            8,
            &[
                ([19, 1], [10, 1], [9, 1]),
                ([93, 5], [94, 5], [9, 1]),
                ([18, 1], [18, 1], [999, 100]),
                ([10, 1], [10, 1], [5, 1]),
            ],
            &[
                ([199, 10], [10, 1], [99, 10]),
                ([1917, 100], [1956, 100], [199, 20]),
                ([199, 10], [199, 10], [5, 1]),
            ],
        ),
        Case::new(
            "pocket 3/3",
            pocket("3", "3"),
            "floor",
            8,
            &[
                ([149, 10], [10, 1], [41, 10]),
                ([1377, 100], [1436, 100], [81, 20]),
                ([2, 1], [2, 1], [2, 1]),
            ],
            &[
                ([14, 1], [10, 1], [5, 1]),
                ([66, 5], [68, 5], [5, 1]),
                ([27, 2], [14, 1], [6, 1]),
                ([12, 1], [12, 1], [401, 100]),
                ([10, 1], [10, 1], [7, 1]),
            ],
        ),
        Case::new(
            "pocket 2/2",
            pocket("2", "2"),
            "floor",
            8,
            &[
                ([149, 10], [10, 1], [41, 10]),
                ([1417, 100], [1456, 100], [81, 20]),
                ([2, 1], [2, 1], [2, 1]),
            ],
            &[
                ([14, 1], [10, 1], [5, 1]),
                ([68, 5], [69, 5], [5, 1]),
                ([13, 1], [13, 1], [401, 100]),
                ([10, 1], [10, 1], [7, 1]),
            ],
        ),
        // The corner's meridian is (-3/5, 4/5) in (y, z) about an axis
        // along x.
        right(
            "corner 3/3",
            corner("3", "3", false),
            &[
                ([18, 1], [15, 1], [8, 1]),
                ([19, 1], [9, 5], [43, 5]),
                ([18, 1], [2, 1], [3, 1]),
                ([1999, 100], [3, 1], [7, 1]),
                ([10, 1], [10, 1], [5, 1]),
            ],
            &[
                ([199, 10], [15, 1], [99, 10]),
                ([399, 20], [123, 100], [936, 100]),
                ([199, 10], [1, 10], [3, 1]),
            ],
        ),
        right(
            "corner 2/2",
            corner("2", "2", false),
            CORNER_22_IN,
            CORNER_22_OUT,
        ),
        right(
            "corner 2/2 from the arc",
            corner("2", "2", true),
            CORNER_22_IN,
            CORNER_22_OUT,
        ),
    ]
}
/// A ball wider than the round its run reaches (RECIPE 4.5.1, the cliff
/// rule): the chain stops there, and its one straight strip keeps both
/// planes, runs on to the other rail's end and the round cuts it. Derived
/// points: kept inside the strip's ball and cut or added outside it, along
/// the run and on its extension over the round, where the round's material
/// (or, in the pocket, its air) meets the strip; the edges the chain does not
/// reach stay sharp. The corner's points past the front round are the test
/// file's.
fn cliffs() -> Vec<Case> {
    vec![
        Case::new(
            "boss 1/3",
            boss(4, "1", "3", false),
            "top",
            1,
            &[
                ([18, 1], [10, 1], [8, 1]),
                ([39, 2], [39, 2], [8, 1]),
                ([10, 1], [199, 10], [99, 10]),
                ([10, 1], [10, 1], [5, 1]),
            ],
            &[([199, 10], [10, 1], [99, 10]), ([39, 2], [39, 2], [99, 10])],
        ),
        Case::new(
            "pocket 1/3",
            pocket("1", "3"),
            "floor",
            1,
            &[
                ([149, 10], [10, 1], [41, 10]),
                ([29, 2], [29, 2], [41, 10]),
                ([2, 1], [2, 1], [2, 1]),
            ],
            &[
                ([14, 1], [10, 1], [5, 1]),
                ([29, 2], [29, 2], [6, 1]),
                ([10, 1], [10, 1], [7, 1]),
            ],
        ),
        Case::new(
            "corner 1/2",
            corner("1", "2", false),
            "right",
            1,
            &[
                ([19, 1], [10, 1], [9, 1]),
                ([19, 1], [1, 2], [9, 1]),
                ([10, 1], [1, 2], [19, 2]),
                ([10, 1], [10, 1], [5, 1]),
            ],
            &[
                ([199, 10], [10, 1], [99, 10]),
                ([199, 10], [1, 2], [93, 10]),
                ([19, 1], [10, 1], [99, 10]),
            ],
        ),
    ]
}
/// An upright whose ball swallows an earlier mitre of two rounds as wide
/// as the ball (RECIPE 4.5.1): its straight strip runs up to the rounds,
/// and the trihedral corner of its two planes and the plane the rounds roll
/// on, less the ball, closes it. Derived points: kept inside the strip's
/// ball and the corner ball, cut or added outside them, in the strip and
/// where the mitre stood (inside both rounds, outside the corner ball); a
/// round away from the corner untouched; with the test file's points.
fn mitre_corners() -> Vec<Case> {
    vec![
        // The ball about (17, 17, 7).
        Case::new(
            "mitred 3/3",
            mitred(4, "3", "3", false),
            "upright",
            1,
            &[
                ([37, 2], [37, 2], [17, 2]),
                ([19, 1], [18, 1], [3, 1]),
                ([19, 1], [18, 1], [8, 1]),
                ([19, 1], [10, 1], [9, 1]),
                ([10, 1], [10, 1], [5, 1]),
            ],
            &[
                ([199, 10], [199, 10], [9, 1]),
                ([19, 1], [19, 1], [99, 10]),
                ([199, 10], [199, 10], [3, 1]),
                ([39, 2], [19, 1], [15, 2]),
            ],
        ),
        // The ball about (18, 18, 8).
        Case::new(
            "mitred 2/2",
            mitred(4, "2", "2", false),
            "upright",
            1,
            &[
                ([19, 1], [19, 1], [9, 1]),
                ([39, 2], [37, 2], [3, 1]),
                ([10, 1], [10, 1], [5, 1]),
            ],
            &[
                ([199, 10], [199, 10], [9, 1]),
                ([39, 2], [39, 2], [17, 2]),
                ([199, 10], [199, 10], [4, 1]),
            ],
        ),
        // Balls about (3, 3, 7), (17, 3, 7), (17, 17, 7) and (3, 17, 7).
        Case::new(
            "mitred loop 3/3",
            mitred(4, "3", "3", true),
            "upright",
            4,
            &[
                ([37, 2], [37, 2], [17, 2]),
                ([3, 2], [3, 2], [17, 2]),
                ([10, 1], [10, 1], [5, 1]),
            ],
            &[
                ([1, 10], [1, 10], [9, 1]),
                ([1, 2], [39, 2], [15, 2]),
                ([39, 2], [1, 2], [15, 2]),
            ],
        ),
        // Concave: balls about (6, 6, 5) and (14, 14, 5) among four. Where
        // the mitre's air stood, inside both floor rounds and outside the
        // corner ball, the corner adds material.
        Case::new(
            "mitred pocket 1/1",
            mitred_pocket("1", "1"),
            "upright",
            4,
            &[
                ([21, 4], [21, 4], [49, 10]),
                ([59, 4], [59, 4], [49, 10]),
                ([51, 10], [51, 10], [7, 1]),
                ([2, 1], [2, 1], [2, 1]),
            ],
            &[
                ([11, 2], [11, 2], [9, 2]),
                ([11, 2], [11, 2], [7, 1]),
                ([10, 1], [10, 1], [7, 1]),
            ],
        ),
    ]
}
/// An upright wider than both rounds of an earlier mitre (RECIPE 4.5.1,
/// the cliff rule at both feet): its strip keeps both planes and runs on
/// past the cliff to the plane the rounds roll on, and both rounds cut it.
/// Derived points: cut or added outside the strip's disk below the cliff
/// and above it inside a round (in the pocket, inside a floor round's air),
/// kept inside both the disk and the round (where the mitre's outer part
/// survives); uncut material.
fn mitre_cliffs() -> Vec<Case> {
    vec![
        // The upright's disk about (16, 16); the rounds' axes at z = 7.
        Case::new(
            "mitred 3/4",
            mitred(4, "3", "4", false),
            "upright",
            1,
            &[
                ([18, 1], [17, 1], [8, 1]),
                ([37, 2], [37, 2], [9, 1]),
                ([10, 1], [10, 1], [5, 1]),
            ],
            &[([39, 2], [39, 2], [3, 1]), ([39, 2], [19, 1], [15, 2])],
        ),
        // Concave: the disk about (13/2, 13/2), the floor rounds' axes at
        // z = 5.
        Case::new(
            "mitred pocket 1/(3/2)",
            mitred_pocket("1", "3/2"),
            "upright",
            4,
            &[
                ([51, 10], [51, 10], [7, 1]),
                ([51, 10], [11, 2], [23, 5]),
                ([2, 1], [2, 1], [2, 1]),
            ],
            &[([6, 1], [6, 1], [9, 2]), ([10, 1], [10, 1], [7, 1])],
        ),
    ]
}
const CORNER_22_IN: &Points = &[
    ([19, 1], [15, 1], [9, 1]),
    ([19, 1], [7, 5], [44, 5]),
    ([19, 1], [11, 10], [46, 5]),
    ([19, 1], [1, 1], [3, 1]),
    ([1999, 100], [2, 1], [8, 1]),
    ([10, 1], [10, 1], [5, 1]),
];
const CORNER_22_OUT: &Points = &[
    ([199, 10], [15, 1], [99, 10]),
    ([399, 20], [83, 100], [956, 100]),
    ([199, 10], [43, 50], [238, 25]),
    ([199, 10], [1, 10], [3, 1]),
];
/// Rational points of each sphere corner's face, `c + (1, 2, 2)` and its
/// permutations within the sector. The B-rep puts them on its boundary.
/// In the field exactly one rolling piece, the sector, holds each on its
/// region's boundary, and its distance range starts at 0. The plan puts them
/// on its boundary too: the sector's sphere is a registered face, whose
/// equation settles the side of each ray from the point.
#[test]
fn a_sphere_corner_face_holds_its_rational_points() {
    let points: [(&str, Recipe, [[i64; 3]; 3]); 3] = [
        (
            "boss 3/3",
            boss(4, "3", "3", false),
            [[18, 19, 9], [19, 18, 9], [19, 19, 8]],
        ),
        (
            "pocket 3/3",
            pocket("3", "3"),
            [[13, 14, 5], [14, 13, 5], [14, 14, 6]],
        ),
        (
            "corner 3/3",
            corner("3", "3", false),
            [[19, 1, 8], [18, 1, 9], [19, 2, 9]],
        ),
    ];
    let width = ratio(1, 1_000_000_000_000);
    for (name, recipe, list) in points {
        let plan = Plan::compile_checked(&recipe, &mut decider()).unwrap();
        let mut built = evaluate(&recipe, &mut decider()).unwrap();
        for p in list {
            let p = p.map(Rational::from_i64);
            assert_eq!(
                brep_location(&mut built, &p),
                Location::Boundary,
                "{name}: the B-rep at {p:?}"
            );
            assert_eq!(location(&plan, &p), Location::Boundary, "{name} at {p:?}");
            let point = Point3::new(
                p[0].clone().into(),
                p[1].clone().into(),
                p[2].clone().into(),
            );
            let mut on = Vec::new();
            for piece in plan.rolling_pieces() {
                match piece.classify(&point, &mut decider()).unwrap() {
                    Relation::Boundary => on.push(piece),
                    Relation::Exterior => {}
                    Relation::Interior => panic!("{name}: {p:?} inside a piece's region"),
                }
            }
            let [piece] = on[..] else {
                panic!("{name}: {p:?} on {} pieces' boundaries", on.len());
            };
            let id = plan
                .select(piece.origin())
                .expect("the sector's face is registered");
            assert!(matches!(
                plan.faces()[id.index()].surface,
                Surface::Sphere(_)
            ));
            let query = plan
                .face_query_checked(
                    id,
                    &RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone()),
                    &mut Workspace::default(),
                    &mut decider(),
                )
                .unwrap();
            assert_eq!(query.state, FaceState::OnFace, "{name} at {p:?}");
            assert!(
                matches!(query.direction, Direction::Unique(_)),
                "{name} at {p:?}: {query:?}"
            );
            let range = piece.distance(&point, &width, &mut decider()).unwrap();
            assert_eq!(range.lower, Rational::zero(), "{name} at {p:?}");
            assert!(
                range.upper.as_ref().is_some_and(|u| *u <= width),
                "{name} at {p:?}: {range:?}"
            );
        }
    }
}

/// The B-rep's membership of a point, from its exact point classifier.
fn brep_location(built: &mut catallax_features::Built, p: &[Rational; 3]) -> Location {
    use catallax_boolean::{Containment, classify_point};
    use catallax_numbers::Outcome;
    let point = RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone());
    match classify_point(&mut built.model, built.result, &point, &mut decider()) {
        Outcome::Proven(Containment::Inside) => Location::Inside,
        Outcome::Proven(Containment::Outside) => Location::Outside,
        Outcome::ProvenCoincident(_) => Location::Boundary,
        other => panic!("the B-rep does not settle {p:?}: {other:?}"),
    }
}

fn location(plan: &Plan, p: &[Rational; 3]) -> Location {
    let point = RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone());
    plan.proven_checked(&point, &mut Workspace::default(), &mut decider())
        .unwrap()
        .location
}

/// Seam control at the junctions where a chain meets the ring sector over an
/// earlier round (reviewer note on #738): a point on the junction plane and
/// its neighbours 1/1000 either side along the plane's normal take the same
/// answer, in the material the chain removes or adds and inside its ball,
/// where the plane is no face. Each point is checked against the B-rep
/// first.
#[test]
fn junctions_over_an_earlier_round_add_no_seam_face() {
    // (case, junction point, axis of the junction normal, answer)
    let junctions: [(&str, [Rational; 3], usize, Location); 8] = [
        // `tangent-boss`: the right strip meets the sector over (18, 18) in
        // the plane y = 18.
        (
            "tangent-boss",
            [ratio(199, 10), ratio(18, 1), ratio(99, 10)],
            1,
            Location::Outside,
        ),
        (
            "tangent-boss",
            [ratio(39, 2), ratio(18, 1), ratio(19, 2)],
            1,
            Location::Inside,
        ),
        // `tangent-pocket`: the floor strip along x = 15 meets the sector
        // over (13, 13) in the plane y = 13.
        (
            "tangent-pocket",
            [ratio(299, 20), ratio(13, 1), ratio(81, 20)],
            1,
            Location::Inside,
        ),
        (
            "tangent-pocket",
            [ratio(29, 2), ratio(13, 1), ratio(9, 2)],
            1,
            Location::Outside,
        ),
        // `unequal-ordered`: the top strip meets the sector over the front
        // round in the plane y = 2, and the sector meets the front strip
        // in the plane z = 8.
        (
            "unequal-ordered",
            [ratio(199, 10), ratio(2, 1), ratio(99, 10)],
            1,
            Location::Outside,
        ),
        (
            "unequal-ordered",
            [ratio(39, 2), ratio(2, 1), ratio(19, 2)],
            1,
            Location::Inside,
        ),
        (
            "unequal-ordered",
            [ratio(199, 10), ratio(1, 10), ratio(8, 1)],
            2,
            Location::Outside,
        ),
        (
            "unequal-ordered",
            [ratio(39, 2), ratio(1, 2), ratio(8, 1)],
            2,
            Location::Inside,
        ),
    ];
    let mut plans: BTreeMap<&str, (Plan, catallax_features::Built)> = BTreeMap::new();
    for (name, point, axis, want) in junctions {
        if !plans.contains_key(name) {
            let recipe = parse(&corpus_row(name)["recipe"]);
            let plan = Plan::compile_checked(&recipe, &mut decider()).unwrap();
            plans.insert(name, (plan, evaluate(&recipe, &mut decider()).unwrap()));
        }
        let (plan, built) = plans.get_mut(name).unwrap();
        for step in [-1, 0, 1] {
            let mut p = point.clone();
            p[axis] = &p[axis] + &ratio(step, 1000);
            assert_eq!(brep_location(built, &p), want, "{name}: the B-rep at {p:?}");
            assert_eq!(location(plan, &p), want, "{name} at {p:?}");
        }
    }
}

#[test]
fn the_field_rolls_each_chain_over_its_earlier_round_and_agrees_on_membership() {
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
        let mut built = case
            .derived
            .then(|| evaluate(&case.recipe, &mut decider()).unwrap());
        for (points, want) in [
            (&case.inside, Location::Inside),
            (&case.outside, Location::Outside),
        ] {
            for p in points {
                if let Some(built) = built.as_mut() {
                    let brep = brep_location(built, p);
                    assert_eq!(brep, want, "{}: the B-rep at {p:?}", case.name);
                }
                let found = location(&plan, p);
                if found != want {
                    failures.push(format!("{} at {p:?}: {found:?}", case.name));
                }
            }
        }
        eprintln!(
            "{}: {} witnesses in {:?}",
            case.name,
            case.inside.len() + case.outside.len(),
            start.elapsed()
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// The cases whose B-rep blend vertices are checked, one test per group so
/// they run side by side.
const VERTEX_GROUPS: &[&[&str]] = &[
    &["tangent-boss", "boss 3/1"],
    &["tangent-pocket", "pocket 3/1"],
    &["unequal-ordered", "corner 3/1", "corner 3/1 from the arc"],
    &["boss 3/2", "corner 3/2", "boss 3/1 named loop"],
    &["placed boss 3/1", "placed corner 3/1"],
    &["boss 3/3", "boss 2/2"],
    &["pocket 3/3", "pocket 2/2"],
    &["corner 3/3", "corner 2/2", "corner 2/2 from the arc"],
    &["placed boss 3/3"],
    &["boss 1/3", "pocket 1/3", "corner 1/2"],
    &["roll-over-mitres", "mitred 3/3", "mitred 2/2"],
    &["mitred loop 3/3", "mitred pocket 1/1", "placed mitred 3/3"],
    &["unequal-mitres", "mitred 3/4", "mitred pocket 1/(3/2)"],
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
    boss_blend_vertices_lie_on_the_field_face = 0;
    pocket_blend_vertices_lie_on_the_field_face = 1;
    corner_blend_vertices_lie_on_the_field_face = 2;
    spindle_blend_vertices_lie_on_the_field_face = 3;
    placed_blend_vertices_lie_on_the_field_face = 4;
    sphere_boss_blend_vertices_lie_on_the_field_face = 5;
    sphere_pocket_blend_vertices_lie_on_the_field_face = 6;
    sphere_corner_blend_vertices_lie_on_the_field_face = 7;
    placed_sphere_blend_vertices_lie_on_the_field_face = 8;
    cliff_blend_vertices_lie_on_the_field_face = 9;
    mitre_blend_vertices_lie_on_the_field_face = 10;
    mitre_loop_blend_vertices_lie_on_the_field_face = 11;
    mitre_cliff_blend_vertices_lie_on_the_field_face = 12;
}

#[test]
fn every_case_has_its_vertices_checked() {
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

/// The display vertices of the B-rep's blend faces made by `step`. Every
/// such face here is an exact cylinder or torus, so the part's tolerance at
/// each vertex is 1 nm.
fn blend_vertices(recipe: &Recipe, step: &str, chord: &str) -> Vec<Point3> {
    step_vertices(recipe, step, chord, false)
}

/// The display vertices of the B-rep's blend faces made by `step`, or of
/// its corner faces when `corners` is set.
fn step_vertices(recipe: &Recipe, step: &str, chord: &str, corners: bool) -> Vec<Point3> {
    let built = evaluate(recipe, &mut decider()).unwrap();
    let mesh = mesh_body(
        &built.model,
        built.result,
        &MeshOptions {
            chord_tolerance: Rational::parse(chord).unwrap(),
            ..MeshOptions::default()
        },
        &mut decider(),
    )
    .unwrap()
    .value;
    let mut indices = std::collections::BTreeSet::new();
    for face in &mesh.faces {
        if !built.model.face_origins(face.face).iter().any(|o| {
            let mut origin = o;
            while let Role::Copy { of, .. } = &origin.role {
                origin = of;
            }
            origin.step.as_str() == step
                && if corners {
                    matches!(origin.role, Role::Corner { .. })
                } else {
                    matches!(origin.role, Role::Blend { .. })
                }
        }) {
            continue;
        }
        for triangle in &mesh.triangles[face.triangles.clone()] {
            indices.extend(triangle.vertices);
        }
    }
    indices
        .into_iter()
        .map(|v| {
            let [x, y, z] = mesh.positions[v].map(|c| Rational::from_f64_exact(c).unwrap());
            Point3::new(x.into(), y.into(), z.into())
        })
        .collect()
}

/// Upper bounds on the distance from each point to the nearest of the
/// pieces' faces, from guaranteed ranges at most `width` wide, and the most
/// work one query took. The pieces are compiled once and every query reuses
/// them. Every 128th point is also queried at 1, 2 and 4 workers, whose
/// records must agree.
fn nearest_witnesses(
    name: &str,
    pieces: &[&Piece],
    points: &[Point3],
    width: &Rational,
) -> (Vec<Option<Rational>>, u64) {
    const CAP: u64 = 20_000;
    let allowance = |n: usize| OperationAllowance::reserve(n as u64 * CAP).unwrap();
    let mut nearest = vec![None::<Rational>; points.len()];
    let mut most = 0;
    for &piece in pieces {
        let index: Vec<usize> = (0..points.len())
            .filter(|&i| nearest[i].as_ref().is_none_or(|u| u > width))
            .collect();
        if index.is_empty() {
            break;
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
            most = most.max(record.work);
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
    (nearest, most)
}

/// Every B-rep display vertex of the rolling step's blend faces lies within
/// the field's guaranteed range of a rolling piece's face, at 1 nm.
fn blend_vertices_lie_on_the_field_face(names: &[&str]) {
    let tolerance = ratio(1, 1_000_000);
    for case in cases()
        .into_iter()
        .filter(|c| names.contains(&c.name.as_str()))
    {
        let start = std::time::Instant::now();
        let plan = Plan::compile_checked(&case.recipe, &mut decider()).unwrap();
        let points = blend_vertices(&case.recipe, case.step, "1/20");
        assert!(!points.is_empty(), "{}: no blend vertices", case.name);
        let pieces: Vec<&Piece> = plan.rolling_pieces().collect();
        let (nearest, most) = nearest_witnesses(&case.name, &pieces, &points, &tolerance);
        let mut worst = Rational::zero();
        for found in &nearest {
            let found = found
                .clone()
                .unwrap_or_else(|| panic!("{}: no witness", case.name));
            assert!(found <= tolerance, "{}: {found:?}", case.name);
            worst = worst.max(found);
        }
        eprintln!(
            "{}: {} blend vertices in {:?}, worst field distance {}, at most {most} work a query",
            case.name,
            nearest.len(),
            start.elapsed(),
            worst.to_f64_nearest()
        );
    }
}

/// One mitre corner derived here: its ball's centre, the signs of the
/// cube's directions along x, y and z, and its radius.
type Corner = ([i64; 3], [i64; 3], i64);

/// Each mitre corner registers its sphere face under a corner name of the
/// rolling step, about the centre and with the radius derived here, and
/// every display vertex of the B-rep's corner faces of that step lies
/// within 1 nm of one of those spheres and at most 1 nm outside the cone of
/// its cell, whose directions are derived here: the radial projection of a
/// point of that cone lies on the face. The test reads the field's
/// registered faces directly; the distance to a sphere is exact.
#[test]
fn mitre_corner_vertices_lie_on_the_field_corner_face() {
    let tolerance = ratio(1, 1_000_000);
    let corners: [(&str, Vec<Corner>); 6] = [
        (
            "roll-over-mitres",
            vec![
                ([1, 1, 9], [-1, -1, 1], 1),
                ([19, 1, 9], [1, -1, 1], 1),
                ([19, 19, 9], [1, 1, 1], 1),
                ([1, 19, 9], [-1, 1, 1], 1),
            ],
        ),
        ("mitred 3/3", vec![([17, 17, 7], [1, 1, 1], 3)]),
        ("mitred 2/2", vec![([18, 18, 8], [1, 1, 1], 2)]),
        (
            "mitred loop 3/3",
            vec![
                ([3, 3, 7], [-1, -1, 1], 3),
                ([17, 3, 7], [1, -1, 1], 3),
                ([17, 17, 7], [1, 1, 1], 3),
                ([3, 17, 7], [-1, 1, 1], 3),
            ],
        ),
        (
            "mitred pocket 1/1",
            vec![
                ([6, 6, 5], [-1, -1, -1], 1),
                ([14, 6, 5], [1, -1, -1], 1),
                ([14, 14, 5], [1, 1, -1], 1),
                ([6, 14, 5], [-1, 1, -1], 1),
            ],
        ),
        ("placed mitred 3/3", vec![([17, 17, 7], [1, 1, 1], 3)]),
    ];
    let cases = cases();
    for (name, expected) in corners {
        let case = cases.iter().find(|c| c.name == name).unwrap();
        let placed = name.starts_with("placed");
        // Each corner's centre and its three unit directions, moved as the
        // case's blocks are.
        let expected: Vec<([Rational; 3], [[Rational; 3]; 3], Rational)> = expected
            .iter()
            .map(|(c, signs, r)| {
                let c = c.map(Rational::from_i64);
                let mut directions: [[Rational; 3]; 3] = std::array::from_fn(|i| {
                    std::array::from_fn(|j| Rational::from_i64(if i == j { signs[i] } else { 0 }))
                });
                if placed {
                    let origin = place(&[Rational::zero(), Rational::zero(), Rational::zero()]);
                    for d in &mut directions {
                        let moved = place(d);
                        *d = std::array::from_fn(|i| moved[i].clone() - origin[i].clone());
                    }
                }
                let c = if placed { place(&c) } else { c };
                (c, directions, Rational::from_i64(*r))
            })
            .collect();
        let plan = Plan::compile_checked(&case.recipe, &mut decider()).unwrap();
        let mut spheres: Vec<([Rational; 3], Rational)> = plan
            .faces()
            .iter()
            .filter(|f| {
                f.origin.step.as_str() == case.step && matches!(f.origin.role, Role::Corner { .. })
            })
            .map(|f| {
                let Surface::Sphere(s) = &f.surface else {
                    panic!("{name}: a corner face that is no sphere");
                };
                let c = [&s.center.x, &s.center.y, &s.center.z]
                    .map(|x| x.as_rational().expect("a rational centre"));
                (c, s.radius.as_rational().expect("a rational radius"))
            })
            .collect();
        spheres.sort();
        let mut want: Vec<_> = expected
            .iter()
            .map(|(c, _, r)| (c.clone(), r.clone()))
            .collect();
        want.sort();
        assert_eq!(spheres, want, "{name}: corner faces");
        let points = step_vertices(&case.recipe, case.step, "1/20", true);
        assert!(!points.is_empty(), "{name}: no corner vertices");
        for p in &points {
            let p = [&p.x, &p.y, &p.z].map(|x| x.as_rational().unwrap());
            let near = expected.iter().any(|(c, directions, r)| {
                let u: [Rational; 3] = std::array::from_fn(|i| p[i].clone() - c[i].clone());
                let squared = u
                    .iter()
                    .fold(Rational::zero(), |sum, x| sum + x.clone() * x.clone());
                let lo = r.clone() - tolerance.clone();
                let hi = r.clone() + tolerance.clone();
                let within = lo.clone() * lo <= squared && squared <= hi.clone() * hi;
                let cone = directions.iter().all(|d| {
                    let along =
                        (0..3).fold(Rational::zero(), |sum, i| sum + u[i].clone() * d[i].clone());
                    along >= -tolerance.clone()
                });
                within && cone
            });
            assert!(near, "{name}: {p:?} is not on a field corner face");
        }
        eprintln!("{name}: {} corner vertices within 1 nm", points.len());
    }
}

/// Rational points of a mitre corner's sphere face, `c + (1, 2, 2)` and its
/// permutations, which lie in the corner's cone. The B-rep puts them on its
/// boundary, and the field's plan does too, on the registered corner face.
#[test]
fn a_mitre_corner_face_holds_its_rational_points() {
    let recipe = mitred(4, "3", "3", false);
    let plan = Plan::compile_checked(&recipe, &mut decider()).unwrap();
    let mut built = evaluate(&recipe, &mut decider()).unwrap();
    for p in [[18, 19, 9], [19, 18, 9], [19, 19, 8]] {
        let p = p.map(Rational::from_i64);
        assert_eq!(
            brep_location(&mut built, &p),
            Location::Boundary,
            "the B-rep at {p:?}"
        );
        assert_eq!(
            location(&plan, &p),
            Location::Boundary,
            "the field at {p:?}"
        );
    }
}

/// The face registry (4.2a): a straight piece's face is the cylinder about
/// its centre line, registered under the piece's name. The upright under the
/// mitred corner at 3/3 has its centre line at `x = y = 17`; the rational
/// points `c + 3 (3/5, 4/5, 0)` and their mirrors about the diagonal lie on
/// the B-rep's boundary. The field's plan puts them on its boundary, the
/// registered face holds them with the unique outward normal away from the
/// centre line, and points 1/1000 off the face along that normal are cut
/// outside the ball and kept inside it, as in the B-rep.
#[test]
fn a_straight_piece_face_is_a_registered_face() {
    let recipe = mitred(4, "3", "3", false);
    let plan = Plan::compile_checked(&recipe, &mut decider()).unwrap();
    let mut built = evaluate(&recipe, &mut decider()).unwrap();
    let registered: Vec<_> = (0..plan.faces().len())
        .filter(|&i| {
            plan.result_face(i).is_some()
                && matches!(plan.faces()[i].surface, Surface::Cylinder(_))
                && plan
                    .rolling_pieces()
                    .any(|piece| piece.origin() == &plan.faces()[i].origin)
        })
        .collect();
    let [index] = registered[..] else {
        panic!("{} registered straight faces", registered.len());
    };
    let id = plan.result_face(index).unwrap();
    for (x, y, normal) in [(94, 97, [3, 4]), (97, 94, [4, 3])] {
        for z in [1, 3, 6] {
            let p = [ratio(x, 5), ratio(y, 5), ratio(z, 1)];
            assert_eq!(brep_location(&mut built, &p), Location::Boundary, "{p:?}");
            assert_eq!(location(&plan, &p), Location::Boundary, "{p:?}");
            let point = RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone());
            let query = plan
                .face_query_checked(id, &point, &mut Workspace::default(), &mut decider())
                .unwrap();
            assert_eq!(query.state, FaceState::OnFace, "{p:?}");
            let Direction::Unique(direction) = query.direction else {
                panic!("{p:?}: {query:?}");
            };
            let expected = [normal[0] as f64 / 5.0, normal[1] as f64 / 5.0, 0.0];
            for (component, want) in direction.iter().zip(expected) {
                assert!(
                    (component.lo() - want).abs() < 1e-9 && (component.hi() - want).abs() < 1e-9,
                    "{p:?}: {direction:?}"
                );
            }
            for (sign, want) in [(1, Location::Outside), (-1, Location::Inside)] {
                let off = [
                    &p[0] + &ratio(sign * normal[0], 5000),
                    &p[1] + &ratio(sign * normal[1], 5000),
                    p[2].clone(),
                ];
                assert_eq!(
                    brep_location(&mut built, &off),
                    want,
                    "the B-rep at {off:?}"
                );
                assert_eq!(location(&plan, &off), want, "{off:?}");
            }
        }
    }
}

// The other recipes of `fillet_over_fillet.rs`.
/// The block's top edges along x_max and y_max, or its whole top loop,
/// rounded together with mitred corners, then the upright under each mitre.
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
/// The pocket's floor loop rounded with mitred corners, then its uprights.
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
/// An L whose inside corner is rounded, then the right edge of the base's
/// top, which rolls over the inner round's end and up the wall.
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
/// A boss on a plate with its uprights rounded, then its base edge: a
/// concave loop over convex rounds.
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
/// The block after `first`, then every edge at `rest`.
fn every(first: Value, rest: &str) -> Recipe {
    recipe(
        4,
        vec![
            block("block", ["20", "20", "10"], ["0", "0", "0"]),
            first,
            json!({"id":"rest","op":"fillet","body":"first","radius":rest,"edges":"all"}),
        ],
    )
}
/// A third round that rolls over the right chain of `corner("3", "2")`.
fn corner_bottom() -> Recipe {
    let mut value = serde_json::to_value(corner("3", "2", false)).unwrap();
    value["steps"].as_array_mut().unwrap().push(fillet(
        "bottom",
        "right",
        "1",
        vec![edge(cap("block", "start"), side("block", "x_max"))],
    ));
    value["result"] = json!("bottom");
    parse(&value)
}

/// The step 5 recipes left for later parts of the step, each with the
/// field's code today. The B-rep builds every one, so each is a field
/// refusal still to clear: a ball outside the round's cylinder (a convex
/// chain over a concave round, or the reverse) crosses the round's tangent
/// faces; an upright under an earlier mitre of rounds wider than its ball,
/// which rolls over it;
/// every edge at once; and a round over a chain of rolling pieces, whose
/// faces the field's fillet routes do not take as supports to roll on or to
/// clear (`field.fillet-contact`).
fn later() -> Vec<(&'static str, Recipe, &'static str)> {
    let overflow = "fillet.unsupported-overflow";
    let contact = "fillet.unsupported-contact";
    let front = |radius: &str| {
        fillet(
            "first",
            "block",
            radius,
            vec![edge(cap("block", "end"), side("block", "y_min"))],
        )
    };
    vec![
        (
            "mixed-ordered",
            parse(&corpus_row("mixed-ordered")["recipe"]),
            overflow,
        ),
        ("bracket 3/1", bracket("3", "1"), overflow),
        ("bracket 3/4", bracket("3", "4"), overflow),
        ("bracket 2/(3/2)", bracket("2", "3/2"), overflow),
        ("based 2/1", based("2", "1"), overflow),
        ("based 2/3", based("2", "3"), overflow),
        ("based 3/1", based("3", "1"), overflow),
        ("mitred loop 1/(1/2)", mitred(4, "1", "1/2", true), overflow),
        ("mitred 3/2", mitred(4, "3", "2", false), overflow),
        ("mitred pocket 1/(1/2)", mitred_pocket("1", "1/2"), overflow),
        (
            "every after a front round 3/1",
            every(front("3"), "1"),
            overflow,
        ),
        (
            "every after a front round 2/2",
            every(front("2"), "2"),
            contact,
        ),
        (
            "every after uprights",
            every(fillet("first", "block", "3", uprights("block")), "1"),
            overflow,
        ),
        (
            "every after a mitred loop",
            every(
                fillet(
                    "first",
                    "block",
                    "2",
                    SIDES
                        .iter()
                        .map(|s| edge(cap("block", "end"), side("block", s)))
                        .collect(),
                ),
                "2",
            ),
            contact,
        ),
        (
            "corner bottom over the right chain",
            corner_bottom(),
            "field.fillet-contact",
        ),
    ]
}

/// The recipes of `fillet_over_fillet.rs` that both kernels refuse, with
/// the field's code: a ball wider than a round chosen by name, and, in
/// version 3, a strip that does not roll over and whose carrier passes the
/// earlier round or mitre.
fn both_refuse() -> Vec<(&'static str, Recipe, &'static str)> {
    vec![
        (
            "boss named loop 1/3",
            boss(4, "1", "3", true),
            "field.rolling-path-reach",
        ),
        (
            "corner 1/2 from the arc",
            corner("1", "2", true),
            "field.rolling-path-reach",
        ),
        (
            "boss version 3 3/1",
            boss(3, "3", "1", false),
            "fillet.unsupported-overflow",
        ),
        (
            "mitred version 3 3/3",
            mitred(3, "3", "3", false),
            "field.fillet-contact",
        ),
    ]
}

#[test]
fn the_field_names_what_it_has_not_rolled_over_yet() {
    for list in [later(), both_refuse()] {
        let mut found = BTreeMap::new();
        let mut want = BTreeMap::new();
        for (name, recipe, code) in list {
            let answer = match Plan::compile_checked(&recipe, &mut decider()) {
                Ok(_) => "built".to_string(),
                Err(refusal) => refusal.code.to_string(),
            };
            found.insert(name, answer);
            want.insert(name, code.to_string());
        }
        assert_eq!(found, want);
    }
}
