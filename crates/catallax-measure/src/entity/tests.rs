//! Closed-form answers for each carrier pair, trimmed faces whose minimum lies
//! on their boundaries, and the two hole walls from the real-use pass.

use crate::{Accuracy, Enclosure, PairEntity, PairRequest, entity_pair};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{Curve, RPoint3, RVector3, Surface};
use catallax_topology::{BodyId, EdgeId, FaceId, Model};
use serde_json::{Value, json};

use PairEntity::{Edge, Face};
use PairRequest::{AngleDegrees, CentreDistance, MinimumDistance};

/// The SolidDraft application allowance for one pair readout.
const WORK: u64 = 20_000_000;

fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn p(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(q(x), q(y), q(z))
}
fn v(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(q(x), q(y), q(z))
}

fn pose(x: i64, y: i64, z: i64) -> Value {
    json!({"rotation": ["1", "0", "0", "0"],
           "translation": [x.to_string(), y.to_string(), z.to_string()]})
}

/// A quarter turn about the x axis takes local z to -y.
fn turned(x: i64, y: i64, z: i64) -> Value {
    json!({"rotation": ["1", "1", "0", "0"],
           "translation": [x.to_string(), y.to_string(), z.to_string()]})
}

fn build(steps: Vec<Value>, result: &str) -> (Model, BodyId) {
    let recipe = json!({"format": "catallax-recipe/3", "name": "pair readouts",
                        "steps": steps, "result": result});
    let recipe = Recipe::parse(&recipe.to_string()).unwrap();
    let body =
        catallax_app::build_incremental_without_display(&recipe, None, &mut Decider::default())
            .unwrap();
    let built = body.checkpoint.built();
    (built.model.clone(), built.result)
}

fn block() -> Value {
    json!({"id": "block", "op": "block", "size": ["60", "40", "20"], "pose": pose(-30, -20, 0)})
}

fn hole(
    id: &str,
    target: &str,
    radius: i64,
    point: RPoint3,
    direction: RVector3,
    depth: Value,
) -> Value {
    let s = |r: &Rational| r.to_decimal_exact().unwrap();
    json!({"id": id, "op": "hole", "target": target, "diameter": (2 * radius).to_string(),
           "depth": depth, "tip": "flat",
           "at": [{"id": "h", "point": [s(&point.x), s(&point.y), s(&point.z)],
                   "direction": [s(&direction.x), s(&direction.y), s(&direction.z)]}]})
}

fn union(id: &str, target: &str, tool: &str) -> Value {
    json!({"id": id, "op": "boolean", "kind": "union", "target": target, "tool": tool})
}

/// The unique face of `body` whose carrier satisfies `want`.
fn face(model: &Model, body: BodyId, want: impl Fn(&Surface) -> bool) -> FaceId {
    let found: Vec<_> = model
        .faces_of_body(body)
        .into_iter()
        .filter(|f| want(model.geom().surface(model.face(*f).surface)))
        .collect();
    assert_eq!(found.len(), 1, "expected one matching face");
    found[0]
}

fn rational(p: &catallax_shapes::Point3) -> RPoint3 {
    p.as_rational().expect("rational point")
}
fn rvec(v: &catallax_shapes::Vector3) -> RVector3 {
    v.as_rational().expect("rational vector")
}

/// A cylinder whose axis passes through `at` along `dir`.
fn wall(at: RPoint3, dir: RVector3) -> impl Fn(&Surface) -> bool {
    move |s| match s {
        Surface::Cylinder(c) => {
            let axis = rvec(&c.axis);
            axis.cross(&dir).is_zero() && rational(&c.origin).minus(&at).cross(&dir).is_zero()
        }
        _ => false,
    }
}

/// A plane through `at` with normal parallel to `normal`.
fn plane(at: RPoint3, normal: RVector3) -> impl Fn(&Surface) -> bool {
    move |s| match s {
        Surface::Plane(pl) => {
            rvec(&pl.normal).cross(&normal).is_zero()
                && rational(&pl.origin).minus(&at).dot(&normal).is_zero()
        }
        _ => false,
    }
}

fn kind(s: &Surface) -> &'static str {
    match s {
        Surface::Sphere(_) => "sphere",
        Surface::Torus(_) => "torus",
        Surface::Cone(_) => "cone",
        _ => "other",
    }
}

/// The unique edge of `body` whose curve satisfies `want`.
fn edge(model: &Model, body: BodyId, want: impl Fn(&Model, EdgeId) -> bool) -> EdgeId {
    let found: Vec<_> = model
        .edges_of_body(body)
        .into_iter()
        .filter(|e| want(model, *e))
        .collect();
    assert_eq!(found.len(), 1, "expected one matching edge");
    found[0]
}

fn ask(
    model: &Model,
    a: PairEntity,
    b: PairEntity,
    request: PairRequest,
) -> Result<Enclosure, catallax_numbers::Refusal> {
    let (answer, work) = with_operation_budget(WORK, || {
        entity_pair(
            model,
            a,
            b,
            request,
            &Accuracy::default(),
            &mut Decider::default(),
        )
    });
    assert!(
        !work.exhausted,
        "{request:?} {a:?} {b:?} exceeded the allowance"
    );
    let answer = answer?;
    assert!(answer.marks.is_empty(), "{:?}", answer.marks);
    // Swapping the entities changes nothing.
    let swapped = entity_pair(
        model,
        b,
        a,
        request,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(swapped.value, answer.value);
    Ok(answer.value)
}

fn exact(model: &Model, a: PairEntity, b: PairEntity, request: PairRequest, want: Rational) {
    let got = ask(model, a, b, request).unwrap();
    assert!(
        got.is_exact() && got.lo == want,
        "{request:?}: {got:?}, want {want:?}"
    );
}

/// Checks `lo <= sqrt(square) - minus <= hi` exactly, with a narrow nonzero width.
fn root_minus(model: &Model, a: PairEntity, b: PairEntity, square: Rational, minus: Rational) {
    root_of(model, a, b, MinimumDistance, square, minus);
}

fn root_of(
    model: &Model,
    a: PairEntity,
    b: PairEntity,
    request: PairRequest,
    square: Rational,
    minus: Rational,
) {
    let got = ask(model, a, b, request).unwrap();
    let lo = &got.lo + &minus;
    let hi = &got.hi + &minus;
    assert!(lo.sign() != catallax_numbers::Sign::Negative, "{got:?}");
    assert!(
        &lo * &lo <= square && square <= &hi * &hi,
        "{got:?} vs sqrt({square:?}) - {minus:?}"
    );
    assert!(got.width() <= Rational::parse("1e-12").unwrap(), "{got:?}");
}

fn refused(model: &Model, a: PairEntity, b: PairEntity, request: PairRequest, code: &str) {
    let refusal = ask(model, a, b, request).unwrap_err();
    assert_eq!(refusal.code, code, "{refusal:?}");
}

#[test]
fn two_parallel_hole_walls_from_the_real_use_pass() {
    // Through holes of radius 4 at (-10, 0) and radius 3 at (12, 0).
    let (model, body) = build(
        vec![
            block(),
            hole(
                "a",
                "block",
                4,
                p(-10, 0, 20),
                v(0, 0, -1),
                json!("through"),
            ),
            hole("b", "a", 3, p(12, 0, 20), v(0, 0, -1), json!("through")),
        ],
        "b",
    );
    let a = Face(face(&model, body, wall(p(-10, 0, 0), v(0, 0, 1))));
    let b = Face(face(&model, body, wall(p(12, 0, 0), v(0, 0, 1))));
    exact(&model, a, b, MinimumDistance, q(15));
    exact(&model, a, b, CentreDistance, q(22));
    exact(&model, a, b, AngleDegrees, q(0));
}

#[test]
fn offset_parallel_hole_walls_enclose_an_irrational_gap() {
    let (model, body) = build(
        vec![
            block(),
            hole(
                "a",
                "block",
                4,
                p(-10, -5, 20),
                v(0, 0, -1),
                json!("through"),
            ),
            hole("b", "a", 3, p(12, 6, 20), v(0, 0, -1), json!("through")),
        ],
        "b",
    );
    let a = Face(face(&model, body, wall(p(-10, -5, 0), v(0, 0, 1))));
    let b = Face(face(&model, body, wall(p(12, 6, 0), v(0, 0, 1))));
    // Axes 11 sqrt(5) apart.
    root_minus(&model, a, b, q(605), q(7));
}

#[test]
fn skew_hole_walls_meet_along_the_common_perpendicular() {
    // One hole along z at (-10, 0); one along x at y = 10, z = 10.
    let (model, body) = build(
        vec![
            block(),
            hole(
                "a",
                "block",
                4,
                p(-10, 0, 20),
                v(0, 0, -1),
                json!("through"),
            ),
            hole("b", "a", 3, p(30, 10, 10), v(-1, 0, 0), json!("through")),
        ],
        "b",
    );
    let a = Face(face(&model, body, wall(p(-10, 0, 0), v(0, 0, 1))));
    let b = Face(face(&model, body, wall(p(0, 10, 10), v(1, 0, 0))));
    exact(&model, a, b, MinimumDistance, q(3));
    exact(&model, a, b, CentreDistance, q(10));
    exact(&model, a, b, AngleDegrees, q(90));
}

#[test]
fn hole_wall_against_planes_and_edges() {
    let (model, body) = build(
        vec![
            block(),
            hole(
                "a",
                "block",
                4,
                p(-10, 0, 20),
                v(0, 0, -1),
                json!("through"),
            ),
        ],
        "a",
    );
    let bore = Face(face(&model, body, wall(p(-10, 0, 0), v(0, 0, 1))));
    let side = Face(face(&model, body, plane(p(0, -20, 0), v(0, 1, 0))));
    let end = Face(face(&model, body, plane(p(30, 0, 0), v(1, 0, 0))));
    let top = Face(face(&model, body, plane(p(0, 0, 20), v(0, 0, 1))));
    exact(&model, bore, side, MinimumDistance, q(16));
    exact(&model, bore, end, MinimumDistance, q(36));
    // Parallel planar faces keep the planar witness search.
    let bottom = Face(face(&model, body, plane(p(0, 0, 0), v(0, 0, 1))));
    exact(&model, top, bottom, MinimumDistance, q(20));
    // The top face holds the bore's rim.
    exact(&model, bore, top, MinimumDistance, q(0));
    exact(&model, bore, side, CentreDistance, q(20));
    exact(&model, bore, side, AngleDegrees, q(0));
    exact(&model, bore, top, AngleDegrees, q(90));
    refused(&model, bore, top, CentreDistance, "measure.centre-domain");

    // The bottom front edge runs along x at y = -20, z = 0.
    let front = edge(&model, body, |m, e| {
        let e = m.edge(e);
        matches!(m.geom().curve(e.curve), Curve::Line(_))
            && [e.start, e.end].iter().all(|v| {
                v.is_some_and(|v| {
                    let at = rational(&m.geom().position(m.vertex(v).point));
                    at.y == q(-20) && at.z == q(0)
                })
            })
    });
    exact(&model, bore, Edge(front), MinimumDistance, q(16));
    exact(&model, bore, Edge(front), AngleDegrees, q(90));
    // The bore's top rim against the side face and against the edge.
    let rim = edge(&model, body, |m, e| match m.geom().curve(m.edge(e).curve) {
        Curve::Circle(c) => rational(&c.center) == p(-10, 0, 20),
        _ => false,
    });
    exact(&model, Edge(rim), side, MinimumDistance, q(16));
    exact(&model, Edge(rim), side, CentreDistance, q(20));
    root_of(&model, Edge(rim), Edge(front), CentreDistance, q(800), q(0));
    // The rim lies in the top face.
    exact(&model, Edge(rim), top, AngleDegrees, q(0));
    exact(&model, Edge(rim), bore, AngleDegrees, q(90));
}

#[test]
fn blind_holes_whose_nearest_points_lie_on_their_rims() {
    // Radius 4 from the top down to z = 15 and radius 3 from the bottom up to
    // z = 5, 15 apart in x. Planar gap 15 - 7 = 8, axial gap 10.
    let (model, body) = build(
        vec![
            block(),
            hole(
                "a",
                "block",
                4,
                p(-10, 0, 20),
                v(0, 0, -1),
                json!({"blind": "5"}),
            ),
            hole("b", "a", 3, p(5, 0, 0), v(0, 0, 1), json!({"blind": "5"})),
        ],
        "b",
    );
    let a = Face(face(&model, body, wall(p(-10, 0, 0), v(0, 0, 1))));
    let b = Face(face(&model, body, wall(p(5, 0, 0), v(0, 0, 1))));
    root_minus(&model, a, b, q(164), q(0));
    exact(&model, a, b, CentreDistance, q(15));
}

#[test]
fn sphere_cylinder_and_their_caps() {
    // A radius 2 cylinder on z in [0, 10] and two separate balls of radius 3.
    let (model, body) = build(
        vec![
            json!({"id": "post", "op": "cylinder", "radius": "2", "height": "10", "pose": pose(0, 0, 0)}),
            json!({"id": "side", "op": "sphere", "radius": "3", "pose": pose(10, 0, 5)}),
            json!({"id": "over", "op": "sphere", "radius": "1", "pose": pose(10, 0, 15)}),
            union("u1", "post", "side"),
            union("u2", "u1", "over"),
        ],
        "u2",
    );
    let post = Face(face(&model, body, wall(p(0, 0, 0), v(0, 0, 1))));
    let side = Face(face(
        &model,
        body,
        |s| matches!(s, Surface::Sphere(b) if rational(&b.center) == p(10, 0, 5)),
    ));
    let over = Face(face(
        &model,
        body,
        |s| matches!(s, Surface::Sphere(b) if rational(&b.center) == p(10, 0, 15)),
    ));
    let cap = Face(face(&model, body, plane(p(0, 0, 10), v(0, 0, 1))));
    exact(&model, post, side, MinimumDistance, q(5));
    exact(&model, post, side, CentreDistance, q(10));
    refused(&model, post, side, AngleDegrees, "measure.angle-domain");
    // Above the band its rim is nearest: sqrt(8^2 + 5^2) - 1.
    root_minus(&model, post, over, q(89), q(1));
    // The cap's foot misses the disc; its rim is nearest as well.
    root_minus(&model, cap, over, q(89), q(1));
    // Sphere to sphere: centres 10 apart.
    exact(&model, side, over, MinimumDistance, q(6));
    exact(&model, side, over, CentreDistance, q(10));
}

#[test]
fn torus_against_plane_sphere_and_coaxial_cylinder() {
    let (model, body) = build(
        vec![
            json!({"id": "ring", "op": "torus", "major": "5", "minor": "1", "pose": pose(0, 0, 0)}),
            json!({"id": "plate", "op": "block", "size": ["20", "20", "2"], "pose": pose(-10, -10, -5)}),
            json!({"id": "ball", "op": "sphere", "radius": "1", "pose": pose(0, 0, 4)}),
            json!({"id": "pin", "op": "cylinder", "radius": "1", "height": "6", "pose": pose(0, 0, -3)}),
            union("u1", "ring", "plate"),
            union("u2", "u1", "ball"),
            union("u3", "u2", "pin"),
        ],
        "u3",
    );
    let ring = Face(face(&model, body, |s| kind(s) == "torus"));
    let top = Face(face(&model, body, plane(p(0, 0, -3), v(0, 0, 1))));
    let ball = Face(face(&model, body, |s| kind(s) == "sphere"));
    let pin = Face(face(&model, body, wall(p(0, 0, 0), v(0, 0, 1))));
    // The core circle sits 3 above the plate's top.
    exact(&model, ring, top, MinimumDistance, q(2));
    // The ball's centre is on the axis, sqrt(41) from the core circle.
    root_minus(&model, ring, ball, q(41), q(2));
    // The pin runs through the hole: 5 - 1 - 1.
    exact(&model, ring, pin, MinimumDistance, q(3));
    exact(&model, ring, pin, CentreDistance, q(0));
    exact(&model, ring, top, AngleDegrees, q(0));
}

#[test]
fn stacked_coaxial_tori() {
    let (model, body) = build(
        vec![
            json!({"id": "low", "op": "torus", "major": "5", "minor": "1", "pose": pose(0, 0, 0)}),
            json!({"id": "high", "op": "torus", "major": "4", "minor": "1", "pose": pose(0, 0, 5)}),
            union("u", "low", "high"),
        ],
        "u",
    );
    let torus = |z: i64| {
        Face(face(
            &model,
            body,
            move |s| matches!(s, Surface::Torus(t) if rational(&t.center) == p(0, 0, z)),
        ))
    };
    // Core circles sqrt(5^2 + 1^2) apart at every azimuth.
    root_of(&model, torus(0), torus(5), MinimumDistance, q(26), q(2));
    exact(&model, torus(0), torus(5), CentreDistance, q(5));
    exact(&model, torus(0), torus(5), AngleDegrees, q(0));
}

#[test]
fn frustum_wall_against_plane_and_spheres() {
    // Base radius 3 at z = 0, top radius 1 at z = 4: apex at z = 6, slope 1/2.
    let (model, body) = build(
        vec![
            json!({"id": "cone", "op": "cone", "base_radius": "3", "top_radius": "1", "height": "4",
                   "pose": pose(0, 0, 0)}),
            json!({"id": "roof", "op": "block", "size": ["20", "20", "2"], "pose": pose(-10, -10, 7)}),
            json!({"id": "near", "op": "sphere", "radius": "1", "pose": pose(6, 0, 4)}),
            json!({"id": "low", "op": "sphere", "radius": "1", "pose": pose(-10, 0, 2)}),
            union("u1", "cone", "roof"),
            union("u2", "u1", "near"),
            union("u3", "u2", "low"),
        ],
        "u3",
    );
    let cone = Face(face(&model, body, |s| kind(s) == "cone"));
    let roof = Face(face(&model, body, plane(p(0, 0, 7), v(0, 0, 1))));
    let near = Face(face(
        &model,
        body,
        |s| matches!(s, Surface::Sphere(b) if rational(&b.center) == p(6, 0, 4)),
    ));
    let low = Face(face(
        &model,
        body,
        |s| matches!(s, Surface::Sphere(b) if rational(&b.center) == p(-10, 0, 2)),
    ));
    // The top rim is 3 below the roof.
    exact(&model, cone, roof, MinimumDistance, q(3));
    // The foot (2, 0, 2) lies on the wall: 2 sqrt(5) - 1.
    root_minus(&model, cone, near, q(20), q(1));
    // The foot falls below the base, so the base rim is nearest: sqrt(53) - 1.
    root_minus(&model, cone, low, q(53), q(1));
    exact(&model, cone, roof, AngleDegrees, q(90));
    exact(&model, cone, near, CentreDistance, q(6));
}

#[test]
fn trimmed_wall_breaking_out_of_a_side_is_measured_on_its_trim() {
    // A radius 4 hole at y = -18 breaks out of the y = -20 face. A ball of
    // radius 1 at (0, -30, 10) faces the removed side; the nearest remaining
    // wall points are the breakout lines at x = +-2 sqrt(3), y = -20.
    let (model, body) = build(
        vec![
            block(),
            hole(
                "a",
                "block",
                4,
                p(0, -18, 20),
                v(0, 0, -1),
                json!("through"),
            ),
            json!({"id": "ball", "op": "sphere", "radius": "1", "pose": pose(0, -30, 10)}),
            union("u", "a", "ball"),
        ],
        "u",
    );
    let bore = Face(face(&model, body, wall(p(0, -18, 0), v(0, 0, 1))));
    let ball = Face(face(&model, body, |s| kind(s) == "sphere"));
    let back = Face(face(&model, body, plane(p(0, 20, 0), v(0, 1, 0))));
    root_minus(&model, bore, ball, q(112), q(1));
    // The retained far side faces the back: 20 - (-18) - 4.
    exact(&model, bore, back, MinimumDistance, q(34));
}

/// Bare finite line edges in a fresh model.
fn segments(pairs: &[[RPoint3; 2]]) -> (Model, Vec<EdgeId>) {
    use catallax_shapes::{Line, PointDef};
    use catallax_topology::{Edge as TopoEdge, Sense, Vertex};
    let mut model = Model::new();
    let mut edges = Vec::new();
    for [a, b] in pairs {
        let line = Line::new(a.clone().into(), b.minus(a).into(), &mut Decider::default()).unwrap();
        let curve = model.geom_mut().add_curve(Curve::Line(line));
        let [a, b] = [a, b].map(|q| {
            let point = model
                .geom_mut()
                .add_point(PointDef::Explicit(q.clone().into()));
            model.add_vertex(Vertex { point })
        });
        edges.push(model.add_edge(TopoEdge {
            curve,
            sense: Sense::Same,
            start: Some(a),
            end: Some(b),
        }));
    }
    (model, edges)
}

#[test]
fn straight_edges_answer_exact_oblique_angles_and_clamped_gaps() {
    // Oblique holes do not build yet, so bare edges carry the exact table.
    let (model, e) = segments(&[
        [p(0, 0, 0), p(2, 0, 0)],
        [p(5, 0, 0), p(6, 1, 0)],
        [p(0, 5, 0), p(1, 6, 0)],
        [p(0, 9, 0), p(1, 9, 1)],
        [p(1, 1, 1), p(1, 3, 1)],
        [p(5, 4, 3), p(5, 5, 3)],
    ]);
    exact(&model, Edge(e[0]), Edge(e[1]), AngleDegrees, q(45));
    exact(&model, Edge(e[2]), Edge(e[3]), AngleDegrees, q(60));
    // Skew segments meet at an interior pair sqrt(2) apart; the clamped
    // pair joins (2, 0, 0) to (5, 4, 3).
    root_of(&model, Edge(e[0]), Edge(e[4]), MinimumDistance, q(2), q(0));
    root_of(&model, Edge(e[0]), Edge(e[5]), MinimumDistance, q(34), q(0));
    for (a, b, square) in [(0, 4, 2), (0, 5, 34)] {
        let native = crate::edge_pair(
            &model,
            e[a],
            e[b],
            crate::EdgePairRequest::MinimumDistance,
            &Accuracy::default(),
            &mut Decider::default(),
        )
        .unwrap()
        .value;
        assert!(&native.lo * &native.lo <= q(square) && q(square) <= &native.hi * &native.hi);
    }
}

#[test]
fn a_rotated_bar_has_an_enclosed_angle() {
    // The rational quaternion (2, 0, 1, 0) turns z to (4/5, 0, 3/5).
    let (model, body) = build(
        vec![
            json!({"id": "post", "op": "cylinder", "radius": "1", "height": "10", "pose": pose(0, 0, 0)}),
            json!({"id": "bar", "op": "cylinder", "radius": "1", "height": "10",
                   "pose": {"rotation": ["2", "0", "1", "0"], "translation": ["10", "0", "0"]}}),
            union("u", "post", "bar"),
        ],
        "u",
    );
    let post = Face(face(&model, body, wall(p(0, 0, 0), v(0, 0, 1))));
    let bar = Face(face(&model, body, wall(p(10, 0, 0), v(4, 0, 3))));
    // acos(3/5) in degrees = 53.1301023541559787031443874409065893...
    let got = ask(&model, post, bar, AngleDegrees).unwrap();
    assert!(
        got.lo <= Rational::parse("53.130102354155978703144387440906590").unwrap()
            && Rational::parse("53.130102354155978703144387440906589").unwrap() <= got.hi
            && got.width() <= Rational::parse("1e-12").unwrap(),
        "{got:?}"
    );
}

#[test]
fn uncovered_pairs_refuse_with_a_type() {
    // A cone and a cylinder on skew axes have no closed form here.
    let (model, body) = build(
        vec![
            json!({"id": "cone", "op": "cone", "base_radius": "3", "top_radius": "1", "height": "4",
                   "pose": pose(0, 0, 0)}),
            json!({"id": "bar", "op": "cylinder", "radius": "1", "height": "20", "pose": turned(-10, 10, 2)}),
            union("u", "cone", "bar"),
        ],
        "u",
    );
    let cone = Face(face(&model, body, |s| kind(s) == "cone"));
    let bar = Face(face(&model, body, wall(p(-10, 0, 2), v(0, 1, 0))));
    refused(
        &model,
        cone,
        bar,
        MinimumDistance,
        "measure.distance-uncovered",
    );
    exact(&model, cone, bar, AngleDegrees, q(90));
}
