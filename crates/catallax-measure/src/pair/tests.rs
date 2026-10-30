//! V2-K contract: the V2 batch 2 designer's K1-K3 oracles, fixed before probes.
//! Physical aliases name entities by their exact endpoints, never by handle order.

use super::*;
use crate::{distance::PreparedFace, edge_pair};
use catallax_numbers::{BudgetPolicy, OperationWork, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{
    Line, Point3, PointDef, Pose, RVector3, Rotation, Surface, distance::Options,
};
use catallax_topology::{BodyId, Edge, FaceId, Sense, Vertex};
use serde_json::json;

const WORK: u64 = 20_000_000;
const ALL: [EdgePairRequest; 3] = [
    EdgePairRequest::MinimumDistance,
    EdgePairRequest::MidpointDistance,
    EdgePairRequest::DirectionAngleDegrees,
];

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn p(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(r(x), r(y), r(z))
}
fn accuracy() -> Accuracy {
    Accuracy {
        relative: Rational::zero(),
        absolute: Rational::parse("1e-8").unwrap(),
    }
}

/// The batch 1 prism `[-20,20] x [-15,15] x [0,depth]` with its radius 6 bore
/// about Z, optionally after a radius 2 fillet of E0.
fn part(depth: i64, fillet: bool) -> (Model, BodyId) {
    let pose = |x: i64, y: i64, z: i64| {
        json!({"rotation": ["1", "0", "0", "0"],
               "translation": [x.to_string(), y.to_string(), z.to_string()]})
    };
    let mut steps = vec![
        json!({"id": "block", "op": "block", "size": ["40", "30", depth.to_string()],
               "pose": pose(-20, -15, 0)}),
        json!({"id": "bore", "op": "cylinder", "radius": "6", "height": "100",
               "pose": pose(0, 0, depth - 100)}),
        json!({"id": "cut", "op": "boolean", "kind": "subtract", "target": "block",
               "tool": "bore"}),
    ];
    if fillet {
        steps.push(
            json!({"id": "round", "op": "fillet", "body": "cut", "radius": "2",
            "edges": [{"between": [{"step": "block", "role": "start"},
                {"step": "block", "role": {"side": {"segment": "y_min"}}}]}]}),
        );
    }
    let result = if fillet { "round" } else { "cut" };
    let recipe = json!({"format": "catallax-recipe/3", "name": "V2 batch 2 block",
                        "steps": steps, "result": result});
    let recipe = Recipe::parse(&recipe.to_string()).unwrap();
    let body =
        catallax_app::build_incremental_without_display(&recipe, None, &mut Decider::default())
            .unwrap();
    let built = body.checkpoint.built();
    (built.model.clone(), built.result)
}

/// E0-E3 at one depth, by their exact physical endpoints.
fn aliases(depth: i64) -> [[RPoint3; 2]; 4] {
    [
        [p(-20, -15, 0), p(20, -15, 0)],
        [p(-20, 15, depth), p(20, 15, depth)],
        [p(20, -15, depth), p(20, 15, depth)],
        [p(20, -15, 0), p(20, -15, depth)],
    ]
}

fn ends(model: &Model, edge: EdgeId) -> Option<[RPoint3; 2]> {
    let e = model.edge(edge);
    let [a, b] =
        [e.start?, e.end?].map(|v| model.geom().position(model.vertex(v).point).as_rational());
    Some([a?, b?])
}

/// The unique straight body edge joining both points, in either order.
fn edge(model: &Model, body: BodyId, [a, b]: &[RPoint3; 2]) -> EdgeId {
    let found: Vec<_> = model
        .edges_of_body(body)
        .into_iter()
        .filter(|e| matches!(model.geom().curve(model.edge(*e).curve), Curve::Line(_)))
        .filter(|e| ends(model, *e).is_some_and(|[s, t]| (&s, &t) == (a, b) || (&s, &t) == (b, a)))
        .collect();
    assert_eq!(found.len(), 1, "physical edge {a:?}..{b:?}");
    found[0]
}

/// The planar face lying wholly at height `z`.
fn cap(model: &Model, body: BodyId, z: i64) -> FaceId {
    let found: Vec<_> = model
        .faces_of_body(body)
        .into_iter()
        .filter(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::Plane(_)
            )
        })
        .filter(|f| {
            let bounds = crate::face_bounding_box(model, *f, &accuracy(), &mut Decider::default())
                .unwrap()
                .value;
            bounds.min[2] == r(z) && bounds.max[2] == r(z)
        })
        .collect();
    assert_eq!(found.len(), 1, "planar cap at z={z}");
    found[0]
}

/// Bare finite line edges in a fresh model. A reversed edge keeps the carrier
/// and travels it from its second end to its first.
fn segments(pairs: &[[RPoint3; 2]], reversed: bool) -> (Model, Vec<EdgeId>) {
    let mut model = Model::new();
    let mut edges = Vec::new();
    for [a, b] in pairs {
        let direction = if a == b {
            RVector3::new(r(1), r(0), r(0))
        } else {
            b.minus(a)
        };
        let line = Line::new(a.clone().into(), direction.into(), &mut Decider::default()).unwrap();
        let curve = model.geom_mut().add_curve(Curve::Line(line));
        let [a, b] = [a, b].map(|q| {
            let point = model
                .geom_mut()
                .add_point(PointDef::Explicit(q.clone().into()));
            model.add_vertex(Vertex { point })
        });
        let (sense, start, end) = if reversed {
            (Sense::Reversed, b, a)
        } else {
            (Sense::Same, a, b)
        };
        edges.push(model.add_edge(Edge {
            curve,
            sense,
            start: Some(start),
            end: Some(end),
        }));
    }
    (model, edges)
}

fn shifted([a, b]: &[RPoint3; 2], by: &RVector3) -> [RPoint3; 2] {
    [a.plus(by), b.plus(by)]
}

/// One request inside the application operation allowance.
fn ask(model: &Model, a: EdgeId, b: EdgeId, request: EdgePairRequest) -> Answered<Enclosure> {
    let (answer, work) = with_operation_budget(WORK, || {
        edge_pair(model, a, b, request, &accuracy(), &mut Decider::default())
    });
    assert!(!work.exhausted, "{request:?} exhausted {WORK}: {work:?}");
    answer
}

/// An unmarked answer's enclosure.
fn value(model: &Model, a: EdgeId, b: EdgeId, request: EdgePairRequest) -> Enclosure {
    let answer = ask(model, a, b, request).unwrap();
    assert!(answer.marks.is_empty(), "{answer:?}");
    answer.value
}

fn code(answer: Answered<Enclosure>) -> &'static str {
    answer.unwrap_err().code
}

#[derive(Clone, Copy, Debug)]
enum Expected {
    /// An exact rational value: `lo == hi`.
    Exact(i64),
    /// The irrational root of this nonsquare integer.
    Root(i64),
}

/// Exact rational containment, without decimal formatting or floating point.
fn check(answer: &Answer<Enclosure>, expected: Expected) {
    let value = &answer.value;
    assert!(answer.marks.is_empty(), "{answer:?}");
    assert!(value.lo.sign() != Sign::Negative, "{value:?}");
    assert!(value.width() <= accuracy().absolute, "{value:?}");
    match expected {
        Expected::Exact(n) => assert_eq!(value, &Enclosure::exact(r(n))),
        Expected::Root(n) => {
            assert!(!value.is_exact(), "sqrt({n}) is irrational: {value:?}");
            assert!(
                &value.lo * &value.lo <= r(n) && r(n) <= &value.hi * &value.hi,
                "{value:?} must enclose sqrt({n})"
            );
        }
    }
}

#[test]
fn finite_edge_pairs_and_points_have_proven_bounds() {
    use EdgePairRequest::{DirectionAngleDegrees as Angle, MidpointDistance, MinimumDistance};
    let offset = RVector3::new(
        Rational::from_ratio(1, 3).unwrap(),
        Rational::from_ratio(2, 3).unwrap(),
        Rational::from_ratio(7, 3).unwrap(),
    );
    let mut most = 0;
    for depth in [10, 15] {
        let (model, body) = part(depth, false);
        let physical = aliases(depth);
        let [e0, e1, e2, e3] = physical.each_ref().map(|ends| edge(&model, body, ends));
        // M0 and M2 are native half-arc-length points, exact on rational lines.
        for (edge, want) in [(e0, p(0, -15, 0)), (e2, p(20, 0, depth))] {
            let midpoint =
                crate::edge_midpoint(&model, edge, &accuracy(), &mut Decider::default()).unwrap();
            assert!(midpoint.marks.is_empty());
            assert_eq!(
                midpoint.value,
                [want.x, want.y, want.z].map(Enclosure::exact)
            );
        }
        let d2 = depth * depth;
        let rows = [
            (0, 1, MinimumDistance, Expected::Root(900 + d2)),
            (0, 2, MinimumDistance, Expected::Exact(depth)),
            (0, 3, MinimumDistance, Expected::Exact(0)),
            (0, 2, MidpointDistance, Expected::Root(625 + d2)),
            (0, 1, Angle, Expected::Exact(0)),
            (0, 2, Angle, Expected::Exact(90)),
            (0, 0, MinimumDistance, Expected::Exact(0)),
            (0, 0, MidpointDistance, Expected::Exact(0)),
            (0, 0, Angle, Expected::Exact(0)),
        ];
        let handles = [e0, e1, e2, e3];
        for (i, j, request, expected) in rows {
            let (a, b) = (handles[i], handles[j]);
            let (answer, work) = with_operation_budget(WORK, || {
                edge_pair(&model, a, b, request, &accuracy(), &mut Decider::default())
            });
            assert!(!work.exhausted);
            most = most.max(work.used);
            let answer = answer.unwrap();
            check(&answer, expected);
            eprintln!(
                "D{depth} E{i}/E{j} {request:?}: [{}, {}] work {}",
                answer.value.lo, answer.value.hi, work.used
            );
            assert_eq!(value(&model, b, a, request), answer.value, "pair order");
            // The same physical segments, reversed and moved off the origin.
            let pair = [physical[i].clone(), physical[j].clone()];
            for shift in [RVector3::zero(), offset.clone()] {
                for reversed in [false, true] {
                    let moved = pair.each_ref().map(|ends| shifted(ends, &shift));
                    let (bare, ids) = segments(&moved, reversed);
                    for (a, b) in [(ids[0], ids[1]), (ids[1], ids[0])] {
                        assert_eq!(
                            value(&bare, a, b, request),
                            answer.value,
                            "D{depth} E{i}/E{j} {request:?} reversed={reversed} shift={shift:?}"
                        );
                    }
                }
            }
        }
    }
    // Carriers meet at (2,0,0), outside both segments. The finite minimum is at
    // the clamped ends (1,0,0) and (2,1,0), never the carrier crossing.
    let neighbour = [[p(0, 0, 0), p(1, 0, 0)], [p(2, 1, 0), p(2, 2, 0)]];
    for shift in [RVector3::zero(), offset] {
        for reversed in [false, true] {
            let moved = neighbour.each_ref().map(|ends| shifted(ends, &shift));
            let (bare, ids) = segments(&moved, reversed);
            let answer = ask(&bare, ids[0], ids[1], MinimumDistance).unwrap();
            check(&answer, Expected::Root(2));
            assert!(answer.value.lo > r(1), "{answer:?}");
            assert_eq!(value(&bare, ids[1], ids[0], MinimumDistance), answer.value);
        }
    }
    eprintln!("most edge pair work: {most}");
}

#[test]
fn query_domains_and_work_are_explicit() {
    use EdgePairRequest::{DirectionAngleDegrees as Angle, MinimumDistance};
    let (model, body) = part(10, false);
    let [e0, e1, ..] = aliases(10).each_ref().map(|ends| edge(&model, body, ends));
    let reference = ask(&model, e0, e1, MinimumDistance).unwrap();
    check(&reference, Expected::Root(1000));
    let still_answers = || assert_eq!(value(&model, e0, e1, MinimumDistance), reference.value);

    // E0 against the bottom bore rim, a closed circle with no unique midpoint.
    let rim: Vec<_> = model
        .edges_of_body(body)
        .into_iter()
        .filter(|e| match model.geom().curve(model.edge(*e).curve) {
            Curve::Circle(c) => c.center.as_rational() == Some(p(0, 0, 0)),
            _ => false,
        })
        .collect();
    let [rim] = rim.as_slice() else {
        panic!("one bottom bore rim: {rim:?}");
    };
    for request in ALL {
        for (a, b) in [(e0, *rim), (*rim, e0)] {
            assert_eq!(
                code(ask(&model, a, b, request)),
                "measure.unsupported-curve"
            );
            still_answers();
        }
    }

    // Direction (1,1,0) against X is oblique: the angle refuses, distance does not.
    let (bare, ids) = segments(
        &[[p(0, 0, 0), p(1, 1, 0)], [p(0, -1, 0), p(1, -1, 0)]],
        false,
    );
    assert_eq!(
        code(ask(&bare, ids[0], ids[1], Angle)),
        "measure.angle-domain"
    );
    check(
        &ask(&bare, ids[0], ids[1], MinimumDistance).unwrap(),
        Expected::Exact(1),
    );
    still_answers();

    // A zero-length segment and an irrational end refuse every request.
    let (bare, ids) = segments(&[[p(3, 4, 5), p(3, 4, 5)], [p(0, 0, 0), p(1, 0, 0)]], false);
    let mut irrational = Model::new();
    let root = Point3::new(Real::from(2).sqrt(), Real::from(0), Real::from(0));
    let line = Line::new(
        root.clone(),
        RVector3::new(r(0), r(1), r(0)).into(),
        &mut Decider::default(),
    )
    .unwrap();
    let curve = irrational.geom_mut().add_curve(Curve::Line(line));
    let [a, b] = [
        root.clone(),
        root.plus(&RVector3::new(r(0), r(1), r(0)).into()),
    ]
    .map(|q| {
        let point = irrational.geom_mut().add_point(PointDef::Explicit(q));
        irrational.add_vertex(Vertex { point })
    });
    let tilted = irrational.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: Some(a),
        end: Some(b),
    });
    for request in ALL {
        for (a, b) in [(ids[0], ids[1]), (ids[1], ids[0])] {
            assert_eq!(code(ask(&bare, a, b, request)), "measure.degenerate-edge");
        }
        assert_eq!(
            code(ask(&irrational, tilted, tilted, request)),
            "measure.irrational-endpoint"
        );
        still_answers();
    }

    // Face pairs outside the parallel planar domain refuse before any search.
    let mut faces = model.clone();
    let f0 = cap(&faces, body, 0);
    let wall = faces
        .faces_of_body(body)
        .into_iter()
        .find(|f| {
            matches!(
                faces.geom().surface(faces.face(*f).surface),
                Surface::Cylinder(_)
            )
        })
        .unwrap();
    let front = faces
        .faces_of_body(body)
        .into_iter()
        .find(|f| {
            let bounds = crate::face_bounding_box(&faces, *f, &accuracy(), &mut Decider::default())
                .unwrap()
                .value;
            bounds.min[1] == r(-15) && bounds.max[1] == r(-15)
        })
        .unwrap();
    let options = Options {
        width: accuracy().absolute,
        max_cells: 4096,
    };
    for (a, b, want) in [
        (wall, f0, "measure.nonplanar-face"),
        (f0, wall, "measure.nonplanar-face"),
        (f0, front, "measure.nonparallel-faces"),
    ] {
        let (answer, work) = with_operation_budget(WORK, || {
            let mut dec = Decider::default();
            let a = PreparedFace::new(&faces, a, &mut dec)?;
            let b = PreparedFace::new(&faces, b, &mut dec)?;
            a.parallel_face_to_face(&b, &mut faces, &options, &mut dec)
        });
        assert!(!work.exhausted);
        assert_eq!(answer.unwrap_err().code, want);
        still_answers();
    }

    // An unknown handle is a programming error: it panics, never answers.
    let unknown = model.edges_of_body(body).into_iter().max().unwrap();
    let small = segments(&[[p(0, 0, 0), p(1, 0, 0)]], false);
    assert!(unknown.index() >= 1);
    let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        edge_pair(
            &small.0,
            small.1[0],
            unknown,
            MinimumDistance,
            &accuracy(),
            &mut Decider::default(),
        )
    }));
    assert!(panicked.is_err());
    still_answers();

    // Real numerical work under zero allowance refuses with its native code.
    for request in ALL {
        let (answer, work): (_, OperationWork) = with_operation_budget(0, || {
            edge_pair(
                &model,
                e0,
                e1,
                request,
                &accuracy(),
                &mut Decider::default(),
            )
        });
        assert!(work.exhausted);
        assert_eq!(code(answer), "numbers.operation-budget");
        still_answers();
    }
    let mut idle = Decider::new(BudgetPolicy {
        max_work: 0,
        ..BudgetPolicy::default()
    });
    let answer = edge_pair(&model, e0, e1, MinimumDistance, &accuracy(), &mut idle);
    assert_eq!(code(answer), "numbers.budget-exhausted");
    still_answers();

    // A width above the request is marked; the outward root stays enclosed.
    let mut coarse = Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        ..BudgetPolicy::default()
    });
    let tight = Accuracy {
        relative: Rational::zero(),
        absolute: Rational::parse("1e-30").unwrap(),
    };
    let answer = edge_pair(&model, e0, e1, MinimumDistance, &tight, &mut coarse).unwrap();
    let value = &answer.value;
    assert!(&value.lo * &value.lo <= r(1000) && r(1000) <= &value.hi * &value.hi);
    assert!(value.width() > tight.absolute, "{value:?}");
    let codes: Vec<_> = answer.marks.iter().map(|m| m.code).collect();
    assert_eq!(codes, ["measure.accuracy-not-reached"]);
    still_answers();
}

#[test]
fn parallel_face_distance_respects_trims_and_completion() {
    let options = Options {
        width: Rational::parse("1e-8").unwrap(),
        max_cells: 4096,
    };
    let distance = |model: &mut Model, a: FaceId, b: FaceId, options: &Options, limit: u64| {
        with_operation_budget(limit, || {
            let mut dec = Decider::default();
            let a = PreparedFace::new(model, a, &mut dec)?;
            let b = PreparedFace::new(model, b, &mut dec)?;
            a.parallel_face_to_face(&b, model, options, &mut dec)
        })
    };
    // Both caps keep (10,0) beside the bore, so the supporting gap is attained.
    // Neither cap centre belongs to its face.
    for (depth, fillet) in [(10, false), (15, false), (15, true)] {
        let (mut model, body) = part(depth, fillet);
        let caps = [cap(&model, body, 0), cap(&model, body, depth)];
        for [a, b] in [caps, [caps[1], caps[0]]] {
            let (bounds, work) = distance(&mut model, a, b, &options, WORK);
            assert!(!work.exhausted, "{work:?}");
            let bounds = bounds.unwrap();
            eprintln!(
                "D{depth} fillet={fillet}: lower {} upper {:?} converged {} cells {} work {}",
                bounds.lower, bounds.upper, bounds.converged, bounds.cells, work.used
            );
            assert!(bounds.converged, "{bounds:?}");
            assert_eq!(bounds.lower, r(depth));
            assert_eq!(bounds.upper, Some(r(depth)));
            let [s, t] = bounds.witnesses.unwrap();
            let mut dec = Decider::default();
            assert!(catallax_boolean::face_contains(&mut model, a, &s, &mut dec).unwrap());
            assert!(catallax_boolean::face_contains(&mut model, b, &t, &mut dec).unwrap());
        }
    }

    // z=0 over [0,1]^2 and z=3 over [2,3]x[0,1]: the trimmed minimum is
    // sqrt(10), never the supporting-plane gap 3.
    let mut model = Model::new();
    let mut dec = Decider::default();
    let mut patch = |x: i64, z: i64| {
        let pose = Pose::new(
            Rotation::identity(),
            RVector3::new(r(x), r(0), r(z)).into(),
            &mut dec,
        )
        .unwrap();
        catallax_topology::block(&mut model, [r(1), r(1), r(1)], &pose).unwrap()
    };
    let (low, high) = (patch(0, -1), patch(2, 3));
    let faces = [cap(&model, low, 0), cap(&model, high, 3)];
    let (bounds, work) = distance(&mut model, faces[0], faces[1], &options, WORK);
    assert!(!work.exhausted, "{work:?}");
    let bounds = bounds.unwrap();
    eprintln!(
        "disjoint patches: lower {} upper {:?} converged {} cells {} work {}",
        bounds.lower, bounds.upper, bounds.converged, bounds.cells, work.used
    );
    let upper = bounds.upper.clone().unwrap();
    assert!(bounds.converged, "{bounds:?}");
    assert!(&bounds.lower * &bounds.lower <= r(10) && r(10) <= &upper * &upper);
    assert!(bounds.lower > r(3), "{bounds:?}");
    assert!(&upper - &bounds.lower <= options.width, "{bounds:?}");

    // A stopped search keeps its incompleteness flag; no allowance, no answer.
    let starved = Options {
        width: options.width.clone(),
        max_cells: 1,
    };
    let (bounds, work) = distance(&mut model, faces[0], faces[1], &starved, WORK);
    assert!(!work.exhausted, "{work:?}");
    assert!(!bounds.unwrap().converged);
    let (bounds, work) = distance(&mut model, faces[0], faces[1], &options, 0);
    assert!(work.exhausted);
    assert_eq!(bounds.unwrap_err().code, "numbers.operation-budget");
}
