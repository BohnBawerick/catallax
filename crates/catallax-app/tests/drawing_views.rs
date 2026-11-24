//! Contour generators on release parts: the P1 bracket, the P3 filleted boss
//! and an elliptic extrusion; on spline faces: an extruded Bezier arch,
//! rational patch roofs and a lofted barrel; and on faces trimmed by
//! intersection branches: a ball with a cylindrical hole.

use std::num::NonZeroUsize;

use catallax_app::DEFAULT_BUILD_WORK;
use catallax_app::views::{
    CONTOUR_FIT_MARK, ContourArc, Contours, DrawingCurve, contour_generators, drawing_curve,
};
use catallax_boolean::{BoolOp, boolean};
use catallax_features::evaluate;
use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, Rational, Real, set_threads, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{Curve, PointId, Pose, RVector3, Rotation, Surface};
use catallax_topology::{BodyId, Model, cylinder, sphere};

#[allow(dead_code)]
mod rational_patch_cases;
use rational_patch_cases::{Roof, solid};

const TARGET: u64 = 20_000_000;

fn vector(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(x.into(), y.into(), z.into())
}

fn zero(value: &Real) -> bool {
    matches!(
        Decider::default().sign(value),
        catallax_numbers::Outcome::ProvenCoincident(_)
    )
}

/// Arc count, line count, arcs running along an edge, and work.
fn summary(model: &Model, body: BodyId, d: &RVector3) -> (usize, usize, usize, u64, Contours) {
    let (answer, work) = with_operation_budget(DEFAULT_BUILD_WORK, || {
        contour_generators(model, body, d, &mut Decider::default())
    });
    let contours = answer.unwrap().value;
    assert!(!work.exhausted);
    let lines = contours
        .arcs
        .iter()
        .filter(|arc| {
            matches!(
                curve_form(contours.model.geom(), arc.curve, &mut Decider::default()).unwrap(),
                CurveForm::Line { .. }
            )
        })
        .count();
    let on_edges = contours
        .arcs
        .iter()
        .filter(|a| a.on_faces.len() > 1)
        .count();
    (contours.arcs.len(), lines, on_edges, work.used, contours)
}

#[test]
fn p1_bracket_generators_are_exact_cylinder_lines_inside_the_target() {
    let recipe = Recipe::parse(include_str!("../../../parts/P1-bracket.json")).unwrap();
    let (built, _) = with_operation_budget(DEFAULT_BUILD_WORK, || {
        evaluate(&recipe, &mut Decider::default())
    });
    let built = built.unwrap();
    // Front: boss and bore share the plane y = 25, the bolt pairs y = 12 and
    // y = 32; the slot ends are seen along their axes. Top: one line on each
    // slot end at z = 18. Right: as front, plus the slot ends' generators on
    // their tangent edges with the slot's flat walls. Isometric: five planes
    // for the vertical holes and two for the slot ends.
    for (name, d, arcs, on_edges) in [
        ("front", vector(0, -1, 0), 12, 0),
        ("top", vector(0, 0, 1), 2, 0),
        ("right", vector(1, 0, 0), 16, 4),
        ("iso", vector(1, -1, 1), 14, 0),
    ] {
        let mut reference = None;
        for threads in [1, 2, 4] {
            set_threads(NonZeroUsize::new(threads));
            let (count, lines, edges, work, contours) = summary(&built.model, built.result, &d);
            eprintln!(
                "P1 {name} at {threads} threads: {count} arcs, {edges} on edges, {work} work"
            );
            assert_eq!((count, lines, edges), (arcs, arcs, on_edges), "{name}");
            assert!(work < TARGET, "{name}: {work} work units");
            for arc in &contours.arcs {
                let Surface::Cylinder(c) = built
                    .model
                    .geom()
                    .surface(built.model.face(arc.face).surface)
                else {
                    panic!("{name}: generator off a cylinder face");
                };
                // Each end lies on the cylinder with its radial normal
                // perpendicular to the view direction.
                for end in [arc.start.unwrap(), arc.end.unwrap()] {
                    let r = contours.model.geom().position(end).minus(&c.origin);
                    let axis2 = c.axis.norm_squared();
                    let radial = r.scaled(&axis2).minus(&c.axis.scaled(&r.dot(&c.axis)));
                    assert!(zero(&radial.dot(&d.to_real())), "{name}");
                    let off_axis = &(&r.norm_squared() * &axis2) - &r.dot(&c.axis).square();
                    assert!(zero(&(&off_axis - &(&c.radius_sq * &axis2))), "{name}");
                }
            }
            match &reference {
                None => reference = Some(work),
                Some(first) => assert_eq!(*first, work, "{name} at {threads} threads"),
            }
        }
        set_threads(None);
    }
}

#[test]
fn elliptic_extrusion_generators_touch_the_ellipse_where_its_normal_turns_away() {
    let source = serde_json::json!({"format":"soliddraft-recipe/1", "name":"ellipse", "steps":[
        {"id":"profile","op":"profile","plane":"XY","loops":{"outer":[
            {"type":"ellipse","center":[1,2],"axis":[4,0],"ratio":0.5}
        ]}},
        {"id":"body","op":"extrude","profile":"profile","distance":2.0}
    ], "result":"body"});
    let app = catallax_app::build_soliddraft(
        &source.to_string(),
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    let (model, body) = (&app.built.model, app.built.result);
    let (count, lines, _, _, front) = summary(model, body, &vector(0, -1, 0));
    assert_eq!((count, lines), (2, 2));
    for arc in &front.arcs {
        for end in [arc.start.unwrap(), arc.end.unwrap()] {
            let p = front.model.geom().position(end);
            assert!(zero(&(&p.y - &Real::from(2))));
            assert!(zero(&(&(&p.x - &Real::from(1)).square() - &Real::from(16))));
        }
    }
    assert!(summary(model, body, &vector(0, 0, 1)).4.arcs.is_empty());
    // The ellipse (x - 1)^2 / 16 + (y - 2)^2 / 4 = 1 has normal
    // ((x - 1) / 16, (y - 2) / 4); seen along (1, 1, 0) the generators satisfy
    // (x - 1) / 16 + (y - 2) / 4 = 0.
    let (count, _, _, _, oblique) = summary(model, body, &vector(1, 1, 0));
    assert_eq!(count, 2);
    for arc in &oblique.arcs {
        for end in [arc.start.unwrap(), arc.end.unwrap()] {
            let p = oblique.model.geom().position(end);
            let along = &(&(&p.x - &Real::from(1)) / &Real::from(16))
                + &(&(&p.y - &Real::from(2)) / &Real::from(4));
            assert!(zero(&along));
        }
    }
}

#[test]
fn p3_fillet_torus_generators_are_exact_along_and_across_its_axis_and_fitted_otherwise() {
    let recipe = Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    let (built, _) = with_operation_budget(DEFAULT_BUILD_WORK, || {
        evaluate(&recipe, &mut Decider::default())
    });
    let built = built.unwrap();
    let model = &built.model;
    let is_torus = |face| {
        matches!(
            model.geom().surface(model.face(face).surface),
            Surface::Torus(_)
        )
    };
    // Expected (arcs, torus arcs, fitted arcs, arcs on edges). The junction
    // fillet is the inner lower quarter of a torus about the boss axis.
    // Front: the two quarter meridians and its contact circle with the block
    // top, plus the boss lines and the two r1 fillets along x. Top: the
    // contact circle with the boss, plus the four r1 fillets' outer lines.
    // Isometric: one fitted arc between the two ends where the branch meets
    // the contact circle with the boss.
    for (name, d, expected) in [
        ("front", vector(0, -1, 0), (7, 3, 0, 3)),
        ("top", vector(0, 0, 1), (5, 1, 0, 5)),
        ("iso", vector(1, -1, 1), (5, 1, 1, 0)),
    ] {
        let mut reference = None;
        for threads in [1, 2, 4] {
            set_threads(NonZeroUsize::new(threads));
            let (answer, work) = with_operation_budget(DEFAULT_BUILD_WORK, || {
                contour_generators(model, built.result, &d, &mut Decider::default())
            });
            assert!(!work.exhausted);
            let answer = answer.unwrap();
            let contours = &answer.value;
            let torus = contours.arcs.iter().filter(|a| is_torus(a.face)).count();
            let fits: Vec<_> = contours.arcs.iter().filter_map(|a| a.fit.clone()).collect();
            for fit in &fits {
                assert!(fit.deviation_mm <= fit.target_mm, "{name}");
            }
            assert_eq!(
                answer.marks.iter().any(|m| m.code == CONTOUR_FIT_MARK),
                !fits.is_empty(),
                "{name}"
            );
            let on_edges = contours
                .arcs
                .iter()
                .filter(|a| a.on_faces.len() > 1)
                .count();
            let counts = (contours.arcs.len(), torus, fits.len(), on_edges);
            eprintln!(
                "P3 {name} at {threads} threads: {counts:?} (arcs, torus, fitted, on edges), \
                 {} work, fits {:?}",
                work.used,
                fits.iter()
                    .map(|f| (f.spans, f.deviation_mm))
                    .collect::<Vec<_>>()
            );
            assert_eq!(counts, expected, "{name}");
            if name == "iso" {
                // Binary64 and libm only: the same bits on every machine.
                // Recorded on arm64 macOS; CI runs x86_64 Linux.
                assert_eq!(fits[0].spans, 82);
                assert_eq!(
                    fits[0].deviation_mm.to_bits(),
                    9.276640353977601e-7_f64.to_bits()
                );
            }
            assert!(work.used < TARGET, "{name}: {} work units", work.used);
            let record = (counts, work.used, fits);
            match &reference {
                None => reference = Some(record),
                Some(first) => assert_eq!(*first, record, "{name} at {threads} threads"),
            }
        }
        set_threads(None);
    }
}

/// The arch profile of a spline extrusion: a cubic Bezier from (0, 0) over
/// (0, 4) and (4, 4) to (4, 0), closed by the floor, extruded 2 along z.
fn arch() -> catallax_app::ApplicationBody {
    let source = serde_json::json!({"format":"soliddraft-recipe/1","name":"spline-profile","steps":[
        {"id":"profile","op":"profile","plane":"XY","loops":{"outer":[
            {"type":"bezier","from":[0,0],"controls":[[0,4],[4,4]],"to":[4,0]},
            {"type":"line","from":[4,0],"to":[0,0]}
        ]}},
        {"id":"body","op":"extrude","profile":"profile","distance":2}
    ],"result":"body"});
    catallax_app::build_soliddraft(
        &source.to_string(),
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .body
}

#[test]
fn extruded_bezier_profiles_give_exact_rulings() {
    let app = arch();
    let (model, body) = (&app.built.model, app.built.result);
    // Seen along (1, 1, 0): one ruling where the arch's tangent is parallel
    // to (1, 1), at u = 1 - sqrt(2) / 2 with x = 12 u^2 - 8 u^3 and
    // y = 12 u (1 - u).
    let answer =
        contour_generators(model, body, &vector(1, 1, 0), &mut Decider::default()).unwrap();
    assert!(answer.is_exact());
    let c = &answer.value;
    assert_eq!(c.arcs.len(), 1);
    let arc = &c.arcs[0];
    assert!(arc.fit.is_none());
    assert_eq!(arc.on_faces, vec![arc.face]);
    assert!(matches!(
        curve_form(c.model.geom(), arc.curve, &mut Decider::default()).unwrap(),
        CurveForm::Line { .. }
    ));
    let u = 1.0 - 0.5_f64.sqrt();
    let (x, y) = (12.0 * u * u - 8.0 * u * u * u, 12.0 * u * (1.0 - u));
    let mut heights = Vec::new();
    for end in [arc.start.unwrap(), arc.end.unwrap()] {
        let p = c.model.geom().position(end);
        assert!((p.x.to_f64() - x).abs() < 1e-12 && (p.y.to_f64() - y).abs() < 1e-12);
        heights.push(p.z.to_f64());
    }
    heights.sort_by(f64::total_cmp);
    assert_eq!(heights, [0.0, 2.0]);
    // Seen along (0, 1, 1): the rulings at both ends of the profile, where
    // its tangent is vertical; each runs along the arch's edge with the floor.
    let c = contour_generators(model, body, &vector(0, 1, 1), &mut Decider::default())
        .unwrap()
        .value;
    assert_eq!(c.arcs.len(), 2);
    for arc in &c.arcs {
        assert!(arc.fit.is_none());
        assert_eq!(arc.on_faces.len(), 2);
    }
}

/// Points of a fitted arc at `count + 1` even parameters, in binary64.
fn samples(c: &Contours, arc: &ContourArc, count: usize) -> Vec<[f64; 3]> {
    let Curve::BSpline(curve) = c.model.geom().curve(arc.curve) else {
        panic!("a fitted arc is a B-spline");
    };
    let [from, to] = curve.domain().map(Rational::to_f64_nearest);
    let mut dec = Decider::default();
    (0..=count)
        .map(|k| {
            let at = from + (to - from) * k as f64 / count as f64;
            let at = Real::from(Rational::from_f64_exact(at).unwrap());
            let p = curve.point_at(&at, &mut dec).unwrap();
            [p.x.to_f64(), p.y.to_f64(), p.z.to_f64()]
        })
        .collect()
}

fn position(c: &Contours, point: Option<PointId>) -> [f64; 3] {
    let p = c.model.geom().position(point.unwrap());
    [p.x.to_f64(), p.y.to_f64(), p.z.to_f64()]
}

/// Contours at 1, 2 and 4 threads: the same arcs, fits and work each time.
fn fitted_contours(model: &Model, body: BodyId, d: &RVector3) -> Contours {
    let mut reference = None;
    let mut first = None;
    for threads in [1, 2, 4] {
        set_threads(NonZeroUsize::new(threads));
        let (answer, work) = with_operation_budget(DEFAULT_BUILD_WORK, || {
            contour_generators(model, body, d, &mut Decider::default())
        });
        assert!(!work.exhausted);
        let answer = answer.unwrap();
        assert!(answer.marks.iter().any(|m| m.code == CONTOUR_FIT_MARK));
        assert!(work.used < TARGET, "{} work units", work.used);
        let fits: Vec<_> = answer.value.arcs.iter().map(|a| a.fit.clone()).collect();
        for fit in &fits {
            let fit = fit.as_ref().expect("every spline arc is fitted");
            assert!(fit.deviation_mm <= fit.target_mm);
        }
        eprintln!(
            "{} arcs at {threads} threads, {} work, fits {:?}",
            fits.len(),
            work.used,
            fits.iter()
                .flatten()
                .map(|f| (f.spans, f.deviation_mm))
                .collect::<Vec<_>>()
        );
        let record = (work.used, fits);
        match &reference {
            None => {
                reference = Some(record);
                first = Some(answer.value);
            }
            Some(r) => assert_eq!(*r, record, "at {threads} threads"),
        }
    }
    set_threads(None);
    first.unwrap()
}

#[test]
fn rational_patch_roofs_give_fitted_generators_on_their_known_curves() {
    // The roof z = 1 + x y over the unit square has normal (-y, -x, 1). Seen
    // along (2, 2, 1) its generator is x + y = 1/2, a parabola from
    // (1/2, 0, 1) to (0, 1/2, 1). The joined roofs meet at x = 1/2; seen
    // along (4, 4, 3) the generator x + y = 3/4 crosses that joint at
    // (1/2, 1/4), one arc on each roof.
    for (kind, d, sum, arcs) in [
        (Roof::Bilinear, vector(2, 2, 1), 0.5, 1),
        (Roof::Rational, vector(2, 2, 1), 0.5, 1),
        (Roof::JoinedDegreeElevated, vector(4, 4, 3), 0.75, 2),
        (Roof::JoinedNonlinear, vector(4, 4, 3), 0.75, 2),
    ] {
        let fixture = solid(kind);
        let c = fitted_contours(&fixture.model, fixture.body, &d);
        assert_eq!(c.arcs.len(), arcs, "{kind:?}");
        let mut ends = Vec::new();
        for arc in &c.arcs {
            assert!(fixture.roofs.contains(&arc.face));
            assert_eq!(arc.on_faces, vec![arc.face]);
            // Within the deviation of a true point, the sample is within
            // sqrt(2) of it from the line and three of it from the roof.
            let bound = 4.0 * arc.fit.as_ref().unwrap().deviation_mm;
            for p in samples(&c, arc, 64) {
                assert!((p[0] + p[1] - sum).abs() <= bound, "{kind:?} {p:?}");
                assert!((p[2] - 1.0 - p[0] * p[1]).abs() <= bound, "{kind:?} {p:?}");
                assert!(p[0] >= -bound && p[1] >= -bound, "{kind:?} {p:?}");
            }
            ends.push(position(&c, arc.start));
            ends.push(position(&c, arc.end));
        }
        // The ends lie on the roof's boundary: y = 0, x = 0 or the joint.
        let mut xs: Vec<f64> = ends.iter().map(|p| p[0]).collect();
        xs.sort_by(f64::total_cmp);
        let expected: &[f64] = if arcs == 1 {
            &[0.0, 0.5]
        } else {
            &[0.0, 0.5, 0.5, 0.75]
        };
        for (x, e) in xs.iter().zip(expected) {
            assert!((x - e).abs() < 1e-9, "{kind:?} ends {ends:?}");
        }
    }
}

#[test]
fn a_parabolic_roof_gives_its_straight_ruling() {
    // z = 1 + x (1 - x) with a nonuniform knot in x has normal
    // (2 x - 1, 0, 1). Seen along (2, 0, 1) its generator is the line
    // x = 1/4, z = 19/16 across the roof.
    let fixture = solid(Roof::NonuniformBump);
    let c = fitted_contours(&fixture.model, fixture.body, &vector(2, 0, 1));
    assert_eq!(c.arcs.len(), 1);
    let arc = &c.arcs[0];
    let bound = 2.0 * arc.fit.as_ref().unwrap().deviation_mm;
    for p in samples(&c, arc, 64) {
        assert!((p[0] - 0.25).abs() <= bound && (p[2] - 19.0 / 16.0).abs() <= bound);
    }
    let mut ys = [position(&c, arc.start)[1], position(&c, arc.end)[1]];
    ys.sort_by(f64::total_cmp);
    assert!(ys[0].abs() < 1e-9 && (ys[1] - 1.0).abs() < 1e-9, "{ys:?}");
}

#[test]
fn a_lofted_barrel_gives_fitted_generators_in_the_plane_through_its_axis() {
    // Radii 1, 2, 1 at heights 0, 1, 2 about z: a solid of revolution, whose
    // generators seen along (1, 2, 0) lie in the plane x + 2 y = 0.
    let recipe = Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/loft-natural-circles.json"
    ))
    .unwrap();
    let (built, _) = with_operation_budget(DEFAULT_BUILD_WORK, || {
        evaluate(&recipe, &mut Decider::default())
    });
    let built = built.unwrap();
    let c = fitted_contours(&built.model, built.result, &vector(1, 2, 0));
    assert!(c.arcs.len() >= 2);
    let mut heights = Vec::new();
    for arc in &c.arcs {
        let bound = 2.0 * arc.fit.as_ref().unwrap().deviation_mm;
        for p in samples(&c, arc, 32) {
            assert!((p[0] + 2.0 * p[1]).abs() / 5.0_f64.sqrt() <= bound, "{p:?}");
        }
        heights.push(position(&c, arc.start)[2]);
        heights.push(position(&c, arc.end)[2]);
    }
    heights.sort_by(f64::total_cmp);
    // Two generators, one on each side, each from the bottom to the top.
    assert!(
        heights[0].abs() < 1e-9 && heights[1].abs() < 1e-9,
        "{heights:?}"
    );
    let n = heights.len();
    assert!((heights[n - 1] - 2.0).abs() < 1e-9 && (heights[n - 2] - 2.0).abs() < 1e-9);
    // The arcs of one generator meet at the chart seams it crosses: they
    // share their end points there, so each generator is one chain from the
    // bottom to the top with no loose end inside the face.
    let mut uses: std::collections::BTreeMap<PointId, Vec<usize>> = Default::default();
    for (i, arc) in c.arcs.iter().enumerate() {
        for end in [arc.start.unwrap(), arc.end.unwrap()] {
            uses.entry(end).or_default().push(i);
        }
    }
    let mut loose = Vec::new();
    for (point, arcs) in &uses {
        match arcs.len() {
            1 => loose.push(position(&c, Some(*point))[2]),
            2 => assert_ne!(arcs[0], arcs[1], "an arc closes on itself"),
            k => panic!("{k} arcs meet at one point"),
        }
    }
    loose.sort_by(f64::total_cmp);
    assert_eq!(loose.len(), 4, "loose ends at heights {loose:?}");
    assert!(loose[0].abs() < 1e-9 && loose[1].abs() < 1e-9, "{loose:?}");
    assert!(
        (loose[2] - 2.0).abs() < 1e-9 && (loose[3] - 2.0).abs() < 1e-9,
        "{loose:?}"
    );
    // Two chains: follow shared ends from each arc.
    let mut chain = vec![usize::MAX; c.arcs.len()];
    let mut chains = 0;
    for first in 0..c.arcs.len() {
        if chain[first] != usize::MAX {
            continue;
        }
        let mut stack = vec![first];
        while let Some(i) = stack.pop() {
            if chain[i] != usize::MAX {
                continue;
            }
            chain[i] = chains;
            for end in [c.arcs[i].start.unwrap(), c.arcs[i].end.unwrap()] {
                stack.extend(uses[&end].iter().copied());
            }
        }
        chains += 1;
    }
    assert_eq!(chains, 2, "{} arcs", c.arcs.len());
}

/// The ball of radius 5 about the origin less the cylinder of radius 2
/// along z through (1, 0): the hole's wall meets the sphere in two closed
/// intersection branches, which trim the generators of both faces.
fn holed_ball() -> (Model, BodyId) {
    let mut model = Model::new();
    let ball = sphere(&mut model, 5.into(), &Pose::identity()).unwrap();
    let hole = cylinder(
        &mut model,
        2.into(),
        20.into(),
        &Pose::from_rational(Rotation::identity(), vector(1, 0, -10)),
    )
    .unwrap();
    let body = boolean(
        &mut model,
        ball,
        hole,
        BoolOp::Subtract,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    (model, body)
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn unit3(v: [f64; 3]) -> [f64; 3] {
    let n = dot3(v, v).sqrt();
    v.map(|x| x / n)
}

/// Distance from the hole's axis, squared.
fn from_axis(p: [f64; 3]) -> f64 {
    (p[0] - 1.0).powi(2) + p[1] * p[1]
}

/// A binary64 image of exact coordinates, good to the last few bits.
fn close(v: [&Real; 3]) -> [f64; 3] {
    let floor = Rational::from_f64_exact(1e-12).unwrap();
    let mut dec = Decider::default();
    v.map(|x| dec.enclose(x, &floor).unwrap().mid().to_f64_nearest())
}

fn point_image(c: &Contours, point: Option<PointId>) -> [f64; 3] {
    let p = c.model.geom().position(point.unwrap());
    close([&p.x, &p.y, &p.z])
}

/// The binary64 middle and length of a line or circle arc, read from its
/// drawing, whose direction along the carrier is decided exactly.
fn middle_and_length(c: &Contours, arc: &ContourArc) -> ([f64; 3], f64) {
    let drawn = drawing_curve(
        &c.model,
        arc.curve,
        arc.start,
        arc.end,
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert!(drawn.fit.is_none());
    match drawn.curve {
        DrawingCurve::Line { a, b } => {
            let ab: [f64; 3] = std::array::from_fn(|k| b[k] - a[k]);
            (
                std::array::from_fn(|k| (a[k] + b[k]) / 2.0),
                dot3(ab, ab).sqrt(),
            )
        }
        DrawingCurve::Arc {
            center,
            normal,
            x_axis,
            radius,
            start_deg,
            end_deg,
        } => {
            let y = cross3(normal, x_axis);
            let (sin, cos) = ((start_deg + end_deg) / 2.0).to_radians().sin_cos();
            (
                std::array::from_fn(|k| center[k] + radius * (cos * x_axis[k] + sin * y[k])),
                radius * (end_deg - start_deg).to_radians(),
            )
        }
        DrawingCurve::Circle {
            center,
            normal,
            radius,
        } => {
            let axis = if normal[0].abs() < 0.9 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 1.0, 0.0]
            };
            let x = unit3(cross3(normal, axis));
            (
                std::array::from_fn(|k| center[k] + radius * x[k]),
                std::f64::consts::TAU * radius,
            )
        }
        other => panic!("unexpected drawing {other:?}"),
    }
}

#[test]
fn generators_trimmed_by_intersection_branches_stay_exact_and_inside_the_target() {
    let (model, body) = holed_ball();
    // From (2, -3, 1) the sphere's generator, a great circle, comes within
    // 0.07 of the hole's wall where it does not cross it, and the cuts
    // where the branches cross it are algebraic points of high degree.
    for d in [
        vector(2, -3, 1),
        vector(1, 0, 1),
        vector(1, -2, 2),
        vector(-1, 1, -2),
        vector(0, 0, 1),
        vector(0, -1, 0),
        vector(1, 3, 0),
    ] {
        let d64 = unit3([
            d.x.to_f64_nearest(),
            d.y.to_f64_nearest(),
            d.z.to_f64_nearest(),
        ]);
        let mut reference = None;
        for threads in [1, 2, 4] {
            set_threads(NonZeroUsize::new(threads));
            let (answer, work) = with_operation_budget(DEFAULT_BUILD_WORK, || {
                contour_generators(&model, body, &d, &mut Decider::default())
            });
            assert!(!work.exhausted, "{d:?}");
            let answer = answer.unwrap();
            assert!(answer.is_exact(), "{d:?}");
            eprintln!(
                "holed ball {d:?} at {threads} threads: {} arcs, {} work",
                answer.value.arcs.len(),
                work.used
            );
            assert!(work.used < TARGET, "{d:?}: {} work units", work.used);
            let c = &answer.value;
            let (mut sphere_length, mut rulings) = (0.0, Vec::new());
            for arc in &c.arcs {
                assert!(arc.fit.is_none(), "{d:?}");
                for end in [arc.start, arc.end].into_iter().flatten() {
                    let p = point_image(c, Some(end));
                    assert!((dot3(p, p) - 25.0).abs() < 1e-9, "{d:?}: end {p:?}");
                    assert!((from_axis(p) - 4.0).abs() < 1e-9, "{d:?}: end {p:?}");
                }
                let (m, length) = middle_and_length(c, arc);
                if matches!(
                    model.geom().surface(model.face(arc.face).surface),
                    Surface::Sphere(_)
                ) {
                    assert!((dot3(m, m) - 25.0).abs() < 1e-9, "{d:?}: {m:?}");
                    assert!(from_axis(m) > 4.0, "{d:?}: {m:?} lies in the hole");
                    assert!(dot3(unit3(m), d64).abs() < 1e-12, "{d:?}: {m:?}");
                    sphere_length += length;
                } else {
                    assert!((from_axis(m) - 4.0).abs() < 1e-9, "{d:?}: {m:?}");
                    assert!(dot3(m, m) < 25.0, "{d:?}: {m:?} lies outside the ball");
                    let normal = unit3([m[0] - 1.0, m[1], 0.0]);
                    assert!(dot3(normal, d64).abs() < 1e-12, "{d:?}: {m:?}");
                    // A whole ruling inside the ball, from branch to branch.
                    let reach = 2.0 * (25.0 - m[0] * m[0] - m[1] * m[1]).sqrt();
                    assert!((length - reach).abs() < 1e-9, "{d:?}: {length} {reach}");
                    rulings.push(m);
                }
            }
            // The great circle across the view, sampled: the part outside
            // the hole is what the sphere's arcs cover.
            let e1 = unit3(cross3(
                d64,
                if d64[2].abs() < 0.9 {
                    [0.0, 0.0, 1.0]
                } else {
                    [1.0, 0.0, 0.0]
                },
            ));
            let e2 = cross3(d64, e1);
            let samples = 1 << 16;
            let outside = (0..samples)
                .filter(|&k| {
                    let (sin, cos) = (std::f64::consts::TAU * k as f64 / samples as f64).sin_cos();
                    from_axis(std::array::from_fn(|i| 5.0 * (cos * e1[i] + sin * e2[i]))) > 4.0
                })
                .count();
            let expected = std::f64::consts::TAU * 5.0 * outside as f64 / samples as f64;
            let step = std::f64::consts::TAU * 5.0 / samples as f64;
            assert!(
                (sphere_length - expected).abs() <= 2.0 * step,
                "{d:?}: {sphere_length} against {expected}"
            );
            // The wall's rulings where its normal is across the view: two,
            // unless the view runs along the axis.
            let across = d64[0].abs() > 0.0 || d64[1].abs() > 0.0;
            assert_eq!(rulings.len(), if across { 2 } else { 0 }, "{d:?}");
            let record = (
                c.arcs
                    .iter()
                    .map(|a| {
                        (
                            a.face,
                            a.start.map(|p| point_image(c, Some(p))),
                            a.end.map(|p| point_image(c, Some(p))),
                        )
                    })
                    .collect::<Vec<_>>(),
                work.used,
            );
            match &reference {
                None => reference = Some(record),
                Some(first) => assert_eq!(*first, record, "{d:?} at {threads} threads"),
            }
        }
        set_threads(None);
    }
}
