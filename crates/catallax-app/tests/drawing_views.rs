//! Contour generators on release parts: the P1 bracket and an elliptic extrusion.

use std::num::NonZeroUsize;

use catallax_app::DEFAULT_BUILD_WORK;
use catallax_app::views::{Contours, contour_generators};
use catallax_features::evaluate;
use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, Real, set_threads, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{RVector3, Surface};
use catallax_topology::{BodyId, Model};

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
