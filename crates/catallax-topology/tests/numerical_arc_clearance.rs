//! Circular holes beside finite arcs retain continuous clearance and winding.
#[path = "common/periodic.rs"]
mod fixture;

use catallax_numbers::{Ball, Decider, NumericalCheck};
use catallax_shapes::{Circle, Curve, Cylinder, RPoint3, RVector3, Surface};
use catallax_topology::{
    Sense,
    numerical::{self as n, NumericalDefinition, NumericalOptions},
};

fn arched_bore(inward: bool, center: [&str; 2], radius: &str) -> NumericalDefinition {
    let mut definition = fixture::bores(&[center[0]], radius);
    let shift = fixture::q(center[1]);
    let real_shift = catallax_numbers::Real::from(shift.clone());
    for vertex in &mut definition.vertices[8..] {
        vertex.position[1] = &vertex.position[1] + &shift;
    }
    for edge in &mut definition.edges[12..] {
        let origin = match &mut edge.curve {
            Curve::Circle(circle) => &mut circle.center,
            Curve::Line(line) => &mut line.origin,
            _ => unreachable!(),
        };
        origin.y = &origin.y + &real_shift;
    }
    for cap in &mut definition.faces[..2] {
        for uv in &mut cap.loops[1][0].uv {
            uv[1] = &uv[1] + &shift;
        }
    }
    let Surface::Cylinder(bore) = &mut definition.faces.last_mut().unwrap().surface else {
        unreachable!()
    };
    bore.origin.y = &bore.origin.y + &real_shift;

    let pi = Ball::pi(128).mid();
    let left = if inward { pi } else { -pi };
    for edge_index in [
        definition.faces[2].loops[0][0].edge,
        definition.faces[2].loops[0][2].edge,
    ] {
        let edge = &mut definition.edges[edge_index];
        let start = &definition.vertices[edge.vertices[0]].position;
        edge.curve = Curve::Circle(Circle {
            center: RPoint3::new(fixture::q("0"), fixture::q("-2"), start[2].clone()).to_real(),
            normal: RVector3::new(fixture::q("0"), fixture::q("0"), fixture::q("1")).to_real(),
            radius_sq: 4.into(),
            ref_dir: RVector3::new(fixture::q("1"), fixture::q("0"), fixture::q("0")).to_real(),
        });
        edge.parameters = if start[0] > fixture::q("0") {
            [fixture::q("0"), left.clone()]
        } else {
            [left.clone(), fixture::q("0")]
        };
    }
    let wall = &mut definition.faces[2];
    wall.surface = Surface::Cylinder(Cylinder {
        origin: RPoint3::new(fixture::q("0"), fixture::q("-2"), fixture::q("0")).to_real(),
        axis: RVector3::new(fixture::q("0"), fixture::q("0"), fixture::q("1")).to_real(),
        radius_sq: 4.into(),
        ref_dir: RVector3::new(fixture::q("1"), fixture::q("0"), fixture::q("0")).to_real(),
    });
    wall.sense = if inward { Sense::Reversed } else { Sense::Same };
    let corners = [
        [left.clone(), fixture::q("0")],
        [fixture::q("0"), fixture::q("0")],
        [fixture::q("0"), fixture::q("1")],
        [left, fixture::q("1")],
    ];
    for (i, trim) in wall.loops[0].iter_mut().enumerate() {
        trim.uv = [corners[i].clone(), corners[(i + 1) % 4].clone()];
    }
    definition
}

#[test]
fn a_finite_arc_can_clear_a_hole_when_its_complete_carrier_does_not() {
    let body = n::build(
        arched_bore(false, ["0", "0"], "1"),
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!((body.counts().shells, body.counts().lumps), (1, 1));
    for cap in ["f0", "f1"] {
        assert!(body.evidence().checks().iter().any(|check| matches!(check,
            NumericalCheck::Passed { name, scope }
                if name == "trim.continuous-circle-hole-clearance" && scope == cap)));
    }
    let arc = &body.definition().edges[3].curve;
    let at = Ball::from_rational(&(-Ball::pi(128).mid()).mul_pow2(-1), 128);
    let point = arc.evaluate(&at, &Decider::default()).unwrap().unwrap();
    assert!(point[0].contains(&fixture::q("0")));
    assert!(point[1].contains(&fixture::q("-4")));
}

#[test]
fn circle_clearance_cannot_transfer_winding_through_an_unsampled_arc_bulge() {
    // Both centers lie between polygon stations on the inward upper arc.
    // The first hole crosses the true arc. The second lies wholly on its
    // wrong side, in the lens between the arc and the sampled chord.
    for center in [["6/5", "-2/5"], ["1.1976", "-0.4032"]] {
        let error = n::build(
            arched_bore(true, center, "1/1000"),
            NumericalOptions {
                tolerance_mm: fixture::q("1e-6"),
            },
            &mut Decider::default(),
        )
        .unwrap_err();
        assert_eq!(error.code, "numerical.invalid-trim");
    }
}

#[test]
fn finite_arc_clearance_reserves_both_boundaries_requested_tolerance() {
    for (tolerance, passes) in [("1/5", false), ("1/100", true)] {
        let result = n::build(
            arched_bore(true, ["0", "4/5"], "1/2"),
            NumericalOptions {
                tolerance_mm: fixture::q(tolerance),
            },
            &mut Decider::default(),
        );
        if passes {
            assert_eq!(
                result.unwrap().evidence().tolerance_mm(),
                &fixture::q(tolerance)
            );
        } else {
            assert_eq!(result.unwrap_err().code, "numerical.invalid-trim");
        }
    }
}
