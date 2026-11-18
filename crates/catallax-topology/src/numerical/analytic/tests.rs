use super::*;
use catallax_shapes::{Circle, Plane, RPoint3, RVector3};

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn split_annulus(center: &str, split_outer: bool) -> (NumericalFace, Vec<NumericalEdge>) {
    let circle = |x: Rational, radius: Rational| Circle {
        center: RPoint3::new(x, q("0"), q("0")).to_real(),
        normal: RVector3::new(q("0"), q("0"), q("1")).to_real(),
        radius_sq: (&radius * &radius).into(),
        ref_dir: RVector3::new(q("1"), q("0"), q("0")).to_real(),
    };
    let turn = period();
    let mut edges = vec![
        NumericalEdge {
            key: "outer".into(),
            curve: Curve::Circle(circle(q("0"), q("3"))),
            parameters: [q("0"), turn.clone()],
            vertices: [0, 0],
        },
        NumericalEdge {
            key: "inner".into(),
            curve: Curve::Circle(circle(q(center), q("1"))),
            parameters: [q("0"), turn.clone()],
            vertices: [1, 1],
        },
    ];
    let split = usize::from(!split_outer);
    let curve = edges[split].curve.clone();
    for i in 0..2 {
        edges.push(NumericalEdge {
            key: format!("arc-{i}"),
            curve: curve.clone(),
            parameters: [
                &turn * &Rational::from_ratio(2 * i + 1, 4).unwrap(),
                &turn * &Rational::from_ratio(2 * i + 3, 4).unwrap(),
            ],
            vertices: if i == 0 { [2, 3] } else { [3, 2] },
        });
    }
    let loops = [q("0"), q(center)]
        .into_iter()
        .enumerate()
        .map(|(i, x)| {
            let radius = q(if i == 0 { "3" } else { "1" });
            let sense = if i == 0 { Sense::Same } else { Sense::Reversed };
            if i != split {
                let at = [&x + &radius, q("0")];
                vec![NumericalTrim {
                    edge: i,
                    sense,
                    uv: [at.clone(), at],
                }]
            } else {
                let a = [x.clone(), radius.clone()];
                let b = [x, -radius];
                let mut trims = vec![
                    NumericalTrim {
                        edge: 2,
                        sense,
                        uv: [a.clone(), b.clone()],
                    },
                    NumericalTrim {
                        edge: 3,
                        sense,
                        uv: [b, a],
                    },
                ];
                if sense == Sense::Reversed {
                    trims.reverse();
                    for trim in &mut trims {
                        trim.uv.swap(0, 1);
                    }
                }
                trims
            }
        })
        .collect();
    (
        NumericalFace {
            key: "annulus".into(),
            surface: Surface::Plane(Plane {
                origin: RPoint3::new(q("0"), q("0"), q("0")).to_real(),
                normal: RVector3::new(q("0"), q("0"), q("1")).to_real(),
                ref_dir: RVector3::new(q("1"), q("0"), q("0")).to_real(),
            }),
            sense: Sense::Same,
            loops,
        },
        edges,
    )
}

#[test]
fn complete_rims_clear_finite_arc_neighbours_in_both_nesting_directions() {
    for split_outer in [false, true] {
        let (face, edges) = split_annulus("1", split_outer);
        circle_clearance(&face, &edges, &q("1/1000"), &mut Decider::default()).unwrap();
    }
}

#[test]
fn arc_bulges_and_sub_tolerance_rim_gaps_do_not_pass_endpoint_clearance() {
    for (center, split_outer) in [("5/2", false), ("5/2", true), ("1999/1000", false)] {
        let (face, edges) = split_annulus(center, split_outer);
        let error =
            circle_clearance(&face, &edges, &q("1/1000"), &mut Decider::default()).unwrap_err();
        assert_eq!(error.code, "numerical.invalid-trim");
    }
}
