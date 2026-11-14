//! A rational diagonal torus boundary, independently of fillet construction.

use catallax_boolean::face_contains;
use catallax_numbers::{Decider, Rational, Real};
use catallax_shapes::{BSplineCurve, Curve, Point3, PointDef, Surface, Torus, Vector3};
use catallax_topology::{Coedge, Edge, Face, Loop, Model, Sense, Vertex};

fn r(value: &str) -> Real {
    Rational::parse(value).unwrap().into()
}

#[test]
fn noncircular_torus_loop_keeps_cells_across_the_seam_and_chart_pole() {
    let mut model = Model::new();
    let mut dec = Decider::default();
    let surface = model.geom_mut().add_surface(Surface::Torus(Torus {
        center: Point3::new(r("0"), r("0"), r("0")),
        axis: Vector3::new(r("0"), r("0"), r("1")),
        ref_dir: Vector3::new(r("1"), r("0"), r("0")),
        major: r("3"),
        minor: r("1"),
    }));
    let vertices = [["4", "0", "0"], ["0", "4", "0"], ["3", "0", "1"]].map(|xyz| {
        let [x, y, z] = xyz.map(r);
        let point = model
            .geom_mut()
            .add_point(PointDef::Explicit(Point3::new(x, y, z)));
        model.add_vertex(Vertex { point })
    });
    // Rational half-angle parameters (u,v) map to
    // ((3+cos(v))*cos(u),(3+cos(v))*sin(u),sin(v)).
    // The triangle (0,0),(1,0),(0,1) has one diagonal quartic edge.
    // These are its exact homogeneous Bernstein coefficients, not samples.
    let nets = [
        vec![
            ["4", "0", "0", "1"],
            ["4", "4", "0", "1"],
            ["0", "8", "0", "2"],
        ],
        vec![
            ["0", "8", "0", "2"],
            ["2", "6", "1", "3/2"],
            ["10/3", "14/3", "4/3", "3/2"],
            ["5", "3", "3/2", "3/2"],
            ["6", "0", "2", "2"],
        ],
        vec![
            ["6", "0", "2", "2"],
            ["4", "0", "1", "1"],
            ["4", "0", "0", "1"],
        ],
    ];
    let mut coedges = Vec::new();
    for (i, net) in nets.into_iter().enumerate() {
        let degree = net.len() - 1;
        let mut knots = vec![Rational::zero(); degree + 1];
        knots.extend(vec![Rational::one(); degree + 1]);
        let curve = BSplineCurve::from_homogeneous(
            degree,
            knots,
            net.into_iter().map(|row| row.map(r)).collect(),
            &mut dec,
        )
        .unwrap();
        let curve = model.geom_mut().add_curve(Curve::BSpline(curve));
        let edge = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(vertices[i]),
            end: Some(vertices[(i + 1) % 3]),
        });
        coedges.push(model.add_coedge(Coedge {
            edge,
            sense: Sense::Same,
        }));
    }
    let lp = model.add_loop(Loop::Coedges(coedges.clone()));
    let face = model.add_face(Face {
        surface,
        sense: Sense::Same,
        loops: vec![lp],
    });
    let mut faces = vec![(face, false)];
    for (sense, complement) in [(Sense::Same, true), (Sense::Reversed, false)] {
        let reversed = coedges
            .iter()
            .rev()
            .map(|coedge| {
                let edge = model.coedge(*coedge).edge;
                model.add_coedge(Coedge {
                    edge,
                    sense: Sense::Reversed,
                })
            })
            .collect();
        let lp = model.add_loop(Loop::Coedges(reversed));
        faces.push((
            model.add_face(Face {
                surface,
                sense,
                loops: vec![lp],
            }),
            complement,
        ));
    }
    for (face, complement) in faces {
        for (xyz, expected) in [
            (["574/625", "1968/625", "24/25"], false), // u=v=3/4
            (["990/289", "528/289", "8/17"], true),    // u=v=1/4
            (["3", "0", "-1"], false),                 // omitted south parallel
            (["990/289", "-528/289", "8/17"], false),  // across the azimuth seam
            (["990/289", "528/289", "8/17"], true),
        ] {
            let [x, y, z] = xyz.map(r);
            assert_eq!(
                face_contains(&mut model, face, &Point3::new(x, y, z), &mut dec).unwrap(),
                expected ^ complement,
                "{xyz:?}"
            );
        }
    }
}
