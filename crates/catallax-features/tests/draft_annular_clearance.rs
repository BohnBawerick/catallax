//! A steep thin annulus whose closest wall points have different axial heights.
use catallax_numbers::{Decider, Fact, FactKind, FactOrigin, Rational, Real};
use catallax_shapes::{
    BSplineCurve, Curve, Cylinder, GeomKey, ObliqueDraft, Plane, Point3, SplineBoundary,
    SplineSupport, Surface, Vector3,
};
use catallax_topology::{
    Body, BodyKind, Coedge, DraftedCylinder, Edge, Face, Loop, Lump, Model, Sense, Shell,
};

fn q(value: &str) -> Real {
    Real::from_rational(Rational::parse(value).unwrap())
}

fn annulus(rotated: bool) -> (Model, catallax_topology::BodyId) {
    let mut dec = Decider::default();
    let origin = Point3::new(q("1/3"), q("2/7"), q("3/11"));
    let axis = if rotated {
        Vector3::new(q("0"), q("3/5"), q("4/5"))
    } else {
        Vector3::new(q("0"), q("0"), q("1"))
    };
    let across = Vector3::new(q("1"), q("0"), q("0"));
    let neutral = Plane {
        origin: origin.clone(),
        normal: across.plus(&axis.scaled(&q("10000"))),
        ref_dir: axis.cross(&across),
    };
    let supports = ["81/20", "4"].map(|radius| {
        ObliqueDraft::new(
            &Cylinder {
                origin: origin.clone(),
                axis: axis.clone(),
                radius_sq: q(radius).square(),
                ref_dir: across.clone(),
            },
            axis.clone(),
            &neutral,
            q("100"),
            [q("0"), q("1/1000")],
            &mut dec,
        )
        .unwrap()
    });
    let mut model = Model::new();
    let rings = supports.each_ref().map(|support| {
        [SplineBoundary::VMin, SplineBoundary::VMax].map(|boundary| {
            let curve = model
                .geom_mut()
                .add_curve(Curve::BSpline(BSplineCurve::boundary(
                    support.carrier(),
                    boundary,
                )));
            model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: None,
                end: None,
                tolerance: None,
            })
        })
    });
    let mut add_face = |surface, sense, paths: Vec<Vec<_>>| {
        let surface = model.geom_mut().add_surface(surface);
        let loops = paths
            .into_iter()
            .map(|path| {
                let coedges = path
                    .into_iter()
                    .map(|(edge, sense)| {
                        let curve = model.edge(edge).curve;
                        model.geom_mut().record_fact(Fact {
                            kind: FactKind::Incident,
                            a: GeomKey::Curve(curve),
                            b: GeomKey::Surface(surface),
                            at: None,
                            origin: FactOrigin::Construction("test.shared-oblique-row"),
                        });
                        model.add_coedge(Coedge { edge, sense })
                    })
                    .collect();
                model.add_loop(Loop::Coedges(coedges))
            })
            .collect();
        model.add_face(Face {
            surface,
            sense,
            loops,
            tolerance: None,
        })
    };
    let mut faces = Vec::new();
    for i in 0..2 {
        let cap = Plane {
            origin: origin.plus(&axis.scaled(&q(if i == 0 { "0" } else { "1/1000" }))),
            normal: axis.scaled(&q(if i == 0 { "-1" } else { "1" })),
            ref_dir: across.clone(),
        };
        let outer = if i == 0 { Sense::Reversed } else { Sense::Same };
        let inner = if i == 0 { Sense::Same } else { Sense::Reversed };
        faces.push(add_face(
            Surface::Plane(cap),
            Sense::Same,
            vec![vec![(rings[0][i], outer)], vec![(rings[1][i], inner)]],
        ));
    }
    for (i, support) in supports.into_iter().enumerate() {
        let lower = if i == 0 { Sense::Same } else { Sense::Reversed };
        let upper = if i == 0 { Sense::Reversed } else { Sense::Same };
        faces.push(add_face(
            Surface::BSpline(SplineSupport::ObliqueDraft(support)),
            if i == 0 { Sense::Same } else { Sense::Reversed },
            vec![vec![(rings[i][0], lower)], vec![(rings[i][1], upper)]],
        ));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        wires: Vec::new(),
        marks: Default::default(),
    });
    (model, body)
}

#[test]
fn steep_narrow_annulus_bounds_cross_height_separation_before_xyz_rounding() {
    for rotated in [false, true] {
        let (model, body) = annulus(rotated);
        let mut dec = Decider::default();
        let checked = DraftedCylinder::check(&model, body, &mut dec).unwrap();
        let inner = checked.hole().unwrap().0;
        let outer_point = checked.support().point([&q("1"), &q("0")], &q("1/2000"));
        let inner_point = inner.point([&q("1"), &q("0")], &q("1/200000"));
        let witness = dec
            .enclose_at(&outer_point.minus(&inner_point).norm_squared(), 128)
            .unwrap();
        let clearance = checked.wall_clearance(&mut dec).unwrap().unwrap();
        let bound = dec.enclose_at(&clearance, 128).unwrap();
        assert!(
            bound.hi() * bound.hi() <= witness.lo(),
            "clearance {bound:?} exceeds an actual cross-height wall distance {witness:?}"
        );

        let mesh = catallax_mesh::mesh_body(&model, body, &Default::default(), &mut dec).unwrap();
        assert_eq!(mesh.value.faces.len(), 4);
        assert!(
            mesh.value
                .faces
                .iter()
                .all(|face| face.max_chord_error.mul_pow2(1) < bound.lo())
        );
        let mut bytes = Vec::new();
        catallax_step::write_step(&model, &[body], &Default::default(), &mut dec, &mut bytes)
            .unwrap();
        assert!(bytes.starts_with(b"ISO-10303-21;"));
    }
}
