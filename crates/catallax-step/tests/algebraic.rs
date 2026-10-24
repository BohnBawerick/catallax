//! Oblique swept-circle definitions, independent mass integrals and STEP readback.
mod common;
use catallax_intersect::{Component, CurveForm, CurveSurface, SurfacePair};
use catallax_numbers::{
    Ball, Decider, Fact, FactKind, FactOrigin, MarkLog, Outcome, Rational, Real,
};
use catallax_shapes::{
    Circle, Curve, EllipticCylinder, GeomKey, Line, Plane, Point3, RPoint3, RVector3, Surface,
};
use catallax_topology::{Body, BodyKind, Coedge, Edge, Face, Loop, Lump, Model, Sense, Shell};
fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn p([x, y, z]: [i64; 3]) -> Point3 {
    RPoint3::new(r(x), r(y), r(z)).into()
}
fn v([x, y, z]: [i64; 3]) -> catallax_shapes::Vector3 {
    RVector3::new(r(x), r(y), r(z)).into()
}
fn loop_for(
    model: &mut Model,
    edge: catallax_topology::EdgeId,
    sense: Sense,
) -> catallax_topology::LoopId {
    let coedge = model.add_coedge(Coedge { edge, sense });
    model.add_loop(Loop::Coedges(vec![coedge]))
}
fn ball_range(lo: Rational, hi: Rational) -> Ball {
    Ball::from_rational(&lo, 128).hull(&Ball::from_rational(&hi, 128))
}
fn elliptic_e_half() -> Ball {
    // Complete E(m) = pi/2 * (1-sum c_n*m^n/(2n-1)), m=1/2.
    // Positive successive terms have ratio below 1/2, so the remaining tail
    // is bounded by twice the first omitted term. No boundary quadrature is reused.
    let mut term = Rational::from_ratio(1, 8).unwrap();
    let mut sum = Rational::zero();
    for n in 1..=80_i64 {
        sum = &sum + &term;
        term = &term
            * &Rational::from_ratio((2 * n - 1) * (2 * n + 1), 8 * (n + 1) * (n + 1)).unwrap();
    }
    let lo = &(&Rational::one() - &sum) - &term.mul_pow2(1);
    let hi = &Rational::one() - &sum;
    Ball::pi(128) * ball_range(lo.mul_pow2(-1), hi.mul_pow2(-1))
}
#[test]
fn oblique_circle_extrusion_has_exact_measurements_and_step_geometry() {
    let mut dec = Decider::default();
    let mut model = Model::new();
    // Circle radius 2, normal n=(0,0,1), swept by (2,0,2).
    // e=(1,0,1)/sqrt(2), a^2=r^2=4, b^2=r^2*(n.e)^2=2.
    let analytic = Surface::EllipticCylinder(
        EllipticCylinder::new(
            p([0, 0, 0]),
            v([1, 0, 1]),
            v([0, 1, 0]),
            Real::from_i64(4),
            Real::from_i64(2),
            &mut dec,
        )
        .unwrap(),
    );
    assert_eq!(
        analytic.domain(),
        [
            catallax_shapes::ParamKind::Periodic,
            catallax_shapes::ParamKind::Line
        ]
    );
    let side = model.geom_mut().add_surface(analytic.clone());
    let mut caps = Vec::new();
    let mut edges = Vec::new();
    for (center, normal) in [([0, 0, 0], [0, 0, -1]), ([2, 0, 2], [0, 0, 1])] {
        let cap = model.geom_mut().add_surface(Surface::Plane(
            Plane::new(p(center), v(normal), v([1, 0, 0]), &mut dec).unwrap(),
        ));
        let curve = model.geom_mut().add_curve(Curve::Circle(
            Circle::new(
                p(center),
                v([0, 0, 1]),
                Real::from_i64(4),
                v([1, 0, 0]),
                &mut dec,
            )
            .unwrap(),
        ));
        for surface in [cap, side] {
            model.geom_mut().record_fact(Fact {
                kind: FactKind::Incident,
                a: GeomKey::Curve(curve),
                b: GeomKey::Surface(surface),
                at: None,
                origin: FactOrigin::Construction("oblique-circle-extrusion"),
            });
        }
        caps.push(cap);
        edges.push(model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: None,
            end: None,
        }));
    }
    let bottom_loop = loop_for(&mut model, edges[0], Sense::Reversed);
    let top_loop = loop_for(&mut model, edges[1], Sense::Same);
    let side_bottom = loop_for(&mut model, edges[0], Sense::Same);
    let side_top = loop_for(&mut model, edges[1], Sense::Reversed);
    let faces = vec![
        model.add_face(Face {
            surface: caps[0],
            sense: Sense::Same,
            loops: vec![bottom_loop],
        }),
        model.add_face(Face {
            surface: caps[1],
            sense: Sense::Same,
            loops: vec![top_loop],
        }),
        model.add_face(Face {
            surface: side,
            sense: Sense::Same,
            loops: vec![side_bottom, side_top],
        }),
    ];
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    let zero = Ball::from_rational(&r(0), 128);
    let position = analytic.evaluate(&zero, &zero, &dec).unwrap();
    for (value, expected) in position.iter().zip([0, 2, 0]) {
        assert!(value.contains(&r(expected)));
    }
    let derivatives = analytic.derivatives(&zero, &zero, &dec).unwrap();
    for (value, expected) in derivatives[0].iter().zip([-1, 0, 1]) {
        assert!(value.contains(&r(expected)));
    }
    let gradient = analytic
        .gradient(&p([0, 2, 0]))
        .expect("analytic fixture gradient");
    for (value, expected) in [&gradient.x, &gradient.y, &gradient.z]
        .into_iter()
        .zip([0, 1, 0])
    {
        assert!(
            dec.compare(value, &Real::from_i64(expected))
                .is_coincident()
        );
    }
    let section =
        catallax_intersect::intersect_surfaces(model.geom_mut(), caps[0], side, &mut dec).unwrap();
    let SurfacePair::Meets(components) = section.value else {
        panic!("missing cap section")
    };
    assert_eq!(components.len(), 1);
    let Component::Curve { curve, .. } = components[0] else {
        panic!("cap section is not a curve")
    };
    let CurveForm::Circle { radius_sq, .. } =
        catallax_intersect::curve_form(model.geom(), curve, &mut dec).unwrap()
    else {
        panic!("the original circular cap became an ellipse")
    };
    assert!(matches!(
        dec.sign(&(&radius_sq - &Real::from_i64(4))),
        Outcome::ProvenCoincident(_)
    ));
    let line = model.geom_mut().add_curve(Curve::Line(
        Line::new(p([-3, 0, 0]), v([1, 0, 0]), &mut dec).unwrap(),
    ));
    let CurveSurface::Meets(hits) =
        catallax_intersect::intersect_curve_surface(model.geom_mut(), line, side, &mut dec)
            .unwrap()
            .value
    else {
        panic!("missing elliptic-cylinder line hits")
    };
    assert_eq!(hits.len(), 2);
    for (hit, expected) in hits.iter().zip([-2, 2]) {
        let position = model.geom().position(hit.point);
        assert!(matches!(
            dec.sign(&(&position.x - &Real::from_i64(expected))),
            Outcome::ProvenCoincident(_)
        ));
        assert!(matches!(
            dec.sign(&position.y),
            Outcome::ProvenCoincident(_)
        ));
        assert!(matches!(
            dec.sign(&position.z),
            Outcome::ProvenCoincident(_)
        ));
    }
    let accuracy = catallax_measure::Accuracy {
        relative: Rational::parse("1e-6").unwrap(),
        absolute: Rational::parse("1e-6").unwrap(),
    };
    let answer = catallax_measure::mass_properties(&model, body, &accuracy, &mut dec).unwrap();
    assert!(answer.marks.is_empty(), "{:?}", answer.marks);
    let pi = Ball::pi(128);
    let expected_volume = &pi * &Ball::from_rational(&r(8), 128);
    let expected_area = &(&pi * &Ball::from_rational(&r(8), 128))
        + &(&(&Ball::from_rational(&r(16), 128)
            * &Ball::from_rational(&r(2), 128).sqrt().unwrap())
            * &elliptic_e_half());
    assert!(
        answer.value.volume.lo <= expected_volume.lo()
            && answer.value.volume.hi >= expected_volume.hi(),
        "{:?} vs {:?}",
        answer.value.volume,
        expected_volume
    );
    assert!(
        answer.value.area.lo <= expected_area.lo() && answer.value.area.hi >= expected_area.hi(),
        "{:?} vs {:?}",
        answer.value.area,
        expected_area
    );
    for (value, expected) in answer
        .value
        .centroid
        .as_ref()
        .unwrap()
        .iter()
        .zip([1, 0, 1])
    {
        assert!(value.contains(&r(expected)));
    }
    let options = catallax_step::StepOptions::default();
    let mut step = Vec::new();
    let report = catallax_step::write_step(&model, &[body], &options, &mut dec, &mut step).unwrap();
    assert_eq!(report.value.faces, 3);
    assert_eq!(report.value.solids, 1);
    assert!(report.value.max_deviation <= options.curve_tolerance);
    let written = common::Step::parse(&step, &options.product_name);
    let extrusion = written.ids("SURFACE_OF_LINEAR_EXTRUSION");
    assert_eq!(extrusion.len(), 1);
    let surface = written.e(extrusion[0], "SURFACE_OF_LINEAR_EXTRUSION");
    let ellipse = written.e(common::reference(&surface.args[1]), "ELLIPSE");
    let (origin, normal, major) = written.placement(common::reference(&ellipse.args[1]));
    let vector = written.e(common::reference(&surface.args[2]), "VECTOR");
    let direction = written.direction(common::reference(&vector.args[1]));
    let magnitude = common::number(&vector.args[2]);
    assert!(magnitude.sign() == catallax_numbers::Sign::Positive);
    let direction = direction.map(|x| &x * &magnitude);
    let major_radius = common::number(&ellipse.args[2]);
    let minor_radius = common::number(&ellipse.args[3]);
    assert!(major_radius > minor_radius && minor_radius > r(0));
    let minor = std::array::from_fn(|i| {
        &normal[(i + 1) % 3] * &major[(i + 2) % 3] - &normal[(i + 2) % 3] * &major[(i + 1) % 3]
    });
    assert_eq!(common::dot(&normal, &major), r(0));
    let epsilon = &report.value.max_deviation / &major_radius;
    let lo = (&r(1) - &epsilon).pow(2).unwrap();
    let hi = (&r(1) + &epsilon).pow(2).unwrap();
    // Project each original swept-circle sample along the written generator
    // onto its written ellipse. Radial scaling gives a certified nearby point
    // on that STEP surface, without calling the kernel's chart or implicit.
    for (cosine, sine) in [
        (r(1), r(0)),
        (r(0), r(1)),
        (
            Rational::from_ratio(3, 5).unwrap(),
            Rational::from_ratio(4, 5).unwrap(),
        ),
        (r(-1), r(0)),
    ] {
        let radial = [cosine.mul_pow2(1), sine.mul_pow2(1)];
        for t in [r(0), Rational::from_ratio(1, 2).unwrap(), r(1)] {
            let shift = t.mul_pow2(1);
            let sample = [&radial[0] + &shift, radial[1].clone(), shift];
            let d = common::sub(&sample, &origin);
            let axial = common::dot(&d, &normal) / common::dot(&direction, &normal);
            let transverse = std::array::from_fn(|i| &d[i] - &(&direction[i] * &axial));
            let radius2 = common::sq(&common::dot(&transverse, &major))
                / (common::dot(&major, &major) * common::sq(&major_radius))
                + common::sq(&common::dot(&transverse, &minor))
                    / (common::dot(&minor, &minor) * common::sq(&minor_radius));
            assert!(
                lo <= radius2 && radius2 <= hi,
                "written extrusion misses the defining swept circle"
            );
        }
    }
    let circles = written.ids("CIRCLE");
    assert_eq!(circles.len(), 2);
    let mut cap_levels = Vec::new();
    for id in circles {
        let circle = written.e(id, "CIRCLE");
        let (center, normal, _) = written.placement(common::reference(&circle.args[1]));
        assert_eq!(common::number(&circle.args[2]), r(2));
        assert_eq!(normal, [r(0), r(0), r(1)]);
        assert_eq!(center[0], center[2]);
        assert_eq!(center[1], r(0));
        cap_levels.push(center[2].clone());
    }
    cap_levels.sort();
    assert_eq!(cap_levels, vec![r(0), r(2)]);
    let mut plane_levels = Vec::new();
    for id in written.ids("PLANE") {
        let plane = written.e(id, "PLANE");
        let (center, normal, _) = written.placement(common::reference(&plane.args[1]));
        let outward = if center[2] == r(0) {
            -1
        } else {
            assert_eq!(center[2], r(2));
            1
        };
        assert_eq!(normal, [r(0), r(0), r(outward)]);
        assert_eq!(center[0], center[2]);
        assert_eq!(center[1], r(0));
        plane_levels.push(center[2].clone());
    }
    plane_levels.sort();
    assert_eq!(plane_levels, vec![r(0), r(2)]);
}
