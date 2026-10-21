use super::*;
fn sense(v: Sense) -> &'static str {
    if v == Sense::Same { ".T." } else { ".F." }
}
fn bind<K: Ord + Copy>(map: &mut BTreeMap<K, usize>, k: K, id: usize) {
    if let Some(old) = map.insert(k, id) {
        assert_eq!(old, id, "shared model entity duplicated in STEP");
    }
}
fn coords(s: &Step, id: usize, p: &[Rational; 3], digits: u32) {
    for (text, exact) in s.xyz_text(id).iter().zip(p) {
        rounded(text, exact, digits);
    }
}
fn direction(s: &Step, id: usize, v: &RVector3, digits: u32) {
    let a = vec_xyz(v);
    let norm2 = dot(&a, &a);
    let texts = list(&s.e(id, "DIRECTION").args[1]);
    if norm2 == n(1) {
        for (text, exact) in texts.iter().zip(a) {
            rounded(text, &exact, digits);
        }
    } else {
        let length = Ball::from_rational(&norm2, 512).sqrt().unwrap();
        for (text, exact) in texts.iter().zip(a) {
            let q = number(text);
            // A written unit component must lie between directed decimal bounds
            // of the exact ratio; use outward rational bounds for sqrt(norm2).
            let a = &exact / &length.lo();
            let b = &exact / &length.hi();
            let lo = number(&decimal(&a.clone().min(b.clone()), digits, false));
            let hi = number(&decimal(&a.max(b), digits, true));
            assert!(lo <= q && q <= hi);
        }
    }
}
fn placement(
    s: &Step,
    id: usize,
    p: &RPoint3,
    axis: &RVector3,
    reference_dir: &RVector3,
    digits: u32,
) {
    let a = &s.e(id, "AXIS2_PLACEMENT_3D").args;
    coords(s, reference(&a[1]), &xyz(p), digits);
    direction(s, reference(&a[2]), axis, digits);
    direction(s, reference(&a[3]), reference_dir, digits);
}
fn surface(s: &Step, id: usize, def: &Surface, digits: u32) {
    let e = &s.entities[&id];
    let a = reference(&e.args[1]);
    match def {
        Surface::Plane(p) => {
            s.e(id, "PLANE");
            placement(s, a, &p.origin, &p.normal, &p.ref_dir, digits);
        }
        Surface::Cylinder(p) => {
            s.e(id, "CYLINDRICAL_SURFACE");
            placement(s, a, &p.origin, &p.axis, &p.ref_dir, digits);
            rounded(&e.args[2], &p.radius, digits);
        }
        Surface::Sphere(p) => {
            s.e(id, "SPHERICAL_SURFACE");
            placement(s, a, &p.center, &p.axis, &p.ref_dir, digits);
            rounded(&e.args[2], &p.radius, digits);
        }
        Surface::Torus(p) => {
            s.e(id, "TOROIDAL_SURFACE");
            placement(s, a, &p.center, &p.axis, &p.ref_dir, digits);
            rounded(&e.args[2], &p.major, digits);
            rounded(&e.args[3], &p.minor, digits);
        }
        Surface::Cone(p) => {
            s.e(id, "CONICAL_SURFACE");
            placement(s, a, &p.apex, &p.axis, &p.ref_dir, digits);
            assert_eq!(e.args[2], "0.");
            let angle = Ball::from_rational(&p.slope, 512).atan();
            let lo = number(&decimal(&angle.lo(), digits, false));
            let hi = number(&decimal(&angle.hi(), digits, true));
            let written = number(&e.args[3]);
            assert!(lo <= written && written <= hi);
        }
        _ => panic!("new surface kind needs a test oracle"),
    }
}
fn curve(s: &Step, id: usize, def: &Curve, digits: u32) {
    let e = &s.entities[&id];
    match def {
        Curve::Line(p) => {
            s.e(id, "LINE");
            coords(s, reference(&e.args[1]), &xyz(&p.origin), digits);
            let v = &s.e(reference(&e.args[2]), "VECTOR").args;
            assert_eq!(v[2], "1.");
            direction(s, reference(&v[1]), &p.dir, digits);
        }
        Curve::Circle(p) => {
            s.e(id, "CIRCLE");
            placement(
                s,
                reference(&e.args[1]),
                &p.center,
                &p.normal,
                &p.ref_dir,
                digits,
            );
            let root = Ball::from_rational(&p.radius_sq, 512).sqrt().unwrap();
            let q = number(&e.args[2]);
            assert!(
                q >= number(&decimal(&root.lo(), digits, false))
                    && q <= number(&decimal(&root.hi(), digits, true))
            );
        }
        Curve::Ellipse(p) => {
            s.e(id, "ELLIPSE");
            placement(
                s,
                reference(&e.args[1]),
                &p.center,
                &p.normal,
                &p.major_dir,
                digits,
            );
            for (text, square) in [(&e.args[2], &p.major_sq), (&e.args[3], &p.minor_sq)] {
                let root = Ball::from_rational(square, 512).sqrt().unwrap();
                let q = number(text);
                assert!(
                    q >= number(&decimal(&root.lo(), digits, false))
                        && q <= number(&decimal(&root.hi(), digits, true))
                );
            }
        }
        Curve::Intersection(_) => {} // Its independent support equations are tested in intersections.rs.
        _ => panic!("new curve kind needs a test oracle"),
    }
}
#[derive(Default)]
struct Maps {
    vertices: BTreeMap<VertexId, usize>,
    edges: BTreeMap<EdgeId, usize>,
    surfaces: BTreeMap<SurfaceId, usize>,
    curves: BTreeMap<CurveId, usize>,
}
fn vertex(
    m: &Model,
    s: &Step,
    model: VertexId,
    written: usize,
    o: &StepOptions,
    d: &Rational,
    maps: &mut Maps,
) {
    bind(&mut maps.vertices, model, written);
    let id = s.vertex(written);
    if let PointDef::Explicit(p) = m.geom().point(m.vertex(model).point) {
        coords(s, id, &xyz(p), o.significant_digits);
        assert!(
            distance_sq(&s.xyz(id), &xyz(p)) <= sq(d),
            "underreported vertex error"
        );
    } else {
        exact_point(s, id, &m.geom().position(m.vertex(model).point), o, d);
    }
}
fn exact_point(s: &Step, id: usize, p: &Point3, o: &StepOptions, d: &Rational) {
    let texts = s.xyz_text(id);
    let written = s.xyz(id);
    let mut error2 = n(0);
    for (i, value) in [&p.x, &p.y, &p.z].into_iter().enumerate() {
        let (lo, hi) = if let Some(q) = value.as_rational() {
            rounded(texts[i], &q, o.significant_digits);
            (q.clone(), q)
        } else {
            let ball = Decider::default().enclose(value, &r("1e-100")).unwrap();
            let lo = ball.lo();
            let hi = ball.hi();
            assert!(written[i] >= number(&decimal(&lo, o.significant_digits, false)));
            assert!(written[i] <= number(&decimal(&hi, o.significant_digits, true)));
            (lo, hi)
        };
        error2 = error2 + sq(&(&written[i] - &lo).abs().max((&written[i] - &hi).abs()));
    }
    assert!(error2 <= sq(d), "derived or added vertex outside deviation");
}

fn loopless_seams(
    s: &Step,
    bound: usize,
    surface_id: usize,
    def: &Surface,
    void: bool,
    o: &StepOptions,
    d: &Rational,
) {
    let bound = s.e(bound, "FACE_BOUND");
    assert_eq!(bound.args[2], if void { ".F." } else { ".T." });
    let loop_ = s.e(reference(&bound.args[1]), "EDGE_LOOP");
    let uses = list(&loop_.args[1])
        .into_iter()
        .map(|id| {
            let e = s.e(reference(id), "ORIENTED_EDGE");
            (reference(&e.args[3]), e.args[4].as_str())
        })
        .collect::<Vec<_>>();
    let (center, e3, e1, radii) = match def {
        Surface::Sphere(p) => (&p.center, &p.axis, &p.ref_dir, vec![p.radius.clone()]),
        Surface::Torus(p) => (
            &p.center,
            &p.axis,
            &p.ref_dir,
            vec![p.minor.clone(), &p.major + &p.minor],
        ),
        _ => panic!("unexpected loopless fixture"),
    };
    // Primitive fixtures use unit frame vectors, including after rational poses.
    assert_eq!(e3.norm_squared(), n(1));
    assert_eq!(e1.norm_squared(), n(1));
    let minus_e2 = e3.cross(e1).scaled(&n(-1));
    let mut edges = Vec::new();
    for (i, radius) in radii.iter().enumerate() {
        let edge_id = uses[i].0;
        let edge = s.e(edge_id, "EDGE_CURVE");
        let (origin, axis, start, end) = match def {
            Surface::Sphere(p) => (
                center.clone(),
                minus_e2.clone(),
                center.plus(&e3.scaled(&(-&p.radius))),
                center.plus(&e3.scaled(&p.radius)),
            ),
            Surface::Torus(p) => {
                let vertex = center.plus(&e1.scaled(&(&p.major + &p.minor)));
                let origin = if i == 0 {
                    center.plus(&e1.scaled(&p.major))
                } else {
                    center.clone()
                };
                (
                    origin,
                    if i == 0 { minus_e2.clone() } else { e3.clone() },
                    vertex.clone(),
                    vertex,
                )
            }
            _ => unreachable!(),
        };
        curve(
            s,
            reference(&edge.args[3]),
            &Curve::Circle(Circle {
                center: origin,
                normal: axis,
                ref_dir: e1.clone(),
                radius_sq: sq(radius),
            }),
            o.significant_digits,
        );
        for (vertex, point) in [(&edge.args[1], start), (&edge.args[2], end)] {
            let id = s.vertex(reference(vertex));
            coords(s, id, &xyz(&point), o.significant_digits);
            assert!(distance_sq(&s.xyz(id), &xyz(&point)) <= sq(d));
            on_written_surface(s, surface_id, &xyz(&point), d);
        }
        edges.push(edge_id);
    }
    let expected = edges
        .iter()
        .map(|&id| (id, ".T."))
        .chain(edges.iter().map(|&id| (id, ".F.")))
        .collect::<Vec<_>>();
    assert_eq!(uses, expected);
}
pub fn assert_topology(m: &Model, bs: &[BodyId], s: &Step, o: &StepOptions, d: &Rational) {
    let mut maps = Maps::default();
    let items = s.items();
    assert_eq!(items[0], 17);
    let mut solids = items[1..].iter();
    for &b in bs {
        for &lump in &m.body(b).lumps {
            let solid = *solids.next().unwrap();
            let shells = &m.lump(lump).shells;
            let e = s.e(
                solid,
                if shells.len() == 1 {
                    "MANIFOLD_SOLID_BREP"
                } else {
                    "BREP_WITH_VOIDS"
                },
            );
            let mut written_shells = vec![reference(&e.args[1])];
            if shells.len() > 1 {
                for oriented in list(&e.args[2]) {
                    let v = &s.e(reference(oriented), "ORIENTED_CLOSED_SHELL").args;
                    assert_eq!(v[1], "*");
                    assert_eq!(v[3], ".F.");
                    written_shells.push(reference(&v[2]));
                }
            }
            assert_eq!(written_shells.len(), shells.len());
            for (index, (&shell, &written)) in shells.iter().zip(&written_shells).enumerate() {
                let written_faces = list(&s.e(written, "CLOSED_SHELL").args[1]);
                let faces = &m.shell(shell).faces;
                assert_eq!(faces.len(), written_faces.len());
                for (&fid, face_ref) in faces.iter().zip(written_faces) {
                    let face = m.face(fid);
                    let f = s.e(reference(face_ref), "ADVANCED_FACE");
                    let surface_id = reference(&f.args[2]);
                    bind(&mut maps.surfaces, face.surface, surface_id);
                    surface(
                        s,
                        surface_id,
                        m.geom().surface(face.surface),
                        o.significant_digits,
                    );
                    let expected_sense = if index == 0 {
                        face.sense
                    } else {
                        flip(face.sense)
                    };
                    assert_eq!(f.args[3], sense(expected_sense));
                    let bounds = list(&f.args[1]);
                    if face.loops.is_empty() {
                        assert_eq!(bounds.len(), 1);
                        loopless_seams(
                            s,
                            reference(bounds[0]),
                            surface_id,
                            m.geom().surface(face.surface),
                            index != 0,
                            o,
                            d,
                        );
                        continue;
                    }
                    assert_eq!(bounds.len(), face.loops.len());
                    for (&lid, bound) in face.loops.iter().zip(bounds) {
                        let bound = s.e(reference(bound), "FACE_BOUND");
                        assert_eq!(bound.args[2], if index == 0 { ".T." } else { ".F." });
                        let loop_id = reference(&bound.args[1]);
                        match m.loop_(lid) {
                            Loop::Vertex(v) => {
                                let lv = reference(&s.e(loop_id, "VERTEX_LOOP").args[1]);
                                vertex(m, s, *v, lv, o, d, &mut maps);
                                if let PointDef::Explicit(p) = m.geom().point(m.vertex(*v).point) {
                                    on_written_surface(s, surface_id, &xyz(p), d);
                                }
                            }
                            Loop::Coedges(coedges) => {
                                let oes = list(&s.e(loop_id, "EDGE_LOOP").args[1]);
                                assert_eq!(oes.len(), coedges.len());
                                for (&cid, oe) in coedges.iter().zip(oes) {
                                    let co = m.coedge(cid);
                                    let oe = s.e(reference(oe), "ORIENTED_EDGE");
                                    assert_eq!(&oe.args[1..3], ["*", "*"]);
                                    assert_eq!(oe.args[4], sense(co.sense));
                                    let eid = reference(&oe.args[3]);
                                    bind(&mut maps.edges, co.edge, eid);
                                    let edge = m.edge(co.edge);
                                    let e = s.e(eid, "EDGE_CURVE");
                                    let curve_id = reference(&e.args[3]);
                                    let spline =
                                        s.entities[&curve_id].kind == "B_SPLINE_CURVE_WITH_KNOTS";
                                    if !spline {
                                        bind(&mut maps.curves, edge.curve, curve_id);
                                    }
                                    assert_eq!(
                                        e.args[4],
                                        if spline { ".T." } else { sense(edge.sense) }
                                    );
                                    curve(
                                        s,
                                        curve_id,
                                        m.geom().curve(edge.curve),
                                        o.significant_digits,
                                    );
                                    for (v, vr) in
                                        [(edge.start, &e.args[1]), (edge.end, &e.args[2])]
                                    {
                                        if let Some(v) = v {
                                            vertex(m, s, v, reference(vr), o, d, &mut maps);
                                            if let PointDef::Explicit(p) =
                                                m.geom().point(m.vertex(v).point)
                                            {
                                                on_written_surface(s, surface_id, &xyz(p), d);
                                            }
                                        }
                                    }
                                    if edge.start.is_none() {
                                        assert_eq!(e.args[1], e.args[2], "ring endpoints");
                                        exact_point(
                                            s,
                                            s.vertex(reference(&e.args[1])),
                                            &m.geom().curve(edge.curve).parameter_origin(m.geom()),
                                            o,
                                            d,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    assert!(solids.next().is_none());
    // Distinct model topology and definitions must not alias merely by coordinates.
    for values in [
        maps.vertices.values().copied().collect::<Vec<_>>(),
        maps.edges.values().copied().collect(),
        maps.surfaces.values().copied().collect(),
        maps.curves.values().copied().collect(),
    ] {
        assert_eq!(values.len(), values.iter().collect::<BTreeSet<_>>().len());
    }
    let mut supports = BTreeMap::new();
    for (&model_edge, &written_edge) in &maps.edges {
        let curve_id = reference(&s.e(written_edge, "EDGE_CURVE").args[3]);
        if s.entities[&curve_id].kind == "B_SPLINE_CURVE_WITH_KNOTS" {
            let Curve::Intersection(curve) = m.geom().curve(m.edge(model_edge).curve) else {
                panic!("analytic model curve was approximated");
            };
            assert!(supports.insert(curve_id, curve.surfaces).is_none());
        }
    }
    check_splines_by_curve(s, d, |curve, p, d| {
        for surface in supports[&curve] {
            on_exact_support(m.geom().surface(surface), p, d);
        }
    });
    let mut owned = BTreeSet::new();
    for id in s.ids("VERTEX_POINT") {
        assert!(
            owned.insert(s.vertex(id)),
            "vertices share coordinate entity"
        );
    }
    for id in s.ids("AXIS2_PLACEMENT_3D") {
        for child in refs(&s.entities[&id].text) {
            assert!(owned.insert(child), "placement children shared");
        }
    }
    for id in s.ids("LINE") {
        let e = &s.entities[&id];
        assert!(owned.insert(reference(&e.args[1])));
        let v = reference(&e.args[2]);
        assert!(owned.insert(v));
        assert!(owned.insert(reference(&s.e(v, "VECTOR").args[1])));
    }
}
