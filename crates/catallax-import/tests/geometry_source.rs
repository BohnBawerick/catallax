//! Source field retention through actual Part 21 bytes.
mod common;
use catallax_import::{
    ImportRefusal,
    geometry_source::{Curve, Geometry, KnotData, ParsedStep, Surface, TrimmingPreference},
};
use catallax_numbers::{Decider, Rational, with_operation_budget};

fn read(data: &str) -> ParsedStep {
    let text = format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('geometry'),'2;1');\nFILE_NAME('geometry.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n{data}\nENDSEC;\nEND-ISO-10303-21;\n"
    );
    ParsedStep::read(text.as_bytes(), &mut Decider::default()).unwrap()
}
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn geometry(s: &ParsedStep, id: u64) -> Geometry {
    s.geometry(id, &mut Decider::default()).unwrap()
}

#[test]
fn simple_and_complex_geometry_have_the_same_declared_fields() {
    let s=read("
#1=CARTESIAN_POINT('p',(1.,2.,3.));
#2=(CARTESIAN_POINT((1.,2.,3.)) GEOMETRIC_REPRESENTATION_ITEM() POINT() REPRESENTATION_ITEM('p'));
#3=CIRCLE('c',#80,2.);
#4=(CIRCLE(2.) CONIC(#80) CURVE() GEOMETRIC_REPRESENTATION_ITEM() REPRESENTATION_ITEM('c'));
#5=CONICAL_SURFACE('cone',#80,3.,0.4);
#6=(CONICAL_SURFACE(3.,0.4) ELEMENTARY_SURFACE(#80) GEOMETRIC_REPRESENTATION_ITEM() REPRESENTATION_ITEM('cone') SURFACE());
#7=AXIS2_PLACEMENT_3D('frame',#1,#20,#21);
#8=(AXIS2_PLACEMENT_3D(#20,#21) GEOMETRIC_REPRESENTATION_ITEM() PLACEMENT(#1) REPRESENTATION_ITEM('frame'));
#9=TRIMMED_CURVE('trim',#3,(PARAMETER_VALUE(0.),#30),(#31,PARAMETER_VALUE(1.)),.F.,.PARAMETER.);
#10=(BOUNDED_CURVE() CURVE() GEOMETRIC_REPRESENTATION_ITEM() REPRESENTATION_ITEM('trim') TRIMMED_CURVE(#3,(PARAMETER_VALUE(0.),#30),(#31,PARAMETER_VALUE(1.)),.F.,.PARAMETER.));
#20=DIRECTION('z',(0.,0.,1.));
#21=DIRECTION('x',(1.,0.,0.));
#30=CARTESIAN_POINT('start',(3.,2.,3.));
#31=CARTESIAN_POINT('end',(2.08060461173628,3.68294196961579,3.));
#80=AXIS2_PLACEMENT_3D('frame',#1,#20,#21);
");
    for (a, b) in [(1, 2), (3, 4), (5, 6), (7, 8), (9, 10)] {
        assert_eq!(geometry(&s, a), geometry(&s, b));
    }
    let Geometry::Curve(Curve::Trimmed(basis, trims, sense, preference)) = geometry(&s, 9) else {
        panic!("trim")
    };
    assert_eq!(basis, 3);
    assert!(!sense);
    assert_eq!(preference, TrimmingPreference::Parameter);
    assert_eq!(trims[0].len(), 2);
    assert_eq!(trims[1].len(), 2);
}

#[test]
fn legacy_transformation_names_do_not_shift_geometric_fields() {
    let s=read("
#1=CARTESIAN_TRANSFORMATION_OPERATOR_3D('name',#20,#21,#30,2.,#22);
#2=CARTESIAN_TRANSFORMATION_OPERATOR_3D('function','description','name',#20,#21,#30,2.,#22);
#3=(CARTESIAN_TRANSFORMATION_OPERATOR(#20,#21,#30,2.) CARTESIAN_TRANSFORMATION_OPERATOR_3D(#22) FUNCTIONALLY_DEFINED_TRANSFORMATION('function','description') GEOMETRIC_REPRESENTATION_ITEM() REPRESENTATION_ITEM('name'));
");
    assert_eq!(geometry(&s, 1), geometry(&s, 2));
    assert_eq!(geometry(&s, 1), geometry(&s, 3));
}

#[test]
fn profiles_offsets_replicas_and_surface_trims_keep_their_wrappers() {
    let s = read(
        "
#1=HYPERBOLA('hyperbola',#90,2.,3.);
#2=PARABOLA('parabola',#90,4.);
#3=OFFSET_CURVE_3D('offset',#1,-0.25,.U.,#91);
#4=CURVE_REPLICA('replica',#3,#92);
#5=SURFACE_OF_LINEAR_EXTRUSION('extrusion',#4,#93);
#6=SURFACE_OF_REVOLUTION('revolution',#2,#94);
#7=OFFSET_SURFACE('offset',#5,-1.,.F.);
#8=RECTANGULAR_TRIMMED_SURFACE('trim',#7,-2.,3.,4.,5.,.F.,.T.);
#9=SURFACE_REPLICA('replica',#8,#92);
#10=PCURVE('pcurve',#9,#95);
#11=SEAM_CURVE('seam',#4,(#10,#10),.PCURVE_S1.);
",
    );
    assert_eq!(
        geometry(&s, 1),
        Geometry::Curve(Curve::Hyperbola(90, q("2"), q("3")))
    );
    assert_eq!(
        geometry(&s, 3),
        Geometry::Curve(Curve::Offset3d(1, q("-0.25"), None, 91))
    );
    assert_eq!(
        geometry(&s, 5),
        Geometry::Surface(Surface::Extrusion(4, 93))
    );
    assert_eq!(
        geometry(&s, 8),
        Geometry::Surface(Surface::RectangularTrim(
            7,
            [q("-2"), q("3"), q("4"), q("5")],
            [false, true]
        ))
    );
    assert_eq!(
        geometry(&s, 11),
        Geometry::Curve(Curve::SurfaceCurve(
            4,
            vec![10, 10],
            "PCURVE_S1".into(),
            true
        ))
    );
}

#[test]
fn rational_splines_keep_control_ids_weights_flags_and_knot_units() {
    let s=read("
#1=(BOUNDED_CURVE() B_SPLINE_CURVE(2,(#10,#11,#12),.UNSPECIFIED.,.F.,.U.) B_SPLINE_CURVE_WITH_KNOTS((3,3),(-2.,5.),.UNSPECIFIED.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE((1.,2.,1.)) REPRESENTATION_ITEM('rational'));
#2=BEZIER_SURFACE('patch',1,1,((#10,#11),(#12,#13)),.UNSPECIFIED.,.F.,.F.,.U.);
");
    let Geometry::Curve(Curve::BSpline(c)) = geometry(&s, 1) else {
        panic!("spline")
    };
    assert_eq!(c.controls, vec![10, 11, 12]);
    assert_eq!(c.weights, Some(vec![q("1"), q("2"), q("1")]));
    assert_eq!(
        c.knots,
        KnotData::Explicit(vec![3, 3], vec![q("-2"), q("5")], "UNSPECIFIED".into())
    );
    assert_eq!(c.self_intersect, None);
    let Geometry::Surface(Surface::BSpline(s)) = geometry(&s, 2) else {
        panic!("surface")
    };
    assert_eq!(s.degree, [1, 1]);
    assert_eq!(s.controls, vec![vec![10, 11], vec![12, 13]]);
    assert_eq!(s.knots, [KnotData::Bezier, KnotData::Bezier]);
}

#[test]
fn certified_view_preserves_the_existing_box_reader() {
    let bytes = common::box_step(common::BoxFixture::new(["2", "3", "4"]));
    let s = ParsedStep::read(bytes.as_bytes(), &mut Decider::default()).unwrap();
    assert_eq!(
        s.certified_view(&Decider::default()).unwrap(),
        catallax_import::read_step(bytes.as_bytes()).unwrap()
    );
}

#[test]
fn source_byte_budget_is_checked_before_parsing() {
    let (result, work) =
        with_operation_budget(2, || ParsedStep::read(b"not STEP", &mut Decider::default()));
    assert!(work.exhausted);
    assert!(matches!(result, Err(ImportRefusal::ResourceLimit)));
}

#[test]
fn representation_contexts_retain_separate_unit_assignments() {
    let s=read("
#1=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#10,#11)) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#12)) REPRESENTATION_CONTEXT('model',''));
#2=GEOMETRIC_REPRESENTATION_CONTEXT('parameters','',2);
#3=GLOBAL_UNIT_ASSIGNED_CONTEXT('second model','',(#20,#21));
#4=CARTESIAN_POINT('point',(1.,2.,3.));
");
    let first = s.representation_context(1).unwrap();
    assert_eq!(first.dimension, Some(3));
    assert_eq!(first.units, vec![10, 11]);
    assert_eq!(first.uncertainties, vec![12]);
    assert_eq!(s.representation_context(2).unwrap().dimension, Some(2));
    assert_eq!(s.representation_context(3).unwrap().units, vec![20, 21]);
    assert!(matches!(
        s.representation_context(4),
        Err(ImportRefusal::Malformed { .. })
    ));
}

#[test]
fn topology_retains_void_orientation_and_surface_intent() {
    use catallax_import::geometry_source::Topology;
    let s=read("
#1=BREP_WITH_VOIDS('solid',#10,(#12));
#2=(BREP_WITH_VOIDS((#12)) FACETED_BREP() GEOMETRIC_REPRESENTATION_ITEM() MANIFOLD_SOLID_BREP(#10) REPRESENTATION_ITEM('solid') SOLID_MODEL());
#3=SHELL_BASED_SURFACE_MODEL('surface',(#10,#11));
#4=FACE_BASED_SURFACE_MODEL('faces',(#10));
#10=CLOSED_SHELL('outer',(#20));
#11=OPEN_SHELL('open',(#20));
#12=ORIENTED_CLOSED_SHELL('inner',*,#10,.F.);
#13=(CLOSED_SHELL() CONNECTED_FACE_SET((#20)) GEOMETRIC_REPRESENTATION_ITEM() REPRESENTATION_ITEM('outer') TOPOLOGICAL_REPRESENTATION_ITEM());
#20=ADVANCED_FACE('face',(#30),#90,.F.);
#21=(ADVANCED_FACE() FACE((#30)) FACE_SURFACE(#90,.F.) GEOMETRIC_REPRESENTATION_ITEM() REPRESENTATION_ITEM('face') TOPOLOGICAL_REPRESENTATION_ITEM());
#22=ORIENTED_FACE('face',*,#20,.F.);
#30=FACE_OUTER_BOUND('bound',#40,.T.);
#31=(FACE_BOUND(#40,.T.) FACE_OUTER_BOUND() REPRESENTATION_ITEM('bound') TOPOLOGICAL_REPRESENTATION_ITEM());
#40=EDGE_LOOP('loop',(#50));
#41=VERTEX_LOOP('loop',#70);
#42=POLY_LOOP('poly',(#80,#81,#82));
#50=ORIENTED_EDGE('edge',*,*,#60,.F.);
#60=EDGE_CURVE('curve',#70,#71,#91,.T.);
#61=(EDGE(#70,#71) EDGE_CURVE(#91,.T.) GEOMETRIC_REPRESENTATION_ITEM() REPRESENTATION_ITEM('curve') TOPOLOGICAL_REPRESENTATION_ITEM());
#70=VERTEX_POINT('vertex',#80);
");
    let t = |id| s.topology(id, &mut Decider::default()).unwrap();
    for (a, b) in [(1, 2), (10, 13), (20, 21), (30, 31), (60, 61)] {
        assert_eq!(t(a), t(b));
    }
    assert_eq!(
        t(1),
        Topology::Solid {
            outer: 10,
            voids: vec![12]
        }
    );
    assert_eq!(t(12), Topology::OrientedShell(10, false));
    assert_eq!(t(3), Topology::SurfaceModel(vec![10, 11]));
    assert_eq!(t(4), Topology::FaceModel(vec![10]));
    assert_eq!(
        t(11),
        Topology::Shell {
            faces: vec![20],
            closed: Some(false)
        }
    );
    assert_eq!(
        t(20),
        Topology::Face {
            bounds: vec![30],
            surface: Some(90),
            same_sense: Some(false)
        }
    );
    assert_eq!(t(22), Topology::OrientedFace(20, false));
    assert_eq!(t(41), Topology::VertexLoop(70));
    assert_eq!(t(42), Topology::PolyLoop(vec![80, 81, 82]));
    assert_eq!(t(50), Topology::OrientedEdge(60, false));
    assert_eq!(
        t(60),
        Topology::Edge {
            vertices: [70, 71],
            curve: Some(91),
            same_sense: Some(true)
        }
    );
    assert_eq!(t(70), Topology::Vertex(80));
}

#[test]
fn representation_maps_retain_both_frames_and_transform_direction() {
    use catallax_import::geometry_source::PlacementLink;
    let s=read("
#1=REPRESENTATION_MAP(#10,#11);
#2=MAPPED_ITEM('mapped',#1,#12);
#3=ITEM_DEFINED_TRANSFORMATION('placement',$,#20,#21);
#4=REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION('link',$,#11,#13,#3);
#5=(REPRESENTATION_RELATIONSHIP('link',$,#11,#13) REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#3) SHAPE_REPRESENTATION_RELATIONSHIP());
");
    let link = |id| s.placement_link(id, &mut Decider::default()).unwrap();
    assert_eq!(link(1), PlacementLink::RepresentationMap(10, 11));
    assert_eq!(link(2), PlacementLink::MappedItem(1, 12));
    assert_eq!(link(3), PlacementLink::ItemDefinedTransformation([20, 21]));
    assert_eq!(link(4), PlacementLink::Relationship([11, 13], 3));
    assert_eq!(link(4), link(5));
}

#[test]
fn exact_mapping_keeps_line_magnitude_and_polyline_index_parameters() {
    use catallax_shapes::Curve as NativeCurve;
    let s = read(
        "
#1=CARTESIAN_POINT('p',(1.,2.,3.));
#2=DIRECTION('direction',(0.,3.,4.));
#3=VECTOR('vector',#2,5.);
#4=LINE('line',#1,#3);
#5=CARTESIAN_POINT('q',(1.,2.,103.));
#6=CARTESIAN_POINT('r',(2.,2.,103.));
#7=POLYLINE('polyline',(#1,#5,#6));
#8=CARTESIAN_POINT('uv',(2.,3.));
#9=POLYLINE('uv',(#8,#8));
",
    );
    let mut dec = Decider::default();
    let NativeCurve::Line(line) = s.native_curve(4, &q("2"), &mut dec).unwrap() else {
        panic!("line")
    };
    assert_eq!(line.origin.y.rational_ref(), Some(&q("4")));
    for (actual, expected) in [(&line.dir.y, q("6")), (&line.dir.z, q("8"))] {
        assert!(matches!(
            dec.sign(&(actual - &catallax_numbers::Real::from_rational(expected))),
            catallax_numbers::Outcome::ProvenCoincident(_)
        ));
    }
    let poly = s
        .spline_curve(7, &[q("1"), q("1"), q("1")], &mut dec)
        .unwrap();
    assert_eq!(poly.knots(), &[q("0"), q("0"), q("1"), q("2"), q("2")]);
    let uv = s
        .spline_curve(9, &[q("10"), q("100"), q("1")], &mut dec)
        .unwrap();
    assert_eq!(uv.controls()[0], [q("20"), q("300"), q("0")]);
    assert_eq!(uv.controls()[0], uv.controls()[1]);
}

#[test]
fn native_mapping_preserves_wrappers_and_complete_torus_definition() {
    let s = read(
        "
#1=CARTESIAN_POINT('p',(0.,0.,0.));
#2=AXIS2_PLACEMENT_3D('frame',#1,$,$);
#3=TOROIDAL_SURFACE('spindle',#2,1.,2.);
#4=PLANE('plane',#2);
#5=OFFSET_SURFACE('offset',#4,2.,.F.);
#6=CIRCLE('circle',#2,1.);
#7=TRIMMED_CURVE('trim',#6,(PARAMETER_VALUE(0.)),(PARAMETER_VALUE(1.)),.T.,.PARAMETER.);
",
    );
    let mut dec = Decider::default();
    assert!(s.native_surface(4, &q("1"), &mut dec).is_ok());
    assert!(s.native_surface(3, &q("1"), &mut dec).is_err());
    assert!(s.native_surface(5, &q("1"), &mut dec).is_err());
    assert!(s.native_curve(7, &q("1"), &mut dec).is_err());
    assert_eq!(
        geometry(&s, 3),
        Geometry::Surface(Surface::Torus(2, q("1"), q("2")))
    );
}

#[test]
fn parameter_charts_evaluate_conics_sweeps_offsets_replicas_and_spindle_sheet() {
    use catallax_import::geometry_source::AngleUnit;
    use catallax_numbers::Ball;
    let s = read(
        "
#1=CARTESIAN_POINT('origin',(1.,2.,3.));
#2=AXIS2_PLACEMENT_3D('frame',#1,$,$);
#3=DIRECTION('z',(0.,0.,1.));
#4=HYPERBOLA('hyperbola',#2,2.,3.);
#5=PARABOLA('parabola',#2,2.);
#6=CIRCLE('circle',#2,2.);
#7=OFFSET_CURVE_3D('offset',#6,0.25,.F.,#3);
#8=VECTOR('extrusion',#3,7.);
#9=SURFACE_OF_LINEAR_EXTRUSION('extrusion',#4,#8);
#10=CONICAL_SURFACE('cone',#2,3.,0.4);
#11=TOROIDAL_SURFACE('spindle',#2,1.,2.);
#12=CYLINDRICAL_SURFACE('cylinder',#2,2.);
#13=OFFSET_SURFACE('offset',#12,0.25,.F.);
#14=DIRECTION('y',(0.,1.,0.));
#15=VECTOR('line direction',#14,1.);
#16=LINE('line',#1,#15);
#17=CARTESIAN_POINT('translation',(10.,20.,30.));
#18=CARTESIAN_TRANSFORMATION_OPERATOR_3D('transform',#14,$,#17,2.,#3);
#19=CURVE_REPLICA('replica',#16,#18);
#20=CARTESIAN_POINT('zero',(0.,0.,0.));
#21=AXIS2_PLACEMENT_3D('ellipse frame',#20,$,$);
#22=ELLIPSE('profile',#21,3.,2.);
#23=AXIS1_PLACEMENT('axis',#20,#14);
#24=SURFACE_OF_REVOLUTION('revolution',#22,#23);
#25=RECTANGULAR_TRIMMED_SURFACE('trim',#13,0.,1.,0.,2.,.F.,.T.);
#26=CURVE_REPLICA('cycle',#27,#18);
#27=CURVE_REPLICA('cycle',#26,#18);
#28=CARTESIAN_POINT('uv',(2.,3.));
#29=CARTESIAN_POINT('uv',(4.,5.));
#30=POLYLINE('uv',(#28,#29));
#31=DEFINITIONAL_REPRESENTATION('uv',(#30),#90);
#32=PCURVE('pcurve',#12,#31);
",
    );
    let cases = [
        (4, false, vec!["0.4"]),
        (5, false, vec!["3"]),
        (7, false, vec!["0.4"]),
        (9, true, vec!["0.4", "2"]),
        (10, true, vec!["0.3", "2"]),
        (11, true, vec!["0.3", "2.4"]),
        (13, true, vec!["0.3", "2"]),
        (19, false, vec!["2"]),
        (24, true, vec!["0.3", "0.4"]),
        (25, true, vec!["0.3", "2"]),
        (32, false, vec!["0.5"]),
    ];
    let mut observations = Vec::new();
    for (id, surface, parameters) in cases {
        let mut dec = Decider::default();
        let charts = if surface {
            s.surface_charts(id, &q("1"), AngleUnit::Radians, &mut dec)
        } else {
            s.curve_charts(id, &q("1"), AngleUnit::Radians, &mut dec)
        }
        .unwrap();
        assert_eq!(charts.len(), 1);
        let parameters = parameters
            .iter()
            .map(|v| Ball::from_rational(&q(v), 128))
            .collect::<Vec<_>>();
        let value = charts[0].evaluate(&parameters, &mut dec).unwrap();
        for coordinate in &value {
            assert!((coordinate.hi() - coordinate.lo()) < q("0.000000000001"));
        }
        if id == 19 {
            for (value, expected) in value.iter().zip(["3", "10.5", "16.5"]) {
                assert!(value.lo() <= q(expected) && value.hi() >= q(expected));
            }
        }
        if id == 25 {
            assert_eq!(charts[0].wrappers(), &[13, 25]);
        }
        if id == 32 {
            for (value, expected) in value.iter().zip(["3", "4", "0"]) {
                assert!(value.lo() <= q(expected) && value.hi() >= q(expected));
            }
        }
        observations.push(serde_json::json!({"entity":id,"bounds":value.map(|v|[v.lo().to_string(),v.hi().to_string()])}));
    }
    assert!(
        s.curve_charts(26, &q("1"), AngleUnit::Radians, &mut Decider::default())
            .is_err()
    );
    if let Ok(path) = std::env::var("CATALLAX_GEOMETRY_OBSERVATIONS") {
        std::fs::write(path, serde_json::to_vec_pretty(&observations).unwrap()).unwrap();
    }
}

#[test]
fn units_are_resolved_per_context_without_rounding_degree_conversion() {
    use catallax_import::geometry_source::AngleUnit;
    let s=read("
#1=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#10,#11)) REPRESENTATION_CONTEXT('mm',''));
#2=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#12,#13)) REPRESENTATION_CONTEXT('metres',''));
#3=(GLOBAL_UNIT_ASSIGNED_CONTEXT((#10,#12)) REPRESENTATION_CONTEXT('conflict',''));
#10=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));
#11=(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));
#12=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.));
#13=(CONVERSION_BASED_UNIT('DEGREE',#14) NAMED_UNIT(*) PLANE_ANGLE_UNIT());
#14=MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.0174532925199433),#11);
");
    let a = s.context_units(1, &mut Decider::default()).unwrap();
    let b = s.context_units(2, &mut Decider::default()).unwrap();
    assert_eq!(a.length_mm, Some(q("1")));
    assert_eq!(b.length_mm, Some(q("1000")));
    assert_eq!(a.plane_angle, Some(AngleUnit::ScaledRadians(q("1"))));
    assert_eq!(
        b.plane_angle,
        Some(AngleUnit::ScaledRadians(q("0.0174532925199433")))
    );
    assert!(s.context_units(3, &mut Decider::default()).is_err());
}

#[test]
fn geometry_sets_and_model_declarations_do_not_infer_closed_solids() {
    use catallax_import::geometry_source::{ModelDeclarationKind, Topology};
    let s=read("
#1=GEOMETRIC_SET('geometry',(#10,#11));
#2=GEOMETRIC_CURVE_SET('curves',(#11));
#3=(GEOMETRIC_CURVE_SET() GEOMETRIC_REPRESENTATION_ITEM() GEOMETRIC_SET((#11)) REPRESENTATION_ITEM('curves'));
#4=SHELL_BASED_SURFACE_MODEL('surface',(#20));
#5=BREP_WITH_VOIDS('solid',#21,(#22));
");
    assert_eq!(
        s.topology(1, &mut Decider::default()).unwrap(),
        Topology::GeometricSet {
            elements: vec![10, 11],
            curves_only: false
        }
    );
    assert_eq!(
        s.topology(2, &mut Decider::default()).unwrap(),
        s.topology(3, &mut Decider::default()).unwrap()
    );
    let models = s.model_declarations(&mut Decider::default()).unwrap();
    assert_eq!(
        models.iter().map(|m| m.entity).collect::<Vec<_>>(),
        vec![1, 2, 3, 4, 5]
    );
    assert_eq!(models[0].kind, ModelDeclarationKind::Geometry);
    assert_eq!(models[3].kind, ModelDeclarationKind::Shells);
    assert_eq!(models[4].kind, ModelDeclarationKind::Solid);
}

#[test]
fn part21_sections_retain_all_models_and_ignore_quoted_delimiters() {
    let text = include_str!("fixtures/surface-square.step");
    let model = text
        .lines()
        .find(|line| line.contains("=SHELL_BASED_SURFACE_MODEL("))
        .unwrap()
        .split_once('=')
        .unwrap()
        .1
        .trim_end_matches(';');
    for declaration in ["DATA;", "DATA('second',('AUTOMOTIVE_DESIGN'));", "data ;"] {
        let source = text.replace(
            "END-ISO-10303-21;",
            &format!("{declaration}\n#100000={model};\nENDSEC;\nEND-ISO-10303-21;"),
        );
        let file = catallax_import::read_step(source.as_bytes()).unwrap();
        assert_eq!(file.surface_models.len(), 2);
        let parsed = ParsedStep::read(source.as_bytes(), &mut Decider::default()).unwrap();
        assert_eq!(
            parsed
                .model_declarations(&mut Decider::default())
                .unwrap()
                .len(),
            2
        );
    }
    let source = text.replace("'surface-square.step'", "'DATA; ENDSEC;'");
    assert_eq!(
        catallax_import::read_step(source.as_bytes())
            .unwrap()
            .surface_models
            .len(),
        1
    );
    let duplicate = text.replace(
        "END-ISO-10303-21;",
        &format!("DATA;\n#32={model};\nENDSEC;\nEND-ISO-10303-21;"),
    );
    assert!(ParsedStep::read(duplicate.as_bytes(), &mut Decider::default()).is_err());
    assert!(ParsedStep::read(b"DATA; #1=POINT();", &mut Decider::default()).is_err());
}

#[test]
fn two_dimensional_replicas_refuse_unapplied_scale_and_second_axis() {
    use catallax_import::geometry_source::AngleUnit;
    use catallax_numbers::Ball;
    let s = read(
        "
#1=CARTESIAN_POINT('',(0.,0.));
#2=DIRECTION('',(1.,0.));
#3=VECTOR('',#2,1.);
#4=LINE('',#1,#3);
#5=CARTESIAN_POINT('',(2.,3.));
#6=CARTESIAN_TRANSFORMATION_OPERATOR_2D('',#2,$,#5,1.);
#7=DIRECTION('',(0.,-1.));
#8=CARTESIAN_TRANSFORMATION_OPERATOR_2D('',#2,$,#5,2.);
#9=CARTESIAN_TRANSFORMATION_OPERATOR_2D('',#2,#7,#5,1.);
#10=CURVE_REPLICA('',#4,#6);
#11=CURVE_REPLICA('',#4,#8);
#12=CURVE_REPLICA('',#4,#9);
",
    );
    let mut dec = Decider::default();
    let charts = s
        .curve_charts(10, &q("1"), AngleUnit::Radians, &mut dec)
        .unwrap();
    let point = charts[0]
        .evaluate(&[Ball::from_rational(&q("2"), 128)], &mut dec)
        .unwrap();
    for (value, expected) in point.iter().zip([q("4"), q("3"), q("0")]) {
        assert!(value.lo() <= expected && value.hi() >= expected);
    }
    for id in [11, 12] {
        let error = s
            .curve_charts(id, &q("1"), AngleUnit::Radians, &mut dec)
            .unwrap_err();
        assert!(
            matches!(error,catallax_import::geometry_source::MappingError::Geometry(r) if r.code=="import.geometry-parameter-mapping")
        );
    }
}

#[test]
fn unicode_whitespace_in_named_data_never_splits_a_character() {
    for whitespace in ["\u{2003}", "\u{00a0}", "\u{3000}"] {
        let text = format!(
            "{whitespace}DATA({whitespace}'section',({whitespace}'AP242')); #1=CARTESIAN_POINT('',(1.,2.,3.)); ENDSEC;"
        );
        let parsed = ParsedStep::read(text.as_bytes(), &mut Decider::default()).unwrap();
        assert_eq!(
            geometry(&parsed, 1),
            Geometry::Point(vec![q("1"), q("2"), q("3")])
        );
    }
}
