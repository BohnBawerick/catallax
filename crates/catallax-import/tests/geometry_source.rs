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
