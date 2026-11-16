use catallax_features::fillet::numerical as fillet;
use catallax_numbers::{Decider, Rational, Real, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{BSplineCurve, BSplineSurface, Curve, Point3, RationalAtlas, SplineSupport, Surface};
use catallax_topology::{Sense, numerical::{self, NumericalBody, NumericalNativeTrim, NumericalOptions}};
use serde_json::json;

fn q(value: &str) -> Rational { Rational::parse(value).unwrap() }

fn linear_uv(edge: &numerical::NumericalEdge, trim: &numerical::NumericalTrim, dec: &mut Decider) -> BSplineCurve {
    let mut uv = trim.uv.clone();
    if trim.sense == Sense::Reversed { uv.reverse(); }
    let mut range = edge.parameters.clone();
    if range[0] > range[1] { range.reverse(); uv.reverse(); }
    BSplineCurve::new(1, vec![range[0].clone(), range[0].clone(), range[1].clone(), range[1].clone()],
        uv.map(|p| Point3::new(p[0].clone().into(), p[1].clone().into(), Real::zero())).to_vec(),
        vec![Real::one(); 2], dec).unwrap()
}

fn source(bore: bool, gap: bool, dec: &mut Decider) -> (NumericalBody, String) {
    let mut value = json!({"format":"catallax-recipe/5","name":"fillet admission controls","steps":[
        {"id":"stock","op":"block","size":["3","3","3"]}
    ],"result":"stock"});
    if bore {
        value["steps"].as_array_mut().unwrap().extend([
            json!({"id":"tool","op":"cylinder","radius":"1/2","height":"5","pose":{"rotation":["1","0","0","0"],"translation":["3/2","3/2","-1"]}}),
            json!({"id":"part","op":"boolean","kind":"subtract","target":"stock","tool":"tool"}),
        ]);
        value["result"] = json!("part");
    }
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let built = catallax_features::evaluate(&recipe, dec).unwrap();
    let tolerance = if gap { q("1/1000") } else { q("1/10000000") };
    let (mut definition, mut bindings) = catallax_step::numerical::definition_from_exact_model(
        &built.model, built.result, recipe.name.clone(), &tolerance, dec,
    ).unwrap();
    let selected = definition.edges.iter().position(|edge| {
        let [a,b] = edge.vertices.map(|i| &definition.vertices[i].position);
        a[0] == q("0") && b[0] == q("0") && a[1] == q("0") && b[1] == q("0") && a[2] != b[2]
    }).unwrap();
    let key = definition.edges[selected].key.clone();
    for face in &mut definition.faces {
        let side = face.loops.iter().flatten().any(|t| t.edge == selected);
        if side {
            let Surface::Plane(plane) = face.surface.clone() else { panic!("box sides are planar") };
            let points: Vec<_> = face.loops.iter().flatten().flat_map(|t| t.uv.iter()).collect();
            let bounds: [[Rational; 2]; 2] = std::array::from_fn(|axis| [
                points.iter().map(|p| p[axis].clone()).min().unwrap(),
                points.iter().map(|p| p[axis].clone()).max().unwrap(),
            ]);
            let vertical = plane.normal.cross(&plane.ref_dir);
            let plane = &plane;
            let vertical = &vertical;
            let controls = bounds[0].iter().flat_map(|u| bounds[1].iter().map(move |v| {
                plane.origin.plus(&plane.ref_dir.scaled(&Real::from(u.clone())))
                    .plus(&vertical.scaled(&Real::from(v.clone())))
            })).collect();
            let knots = bounds.each_ref().map(|b| vec![b[0].clone(),b[0].clone(),b[1].clone(),b[1].clone()]);
            let carrier = BSplineSurface::new([1,1],knots,[2,2],controls,vec![Real::one();4],dec).unwrap();
            let [atlas]: [_;1] = RationalAtlas::partition(carrier,[false;2],dec).unwrap().try_into().unwrap();
            face.surface = Surface::BSpline(SplineSupport::RationalAtlas(atlas));
        }
        if side || matches!(face.surface, Surface::Cylinder(_)) {
            bindings.retain(|binding| binding.face != face.key);
            for (boundary,row) in face.loops.iter().enumerate() {
                for (use_index,trim) in row.iter().enumerate() {
                    let edge = &definition.edges[trim.edge];
                    bindings.push(NumericalNativeTrim { face:face.key.clone(), boundary, use_index,
                        edge:edge.key.clone(), curve:linear_uv(edge,trim,dec) });
                }
            }
        }
    }
    if gap {
        let border = definition.edges.iter_mut().find(|edge| {
            let [a,b] = edge.vertices.map(|i| &definition.vertices[i].position);
            a[1] == q("0") && b[1] == q("0") && a[2] == q("0") && b[2] == q("0") && a[0] != b[0]
        }).unwrap();
        let Curve::Line(line) = &mut border.curve else { panic!("box border is linear") };
        line.origin.z = &line.origin.z + &Real::from(q("1/10000"));
    }
    let body = numerical::build_with_native_trims(definition,bindings,NumericalOptions{tolerance_mm:tolerance},dec).unwrap();
    (body,key)
}

#[test]
fn a_stricter_fillet_request_refuses_a_loose_source_border() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let (exact,key) = source(false,false,&mut dec);
        let loose = exact.rebuild(
            exact.definition().clone(), NumericalOptions { tolerance_mm:q("1/1000") }, &mut dec,
        ).unwrap();
        let tight = fillet::build(loose,&key,"round",q("1/10"),q("1/10000000"),&mut dec).unwrap();
        assert!(tight.body().exact_topology_faults().is_empty());
        assert_eq!(tight.body().options().tolerance_mm,q("1/10000000"));
        let (source,key) = source(false,true,&mut dec);
        let refused = fillet::build(source,&key,"round",q("1/10"),q("1/10000000"),&mut dec).unwrap_err();
        assert_eq!(refused.code,"fillet.numerical-body");
    });
    assert!(!work.exhausted,"{work:?}");
}

#[test]
fn rounding_an_outer_edge_preserves_both_untouched_cylinder_seam_uses() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let (source,key) = source(true,false,&mut dec);
        let seam = source.native_trims().iter().find(|a| source.native_trims().iter().any(|b|
            a.face == b.face && a.edge == b.edge && (a.boundary,a.use_index) != (b.boundary,b.use_index)
        )).unwrap();
        let face = seam.face.clone();
        let edge = seam.edge.clone();
        let sample = |body: &NumericalBody, dec: &mut Decider| {
            body.native_trims().iter().filter(|t| t.face == face && t.edge == edge).map(|t| {
                let range = t.curve.domain();
                let at = Real::from((range[0] + range[1]).mul_pow2(-1));
                let p = t.curve.point_at(&at, dec).unwrap();
                (t.boundary,t.use_index,p.x.as_rational().unwrap(),p.y.as_rational().unwrap())
            }).collect::<Vec<_>>()
        };
        let before = sample(&source,&mut dec);
        assert_eq!(before.len(),2);
        assert_ne!((&before[0].2,&before[0].3),(&before[1].2,&before[1].3));
        let rounded = fillet::build(source,&key,"round",q("1/10"),q("1/10000000"),&mut dec).unwrap();
        assert_eq!(sample(rounded.body(),&mut dec),before);
        assert!(rounded.body().exact_topology_faults().is_empty());
    });
    assert!(!work.exhausted,"{work:?}");
}
