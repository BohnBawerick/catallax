//! Frozen projector controls for the public P3 plane/cylinder/torus family.
use catallax_numbers::{
    BudgetPolicy, Decider, Integer, Outcome, Rational, Real, Refusal, Sign, StrictSign,
    with_operation_budget,
};
use catallax_recipe::{Origin, Recipe};
use catallax_shapes::{Isometry, Point3, RPoint3, RVector3, Surface};
use catallax_topology::Sense;
use std::collections::{BTreeMap, BTreeSet};
fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn reference_sign(x: &Real, dec: &mut Decider, code: &'static str) -> Result<Sign, Refusal> {
    match dec.sign(x) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(e) => Err(Refusal::new(code, e.detail)),
        Outcome::Refused(e) => Err(e),
    }
}
fn integer_sqrt(n: &Integer) -> Option<Integer> {
    if n.sign() == Sign::Negative {
        return None;
    }
    if n.is_zero() {
        return Some(Integer::zero());
    }
    let mut x = Integer::one().mul_pow2(n.bit_length().div_ceil(2));
    loop {
        let q = n.div_rem_floor(&x)?.0;
        let y = (&x + &q).div_rem_floor(&Integer::from_i64(2))?.0;
        if y >= x {
            return (&x * &x == *n).then_some(x);
        }
        x = y;
    }
}

fn sqrt_rational(n: &Rational) -> Option<Rational> {
    Rational::new(integer_sqrt(n.numerator())?, integer_sqrt(n.denominator())?)
}
fn reference_probe(
    surface: &Surface,
    target: &Point3,
    width: &Rational,
    dec: &mut Decider,
) -> Result<Option<Point3>, Refusal> {
    match surface {
        Surface::Plane(s) => {
            let p = RPoint3::new(
                dec.enclose(&target.x, width)?.mid(),
                dec.enclose(&target.y, width)?.mid(),
                dec.enclose(&target.z, width)?.mid(),
            )
            .to_real();
            let displacement = p.minus(&s.origin);
            Ok(Some(p.plus(&s.normal.scaled(
                &(-displacement.dot(&s.normal) / s.normal.norm_squared()),
            ))))
        }
        Surface::Cylinder(s) => {
            let base = s.radial_reference();
            let tangent = s.axis.cross(&base);
            let squared = base.norm_squared();
            let tangent_squared = tangent.norm_squared();
            let d = target.minus(&s.origin);
            let h = d.dot(&s.axis) / s.axis.norm_squared();
            let along = d.dot(&base);
            // Keep the stereographic denominator away from the opposite pole,
            // including when the target is a rounded point just off the surface.
            let senses =
                if reference_sign(&along, dec, "boolean.imprint-ambiguous")? == Sign::Negative {
                    [Sense::Reversed, Sense::Same]
                } else {
                    [Sense::Same, Sense::Reversed]
                };
            for sense in senses {
                let b = if sense == Sense::Same {
                    base.clone()
                } else {
                    base.scaled(&Real::from_i64(-1))
                };
                let v = s.axis.cross(&b);
                let denominator = if sense == Sense::Same {
                    &squared + &along
                } else {
                    &squared - &along
                };
                match dec.sign(&denominator) {
                    Outcome::ProvenCoincident(_) | Outcome::NotExact(_) => continue,
                    Outcome::Refused(e) => return Err(e),
                    Outcome::Proven(_) => {}
                }
                let t = (d.dot(&v) / denominator) * (&squared / &tangent_squared);
                let t = Real::from(dec.enclose(&t, width)?.mid());
                let h = Real::from(dec.enclose(&h, width)?.mid());
                let z = &tangent_squared * &t.square();
                let denominator = &squared + &z;
                let radial = b
                    .scaled(&(&(&squared - &z) / &denominator))
                    .plus(&v.scaled(&((&squared * &t) * Real::from_i64(2) / denominator)));
                return Ok(Some(s.origin.plus(&radial).plus(&s.axis.scaled(&h))));
            }
            Ok(None)
        }
        _ => Ok(None),
    }
}
fn reference(
    surface: &Surface,
    target: &Point3,
    width: &Rational,
    dec: &mut Decider,
) -> Result<Option<RPoint3>, Refusal> {
    match surface {
        Surface::Torus(s) => {
            let (Some(center), Some(axis), Some(reference), Some(major), Some(minor)) = (
                s.center.as_rational(),
                s.axis.as_rational(),
                s.ref_dir.as_rational(),
                s.major.as_rational(),
                s.minor.as_rational(),
            ) else {
                return Ok(None);
            };
            let Some(axis_length) = sqrt_rational(&axis.norm_squared()) else {
                return Ok(None);
            };
            let axis = axis.scaled(&axis_length.recip().unwrap());
            let d = target.minus(&center.to_real());
            let Some(radial) = circle_direction(&axis, &reference, &d, width, dec)? else {
                return Ok(None);
            };
            let height = d.dot(&axis.to_real());
            let planar = d.minus(&axis.to_real().scaled(&height));
            let x = planar.norm_squared().sqrt() - Real::from(major.clone());
            let tube = catallax_shapes::Vector3::new(x, height, Real::zero());
            let z = catallax_shapes::RVector3::new(r(0), r(0), r(1));
            let x = catallax_shapes::RVector3::new(r(1), r(0), r(0));
            let Some(tube) = circle_direction(&z, &x, &tube, width, dec)? else {
                return Ok(None);
            };
            Ok(Some(
                center
                    .plus(&radial.scaled(&(&major + &(&minor * &tube.x))))
                    .plus(&axis.scaled(&(&minor * &tube.y))),
            ))
        }
        _ => Ok(reference_probe(surface, target, width, dec)?.and_then(|p| p.as_rational())),
    }
}
fn circle_direction(
    axis: &catallax_shapes::RVector3,
    reference: &catallax_shapes::RVector3,
    target: &catallax_shapes::Vector3,
    width: &Rational,
    dec: &mut Decider,
) -> Result<Option<catallax_shapes::RVector3>, Refusal> {
    let Some(a) = sqrt_rational(&axis.norm_squared()) else {
        return Ok(None);
    };
    let Some(b) = sqrt_rational(&reference.norm_squared()) else {
        return Ok(None);
    };
    let axis = axis.scaled(&a.recip().unwrap());
    let u = reference.scaled(&b.recip().unwrap());
    let v = axis.cross(&u);
    let x = target.dot(&u.to_real());
    let y = target.dot(&v.to_real());
    let length = (&x.square() + &y.square()).sqrt();
    for sign in [1, -1] {
        let denominator = &length + &(&x * &Real::from_i64(sign));
        if !matches!(dec.sign(&denominator), Outcome::Proven(_)) {
            continue;
        }
        let t = dec.enclose(&(&y / &denominator), width)?.mid();
        let squared = &t * &t;
        let denominator = &Rational::one() + &squared;
        let numerator = Rational::one() - squared;
        let cosine = (&numerator / &denominator) * r(sign);
        let sine = r(2) * t / denominator;
        return Ok(Some(u.scaled(&cosine).plus(&v.scaled(&sine))));
    }
    Ok(None)
}

fn coordinates(p: Option<RPoint3>) -> Option<[Rational; 3]> {
    p.map(|p| [p.x, p.y, p.z])
}
fn check(surface: &Surface, target: &RPoint3, width: &Rational) {
    let actual = catallax_boolean::surface_witness(
        surface,
        &target.to_real(),
        width,
        &mut Decider::default(),
    )
    .unwrap();
    let expected = reference(surface, &target.to_real(), width, &mut Decider::default()).unwrap();
    if let Some(p) = &actual {
        let polynomial = surface.implicit().unwrap();
        let polynomial = polynomial.rational().unwrap();
        assert!(
            polynomial
                .eval(&[p.x.clone(), p.y.clone(), p.z.clone()])
                .is_zero(),
            "not on carrier: {surface:?}"
        );
    }
    assert_eq!(
        coordinates(actual),
        coordinates(expected),
        "surface={surface:?}, target={target:?}"
    );
}
#[test]
fn p3_forward_witnesses_and_rigid_images_match_the_frozen_projector() {
    let recipe = Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    let (checked, work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let plan = catallax_field::Plan::compile_checked(&recipe, &mut dec).unwrap();
        let built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
        let options = catallax_compare::Options::default();
        let mesh = catallax_mesh::mesh_body(&built.model, built.result, &options.mesh, &mut dec)
            .unwrap()
            .value;
        let mut registry = BTreeMap::<Origin, Vec<usize>>::new();
        for (i, f) in plan.faces().iter().enumerate() {
            if plan.result_face(i).is_some() {
                registry.entry(f.origin.clone()).or_default().push(i);
            }
        }
        let mut incident = vec![BTreeSet::new(); mesh.positions.len()];
        for (i, f) in mesh.faces.iter().enumerate() {
            for t in &mesh.triangles[f.triangles.clone()] {
                for &v in &t.vertices {
                    incident[v].insert(i);
                }
            }
        }
        let reflection = Isometry::reflection(
            &RPoint3::new(r(1), r(2), r(3)),
            &RVector3::new(r(1), r(2), r(2)),
        )
        .unwrap();
        let mut checked = 0;
        for (i, p) in mesh.positions.iter().enumerate() {
            let [x, y, z] = p.map(|x| Rational::from_f64_exact(x).unwrap());
            let p = RPoint3::new(x, y, z);
            let error = incident[i]
                .iter()
                .map(|&f| &mesh.faces[f].max_chord_error)
                .max()
                .cloned()
                .unwrap();
            let width = (&options.accuracy - &error).mul_pow2(-1).mul_pow2(-20);
            let mut selected = BTreeSet::new();
            for &f in &incident[i] {
                for o in built.model.face_origins(mesh.faces[f].face) {
                    if let Some(ids) = registry.get(o) {
                        selected.extend(ids.iter().copied());
                    }
                }
            }
            for id in selected {
                let surface = &plan.faces()[id].surface;
                check(surface, &p, &width);
                checked += 1;
                if i % 97 == 0 {
                    let image = reflection.surface(surface);
                    let image_p = reflection.apply_point(&p.to_real()).as_rational().unwrap();
                    check(&image, &image_p, &width);
                    checked += 1;
                }
            }
        }
        checked
    });
    assert!(checked >= 10_098);
    assert!(!work.exhausted);
    eprintln!("frozen projection controls: {checked}, work={work:?}");
}
#[test]
fn closed_contact_offsets_keep_the_same_points_and_zero_work_refuses() {
    let recipe = Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    let plan = catallax_field::Plan::compile(&recipe).unwrap();
    let width = Rational::parse("1/419430400").unwrap();
    let step = Rational::parse("1/3200").unwrap();
    for surface in plan
        .faces()
        .iter()
        .enumerate()
        .filter(|(i, _)| plan.result_face(*i).is_some())
        .map(|(_, f)| &f.surface)
    {
        for p in [[28, 15, 12], [20, 15, 22], [1, 1, 10], [10, 15, 12]] {
            for z in -1..=1 {
                for y in -1..=1 {
                    for x in -1..=1 {
                        let point = RPoint3::new(
                            r(p[0]) + (&r(x) * &step),
                            r(p[1]) + (&r(y) * &step),
                            r(p[2]) + (&r(z) * &step),
                        );
                        check(surface, &point, &width);
                    }
                }
            }
        }
        let point = RPoint3::new(r(23), r(17), r(13));
        let (answer, work) = with_operation_budget(0, || {
            catallax_boolean::surface_witness(
                surface,
                &point.to_real(),
                &width,
                &mut Decider::default(),
            )
        });
        assert!(work.exhausted);
        assert!(answer.is_err());
        for bits in [1, 8, 32] {
            let policy = BudgetPolicy {
                max_precision_bits: bits,
                ..Default::default()
            };
            let mut a = Decider::new(policy);
            let mut b = Decider::new(policy);
            let a = catallax_boolean::surface_witness(surface, &point.to_real(), &width, &mut a)
                .map(coordinates)
                .map_err(|e| (e.code, e.detail));
            let b = reference(surface, &point.to_real(), &width, &mut b)
                .map(coordinates)
                .map_err(|e| (e.code, e.detail));
            assert_eq!(a, b);
        }
    }
}
