//! Rational supporting-surface candidates for analytic intersection cells.
//!
//! Exact halfspaces bound planar candidate polygons. Analytic rational charts
//! supply the curved candidates. Neither construction asserts surviving trim
//! membership, absence, nearest distance or complete geometric coverage.

use super::{Bounds, circle, fraction, placement_error};
use crate::{FaceId, Plan, geometry::Leaf};
use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::{RPoint3, RVector3, Surface, Vector3};
use std::collections::BTreeSet;

fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn rational(x: &Real) -> Result<Rational, Refusal> {
    x.as_rational().ok_or_else(placement_error)
}
fn vector(v: &Vector3) -> Result<RVector3, Refusal> {
    Ok(RVector3::new(
        rational(&v.x)?,
        rational(&v.y)?,
        rational(&v.z)?,
    ))
}
fn unit(v: &Vector3) -> Result<RVector3, Refusal> {
    let v = vector(v)?;
    let length = crate::profile::sqrt(&v.norm_squared())
        .filter(|n| !n.is_zero())
        .ok_or_else(placement_error)?;
    Ok(v.scaled(&(q(1) / length)))
}

/// An exact inequality n.p <= d; normalization never changes its material side.
#[derive(Clone)]
struct Halfspace {
    n: RVector3,
    d: Rational,
}
impl Halfspace {
    fn from_poly(poly: &catallax_numbers::RealPoly) -> Option<Self> {
        let poly = poly.rational()?;
        let mut n = [q(0), q(0), q(0)];
        let mut constant = q(0);
        for (m, c) in poly.terms() {
            match m.exponents().iter().sum::<u32>() {
                0 => constant = c,
                1 => n[m.exponents().iter().position(|&e| e == 1)?] = c,
                _ => return None,
            }
        }
        let [x, y, z] = n;
        Some(Self {
            n: RVector3::new(x, y, z),
            d: -constant,
        })
    }
    fn contains(&self, p: &RPoint3) -> bool {
        self.n.dot(&p.minus(&RPoint3::origin())) <= self.d
    }
}

struct Emitter<'a> {
    bounds: &'a Bounds,
    seen: BTreeSet<[Rational; 3]>,
    callback: &'a mut dyn FnMut(RPoint3) -> bool,
    dec: Decider,
}
impl Emitter<'_> {
    fn point(&mut self, p: RPoint3) -> Result<bool, Refusal> {
        self.dec.check_deadline()?;
        let xyz = [p.x.clone(), p.y.clone(), p.z.clone()];
        if xyz
            .iter()
            .zip(self.bounds)
            .any(|(x, [lo, hi])| x < lo || x > hi)
            || !self.seen.insert(xyz)
        {
            return Ok(true);
        }
        Ok((self.callback)(p))
    }
}

fn box_planes(bounds: &Bounds) -> Vec<Halfspace> {
    let mut planes = Vec::new();
    for (axis, [lo, hi]) in bounds.iter().enumerate() {
        for (sign, d) in [(-1, -lo), (1, hi.clone())] {
            let [x, y, z] = std::array::from_fn(|i| if i == axis { q(sign) } else { q(0) });
            planes.push(Halfspace {
                n: RVector3::new(x, y, z),
                d,
            });
        }
    }
    planes
}

fn planar(
    face: &Halfspace,
    constraints: &[Halfspace],
    seeds: &[RPoint3],
    bounds: &Bounds,
    divisions: u32,
    emit: &mut Emitter<'_>,
) -> Result<bool, Refusal> {
    let mut planes = constraints.to_vec();
    planes.extend(box_planes(bounds));
    let mut vertices = BTreeSet::new();
    for (i, a) in planes.iter().enumerate() {
        for b in &planes[i + 1..] {
            emit.dec.check_deadline()?;
            let determinant = face.n.dot(&a.n.cross(&b.n));
            if determinant.is_zero() {
                continue;
            }
            let v =
                a.n.cross(&b.n)
                    .scaled(&face.d)
                    .plus(&b.n.cross(&face.n).scaled(&a.d))
                    .plus(&face.n.cross(&a.n).scaled(&b.d))
                    .scaled(&(q(1) / determinant));
            let point = RPoint3::origin().plus(&v);
            if planes.iter().all(|plane| plane.contains(&point)) {
                vertices.insert([point.x, point.y, point.z]);
            }
        }
    }
    // A lower-dimensional linear section supplies no open face sample.
    // The comparator retains its ordinary missing-coverage outcome.
    if vertices.len() < 3 {
        return Ok(true);
    }
    let count =
        Rational::from_i64(i64::try_from(vertices.len()).expect("bounded candidate polygon"));
    let center: [Rational; 3] =
        std::array::from_fn(|i| &vertices.iter().fold(q(0), |sum, p| &sum + &p[i]) / &count);
    let point = |p: [Rational; 3]| {
        let [x, y, z] = p;
        RPoint3::new(x, y, z)
    };
    if !emit.point(point(center.clone()))? {
        return Ok(false);
    }
    // Project analytic reference points onto the plane. In particular, these
    // seed offset circular caps without guessing a centre from a display mesh.
    // A seed is still only a candidate and may fail the nonlinear trim query.
    for seed in seeds {
        emit.dec.check_deadline()?;
        let t = (&face.d - &face.n.dot(&seed.minus(&RPoint3::origin()))) / face.n.norm_squared();
        let p = seed.plus(&face.n.scaled(&t));
        let interior = planes.iter().all(|plane| {
            let value = plane.n.dot(&p.minus(&RPoint3::origin()));
            value < plane.d || (value == plane.d && plane.n.cross(&face.n).is_zero())
        });
        if interior && !emit.point(p)? {
            return Ok(false);
        }
    }
    // Strict convex combinations avoid polygon edges whenever the section
    // has relative interior. No polygon triangulation or boundary claim is used.
    for vertex in vertices {
        for i in 0..divisions {
            let t = fraction(i, divisions);
            let p = std::array::from_fn(|j| &center[j] + &(&(&vertex[j] - &center[j]) * &t));
            if !emit.point(point(p))? {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn angles(divisions: u32) -> Vec<[Rational; 2]> {
    let mut values = vec![[q(1), q(0)], [q(0), q(1)], [q(-1), q(0)], [q(0), q(-1)]];
    values.extend(circle(divisions));
    values
}

fn axial_samples(
    origin: &RPoint3,
    axis: &RVector3,
    bounds: &Bounds,
    anchors: &[Halfspace],
    divisions: u32,
    ray: bool,
) -> Vec<Rational> {
    let xyz = [&origin.x, &origin.y, &origin.z];
    let direction = [&axis.x, &axis.y, &axis.z];
    let mut range = [q(0), q(0)];
    for i in 0..3 {
        let a = &(&bounds[i][0] - xyz[i]) * direction[i];
        let b = &(&bounds[i][1] - xyz[i]) * direction[i];
        range[0] = &range[0] + &a.clone().min(b.clone());
        range[1] = &range[1] + &a.max(b);
    }
    if ray {
        range[0] = range[0].clone().max(q(0));
    }
    if range[0] >= range[1] {
        return Vec::new();
    }
    let mut knots = BTreeSet::from(range.clone());
    for plane in anchors {
        let denominator = plane.n.dot(axis);
        if !denominator.is_zero() {
            // These are only sampling seeds. Other operands' planes are never
            // imposed as material constraints on this cell or on the result.
            let t = (&plane.d - &plane.n.dot(&origin.minus(&RPoint3::origin()))) / denominator;
            if t > range[0] && t < range[1] {
                knots.insert(t);
            }
        }
    }
    let knots: Vec<_> = knots.into_iter().collect();
    let mut values = BTreeSet::new();
    for pair in knots.windows(2) {
        values.insert((&pair[0] + &pair[1]) / q(2));
        for i in 0..divisions {
            values.insert(&pair[0] + &((&pair[1] - &pair[0]) * fraction(i, divisions)));
        }
    }
    values.into_iter().collect()
}

fn curved(
    surface: &Surface,
    bounds: &Bounds,
    anchors: &[Halfspace],
    divisions: u32,
    emit: &mut Emitter<'_>,
) -> Result<bool, Refusal> {
    let angular = angles(divisions);
    if let Surface::Sphere(s) = surface {
        let center = s.center.as_rational().ok_or_else(placement_error)?;
        let radius = rational(&s.radius)?;
        for [c, s] in &angular {
            if c < &q(0) {
                continue;
            }
            for [x, y] in &angular {
                let p = center.plus(&RVector3::new(
                    &(&radius * c) * x,
                    &(&radius * c) * y,
                    &radius * s,
                ));
                if !emit.point(p)? {
                    return Ok(false);
                }
            }
        }
        return Ok(true);
    }
    let (origin, axis, reference, major, minor, slope, torus) = match surface {
        Surface::Cylinder(s) => (
            &s.origin,
            &s.axis,
            &s.ref_dir,
            rational(&s.radius())?,
            None,
            None,
            false,
        ),
        Surface::Cone(s) => (
            &s.apex,
            &s.axis,
            &s.ref_dir,
            q(0),
            None,
            Some(rational(&s.slope)?),
            false,
        ),
        Surface::EllipticCylinder(s) => (
            &s.origin,
            &s.axis,
            &s.major_dir,
            crate::profile::sqrt(&rational(&s.major_sq)?).ok_or_else(placement_error)?,
            Some(crate::profile::sqrt(&rational(&s.minor_sq)?).ok_or_else(placement_error)?),
            None,
            false,
        ),
        Surface::Torus(s) => (
            &s.center,
            &s.axis,
            &s.ref_dir,
            rational(&s.major)?,
            Some(rational(&s.minor)?),
            None,
            true,
        ),
        _ => {
            return Err(Refusal::new(
                "field.sample-unsupported",
                "cell support has no rational candidate chart",
            ));
        }
    };
    let origin = origin.as_rational().ok_or_else(placement_error)?;
    let axis = unit(axis)?;
    let u = unit(reference)?;
    let v = axis.cross(&u);
    if torus {
        let minor = minor.expect("torus minor radius");
        for [c, s] in &angular {
            let rho = &major + &(&minor * c);
            for [x, y] in &angular {
                let radial = u.scaled(x).plus(&v.scaled(y));
                let p = origin
                    .plus(&radial.scaled(&rho))
                    .plus(&axis.scaled(&(&minor * s)));
                if !emit.point(p)? {
                    return Ok(false);
                }
            }
        }
    } else {
        for t in axial_samples(&origin, &axis, bounds, anchors, divisions, slope.is_some()) {
            let radius = slope.as_ref().map_or_else(|| major.clone(), |s| s * &t);
            let other = minor.as_ref().unwrap_or(&radius);
            for [x, y] in &angular {
                let radial = u.scaled(&(&radius * x)).plus(&v.scaled(&(other * y)));
                let p = origin.plus(&axis.scaled(&t)).plus(&radial);
                if !emit.point(p)? {
                    return Ok(false);
                }
            }
        }
    }
    Ok(true)
}

pub(super) fn visit(
    plan: &Plan,
    divisions: u32,
    bounds: &Bounds,
    visit: &mut dyn FnMut(RPoint3, FaceId) -> bool,
) -> Result<bool, Refusal> {
    if !plan.faces.iter().any(|face| {
        face.in_result
            && plan.leaf_registers[face.leaf].is_some()
            && matches!(&plan.leaves[face.leaf], Leaf::Cell(cell) if cell.capped_chart().is_none())
    }) {
        return Ok(true);
    }
    let anchors: Vec<_> = plan
        .leaves
        .iter()
        .flat_map(|leaf| leaf.constraints())
        .filter_map(|c| Halfspace::from_poly(&c.polynomial))
        .collect();
    let seeds: Vec<_> = plan
        .leaves
        .iter()
        .flat_map(|leaf| leaf.constraints())
        .filter_map(|c| c.reference_point()?.as_rational())
        .collect();
    let conical = plan.faces.iter().any(|face| {
        face.in_result
            && plan.leaf_registers[face.leaf].is_some()
            && matches!(face.surface, Surface::Cone(_))
    });
    for (index, face) in plan.faces.iter().enumerate() {
        let Some(id) = plan.result_face(index) else {
            continue;
        };
        // Native drafted cap charts own their candidate traversal. This pass
        // supplements the other analytic conjunctions, without duplicating it.
        if !matches!(&plan.leaves[face.leaf], Leaf::Cell(cell) if cell.capped_chart().is_none()) {
            continue;
        }
        if conical && let Surface::Torus(torus) = &face.surface {
            if !crate::distance_samples::torus_candidates(
                plan,
                torus,
                bounds,
                divisions,
                &mut Decider::default(),
                &mut |p| visit(p, id),
            )? {
                return Ok(false);
            }
            continue;
        }
        if conical
            && let Surface::Cone(cone) = &face.surface
            && cone.apex.as_rational().is_none()
        {
            if !crate::distance_samples::cone_candidates(
                cone,
                bounds,
                &plan.leaves,
                &plan.leaf_registers,
                divisions,
                &mut Decider::default(),
                &mut |p| visit(p, id),
            )? {
                return Ok(false);
            }
            continue;
        }
        let mut callback = |p| visit(p, id);
        let mut emit = Emitter {
            bounds,
            seen: BTreeSet::new(),
            callback: &mut callback,
            dec: Decider::default(),
        };
        emit.dec.check_deadline()?;
        let complete = if matches!(face.surface, Surface::Plane(_)) {
            let plane = Halfspace::from_poly(&face.polynomial).ok_or_else(placement_error)?;
            let constraints: Vec<_> = plan.leaves[face.leaf]
                .constraints()
                .iter()
                .filter_map(|c| Halfspace::from_poly(&c.polynomial))
                .collect();
            planar(&plane, &constraints, &seeds, bounds, divisions, &mut emit)?
        } else {
            curved(&face.surface, bounds, &anchors, divisions, &mut emit)?
        };
        if !complete {
            return Ok(false);
        }
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::{Cone, Cylinder, EllipticCylinder, Sphere, Torus};

    fn vector(x: i64, y: i64, z: i64) -> Vector3 {
        RVector3::new(q(x), q(y), q(z)).to_real()
    }
    fn point(x: i64, y: i64, z: i64) -> RPoint3 {
        RPoint3::new(q(x), q(y), q(z))
    }

    #[test]
    fn oblique_planar_candidates_obey_exact_support_and_all_linear_trims() {
        let bounds = [[q(0), q(3)], [q(0), q(3)], [q(0), q(3)]];
        let plane = Halfspace {
            n: RVector3::new(q(1), q(1), q(1)),
            d: q(3),
        };
        let mut points = Vec::new();
        let mut callback = |p| {
            points.push(p);
            true
        };
        let mut emit = Emitter {
            bounds: &bounds,
            seen: BTreeSet::new(),
            callback: &mut callback,
            dec: Decider::default(),
        };
        assert!(
            planar(
                &plane,
                &[],
                &[point(3, 0, 0), point(1, 1, 1)],
                &bounds,
                4,
                &mut emit
            )
            .unwrap()
        );
        assert!(points.len() > 3);
        for p in points {
            assert_eq!(&(&p.x + &p.y) + &p.z, q(3));
            assert!([p.x, p.y, p.z].iter().all(|x| x > &q(0) && x < &q(3)));
        }
    }

    #[test]
    fn every_native_curved_chart_emits_exact_rational_support_points() {
        let mut dec = Decider::default();
        let surfaces = [
            Surface::Sphere(
                Sphere::new(
                    point(3, -2, 5).to_real(),
                    q(2).into(),
                    vector(0, 0, 2),
                    vector(2, 0, 0),
                    &mut dec,
                )
                .unwrap(),
            ),
            Surface::Cylinder(
                Cylinder::new(
                    point(3, -2, 5).to_real(),
                    vector(0, 0, 2),
                    q(2).into(),
                    vector(2, 0, 0),
                    &mut dec,
                )
                .unwrap(),
            ),
            Surface::Cone(
                Cone::new(
                    point(3, -2, -5).to_real(),
                    vector(0, 0, 2),
                    Rational::from_ratio(3, 2).unwrap().into(),
                    vector(2, 0, 0),
                    &mut dec,
                )
                .unwrap(),
            ),
            Surface::Torus(
                Torus::new(
                    point(1, 2, 3).to_real(),
                    vector(0, 0, 2),
                    q(3).into(),
                    q(1).into(),
                    vector(2, 0, 0),
                    &mut dec,
                )
                .unwrap(),
            ),
            Surface::EllipticCylinder(
                EllipticCylinder::new(
                    point(1, 2, 3).to_real(),
                    vector(0, 0, 3),
                    vector(2, 0, 0),
                    q(4).into(),
                    q(1).into(),
                    &mut dec,
                )
                .unwrap(),
            ),
        ];
        let bounds = std::array::from_fn(|_| [q(-10), q(10)]);
        for surface in surfaces {
            let mut points = Vec::new();
            let mut callback = |p| {
                points.push(p);
                true
            };
            let mut emit = Emitter {
                bounds: &bounds,
                seen: BTreeSet::new(),
                callback: &mut callback,
                dec: Decider::default(),
            };
            assert!(curved(&surface, &bounds, &[], 3, &mut emit).unwrap());
            assert!(!points.is_empty(), "{:?}", surface.kind());
            let implicit = surface.implicit().expect("analytic sample carrier");
            let poly = implicit.rational().unwrap();
            for p in points {
                assert!(
                    poly.eval(&[p.x, p.y, p.z]).is_zero(),
                    "{:?}",
                    surface.kind()
                );
            }
        }
    }

    #[test]
    fn cancellation_and_nonrational_charts_cannot_claim_complete_sampling() {
        let mut dec = Decider::default();
        let bounds = std::array::from_fn(|_| [q(-3), q(3)]);
        let sphere = Surface::Sphere(
            Sphere::new(
                point(0, 0, 0).to_real(),
                q(2).into(),
                vector(0, 0, 1),
                vector(1, 0, 0),
                &mut dec,
            )
            .unwrap(),
        );
        let mut calls = 0;
        let mut callback = |_| {
            calls += 1;
            false
        };
        let mut emit = Emitter {
            bounds: &bounds,
            seen: BTreeSet::new(),
            callback: &mut callback,
            dec: Decider::default(),
        };
        assert!(!curved(&sphere, &bounds, &[], 4, &mut emit).unwrap());
        assert_eq!(calls, 1);
        let cylinder = Surface::Cylinder(
            Cylinder::from_radius_sq(
                point(0, 0, 0).to_real(),
                vector(0, 0, 1),
                q(2).into(),
                vector(1, 0, 0),
                &mut dec,
            )
            .unwrap(),
        );
        let mut callback = |_| panic!("an irrational chart must not round a candidate");
        let mut emit = Emitter {
            bounds: &bounds,
            seen: BTreeSet::new(),
            callback: &mut callback,
            dec: Decider::default(),
        };
        assert_eq!(
            curved(&cylinder, &bounds, &[], 4, &mut emit)
                .unwrap_err()
                .code,
            "field.sample-algebraic-placement"
        );
    }

    fn drilled(mirrored: bool) -> Plan {
        let mut value = serde_json::json!({"format":"catallax-recipe/3","name":"candidate cells","steps":[
            {"id":"stock","op":"block","size":["10","10","10"],"pose":{"rotation":["1","0","0","0"],"translation":["-5","-5","0"]}},
            {"id":"hole","op":"hole","target":"stock","at":[{"id":"centre","point":["2","0","0"],"direction":["0","0","1"]}],"diameter":"2","depth":{"blind":"1"},"tip":"flat"}
        ],"result":"hole"});
        if mirrored {
            value["steps"].as_array_mut().unwrap().push(serde_json::json!({"id":"image","op":"mirror","body":"hole","plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}));
            value["result"] = "image".into();
        }
        Plan::compile(&catallax_recipe::Recipe::parse(&value.to_string()).unwrap()).unwrap()
    }

    #[test]
    fn real_cell_visitor_keeps_trim_queries_and_transformed_origins_authoritative() {
        for mirrored in [false, true] {
            let plan = drilled(mirrored);
            let mut samples = Vec::new();
            assert!(
                plan.visit_candidates(4, |p, id| {
                    samples.push((p, id));
                    true
                })
                .unwrap()
            );
            assert!(!samples.is_empty());
            let mut workspace = crate::Workspace::default();
            let mut surviving = 0;
            let mut discarded = 0;
            let mut floor = false;
            let mut wall = false;
            for (p, id) in samples {
                let origin = &plan.faces()[id.index()].origin;
                if mirrored {
                    assert!(matches!(origin.role, catallax_recipe::Role::Copy { .. }));
                }
                match plan.face_query(id, &p, &mut workspace).unwrap().state {
                    crate::FaceState::OnFace => {
                        surviving += 1;
                        let mut role = &origin.role;
                        while let catallax_recipe::Role::Copy { of, .. } = role {
                            role = &of.role;
                        }
                        if let catallax_recipe::Role::Hole { part, .. } = role {
                            floor |= *part == catallax_recipe::HolePart::Floor;
                            wall |= *part == catallax_recipe::HolePart::Wall;
                        }
                    }
                    crate::FaceState::NotOnFace => discarded += 1,
                    other => panic!("unexpected face state {other:?}"),
                }
            }
            assert!(surviving > 0);
            assert!(
                floor && wall,
                "both the offset planar cap and wall need surviving samples"
            );
            assert!(discarded > 0, "candidate membership must still be checked");
            let mut calls = 0;
            assert!(
                !plan
                    .visit_candidates(4, |_, _| {
                        calls += 1;
                        false
                    })
                    .unwrap()
            );
            assert_eq!(calls, 1);
        }
    }

    #[test]
    fn zero_allowance_refuses_even_when_the_cell_was_compiled_beforehand() {
        let plan = drilled(false);
        let (result, work) =
            catallax_numbers::with_operation_budget(0, || plan.visit_candidates(4, |_, _| true));
        assert!(work.exhausted);
        assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    }

    #[test]
    fn noncell_feature_order_and_cell_stock_exclusion_remain_unchanged() {
        let recipe = catallax_recipe::Recipe::parse(&serde_json::json!({
            "format":"catallax-recipe/3","name":"unchanged feature candidates","steps":[
                {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[
                    {"kind":"line","id":"bottom","to":["4","0"]},
                    {"kind":"line","id":"right","to":["4","3"]},
                    {"kind":"line","id":"top","to":["0","3"]},
                    {"kind":"line","id":"left","to":["0","0"]}
                ]]},
                {"id":"body","op":"extrude","profile":"profile","distance":"2"}
            ],"result":"body"
        }).to_string()).unwrap();
        let plan = Plan::compile(&recipe).unwrap();
        let mut samples = Vec::new();
        assert!(
            plan.visit_candidates(4, |p, id| {
                samples.push((p, id));
                true
            })
            .unwrap()
        );
        assert_eq!(samples.len(), 96);
        let first = &samples[0];
        assert_eq!(
            first.0,
            RPoint3::new(
                Rational::from_ratio(1, 2).unwrap(),
                q(0),
                Rational::from_ratio(1, 4).unwrap()
            )
        );
        assert!(
            matches!(&plan.faces()[first.1.index()].origin.role, catallax_recipe::Role::Side{segment,..} if segment.as_str()=="bottom")
        );
        let mut calls = 0;
        assert!(
            !plan
                .visit_candidates(4, |p, id| {
                    assert_eq!((&p, id), (&first.0, first.1));
                    calls += 1;
                    false
                })
                .unwrap()
        );
        assert_eq!(calls, 1);
        let stock = catallax_recipe::Recipe::parse(r#"{"format":"catallax-recipe/3","name":"stock","steps":[{"id":"body","op":"sphere","radius":"2","pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"body"}"#).unwrap();
        let plan = Plan::compile(&stock).unwrap();
        let bounds = [[q(-2), q(2)], [q(-2), q(2)], [q(-2), q(2)]];
        assert!(
            visit(&plan, 4, &bounds, &mut |_, _| {
                panic!("the analytic-cell pass must not duplicate stock sampling")
            })
            .unwrap()
        );
        let mut calls = 0;
        let mut workspace = crate::Workspace::default();
        assert!(
            plan.visit_candidates(4, |p, id| {
                assert_eq!(plan.result_face(id.index()), Some(id));
                assert_eq!(
                    plan.face_query(id, &p, &mut workspace).unwrap().state,
                    crate::FaceState::OnFace
                );
                calls += 1;
                true
            })
            .unwrap()
        );
        assert!(calls > 0);
    }
}
