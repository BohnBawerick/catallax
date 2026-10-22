use crate::{Completion, Finding, NotComparable, Options, Pass, Reason, Report, Revision, Verdict};
use catallax_features::{Built, evaluate};
use catallax_field::{Direction, FaceState, Location, Plan, Workspace};
use catallax_numbers::{Ball, Decider, Interval, Rational};
use catallax_recipe::{Operation, Origin, Recipe};
use catallax_shapes::{RPoint3, RVector3, Surface};
use catallax_topology::Sense;
use std::collections::{BTreeMap, BTreeSet};

const BITS: u32 = 192;
fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn ball(x: &Rational) -> Ball {
    Ball::from_rational(x, BITS)
}
fn norm(v: &RVector3) -> Ball {
    ball(&v.norm_squared()).sqrt().unwrap()
}
fn absolute(x: Ball) -> Ball {
    if x.lo() >= r(0) {
        x
    } else if x.hi() <= r(0) {
        -x
    } else {
        ball(&r(0)).hull(&ball(&x.lo().abs().max(x.hi().abs())))
    }
}
fn support(surface: &Surface) -> catallax_numbers::MPoly {
    let polynomial = surface.implicit();
    let coefficient = polynomial.terms().into_iter().next().unwrap().1;
    polynomial.scale(&coefficient.recip().unwrap())
}
fn witness(
    surface: &Surface,
    polynomial: &catallax_numbers::MPoly,
    p: &RPoint3,
    width: &Rational,
    dec: &mut Decider,
) -> Result<RPoint3, String> {
    if polynomial
        .eval(&[p.x.clone(), p.y.clone(), p.z.clone()])
        .is_zero()
    {
        return Ok(p.clone());
    }
    catallax_boolean::surface_witness(surface, &p.to_real(), width, dec)
        .map_err(|e| format!("{}: {}", e.code, e.detail))?
        .ok_or_else(|| "no rational projection for this surface frame".into())
}
fn support_distance(surface: &Surface, p: &RPoint3) -> Ball {
    match surface {
        Surface::Plane(s) => absolute(
            ball(&p.minus(&s.origin).dot(&s.normal))
                .checked_div(&norm(&s.normal))
                .unwrap(),
        ),
        Surface::Sphere(s) => absolute(norm(&p.minus(&s.center)) - ball(&s.radius)),
        Surface::Cylinder(s) => {
            let d = p.minus(&s.origin);
            let h = d.dot(&s.axis) / s.axis.norm_squared();
            absolute(norm(&d.minus(&s.axis.scaled(&h))) - ball(&s.radius))
        }
        Surface::Cone(s) => {
            let d = p.minus(&s.apex);
            let t = d.dot(&s.axis) / s.axis.norm_squared();
            let radial = d.minus(&s.axis.scaled(&t));
            let height = ball(&t) * norm(&s.axis);
            absolute(norm(&radial) - ball(&s.slope) * absolute(height))
                .checked_div(&(ball(&r(1)) + ball(&(&s.slope * &s.slope))).sqrt().unwrap())
                .unwrap()
        }
        Surface::Torus(s) => {
            let d = p.minus(&s.center);
            let t = d.dot(&s.axis) / s.axis.norm_squared();
            let radial = d.minus(&s.axis.scaled(&t));
            let rho = norm(&radial) - ball(&s.major);
            let h = ball(&t) * norm(&s.axis);
            absolute((rho.square() + h.square()).sqrt().unwrap() - ball(&s.minor))
        }
        // Every distance is nonnegative. An unknown support supplies no tighter
        // lower bound; the missing projection still prevents agreement.
        _ => ball(&r(0)),
    }
}
fn direction(
    gradient: &[catallax_numbers::MPoly; 3],
    p: &RPoint3,
    sense: Sense,
) -> Option<[Ball; 3]> {
    let xyz = [p.x.clone(), p.y.clone(), p.z.clone()];
    let n: [Rational; 3] = std::array::from_fn(|i| gradient[i].eval(&xyz));
    let length = norm(&RVector3::new(n[0].clone(), n[1].clone(), n[2].clone()));
    if length.lo() <= r(0) {
        return None;
    }
    let sign = if sense == Sense::Same { 1 } else { -1 };
    let inverse = length.recip()?;
    Some(n.map(|x| &ball(&(x * r(sign))) * &inverse))
}
fn interval_ball(x: Interval) -> Option<Ball> {
    Some(ball(&Rational::from_f64_exact(x.lo())?).hull(&ball(&Rational::from_f64_exact(x.hi())?)))
}
fn coverage(recipe: &Recipe) -> Vec<NotComparable> {
    recipe.steps.iter().filter_map(|step| {
        let reason=match &step.op {
            Operation::Boolean {target,tool,..} if target==tool =>
                "B-rep recipe evaluation asserts distinct bodies when both operands name the same step",
            Operation::Block{..}|Operation::Cylinder{..}|Operation::Sphere{..}|
            Operation::Cone{..}|Operation::Torus{..}|Operation::Boolean{..}|
            Operation::Sketch{..}|Operation::Extrude{..}|Operation::ExtrudeThrough{..}|
            Operation::Revolve{..}=>return None,
            Operation::Path{..}=>"neither kernel evaluates recipe paths",
            Operation::Hole{..}=>"neither kernel evaluates recipe holes",
            Operation::Fillet{..}=>"neither kernel evaluates recipe fillets",
            Operation::Chamfer{..}=>"neither kernel evaluates recipe chamfers",
            Operation::Shell{..}=>"neither kernel evaluates recipe shelling",
            Operation::Offset{..}=>"neither kernel evaluates recipe offsets",
            Operation::Thicken{..}=>"neither kernel evaluates recipe thickening",
            Operation::Draft{..}=>"neither kernel evaluates recipe draft",
            Operation::Sweep{..}=>"neither kernel evaluates recipe sweeps",
            Operation::Loft{..}=>"neither kernel evaluates recipe lofts",
            Operation::Pattern{..}=>"neither kernel evaluates recipe patterns",
            Operation::Mirror{..}=>"neither kernel evaluates recipe mirroring",
            Operation::Lattice{..}=>"field lattice evaluation is not implemented; B-rep lattice evaluation is not planned",
            Operation::Gyroid{..}=>"field gyroid evaluation is not implemented; B-rep gyroid evaluation is not planned",
            Operation::SmoothBlend{..}=>"field smooth blending is not implemented; B-rep smooth blending is not planned",
            Operation::ImportStep{..}=>"neither kernel imports recipe STEP geometry",
            Operation::ImportMesh{..}=>"field mesh import is not implemented; B-rep mesh import is not planned",
            _=>"operation is not implemented by both kernel adapters",
        };
        Some(NotComparable {step:step.id.clone(),operation:step.op.name().into(),reason:reason.into()})
    }).collect()
}
fn unresolved(
    report: &mut Report,
    pass: Pass,
    p: &RPoint3,
    names: Vec<Vec<Origin>>,
    field: Vec<Origin>,
    reason: impl Into<String>,
) {
    report.findings.push(Finding {
        pass,
        point: p.clone(),
        brep_faces: names,
        field_faces: field,
        reason: Reason::Unresolved(reason.into()),
        gap: None,
    });
}
fn whole_failure(report: &mut Report, reason: impl Into<String>) {
    unresolved(
        report,
        Pass::Forward,
        &RPoint3::origin(),
        vec![],
        vec![],
        reason,
    );
}
fn distance_finding(
    report: &mut Report,
    pass: Pass,
    p: &RPoint3,
    names: Vec<Vec<Origin>>,
    field: Vec<Origin>,
    bounds: (Rational, Option<Rational>),
    options: &Options,
) {
    let (lo, hi) = bounds;
    if lo > options.accuracy {
        let low = Interval::from_rational(&lo).lo();
        let high = hi
            .as_ref()
            .map_or(f64::INFINITY, |h| Interval::from_rational(h).hi());
        report.findings.push(Finding {
            pass,
            point: p.clone(),
            brep_faces: names,
            field_faces: field,
            reason: Reason::Gap,
            gap: Some(Interval::new(low, high).unwrap()),
        });
    } else if hi.is_none_or(|h| h > options.accuracy) {
        unresolved(
            report,
            pass,
            p,
            names,
            field,
            "distance upper bound is unbounded or exceeds display accuracy",
        );
    }
}
pub(crate) fn run(
    revision: &Revision,
    recipe: &Recipe,
    built: Option<&Built>,
    options: &Options,
) -> Completion {
    if revision.is_cancelled() {
        return Completion::Cancelled;
    }
    let mut report = Report {
        revision: revision.id(),
        region: options.region.clone(),
        verdict: Verdict::Unresolved,
        forward_vertices: 0,
        reverse_samples: 0,
        findings: vec![],
        not_comparable: coverage(recipe),
    };
    let valid = options.accuracy > r(0)
        && (1..=64).contains(&options.divisions)
        && options.normal_cosine >= r(0)
        && options.normal_cosine < r(1)
        && options
            .region
            .as_ref()
            .is_none_or(|b| (0..3).all(|i| b.min[i] <= b.max[i]));
    if !valid {
        whole_failure(&mut report, "invalid comparison options");
        return revision.publish(report);
    }
    if !report.not_comparable.is_empty() {
        return revision.publish(report);
    }
    let field = match Plan::compile(recipe) {
        Ok(f) => f,
        Err(e) => {
            whole_failure(&mut report, format!("field build {}: {}", e.code, e.detail));
            return revision.publish(report);
        }
    };
    let supports: Vec<_> = field
        .faces()
        .iter()
        .map(|face| support(&face.surface))
        .collect();
    let registry: BTreeMap<_, _> = field
        .faces()
        .iter()
        .enumerate()
        .map(|(i, face)| {
            (
                face.origin.clone(),
                (field.select(&face.origin).unwrap(), i),
            )
        })
        .collect();
    let groups: Vec<Vec<_>> = supports
        .iter()
        .map(|source| {
            supports
                .iter()
                .enumerate()
                .filter(|(_, other)| source == *other)
                .map(|(i, _)| i)
                .collect()
        })
        .collect();
    let (mut model, result) = if let Some(b) = built {
        (b.model.clone(), b.result)
    } else {
        let built = match evaluate(recipe, &mut Decider::default()) {
            Ok(b) => b,
            Err(e) => {
                whole_failure(&mut report, format!("B-rep build: {e}"));
                return revision.publish(report);
            }
        };
        (built.model, built.result)
    };
    if revision.is_cancelled() {
        return Completion::Cancelled;
    }
    let mesh =
        match catallax_mesh::mesh_body(&model, result, &options.mesh, &mut Decider::default()) {
            Ok(m) => {
                if !m.marks.is_empty() {
                    whole_failure(&mut report, "mesh inherits marked geometry");
                }
                m.value
            }
            Err(e) => {
                whole_failure(&mut report, format!("mesh {}: {}", e.code, e.detail));
                return revision.publish(report);
            }
        };
    if model.body(result).kind == catallax_topology::BodyKind::Empty
        && field.is_proven_empty()
        && report.findings.is_empty()
    {
        report.verdict = Verdict::Agree;
        return revision.publish(report);
    }
    let face_origins: Vec<_> = mesh
        .faces
        .iter()
        .map(|f| model.face_origins(f.face).to_vec())
        .collect();
    let surfaces: Vec<_> = mesh
        .faces
        .iter()
        .map(|f| model.geom().surface(model.face(f.face).surface).clone())
        .collect();
    let polynomials: Vec<_> = surfaces.iter().map(Surface::implicit).collect();
    let gradients: Vec<[_; 3]> = polynomials
        .iter()
        .map(|p| std::array::from_fn(|i| p.partial(i)))
        .collect();
    let mut incident: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); mesh.positions.len()];
    for (i, face) in mesh.faces.iter().enumerate() {
        for triangle in &mesh.triangles[face.triangles.clone()] {
            for &v in &triangle.vertices {
                incident[v].insert(i);
            }
        }
    }
    let mut selected = BTreeSet::new();
    let mut sample_points = vec![None; mesh.faces.len()];
    let mut storage = Workspace::default();
    for (i, position) in mesh.positions.iter().enumerate() {
        if revision.is_cancelled() {
            return Completion::Cancelled;
        }
        let [x, y, z] =
            position.map(|x| Rational::from_f64_exact(x).expect("finite certified mesh"));
        let p = RPoint3::new(x, y, z);
        let error = incident[i]
            .iter()
            .map(|&f| &mesh.faces[f].max_chord_error)
            .max()
            .cloned()
            .unwrap_or_else(|| options.mesh.chord_tolerance.clone());
        if options
            .region
            .as_ref()
            .is_some_and(|b| !b.contains(&p, &error))
        {
            continue;
        }
        report.forward_vertices += 1;
        selected.extend(incident[i].iter().copied());
        for &face in &incident[i] {
            sample_points[face].get_or_insert_with(|| p.clone());
        }
        let query = field.proven(&p, &mut storage);
        if incident[i].is_empty() {
            unresolved(
                &mut report,
                Pass::Forward,
                &p,
                vec![],
                vec![],
                "mesh vertex has no incident face deviation certificate",
            );
            continue;
        }
        let d = query.distance.millimetres;
        let upper = Rational::from_f64_exact(d.hi()).map(|h| &h + &error);
        if query.location != Location::Indeterminate
            && upper.as_ref().is_some_and(|h| h <= &options.accuracy)
        {
            continue;
        }
        let names: Vec<Vec<Origin>> = incident[i]
            .iter()
            .map(|&f| face_origins[f].clone())
            .collect();
        let field_names: Vec<_> = names
            .iter()
            .flatten()
            .filter(|o| registry.contains_key(*o))
            .cloned()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        if query.location == Location::Indeterminate {
            unresolved(
                &mut report,
                Pass::Forward,
                &p,
                names,
                field_names,
                "field point classification is indeterminate",
            );
            continue;
        }
        let lo = Rational::from_f64_exact(d.lo());
        if let Some(lo) = lo {
            let lower = (&lo - &error).max(r(0));
            distance_finding(
                &mut report,
                Pass::Forward,
                &p,
                names,
                field_names,
                (lower, upper),
                options,
            );
        } else if d.lo() == f64::INFINITY && field.is_proven_empty() {
            report.findings.push(Finding {
                pass: Pass::Forward,
                point: p,
                brep_faces: names,
                field_faces: field_names,
                reason: Reason::Gap,
                gap: Some(d),
            });
        } else {
            unresolved(
                &mut report,
                Pass::Forward,
                &p,
                names,
                field_names,
                "field boundary distance is not finite",
            );
        }
    }
    #[cfg(test)]
    if let Some(barriers) = &revision.checkpoint {
        barriers.0.wait();
        barriers.1.wait();
    }
    if options.region.is_none() {
        selected.extend(0..mesh.faces.len());
    }
    let mut dec = Decider::default();
    let width = options.accuracy.mul_pow2(-20);
    let mut observed: BTreeMap<usize, BTreeSet<Origin>> = BTreeMap::new();
    let mut oriented = BTreeSet::new();
    let mut check_sample = |p: RPoint3, origin: &Origin| {
        if options
            .region
            .as_ref()
            .is_some_and(|b| !b.contains(&p, &r(0)))
        {
            return;
        }
        let (id, source) = registry[origin];
        let query = field
            .face_query(id, &p, &mut storage)
            .expect("selected field face");
        match query.state {
            FaceState::NotOnFace => return,
            FaceState::OnFace => (),
            _ => {
                unresolved(
                    &mut report,
                    Pass::Reverse,
                    &p,
                    vec![],
                    vec![origin.clone()],
                    "field trim membership is indeterminate",
                );
                return;
            }
        }
        report.reverse_samples += 1;
        let mut origins = BTreeSet::new();
        for &index in &groups[source] {
            let face = &field.faces()[index];
            let coactive = if index == source {
                query
            } else {
                field
                    .face_query(registry[&face.origin].0, &p, &mut storage)
                    .unwrap()
            };
            match coactive.state {
                FaceState::OnFace => {
                    origins.insert(face.origin.clone());
                }
                FaceState::NotOnFace => (),
                _ => {
                    unresolved(
                        &mut report,
                        Pass::Reverse,
                        &p,
                        vec![],
                        vec![face.origin.clone()],
                        "coactive field origin is indeterminate",
                    );
                }
            }
        }
        let origins: Vec<_> = origins.into_iter().collect();
        let candidates: Vec<usize> = mesh
            .faces
            .iter()
            .enumerate()
            .filter(|(i, _)| face_origins[*i].iter().any(|o| origins.contains(o)))
            .map(|(i, _)| i)
            .collect();
        let candidates = if candidates.is_empty() {
            (0..mesh.faces.len()).collect()
        } else {
            candidates
        };
        let names = candidates
            .iter()
            .map(|&i| face_origins[i].clone())
            .collect::<Vec<_>>();
        let mut lower: Option<Rational> = None;
        let mut closest: Option<(Rational, usize, RPoint3)> = None;
        let mut projection_errors = Vec::new();
        for &i in &candidates {
            if revision.is_cancelled() {
                return;
            }
            let face = mesh.faces[i].face;
            let surface = &surfaces[i];
            let d = support_distance(surface, &p).lo().max(r(0));
            lower = Some(lower.map_or(d.clone(), |old| old.min(d)));
            let q = match witness(surface, &polynomials[i], &p, &width, &mut dec) {
                Ok(q) => q,
                Err(e) => {
                    projection_errors.push(e);
                    continue;
                }
            };
            match catallax_boolean::face_contains(&mut model, face, &q, &mut dec) {
                Ok(true) => (),
                Ok(false) => continue,
                Err(e) => {
                    projection_errors.push(format!("{}: {}", e.code, e.detail));
                    continue;
                }
            }
            let upper = norm(&p.minus(&q)).hi();
            if closest.as_ref().is_none_or(|(best, _, _)| upper < *best) {
                closest = Some((upper, i, q));
            }
        }
        let upper = closest.as_ref().map(|(d, _, _)| d.clone());
        if let Some(lower) = lower {
            distance_finding(
                &mut report,
                Pass::Reverse,
                &p,
                names.clone(),
                origins.clone(),
                (lower, upper.clone()),
                options,
            );
        } else {
            report.findings.push(Finding {
                pass: Pass::Reverse,
                point: p.clone(),
                brep_faces: names.clone(),
                field_faces: origins.clone(),
                reason: Reason::Gap,
                gap: Some(Interval::point(f64::INFINITY).unwrap()),
            });
        }
        let Some((upper, i, q)) = closest else {
            if !projection_errors.is_empty() {
                unresolved(
                    &mut report,
                    Pass::Reverse,
                    &p,
                    names,
                    origins,
                    projection_errors.join("; "),
                );
            }
            return;
        };
        if upper > options.accuracy {
            return;
        }
        selected.insert(i);
        sample_points[i].get_or_insert_with(|| p.clone());
        let face = mesh.faces[i].face;
        let expected = face_origins[i].as_slice();
        if origins.iter().any(|o| !expected.contains(o)) {
            report.findings.push(Finding {
                pass: Pass::Reverse,
                point: p.clone(),
                brep_faces: vec![expected.to_vec()],
                field_faces: origins.clone(),
                reason: Reason::Origins,
                gap: None,
            });
        }
        observed.entry(i).or_default().extend(origins);
        let brep = model.face(face);
        if let (Direction::Unique(field_normal), Some(brep_normal)) =
            (query.direction, direction(&gradients[i], &q, brep.sense))
        {
            let mut dot = ball(&r(0));
            for j in 0..3 {
                let Some(n) = interval_ball(field_normal[j]) else {
                    unresolved(
                        &mut report,
                        Pass::Reverse,
                        &p,
                        vec![expected.to_vec()],
                        vec![origin.clone()],
                        "normal enclosure is unbounded",
                    );
                    return;
                };
                dot = dot + (&brep_normal[j] * &n);
            }
            if dot.hi() < options.normal_cosine {
                report.findings.push(Finding {
                    pass: Pass::Reverse,
                    point: p,
                    brep_faces: vec![expected.to_vec()],
                    field_faces: vec![origin.clone()],
                    reason: Reason::Direction,
                    gap: None,
                });
            } else if dot.lo() >= options.normal_cosine {
                oriented.insert(i);
            } else {
                unresolved(
                    &mut report,
                    Pass::Reverse,
                    &p,
                    vec![expected.to_vec()],
                    vec![origin.clone()],
                    "normal cosine enclosure crosses the comparison threshold",
                );
            }
        } else {
            unresolved(
                &mut report,
                Pass::Reverse,
                &p,
                vec![expected.to_vec()],
                vec![origin.clone()],
                "surface direction is not unique and certified",
            );
        }
    };
    crate::sample::visit(recipe, options.divisions, revision, &mut check_sample);
    let feature_sampling = field.visit_feature_candidates(options.divisions, |p, origin| {
        if revision.is_cancelled() {
            return false;
        }
        check_sample(p, origin);
        !revision.is_cancelled()
    });
    if let Err(error) = feature_sampling {
        whole_failure(
            &mut report,
            format!("field boundary sampling {}: {}", error.code, error.detail),
        );
    }
    if revision.is_cancelled() {
        return Completion::Cancelled;
    }
    for i in selected {
        let expected = face_origins[i].as_slice();
        let names = observed.get(&i);
        let point = sample_points[i].clone().unwrap_or_else(RPoint3::origin);
        if expected.iter().any(|o| !registry.contains_key(o)) {
            report.findings.push(Finding {
                pass: Pass::Reverse,
                point: point.clone(),
                brep_faces: vec![expected.to_vec()],
                field_faces: names.map_or_else(Vec::new, |n| n.iter().cloned().collect()),
                reason: if options.region.is_some() {
                    Reason::Unresolved("a merged-face origin absent from the field registry may belong outside the checked region".into())
                } else {
                    Reason::Origins
                },
                gap: None,
            });
        }
        let missing = options.region.is_none()
            && expected
                .iter()
                .any(|o| names.is_none_or(|n| !n.contains(o)));
        if expected.is_empty() || missing || !oriented.contains(&i) {
            unresolved(
                &mut report,
                Pass::Reverse,
                &point,
                vec![expected.to_vec()],
                vec![],
                "selected B-rep face lacks complete sampled origin/direction coverage",
            );
        }
    }
    if report.forward_vertices == 0 || report.reverse_samples == 0 {
        whole_failure(
            &mut report,
            "no two-way boundary samples in the requested scope",
        );
    }
    report.verdict = if report
        .findings
        .iter()
        .any(|f| !matches!(f.reason, Reason::Unresolved(_)))
    {
        Verdict::Disagree
    } else if report.findings.is_empty() {
        Verdict::Agree
    } else {
        Verdict::Unresolved
    };
    revision.publish(report)
}
