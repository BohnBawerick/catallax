use crate::{Completion, Finding, NotComparable, Options, Pass, Reason, Report, Revision, Verdict};
use catallax_features::{Built, evaluate};
use catallax_field::{Direction, FaceState, Location, Plan, Workspace};
use catallax_numbers::{
    Ball, Decider, Interval, Outcome, Rational, Real, RealPoly, Refusal, Sign, StrictSign,
};
use catallax_recipe::{Operation, Origin, Recipe};
use catallax_shapes::{Point3, RPoint3, RVector3, Surface, Vector3};
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
fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    if let Some(q) = value.rational_ref() {
        return Ok(q.sign());
    }
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::Refused(e) => Err(e),
        Outcome::NotExact(e) => Err(Refusal::new(e.code, e.detail)),
    }
}
fn real_ball(value: &Real, dec: &mut Decider) -> Result<Ball, Refusal> {
    if let Some(q) = value.rational_ref() {
        Ok(ball(q))
    } else {
        dec.enclose_budgeted_at(value, BITS)
    }
}
fn real_norm(vector: &Vector3, dec: &mut Decider) -> Result<Ball, Refusal> {
    Ok(real_ball(&vector.norm_squared(), dec)?
        .sqrt()
        .expect("a squared norm is nonnegative"))
}
fn divide(a: Ball, b: Ball) -> Result<Ball, Refusal> {
    a.checked_div(&b).ok_or_else(|| {
        Refusal::new(
            "numbers.divisor-not-separated",
            "comparison support norm is not separated from zero",
        )
    })
}
fn support(surface: &Surface, dec: &mut Decider) -> Result<RealPoly, Refusal> {
    let polynomial = surface.implicit().ok_or_else(|| {
        Refusal::new(
            "compare.parametric-support",
            "support comparison needs a finite parametric correspondence solver",
        )
    })?;
    for (_, coefficient) in polynomial.terms() {
        if sign(&coefficient, dec)? != Sign::Zero {
            return Ok(polynomial.scale(&(&Real::one() / &coefficient)));
        }
    }
    Err(Refusal::new(
        "numbers.zero-polynomial",
        "surface support vanishes identically",
    ))
}
fn witness(
    surface: &Surface,
    polynomial: &RealPoly,
    p: &RPoint3,
    width: &Rational,
    dec: &mut Decider,
) -> Result<RPoint3, String> {
    let xyz = [p.x.clone(), p.y.clone(), p.z.clone()];
    let on_support = if let Some(q) = polynomial.rational() {
        q.eval(&xyz).is_zero()
    } else {
        sign(&polynomial.eval_real(&xyz.map(Real::from)), dec)
            .map_err(|e| format!("{}: {}", e.code, e.detail))?
            == Sign::Zero
    };
    if on_support {
        return Ok(p.clone());
    }
    catallax_boolean::surface_witness(surface, &p.to_real(), width, dec)
        .map_err(|e| format!("{}: {}", e.code, e.detail))?
        .ok_or_else(|| "no rational projection for this surface frame".into())
}
fn support_distance(surface: &Surface, p: &RPoint3, dec: &mut Decider) -> Result<Ball, Refusal> {
    let p = p.to_real();
    Ok(match surface {
        Surface::Plane(s) => absolute(divide(
            real_ball(&p.minus(&s.origin).dot(&s.normal), dec)?,
            real_norm(&s.normal, dec)?,
        )?),
        Surface::Sphere(s) => {
            absolute(real_norm(&p.minus(&s.center), dec)? - real_ball(&s.radius, dec)?)
        }
        Surface::Cylinder(s) => {
            let d = p.minus(&s.origin);
            let h = d.dot(&s.axis) / s.axis.norm_squared();
            absolute(real_norm(&d.minus(&s.axis.scaled(&h)), dec)? - real_ball(&s.radius(), dec)?)
        }
        Surface::Cone(s) => {
            let d = p.minus(&s.apex);
            let t = d.dot(&s.axis) / s.axis.norm_squared();
            let radial = d.minus(&s.axis.scaled(&t));
            let height = real_ball(&t, dec)? * real_norm(&s.axis, dec)?;
            let residual =
                absolute(real_norm(&radial, dec)? - real_ball(&s.slope, dec)? * absolute(height));
            let length = (ball(&r(1)) + real_ball(&(&s.slope * &s.slope), dec)?)
                .sqrt()
                .expect("one plus slope squared is positive");
            divide(residual, length)?
        }
        Surface::Torus(s) => {
            let d = p.minus(&s.center);
            let t = d.dot(&s.axis) / s.axis.norm_squared();
            let radial = d.minus(&s.axis.scaled(&t));
            let rho = real_norm(&radial, dec)? - real_ball(&s.major, dec)?;
            let h = real_ball(&t, dec)? * real_norm(&s.axis, dec)?;
            absolute(
                (rho.square() + h.square())
                    .sqrt()
                    .expect("a squared norm is nonnegative")
                    - real_ball(&s.minor, dec)?,
            )
        }
        // An unknown support supplies only the universal nonnegative distance
        // lower bound. A missing projection still prevents sampled agreement.
        _ => ball(&r(0)),
    })
}
fn direction(
    gradient: &[RealPoly; 3],
    p: &Point3,
    sense: Sense,
    dec: &mut Decider,
) -> Result<Option<[Ball; 3]>, Refusal> {
    let xyz = [p.x.clone(), p.y.clone(), p.z.clone()];
    let orientation = if sense == Sense::Same { 1 } else { -1 };
    if let [Some(x), Some(y), Some(z)] = gradient.each_ref().map(RealPoly::rational)
        && let Some(point) = p.as_rational()
    {
        let xyz = [point.x, point.y, point.z];
        let n = [x.eval(&xyz), y.eval(&xyz), z.eval(&xyz)];
        let length = norm(&RVector3::new(n[0].clone(), n[1].clone(), n[2].clone()));
        if length.lo() <= r(0) {
            return Ok(None);
        }
        let Some(inverse) = length.recip() else {
            return Ok(None);
        };
        return Ok(Some(n.map(|x| &ball(&(x * r(orientation))) * &inverse)));
    }
    let n = gradient.each_ref().map(|g| g.eval_real(&xyz));
    let length = real_norm(&Vector3::new(n[0].clone(), n[1].clone(), n[2].clone()), dec)?;
    if length.lo() <= r(0) {
        return Ok(None);
    }
    let Some(inverse) = length.recip() else {
        return Ok(None);
    };
    let orientation = Real::from_i64(orientation);
    Ok(Some([
        &real_ball(&(&n[0] * &orientation), dec)? * &inverse,
        &real_ball(&(&n[1] * &orientation), dec)? * &inverse,
        &real_ball(&(&n[2] * &orientation), dec)? * &inverse,
    ]))
}
fn interval_ball(x: Interval) -> Option<Ball> {
    Some(ball(&Rational::from_f64_exact(x.lo())?).hull(&ball(&Rational::from_f64_exact(x.hi())?)))
}
fn coverage(recipe: &Recipe) -> Vec<NotComparable> {
    recipe.steps.iter().filter_map(|step| {
        let reason=match &step.op {
            Operation::Block{..}|Operation::Cylinder{..}|Operation::Sphere{..}|
            Operation::Cone{..}|Operation::Torus{..}|Operation::Boolean{..}|
            Operation::Sketch{..}|Operation::Extrude{..}|Operation::ExtrudeThrough{..}|
            Operation::ExtrudeTo{..}|Operation::ExtrudeNext{..}|
            Operation::Revolve{..}|Operation::Hole{..}|Operation::Pattern{..}|Operation::Mirror{..}|
            Operation::Fillet{..}|Operation::Draft{..}|Operation::Chamfer{..}|Operation::Shell{..}|Operation::Offset{..}|Operation::Thicken{..}=>return None,
            Operation::Path{..}=>"neither kernel evaluates recipe paths",
            Operation::Sweep{..}=>"neither kernel evaluates recipe sweeps",
            Operation::Loft{..}=>"neither kernel evaluates recipe lofts",
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
struct FieldWitness {
    point: Point3,
    remaining: Rational,
    controls: catallax_shapes::distance::Options,
    upper: Option<Rational>,
}
impl FieldWitness {
    fn seed(
        prepared: &catallax_field::distance::Prepared<'_>,
        preferred: &[catallax_field::FaceId],
        point: &RPoint3,
        error: &Rational,
        options: &Options,
        dec: &mut Decider,
    ) -> Self {
        let point = point.to_real();
        let remaining = &options.accuracy - error;
        let mut controls = catallax_shapes::distance::Options {
            width: if remaining > r(0) {
                remaining.mul_pow2(-1)
            } else {
                options.accuracy.mul_pow2(-2)
            },
            ..Default::default()
        };
        let cell_limit = controls.max_cells;
        controls.max_cells = 1;
        let mut upper = None;
        for &id in preferred {
            if let Ok(distance) = prepared.point_to_face(id, &point, &controls, dec)
                && let Some(bound) = distance.upper
                && upper.as_ref().is_none_or(|old| bound < *old)
            {
                upper = Some(bound);
            }
            if upper.as_ref().is_some_and(|bound| bound <= &remaining) {
                break;
            }
        }
        controls.max_cells = cell_limit;
        Self {
            point,
            remaining,
            controls,
            upper,
        }
    }
    fn near_upper(&self) -> Option<&Rational> {
        self.upper
            .as_ref()
            .filter(|bound| *bound <= &self.remaining)
    }
    fn refine(
        &mut self,
        prepared: &catallax_field::distance::Prepared<'_>,
        proven_gap: bool,
        dec: &mut Decider,
    ) {
        if (self.upper.is_none()
            || (!proven_gap
                && self
                    .upper
                    .as_ref()
                    .is_some_and(|bound| bound > &self.remaining)))
            && let Ok(distance) = prepared.point_to_boundary(&self.point, &self.controls, dec)
            && let Some(bound) = distance.upper
            && self.upper.as_ref().is_none_or(|old| bound < *old)
        {
            self.upper = Some(bound);
        }
    }
}
fn retain_distance_upper(report: &mut Report, upper: Option<&Rational>) {
    if let Some(upper) = upper {
        if let Some(old) = &mut report.sampled_distance_upper
            && upper > &*old
        {
            *old = upper.clone();
        }
    } else {
        report.sampled_distance_upper = None;
    }
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
    retain_distance_upper(report, hi.as_ref());
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
// A rounded display corner need not lie on any finite sweep support patch.
// Its exact topological endpoint is only a candidate: the independent field
// must prove it is on the result boundary before it supplies a distance upper bound.
fn endpoint_distance(
    model: &catallax_topology::Model,
    field: &Plan,
    vertices: &BTreeSet<catallax_topology::VertexId>,
    point: &RPoint3,
    storage: &mut Workspace,
    dec: &mut Decider,
) -> Result<Option<Rational>, Refusal> {
    let mut upper: Option<Rational> = None;
    for vertex in vertices {
        let Some(candidate) = model
            .geom()
            .position(model.vertex(*vertex).point)
            .as_rational()
        else {
            continue;
        };
        if field.proven_checked(&candidate, storage, dec)?.location == Location::Boundary {
            let distance = norm(&point.minus(&candidate)).hi();
            upper = Some(match upper {
                Some(old) => old.min(distance),
                None => distance,
            });
        }
    }
    Ok(upper)
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
        sampled_distance_upper: Some(r(0)),
        findings: vec![],
        not_comparable: coverage(recipe),
        refusal: None,
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
    let mut dec = Decider::default();
    let field = match Plan::compile_checked(recipe, &mut dec) {
        Ok(f) => f,
        Err(e) => {
            whole_failure(&mut report, format!("field build {}: {}", e.code, e.detail));
            return revision.publish(report);
        }
    };
    let supports = match field
        .faces()
        .iter()
        .map(|face| support(&face.surface, &mut dec))
        .collect::<Result<Vec<_>, _>>()
    {
        Ok(supports) => supports,
        Err(e) => {
            whole_failure(
                &mut report,
                format!("support normalisation {}: {}", e.code, e.detail),
            );
            return revision.publish(report);
        }
    };
    let mut registry = BTreeMap::<Origin, Vec<catallax_field::FaceId>>::new();
    for (index, face) in field.faces().iter().enumerate() {
        if let Some(id) = field.result_face(index) {
            registry.entry(face.origin.clone()).or_default().push(id);
        }
    }
    let mut groups = Vec::with_capacity(supports.len());
    for source in &supports {
        let mut group = Vec::new();
        for (i, other) in supports.iter().enumerate() {
            match source.equals(other, &mut dec) {
                Ok(true) => group.push(i),
                Ok(false) => (),
                Err(e) => {
                    whole_failure(
                        &mut report,
                        format!("support comparison {}: {}", e.code, e.detail),
                    );
                    return revision.publish(report);
                }
            }
        }
        groups.push(group);
    }
    let (mut model, result) = if let Some(b) = built {
        (b.model.clone(), b.result)
    } else {
        let built = match evaluate(recipe, &mut dec) {
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
    let mesh = match catallax_mesh::mesh_body(&model, result, &options.mesh, &mut dec) {
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
    let Some(polynomials) = surfaces
        .iter()
        .map(Surface::implicit)
        .collect::<Option<Vec<_>>>()
    else {
        whole_failure(
            &mut report,
            "surface comparison needs a finite parametric correspondence solver",
        );
        return revision.publish(report);
    };
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
    let mut endpoints = BTreeMap::<usize, BTreeSet<catallax_topology::VertexId>>::new();
    for boundary in &mesh.edges {
        let edge = model.edge(boundary.edge);
        for (index, vertex) in boundary
            .vertices
            .first()
            .copied()
            .zip(edge.start)
            .into_iter()
            .chain(boundary.vertices.last().copied().zip(edge.end))
        {
            endpoints.entry(index).or_default().insert(vertex);
        }
    }
    let mut selected = BTreeSet::new();
    let mut sample_points = vec![None; mesh.faces.len()];
    let mut storage = Workspace::default();
    let width = options.accuracy.mul_pow2(-20);
    let mut distance_domains = None;
    let mut preferred = Vec::new();
    let mut prefer_witness = false;
    let choose_preferred = |i: usize, preferred: &mut Vec<catallax_field::FaceId>| {
        preferred.clear();
        for &index in &incident[i] {
            for origin in &face_origins[index] {
                if let Some(candidates) = registry.get(origin) {
                    preferred.extend(candidates.iter().copied());
                }
            }
        }
        // Names order candidates, never certify membership or exclude faces.
        // Incident planes have cheap linear projections; curved candidates
        // remain available when a cap projection misses its trim.
        preferred.sort_unstable_by_key(|id| {
            (
                !matches!(&field.faces()[id.index()].surface, Surface::Plane(_)),
                id.index(),
            )
        });
        preferred.dedup();
    };
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
        let mut witness = None;
        if prefer_witness {
            let prepared =
                distance_domains.get_or_insert_with(|| field.prepare_distances(&mut dec));
            if let Ok(prepared) = prepared {
                choose_preferred(i, &mut preferred);
                if !preferred.is_empty() {
                    let seed =
                        FieldWitness::seed(prepared, &preferred, &p, &error, options, &mut dec);
                    if let Some(bound) = seed.near_upper() {
                        // An admitted boundary point already proves proximity.
                        // Signed classification at the mesh point adds no
                        // distance evidence; reverse names/normals still run.
                        let upper = bound + &error;
                        retain_distance_upper(&mut report, Some(&upper));
                        continue;
                    }
                    witness = Some(seed);
                }
            }
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
        let mut upper = Rational::from_f64_exact(d.hi()).map(|h| &h + &error);
        if query.location != Location::Indeterminate
            && upper.as_ref().is_none_or(|h| h > &options.accuracy)
        {
            // A residual has no metric upper bound. Seek an independent,
            // trim-checked field boundary point instead of treating it as distance.
            'witness: for &face in &incident[i] {
                for origin in &face_origins[face] {
                    let Some(candidates) = registry.get(origin) else {
                        continue;
                    };
                    for &id in candidates {
                        let source = id.index();
                        let Ok(q) = self::witness(
                            &field.faces()[source].surface,
                            &supports[source],
                            &p,
                            &width,
                            &mut dec,
                        ) else {
                            continue;
                        };
                        let candidate = norm(&q.minus(&p)).hi() + error.clone();
                        if candidate <= options.accuracy
                            && field
                                .face_query_checked(id, &q, &mut storage, &mut dec)
                                .is_ok_and(|q| q.state == FaceState::OnFace)
                        {
                            upper = Some(candidate);
                            break 'witness;
                        }
                    }
                }
            }
        }
        let proven_gap =
            Rational::from_f64_exact(d.lo()).is_some_and(|lo| &lo - &error > options.accuracy);
        if query.location != Location::Indeterminate
            && (upper.is_none()
                || (!proven_gap && upper.as_ref().is_some_and(|h| h > &options.accuracy)))
        {
            let prepared =
                distance_domains.get_or_insert_with(|| field.prepare_distances(&mut dec));
            if let Ok(prepared) = prepared {
                let witness = witness.get_or_insert_with(|| {
                    choose_preferred(i, &mut preferred);
                    FieldWitness::seed(prepared, &preferred, &p, &error, options, &mut dec)
                });
                witness.refine(prepared, proven_gap, &mut dec);
                prefer_witness |= witness.near_upper().is_some();
                if let Some(bound) = &witness.upper {
                    let bound = bound + &error;
                    upper = Some(match upper {
                        Some(old) => old.min(bound),
                        None => bound,
                    });
                }
            }
        }
        if query.location != Location::Indeterminate
            && upper.as_ref().is_some_and(|h| h <= &options.accuracy)
        {
            retain_distance_upper(&mut report, upper.as_ref());
            continue;
        }
        if let Some(vertices) = endpoints.get(&i) {
            match endpoint_distance(&model, &field, vertices, &p, &mut storage, &mut dec) {
                Ok(Some(distance)) if &distance + &error <= options.accuracy => {
                    let upper = &distance + &error;
                    retain_distance_upper(&mut report, Some(&upper));
                    continue;
                }
                Err(error) => whole_failure(
                    &mut report,
                    format!("field endpoint witness {}: {}", error.code, error.detail),
                ),
                _ => (),
            }
        }
        retain_distance_upper(&mut report, upper.as_ref());
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
    let mut observed: BTreeMap<usize, BTreeSet<Origin>> = BTreeMap::new();
    let mut oriented = BTreeSet::new();
    let mut brep_distance_domains: Vec<
        Option<Result<catallax_measure::distance::PreparedFace, Refusal>>,
    > = (0..mesh.faces.len()).map(|_| None).collect();
    let sampling_failed = std::cell::Cell::new(false);
    let mut sample_refusal = None;
    let mut check_sample = |p: RPoint3, id: catallax_field::FaceId| {
        if sampling_failed.get() {
            return;
        }
        if options
            .region
            .as_ref()
            .is_some_and(|b| !b.contains(&p, &r(0)))
        {
            return;
        }
        let source = id.index();
        let origin = &field.faces()[source].origin;
        let query = match field.closed_face_query_checked(id, &p, &mut storage, &mut dec) {
            Ok(query) => query,
            Err(error) => {
                sampling_failed.set(true);
                unresolved(
                    &mut report,
                    Pass::Reverse,
                    &p,
                    vec![],
                    vec![origin.clone()],
                    format!("field face query {}: {}", error.code, error.detail),
                );
                sample_refusal = Some(error);
                return;
            }
        };
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
            let Some(candidate) = field.result_face(index) else {
                continue;
            };
            let face = &field.faces()[index];
            let coactive = if index == source {
                query
            } else {
                match field.closed_face_query_checked(candidate, &p, &mut storage, &mut dec) {
                    Ok(query) => query,
                    Err(error) => {
                        sampling_failed.set(true);
                        unresolved(
                            &mut report,
                            Pass::Reverse,
                            &p,
                            vec![],
                            vec![face.origin.clone()],
                            format!("coactive field face query {}: {}", error.code, error.detail),
                        );
                        sample_refusal = Some(error);
                        return;
                    }
                }
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
        let mut closest: Option<(Rational, usize, Point3)> = None;
        let mut projection_errors = Vec::new();
        for &i in &candidates {
            if revision.is_cancelled() {
                return;
            }
            let face = mesh.faces[i].face;
            let surface = &surfaces[i];
            let d = match support_distance(surface, &p, &mut dec) {
                Ok(d) => d.lo().max(r(0)),
                Err(e) => {
                    projection_errors.push(format!("{}: {}", e.code, e.detail));
                    r(0) // The universal lower bound remains valid.
                }
            };
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
                closest = Some((upper, i, q.to_real()));
            }
        }
        if closest.is_none()
            || (lower.as_ref().is_none_or(|lo| lo <= &options.accuracy)
                && closest
                    .as_ref()
                    .is_some_and(|(upper, _, _)| upper > &options.accuracy))
        {
            let controls = catallax_shapes::distance::Options {
                width: options.accuracy.mul_pow2(-2),
                ..Default::default()
            };
            let mut trimmed_lower: Option<Rational> = None;
            for &i in &candidates {
                let prepared = brep_distance_domains[i].get_or_insert_with(|| {
                    catallax_measure::distance::PreparedFace::new(
                        &model,
                        mesh.faces[i].face,
                        &mut dec,
                    )
                });
                let distance = match prepared {
                    Ok(prepared) => {
                        prepared.point_to_face(&mut model, &p.to_real(), &controls, &mut dec)
                    }
                    Err(e) => Err(e.clone()),
                };
                match distance {
                    Ok(distance) => {
                        trimmed_lower = Some(
                            trimmed_lower
                                .map_or(distance.lower.clone(), |old| old.min(distance.lower)),
                        );
                        if let (Some(upper), Some([_, q])) = (distance.upper, distance.witnesses)
                            && closest.as_ref().is_none_or(|(best, _, _)| upper < *best)
                        {
                            closest = Some((upper, i, q));
                        }
                    }
                    Err(e) => {
                        trimmed_lower = Some(r(0));
                        projection_errors.push(format!("{}: {}", e.code, e.detail));
                    }
                }
            }
            if let Some(trimmed) = trimmed_lower {
                lower = Some(lower.map_or(trimmed.clone(), |old| old.max(trimmed)));
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
            retain_distance_upper(&mut report, None);
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
        let brep_normal = match direction(&gradients[i], &q, brep.sense, &mut dec) {
            Ok(n) => n,
            Err(e) => {
                unresolved(
                    &mut report,
                    Pass::Reverse,
                    &p,
                    vec![expected.to_vec()],
                    vec![origin.clone()],
                    format!("{}: {}", e.code, e.detail),
                );
                return;
            }
        };
        if let (Direction::Unique(field_normal), Some(brep_normal)) = (query.direction, brep_normal)
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
    let sampling = field.visit_candidates(options.divisions, |p, id| {
        if revision.is_cancelled() {
            return false;
        }
        check_sample(p, id);
        !revision.is_cancelled() && !sampling_failed.get()
    });
    let distance_sampling = if recipe.steps.iter().any(|step| matches!(step.op,
        Operation::Offset { .. } | Operation::Shell { .. } | Operation::Thicken { .. })) {
        field.visit_distance_candidates(options.divisions, |p, origin| {
            if revision.is_cancelled() { return false; }
            check_sample(p, origin);
            !revision.is_cancelled() && !sampling_failed.get()
        })
    } else {
        Ok(true)
    };
    if let Some(refusal) = sample_refusal {
        whole_failure(
            &mut report,
            format!("field face query {}: {}", refusal.code, refusal.detail),
        );
        report.refusal = Some(refusal);
        return revision.publish(report);
    }
    if let Err(error) = sampling {
        whole_failure(
            &mut report,
            format!("field boundary sampling {}: {}", error.code, error.detail),
        );
    }
    if let Err(error) = distance_sampling {
        whole_failure(
            &mut report,
            format!("distance-feature boundary sampling {}: {}", error.code, error.detail),
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
