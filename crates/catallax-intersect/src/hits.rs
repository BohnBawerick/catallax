use crate::{
    Contact, CurveForm, CurveHit, CurveSurface,
    chart::{self, Chart},
    math::*,
    pair::Plan,
    query,
    solve::{self, Solutions},
};
use catallax_numbers::{
    Answer, Answered, BoundKind, Decider, FactKind, Mark, MarkLog, Outcome, Rational, Real,
    RealPoly, Refusal, Sign,
};
use catallax_shapes::{
    Circle, Curve, CurveId, Ellipse, GeomKey, GeomStore, Line, Meet, MeetOf, Point3, PointDef,
    Position, SurfaceId,
};
use std::sync::Arc;

#[derive(Clone)]
struct Candidate {
    position: Point3,
    contact: Contact,
    uncertain_branch: bool,
}
fn distance(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<Rational, Refusal> {
    let d = a.minus(b);
    let mut sum = r(0);
    for x in [&d.x, &d.y, &d.z] {
        let b = dec.enclose(x, &r(1))?;
        let m = b.lo().abs().max(b.hi().abs());
        sum = sum + m;
    }
    Ok(sum)
}
fn mark(
    marks: &mut MarkLog,
    n: catallax_numbers::NotExact,
    kind: BoundKind,
    bound: Rational,
    detail: &str,
) {
    marks.push(Mark::from_not_exact(&n, kind, bound, detail));
}
fn carrier_equations(
    geom: &GeomStore,
    c: CurveId,
    f: &CurveForm,
    dec: &mut Decider,
) -> Result<Vec<RealPoly>, Refusal> {
    let curve = match f {
        CurveForm::Line { origin, dir } => {
            let mut rational_origin = origin.as_rational();
            if rational_origin.is_none()
                && let Some(plan) = query::validate(geom, c, dec)?
            {
                for spec in &plan.curves {
                    if let CurveForm::Line {
                        origin,
                        dir: canonical,
                    } = &spec.form
                        && sign(&canonical.dot(dir), dec)? == Sign::Positive
                    {
                        rational_origin = origin.as_rational();
                        if rational_origin.is_some() {
                            break;
                        }
                    }
                }
            }
            rational_origin.map(|origin| {
                Curve::Line(Line {
                    origin: origin.to_real(),
                    dir: dir.clone(),
                })
            })
        }
        CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ref_dir,
        } if center.as_rational().is_some() && radius_sq.rational_ref().is_some() => {
            Some(Curve::Circle(Circle {
                center: center.clone(),
                normal: normal.clone(),
                radius_sq: radius_sq.clone(),
                ref_dir: ref_dir.clone(),
            }))
        }
        CurveForm::Ellipse {
            center,
            normal,
            major_dir,
            major_sq,
            minor_sq,
        } if center.as_rational().is_some()
            && major_sq.rational_ref().is_some()
            && minor_sq.rational_ref().is_some() =>
        {
            Some(Curve::Ellipse(Ellipse {
                center: center.clone(),
                normal: normal.clone(),
                major_dir: major_dir.clone(),
                major_sq: major_sq.clone(),
                minor_sq: minor_sq.clone(),
            }))
        }
        _ => None,
    };
    Ok(curve.as_ref().unwrap_or(geom.curve(c)).equations(geom))
}
fn keep_branch(
    geom: &GeomStore,
    c: CurveId,
    s: SurfaceId,
    p: &Point3,
    plan: &Option<Arc<Plan>>,
    dec: &mut Decider,
    marks: &mut MarkLog,
) -> Result<Option<bool>, Refusal> {
    let supports: &[SurfaceId] = match geom.curve(c) {
        Curve::Intersection(i) => &i.surfaces,
        _ => &[],
    };
    for surface in std::iter::once(s).chain(supports.iter().copied()) {
        match query::surface_domain(geom.surface(surface), p, dec) {
            Ok(true) => {}
            Ok(false) => return Ok(None),
            Err(e) => return Err(Refusal::from(e)),
        }
    }
    let Some(d) = plan.as_ref().and_then(|p| p.divider.as_ref()) else {
        return Ok(Some(false));
    };
    let Curve::Intersection(i) = geom.curve(c) else {
        unreachable!()
    };
    let witness = d.value(&geom.position(i.witness));
    let value = d.value(p);
    let ps = match dec.sign(&value) {
        Outcome::Proven(s) => Some(s),
        Outcome::ProvenCoincident(_) => return Ok(Some(false)),
        Outcome::Refused(e) => return Err(e),
        // Only legacy divider plans permit bounded, marked component choices.
        Outcome::NotExact(n) if !d.canonical_direction => {
            return Err(Refusal::from(Failure::Uncertain(n)));
        }
        Outcome::NotExact(n) => {
            let bound = &n.bound / &d.length_lower;
            mark(
                marks,
                n,
                BoundKind::Distance,
                bound,
                "candidate kept on unresolved branch",
            );
            return Ok(Some(true));
        }
    };
    match dec.sign(&witness) {
        Outcome::Proven(s) => Ok((Some(s) == ps).then_some(false)),
        Outcome::ProvenCoincident(_) => panic!("regular witness on divider"),
        Outcome::Refused(e) => Err(e),
        Outcome::NotExact(n) if !d.canonical_direction => Err(Refusal::from(Failure::Uncertain(n))),
        Outcome::NotExact(n) => {
            let b = dec.enclose(&value, &r(1))?;
            let bound = &b.lo().abs().max(b.hi().abs()) / &d.length_lower;
            mark(
                marks,
                n,
                BoundKind::Distance,
                bound,
                "candidate kept on unresolved witness branch",
            );
            Ok(Some(true))
        }
    }
}
fn convolution(a: &[Real], b: &[Real], n: usize) -> Vec<Real> {
    (0..=n)
        .map(|k| (0..=k).fold(Real::zero(), |v, i| v + &a[i] * &b[k - i]))
        .collect()
}
fn polynomial_jet(poly: &RealPoly, values: &[[Real; 3]], n: usize) -> Vec<Real> {
    let mut out = vec![Real::zero(); n + 1];
    for (m, c) in poly.terms() {
        let mut term = vec![Real::zero(); n + 1];
        term[0] = real(&c);
        for (axis, exp) in m.exponents().iter().enumerate() {
            let v: Vec<_> = values.iter().map(|v| v[axis].clone()).collect();
            for _ in 0..*exp {
                term = convolution(&term, &v, n);
            }
        }
        for k in 0..=n {
            out[k] = &out[k] + &term[k];
        }
    }
    out
}
fn local_sign(
    ch: &Chart,
    t: &Real,
    incoming: bool,
    poly: &RealPoly,
    dec: &mut Decider,
) -> Result<Sign, Failure> {
    // Higher Taylor coefficients cannot affect the first nonzero one. Most
    // contacts need only order one; construct later coefficients only if zero
    // has actually been proved, rather than expanding all eight in advance.
    for k in 1..=8 {
        let j = polynomial_jet(poly, &ch.exact_jets(t, k), k);
        let s = sign(&j[k], dec)?;
        if s != Sign::Zero {
            return Ok(if incoming && k % 2 == 1 {
                if s == Sign::Positive {
                    Sign::Negative
                } else {
                    Sign::Positive
                }
            } else {
                s
            });
        }
    }
    panic!("a non-contained quartic on a conic has contact order at most eight")
}
fn contact(
    charts: &[Chart],
    p: &Point3,
    poly: &RealPoly,
    dec: &mut Decider,
) -> Result<Contact, Failure> {
    let mut incoming = None;
    let mut leaving = None;
    for ch in charts {
        if let Some(t) = ch.locate(p, dec)? {
            if sign(&t, dec)? != Sign::Zero && incoming.is_none() {
                incoming = Some(local_sign(ch, &t, true, poly, dec)?);
            }
            if sign(&(&Real::one() - &t), dec)? != Sign::Zero && leaving.is_none() {
                leaving = Some(local_sign(ch, &t, false, poly, dec)?);
            }
            if incoming.is_some() && leaving.is_some() {
                break;
            }
        }
    }
    Ok(
        if incoming.expect("arriving chart") == leaving.expect("leaving chart") {
            Contact::Tangent
        } else {
            Contact::Crossing
        },
    )
}
fn compare_locations(
    a: &(usize, Real),
    b: &(usize, Real),
    dec: &mut Decider,
) -> Result<std::cmp::Ordering, Failure> {
    if a.0 != b.0 {
        return Ok(a.0.cmp(&b.0));
    }
    Ok(match sign(&(&a.1 - &b.1), dec)? {
        Sign::Negative => std::cmp::Ordering::Less,
        Sign::Positive => std::cmp::Ordering::Greater,
        Sign::Zero => std::cmp::Ordering::Equal,
    })
}
fn compare(
    a: &Candidate,
    b: &Candidate,
    form: &CurveForm,
    charts: &[Chart],
    origin: &Option<(usize, Real)>,
    dec: &mut Decider,
) -> Result<std::cmp::Ordering, Failure> {
    if let CurveForm::Line { dir, .. } = form {
        return Ok(match sign(&a.position.minus(&b.position).dot(dir), dec)? {
            Sign::Negative => std::cmp::Ordering::Less,
            Sign::Positive => std::cmp::Ordering::Greater,
            Sign::Zero => std::cmp::Ordering::Equal,
        });
    }
    if let CurveForm::Parabola {
        normal, axis_dir, ..
    }
    | CurveForm::Hyperbola {
        normal, axis_dir, ..
    } = form
    {
        return Ok(
            match sign(
                &a.position.minus(&b.position).dot(&normal.cross(axis_dir)),
                dec,
            )? {
                Sign::Negative => std::cmp::Ordering::Less,
                Sign::Positive => std::cmp::Ordering::Greater,
                Sign::Zero => std::cmp::Ordering::Equal,
            },
        );
    }
    let la = chart::locate(charts, &a.position, dec)?;
    let lb = chart::locate(charts, &b.position, dec)?;
    let o = origin.as_ref().unwrap();
    let wrap_a = compare_locations(&la, o, dec)? == std::cmp::Ordering::Less;
    let wrap_b = compare_locations(&lb, o, dec)? == std::cmp::Ordering::Less;
    if wrap_a != wrap_b {
        return Ok(wrap_a.cmp(&wrap_b));
    }
    compare_locations(&la, &lb, dec)
}
pub(super) fn intersect(
    geom: &mut GeomStore,
    c: CurveId,
    s: SurfaceId,
    dec: &mut Decider,
) -> Answered<CurveSurface> {
    let plan = query::validate(geom, c, dec)?;
    if geom
        .facts()
        .holds(FactKind::Incident, GeomKey::Curve(c), GeomKey::Surface(s))
    {
        dec.note_fact_decision();
        return Ok(Answer::exact(CurveSurface::Contained));
    }
    let form = query::form(geom, c, dec)?;
    let equations = carrier_equations(geom, c, &form, dec)?;
    let implicit = geom.implicit(s);
    let solutions = solve::intersections(&equations, &implicit, dec)?;
    let Solutions::Points(points) = solutions else {
        return Ok(Answer::exact(CurveSurface::Contained));
    };
    let mut marks = MarkLog::new();
    let mut candidates = vec![];
    for position in points {
        if let Some(uncertain_branch) = keep_branch(geom, c, s, &position, &plan, dec, &mut marks)?
        {
            candidates.push(Candidate {
                position,
                contact: Contact::Crossing,
                uncertain_branch,
            });
        }
    }
    if candidates.is_empty() {
        return Ok(Answer::with_marks(CurveSurface::Disjoint, marks));
    }
    let open_conic = matches!(
        form,
        CurveForm::Parabola { .. } | CurveForm::Hyperbola { .. }
    );
    let charts = if matches!(form, CurveForm::Line { .. }) || open_conic {
        vec![]
    } else {
        query::charts(geom, c, dec).map_err(Refusal::from)?
    };
    for candidate in &mut candidates {
        let answer = if let CurveForm::Line { dir, .. } = &form {
            match sign(&geom.surface(s).gradient(&candidate.position).dot(dir), dec) {
                Ok(Sign::Zero) if matches!(geom.surface(s), catallax_shapes::Surface::Torus(_)) => {
                    match chart::line_local(&candidate.position, dir, dec) {
                        Ok(local) => contact(
                            std::slice::from_ref(&local),
                            &candidate.position,
                            &implicit,
                            dec,
                        ),
                        Err(e) => Err(Failure::Refused(e)),
                    }
                }
                Ok(Sign::Zero) => Ok(Contact::Tangent),
                Ok(_) => Ok(Contact::Crossing),
                Err(e) => Err(e),
            }
        } else if open_conic {
            match chart::conic_local(&form, &candidate.position, dec) {
                Ok(local) => contact(
                    std::slice::from_ref(&local),
                    &candidate.position,
                    &implicit,
                    dec,
                ),
                Err(e) => Err(Failure::Refused(e)),
            }
        } else {
            contact(&charts, &candidate.position, &implicit, dec)
        };
        candidate.contact = match answer {
            Ok(c) => c,
            Err(Failure::Refused(e)) => return Err(e),
            Err(Failure::Uncertain(n)) => {
                let bound = n.bound.clone();
                mark(
                    &mut marks,
                    n,
                    BoundKind::Value,
                    bound,
                    "unresolved hit contact",
                );
                Contact::Tangent
            }
        };
    }
    let origin = if charts.is_empty() {
        None
    } else {
        Some(
            chart::locate(&charts, &geom.curve(c).parameter_origin(geom), dec)
                .map_err(Refusal::from)?,
        )
    };
    let mut sorted: Vec<Candidate> = vec![];
    for mut candidate in candidates {
        let mut i = 0;
        loop {
            if i == sorted.len() {
                sorted.push(candidate);
                break;
            }
            match compare(&candidate, &sorted[i], &form, &charts, &origin, dec) {
                Ok(std::cmp::Ordering::Less) => {
                    sorted.insert(i, candidate);
                    break;
                }
                Ok(std::cmp::Ordering::Greater) => i += 1,
                Ok(std::cmp::Ordering::Equal) => {
                    sorted[i].contact = Contact::Tangent;
                    sorted[i].uncertain_branch |= candidate.uncertain_branch;
                    break;
                }
                Err(Failure::Refused(e)) => return Err(e),
                Err(Failure::Uncertain(n)) => {
                    let bound = distance(&candidate.position, &sorted[i].position, dec)?;
                    mark(
                        &mut marks,
                        n,
                        BoundKind::Distance,
                        bound,
                        "merged hits with unresolved order",
                    );
                    candidate.contact = Contact::Tangent;
                    candidate.uncertain_branch |= sorted[i].uncertain_branch;
                    sorted.remove(i);
                }
            }
        }
    }
    let hits = sorted
        .into_iter()
        .map(|candidate| {
            let of = if candidate.uncertain_branch
                && let Curve::Intersection(i) = geom.curve(c)
            {
                let mut surfaces = i.surfaces.to_vec();
                if !surfaces.contains(&s) {
                    surfaces.push(s);
                }
                MeetOf::Surfaces(surfaces)
            } else {
                MeetOf::CurveSurface(c, s)
            };
            let point = geom.add_point(PointDef::Meet(Meet {
                of,
                position: Position::Closed(candidate.position),
            }));
            CurveHit {
                point,
                contact: candidate.contact,
            }
        })
        .collect();
    Ok(Answer::with_marks(CurveSurface::Meets(hits), marks))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::BudgetPolicy;
    use catallax_shapes::{Cone, RPoint3, RVector3, Surface, Torus};

    fn point(x: i64, y: i64, z: i64) -> RPoint3 {
        RPoint3::new(r(x), r(y), r(z))
    }
    fn vector(x: i64, y: i64, z: i64) -> RVector3 {
        RVector3::new(r(x), r(y), r(z))
    }
    fn assert_never_kept(geom: &GeomStore, curve: CurveId, surface: SurfaceId, position: Point3) {
        let plan = query::validate(geom, curve, &mut Decider::default()).unwrap();
        for (mut decider, limited) in [
            (Decider::default(), false),
            (
                Decider::new(BudgetPolicy {
                    max_precision_bits: 64,
                    max_work: 0,
                }),
                true,
            ),
        ] {
            let answer = keep_branch(
                geom,
                curve,
                surface,
                &position,
                &plan,
                &mut decider,
                &mut MarkLog::new(),
            );
            assert!(
                matches!(&answer, Ok(None)) || (limited && answer.is_err()),
                "a budget cannot select a candidate on the wrong component: {answer:?}"
            );
        }
    }

    #[test]
    fn uncertain_nappe_cannot_keep_a_negative_cone_hit() {
        let mut geom = GeomStore::new();
        let curve = geom.add_curve(Curve::Line(
            Line::new(
                point(1, 0, 0).into(),
                vector(0, 0, 1).into(),
                &mut Decider::default(),
            )
            .unwrap(),
        ));
        let surface = geom.add_surface(Surface::Cone(
            Cone::new(
                point(0, 0, 0).into(),
                vector(0, 0, 1).into(),
                r(1).into(),
                vector(1, 0, 0).into(),
                &mut Decider::default(),
            )
            .unwrap(),
        ));
        assert_never_kept(&geom, curve, surface, point(1, 0, -1).to_real());
    }

    #[test]
    fn uncertain_radius_cannot_keep_a_hit_on_the_other_torus_circle() {
        let mut geom = GeomStore::new();
        let plane = geom.add_surface(Surface::Plane(
            catallax_shapes::Plane::new(
                RPoint3::new(r(0), r(0), &r(1) / &r(2)).into(),
                vector(0, 0, 1).into(),
                vector(1, 0, 0).into(),
                &mut Decider::default(),
            )
            .unwrap(),
        ));
        let torus = geom.add_surface(Surface::Torus(
            Torus::new(
                point(0, 0, 0).into(),
                vector(0, 0, 1).into(),
                r(3).into(),
                r(1).into(),
                vector(1, 0, 0).into(),
                &mut Decider::default(),
            )
            .unwrap(),
        ));
        let answer =
            crate::intersect_surfaces(&mut geom, plane, torus, &mut Decider::default()).unwrap();
        let crate::SurfacePair::Meets(components) = answer.value else {
            panic!("expected two circular sections")
        };
        let [
            crate::Component::Curve { curve: first, .. },
            crate::Component::Curve { curve: other, .. },
        ] = components.as_slice()
        else {
            panic!("expected two circular sections")
        };
        let position = geom.curve(*other).parameter_origin(&geom);
        assert_never_kept(&geom, *first, plane, position);
    }
}
