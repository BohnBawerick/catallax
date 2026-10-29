use crate::{
    Contact, CurveForm, CurveHit, CurveSurface,
    chart::{self, Chart},
    math::*,
    pair::{CylinderTorusMode, Frame, Plan},
    query,
    solve::{self, Solutions},
};
use catallax_numbers::{
    Answer, Answered, BoundKind, Decider, Fact, FactKind, FactOrigin, Mark, MarkLog, Outcome,
    Rational, Real, RealPoly, Refusal, Sign, StrictSign,
};
use catallax_shapes::{
    Circle, Curve, CurveId, Cylinder, Ellipse, GeomKey, GeomStore, Line, Meet, MeetOf, Plane,
    Point3, PointDef, Position, Surface, SurfaceId, Vector3,
};
use std::sync::Arc;

#[derive(Clone)]
struct Candidate {
    position: Point3,
    contact: Contact,
    uncertain_branch: bool,
    proven_tangent: bool,
}
fn known_tangent(geom: &GeomStore, c: CurveId, s: SurfaceId, position: &Point3) -> bool {
    let a = GeomKey::Curve(c);
    let b = GeomKey::Surface(s);
    geom.facts().iter().any(|fact| {
        if fact.kind != FactKind::Tangent
            || !((fact.a == a && fact.b == b) || (fact.a == b && fact.b == a))
        {
            return false;
        }
        let Some(GeomKey::Point(at)) = fact.at else {
            return false;
        };
        let point = match geom.point(at) {
            PointDef::Explicit(point)
            | PointDef::Meet(Meet {
                position: Position::Closed(point),
                ..
            }) => point,
            _ => return false,
        };
        point.x.same_definition(&position.x)
            && point.y.same_definition(&position.y)
            && point.z.same_definition(&position.z)
    })
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
    curve
        .as_ref()
        .unwrap_or(geom.curve(c))
        .equations(geom)
        .ok_or_else(|| {
            Refusal::new(
                "intersect.parametric-curve",
                "this intersection needs a finite parametric curve solver",
            )
        })
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
    if let Some(keep) = query::interval_component_contains(geom, c, p, plan.as_deref(), dec)? {
        return Ok(keep.then_some(false));
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
// These residuals are constant around the complete circle. Keeping squared
// radii avoids adjoining an unnecessary square root to the coefficient field.
fn circular_section(
    geom: &GeomStore,
    c: CurveId,
    s: SurfaceId,
    form: &CurveForm,
    dec: &mut Decider,
) -> Result<Option<CurveSurface>, Refusal> {
    let CurveForm::Circle {
        center,
        normal,
        radius_sq,
        ..
    } = form
    else {
        return Ok(None);
    };
    let surface = geom.surface(s);
    if center.as_rational().is_some()
        && normal.as_rational().is_some()
        && radius_sq.rational_ref().is_some()
        && surface.has_rational_definition()
    {
        // Rational sections use the closed-form or general solver below.
        return Ok(None);
    }
    let residual = if let Surface::Plane(plane) = surface {
        if sign(&normal.cross(&plane.normal).norm_squared(), dec)? != Sign::Zero {
            return Ok(None);
        }
        center.minus(&plane.origin).dot(&plane.normal)
    } else {
        if !matches!(surface, Surface::Cylinder(_) | Surface::Cone(_))
            || !geom
                .facts()
                .holds(FactKind::Coaxial, GeomKey::Curve(c), GeomKey::Surface(s))
        {
            return Ok(None);
        }
        dec.note_fact_decision();
        match surface {
            Surface::Cylinder(cylinder) => radius_sq - &cylinder.radius_sq,
            Surface::Cone(cone) => {
                let height = center.minus(&cone.apex).dot(&cone.axis);
                if sign(&height, dec)? == Sign::Negative {
                    return Ok(Some(CurveSurface::Disjoint));
                }
                radius_sq
                    - &(&cone.slope.square() * &(&height.square() / &cone.axis.norm_squared()))
            }
            _ => unreachable!("circular carrier checked above"),
        }
    };
    Ok(Some(if sign(&residual, dec)? == Sign::Zero {
        CurveSurface::Contained
    } else {
        CurveSurface::Disjoint
    }))
}

// A line already supplies a separating parameter. Eliminating x, y and z
// independently loses that correlation and needlessly eliminates coefficients.
fn algebraic_line(
    origin: &Point3,
    direction: &catallax_shapes::Vector3,
    surface: &Surface,
    implicit: &RealPoly,
    dec: &mut Decider,
) -> Result<Solutions, Refusal> {
    // Substitute into a rotational carrier in its defining frame. Expanding a
    // translated cone in world coordinates first obscures its shared axial
    // expressions and then requires eliminating the same identities again.
    let rotational = match surface {
        Surface::Cone(cone) => Some((&cone.apex, &cone.axis, Some(&cone.slope), None)),
        Surface::Cylinder(cylinder) => Some((
            &cylinder.origin,
            &cylinder.axis,
            None,
            Some(&cylinder.radius_sq),
        )),
        _ => None,
    };
    let polynomial = if let Some((center, axis, slope, radius_sq)) = rotational {
        let offset = origin.minus(center);
        let perpendicular = offset.cross(axis);
        let velocity = direction.cross(axis);
        let mut a = velocity.norm_squared();
        let mut b = &Real::from_i64(2) * &perpendicular.dot(&velocity);
        let mut c = perpendicular.norm_squared();
        if let Some(slope) = slope {
            let height = offset.dot(axis);
            let rate = direction.dot(axis);
            let squared = slope.square();
            a = &a - &(&squared * &rate.square());
            b = &b - &(&Real::from_i64(2) * &(&squared * &(&height * &rate)));
            c = &c - &(&squared * &height.square());
        } else if let Some(radius_sq) = radius_sq {
            c = &c - &(radius_sq * &axis.norm_squared());
        }
        RealPoly::from_terms(1, vec![(vec![0], c), (vec![1], b), (vec![2], a)])
    } else {
        let coordinates: [RealPoly; 3] = std::array::from_fn(|axis| {
            let origin = [&origin.x, &origin.y, &origin.z];
            let direction = [&direction.x, &direction.y, &direction.z];
            RealPoly::from_terms(
                1,
                vec![
                    (vec![0], origin[axis].clone()),
                    (vec![1], direction[axis].clone()),
                ],
            )
        });
        let mut polynomial = RealPoly::zero(1);
        for (monomial, coefficient) in implicit.terms() {
            dec.check_deadline()?;
            let mut term = RealPoly::constant(1, coefficient);
            for (axis, power) in monomial.exponents().iter().enumerate() {
                if *power != 0 {
                    term = &term * &coordinates[axis].pow_u32(*power);
                }
            }
            polynomial = &polynomial + &term;
        }
        polynomial
    };
    if polynomial.is_zero(dec)? {
        return Ok(Solutions::Contained);
    }
    let roots = polynomial.real_roots(dec)?;
    let mut points = Vec::with_capacity(roots.len());
    for root in roots {
        dec.check_deadline()?;
        points.push(origin.plus(&direction.scaled(&root)));
    }
    Ok(Solutions::Points(points))
}

fn cylinder_excludes_plane(
    cylinder: &Cylinder,
    origin: &Point3,
    normal: &Vector3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    if !normal
        .dot(&cylinder.axis)
        .rational_ref()
        .is_some_and(Rational::is_zero)
    {
        return Ok(false);
    }
    let height = cylinder.origin.minus(origin).dot(normal);
    let separation = &height.square() - &(&cylinder.radius_sq * &normal.norm_squared());
    Ok(sign(&separation, dec)? == Sign::Positive)
}

fn cylinder_torus_excludes_plane(
    witness: &Point3,
    plane: &Plane,
    plan: &Plan,
    dec: &Decider,
) -> Result<bool, Refusal> {
    let Some(Frame::CylinderTorus {
        origin,
        e1,
        e2,
        e3,
        radius,
        radial_shift,
        radial_slope,
        lo,
        span,
        lower_radius,
        upper_radius,
        mode,
        ..
    }) = &plan.frame
    else {
        return Ok(false);
    };
    // Bound the complete component in its orthonormal frame, then project the
    // plane into that frame. Bounding world coordinates first loses correlations
    // under a rigid placement. Strict exclusion keeps endpoints and tangencies.
    let attempt = || -> Result<bool, Refusal> {
        let enclose = |value: &Real| dec.enclose_at(value, 128);
        let x = enclose(lo)?.hull(&enclose(&(lo + span))?);
        let Some(mut y) = (&enclose(radius)?.square() - &x.square()).sqrt() else {
            return Ok(false);
        };
        let Some(radial) = (&enclose(radial_shift)? + &(&enclose(radial_slope)? * &x)).sqrt()
        else {
            return Ok(false);
        };
        let Some(mut z) =
            (&(&enclose(upper_radius)? - &radial) * &(&radial - &enclose(lower_radius)?)).sqrt()
        else {
            return Ok(false);
        };
        match mode {
            CylinderTorusMode::Single => {
                y = y.hull(&(-&y));
                z = z.hull(&(-&z));
            }
            CylinderTorusMode::Axial | CylinderTorusMode::Transverse => {
                let (selected, both, axis) = if *mode == CylinderTorusMode::Axial {
                    (&mut z, &mut y, e3)
                } else {
                    (&mut y, &mut z, e2)
                };
                match enclose(&witness.minus(origin).dot(axis))?.sign() {
                    Some(StrictSign::Positive) => {}
                    Some(StrictSign::Negative) => *selected = -&*selected,
                    None => return Ok(false),
                }
                *both = both.hull(&(-&*both));
            }
        }
        let mut distance = enclose(&origin.minus(&plane.origin).dot(&plane.normal))?;
        for (axis, bound) in [e1, e2, e3].into_iter().zip([x, y, z]) {
            distance = &distance + &(&enclose(&axis.dot(&plane.normal))? * &bound);
        }
        Ok(!distance.contains_zero())
    };
    match attempt() {
        Err(refusal) if refusal.code != "numbers.deadline" => Ok(false),
        answer => answer,
    }
}

fn torus_excludes_plane(
    torus: &catallax_shapes::Torus,
    plane: &Plane,
    dec: &Decider,
) -> Result<bool, Refusal> {
    // The complete ring torus projects onto [-h,h], where
    // h = R |normal_radial| + r |normal|. Strict separation preserves contacts.
    let radial_sq = plane.normal.cross(&torus.axis).norm_squared() / torus.axis.norm_squared();
    let extent = &torus.major * &root_real(&radial_sq)
        + &torus.minor * &root_real(&plane.normal.norm_squared());
    let separation = torus
        .center
        .minus(&plane.origin)
        .dot(&plane.normal)
        .square()
        - extent.square();
    match dec.enclose_at(&separation, 128) {
        Ok(bound) => Ok(bound.sign() == Some(StrictSign::Positive)),
        Err(refusal) if refusal.code == "numbers.deadline" => Err(refusal),
        Err(_) => Ok(false),
    }
}

fn branch_plane_disjoint(
    geom: &GeomStore,
    curve: CurveId,
    surface: SurfaceId,
    plan: Option<&Plan>,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let (Curve::Intersection(branch), Surface::Plane(plane)) =
        (geom.curve(curve), geom.surface(surface))
    else {
        return Ok(false);
    };
    // The entire carrier lies on each support. This excludes no point at
    // tangency and asserts nothing about trim membership or orientation.
    for support in &branch.surfaces {
        let excluded = match geom.surface(*support) {
            Surface::Cylinder(cylinder) => {
                cylinder_excludes_plane(cylinder, &plane.origin, &plane.normal, dec)?
            }
            Surface::Torus(torus) => torus_excludes_plane(torus, plane, dec)?,
            _ => false,
        };
        if excluded {
            return Ok(true);
        }
    }
    if let Some(plan) = plan {
        return cylinder_torus_excludes_plane(&geom.position(branch.witness), plane, plan, dec);
    }
    Ok(false)
}

fn circular_cylinder(
    form: &CurveForm,
    surface: &Surface,
    dec: &mut Decider,
) -> Result<Option<Solutions>, Refusal> {
    let (
        CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ..
        },
        Surface::Cylinder(cylinder),
    ) = (form, surface)
    else {
        return Ok(None);
    };
    // Parallel normals use the circular section below. A circle plane parallel
    // to the cylinder axis can first be excluded by its distance from that axis.
    if !normal
        .cross(&cylinder.axis)
        .norm_squared()
        .rational_ref()
        .is_some_and(Rational::is_zero)
    {
        if cylinder_excludes_plane(cylinder, center, normal, dec)? {
            return Ok(Some(Solutions::Points(Vec::new())));
        }
        return Ok(None);
    }
    let offset = cylinder.origin.minus(center);
    let radial = offset.minus(&normal.scaled(&(offset.dot(normal) / normal.norm_squared())));
    let distance_sq = radial.norm_squared();
    if sign(&distance_sq, dec)? == Sign::Zero {
        return Ok(Some(
            if sign(&(radius_sq - &cylinder.radius_sq), dec)? == Sign::Zero {
                Solutions::Contained
            } else {
                Solutions::Points(Vec::new())
            },
        ));
    }
    // The radical plane proves separation without eliminating world coordinates.
    // Nonempty intersections retain the original solver's RUR witnesses.
    let fraction =
        (&(&distance_sq + radius_sq) - &cylinder.radius_sq) / (&Real::from_i64(2) * &distance_sq);
    let height_numerator = radius_sq - &(&fraction.square() * &distance_sq);
    Ok(if sign(&height_numerator, dec)? == Sign::Negative {
        Some(Solutions::Points(Vec::new()))
    } else {
        None
    })
}

// Encloses the circle's axis-aligned box relative to the apex and evaluates
// the double cone's equation there. A strictly signed value proves that no
// point of the circle meets either nappe; zero in the enclosure proves nothing.
// Rational pairs keep the direct solver and its accounting. An algebraic pair
// far from the cone would otherwise eliminate its coefficients first.
fn circle_misses_cone(
    form: &CurveForm,
    surface: &Surface,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let (
        CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ..
        },
        Surface::Cone(cone),
    ) = (form, surface)
    else {
        return Ok(false);
    };
    if center.as_rational().is_some()
        && normal.as_rational().is_some()
        && radius_sq.rational_ref().is_some()
        && surface.has_rational_definition()
    {
        return Ok(false);
    }
    let mut local = cone.clone();
    local.apex = Point3::origin();
    let Some(implicit) = Surface::Cone(local).implicit() else {
        return Ok(false);
    };
    let enclosed = (|| {
        let radius_sq = dec.enclose_at(radius_sq, 128)?;
        let normal_sq = normal.norm_squared();
        let offset = center.minus(&cone.apex);
        let mut bounds = Vec::with_capacity(3);
        for (n, x) in [
            (&normal.x, &offset.x),
            (&normal.y, &offset.y),
            (&normal.z, &offset.z),
        ] {
            let factor = dec.enclose_at(&(&Real::one() - &(&n.square() / &normal_sq)), 128)?;
            let Some(width) = (&radius_sq * &factor).sqrt() else {
                return Ok(None);
            };
            let x = dec.enclose_at(x, 128)?;
            bounds.push((&x - &width).hull(&(&x + &width)));
        }
        implicit.eval_ball(&bounds, dec).map(Some)
    })();
    // The exclusion is optional. An unsettled enclosure keeps the solver.
    match enclosed {
        Ok(value) => Ok(value.is_some_and(|value| !value.contains_zero())),
        Err(refusal) if refusal.code == "numbers.deadline" => Err(refusal),
        Err(_) => Ok(false),
    }
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
    if matches!(form, CurveForm::Branch) && branch_plane_disjoint(geom, c, s, plan.as_deref(), dec)?
    {
        return Ok(Answer::exact(CurveSurface::Disjoint));
    }
    if let Some(answer) = circular_section(geom, c, s, &form, dec)? {
        return Ok(Answer::exact(answer));
    }
    if circle_misses_cone(&form, geom.surface(s), dec)? {
        return Ok(Answer::exact(CurveSurface::Disjoint));
    }
    let implicit = geom.implicit(s).ok_or_else(|| {
        Refusal::new(
            "intersect.parametric-surface",
            "this intersection needs a finite parametric surface solver",
        )
    })?;
    let solutions = if let Some(solutions) = circular_cylinder(&form, geom.surface(s), dec)? {
        solutions
    } else if let CurveForm::Line { origin, dir } = &form
        && !geom.surface(s).has_rational_definition()
    {
        algebraic_line(origin, dir, geom.surface(s), implicit, dec)?
    } else {
        let equations = carrier_equations(geom, c, &form, dec)?;
        solve::intersections(&equations, implicit, dec)?
    };
    let Solutions::Points(points) = solutions else {
        if matches!(geom.surface(s), Surface::Cone(_)) && matches!(form, CurveForm::Line { .. }) {
            return Err(Refusal::new(
                "intersect.partial-containment-required",
                "a complete generator line is only partly on the cone's retained nappe",
            ));
        }
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
                proven_tangent: false,
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
    let recorded_contact =
        geom.facts()
            .holds(FactKind::Tangent, GeomKey::Curve(c), GeomKey::Surface(s));
    for candidate in &mut candidates {
        if recorded_contact && known_tangent(geom, c, s, &candidate.position) {
            dec.note_fact_decision();
            candidate.contact = Contact::Tangent;
            candidate.proven_tangent = true;
            continue;
        }
        let answer = if let CurveForm::Line { dir, .. } = &form {
            match sign(
                &geom
                    .surface(s)
                    .normal_at(&candidate.position, dec)?
                    .dot(dir),
                dec,
            ) {
                Ok(Sign::Zero) if matches!(geom.surface(s), catallax_shapes::Surface::Torus(_)) => {
                    match chart::line_local(&candidate.position, dir, dec) {
                        Ok(local) => contact(
                            std::slice::from_ref(&local),
                            &candidate.position,
                            implicit,
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
                    implicit,
                    dec,
                ),
                Err(e) => Err(Failure::Refused(e)),
            }
        } else {
            contact(&charts, &candidate.position, implicit, dec)
        };
        candidate.proven_tangent = matches!(&answer, Ok(Contact::Tangent));
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
                    sorted[i].proven_tangent &= candidate.proven_tangent;
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
                    candidate.proven_tangent = false;
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
            if candidate.proven_tangent && !candidate.uncertain_branch {
                geom.record_fact(Fact {
                    kind: FactKind::Tangent,
                    a: GeomKey::Curve(c),
                    b: GeomKey::Surface(s),
                    at: Some(GeomKey::Point(point)),
                    origin: FactOrigin::Construction("intersect.curve-surface"),
                });
            }
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
    use catallax_shapes::{Cone, Line, RPoint3, RVector3, Surface, Torus};

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

    // A cone with slope sqrt(3) and apex at the offset, and a horizontal circle
    // of the given squared radius centred at (x, 0, z) from that apex.
    fn circle_and_cone(
        offset: i64,
        x: Real,
        z: i64,
        radius_sq: Real,
    ) -> (GeomStore, CurveId, SurfaceId) {
        let mut dec = Decider::default();
        let shift = |p: [Real; 3]| {
            let [x, y, z] = p;
            Point3::new(
                &x + &Real::from_i64(offset),
                &y - &Real::from_i64(offset),
                &z + &Real::from_i64(offset),
            )
        };
        let mut geom = GeomStore::new();
        let surface = geom.add_surface(Surface::Cone(
            Cone::new(
                shift([Real::zero(), Real::zero(), Real::zero()]),
                vector(0, 0, 1).into(),
                Real::from_i64(3).sqrt(),
                vector(1, 0, 0).into(),
                &mut dec,
            )
            .unwrap(),
        ));
        let curve = geom.add_curve(Curve::Circle(
            Circle::new(
                shift([x, Real::zero(), Real::from_i64(z)]),
                vector(0, 0, 1).into(),
                radius_sq,
                vector(1, 0, 0).into(),
                &mut dec,
            )
            .unwrap(),
        ));
        (geom, curve, surface)
    }

    #[test]
    fn a_circle_away_from_both_nappes_of_an_algebraic_cone_is_excluded_exactly() {
        for offset in [0, 1_000_000] {
            // At height one the cone has radius sqrt(3); the unit circle spans x in [9, 11].
            // The lower nappe is the same distance away at height minus one.
            for z in [1, -1] {
                let (mut geom, curve, surface) =
                    circle_and_cone(offset, Real::from_i64(10), z, Real::one());
                let mut dec = Decider::default();
                let form = query::form(&geom, curve, &mut dec).unwrap();
                assert!(circle_misses_cone(&form, geom.surface(surface), &mut dec).unwrap());
                let answer = intersect(&mut geom, curve, surface, &mut dec).unwrap();
                assert!(answer.marks.is_empty());
                assert!(
                    matches!(answer.value, CurveSurface::Disjoint),
                    "offset {offset}, z {z}"
                );
            }
        }
    }

    #[test]
    fn circles_touching_or_crossing_an_algebraic_cone_keep_the_solver() {
        let three = Real::from_i64(3).sqrt();
        for offset in [0, 1_000_000] {
            // Crossing: the unit circle about x = 2 meets the radius sqrt(3) section twice.
            // Touching: the circle about x = 3 has radius 3 - sqrt(3) and touches at x = sqrt(3).
            for (x, radius_sq, expected) in [
                (Real::from_i64(2), Real::one(), 2),
                (Real::from_i64(3), (&Real::from_i64(3) - &three).square(), 1),
            ] {
                let (mut geom, curve, surface) = circle_and_cone(offset, x, 1, radius_sq);
                let mut dec = Decider::default();
                let form = query::form(&geom, curve, &mut dec).unwrap();
                assert!(!circle_misses_cone(&form, geom.surface(surface), &mut dec).unwrap());
                let answer = intersect(&mut geom, curve, surface, &mut dec).unwrap();
                let CurveSurface::Meets(hits) = answer.value else {
                    panic!("offset {offset}: expected {expected} hits")
                };
                assert_eq!(hits.len(), expected, "offset {offset}");
                if expected == 1 {
                    assert_eq!(hits[0].contact, Contact::Tangent);
                }
            }
        }
    }

    #[test]
    fn rational_circle_and_cone_pairs_keep_the_direct_solver() {
        let mut dec = Decider::default();
        let cone = Surface::Cone(
            Cone::new(
                point(0, 0, 0).into(),
                vector(0, 0, 1).into(),
                r(1).into(),
                vector(1, 0, 0).into(),
                &mut dec,
            )
            .unwrap(),
        );
        let form = CurveForm::Circle {
            center: point(10, 0, 1).into(),
            normal: vector(0, 0, 1).into(),
            radius_sq: r(1).into(),
            ref_dir: vector(1, 0, 0).into(),
        };
        assert!(!circle_misses_cone(&form, &cone, &mut dec).unwrap());
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
