//! Exact interior probes of a retained section arc, without a fitted trace.

use super::*;
use catallax_numbers::ZeroProof;
use catallax_shapes::{CurveId, PatchSectionBranch, Plane, SectionHits, SectionSurface};

fn unresolved(detail: &str) -> Refusal {
    refusal("boolean.patch-section-probe", detail)
}

fn distinct(model: &Model, a: PointId, b: PointId, dec: &mut Decider) -> Result<bool, Refusal> {
    match model.geom().same_point(a, b, dec) {
        Outcome::Proven(_) => Ok(true),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => Ok(false),
        Outcome::ProvenCoincident(_) => Err(unresolved("arc endpoints need exact point identity")),
        Outcome::NotExact(error) => Err(Refusal::from(error)),
        Outcome::Refused(error) => Err(error),
    }
}

fn retain(
    model: &mut Model,
    curve: CurveId,
    from: Option<PointId>,
    to: Option<PointId>,
    position: Point3,
    dec: &mut Decider,
) -> Result<Option<PointId>, Refusal> {
    let point = model
        .geom_mut()
        .add_point(PointDef::Explicit(position.clone()));
    for end in [from, to].into_iter().flatten() {
        if !distinct(model, point, end, dec)? {
            return Ok(None);
        }
    }
    if !crate::geometry::on_arc(
        model,
        curve,
        from,
        to,
        point,
        dec,
        "boolean.patch-section-probe",
    )? {
        return Ok(None);
    }
    Ok(Some(model.geom_mut().add_point(PointDef::Meet(Meet {
        of: MeetOf::Curves(curve, curve),
        position: Position::Closed(position),
    }))))
}

pub(super) fn sample(
    model: &mut Model,
    curve: CurveId,
    branch: &PatchSectionBranch,
    from: Option<PointId>,
    to: Option<PointId>,
    dec: &mut Decider,
) -> Result<PointId, Refusal> {
    if (from.is_none() || to.is_none()) && !(from.is_none() && to.is_none() && branch.is_closed()) {
        return Err(unresolved("an open section arc needs two finite endpoints"));
    }
    for end in [from, to].into_iter().flatten() {
        if !branch.contains(&model.geom().position(end), dec)? {
            return Err(unresolved(
                "an arc endpoint is outside its retained section component",
            ));
        }
    }
    for position in branch.witnesses() {
        if let Some(point) = retain(model, curve, from, to, position, dec)? {
            return Ok(point);
        }
    }
    let (Some(from), Some(to)) = (from, to) else {
        return Err(unresolved(
            "the closed section has no interior cell witness",
        ));
    };
    let a = model.geom().position(from);
    let b = model.geom().position(to);
    let a = [&a.x, &a.y, &a.z];
    let b = [&b.x, &b.y, &b.z];
    let limit = dec.policy().max_precision_bits;
    if limit == 0 {
        return Err(unresolved("a section arc probe has no enclosure precision"));
    }
    let mut bits = limit.min(64);
    loop {
        for axis in 0..3 {
            let first = dec.enclose_budgeted_at(a[axis], bits)?;
            let last = dec.enclose_budgeted_at(b[axis], bits)?;
            let gap = if first.hi() < last.lo() {
                Some((first.hi(), last.lo()))
            } else if last.hi() < first.lo() {
                Some((last.hi(), first.lo()))
            } else {
                None
            };
            let Some((lo, hi)) = gap else { continue };
            dec.reserve_exact_arithmetic(&lo, &hi, false)?;
            let sum = &lo + &hi;
            let half = Rational::from_ratio(1, 2).unwrap();
            dec.reserve_exact_arithmetic(&sum, &half, true)?;
            let at = Real::from_rational(sum * half);
            let mut origin = [Real::zero(), Real::zero(), Real::zero()];
            origin[axis] = at;
            let [x, y, z] = origin;
            let cutting = SectionSurface::Plane(Plane {
                origin: Point3::new(x, y, z),
                normal: unit(axis),
                ref_dir: unit((axis + 1) % 3),
            });
            // Continuity forces the retained arc across this strict coordinate
            // gap. The complete query still checks every possible crossing.
            let SectionHits::Points(points) = branch.surface_hits(&cutting, dec)? else {
                return Err(unresolved(
                    "a strict endpoint coordinate cut contained the whole branch",
                ));
            };
            for position in points {
                if let Some(point) = retain(model, curve, Some(from), Some(to), position, dec)? {
                    return Ok(point);
                }
            }
            return Err(unresolved(
                "the complete coordinate cut found no proved interior arc point",
            ));
        }
        if bits == limit {
            break;
        }
        bits = bits.saturating_mul(2).min(limit);
    }
    Err(unresolved(
        "no coordinate enclosure separated the section arc endpoints",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{StrictOrder, StrictSign, with_operation_budget};
    use catallax_shapes::{BSplineSurface, Curve, RationalAtlas};
    fn q(n: i64, d: i64) -> Real {
        Rational::from_ratio(n, d).unwrap().into()
    }
    fn exact_zero(value: &Real, dec: &mut Decider) -> Result<(), Refusal> {
        match dec.sign(value) {
            Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => Ok(()),
            other => Err(Refusal::new("test.exact-zero", format!("{other:?}"))),
        }
    }
    fn fixture(short: bool, dec: &mut Decider) -> Result<(), Refusal> {
        let mut controls = Vec::new();
        let mut weights = Vec::new();
        for (x, y, w) in [
            (q(1, 1), q(0, 1), q(1, 1)),
            (q(1, 1), q(1, 2), q(1, 1)),
            (q(3, 5), q(4, 5), q(5, 4)),
        ] {
            for z in [Real::zero(), Real::one()] {
                controls.push(Point3::new(x.clone(), y.clone(), z));
                weights.push(w.clone());
            }
        }
        let knots = |values: &[i64]| values.iter().map(|n| Rational::from_i64(*n)).collect();
        let source = BSplineSurface::new(
            [2, 1],
            [knots(&[0, 0, 0, 1, 1, 1]), knots(&[0, 0, 1, 1])],
            [3, 2],
            controls,
            weights,
            dec,
        )?;
        let mut atlases = RationalAtlas::partition(source, [false, false], dec)?;
        assert_eq!(atlases.len(), 1);
        let cutting = SectionSurface::Plane(Plane {
            origin: Point3::new(q(0, 1), q(0, 1), q(1, 2)),
            normal: Vector3::new(q(0, 1), q(0, 1), q(1, 1)),
            ref_dir: Vector3::new(q(1, 1), q(0, 1), q(0, 1)),
        });
        let section = atlases.pop().unwrap().sections(&cutting, dec)?;
        assert_eq!(section.branches().len(), 1);
        let branch = section.branches()[0].clone();
        assert!(!branch.is_closed());
        let ends = if short {
            [
                Point3::new(q(15, 17), q(8, 17), q(1, 2)),
                Point3::new(q(4, 5), q(3, 5), q(1, 2)),
            ]
        } else {
            branch.endpoints()
        };
        for p in &ends {
            assert!(branch.contains(p, dec)?);
        }
        let mut model = Model::new();
        let curve = model
            .geom_mut()
            .add_curve(Curve::PatchSection(branch.clone()));
        let mut ends = ends.map(|p| model.geom_mut().add_point(PointDef::Explicit(p)));
        match catallax_intersect::order_on_curve(model.geom(), curve, ends[0], ends[1], dec) {
            Outcome::Proven(StrictOrder::Less) => (),
            Outcome::Proven(StrictOrder::Greater) => ends.swap(0, 1),
            other => return Err(Refusal::new("test.endpoint-order", format!("{other:?}"))),
        }
        let sample =
            super::super::curve_sample(&mut model, curve, Some(ends[0]), Some(ends[1]), dec)?;
        let position = model.geom().position(sample);
        assert!(branch.contains(&position, dec)?);
        exact_zero(
            &(&position.x.square() + &position.y.square() - Real::one()),
            dec,
        )?;
        exact_zero(&(&position.z - &q(1, 2)), dec)?;
        let (lo, hi) = if short {
            (q(4, 5), q(15, 17))
        } else {
            (q(3, 5), q(1, 1))
        };
        assert!(matches!(
            dec.sign(&(&position.x - &lo)),
            Outcome::Proven(StrictSign::Positive)
        ));
        assert!(matches!(
            dec.sign(&(&hi - &position.x)),
            Outcome::Proven(StrictSign::Positive)
        ));
        Ok(())
    }
    #[test]
    fn rational_cylinder_whole_and_short_arcs_have_exact_interior_probes() {
        for short in [false, true] {
            let (answer, work) =
                with_operation_budget(20_000_000, || fixture(short, &mut Decider::default()));
            assert!(answer.is_ok());
            assert!(!work.exhausted);
        }
    }
}
