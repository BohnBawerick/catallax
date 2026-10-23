use crate::{
    CurveForm, CurveTrace, TracePiece,
    chart::{self, Chart},
    math::*,
    query,
};
use catallax_numbers::{Answer, Answered, Decider, Real, Refusal, Sign};
use catallax_shapes::{BSplineCurve, CurveId, GeomStore, Point3, PointId};

fn ends(
    geom: &GeomStore,
    c: CurveId,
    from: Option<PointId>,
    to: Option<PointId>,
) -> (Point3, Point3) {
    match (from, to) {
        (Some(a), Some(b)) => (geom.position(a), geom.position(b)),
        (None, None) => {
            let p = geom.curve(c).parameter_origin(geom);
            (p.clone(), p)
        }
        _ => panic!("arc endpoints must both be present or both absent"),
    }
}
fn segments(
    charts: &[Chart],
    from: &Point3,
    to: &Point3,
    whole: bool,
    dec: &mut Decider,
) -> Result<Vec<(usize, Real, Real)>, Failure> {
    let (i, a) = chart::locate(charts, from, dec)?;
    let (j, b) = if whole {
        (i, a.clone())
    } else {
        chart::locate(charts, to, dec)?
    };
    let full = whole || sign(&from.minus(to).norm_squared(), dec)? == Sign::Zero;
    if !full && i == j && sign(&(&b - &a), dec)? == Sign::Positive {
        return Ok(vec![(i, a, b)]);
    }
    let mut pieces = vec![];
    if sign(&(&Real::one() - &a), dec)? == Sign::Positive {
        pieces.push((i, a.clone(), Real::one()));
    }
    let mut index = (i + 1) % charts.len();
    while index != j {
        pieces.push((index, Real::zero(), Real::one()));
        index = (index + 1) % charts.len();
    }
    if sign(&b, dec)? == Sign::Positive {
        pieces.push((j, Real::zero(), b));
    }
    Ok(pieces)
}
fn open_ends(
    geom: &GeomStore,
    c: CurveId,
    from: Option<PointId>,
    to: Option<PointId>,
    dec: &mut Decider,
) -> Result<(Point3, Point3), Refusal> {
    let (a, b) = (
        from.expect("open curve needs start"),
        to.expect("open curve needs end"),
    );
    match query::order(geom, c, a, b, dec) {
        catallax_numbers::Outcome::Proven(catallax_numbers::StrictOrder::Less) => {
            Ok((geom.position(a), geom.position(b)))
        }
        catallax_numbers::Outcome::Proven(_) | catallax_numbers::Outcome::ProvenCoincident(_) => {
            panic!("open curve endpoints must be ordered")
        }
        catallax_numbers::Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        catallax_numbers::Outcome::Refused(e) => Err(e),
    }
}
fn spline_ranges(
    geom: &GeomStore,
    c: CurveId,
    curve: &BSplineCurve,
    from: Option<PointId>,
    to: Option<PointId>,
    dec: &mut Decider,
) -> Result<Vec<(usize, Real, Real)>, Refusal> {
    let closed = crate::curve_is_closed(geom, c, dec)?;
    let [lo, hi] = curve.domain().map(|q| Real::from_rational(q.clone()));
    let (a, b, whole) = match (from, to) {
        (None, None) if closed => (lo.clone(), hi.clone(), false),
        (Some(a), Some(b)) => {
            let start = query::spline_parameter(curve, &geom.position(a), dec)?
                .value()
                .clone();
            let end = if a == b {
                start.clone()
            } else {
                query::spline_parameter(curve, &geom.position(b), dec)?
                    .value()
                    .clone()
            };
            let equal = sign(&(&end - &start), dec)? == Sign::Zero;
            if equal && !closed {
                return Err(Refusal::new(
                    "intersect.spline-arc",
                    "open spline arc needs distinct ordered endpoints",
                ));
            }
            (start, end, equal)
        }
        _ => {
            return Err(Refusal::new(
                "intersect.spline-arc",
                "an open spline needs both finite arc endpoints",
            ));
        }
    };
    let order = sign(&(&b - &a), dec)?;
    if !closed && order != Sign::Positive {
        return Err(Refusal::new(
            "intersect.spline-arc",
            "open spline endpoints oppose its natural parameter direction",
        ));
    }
    let mut pieces = Vec::new();
    let mut append = |from: &Real, to: &Real| -> Result<(), Refusal> {
        for (span, interval) in curve.spans() {
            dec.check_deadline()?;
            let lower = Real::from_rational(interval[0].clone());
            let upper = Real::from_rational(interval[1].clone());
            let a = if sign(&(from - &lower), dec)? == Sign::Positive {
                from.clone()
            } else {
                lower
            };
            let b = if sign(&(to - &upper), dec)? == Sign::Negative {
                to.clone()
            } else {
                upper
            };
            if sign(&(&b - &a), dec)? == Sign::Positive {
                pieces.push((span, a, b));
            }
        }
        Ok(())
    };
    if whole || order == Sign::Negative {
        append(&a, &hi)?;
        append(&lo, &b)?;
    } else {
        append(&a, &b)?;
    }
    Ok(pieces)
}
pub(super) fn point(
    geom: &GeomStore,
    c: CurveId,
    from: Option<PointId>,
    to: Option<PointId>,
    dec: &mut Decider,
) -> Answered<Point3> {
    let form = query::form(geom, c, dec)?;
    if let CurveForm::BSpline(curve) = &form {
        let ranges = spline_ranges(geom, c, curve, from, to, dec)?;
        let (_, a, b) = ranges.first().ok_or_else(|| {
            Refusal::new(
                "intersect.spline-arc",
                "spline arc has no positive parameter span",
            )
        })?;
        return Ok(Answer::exact(
            curve.point_at(&((a + b) / Real::from_i64(2)), dec)?,
        ));
    }
    if matches!(form, CurveForm::Line { .. }) {
        let (a, b) = open_ends(geom, c, from, to, dec)?;
        return Ok(Answer::exact(midpoint(&a, &b)));
    }
    if matches!(
        form,
        CurveForm::Parabola { .. } | CurveForm::Hyperbola { .. }
    ) {
        let (a, b) = open_ends(geom, c, from, to, dec)?;
        let chart = chart::conic_arc(&form, &a, &b, dec)?;
        return Ok(Answer::exact(
            chart.position(&(Real::one() / Real::from_i64(2))),
        ));
    }
    let charts = query::charts(geom, c, dec).map_err(Refusal::from)?;
    let (a, b) = ends(geom, c, from, to);
    let pieces = segments(&charts, &a, &b, from == to, dec).map_err(Refusal::from)?;
    let (i, a, b) = pieces.into_iter().next().expect("nonempty arc");
    let mid = (a + b) / Real::from_i64(2);
    Ok(Answer::exact(charts[i].position(&mid)))
}
pub(super) fn trace(
    geom: &GeomStore,
    c: CurveId,
    from: Option<PointId>,
    to: Option<PointId>,
    dec: &mut Decider,
) -> Answered<CurveTrace> {
    let form = query::form(geom, c, dec)?;
    if let CurveForm::BSpline(curve) = &form {
        let ranges = spline_ranges(geom, c, curve, from, to, dec)?;
        let mut pieces = Vec::with_capacity(ranges.len());
        for (span, from, to) in ranges {
            pieces.push(TracePiece::new(chart::spline_span(
                curve, span, from, to, dec,
            )?));
        }
        return Ok(Answer::exact(CurveTrace { pieces }));
    }
    if matches!(form, CurveForm::Line { .. }) {
        let (a, b) = open_ends(geom, c, from, to, dec)?;
        return Ok(Answer::exact(CurveTrace {
            pieces: vec![TracePiece::new(chart::segment(&a, &b, dec)?)],
        }));
    }
    if matches!(
        form,
        CurveForm::Parabola { .. } | CurveForm::Hyperbola { .. }
    ) {
        let (a, b) = open_ends(geom, c, from, to, dec)?;
        let chart = chart::conic_arc(&form, &a, &b, dec)?;
        return Ok(Answer::exact(CurveTrace {
            pieces: vec![TracePiece::new(chart.coords)],
        }));
    }
    let charts = query::charts(geom, c, dec).map_err(Refusal::from)?;
    let (a, b) = ends(geom, c, from, to);
    let segments = segments(&charts, &a, &b, from == to, dec).map_err(Refusal::from)?;
    let mut pieces = vec![];
    for (i, a, b) in segments {
        pieces.push(TracePiece::new(charts[i].restricted(a, b, dec)?));
    }
    Ok(Answer::exact(CurveTrace { pieces }))
}
