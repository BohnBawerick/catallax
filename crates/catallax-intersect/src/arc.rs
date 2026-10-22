use crate::{
    CurveForm, CurveTrace, TracePiece,
    chart::{self, Chart},
    math::*,
    query,
};
use catallax_numbers::{Answer, Answered, Decider, Real, Refusal, Sign};
use catallax_shapes::{CurveId, GeomStore, Point3, PointId};

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
pub(super) fn point(
    geom: &GeomStore,
    c: CurveId,
    from: Option<PointId>,
    to: Option<PointId>,
    dec: &mut Decider,
) -> Answered<Point3> {
    let form = query::form(geom, c, dec)?;
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
