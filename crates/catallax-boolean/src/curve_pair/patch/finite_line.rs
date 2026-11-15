//! Finite spline meetings using a proved straight component's two line planes.

use super::*;

fn coordinates(point: &Point3) -> [&Real; 3] {
    [&point.x, &point.y, &point.z]
}

pub(super) fn meetings(
    model: &mut Model,
    a: CurveId,
    b: CurveId,
    dec: &mut Decider,
) -> Result<Option<Meetings>, Refusal> {
    if !matches!(
        (model.geom().curve(a), model.geom().curve(b)),
        (Curve::BSpline(_), Curve::BSpline(_))
    ) {
        return Ok(None);
    }
    let (a, b, ends) = if let Some(ends) =
        catallax_intersect::native_spline_line_segment(model.geom(), a, dec)?
    {
        (a, b, ends)
    } else if let Some(ends) = catallax_intersect::native_spline_line_segment(model.geom(), b, dec)?
    {
        (b, a, ends)
    } else {
        return Ok(None);
    };
    let Curve::BSpline(other) = model.geom().curve(b).clone() else {
        unreachable!()
    };
    let origin = coordinates(&ends[0]);
    let last = coordinates(&ends[1]);
    let direction = [
        spline::exact(last[0], &(-origin[0]), false, dec)?,
        spline::exact(last[1], &(-origin[1]), false, dec)?,
        spline::exact(last[2], &(-origin[2]), false, dec)?,
    ];
    let mut selected = None;
    for (axis, value) in direction.iter().enumerate() {
        let sense = sign(value, dec, "boolean.patch-section-meetings")?;
        if sense != Sign::Zero {
            selected = Some((axis, sense));
            break;
        }
    }
    let (axis, sense) = selected.expect("a certified segment has distinct ends");
    let mut planes = Vec::new();
    for other in (0..3).filter(|other| *other != axis) {
        let mut first = vec![0; 3];
        first[other] = 1;
        let mut second = vec![0; 3];
        second[axis] = 1;
        let a = spline::exact(&direction[other], origin[axis], true, dec)?;
        let b = spline::exact(&direction[axis], origin[other], true, dec)?;
        let offset = spline::exact(&a, &(-b), false, dec)?;
        planes.push(RealPoly::from_terms(
            3,
            vec![
                (first, direction[axis].clone()),
                (second, -&direction[other]),
                (vec![0; 3], offset),
            ],
        ));
    }
    let mut candidates = None;
    for plane in &planes {
        if let Some(points) = spline::hits(&other, plane, dec)? {
            candidates = Some(points);
            break;
        }
    }
    let candidates = candidates
        .ok_or_else(|| unresolved("finite straight spline overlap needs interval ownership"))?;
    let mut points = Vec::new();
    for point in candidates {
        let xyz = [point.x.clone(), point.y.clone(), point.z.clone()];
        let mut on_line = true;
        for plane in &planes {
            if sign(
                &plane.eval_real(&xyz),
                dec,
                "boolean.patch-section-meetings",
            )? != Sign::Zero
            {
                on_line = false;
                break;
            }
        }
        if !on_line {
            continue;
        }
        let at = coordinates(&point)[axis];
        let before = spline::exact(at, &(-origin[axis]), false, dec)?;
        let after = spline::exact(last[axis], &(-at), false, dec)?;
        let mut inside = true;
        for delta in [before, after] {
            let side = sign(&delta, dec, "boolean.patch-section-meetings")?;
            if side != Sign::Zero && side != sense {
                inside = false;
                break;
            }
        }
        if inside {
            store(model, a, b, point, &mut points, dec)?;
        }
    }
    Ok(Some(Meetings::Points(points)))
}
