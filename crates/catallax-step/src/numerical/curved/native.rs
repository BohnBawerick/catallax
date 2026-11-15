//! AP242 native spline and finite linear pcurve entities.

use super::*;
use catallax_shapes::{BSplineCurve, BSplineSurface};

pub(super) fn token(w: &mut Writer<'_>, value: &Rational) -> Result<Field, Refusal> {
    w.dec.charge_work(1)?;
    Ok(text(
        crate::number::rational(value, w.options.significant_digits).text,
    ))
}

fn knots(w: &mut Writer<'_>, source: &[Rational]) -> Result<(Field, Field), Refusal> {
    let mut multiplicities = Vec::new();
    let mut values = Vec::new();
    let mut start = 0;
    while start < source.len() {
        let mut end = start + 1;
        while end < source.len() && source[end] == source[start] {
            end += 1;
        }
        multiplicities.push(text((end - start).to_string()));
        values.push(token(w, &source[start])?);
        start = end;
    }
    Ok((Field::List(multiplicities), Field::List(values)))
}

pub(super) fn curve(w: &mut Writer<'_>, curve: &BSplineCurve) -> Result<Id, Refusal> {
    curve_dimension(w, curve, false)
}

fn curve_dimension(w: &mut Writer<'_>, curve: &BSplineCurve, uv: bool) -> Result<Id, Refusal> {
    let mut points = Vec::new();
    for point in curve.controls() {
        if uv {
            let values = [&point.x, &point.y]
                .map(|q| {
                    q.as_rational()
                        .ok_or_else(|| bad("unsupported", "decimal UV control is not rational"))
                })
                .into_iter()
                .collect::<Result<Vec<_>, _>>()?;
            let fields = values
                .iter()
                .map(|q| token(w, q))
                .collect::<Result<Vec<_>, _>>()?;
            points.push(
                w.graph
                    .add("CARTESIAN_POINT", vec![text("''"), Field::List(fields)]),
            );
        } else {
            points.push(w.point(&rational_point(point)?)?);
        }
    }
    let weights = curve
        .weights()
        .map(|q| {
            token(
                w,
                &q.as_rational()
                    .ok_or_else(|| bad("unsupported", "decimal weight is not rational"))?,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let (multiplicities, knots) = knots(w, curve.knots())?;
    Ok(w.graph.add_complex(vec![
        ("BOUNDED_CURVE", vec![]),
        (
            "B_SPLINE_CURVE",
            vec![
                text(curve.degree().to_string()),
                refs(points),
                text(".UNSPECIFIED."),
                text(".U."),
                text(".U."),
            ],
        ),
        (
            "B_SPLINE_CURVE_WITH_KNOTS",
            vec![multiplicities, knots, text(".UNSPECIFIED.")],
        ),
        ("CURVE", vec![]),
        ("GEOMETRIC_REPRESENTATION_ITEM", vec![]),
        ("RATIONAL_B_SPLINE_CURVE", vec![Field::List(weights)]),
        ("REPRESENTATION_ITEM", vec![text("''")]),
    ]))
}

pub(super) fn finite_curve(
    w: &mut Writer<'_>,
    spline: &BSplineCurve,
    parameters: &[Rational; 2],
) -> Result<Id, Refusal> {
    let basis = curve(w, spline)?;
    finite_basis(w, basis, parameters)
}

fn finite_basis(w: &mut Writer<'_>, basis: Id, parameters: &[Rational; 2]) -> Result<Id, Refusal> {
    let mut interval = parameters.clone();
    interval.sort();
    let trims = interval
        .iter()
        .map(|q| {
            w.dec.charge_work(1)?;
            let token = crate::number::rational(q, w.options.significant_digits).text;
            Ok(Field::List(vec![text(format!("PARAMETER_VALUE({token})"))]))
        })
        .collect::<Result<Vec<_>, Refusal>>()?;
    Ok(w.graph.add(
        "TRIMMED_CURVE",
        vec![
            text("''"),
            basis.into(),
            trims[0].clone(),
            trims[1].clone(),
            boolean(true),
            text(".PARAMETER."),
        ],
    ))
}

pub(super) fn surface(w: &mut Writer<'_>, surface: &BSplineSurface) -> Result<Id, Refusal> {
    let points = surface
        .controls()
        .iter()
        .map(|p| w.point(&rational_point(p)?).map(Field::from))
        .collect::<Result<Vec<_>, _>>()?;
    let weights = surface
        .weights()
        .iter()
        .map(|q| {
            token(
                w,
                &q.as_rational()
                    .ok_or_else(|| bad("unsupported", "decimal weight is not rational"))?,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let [u, v] = surface.knots();
    let (u_m, u_k) = knots(w, u)?;
    let (v_m, v_k) = knots(w, v)?;
    let [du, dv] = surface.degrees();
    let columns = surface.control_count()[1];
    let grid = |values: Vec<Field>| {
        Field::List(
            values
                .chunks(columns)
                .map(|row| Field::List(row.to_vec()))
                .collect(),
        )
    };
    Ok(w.graph.add_complex(vec![
        ("BOUNDED_SURFACE", vec![]),
        (
            "B_SPLINE_SURFACE",
            vec![
                text(du.to_string()),
                text(dv.to_string()),
                grid(points),
                text(".UNSPECIFIED."),
                boolean(false),
                boolean(false),
                boolean(false),
            ],
        ),
        (
            "B_SPLINE_SURFACE_WITH_KNOTS",
            vec![u_m, v_m, u_k, v_k, text(".UNSPECIFIED.")],
        ),
        ("GEOMETRIC_REPRESENTATION_ITEM", vec![]),
        ("RATIONAL_B_SPLINE_SURFACE", vec![grid(weights)]),
        ("REPRESENTATION_ITEM", vec![text("''")]),
        ("SURFACE", vec![]),
    ]))
}

pub(super) fn pcurve(
    w: &mut Writer<'_>,
    surface: Id,
    context: Id,
    parameters: &[Rational; 2],
    mut uv: [[Rational; 2]; 2],
) -> Result<Id, Refusal> {
    let interval = if parameters[0] < parameters[1] {
        parameters.each_ref()
    } else {
        uv.reverse();
        [&parameters[1], &parameters[0]]
    };
    let mut points = Vec::new();
    for point in &uv {
        let coordinates = point
            .iter()
            .map(|q| token(w, q))
            .collect::<Result<Vec<_>, _>>()?;
        points.push(w.graph.add(
            "CARTESIAN_POINT",
            vec![text("''"), Field::List(coordinates)],
        ));
    }
    let values = interval
        .into_iter()
        .map(|q| token(w, q))
        .collect::<Result<Vec<_>, _>>()?;
    let curve = w.graph.add(
        "B_SPLINE_CURVE_WITH_KNOTS",
        vec![
            text("''"),
            text("1"),
            refs(points),
            text(".UNSPECIFIED."),
            boolean(false),
            boolean(false),
            text("(2,2)"),
            Field::List(values),
            text(".UNSPECIFIED."),
        ],
    );
    let representation = w.graph.add(
        "DEFINITIONAL_REPRESENTATION",
        vec![text("''"), refs([curve]), context.into()],
    );
    Ok(w.graph.add(
        "PCURVE",
        vec![text("''"), surface.into(), representation.into()],
    ))
}

pub(super) fn native_pcurve(
    w: &mut Writer<'_>,
    surface: Id,
    context: Id,
    parameters: &[Rational; 2],
    spline: &BSplineCurve,
) -> Result<Id, Refusal> {
    let basis = curve_dimension(w, spline, true)?;
    let finite = finite_basis(w, basis, parameters)?;
    let representation = w.graph.add(
        "DEFINITIONAL_REPRESENTATION",
        vec![text("''"), refs([finite]), context.into()],
    );
    Ok(w.graph.add(
        "PCURVE",
        vec![text("''"), surface.into(), representation.into()],
    ))
}
