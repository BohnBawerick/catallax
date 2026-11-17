use super::*;

pub(super) fn parameters(
    point: &Vec3,
    frame: &crate::Placement,
    major: &Rational,
    dec: &mut Decider,
) -> Result<[Rational; 2], ImportError> {
    let [x, y, z] = coordinates(point, frame, dec)?;
    let radial = (x.square() + y.square())
        .sqrt()
        .ok_or_else(|| unsupported("torus radial coordinate has no nonnegative enclosure"))?;
    Ok([
        angle(&x, &y)?,
        angle(&(radial - Ball::from_rational(major, BITS)), &z)?,
    ])
}

pub(super) fn unwrap_end(
    curve: &n::NumericalEdge,
    face: &crate::Face,
    forward: bool,
    start: &[Rational; 2],
    mut end: [Rational; 2],
    dec: &mut Decider,
) -> Result<[Rational; 2], ImportError> {
    let [a, b] = if forward {
        curve.parameters.each_ref()
    } else {
        [&curve.parameters[1], &curve.parameters[0]]
    };
    let mut previous = start.clone();
    for i in 1..4 {
        dec.charge_work(1)?;
        let fraction = Rational::from_ratio(i, 4).expect("nonzero denominator");
        let t = a + &(&fraction * &(b - a));
        let p = curve
            .curve
            .evaluate(&Ball::from_rational(&t, BITS), dec)?
            .ok_or_else(|| unsupported("torus boundary has no parameter evaluator"))?
            .each_ref()
            .map(Ball::mid);
        let mut uv = super::parameters(&p, face, dec)?;
        for axis in 0..2 {
            uv[axis] = near(uv[axis].clone(), &previous[axis], dec)?;
        }
        previous = uv;
    }
    for axis in 0..2 {
        end[axis] = near(end[axis].clone(), &previous[axis], dec)?;
    }
    Ok(end)
}

pub(super) fn retained_refusal(
    source: &crate::geometry_source::RetainedGeometry,
    dec: &mut Decider,
) -> Result<(), ImportError> {
    use crate::geometry_source::{Geometry, Surface as SourceSurface};
    if matches!(
        source.source().geometry(source.entity(), dec)?,
        Geometry::Surface(SourceSurface::Torus(_, _, _))
    ) {
        return Err(Refusal::new(
            "import.toroidal-sheet-unverified",
            "a retained horn or spindle torus needs a certificate for its original finite parameter sheet; the inner sheet and singular locus are not discarded",
        ).into());
    }
    Ok(())
}

fn rounded_point(p: &catallax_shapes::Point3, dec: &mut Decider) -> Result<Vec3, ImportError> {
    Ok([
        dec.enclose_budgeted_at(&p.x, BITS)?.mid(),
        dec.enclose_budgeted_at(&p.y, BITS)?.mid(),
        dec.enclose_budgeted_at(&p.z, BITS)?.mid(),
    ])
}

pub(super) fn join_rims(
    definition: &mut n::NumericalDefinition,
    loops: &mut Vec<Vec<n::NumericalTrim>>,
    solid: &Solid,
    face: &crate::Face,
    repairs: &mut Vec<NumericalTopologyRepair>,
    dec: &mut Decider,
) -> Result<(), ImportError> {
    use crate::planar_arithmetic as a;
    let crate::Surface::Torus(frame, major, minor) = &face.surface else {
        unreachable!()
    };
    if loops.len() != 2 || loops.iter().any(|row| row.len() != 1) {
        return Ok(());
    }
    let first = loops[0][0].clone();
    let mut second = loops[1][0].clone();
    let rims = [
        &definition.edges[first.edge],
        &definition.edges[second.edge],
    ];
    if rims
        .iter()
        .any(|edge| edge.vertices[0] != edge.vertices[1] || !matches!(edge.curve, Curve::Circle(_)))
    {
        return Ok(());
    }
    let winding = (0..2)
        .find(|axis| {
            (&first.uv[1][*axis] - &first.uv[0][*axis]).abs() == turn()
                && first.uv[1][1 - axis] == first.uv[0][1 - axis]
        })
        .ok_or_else(|| {
            Refusal::new(
                "import.torus-rim-unverified",
                "a complete torus rim must wind in one isoparametric direction",
            )
        })?;
    let other = 1 - winding;
    if &first.uv[1][winding] - &first.uv[0][winding]
        != &second.uv[0][winding] - &second.uv[1][winding]
        || second.uv[0][other] != second.uv[1][other]
        || near(second.uv[0][winding].clone(), &first.uv[0][winding], dec)? != first.uv[0][winding]
    {
        return Err(Refusal::new(
            "import.torus-rim-unverified",
            "torus rims need opposite winding and anchors on the same isoparametric circle",
        )
        .into());
    }
    let shift = &first.uv[1][winding] - &second.uv[0][winding];
    for uv in &mut second.uv {
        uv[winding] = &uv[winding] + &shift;
    }
    let delta = &second.uv[0][other] - &first.uv[0][other];
    if delta.is_zero() {
        return Err(Refusal::new(
            "import.torus-rim-unverified",
            "distinct torus rims cannot bound a zero-width strip",
        )
        .into());
    }
    let positive =
        (face.same_sense == (first.uv[1][winding] > first.uv[0][winding])) == (winding == 0);
    let mut other_shift = Rational::zero();
    if positive && delta.sign() != Sign::Positive {
        other_shift = turn();
    }
    if !positive && delta.sign() != Sign::Negative {
        other_shift = -turn();
    }
    for uv in &mut second.uv {
        uv[other] = &uv[other] + &other_shift;
    }
    let vertices = [rims[0].vertices[0], rims[1].vertices[0]];
    let start = &definition.vertices[vertices[0]].position;
    let circle = if winding == 0 {
        let offset = a::sub(start, &frame.origin, dec)?;
        let height = a::arithmetic(
            &a::dot(&offset, &frame.axis, dec)?,
            &a::dot(&frame.axis, &frame.axis, dec)?,
            '/',
            dec,
        )?;
        let axial = a::offset(&frame.origin, &frame.axis, &height, '+', dec)?;
        let radial = a::sub(start, &axial, dec)?;
        let ray = vector(&radial);
        let ray = ray.scaled(&(&Real::one() / &ray.norm_squared().sqrt()));
        let center = rounded_point(
            &point(&frame.origin).plus(&ray.scaled(&Real::from(major.clone()))),
            dec,
        )?;
        Circle {
            center: point(&center),
            normal: vector(&a::cross(&radial, &frame.axis, dec)?),
            radius_sq: Real::from(minor * minor),
            ref_dir: vector(&radial),
        }
    } else {
        let at = Ball::from_rational(&first.uv[0][other], BITS);
        let height = minor * &at.sin().mid();
        let radius = major + &(minor * &at.cos().mid());
        let axis = vector(&frame.axis);
        let axis = axis.scaled(&(&Real::one() / &axis.norm_squared().sqrt()));
        let center = rounded_point(
            &point(&frame.origin).plus(&axis.scaled(&Real::from(height))),
            dec,
        )?;
        Circle {
            center: point(&center),
            normal: vector(&frame.axis),
            radius_sq: Real::from(&radius * &radius),
            ref_dir: vector(&frame.ref_dir),
        }
    };
    let seam = definition.edges.len();
    let key = format!("derived:{}:torus-iso-seam", face.name);
    repairs.push(NumericalTopologyRepair {
        solid: solid.name.clone(),
        face: face.name.clone(),
        source_rim_edges: rims.map(|edge| edge.key.clone()),
        source_vertices: vertices.map(|vertex| definition.vertices[vertex].key.clone()),
        derived_seam_edge: key.clone(),
        method: "paired_torus_iso_circle/1".into(),
    });
    definition.edges.push(n::NumericalEdge {
        key,
        curve: Curve::Circle(circle),
        vertices,
        parameters: [first.uv[1][other].clone(), second.uv[0][other].clone()],
    });
    *loops = vec![vec![
        first.clone(),
        n::NumericalTrim {
            edge: seam,
            sense: Sense::Same,
            uv: [first.uv[1].clone(), second.uv[0].clone()],
        },
        second.clone(),
        n::NumericalTrim {
            edge: seam,
            sense: Sense::Reversed,
            uv: [second.uv[1].clone(), first.uv[0].clone()],
        },
    ]];
    Ok(())
}
