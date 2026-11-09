use super::*;
use catallax_topology::BezierExtrusionBand;

fn point(
    band: &BezierExtrusionBand<'_>,
    u: &Rational,
    v: &Real,
    dec: &mut Decider,
) -> Result<[Real; 3], Refusal> {
    let base = band
        .support()
        .lower_edge()
        .point_at(&Real::from(u.clone()), dec)?;
    let point = base.plus(&band.support().direction().scaled(v));
    Ok([point.x, point.y, point.z])
}

fn enclosed(point: &[Real; 3], bits: u32, dec: &Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        dec.enclose_at(&point[0], bits)?,
        dec.enclose_at(&point[1], bits)?,
        dec.enclose_at(&point[2], bits)?,
    ])
}

pub(super) fn band(
    model: &Model,
    face: FaceId,
    tolerance: &Rational,
    options: &MeshOptions,
    seeds: &mut BTreeMap<EdgeId, Seed>,
    dec: &mut Decider,
) -> Result<Envelope, Refusal> {
    let checked = BezierExtrusionBand::of_face(model, face, dec)?
        .ok_or_else(|| failure("a Bezier side has no checked complete band"))?;
    let bits = BITS.min(dec.policy().max_precision_bits);
    if bits == 0 {
        return Err(Refusal::new(
            "numbers.precision-limit",
            "no Bezier grid query precision",
        ));
    }
    if tolerance <= &Rational::zero() || options.max_vertices < 4 {
        return Err(Refusal::new(
            "mesh.invalid-options",
            "a Bezier grid needs positive tolerance and four vertices",
        ));
    }
    let controls = checked.support().controls();
    let degree = controls.len() - 1;
    let factor = Real::from_i64((degree * (degree - 1)) as i64);
    let mut acceleration = Rational::zero();
    for triple in controls.windows(3) {
        let first = triple[1].minus(&triple[0]);
        let second = triple[2].minus(&triple[1]);
        let acceleration_vector = second.minus(&first).scaled(&factor);
        let mut bound = Rational::zero();
        for coordinate in [
            &acceleration_vector.x,
            &acceleration_vector.y,
            &acceleration_vector.z,
        ] {
            let enclosure = dec.enclose_at(coordinate, bits)?;
            bound = bound + enclosure.lo().abs().max(enclosure.hi().abs());
        }
        acceleration = acceleration.max(bound);
    }
    let mut steps = 1usize;
    let chord = loop {
        dec.check_deadline()?;
        if steps > i64::MAX as usize
            || steps
                .checked_add(1)
                .and_then(|n| n.checked_mul(2))
                .is_none_or(|n| n > options.max_vertices)
        {
            return Err(Refusal::new(
                "mesh.resource-limit",
                "Bezier grid exceeds the requested vertex limit",
            ));
        }
        let width = Rational::from_ratio(1, steps as i64).unwrap();
        // The axial coordinate is affine. Each triangle therefore differs
        // from the ruled patch only by the curve's linear interpolation error.
        let error = (&(&acceleration * &width) * &width).mul_pow2(-3);
        if error <= *tolerance {
            break error;
        }
        steps = steps
            .checked_mul(2)
            .ok_or_else(|| failure("Bezier grid subdivision overflow"))?;
    };
    let mut points = Vec::with_capacity((steps + 1) * 2);
    for index in 0..=steps {
        let u = Rational::from_ratio(index as i64, steps as i64).unwrap();
        for v in checked.axial() {
            points.push(enclosed(&point(&checked, &u, v, dec)?, bits, dec)?);
        }
    }
    let grid = EnvelopeGrid {
        columns: steps + 1,
        rows: 2,
        points,
        chord,
    };
    let mut rows = Vec::with_capacity(4);
    let mut along_u = None;
    for (edge, row, sense, forward) in checked.rows() {
        let points = row_points(&grid, row);
        let count = points.len();
        let edge_definition = model.edge(edge);
        let endpoints = if forward {
            [edge_definition.start, edge_definition.end]
        } else {
            [edge_definition.end, edge_definition.start]
        };
        let mut error = grid.chord.clone();
        for (vertex, index) in endpoints.into_iter().zip([0, count - 1]) {
            let vertex = vertex.ok_or_else(|| failure("a checked Bezier row lost its endpoint"))?;
            let point = model.geom().position(model.vertex(vertex).point);
            let exact = [point.x, point.y, point.z];
            let rounded = chart::position(&enclosed(&exact, BITS, dec)?)?;
            error = error.max(&grid.chord + &chart::residual(&rounded, &points[index]));
        }
        let mut exact = Vec::with_capacity(count);
        for index in 0..count {
            let t = Rational::from_ratio(index as i64, count as i64 - 1).unwrap();
            let [lo, hi] = checked.axial();
            let (u, v) = match row {
                BoundaryRow::Rail => (t, lo.clone()),
                BoundaryRow::Mouth => (t, hi.clone()),
                BoundaryRow::Negative => (Rational::zero(), lo + &((hi - lo) * Real::from(t))),
                BoundaryRow::Positive => (Rational::one(), lo + &((hi - lo) * Real::from(t))),
            };
            exact.push(point(&checked, &u, &v, dec)?);
        }
        if row == BoundaryRow::Rail {
            along_u = Some((sense == Sense::Same) == forward);
        }
        rows.push((edge, row, forward));
        seeds.insert(
            edge,
            Seed {
                points,
                error,
                forward,
                exact,
            },
        );
    }
    Ok(Envelope {
        grid,
        rows,
        along_u: along_u.ok_or_else(|| failure("a Bezier band has no lower rail"))?,
        ends: Ends::Meridians,
    })
}
