//! Rational ruled draft export using the spline writer's shared chart topology.
use super::*;
use catallax_shapes::ObliqueDraft;
use catallax_topology::DraftedCylinder;

pub(crate) fn build(
    model: &Model,
    body: BodyId,
    options: &CertifiedSplineStepOptions,
    dec: &mut Decider,
) -> Answered<PartDocument> {
    if options.file_tolerance.sign() == Sign::Negative {
        return Err(bad("negative file tolerance"));
    }
    let view = DraftedCylinder::check(model, body, dec)?;
    let (rounded, weights, charts) = round_support(view.support(), options, dec)?;
    let hole = if let Some((support, face, _)) = view.hole() {
        let (inner, inner_weights, inner_charts) = round_support(support, options, dec)?;
        for i in 0..2 {
            if sign(
                &rounded.centers[i].minus(&inner.centers[i]).norm_squared(),
                dec,
            )? != Sign::Zero
            {
                return Err(bad("rounded annular walls have different cap planes"));
            }
        }
        let displacement = Real::from_rational((&rounded.bound + &inner.bound).mul_pow2(1));
        let clearance = view.wall_clearance(dec)?.unwrap();
        if sign(&(clearance - displacement), dec)? != Sign::Positive {
            return Err(bad("rounded annular walls exceed their separation bound"));
        }
        Some(RoundedHole {
            face,
            rounded: inner,
            charts: inner_charts,
            weights: inner_weights,
        })
    } else {
        None
    };
    build_rounded(
        model,
        view.faces(),
        rounded,
        charts,
        4,
        Some(&weights),
        hole,
        options,
        dec,
    )
}

fn round_support(
    support: &ObliqueDraft,
    options: &CertifiedSplineStepOptions,
    dec: &mut Decider,
) -> Result<(Rounded, Vec<Real>, Charts), Refusal> {
    let axis = integer_direction(support.axis())?;
    let reference = integer_direction(support.across())?;
    let minor = integer_direction(&axis.cross(&reference))?;
    let u2 = reference.norm_squared();
    let v2 = minor.norm_squared();
    let n2 = axis.norm_squared();
    let center_grid = Grid::new(options.decimal_places, dec)?;
    let u_grid = Grid::for_direction(options.decimal_places, &reference, dec)?;
    let v_grid = Grid::for_direction(options.decimal_places, &minor, dec)?;
    let n_grid = Grid::for_direction(options.decimal_places, &axis, dec)?;
    let source_centers = support
        .heights()
        .each_ref()
        .map(|h| support.origin().plus(&support.axis().scaled(h)));
    let lower = center_grid.point(&source_centers[0], dec)?;
    let sweep = source_centers[1].minus(&source_centers[0]);
    let height = n_grid.round(&(&sweep.dot(&axis) / &n2), dec)?;
    positive(&height, dec)?;
    let centers = [lower.clone(), lower.plus(&axis.scaled(&height))];
    let mut controls: Vec<[Point3; 2]> = Vec::with_capacity(17);
    let mut weights = Vec::with_capacity(17);
    let carrier = support.carrier();
    let bits = options
        .decimal_places
        .saturating_add(16)
        .saturating_mul(4)
        .max(128)
        .min(dec.policy().max_precision_bits);
    let mut bound = Rational::zero();
    for (a, b) in source_centers.iter().zip(&centers) {
        bound = bound.max(difference_bound(a, b, bits, dec)?);
    }
    for i in 0..17 {
        let mut pair = Vec::with_capacity(2);
        for v in 0..2 {
            let source = &carrier.controls()[i * 2 + v];
            let offset = source.minus(&source_centers[v]);
            let x = u_grid.round(&(&offset.dot(&reference) / &u2), dec)?;
            let y = v_grid.round(&(&offset.dot(&minor) / &v2), dec)?;
            let p = centers[v]
                .plus(&reference.scaled(&x))
                .plus(&minor.scaled(&y));
            pair.push(p);
        }
        let mut pair: [Point3; 2] = pair.try_into().unwrap();
        if i == 16 {
            pair = controls[0].clone();
        }
        for (v, point) in pair.iter().enumerate() {
            bound = bound.max(difference_bound(
                &carrier.controls()[i * 2 + v],
                point,
                bits,
                dec,
            )?);
        }
        controls.push(pair);
        // Squaring the half-angle denominator produces thirds in the
        // Bernstein weights. A common factor three makes every weight an
        // exact terminating integer; no independent weight rounding occurs.
        let weight = &carrier.weights()[i * 2] * &Real::from_i64(3);
        numeral(&weight)?;
        weights.push(weight);
    }
    if bound > options.file_tolerance {
        return Err(Refusal::new(
            "step.tolerance-not-reached",
            format!(
                "shared ruled geometry reaches {bound} mm, above the fixed {} mm tolerance",
                options.file_tolerance
            ),
        ));
    }
    let lean_sq = &support.lean()[0].square() + &support.lean()[1].square();
    let lean = lean_sq.sqrt();
    let margin = [0, 1]
        .into_iter()
        .map(|i| support.amplitude(&support.heights()[i]) - lean.clone())
        .try_fold(None::<Rational>, |old, x| -> Result<_, Refusal> {
            let low = dec.enclose_budgeted_at(&x, bits)?.lo();
            Ok(Some(old.map_or(low.clone(), |v| v.min(low))))
        })?
        .unwrap();
    if bound.mul_pow2(2) >= margin {
        return Err(bad(
            "rounded ruled surface lacks a positive radial correspondence margin",
        ));
    }
    certify_radial(&controls, &weights, &centers, &reference, &minor, dec)?;
    let rounded = Rounded {
        controls,
        centers,
        axis,
        reference,
        bound,
    };
    let charts = open_charts(carrier.knots()[0], 4)?;
    Ok((rounded, weights, charts))
}

pub(super) struct HoleGraph {
    pub(super) rings: [[Id; 2]; 2],
    pub(super) faces: Vec<Field>,
    pub(super) face_map: Vec<StepFaceRegion>,
}

pub(super) fn append_hole(
    graph: &mut Graph,
    hole: &RoundedHole,
    context: Id,
) -> Result<HoleGraph, Refusal> {
    let mut controls: Vec<[Id; 2]> = Vec::with_capacity(hole.rounded.controls.len());
    for (index, pair) in hole.rounded.controls.iter().enumerate() {
        controls.push(if index + 1 == hole.rounded.controls.len() {
            controls[0]
        } else {
            [point(graph, &pair[0])?, point(graph, &pair[1])?]
        });
    }
    let charts = &hole.charts;
    let rows = [&controls[..=charts.cut], &controls[charts.cut..]];
    let weights = [&hole.weights[..=charts.cut], &hole.weights[charts.cut..]];
    let mut sides = Vec::with_capacity(2);
    let mut arcs = Vec::with_capacity(2);
    for i in 0..2 {
        let (multiplicities, knots) = knot_fields(&charts.knots[i])?;
        sides.push(surface(
            graph,
            rows[i],
            4,
            weights[i],
            &multiplicities,
            &knots,
        )?);
        arcs.push([
            curve(graph, rows[i], 0, 4, weights[i], &multiplicities, &knots)?,
            curve(graph, rows[i], 1, 4, weights[i], &multiplicities, &knots)?,
        ]);
    }
    let corners = [controls[0], controls[charts.cut]];
    let vertices =
        corners.map(|pair| pair.map(|row| graph.add("VERTEX_POINT", vec![text("''"), row.into()])));
    let mut rings = [[None; 2]; 2];
    for i in 0..2 {
        for cap in 0..2 {
            let pc = pcurve_line(
                graph,
                sides[i],
                context,
                ["0.".to_owned(), format!("{cap}.")],
                "(1.,0.)",
            );
            rings[i][cap] = Some(edge_curve(
                graph,
                [vertices[i][cap], vertices[1 - i][cap]],
                arcs[i][cap],
                vec![pc],
            ));
        }
    }
    let rings = rings.map(|pair| pair.map(Option::unwrap));
    let mut joins = Vec::with_capacity(2);
    for i in 0..2 {
        let curve = graph.add(
            "B_SPLINE_CURVE_WITH_KNOTS",
            vec![
                text("''"),
                text("1"),
                refs(corners[i]),
                text(".UNSPECIFIED."),
                boolean(false),
                boolean(false),
                text("(2,2)"),
                text("(0.,1.)"),
                text(".UNSPECIFIED."),
            ],
        );
        let parameters = if i == 0 {
            [&charts.domains[0][0], &charts.domains[1][1]]
        } else {
            [&charts.domains[0][1], &charts.domains[1][0]]
        };
        let mut pcurves = Vec::with_capacity(2);
        for (surface, value) in sides.iter().zip(parameters) {
            pcurves.push(pcurve_line(
                graph,
                *surface,
                context,
                [parameter(value)?, "0.".to_owned()],
                "(0.,1.)",
            ));
        }
        joins.push(edge_curve(graph, vertices[i], curve, pcurves));
    }
    let mut faces = Vec::with_capacity(2);
    let mut face_map = Vec::with_capacity(2);
    for (i, domain) in charts.domains.iter().cloned().enumerate() {
        let uses = [
            (joins[i], true),
            (rings[i][1], true),
            (joins[1 - i], false),
            (rings[i][0], false),
        ];
        let name = format!("face_{}_patch_{i}", hole.face.index());
        faces.push(advanced_face(graph, &name, &uses, sides[i], false));
        face_map.push(StepFaceRegion {
            face: hole.face,
            file_region: name,
            file_domain: Some([domain, [Rational::zero(), Rational::one()]]),
        });
    }
    Ok(HoleGraph {
        rings,
        faces,
        face_map,
    })
}

fn choose(n: usize, k: usize) -> i64 {
    (0..k).fold(1, |a, i| a * (n - i) as i64 / (i + 1) as i64)
}
fn product(a: &[Rational], da: [usize; 2], b: &[Rational], db: [usize; 2]) -> Vec<Rational> {
    let d = [da[0] + db[0], da[1] + db[1]];
    let mut result = vec![Rational::zero(); (d[0] + 1) * (d[1] + 1)];
    for i in 0..=da[0] {
        for j in 0..=da[1] {
            for k in 0..=db[0] {
                for l in 0..=db[1] {
                    let scale = Rational::from_ratio(
                        choose(da[0], i) * choose(da[1], j) * choose(db[0], k) * choose(db[1], l),
                        choose(d[0], i + k) * choose(d[1], j + l),
                    )
                    .unwrap();
                    let at = (i + k) * (d[1] + 1) + j + l;
                    result[at] = &result[at]
                        + &(&(&a[i * (da[1] + 1) + j] * &b[k * (db[1] + 1) + l]) * &scale);
                }
            }
        }
    }
    result
}
fn certify_radial(
    controls: &[[Point3; 2]],
    weights: &[Real],
    centers: &[Point3; 2],
    across: &Vector3,
    side: &Vector3,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for span in 0..4 {
        dec.check_deadline()?;
        let mut x = Vec::with_capacity(10);
        let mut y = Vec::with_capacity(10);
        for i in span * 4..=span * 4 + 4 {
            for v in 0..2 {
                let p = controls[i][v].minus(&centers[v]);
                x.push(
                    (&p.dot(across) * &weights[i])
                        .as_rational()
                        .ok_or_else(|| bad("rounded ruled coordinates are not rational"))?,
                );
                y.push(
                    (&p.dot(side) * &weights[i])
                        .as_rational()
                        .ok_or_else(|| bad("rounded ruled coordinates are not rational"))?,
                );
            }
        }
        let derivative = |values: &[Rational]| -> Vec<Rational> {
            (0..8)
                .map(|i| (&values[i + 2] - &values[i]) * Rational::from_i64(4))
                .collect()
        };
        let dx = derivative(&x);
        let dy = derivative(&y);
        let a = product(&x, [4, 1], &dy, [3, 1]);
        let b = product(&y, [4, 1], &dx, [3, 1]);
        if a.iter().zip(b).any(|(a, b)| *a <= b) {
            return Err(bad(
                "rounded ruled chart has no proved strictly increasing polar angle",
            ));
        }
    }
    // Positive angular derivatives on every tensor span, exact closure and
    // the origin-avoiding homotopy prove one simple radial section at every
    // height. Both caps use the same rows. The common positive axial sweep
    // then proves one embedded closed material boundary, despite saddle twist.
    Ok(())
}

pub(super) fn surface(
    graph: &mut Graph,
    rows: &[[Id; 2]],
    degree: usize,
    weights: &[Real],
    multiplicities: &Field,
    knots: &Field,
) -> Result<Id, Refusal> {
    let weighted = weights
        .iter()
        .map(|w| {
            let n = numeral(w)?;
            Ok(Field::List(vec![text(n.clone()), text(n)]))
        })
        .collect::<Result<Vec<_>, Refusal>>()?;
    Ok(graph.add_complex(vec![
        ("BOUNDED_SURFACE", vec![]),
        (
            "B_SPLINE_SURFACE",
            vec![
                text(degree.to_string()),
                text("1"),
                Field::List(rows.iter().map(|pair| refs(*pair)).collect()),
                text(".UNSPECIFIED."),
                boolean(false),
                boolean(false),
                boolean(false),
            ],
        ),
        (
            "B_SPLINE_SURFACE_WITH_KNOTS",
            vec![
                multiplicities.clone(),
                text("(2,2)"),
                knots.clone(),
                text("(0.,1.)"),
                text(".UNSPECIFIED."),
            ],
        ),
        ("GEOMETRIC_REPRESENTATION_ITEM", vec![]),
        ("RATIONAL_B_SPLINE_SURFACE", vec![Field::List(weighted)]),
        ("REPRESENTATION_ITEM", vec![text("''")]),
        ("SURFACE", vec![]),
    ]))
}
pub(super) fn curve(
    graph: &mut Graph,
    rows: &[[Id; 2]],
    row: usize,
    degree: usize,
    weights: &[Real],
    multiplicities: &Field,
    knots: &Field,
) -> Result<Id, Refusal> {
    let weights = weights
        .iter()
        .map(|w| numeral(w).map(text))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(graph.add_complex(vec![
        ("BOUNDED_CURVE", vec![]),
        (
            "B_SPLINE_CURVE",
            vec![
                text(degree.to_string()),
                refs(rows.iter().map(|p| p[row])),
                text(".UNSPECIFIED."),
                boolean(false),
                boolean(false),
            ],
        ),
        (
            "B_SPLINE_CURVE_WITH_KNOTS",
            vec![multiplicities.clone(), knots.clone(), text(".UNSPECIFIED.")],
        ),
        ("CURVE", vec![]),
        ("GEOMETRIC_REPRESENTATION_ITEM", vec![]),
        ("RATIONAL_B_SPLINE_CURVE", vec![Field::List(weights)]),
        ("REPRESENTATION_ITEM", vec![text("''")]),
    ]))
}
