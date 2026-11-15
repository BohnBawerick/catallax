//! Complete finite plane points from the retained graph, without contact labels.

use super::*;
use crate::rational_patch::boundary_rows;
use crate::{Plane, SectionHits};

fn missing(detail: &str) -> Refusal {
    Refusal::new("intersect.patch-pair-plane", detail)
}

fn refinable(error: &Refusal) -> bool {
    error.code == "intersect.patch-pair-trace-jacobian"
        || (error.code == "intersect.patch-pair-trace"
            && error.detail == "graph homogeneous denominator needs a positive enclosure")
}

fn equation(span: &Span, plane: &Plane) -> RealPoly {
    let origin = [&plane.origin.x, &plane.origin.y, &plane.origin.z];
    let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z];
    (0..3).fold(RealPoly::zero(2), |sum, axis| {
        sum + (&span.coordinates[axis] - &span.coordinates[3].scale(origin[axis]))
            .scale(normal[axis])
    })
}

fn value(point: &Point3, plane: &Plane, dec: &mut Decider) -> Result<Sign, Refusal> {
    exact_sign(&point.minus(&plane.origin).dot(&plane.normal), dec)
}

fn store(points: &mut Vec<Point3>, point: Point3, dec: &mut Decider) -> Result<(), Refusal> {
    for previous in points.iter() {
        if same_point(previous, &point, dec)? {
            return Ok(());
        }
    }
    points.push(point);
    Ok(())
}

fn projected(jet: &[Ball; 3], normal: &[Ball; 3]) -> Ball {
    &(&jet[0] * &normal[0]) + &(&(&jet[1] * &normal[1]) + &(&jet[2] * &normal[2]))
}

fn split(
    panel: [Rational; 2],
    depth: usize,
    pending: &mut Vec<([Rational; 2], usize)>,
) -> Result<(), Refusal> {
    if depth == 16 {
        return Err(missing(
            "a complete graph panel needs further finite contact isolation",
        ));
    }
    let middle = (&panel[0] + &panel[1]).mul_pow2(-1);
    pending.push(([middle.clone(), panel[1].clone()], depth + 1));
    pending.push(([panel[0].clone(), middle], depth + 1));
    Ok(())
}

fn root(
    piece: &PairTracePiece,
    cut: &RealPoly,
    panel: &Ball,
    dec: &mut Decider,
) -> Result<Point3, Refusal> {
    let at = piece.parameters(panel, dec)?;
    let arc = piece.arc();
    let domain = [
        arc.first_cell[0].clone(),
        arc.first_cell[1].clone(),
        arc.second_cell[0].clone(),
        arc.second_cell[1].clone(),
    ];
    let panel_width = panel.hi() - panel.lo();
    let bounds: Vec<_> = at
        .iter()
        .zip(&domain)
        .map(|(value, [lo, hi])| {
            let width = hi - lo;
            let padding = &(&width * &panel_width) * &Rational::one().mul_pow2(-24);
            [
                std::cmp::max(lo.clone(), &value.lo() - &padding),
                std::cmp::min(hi.clone(), &value.hi() + &padding),
            ]
        })
        .collect();
    let mut equations: Vec<_> = piece
        .system
        .polynomials
        .iter()
        .map(|equation| {
            RealPoly::from_terms(
                4,
                equation
                    .terms()
                    .into_iter()
                    .map(|(m, c)| (m.exponents()[..4].to_vec(), c))
                    .collect(),
            )
        })
        .collect();
    equations.push(lift(cut, &[0, 1], 4));
    // Opposed exact endpoint signs and a strict projected derivative supply
    // existence in this graph panel. The joint root certificate selects it.
    let at = RealPoly::simple_algebraic_root_in_reconstructed(
        &equations,
        &bounds,
        &[0, 1],
        panel.precision(),
        dec,
    )?;
    Ok(piece.first().point(&[at[0].clone(), at[1].clone()]))
}

impl PatchPairCurve {
    /// Complete finite points on a plane, or whole-component containment.
    /// No two-sided contact label is inferred at a finite endpoint. Partial
    /// containment and contacts without complete panel coverage refuse.
    pub fn plane_hits(&self, plane: &Plane, dec: &mut Decider) -> Result<SectionHits, Refusal> {
        crate::math::validate_point(&plane.origin, dec)?;
        crate::math::validate_vector(&plane.normal, dec)?;
        if exact_sign(&plane.normal.norm_squared(), dec)? != Sign::Positive {
            return Err(missing("a plane needs a nonzero normal"));
        }
        let bits = 128.min(dec.policy().max_precision_bits);
        if bits == 0 {
            return Err(missing("no plane-query enclosure precision is available"));
        }
        let normal = [
            dec.enclose_budgeted_at(&plane.normal.x, bits)?,
            dec.enclose_budgeted_at(&plane.normal.y, bits)?,
            dec.enclose_budgeted_at(&plane.normal.z, bits)?,
        ];
        let offset = dec.enclose_budgeted_at(
            &Vector3::new(
                plane.origin.x.clone(),
                plane.origin.y.clone(),
                plane.origin.z.clone(),
            )
            .dot(&plane.normal),
            bits,
        )?;
        let mut points = Vec::new();
        let mut contained = 0;
        for piece in self.pieces.iter() {
            let cut = equation(piece.first(), plane);
            let arc = piece.arc();
            let second = piece
                .definition
                .second
                .span(arc.second_chart, arc.second_span);
            let second_cut = equation(second, plane);
            if cut.is_zero(dec)? || second_cut.is_zero(dec)? {
                contained += 1;
                continue;
            }
            let mut border_only = None;
            if piece.clipped_ends.is_none() {
                border_only = boundary_rows(&cut, &arc.first_cell, dec)?;
                if border_only.is_none() {
                    border_only = boundary_rows(&second_cut, &arc.second_cell, dec)?;
                }
            }
            if let Some(rows) = border_only {
                // The complete cell-pair certificate accounts for every
                // border meeting as one of this arc's owned endpoints.
                if rows.rows().next().is_some() {
                    for end in 0..2 {
                        let point = endpoint(piece, end);
                        if value(&point, plane, dec)? == Sign::Zero {
                            store(&mut points, point, dec)?;
                        }
                    }
                }
                continue;
            }
            let mut cache = BTreeMap::new();
            let mut pending = vec![([Rational::zero(), Rational::one()], 0)];
            while let Some((panel, depth)) = pending.pop() {
                let interval = Ball::from_rational(&panel[0], bits)
                    .hull(&Ball::from_rational(&panel[1], bits));
                let jet = match piece.taylor(&interval, 0, dec) {
                    Ok(jet) => jet,
                    Err(error) if refinable(&error) => {
                        split(panel, depth, &mut pending)?;
                        continue;
                    }
                    Err(error) => return Err(error),
                };
                if !(&projected(&jet[0], &normal) - &offset).contains_zero() {
                    continue;
                }
                let mut ends = Vec::new();
                for parameter in &panel {
                    if let std::collections::btree_map::Entry::Vacant(entry) =
                        cache.entry(parameter.clone())
                    {
                        entry.insert(piece.point_at(parameter, dec)?);
                    }
                    ends.push(cache[parameter].clone());
                }
                let signs = [value(&ends[0], plane, dec)?, value(&ends[1], plane, dec)?];
                let jets = match piece.taylor(&interval, 1, dec) {
                    Ok(jets) => jets,
                    Err(error) if refinable(&error) => {
                        split(panel, depth, &mut pending)?;
                        continue;
                    }
                    Err(error) => return Err(error),
                };
                if !projected(&jets[1], &normal).contains_zero() {
                    if signs == [Sign::Zero, Sign::Zero] {
                        return Err(missing(
                            "strict monotonicity conflicts with two zero endpoints",
                        ));
                    }
                    for end in 0..2 {
                        if signs[end] == Sign::Zero {
                            store(&mut points, ends[end].clone(), dec)?;
                        }
                    }
                    if signs[0] != Sign::Zero && signs[1] != Sign::Zero && signs[0] != signs[1] {
                        match root(piece, &cut, &interval, dec) {
                            Ok(point) => store(&mut points, point, dec)?,
                            Err(error)
                                if error.code == "numbers.root-not-certified"
                                    || refinable(&error) =>
                            {
                                split(panel, depth, &mut pending)?;
                            }
                            Err(error) => return Err(error),
                        }
                    }
                    continue;
                }
                let jets = match piece.taylor(&interval, 2, dec) {
                    Ok(jets) => jets,
                    Err(error) if refinable(&error) => {
                        split(panel, depth, &mut pending)?;
                        continue;
                    }
                    Err(error) => return Err(error),
                };
                if !projected(&jets[2], &normal).contains_zero() {
                    if signs == [Sign::Zero, Sign::Zero] {
                        for point in ends {
                            store(&mut points, point, dec)?;
                        }
                        continue;
                    }
                    let mut isolated = false;
                    for end in 0..2 {
                        if signs[end] == Sign::Zero {
                            let mut tangent_contact = true;
                            for tangent in self.tangents(&ends[end], dec)? {
                                tangent_contact &=
                                    exact_sign(&tangent.dot(&plane.normal), dec)? == Sign::Zero;
                            }
                            if tangent_contact {
                                store(&mut points, ends[end].clone(), dec)?;
                                isolated = true;
                            }
                        }
                    }
                    if isolated {
                        continue;
                    }
                }
                split(panel, depth, &mut pending)?;
            }
        }
        if contained == self.pieces.len() {
            return Ok(SectionHits::Contained);
        }
        if contained != 0 {
            return Err(missing(
                "partial plane containment needs finite interval ownership",
            ));
        }
        Ok(SectionHits::Points(points))
    }
}
