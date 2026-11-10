//! Tensor grids bounded against complete rational chart cells, not samples.

use super::*;
use crate::EnvelopeGrid;

/// A bounded display grid together with its original knot-coordinate stations.
#[derive(Clone, Debug)]
pub struct RationalPatchGrid {
    /// Enclosed positions and the two-way file-independent chord bound.
    pub grid: EnvelopeGrid,
    /// Original u and v stations, including every knot and required shared row.
    pub parameters: [Vec<Real>; 2],
    /// Whole-cell unit natural normal boxes, in u-major order, when angular control
    /// was requested. The mesh writer must also check its rounded triangles.
    pub normal_boxes: Option<Vec<[Ball; 3]>>,
    /// Projected embedding margin, including an explicit rounding perturbation check.
    pub embedding: RationalGridEmbedding,
}

/// A strong-monotonicity certificate for the piecewise-affine grid map.
#[derive(Clone, Debug)]
pub struct RationalGridEmbedding {
    lower: Rational,
    projection_norms: [Rational; 2],
    step_lower: [Rational; 2],
}

impl RationalGridEmbedding {
    /// Whether perturbing every source vertex by at most this L1 error keeps
    /// the complete piecewise-affine projection strongly monotone and embedded.
    pub fn accepts_rounding(&self, error: &Rational) -> bool {
        if error < &Rational::zero() {
            return false;
        }
        let twice = error * &Rational::from_i64(2);
        let p = self.projection_norms.each_ref().map(|norm| &twice * norm);
        let diagonal = (&p[0] / &self.step_lower[0]).max(&p[1] / &self.step_lower[1]);
        let off = (&p[0] / &self.step_lower[1] + &p[1] / &self.step_lower[0]).mul_pow2(-1);
        self.lower > diagonal + off
    }
}

fn minimum_step(values: &[Real], bits: u32, dec: &mut Decider) -> Result<Rational, Refusal> {
    let mut lower = None::<Rational>;
    for pair in values.windows(2) {
        let gap = &pair[1] - &pair[0];
        let mut precision = bits;
        let mut enclosed = dec.enclose_budgeted_at(&gap, precision)?;
        while enclosed.lo() <= Rational::zero() && precision < dec.policy().max_precision_bits {
            precision = precision
                .saturating_mul(2)
                .min(dec.policy().max_precision_bits);
            enclosed = dec.enclose_budgeted_at(&gap, precision)?;
        }
        if enclosed.lo() <= Rational::zero() {
            return Err(unproven("a display station gap has no positive enclosure"));
        }
        lower = Some(lower.map_or_else(|| enclosed.lo(), |old| old.min(enclosed.lo())));
    }
    lower.ok_or_else(|| unproven("a display axis has no finite interval"))
}

fn l1(value: &[Ball; 3]) -> Rational {
    value
        .iter()
        .fold(Rational::zero(), |sum, x| sum + magnitude(x))
}

fn normal_radius(
    first: &[[Ball; 3]; 2],
    second: &[[Ball; 3]; 3],
    widths: [&Rational; 2],
) -> Option<Rational> {
    let normal = cross(&first[0], &first[1]);
    let lower_squared = normal.iter().fold(Rational::zero(), |sum, component| {
        let (lo, hi) = (component.lo(), component.hi());
        let lower = if lo > Rational::zero() {
            lo
        } else if hi < Rational::zero() {
            -hi
        } else {
            Rational::zero()
        };
        sum + &lower * &lower
    });
    if lower_squared.is_zero() {
        return None;
    }
    let add = |a: [Ball; 3], b: [Ball; 3]| std::array::from_fn(|i| &a[i] + &b[i]);
    let derivatives = [
        add(cross(&second[0], &first[1]), cross(&first[0], &second[1])),
        add(cross(&second[1], &first[1]), cross(&first[0], &second[2])),
    ];
    // |d(n/|n|)| = |n cross dn| / |n|^2 removes parameter-speed
    // variation. Integrate this complete-cell bound from its centre.
    let variation = &l1(&cross(&normal, &derivatives[0])) * widths[0]
        + &l1(&cross(&normal, &derivatives[1])) * widths[1];
    Some((variation / lower_squared).mul_pow2(-1))
}

fn halve_axis(axis: &[Real]) -> Vec<Real> {
    let mut refined = Vec::with_capacity(axis.len() * 2 - 1);
    refined.push(axis[0].clone());
    for pair in axis.windows(2) {
        refined.push(&(&pair[0] + &pair[1]) / &Real::from_i64(2));
        refined.push(pair[1].clone());
    }
    refined
}

fn enclose_axis(
    values: &[Real],
    domain: &[Rational; 2],
    bits: u32,
    dec: &mut Decider,
) -> Result<Vec<[Rational; 2]>, Refusal> {
    values
        .iter()
        .map(|value| {
            if let Some(q) = value.rational_ref() {
                return Ok([q.clone(), q.clone()]);
            }
            let enclosed = dec.enclose_budgeted_at(value, bits)?;
            Ok([
                enclosed.lo().max(domain[0].clone()),
                enclosed.hi().min(domain[1].clone()),
            ])
        })
        .collect()
}

fn owners(
    values: &[Real],
    knots: &[Rational],
    degree: usize,
    dec: &mut Decider,
) -> Result<Vec<usize>, Refusal> {
    let count = knots.len() - degree - 1;
    let mut span = degree;
    let mut indices = Vec::with_capacity(values.len() - 1);
    for pair in values.windows(2) {
        while span + 1 < count
            && exact_sign(
                &(&pair[0] - &Real::from_rational(knots[span + 1].clone())),
                dec,
            )? != Sign::Negative
        {
            span += 1;
        }
        if exact_sign(
            &(&Real::from_rational(knots[span + 1].clone()) - &pair[1]),
            dec,
        )? == Sign::Negative
        {
            return Err(unproven("a display cell crosses an unowned knot boundary"));
        }
        indices.push(span);
    }
    Ok(indices)
}

impl RationalPatch {
    /// Bounds barycentric chord error over complete cells, retaining all knots
    /// and required stations. Algebraic stations allow exact shared boundaries
    /// with independently parameterised neighbouring faces. Angular refinement
    /// bounds normal variation; the returned boxes certify rounded facets too.
    #[allow(clippy::too_many_arguments)]
    pub fn display_grid(
        &self,
        tolerance: &Rational,
        angular: Option<&Rational>,
        max_edge: Option<&Rational>,
        max_vertices: usize,
        required: [&[Real]; 2],
        dec: &mut Decider,
    ) -> Result<RationalPatchGrid, Refusal> {
        if tolerance <= &Rational::zero()
            || angular.is_some_and(|value| value <= &Rational::zero())
            || max_edge.is_some_and(|value| value <= &Rational::zero())
        {
            return Err(Refusal::new(
                "mesh.invalid-options",
                "patch display controls must be positive",
            ));
        }
        let bits = dec.policy().max_precision_bits.min(160);
        if bits == 0 {
            return Err(unproven("no display enclosure precision is available"));
        }
        let mut axes: [Vec<Real>; 2] = std::array::from_fn(|axis| {
            let mut values = self.carrier().knots()[axis].to_vec();
            values.dedup();
            values.into_iter().map(Real::from_rational).collect()
        });
        for axis in 0..2 {
            for value in required[axis] {
                let [lo, hi] = self.domain()[axis]
                    .each_ref()
                    .map(|q| Real::from_rational(q.clone()));
                if exact_sign(&(value - &lo), dec)? == Sign::Negative
                    || exact_sign(&(&hi - value), dec)? == Sign::Negative
                {
                    return Err(Refusal::new(
                        "mesh.rational-patch-station",
                        "a required shared station leaves the active domain",
                    ));
                }
                insert(&mut axes[axis], value.clone(), dec)?;
            }
        }
        let frame = [
            math::budgeted_vector(&self.data.projection[0], bits, dec)?,
            math::budgeted_vector(&self.data.projection[1], bits, dec)?,
        ];
        let half_tolerance = tolerance.mul_pow2(-1);
        let mut workspace = self.carrier().span_workspace();
        let knots = self.carrier().knots();
        let degrees = self.carrier().degrees();
        let scales = self
            .data
            .widths
            .each_ref()
            .map(|width| Ball::from_rational(width, bits));
        let clip = |span: [usize; 2], domain: Domain| -> Domain {
            std::array::from_fn(|axis| {
                [
                    domain[axis][0].clone().max(knots[axis][span[axis]].clone()),
                    domain[axis][1]
                        .clone()
                        .min(knots[axis][span[axis] + 1].clone()),
                ]
            })
        };
        loop {
            dec.check_deadline()?;
            let vertices = axes[0]
                .len()
                .checked_mul(axes[1].len())
                .ok_or_else(|| Refusal::new("mesh.vertex-limit", "patch grid size overflow"))?;
            if vertices > max_vertices {
                return Err(Refusal::new(
                    "mesh.vertex-limit",
                    "the proven rational patch grid exceeds the requested vertex cap",
                ));
            }
            let boxes = [
                enclose_axis(&axes[0], &self.domain()[0], bits, dec)?,
                enclose_axis(&axes[1], &self.domain()[1], bits, dec)?,
            ];
            let owners = [
                owners(&axes[0], knots[0], degrees[0], dec)?,
                owners(&axes[1], knots[1], degrees[1], dec)?,
            ];
            let mut accepted = true;
            let mut worst = Rational::zero();
            let mut lower = None::<Rational>;
            let mut normal_boxes =
                angular.map(|_| Vec::with_capacity((axes[0].len() - 1) * (axes[1].len() - 1)));
            'cells: for (i, u) in boxes[0].windows(2).enumerate() {
                for (j, v) in boxes[1].windows(2).enumerate() {
                    dec.check_deadline()?;
                    let span = [owners[0][i], owners[1][j]];
                    let domain = clip(
                        span,
                        [
                            [u[0][0].clone(), u[1][1].clone()],
                            [v[0][0].clone(), v[1][1].clone()],
                        ],
                    );
                    let jet =
                        self.carrier()
                            .jet_on_span(span, &domain, bits, dec, &mut workspace)?;
                    let du = &domain[0][1] - &domain[0][0];
                    let dv = &domain[1][1] - &domain[1][0];
                    // Taylor's integral remainder against affine barycentric
                    // interpolation bounds both patch-to-mesh directions.
                    let chord = (&l1(&jet.second[0]) * &(&du * &du)
                        + &l1(&jet.second[1]) * &(&du * &dv) * Rational::from_i64(2)
                        + &l1(&jet.second[2]) * &(&dv * &dv))
                        .mul_pow2(-1);
                    let margin = projection_margin(&jet.first, &frame, &scales, bits);
                    if chord > half_tolerance || margin <= Rational::zero() {
                        accepted = false;
                        break 'cells;
                    }
                    lower =
                        Some(lower.map_or_else(|| margin.clone(), |old| old.min(margin.clone())));
                    if let Some(limit) = max_edge {
                        let edge = &l1(&jet.first[0]) * &du + &l1(&jet.first[1]) * &dv;
                        if edge > *limit {
                            accepted = false;
                            break 'cells;
                        }
                    }
                    if let Some(limit) = angular {
                        let Some(radius) = normal_radius(&jet.first, &jet.second, [&du, &dv])
                        else {
                            accepted = false;
                            break 'cells;
                        };
                        if &radius * &Rational::from_i64(8) > *limit {
                            accepted = false;
                            break 'cells;
                        }
                        let centre = midpoint(&domain).map(|value| [value.clone(), value]);
                        let first = self.carrier().derivatives_on_span(
                            span,
                            &centre,
                            bits,
                            dec,
                            &mut workspace,
                        )?;
                        let normal = cross(&first[0], &first[1]);
                        let length = dot(&normal, &normal)
                            .sqrt()
                            .ok_or_else(|| unproven("display normal length was not separated"))?;
                        let divide = |value: &Ball| {
                            value.checked_div(&length).ok_or_else(|| {
                                unproven("display unit normal denominator was not separated")
                            })
                        };
                        let unit = [
                            divide(&normal[0])?,
                            divide(&normal[1])?,
                            divide(&normal[2])?,
                        ];
                        let error = Ball::from_rational(&radius, bits);
                        let error = error.hull(&(-&error));
                        let normal = unit.each_ref().map(|component| component + &error);
                        let normals = normal_boxes.as_mut().expect("angular boxes requested");
                        debug_assert_eq!(normals.len(), i * (axes[1].len() - 1) + j);
                        normals.push(normal);
                    }
                    worst = worst.max(chord);
                }
            }
            if accepted {
                let mut points = Vec::with_capacity(vertices);
                let mut rounding = Rational::zero();
                for (i, u) in boxes[0].iter().enumerate() {
                    for (j, v) in boxes[1].iter().enumerate() {
                        let span = [
                            owners[0][i.min(owners[0].len() - 1)],
                            owners[1][j.min(owners[1].len() - 1)],
                        ];
                        let domain = clip(span, [u.clone(), v.clone()]);
                        let point = self.carrier().value_on_span(
                            span,
                            &domain,
                            bits,
                            dec,
                            &mut workspace,
                        )?;
                        rounding = rounding.max(
                            point
                                .iter()
                                .fold(Rational::zero(), |sum, component| sum + component.rad()),
                        );
                        points.push(point);
                    }
                }
                let chord = worst + rounding;
                if chord > *tolerance {
                    return Err(unproven(
                        "display coordinate uncertainty exceeds the requested chord tolerance",
                    ));
                }
                let embedding = RationalGridEmbedding {
                    lower: lower.expect("a finite tensor grid has cells"),
                    projection_norms: frame.each_ref().map(l1),
                    step_lower: [
                        &minimum_step(&axes[0], bits, dec)? / &self.data.widths[0],
                        &minimum_step(&axes[1], bits, dec)? / &self.data.widths[1],
                    ],
                };
                return Ok(RationalPatchGrid {
                    grid: EnvelopeGrid {
                        columns: axes[0].len(),
                        rows: axes[1].len(),
                        points,
                        chord,
                    },
                    parameters: axes,
                    normal_boxes,
                    embedding,
                });
            }
            let next = axes[0]
                .len()
                .checked_mul(2)
                .and_then(|n| n.checked_sub(1))
                .and_then(|u| {
                    axes[1]
                        .len()
                        .checked_mul(2)
                        .and_then(|n| n.checked_sub(1))
                        .and_then(|v| u.checked_mul(v))
                });
            if next.is_none_or(|count| count > max_vertices) {
                return Err(Refusal::new(
                    "mesh.vertex-limit",
                    "rational patch refinement exceeds the requested vertex cap",
                ));
            }
            axes = axes.each_ref().map(|axis| halve_axis(axis));
        }
    }
}
