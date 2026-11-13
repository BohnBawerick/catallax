//! Whole-cell positional bounds, without an embedding claim.

use super::*;
use catallax_shapes::BSplineSurface;

const BITS: u32 = 64;

fn ball(q: &Rational) -> Ball {
    Ball::from_rational(q, BITS)
}
fn magnitude(x: &Ball) -> Rational {
    x.lo().abs().max(x.hi().abs())
}
fn norm(x: &[Ball; 3]) -> Rational {
    x.iter().fold(Rational::zero(), |s, x| s + magnitude(x))
}
fn blend<const N: usize>(a: &[Ball; N], b: &[Ball; N], t: &Ball) -> [Ball; N] {
    std::array::from_fn(|i| &a[i] + &(&(&b[i] - &a[i]) * t))
}
fn split<const N: usize>(net: &[[Ball; N]], t: &Ball) -> (Vec<[Ball; N]>, Vec<[Ball; N]>) {
    let mut work = net.to_vec();
    let mut left = vec![work[0].clone()];
    let mut right = vec![work.last().unwrap().clone()];
    while work.len() > 1 {
        work = work.windows(2).map(|p| blend(&p[0], &p[1], t)).collect();
        left.push(work[0].clone());
        right.push(work.last().unwrap().clone());
    }
    right.reverse();
    (left, right)
}
fn restrict(net: &[[Ball; 4]], lo: &Rational, hi: &Rational) -> Vec<[Ball; 4]> {
    if hi.is_zero() {
        return vec![net[0].clone(); net.len()];
    }
    let left = split(net, &ball(hi)).0;
    if lo.is_zero() {
        left
    } else {
        split(&left, &ball(&(lo / hi))).1
    }
}
#[derive(Clone)]
struct Span {
    domain: [[Rational; 2]; 2],
    net: Vec<[Ball; 4]>,
}
#[derive(Clone)]
pub(super) struct PositionSurface {
    degrees: [usize; 2],
    domain: [[Rational; 2]; 2],
    spans: Vec<Span>,
    points: std::cell::RefCell<BTreeMap<[Rational; 4], [Ball; 3]>>,
}
pub(super) struct Controls<'a> {
    pub tolerance: &'a Rational,
    pub max_vertices: usize,
    pub displacement: &'a Rational,
}

impl PositionSurface {
    pub(super) fn new(surface: &BSplineSurface, dec: &Decider) -> Result<Self, Refusal> {
        let mut spans = Vec::new();
        for index in surface.span_indices() {
            let domain = std::array::from_fn(|a| {
                [
                    surface.knots()[a][index[a]].clone(),
                    surface.knots()[a][index[a] + 1].clone(),
                ]
            });
            let net = surface
                .span_bezier_controls(index, dec)?
                .iter()
                .map(|p| {
                    let [x, y, z, w] = p.each_ref().map(|v| dec.enclose_at(v, BITS));
                    Ok([x?, y?, z?, w?])
                })
                .collect::<Result<Vec<_>, Refusal>>()?;
            spans.push(Span { domain, net });
        }
        Ok(Self {
            degrees: surface.degrees(),
            domain: surface.domain().map(|a| a.map(Clone::clone)),
            spans,
            points: std::cell::RefCell::new(BTreeMap::new()),
        })
    }
    fn net_on(
        &self,
        span: &Span,
        domain: &[[Rational; 2]; 2],
        dec: &mut Decider,
    ) -> Result<Vec<[Ball; 4]>, Refusal> {
        let [p, q] = self.degrees;
        // Two de Casteljau restrictions per row and column. At 64 bits each
        // four-coordinate blend pays four interval multiplies and eight adds.
        let blends = (p + 1).saturating_mul(q + 1).saturating_mul(p + q);
        dec.charge_work((blends as u64).saturating_mul(64))?;
        let [p, q] = self.degrees;
        let local: [[Rational; 2]; 2] = std::array::from_fn(|a| {
            std::array::from_fn(|b| {
                (&domain[a][b] - &span.domain[a][0]) / (&span.domain[a][1] - &span.domain[a][0])
            })
        });
        let mut rows = Vec::with_capacity(span.net.len());
        for row in span.net.chunks(q + 1) {
            rows.extend(restrict(row, &local[1][0], &local[1][1]));
        }
        let mut result = rows.clone();
        for j in 0..=q {
            let column: Vec<_> = (0..=p).map(|i| rows[i * (q + 1) + j].clone()).collect();
            for (i, value) in restrict(&column, &local[0][0], &local[0][1])
                .into_iter()
                .enumerate()
            {
                result[i * (q + 1) + j] = value;
            }
        }
        Ok(result)
    }
    pub(super) fn point(&self, at: &[Ball; 2], dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
        if at
            .iter()
            .zip(&self.domain)
            .any(|(a, d)| a.lo() < d[0] || a.hi() > d[1])
        {
            return Err(failure(
                "positional parameter enclosure leaves the complete source domain",
            ));
        }
        let key: [Rational; 4] = std::array::from_fn(|i| {
            if i % 2 == 0 {
                at[i / 2].lo()
            } else {
                at[i / 2].hi()
            }
        });
        dec.charge_work(1)?;
        if let Some(point) = self.points.borrow().get(&key) {
            return Ok(point.clone());
        }
        let mut result: Option<[Ball; 3]> = None;
        for span in &self.spans {
            let domain: [[Rational; 2]; 2] = std::array::from_fn(|a| {
                [
                    at[a].lo().max(span.domain[a][0].clone()),
                    at[a].hi().min(span.domain[a][1].clone()),
                ]
            });
            if domain.iter().any(|a| a[0] > a[1]) {
                continue;
            }
            let [p, q] = self.degrees;
            let blends = ((p + 1) * q * (q + 1) + p * (p + 1)) / 2;
            dec.charge_work((blends as u64).saturating_mul(64).saturating_add(128))?;
            let local: [Ball; 2] = std::array::from_fn(|a| {
                let lo = (&domain[a][0] - &span.domain[a][0])
                    / (&span.domain[a][1] - &span.domain[a][0]);
                let hi = (&domain[a][1] - &span.domain[a][0])
                    / (&span.domain[a][1] - &span.domain[a][0]);
                ball(&lo).hull(&ball(&hi))
            });
            let evaluate = |net: &[[Ball; 4]], t: &Ball| {
                let mut work = net.to_vec();
                while work.len() > 1 {
                    work = work.windows(2).map(|p| blend(&p[0], &p[1], t)).collect();
                }
                work.remove(0)
            };
            let rows: Vec<_> = span
                .net
                .chunks(q + 1)
                .map(|row| evaluate(row, &local[1]))
                .collect();
            let homogeneous = evaluate(&rows, &local[0]);
            let weight = homogeneous[3].lo();
            if weight <= Rational::zero() {
                return Err(failure(format!(
                    "positional weight lower bound {weight} is not positive"
                )));
            }
            let point: [Ball; 3] = std::array::from_fn(|i| {
                homogeneous[i]
                    .checked_div(&homogeneous[3])
                    .expect("positive weight")
            });
            result = Some(match result {
                None => point,
                Some(old) => std::array::from_fn(|i| old[i].hull(&point[i])),
            });
        }
        let result =
            result.ok_or_else(|| failure("positional point leaves the source knot cover"))?;
        self.points.borrow_mut().insert(key, result.clone());
        Ok(result)
    }
    /// Bounds derivative residuals against a constant affine reference over
    /// every intersected source span. Translation and affine subtraction take
    /// place before the positive quotient, preserving shared weight factors.
    pub(super) fn derivative_residual(
        &self,
        domain: &[[Rational; 2]; 2],
        origin: &[Rational; 3],
        reference: &[[Rational; 3]; 2],
        dec: &mut Decider,
    ) -> Result<[[Ball; 3]; 2], Refusal> {
        let [p, q] = self.degrees;
        let mut result: Option<[[Ball; 3]; 2]> = None;
        for span in &self.spans {
            let mut part: [[Rational; 2]; 2] = std::array::from_fn(|a| {
                [
                    domain[a][0].clone().max(span.domain[a][0].clone()),
                    domain[a][1].clone().min(span.domain[a][1].clone()),
                ]
            });
            if part.iter().any(|a| a[0] >= a[1]) {
                continue;
            }
            for a in 0..2 {
                if (&part[a][1] - &part[a][0]).mul_pow2(8) < &span.domain[a][1] - &span.domain[a][0]
                {
                    part[a] = span.domain[a].clone();
                }
            }
            let mut net = self.net_on(span, &part, dec)?;
            dec.charge_work((net.len() as u64).saturating_mul(768))?;
            for h in &mut net {
                for i in 0..3 {
                    h[i] = &h[i] - &(&h[3] * &ball(&origin[i]));
                }
            }
            let weight = net
                .iter()
                .map(|h| h[3].clone())
                .reduce(|a, b| a.hull(&b))
                .unwrap();
            if weight.lo() <= Rational::zero() {
                return Err(failure(format!(
                    "evaluator weight lower bound {} is not positive",
                    weight.lo()
                )));
            }
            let position: [Ball; 3] = std::array::from_fn(|i| {
                net.iter()
                    .map(|h| h[i].checked_div(&h[3]).expect("positive weights"))
                    .reduce(|a, b| a.hull(&b))
                    .unwrap()
            });
            let mut residual: [[Ball; 3]; 2] =
                std::array::from_fn(|_| std::array::from_fn(|_| ball(&Rational::zero())));
            for a in 0..2 {
                let degree = [p, q][a];
                let factor = ball(
                    &(Rational::from_i64(degree as i64)
                        * ((&domain[a][1] - &domain[a][0]) / (&part[a][1] - &part[a][0]))),
                );
                let delta = |i: usize, j: usize| -> [Ball; 4] {
                    let first = i * (q + 1) + j;
                    let next = first + if a == 0 { q + 1 } else { 1 };
                    std::array::from_fn(|k| &(&net[next][k] - &net[first][k]) * &factor)
                };
                let mut derivative_weight: Option<Ball> = None;
                let mut numerator: Option<[Ball; 3]> = None;
                for i in 0..=p {
                    for j in 0..=q {
                        let n = [i, j][a];
                        let mut elevated: [Ball; 4] =
                            std::array::from_fn(|_| ball(&Rational::zero()));
                        if n > 0 {
                            let before = if a == 0 {
                                delta(i - 1, j)
                            } else {
                                delta(i, j - 1)
                            };
                            let fraction =
                                ball(&Rational::from_ratio(n as i64, degree as i64).unwrap());
                            for k in 0..4 {
                                elevated[k] = &elevated[k] + &(&before[k] * &fraction);
                            }
                        }
                        if n < degree {
                            let after = delta(i, j);
                            let fraction = ball(
                                &Rational::from_ratio((degree - n) as i64, degree as i64).unwrap(),
                            );
                            for k in 0..4 {
                                elevated[k] = &elevated[k] + &(&after[k] * &fraction);
                            }
                        }
                        derivative_weight = Some(
                            derivative_weight
                                .map_or(elevated[3].clone(), |old| old.hull(&elevated[3])),
                        );
                        let value: [Ball; 3] = std::array::from_fn(|k| {
                            &elevated[k] - &(&ball(&reference[a][k]) * &net[i * (q + 1) + j][3])
                        });
                        numerator = Some(match numerator {
                            None => value,
                            Some(old) => std::array::from_fn(|k| old[k].hull(&value[k])),
                        });
                    }
                }
                let numerator = numerator.unwrap();
                let dw = derivative_weight.unwrap();
                residual[a] = std::array::from_fn(|k| {
                    (&numerator[k] - &(&position[k] * &dw))
                        .checked_div(&weight)
                        .expect("positive weights")
                });
            }
            result = Some(match result {
                None => residual,
                Some(old) => std::array::from_fn(|a| {
                    std::array::from_fn(|k| old[a][k].hull(&residual[a][k]))
                }),
            });
        }
        result.ok_or_else(|| failure("evaluator derivative cover is empty"))
    }

    fn cell(
        &self,
        parameters: &[[Ball; 2]; 4],
        dec: &mut Decider,
    ) -> Result<(Rational, [Rational; 2]), Refusal> {
        let domain: [[Rational; 2]; 2] = std::array::from_fn(|a| {
            [
                parameters.iter().map(|p| p[a].lo()).min().unwrap(),
                parameters.iter().map(|p| p[a].hi()).max().unwrap(),
            ]
        });
        if domain.iter().any(|a| a[0] >= a[1]) {
            return Err(failure("evaluator cell has zero parameter width"));
        }
        let mut points = Vec::new();
        for at in parameters {
            points.push(self.point(at, dec)?);
        }
        let points: [[Ball; 3]; 4] = points.try_into().unwrap();
        dec.charge_work(4096)?;
        let gradients: [[Ball; 3]; 2] = [
            std::array::from_fn(|k| {
                (&points[2][k] - &points[0][k]).hull(&(&points[3][k] - &points[1][k]))
            }),
            std::array::from_fn(|k| {
                (&points[3][k] - &points[2][k]).hull(&(&points[1][k] - &points[0][k]))
            }),
        ];
        let normalized: [[Ball; 2]; 4] = std::array::from_fn(|i| {
            std::array::from_fn(|a| {
                (&parameters[i][a] - &ball(&domain[a][0]))
                    .checked_div(&ball(&(&domain[a][1] - &domain[a][0])))
                    .expect("positive cell widths")
            })
        });
        let xr =
            (&normalized[2][0] - &normalized[0][0]).hull(&(&normalized[3][0] - &normalized[1][0]));
        let xt =
            (&normalized[1][0] - &normalized[0][0]).hull(&(&normalized[3][0] - &normalized[2][0]));
        let yt =
            (&normalized[1][1] - &normalized[0][1]).hull(&(&normalized[3][1] - &normalized[2][1]));
        if xr.lo() <= Rational::zero() || yt.lo() <= Rational::zero() {
            return Err(failure("evaluator cut map needs positive axis speeds"));
        }
        let reference: [[Rational; 3]; 2] = [
            std::array::from_fn(|k| gradients[0][k].mid() / xr.mid()),
            std::array::from_fn(|k| {
                (gradients[1][k].mid() - (gradients[0][k].mid() / xr.mid()) * xt.mid()) / yt.mid()
            }),
        ];
        let origin: [Rational; 3] = std::array::from_fn(|k| {
            points
                .iter()
                .fold(Rational::zero(), |s, p| s + p[k].mid())
                .mul_pow2(-2)
        });
        let raw = self.derivative_residual(&domain, &origin, &reference, dec)?;
        let first: [[Ball; 3]; 2] = [
            std::array::from_fn(|k| {
                &(&(&raw[0][k] * &xr) + &(&ball(&reference[0][k]) * &xr)) - &gradients[0][k]
            }),
            std::array::from_fn(|k| {
                &(&(&(&raw[0][k] * &xt) + &(&raw[1][k] * &yt))
                    + &(&(&ball(&reference[0][k]) * &xt) + &(&ball(&reference[1][k]) * &yt)))
                    - &gradients[1][k]
            }),
        ];
        // A 3x3 evaluator stencil includes the common triangle diagonal.
        // Every parameter point is within 1/4 of a stencil node on each axis.
        let contributions = first.each_ref().map(|v| norm(v).mul_pow2(-2));
        let mut measured = Rational::zero();
        for i in 0..=2 {
            for j in 0..=2 {
                let r = ball(&Rational::from_ratio(i, 2).unwrap());
                let t = ball(&Rational::from_ratio(j, 2).unwrap());
                let at = blend(
                    &blend(&parameters[0], &parameters[2], &r),
                    &blend(&parameters[1], &parameters[3], &r),
                    &t,
                );
                let value = self.point(&at, dec)?;
                let facet: [Ball; 3] = if i >= j {
                    std::array::from_fn(|k| {
                        &(&points[0][k] + &(&(&points[2][k] - &points[0][k]) * &r))
                            + &(&(&points[3][k] - &points[2][k]) * &t)
                    })
                } else {
                    std::array::from_fn(|k| {
                        &(&points[0][k] + &(&(&points[3][k] - &points[1][k]) * &r))
                            + &(&(&points[1][k] - &points[0][k]) * &t)
                    })
                };
                measured = measured.max(norm(&std::array::from_fn(|k| &value[k] - &facet[k])));
            }
        }
        Ok((
            &(&measured + &contributions[0]) + &contributions[1],
            contributions,
        ))
    }
    pub(super) fn grid(
        &self,
        mut axes: [Vec<Real>; 2],
        required: [&[Real]; 2],
        controls: Controls<'_>,
        map: impl Fn(&[Ball; 2]) -> Result<[Ball; 2], Refusal>,
        dec: &mut Decider,
    ) -> Result<super::rational_patch::Grid, Refusal> {
        let Controls {
            tolerance,
            max_vertices,
            displacement,
        } = controls;
        for a in 0..2 {
            for value in required[a] {
                super::rational_patch::insert(&mut axes[a], value.clone(), dec)?;
            }
        }
        let mut accepted = BTreeMap::<[Rational; 8], (Rational, [Rational; 2])>::new();
        let mut last_failed = Rational::zero();
        loop {
            dec.check_deadline()?;
            let count = axes[0]
                .len()
                .checked_mul(axes[1].len())
                .filter(|n| *n <= max_vertices)
                .ok_or_else(|| {
                    Refusal::new("mesh.display-bound-exceeded",format!("complete-cell deviation bound {last_failed} mm cannot meet {tolerance} mm within vertex cap {max_vertices}"))
                })?;
            let boxes: [Vec<Ball>; 2] = [0, 1]
                .map(|a| {
                    axes[a]
                        .iter()
                        .map(|t| dec.enclose_at(t, BITS))
                        .collect::<Result<_, _>>()
                })
                .into_iter()
                .collect::<Result<Vec<_>, _>>()?
                .try_into()
                .unwrap();
            let mut split_cell = None;
            let mut worst = Rational::zero();
            'cells: for i in 0..axes[0].len() - 1 {
                for j in 0..axes[1].len() - 1 {
                    dec.charge_work(1)?;
                    let key = std::array::from_fn(|k| {
                        let a = k / 4;
                        let n = [i, j][a] + (k % 4) / 2;
                        if k % 2 == 0 {
                            boxes[a][n].lo()
                        } else {
                            boxes[a][n].hi()
                        }
                    });
                    let (chord, contributions) = if let Some(chord) = accepted.get(&key) {
                        chord.clone()
                    } else {
                        let mut parameters = Vec::new();
                        for a in 0..2 {
                            for b in 0..2 {
                                parameters.push(map(&[
                                    boxes[0][i + a].clone(),
                                    boxes[1][j + b].clone(),
                                ])?);
                            }
                        }
                        let chord = self.cell(&parameters.try_into().unwrap(), dec)?;
                        accepted.insert(key, chord.clone());
                        chord
                    };
                    if &chord + displacement > tolerance * &Rational::from_ratio(7, 8).unwrap() {
                        last_failed = &chord + displacement;
                        split_cell =
                            Some((usize::from(contributions[1] > contributions[0]), [i, j]));
                        break 'cells;
                    }
                    worst = worst.max(chord);
                }
            }
            if let Some((a, cell)) = split_cell {
                let i = cell[a];
                let midpoint = &(&axes[a][i] + &axes[a][i + 1]) / &Real::from_i64(2);
                axes[a].insert(i + 1, midpoint);
                continue;
            }
            let mut points = Vec::with_capacity(count);
            let mut uncertainty = Rational::zero();
            for u in &boxes[0] {
                for v in &boxes[1] {
                    let point = self.point(&map(&[u.clone(), v.clone()])?, dec)?;
                    uncertainty = uncertainty.max(
                        point
                            .iter()
                            .fold(Rational::zero(), |s, x| s + x.rad().mul_pow2(1)),
                    );
                    points.push(point);
                }
            }
            let chord = &worst + displacement + uncertainty;
            if chord > *tolerance {
                return Err(failure(format!(
                    "positional bound {chord} mm exceeds {tolerance} mm"
                )));
            }
            return Ok(super::rational_patch::Grid {
                grid: EnvelopeGrid {
                    columns: axes[0].len(),
                    rows: axes[1].len(),
                    points,
                    chord,
                },
                parameters: axes,
                normal_boxes: None,
                embedding: super::rational_patch::Embedding::PositionOnly,
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::Point3;
    fn saddle() -> PositionSurface {
        let mut dec = Decider::default();
        let controls = (0..2)
            .flat_map(|u| {
                (0..2).map(move |v| {
                    Point3::new(Real::from_i64(u), Real::from_i64(v), Real::from_i64(u * v))
                })
            })
            .collect();
        let knots = [0, 0, 1, 1].map(Rational::from_i64).to_vec();
        let surface = BSplineSurface::new(
            [1, 1],
            [knots.clone(), knots],
            [2, 2],
            controls,
            vec![Real::from_i64(1); 4],
            &mut dec,
        )
        .unwrap();
        PositionSurface::new(&surface, &dec).unwrap()
    }
    #[test]
    fn bilinear_saddle_bound_covers_the_interior_and_refuses_above_level() {
        let surface = saddle();
        let parameters: [[Ball; 2]; 4] = std::array::from_fn(|i| {
            [
                ball(&Rational::from_i64((i / 2) as i64)),
                ball(&Rational::from_i64((i % 2) as i64)),
            ]
        });
        let measured = surface
            .cell(&parameters, &mut Decider::default())
            .unwrap()
            .0;
        // At the diagonal midpoint the surface has z=1/4 and the facet z=1/2.
        assert!(measured >= Rational::from_ratio(1, 4).unwrap());
        assert!(measured <= Rational::from_i64(1));
        let axes = [
            vec![Real::from_i64(0), Real::from_i64(1)],
            vec![Real::from_i64(0), Real::from_i64(1)],
        ];
        let error = match surface.grid(
            axes,
            [&[], &[]],
            Controls {
                tolerance: &Rational::from_ratio(1, 100).unwrap(),
                max_vertices: 4,
                displacement: &Rational::zero(),
            },
            |at| Ok(at.clone()),
            &mut Decider::default(),
        ) {
            Err(error) => error,
            Ok(_) => panic!("a 0.25 mm deviation cannot receive the 0.01 mm rule"),
        };
        assert_eq!(error.code, "mesh.display-bound-exceeded");
        assert!(error.detail.contains("deviation bound"), "{error:?}");
        assert!(error.detail.contains("1/100"), "{error:?}");
    }
    #[test]
    fn derivatives_reject_a_peak_between_zero_samples() {
        let mut dec = Decider::default();
        let heights = [
            Rational::zero(),
            Rational::from_ratio(-1, 6).unwrap(),
            Rational::from_ratio(1, 6).unwrap(),
            Rational::zero(),
        ];
        let controls = (0..4)
            .flat_map(|i| {
                let height = heights[i].clone();
                (0..2).map(move |j| {
                    Point3::new(
                        Real::from(Rational::from_ratio(i as i64, 3).unwrap()),
                        Real::from_i64(j),
                        Real::from(height.clone()),
                    )
                })
            })
            .collect();
        let knots = [
            vec![Rational::zero(); 4]
                .into_iter()
                .chain(vec![Rational::one(); 4])
                .collect(),
            [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
        ];
        let carrier = BSplineSurface::new(
            [3, 1],
            knots,
            [4, 2],
            controls,
            vec![Real::from_i64(1); 8],
            &mut dec,
        )
        .unwrap();
        let surface = PositionSurface::new(&carrier, &dec).unwrap();
        let parameters = std::array::from_fn(|i| {
            [
                ball(&Rational::from_i64((i / 2) as i64)),
                ball(&Rational::from_i64((i % 2) as i64)),
            ]
        });
        for u in [
            Rational::zero(),
            Rational::from_ratio(1, 2).unwrap(),
            Rational::one(),
        ] {
            let value = surface
                .point(&[ball(&u), ball(&Rational::zero())], &mut dec)
                .unwrap();
            assert!(value[2].lo() <= Rational::zero() && value[2].hi() >= Rational::zero());
        }
        let measured = surface.cell(&parameters, &mut dec).unwrap().0;
        assert!(
            measured >= Rational::from_ratio(3, 64).unwrap(),
            "the unsampled quarter-point peak is covered"
        );
        assert!(measured > Rational::from_ratio(1, 100).unwrap());
    }
}
