//! Whole-domain embedding proofs, independent of display tessellation.

use super::*;
use catallax_numbers::Ball;
use catallax_recipe::LoftShape;
use catallax_shapes::{BSplineSurface, Point3, Vector3};

const MAX_CELLS: usize = 16_384;
const MAX_DEPTH: usize = 24;

pub(super) enum Proof {
    Empty,
    Solid { positive: bool, end_positive: bool },
}

fn vector(values: &[Real]) -> Vector3 {
    Vector3::new(values[0].clone(), values[1].clone(), values[2].clone())
}
fn components(v: &Vector3) -> [&Real; 3] {
    [&v.x, &v.y, &v.z]
}
fn origin(profile: &Profile) -> Point3 {
    let p = &profile.pose.translation;
    Point3::new(p.x.clone(), p.y.clone(), p.z.clone())
}
fn frame(profile: &Profile) -> [Vector3; 3] {
    [
        Vector3::new(Real::one(), Real::zero(), Real::zero()),
        Vector3::new(Real::zero(), Real::one(), Real::zero()),
        Vector3::new(Real::zero(), Real::zero(), Real::one()),
    ]
    .map(|v| profile.pose.rotation.apply(&v))
}
fn pole(h: &[Real; 4]) -> Point3 {
    Point3::new(&h[0] / &h[3], &h[1] / &h[3], &h[2] / &h[3])
}
fn local(h: &[Real; 4], origin: &Point3, frame: &[Vector3; 3]) -> [Real; 2] {
    let v = Vector3::new(
        &h[0] - &(&origin.x * &h[3]),
        &h[1] - &(&origin.y * &h[3]),
        &h[2] - &(&origin.z * &h[3]),
    );
    [v.dot(&frame[0]), v.dot(&frame[1])]
}
fn cross(a: &[Real; 2], b: &[Real; 2]) -> Real {
    &a[0] * &b[1] - &a[1] * &b[0]
}

// Fit from controls only proposes the map. Every homogeneous coefficient of
// every loop must then satisfy that map exactly; samples never admit a body.
fn affine(
    loops: &[Vec<math::Interval>],
    o: &Point3,
    f: &[Vector3; 3],
    dec: &mut Decider,
) -> Result<Option<Vec<[Real; 9]>>, Refusal> {
    let all: Vec<_> = loops
        .iter()
        .flat_map(|parts| parts.iter())
        .flat_map(|part| (0..part.sections[0].len()).map(move |k| (part, k)))
        .collect();
    let (first, k) = all[0];
    let a_h = &first.sections[0][k];
    let a = local(a_h, o, f).map(|x| &x / &a_h[3]);
    let mut basis = None;
    'search: for (bpart, bindex) in &all {
        let bh = &bpart.sections[0][*bindex];
        let b = local(bh, o, f).map(|x| &x / &bh[3]);
        let ab = [&b[0] - &a[0], &b[1] - &a[1]];
        if sign(&(&ab[0].square() + &ab[1].square()), dec)? == Sign::Zero {
            continue;
        }
        for (cpart, cindex) in &all {
            let ch = &cpart.sections[0][*cindex];
            let c = local(ch, o, f).map(|x| &x / &ch[3]);
            let ac = [&c[0] - &a[0], &c[1] - &a[1]];
            let determinant = cross(&ab, &ac);
            if sign(&determinant, dec)? != Sign::Zero {
                basis = Some((*bpart, *bindex, *cpart, *cindex, ab, ac, determinant));
                break 'search;
            }
        }
    }
    let Some((bpart, bindex, cpart, cindex, ab, ac, determinant)) = basis else {
        return Err(unproven(
            "the profile has no independent affine reference controls",
        ));
    };
    let mut maps = Vec::new();
    for station in 0..first.sections.len() {
        let pa = pole(&first.sections[station][k]);
        let db = pole(&bpart.sections[station][bindex]).minus(&pa);
        let dc = pole(&cpart.sections[station][cindex]).minus(&pa);
        let x = db
            .scaled(&ac[1])
            .minus(&dc.scaled(&ab[1]))
            .scaled(&(&Real::one() / &determinant));
        let y = dc
            .scaled(&ab[0])
            .minus(&db.scaled(&ac[0]))
            .scaled(&(&Real::one() / &determinant));
        let translation = Vector3::new(pa.x.clone(), pa.y.clone(), pa.z.clone())
            .minus(&x.scaled(&a[0]))
            .minus(&y.scaled(&a[1]));
        for (part, index) in &all {
            let base = &part.sections[0][*index];
            let target = &part.sections[station][*index];
            let xy = local(base, o, f);
            let expected = x
                .scaled(&xy[0])
                .plus(&y.scaled(&xy[1]))
                .plus(&translation.scaled(&base[3]));
            for (axis, value) in components(&expected).iter().enumerate() {
                if sign(&(&target[axis] - *value), dec)? != Sign::Zero {
                    return Ok(None);
                }
            }
        }
        maps.push([
            x.x,
            x.y,
            x.z,
            y.x,
            y.y,
            y.z,
            translation.x,
            translation.y,
            translation.z,
        ]);
    }
    Ok(Some(maps))
}

fn split<const N: usize>(
    bounds: &[[Rational; 2]; N],
    widths: &[Rational; N],
) -> [[[Rational; 2]; N]; 2] {
    let mut axis = 0;
    for i in 1..N {
        if &(&bounds[i][1] - &bounds[i][0]) * &widths[axis]
            > &(&bounds[axis][1] - &bounds[axis][0]) * &widths[i]
        {
            axis = i;
        }
    }
    let middle = (&bounds[axis][0] + &bounds[axis][1]).mul_pow2(-1);
    let mut low = bounds.clone();
    let mut high = bounds.clone();
    low[axis][1] = middle.clone();
    high[axis][0] = middle;
    [low, high]
}

fn prove<const N: usize>(
    bounds: [[Rational; 2]; N],
    mut accepts: impl FnMut(&[[Rational; 2]; N], &mut Decider) -> Result<bool, Refusal>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let widths = bounds.each_ref().map(|[lo, hi]| hi - lo);
    let mut pending = vec![(bounds, 0)];
    let mut visited = 0;
    while let Some((cell, depth)) = pending.pop() {
        dec.check_deadline()?;
        visited += 1;
        if visited > MAX_CELLS {
            return Err(unproven(
                "whole-domain admissibility exceeded its finite proof cover",
            ));
        }
        if accepts(&cell, dec)? {
            continue;
        }
        if depth == MAX_DEPTH {
            return Err(unproven(
                "whole-domain separation or orientation was not proved",
            ));
        }
        let [low, high] = split(&cell, &widths);
        pending.push((high, depth + 1));
        pending.push((low, depth + 1));
    }
    Ok(())
}

fn balls<const N: usize>(bounds: &[[Rational; 2]; N], bits: u32) -> [Ball; N] {
    bounds
        .each_ref()
        .map(|[lo, hi]| Ball::from_rational(lo, bits).hull(&Ball::from_rational(hi, bits)))
}

fn polynomial(values: &[Real], variable: usize, nvars: usize, offset: usize) -> RealPoly {
    let n = values.len() - 1;
    let mut row = values.to_vec();
    let mut choose = Rational::one();
    let mut terms = Vec::with_capacity(values.len());
    for order in 0..=n {
        let mut powers = vec![0; nvars];
        powers[variable] = order as u32;
        terms.push((powers, &row[0] * &Real::from_rational(choose.clone())));
        for i in 0..n - order {
            row[i] = &row[i + 1] - &row[i];
        }
        if order < n {
            choose = &choose
                * &Rational::from_ratio((n - order) as i64, (order + 1) as i64)
                    .expect("positive degree");
        }
    }
    let poly = RealPoly::from_terms(nvars, terms);
    poly.substitute(
        variable,
        &(RealPoly::var(nvars, variable)
            - RealPoly::constant(nvars, Real::from_i64(offset as i64))),
    )
}

fn positive(
    poly: &RealPoly,
    bounds: [[Rational; 2]; 2],
    bits: u32,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    prove(
        bounds,
        |cell, dec| Ok(poly.eval_ball(&balls(cell, bits), dec)?.lo() > Rational::zero()),
        dec,
    )
}

fn scalar_positive(values: &[Real], bits: u32, dec: &mut Decider) -> Result<(), Refusal> {
    let poly = polynomial(values, 0, 1, 0);
    prove(
        [[Rational::zero(), Rational::one()]],
        |cell, dec| Ok(poly.eval_ball(&balls(cell, bits), dec)?.lo() > Rational::zero()),
        dec,
    )
}

fn bounds(source: &BSplineSurface, span: [usize; 2]) -> [[Rational; 2]; 2] {
    let knots = source.knots();
    std::array::from_fn(|axis| {
        [
            knots[axis][span[axis]].clone(),
            knots[axis][span[axis] + 1].clone(),
        ]
    })
}

fn hull(
    source: &BSplineSurface,
    span: [usize; 2],
    domain: &[[Rational; 2]; 2],
    bits: u32,
    dec: &mut Decider,
) -> Result<[[Rational; 2]; 3], Refusal> {
    let net = source.span_bezier_controls_on(span, domain.each_ref().map(|x| x.each_ref()), dec)?;
    let mut result: Option<[[Rational; 2]; 3]> = None;
    for h in net {
        let p: [Ball; 3] = [
            dec.enclose_budgeted_at(&(&h[0] / &h[3]), bits)?,
            dec.enclose_budgeted_at(&(&h[1] / &h[3]), bits)?,
            dec.enclose_budgeted_at(&(&h[2] / &h[3]), bits)?,
        ];
        if let Some(result) = &mut result {
            for axis in 0..3 {
                result[axis][0] = result[axis][0].clone().min(p[axis].lo());
                result[axis][1] = result[axis][1].clone().max(p[axis].hi());
            }
        } else {
            result = Some(p.map(|x| [x.lo(), x.hi()]));
        }
    }
    Ok(result.expect("positive Bezier net"))
}

fn separate(sources: &[BSplineSurface], bits: u32, dec: &mut Decider) -> Result<(), Refusal> {
    for (index, a) in sources.iter().enumerate() {
        for b in &sources[index + 1..] {
            for sa in a.span_indices() {
                let da = bounds(a, sa);
                for sb in b.span_indices() {
                    let db = bounds(b, sb);
                    if da[1] != db[1] {
                        continue;
                    }
                    prove(
                        [da[0].clone(), db[0].clone(), da[1].clone()],
                        |cell, dec| {
                            let ha = hull(a, sa, &[cell[0].clone(), cell[2].clone()], bits, dec)?;
                            let hb = hull(b, sb, &[cell[1].clone(), cell[2].clone()], bits, dec)?;
                            Ok((0..3)
                                .any(|axis| ha[axis][1] < hb[axis][0] || hb[axis][1] < ha[axis][0]))
                        },
                        dec,
                    )?;
                }
            }
        }
    }
    Ok(())
}

fn polar(
    loops: &[Vec<math::Interval>],
    sources: &[BSplineSurface],
    o: &Point3,
    f: &[Vector3; 3],
    shape: LoftShape,
    bits: u32,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for (parts, source) in loops.iter().zip(sources) {
        let centers: Vec<[Real; 2]> = (0..parts[0].sections.len())
            .map(|station| {
                let mut sum = [Real::zero(), Real::zero()];
                for part in parts {
                    let p = &part.sections[station][0];
                    let xy = local(p, o, f);
                    for axis in 0..2 {
                        sum[axis] = &sum[axis] + &(&xy[axis] / &p[3]);
                    }
                }
                sum.map(|x| x / Real::from_i64(parts.len() as i64))
            })
            .collect();
        let centers = math::stations(&centers, shape);
        let mut orientation = None;
        for span in source.span_indices() {
            let domain = bounds(source, span);
            let band = centers
                .iter()
                .enumerate()
                .find(|(i, _)| domain[1][0] == Rational::from_i64(*i as i64))
                .expect("integer loft station");
            let h = source.span_polynomials(span, dec)?;
            let xy: [RealPoly; 2] = std::array::from_fn(|axis| {
                let direction = components(&f[axis]);
                let origin = [&o.x, &o.y, &o.z];
                (0..3).fold(RealPoly::zero(2), |sum, i| {
                    sum + (&h[i] - &h[3].scale(origin[i])).scale(direction[i])
                })
            });
            let center: [RealPoly; 2] = std::array::from_fn(|axis| {
                polynomial(
                    &band.1.iter().map(|p| p[axis].clone()).collect::<Vec<_>>(),
                    1,
                    2,
                    band.0,
                )
            });
            let radial = [
                &xy[0] - &(&center[0] * &h[3]),
                &xy[1] - &(&center[1] * &h[3]),
            ];
            let tangent = [
                &xy[0].partial(0) * &h[3] - &xy[0] * &h[3].partial(0),
                &xy[1].partial(0) * &h[3] - &xy[1] * &h[3].partial(0),
            ];
            let mut jacobian = &radial[0] * &tangent[1] - &radial[1] * &tangent[0];
            let expected = match orientation {
                Some(value) => value,
                None => {
                    let at = domain
                        .each_ref()
                        .map(|[lo, hi]| Real::from_rational((lo + hi).mul_pow2(-1)));
                    let value = sign(&jacobian.eval_real(&at), dec)?;
                    if value == Sign::Zero {
                        return Err(unproven(
                            "the proposed polar center has no strict angular Jacobian",
                        ));
                    }
                    orientation = Some(value);
                    value
                }
            };
            if expected == Sign::Negative {
                jacobian = -jacobian;
            }
            positive(&jacobian, domain, bits, dec)?;
        }
    }
    separate(sources, bits, dec)
}

fn spatial_affine(
    maps: &[[Real; 9]],
    loops: &[Vec<math::Interval>],
    o: &Point3,
    f: &[Vector3; 3],
    shape: LoftShape,
    bits: u32,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let mut xy: Option<[[Rational; 2]; 2]> = None;
    for part in loops.iter().flatten() {
        for h in &part.sections[0] {
            let p = local(h, o, f);
            let values = [
                dec.enclose_budgeted_at(&(&p[0] / &h[3]), bits)?,
                dec.enclose_budgeted_at(&(&p[1] / &h[3]), bits)?,
            ];
            if let Some(xy) = &mut xy {
                for axis in 0..2 {
                    xy[axis][0] = xy[axis][0].clone().min(values[axis].lo());
                    xy[axis][1] = xy[axis][1].clone().max(values[axis].hi());
                }
            } else {
                xy = Some(values.map(|x| [x.lo(), x.hi()]));
            }
        }
    }
    let xy = xy.expect("closed profile");
    let widths = [
        &xy[0][1] - &xy[0][0],
        &xy[1][1] - &xy[1][0],
        Rational::from_i64((maps.len() - 1) as i64),
    ];
    let bands = math::stations(maps, shape);
    let polynomials: Vec<[RealPoly; 3]> = bands
        .iter()
        .enumerate()
        .map(|(band, controls)| {
            std::array::from_fn(|axis| {
                let x = polynomial(
                    &controls.iter().map(|c| c[axis].clone()).collect::<Vec<_>>(),
                    2,
                    3,
                    band,
                );
                let y = polynomial(
                    &controls
                        .iter()
                        .map(|c| c[axis + 3].clone())
                        .collect::<Vec<_>>(),
                    2,
                    3,
                    band,
                );
                let b = polynomial(
                    &controls
                        .iter()
                        .map(|c| c[axis + 6].clone())
                        .collect::<Vec<_>>(),
                    2,
                    3,
                    band,
                );
                &x * &RealPoly::var(3, 0) + &y * &RealPoly::var(3, 1) + b
            })
        })
        .collect();
    let middle_band = (bands.len() - 1) / 2;
    let at = [
        Real::from_rational((&xy[0][0] + &xy[0][1]).mul_pow2(-1)),
        Real::from_rational((&xy[1][0] + &xy[1][1]).mul_pow2(-1)),
        Real::from_rational(
            Rational::from_i64(middle_band as i64) + Rational::from_ratio(1, 2).expect("nonzero"),
        ),
    ];
    let tangents: [Vector3; 3] = std::array::from_fn(|axis| {
        Vector3::new(
            polynomials[middle_band][0].partial(axis).eval_real(&at),
            polynomials[middle_band][1].partial(axis).eval_real(&at),
            polynomials[middle_band][2].partial(axis).eval_real(&at),
        )
        .scaled(&Real::from_rational(widths[axis].clone()))
    });
    let determinant = tangents[0].dot(&tangents[1].cross(&tangents[2]));
    let orientation = sign(&determinant, dec)?;
    if orientation == Sign::Zero {
        return Err(unproven(
            "the affine volume map has no independent reference Jacobian",
        ));
    }
    let projection = [
        tangents[1].cross(&tangents[2]),
        tangents[2].cross(&tangents[0]),
        tangents[0].cross(&tangents[1]),
    ]
    .map(|v| v.scaled(&(&Real::one() / &determinant)));
    for (band, coordinates) in polynomials.iter().enumerate() {
        let projected: [[RealPoly; 3]; 3] = std::array::from_fn(|row| {
            std::array::from_fn(|column| {
                let p = components(&projection[row]);
                (0..3)
                    .fold(RealPoly::zero(3), |sum, axis| {
                        sum + coordinates[axis].partial(column).scale(p[axis])
                    })
                    .scale(&Real::from_rational(widths[column].clone()))
            })
        });
        let symmetric: [[RealPoly; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (&projected[i][j] + &projected[j][i]).scale(&Real::from_rational(
                    Rational::from_ratio(1, 2).expect("nonzero"),
                ))
            })
        });
        prove(
            [
                xy[0].clone(),
                xy[1].clone(),
                [
                    Rational::from_i64(band as i64),
                    Rational::from_i64((band + 1) as i64),
                ],
            ],
            |cell, dec| {
                let at = balls(cell, bits);
                for (i, row) in symmetric.iter().enumerate() {
                    let mut lower = row[i].eval_ball(&at, dec)?.lo();
                    for (j, entry) in row.iter().enumerate() {
                        if i == j {
                            continue;
                        }
                        let value = entry.eval_ball(&at, dec)?;
                        lower = lower - value.lo().abs().max(value.hi().abs());
                    }
                    if lower <= Rational::zero() {
                        return Ok(false);
                    }
                }
                Ok(true)
            },
            dec,
        )?;
    }
    Ok(orientation == Sign::Positive)
}

pub(super) fn certify(
    loops: &[Vec<math::Interval>],
    sources: &[BSplineSurface],
    profiles: &[&Profile],
    shape: LoftShape,
    dec: &mut Decider,
) -> Result<Proof, Refusal> {
    let f = frame(profiles[0]);
    let o = origin(profiles[0]);
    let last_frame = frame(profiles.last().expect("loft sections"));
    let bits = dec.policy().max_precision_bits.min(128);
    if bits == 0 {
        return Err(unproven(
            "no arithmetic precision is available for loft admissibility",
        ));
    }
    let mut parallel = true;
    for profile in profiles {
        parallel &= sign(&frame(profile)[2].cross(&f[2]).norm_squared(), dec)? == Sign::Zero;
    }
    let heights: Vec<_> = profiles
        .iter()
        .map(|p| [origin(p).minus(&o).dot(&f[2])])
        .collect();
    if parallel {
        let mut flat = true;
        for h in &heights {
            flat &= sign(&h[0], dec)? == Sign::Zero;
        }
        if flat {
            return Ok(Proof::Empty);
        }
        let end = &heights.last().expect("loft sections")[0];
        let direction = sign(end, dec)?;
        if direction == Sign::Zero {
            return Err(unproven(
                "the station planes return to their starting plane",
            ));
        }
        for controls in math::stations(&heights, shape) {
            let degree = Real::from_i64((controls.len() - 1) as i64);
            let derivative: Vec<_> = controls
                .windows(2)
                .map(|p| {
                    let value = &(&p[1][0] - &p[0][0]) * &degree;
                    if direction == Sign::Positive {
                        value
                    } else {
                        -value
                    }
                })
                .collect();
            scalar_positive(&derivative, bits, dec)?;
        }
        if let Some(maps) = affine(loops, &o, &f, dec)? {
            for controls in math::stations(&maps, shape) {
                let entries: [Vec<Real>; 4] = std::array::from_fn(|index| {
                    controls
                        .iter()
                        .map(|c| {
                            vector(&c[(index / 2) * 3..(index / 2 + 1) * 3]).dot(&f[index % 2])
                        })
                        .collect()
                });
                let a = math::product(&entries[0], &entries[3]);
                let b = math::product(&entries[1], &entries[2]);
                scalar_positive(
                    &a.into_iter().zip(b).map(|(a, b)| a - b).collect::<Vec<_>>(),
                    bits,
                    dec,
                )?;
            }
        } else {
            polar(loops, sources, &o, &f, shape, bits, dec)?;
        }
        let positive = direction == Sign::Positive;
        let end_positive = positive == (sign(&last_frame[2].dot(&f[2]), dec)? == Sign::Positive);
        return Ok(Proof::Solid {
            positive,
            end_positive,
        });
    }
    let maps = affine(loops, &o, &f, dec)?.ok_or_else(|| {
        unproven("nonparallel nonaffine sections need a complete spatial embedding proof")
    })?;
    let positive = spatial_affine(&maps, loops, &o, &f, shape, bits, dec)?;
    let last = maps.last().expect("loft sections");
    let normal = vector(&last[..3]).cross(&vector(&last[3..6]));
    let end_positive = positive == (sign(&normal.dot(&last_frame[2]), dec)? == Sign::Positive);
    Ok(Proof::Solid {
        positive,
        end_positive,
    })
}

/// A pointwise translated family on one straight, strictly traversed segment
/// has exactly the existing prism image. Preserve the authored correspondence
/// when proving the reduction; merged profile counts never select it.
pub(super) fn prism(
    loops: &[Vec<math::Interval>],
    profile: &Profile,
    shape: LoftShape,
    dec: &mut Decider,
) -> Result<Option<catallax_shapes::RVector3>, Refusal> {
    // Keep spline prisms on the rational atlas route, including its exchange
    // cover; the older polynomial extrusion carrier has narrower STEP coverage.
    if profile.has_bezier() {
        return Ok(None);
    }
    let first = &loops[0][0].sections;
    let start = pole(&first[0][0]);
    let deltas: Vec<_> = first
        .iter()
        .map(|section| pole(&section[0]).minus(&start))
        .collect();
    let travel = deltas.last().expect("loft sections");
    let length = travel.norm_squared();
    if sign(&length, dec)? == Sign::Zero {
        return Ok(None);
    }
    let normal = frame(profile)[2].clone();
    if sign(&travel.dot(&normal), dec)? == Sign::Zero {
        return Ok(None);
    }
    if profile
        .loops
        .iter()
        .flat_map(|boundary| &boundary.segments)
        .any(|segment| segment.circle.is_some())
        && sign(&travel.cross(&normal).norm_squared(), dec)? != Sign::Zero
    {
        return Ok(None);
    }
    let rational = components(travel).map(|coordinate| coordinate.as_rational());
    let [Some(x), Some(y), Some(z)] = rational else {
        return Ok(None);
    };
    let mut parameters = Vec::with_capacity(deltas.len());
    for delta in &deltas {
        if sign(&delta.cross(travel).norm_squared(), dec)? != Sign::Zero {
            return Ok(None);
        }
        parameters.push([delta.dot(travel) / length.clone()]);
    }
    for part in loops.iter().flatten() {
        for (station, delta) in deltas.iter().enumerate() {
            for (base, next) in part.sections[0].iter().zip(&part.sections[station]) {
                for (axis, component) in components(delta).iter().enumerate() {
                    let difference = &(&next[axis] - &base[axis]) - &(*component * &base[3]);
                    if sign(&difference, dec)? != Sign::Zero {
                        return Ok(None);
                    }
                }
            }
        }
    }
    let bits = dec.policy().max_precision_bits.min(128);
    if bits == 0 {
        return Err(unproven(
            "no precision is available for the prism correspondence proof",
        ));
    }
    for controls in math::stations(&parameters, shape) {
        let degree = Real::from_i64((controls.len() - 1) as i64);
        let derivative: Vec<_> = controls
            .windows(2)
            .map(|pair| &(&pair[1][0] - &pair[0][0]) * &degree)
            .collect();
        scalar_positive(&derivative, bits, dec)?;
    }
    Ok(Some(catallax_shapes::RVector3::new(x, y, z)))
}
