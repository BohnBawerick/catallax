//! Equal crossing cylinders have rational contact curves, unlike a regular
//! saddle. Put q = (u - 1)/2, s = 2q/(1+q²), c = (1-q²)/(1+q²). On one band
//! the centre is P (s, a c, (C a - b)c/S), a,b in {-1,1}.
//!
//! Both contact rays A,B are rational unit vectors perpendicular to the centre
//! tangent. The second intersection of the ray from -A in direction
//! (2-v)A+vB with the unit sphere traces the shorter circular section exactly.
//! Clearing denominators gives degree (6,2), with positive weights on each
//! half of q in [-1,1]. At either end A=B and the whole section is one point.
//!
//! A strict tube-reach check proves interior regularity, including the finite,
//! nonzero limit after dividing its Jacobian by 1-q². It also keeps each band
//! in its two signed axial half-planes, so different bands meet only at poles.
//! No sphere face or positive minimum width is invented at those poles.

use super::*;
use catallax_numbers::{Outcome, ZeroProof};

pub(super) type Rows = [[Vec<Real>; 4]; 3];

#[derive(Debug)]
pub(super) struct Exact {
    pub knots: [Rational; 3],
    // Row zero uses the wide frame; rows one and two use the narrow frame,
    // just as the fitted representation's line equations do.
    pub rows: [Rows; 2],
}

/// An exact-zero test whose undecided answer is kept apart from "nonzero".
pub(super) enum Zero {
    Yes,
    No,
    Undecided(Refusal),
}

/// Whether `value` is exactly zero. A zero proven only by a separation bound
/// is not accepted as an identity here, and neither it nor a `NotExact` sign
/// counts as "nonzero": both are undecided.
pub(super) fn zero(value: &Real, dec: &mut Decider) -> Result<Zero, Refusal> {
    Ok(match dec.sign(value) {
        Outcome::Proven(_) => Zero::No,
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => Zero::Yes,
        Outcome::ProvenCoincident(_) => Zero::Undecided(unresolved(
            "a zero proven only by a separation bound is not an exact identity here",
        )),
        Outcome::NotExact(n) => Zero::Undecided(Refusal::new(n.code, n.detail)),
        Outcome::Refused(error) => return Err(error),
    })
}

/// Whether `value` is exactly zero; an undecided answer is a refusal.
pub(super) fn equal(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    match zero(value, dec)? {
        Zero::Yes => Ok(true),
        Zero::No => Ok(false),
        Zero::Undecided(refusal) => Err(refusal),
    }
}

/// Whether every value is exactly zero, deciding them in order. One proven
/// nonzero value answers no even after an undecided one. Otherwise the first
/// undecided value's refusal is kept.
pub(super) fn all_zero(values: &[Real], dec: &mut Decider) -> Result<Zero, Refusal> {
    let mut undecided = None;
    for value in values {
        match zero(value, dec)? {
            Zero::Yes => {}
            Zero::No => return Ok(Zero::No),
            Zero::Undecided(refusal) => {
                undecided.get_or_insert(refusal);
            }
        }
    }
    Ok(undecided.map_or(Zero::Yes, Zero::Undecided))
}

fn signs(def: &RollingSaddleDefinition) -> [Real; 2] {
    let SaddleRing::Crossing { narrow_positive } = def.ring else {
        unreachable!("a collapsed construction keeps its crossing law")
    };
    [
        Real::from_i64(def.branch),
        Real::from_i64(if narrow_positive { 1 } else { -1 }),
    ]
}

pub(super) fn through(
    mut def: RollingSaddleDefinition,
    point: &Point3,
    dec: &mut Decider,
) -> Result<RollingSaddleDefinition, Refusal> {
    math::positive(&def.radius, "shapes.non-positive-radius", dec)?;
    math::validate_point(point, dec)?;
    for which in 0..2 {
        if !equal(&(&def.axis_gap(which, point) - &def.wide.square()), dec)? {
            return Err(invalid(
                "the crossing witness is not on both source cylinders",
            ));
        }
    }
    let [_, wide, _] = def.local(point);
    let narrow = point.minus(&def.foot).dot(&def.side_axis());
    let [a, b] = [wide, narrow]
        .each_ref()
        .map(|value| math::sign(value, dec));
    let [a, b] = [a?, b?];
    if a == Sign::Zero || b == Sign::Zero {
        return Err(invalid("a crossing pole alone does not select a band"));
    }
    def.branch = if a == Sign::Positive { 1 } else { -1 };
    def.ring = SaddleRing::Crossing {
        narrow_positive: b == Sign::Positive,
    };
    reach(&def, dec)?;
    Ok(def)
}

fn reach(def: &RollingSaddleDefinition, dec: &mut Decider) -> Result<(), Refusal> {
    let [a, b] = signs(def);
    let delta = &(&a * &b) * &def.cosine;
    let offset = def.offset(0);
    // The centre ellipse has minor P and major P sqrt(2/(1+a b C)).
    // Its reach is its minimum curvature radius, P sqrt((1+a b C)/2).
    // This is a geometric sufficient condition, not a fit/resource threshold.
    let twice = Real::from_i64(2);
    if math::sign(
        &(&(&offset.square() * &(&Real::one() + &delta)) - &(&twice * &def.radius.square())),
        dec,
    )? != Sign::Positive
    {
        return Err(invalid(
            "the crossing band is not proved inside the centre ellipse's regular tube",
        ));
    }
    Ok(())
}

pub(super) fn section(def: &RollingSaddleDefinition, u: &Rational) -> Result<[Point3; 3], Refusal> {
    if *u < q(-1) || *u > q(3) {
        return Err(invalid("crossing parameter is outside [-1,3]"));
    }
    let t = (u - &q(1)) / q(2);
    let denominator = Rational::one() + &t * &t;
    let sin = real(&(&t * &q(2) / denominator.clone()));
    let cos = real(&((Rational::one() - &t * &t) / denominator));
    let [a, b] = signs(def);
    let offset = def.offset(0);
    let x = &offset * &sin;
    let y = &(&offset * &a) * &cos;
    let z = &(&(&def.cosine * &y) - &(&(&offset * &b) * &cos)) / &def.sine;
    let centre = def.place(&x, &y, &z);
    let step = &def.wide - &offset;
    Ok([
        centre.clone(),
        def.place(&(&def.wide * &sin), &(&(&def.wide * &a) * &cos), &z),
        centre.plus(
            &def.across
                .scaled(&(&step * &sin))
                .plus(&def.side_axis().scaled(&(&(&step * &b) * &cos))),
        ),
    ])
}

fn scaled(values: &[Real], factor: &Real) -> Vec<Real> {
    values.iter().map(|value| value * factor).collect()
}

fn sum(a: &[Real], b: &[Real]) -> Vec<Real> {
    combine(a, b, &Real::one(), &Real::one())
}

fn difference(a: &[Real], b: &[Real]) -> Vec<Real> {
    combine(a, b, &Real::one(), &-Real::one())
}

fn derivative(values: &[Real]) -> Vec<Real> {
    let degree = Real::from_i64(values.len() as i64 - 1);
    values
        .windows(2)
        .map(|pair| &(&pair[1] - &pair[0]) * &degree)
        .collect()
}

fn zero_coefficients(values: &[Real], dec: &mut Decider) -> Result<(), Refusal> {
    for value in values {
        dec.check_deadline()?;
        if !equal(value, dec)? {
            return Err(unresolved(
                "a crossing carrier failed its whole-polynomial identity",
            ));
        }
    }
    Ok(())
}

// Tensor Bernstein multiplication, with v as the outer index.
fn times(a: &[Vec<Real>], b: &[Vec<Real>]) -> Vec<Vec<Real>> {
    let (m, n) = (a.len() - 1, b.len() - 1);
    let size = a[0].len() + b[0].len() - 1;
    let mut answer = vec![vec![Real::zero(); size]; m + n + 1];
    for (i, a) in a.iter().enumerate() {
        for (j, b) in b.iter().enumerate() {
            let scale = crate::rolling_envelope::binomial(m, i)
                * crate::rolling_envelope::binomial(n, j)
                / crate::rolling_envelope::binomial(m + n, i + j);
            let term = scaled(&product(a, b), &real(&scale));
            answer[i + j] = sum(&answer[i + j], &term);
        }
    }
    answer
}

// Check the actual homogeneous coefficients, not samples or a fit residual.
// For physical E = (d Hxyz - Cnum Hw)/(1,1,S), prove |E|² = r²(d Hw)²
// and E.C' = 0 on every (u,v), after clearing the positive S² denominator.
// Together with the positive circular-section weights and strict reach, these
// prove the complete band's two-sided distance bound is exactly zero.
fn certify(
    def: &RollingSaddleDefinition,
    rows: &Rows,
    centre: &[Vec<Real>; 3],
    denominator: &[Real],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    // The transient z coefficients are S times physical height. Clear S²
    // from both identities rather than eliminating the same frame root in
    // every coefficient. The represented carrier is unscaled below.
    let sine_sq = &Real::one() - &def.cosine.square();
    let weight: Vec<_> = rows.iter().map(|row| row[3].clone()).collect();
    let radial: [Vec<Vec<Real>>; 3] = std::array::from_fn(|k| {
        rows.iter()
            .map(|row| {
                difference(
                    &product(denominator, &row[k]),
                    &product(&centre[k], &row[3]),
                )
            })
            .collect()
    });
    let dw: Vec<_> = weight.iter().map(|w| product(denominator, w)).collect();
    let mut sphere = times(&dw, &dw)
        .iter()
        .map(|row| scaled(row, &-(&def.radius.square() * &sine_sq)))
        .collect::<Vec<_>>();
    for (k, axis) in radial.iter().enumerate() {
        for (total, mut term) in sphere.iter_mut().zip(times(axis, axis)) {
            if k != 2 {
                for value in &mut term {
                    *value = &*value * &sine_sq;
                }
            }
            *total = sum(total, &term);
        }
    }
    for row in &sphere {
        zero_coefficients(row, dec)?;
    }
    let mut normal = vec![vec![Real::zero(); 12]; 3];
    for k in 0..3 {
        let tangent = difference(
            &product(denominator, &derivative(&centre[k])),
            &product(&centre[k], &derivative(denominator)),
        );
        for (total, radial) in normal.iter_mut().zip(&radial[k]) {
            let mut term = product(radial, &tangent);
            if k != 2 {
                for value in &mut term {
                    *value = &*value * &sine_sq;
                }
            }
            *total = sum(total, &term);
        }
    }
    for row in &normal {
        zero_coefficients(row, dec)?;
    }
    Ok(())
}

fn span(def: &RollingSaddleDefinition, lo: i64, dec: &mut Decider) -> Result<Rows, Refusal> {
    let sections =
        catallax_numbers::OperationSpan::begin("shapes", &"equal-crossing", "fillet.sections", dec);
    let one = Real::one();
    let two = Real::from_i64(2);
    let t = [Real::from_i64(lo), Real::from_i64(lo + 1)];
    let square = product(&t, &t);
    let ones = vec![one.clone(); 3];
    let d = sum(&ones, &square);
    let cos = difference(&ones, &square);
    let sin = scaled(&elevate(&t, 2), &two);
    let d2 = product(&d, &d);
    let [a, b] = signs(def);
    let offset = def.offset(0);
    let side = Real::from_i64(def.sides[0]);
    let delta = &(&a * &b) * &def.cosine;
    let h = sum(
        &scaled(&product(&sin, &sin), &two),
        &scaled(&product(&cos, &cos), &(&one + &delta)),
    );
    let centre = [
        scaled(&sin, &offset),
        scaled(&cos, &(&offset * &a)),
        scaled(&cos, &(&offset * &(&(&def.cosine * &a) - &b))),
    ];
    let first = [
        scaled(&sin, &-&side),
        scaled(&cos, &-(&side * &a)),
        vec![Real::zero(); 3],
    ];
    let last = [
        scaled(&sin, &-&side),
        scaled(&cos, &-(&(&side * &b) * &def.cosine)),
        scaled(&cos, &(&(&side * &b) * &(&one - &def.cosine.square()))),
    ];
    let w0 = scaled(&product(&d2, &d), &two);
    let w1 = product(&h, &d);
    let rows = std::array::from_fn(|j| {
        std::array::from_fn(|k| {
            if k == 3 {
                return if j == 0 { w0.clone() } else { w1.clone() };
            }
            match j {
                0 => scaled(
                    &product(&d2, &sum(&centre[k], &scaled(&first[k], &def.radius))),
                    &two,
                ),
                1 => sum(
                    &product(&h, &centre[k]),
                    &scaled(&product(&d2, &sum(&first[k], &last[k])), &def.radius),
                ),
                _ => product(&h, &sum(&centre[k], &scaled(&last[k], &def.radius))),
            }
        })
    });
    sections.finish(dec, None, None);
    let certification = catallax_numbers::OperationSpan::begin(
        "shapes",
        &"equal-crossing",
        "fillet.certification",
        dec,
    );
    certify(def, &rows, &centre, &d, dec)?;
    certification.finish(dec, None, None);
    Ok(rows)
}

fn contact_queries(
    def: &RollingSaddleDefinition,
    carrier: &BSplineSurface,
    dec: &mut Decider,
) -> Result<Arc<[BSplineCurve; 2]>, Refusal> {
    let [_, wide, narrow] = section(def, &q(1))?;
    let ahead = def.across.scaled(&def.wide);
    let back = ahead.scaled(&-Real::one());
    let start = def.origin.plus(&back);
    let end = def.origin.plus(&ahead);
    let rail = |contact: Point3, side, dec: &mut Decider| {
        // Exact half-ellipse in the very same q=(u-1)/2 chart. Cancelling
        // 2(1+q²)² or h from a contact row removes no real point or pole.
        let query = BSplineCurve::new(
            2,
            [vec![q(-1); 3], vec![q(1); 2], vec![q(3); 3]].concat(),
            vec![
                start.clone(),
                contact.plus(&back),
                contact.clone(),
                contact.plus(&ahead),
                end.clone(),
            ],
            [2, 1, 1, 1, 2].into_iter().map(Real::from_i64).collect(),
            dec,
        )?;
        BSplineCurve::boundary(carrier, side).with_query_curve(query, dec)
    };
    Ok(Arc::new([
        rail(wide, SplineBoundary::VMin, dec)?,
        rail(narrow, SplineBoundary::VMax, dec)?,
    ]))
}

pub(super) fn build(
    definition: RollingSaddleDefinition,
    dec: &mut Decider,
) -> Result<RollingSaddle, Refusal> {
    reach(&definition, dec)?;
    let mut spans = [span(&definition, -1, dec)?, span(&definition, 0, dec)?];
    let sections =
        catallax_numbers::OperationSpan::begin("shapes", &"equal-crossing", "fillet.sections", dec);
    for (left, right) in spans[0].iter().zip(&spans[1]) {
        for (a, b) in left.iter().zip(right) {
            zero_coefficients(&[&a[6] - &b[0]], dec)?;
        }
    }
    // Restore physical height only after the frame-root-free coefficient proof.
    for rows in &mut spans {
        for row in rows {
            for z in &mut row[2] {
                *z = &*z / &definition.sine;
            }
        }
    }
    let mut controls = Vec::with_capacity(39);
    let mut weights = Vec::with_capacity(39);
    for (span, rows) in spans.iter().enumerate() {
        for i in usize::from(span != 0)..=6 {
            for row in rows {
                let weight = &row[3][i];
                controls.push(definition.place(
                    &(&row[0][i] / weight),
                    &(&row[1][i] / weight),
                    &(&row[2][i] / weight),
                ));
                weights.push(weight.clone());
            }
        }
    }
    let knot = |at, count| vec![q(at); count];
    let carrier = BSplineSurface::new(
        [6, 2],
        [
            [knot(-1, 7), knot(1, 6), knot(3, 7)].concat(),
            [knot(0, 3), knot(1, 3)].concat(),
        ],
        [13, 3],
        controls,
        weights,
        dec,
    )?;
    let net = Arc::new(Net::new(&definition, &carrier, dec)?);
    for rows in &mut spans {
        for row in &mut rows[1..] {
            let e = difference(
                &scaled(&row[1], &definition.cosine),
                &scaled(&row[2], &definition.sine),
            );
            let h = sum(
                &scaled(&row[1], &definition.sine),
                &scaled(&row[2], &definition.cosine),
            );
            row[1] = e;
            row[2] = h;
        }
    }
    sections.finish(dec, None, None);
    let certification = catallax_numbers::OperationSpan::begin(
        "shapes",
        &"equal-crossing",
        "fillet.certification",
        dec,
    );
    let cached_rails = Some(contact_queries(&definition, &carrier, dec)?);
    certification.finish(dec, None, None);
    Ok(RollingSaddle {
        data: Arc::new(Data {
            definition,
            construction: Construction::Crossing(Arc::new(Exact {
                knots: [q(-1), q(1), q(3)],
                rows: spans,
            })),
            carrier,
            cached_rails,
            net,
            bound: EvidenceBound {
                upper: Rational::zero(),
                unit: EvidenceUnit::Millimetres,
            },
            reversing: false,
        }),
    })
}

pub(super) fn parameter_table(
    saddle: &RollingSaddle,
    rail: usize,
    tolerance: &Rational,
    dec: &Decider,
) -> Result<crate::ParameterTable, Refusal> {
    if rail > 1 || *tolerance <= Rational::zero() {
        return Err(invalid(
            "a crossing table needs rail 0 or 1 and a positive tolerance",
        ));
    }
    let def = saddle.definition();
    let bits = saddle.data.net.bits;
    let [a, b_sign] = signs(def);
    let scale = if rail == 0 {
        &(&(&def.cosine * &a) - &b_sign) * &def.offset(0)
    } else {
        &(&a - &(&def.cosine * &b_sign)) * &def.offset(0)
    };
    let height_scale = dec.enclose_at(&(&scale / &def.sine), bits)?;
    let radius = dec.enclose_at(&def.wide, bits)?.lo();
    if radius <= Rational::zero() {
        return Err(unresolved(
            "a crossing table needs a positive cylinder-radius lower bound",
        ));
    }
    let one = b(&q(1), bits);
    let two = b(&q(2), bits);
    let half = b(&Rational::from_ratio(1, 2).expect("nonzero"), bits);
    let argument = |u: Ball| Jet::linear(&(&u - &one) * &half, half.clone());
    let height = |f: &Jet| {
        let square = f.square();
        Jet::constant(one.clone())
            .sub(&square)
            .div(&Jet::constant(one.clone()).add(&square))
            .map(|value| value.scale(&height_scale))
    };
    let mut pieces = 4;
    loop {
        let h = q(4) / q(pieces);
        let h4 = h.pow(4).expect("positive width");
        let mut angle_bound = Rational::zero();
        let mut height_bound = Rational::zero();
        for i in 0..pieces {
            dec.check_deadline()?;
            let lo = q(-1) + &h * &q(i);
            let hi = &lo + &h;
            let f = argument(b(&lo, bits).hull(&b(&hi, bits)));
            let slope = table::turning(&f, bits)?;
            angle_bound = angle_bound.max(magnitude(&slope.0[3]) * h4.clone() / q(64));
            height_bound = height_bound.max(magnitude(&height(&f)?.0[4]) * h4.clone() / q(16));
        }
        if &angle_bound + &(&height_bound / &radius) <= *tolerance {
            let mut table = crate::ParameterTable {
                knots: Vec::with_capacity(pieces as usize + 1),
                angle: Vec::with_capacity(pieces as usize + 1),
                height: Vec::with_capacity(pieces as usize + 1),
                bound: angle_bound,
                height_bound,
            };
            for i in 0..=pieces {
                dec.check_deadline()?;
                let u = q(-1) + &h * &q(i);
                let f = argument(b(&u, bits));
                let angle = &two * &f.0[0].atan_checked(dec)?;
                let slope = table::turning(&f, bits)?;
                let height = height(&f)?;
                table.bound = &table.bound + &(angle.rad() + slope.0[0].rad());
                table.height_bound = &table.height_bound + &(height.0[0].rad() + height.0[1].rad());
                table.knots.push(u);
                table.angle.push([angle.mid(), slope.0[0].mid()]);
                table.height.push([height.0[0].mid(), height.0[1].mid()]);
            }
            if &table.bound + &(&table.height_bound / &radius) <= *tolerance {
                return Ok(table);
            }
        }
        if pieces == 4096 {
            return Err(unresolved(
                "the crossing rail table misses its angle/height tolerance",
            ));
        }
        pieces *= 2;
    }
}

pub(super) fn same_ellipse(
    def: &RollingSaddleDefinition,
    scale: &Real,
    ellipse: &crate::Ellipse,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let [a, b] = signs(def);
    let major = def
        .along
        .scaled(&a)
        .plus(&def.axis.scaled(&(&(&(&def.cosine * &a) - &b) / &def.sine)));
    let normal = major.cross(&def.across);
    let gaps = [
        ellipse.center.minus(&def.origin).norm_squared(),
        ellipse.normal.cross(&normal).norm_squared(),
        &ellipse.minor_sq - &scale.square(),
        &ellipse.major_sq - &(&scale.square() * &major.norm_squared()),
        ellipse.major_dir.cross(&major).norm_squared(),
    ];
    match all_zero(&gaps, dec)? {
        Zero::Yes => Ok(true),
        Zero::No => Ok(false),
        Zero::Undecided(refusal) => Err(refusal),
    }
}

// Homogeneous power equations of the same antipodal projection. Keeping C
// symbolic here lets the existing three frame-root generators serve a spatial
// crossing too; its common perpendicular need not have rational length.
// Returned z is S times physical height, as in the coefficient proof.
pub(super) fn line_patch(
    def: &RollingSaddleDefinition,
    index: usize,
    nvars: usize,
    t: &catallax_numbers::MPoly,
    v: &catallax_numbers::MPoly,
    cosine: &catallax_numbers::MPoly,
) -> Option<[catallax_numbers::MPoly; 4]> {
    use catallax_numbers::MPoly;
    let radius = def.radius.as_rational()?;
    let offset = def.offset(0).as_rational()?;
    let one = MPoly::constant(nvars, q(1));
    let two = MPoly::constant(nvars, q(2));
    let tau = t + &MPoly::constant(nvars, q(index as i64 - 1));
    let square = &tau * &tau;
    let d = &one + &square;
    let sin = tau.scale(&q(2));
    let cos = &one - &square;
    let SaddleRing::Crossing { narrow_positive } = def.ring else {
        return None;
    };
    let [a, b] = [q(def.branch), q(if narrow_positive { 1 } else { -1 })];
    let side = q(def.sides[0]);
    let centre = [
        sin.scale(&offset),
        cos.scale(&(&offset * &a)),
        &cos * &(&cosine.scale(&a) - &one.scale(&b)).scale(&offset),
    ];
    let first = [
        sin.scale(&-&side),
        cos.scale(&-(&side * &a)),
        MPoly::zero(nvars),
    ];
    let last = [
        first[0].clone(),
        &cos * &cosine.scale(&-(&side * &b)),
        &cos * &(&one - &(cosine * cosine)).scale(&(&side * &b)),
    ];
    let start = &two - v;
    let dot = &(&sin * &sin) + &(&(&cos * &cos) * &cosine.scale(&(&a * &b)));
    let d2 = &d * &d;
    let ray: [MPoly; 3] = std::array::from_fn(|k| &(&start * &first[k]) + &(v * &last[k]));
    let meet = &(&start * &d2) + &(v * &dot);
    let divisor = &(&(&start * &start) + &(v * v)) * &d2;
    let divisor = &divisor + &(&(&start * v) * &dot).scale(&q(2));
    let mut homogeneous: [MPoly; 4] = std::array::from_fn(|_| MPoly::zero(nvars));
    for k in 0..3 {
        let radial = &(&meet * &ray[k]).scale(&q(2)) - &(&divisor * &first[k]);
        homogeneous[k] = &(&centre[k] * &divisor) + &radial.scale(&radius);
    }
    homogeneous[3] = &d * &divisor;
    Some(homogeneous)
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::BudgetPolicy;

    fn vector(x: i64, y: i64, z: i64) -> Vector3 {
        Vector3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
    }
    // A cylinder of radius sqrt(2) about an axis through the origin.
    fn cylinder(axis: Vector3, reference: Vector3) -> Cylinder {
        Cylinder {
            origin: Point3::new(Real::zero(), Real::zero(), Real::zero()),
            axis,
            radius_sq: Real::from_i64(2),
            ref_dir: reference,
        }
    }

    #[test]
    fn undecided_equal_radii_never_select_another_saddle_family() {
        // Equal irrational radii on crossing axes: the pole-to-pole band. The
        // radii's difference is an exact zero that a small allowance cannot
        // prove, and it once read as "unequal".
        let wide = cylinder(vector(0, 0, 1), vector(1, 0, 0));
        let narrow = cylinder(vector(1, 0, 0), vector(0, 1, 0));
        let point = Point3::new(Real::one(), Real::one(), Real::one());
        let radius = Real::from(Rational::from_ratio(1, 2).unwrap());
        let sides = [BallSide::Outside; 2];
        let exact = RollingSaddleDefinition::crossing(
            &wide,
            &narrow,
            &radius,
            sides,
            &point,
            &mut Decider::default(),
        )
        .unwrap();
        assert!(matches!(exact.ring, SaddleRing::Crossing { .. }));
        let mut refused = 0;
        for max_precision_bits in [16, 64, 256] {
            for max_work in (0..2_000).step_by(7) {
                let mut dec = Decider::new(BudgetPolicy {
                    max_precision_bits,
                    max_work,
                });
                match RollingSaddleDefinition::crossing(
                    &wide, &narrow, &radius, sides, &point, &mut dec,
                ) {
                    Ok(def) => assert_eq!(
                        (def.ring, def.branch),
                        (exact.ring, exact.branch),
                        "{max_precision_bits} bits, {max_work} work"
                    ),
                    Err(e) => {
                        assert!(e.code.starts_with("numbers."), "{e:?}");
                        refused += 1;
                    }
                }
            }
        }
        assert!(refused > 0, "the sweep must reach an undecided allowance");
    }

    #[test]
    fn undecided_centre_ellipse_identity_is_not_a_no() {
        let wide = cylinder(vector(0, 0, 1), vector(1, 0, 0));
        let narrow = cylinder(vector(1, 0, 0), vector(0, 1, 0));
        let point = Point3::new(Real::one(), Real::one(), Real::one());
        let radius = Real::from(Rational::from_ratio(1, 2).unwrap());
        let def = RollingSaddleDefinition::crossing(
            &wide,
            &narrow,
            &radius,
            [BallSide::Outside; 2],
            &point,
            &mut Decider::default(),
        )
        .unwrap();
        // The band's own centre ellipse, written out from its definition.
        let [a, b] = signs(&def);
        let major = def
            .along
            .scaled(&a)
            .plus(&def.axis.scaled(&(&(&(&def.cosine * &a) - &b) / &def.sine)));
        let scale = def.offset(0);
        let centres = crate::Ellipse {
            center: def.origin.clone(),
            normal: major.cross(&def.across),
            major_dir: major.clone(),
            major_sq: &scale.square() * &major.norm_squared(),
            minor_sq: scale.square(),
        };
        let shifted = crate::Ellipse {
            center: def.origin.plus(&vector(100, 0, 0)),
            ..centres.clone()
        };
        let mut dec = Decider::default();
        assert!(same_ellipse(&def, &scale, &centres, &mut dec).unwrap());
        assert!(!same_ellipse(&def, &scale, &shifted, &mut dec).unwrap());
        let mut refused = 0;
        for max_work in 0..400 {
            let limited = BudgetPolicy {
                max_precision_bits: 64,
                max_work,
            };
            match same_ellipse(&def, &scale, &centres, &mut Decider::new(limited)) {
                Ok(same) => assert!(same, "{max_work} work: the centre ellipse read as another"),
                Err(e) => {
                    assert!(e.code.starts_with("numbers."), "{e:?}");
                    refused += 1;
                }
            }
            // A proven nonzero gap still answers no after an undecided one.
            if let Ok(same) = same_ellipse(&def, &scale, &shifted, &mut Decider::new(limited)) {
                assert!(!same);
            }
        }
        assert!(refused > 0, "the sweep must reach an undecided allowance");
    }
}
