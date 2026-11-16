//! Circular sections joined by a C1 homogeneous cubic fit.
//!
//! The atlas proves only the carrier's regular charts and joins. Its agreement
//! with the rolling ball and support rails is numerical, checked at the stated
//! finite sample locations and retained separately from that atlas authority.

use super::*;
use crate::BSplineSurface;

const ROWS: usize = 16;
type Packed = [Rational; ROWS];

/// A numerical rolling-ball carrier and its checked native parameter rails.
#[derive(Clone, Debug)]
pub struct NumericalSplineFillet {
    spine: NumericalFilletSpine,
    surface: RationalAtlas,
    rails: [BSplineCurve; 2],
    parameter_rails: [BSplineCurve; 2],
    sampled_error_mm: Rational,
    samples: usize,
}

fn packed(problem: &Problem, node: &NumericalSpineNode, dec: &mut Decider) -> Result<Packed, Refusal> {
    let radius = problem.supports[0].distance.abs().to_f64_nearest();
    let mut points = [[0.0; 3]; 2];
    let mut directions = [[0.0; 3]; 2];
    for side in 0..2 {
        let support = &problem.supports[side];
        let at = node.parameters[side].each_ref().map(Rational::to_f64_nearest);
        let (_, chart) = support.chart(at)?;
        dec.charge_work(chart.evaluation_work())?;
        let (point, normal) = chart.contact_f64(at);
        points[side] = point;
        let sign = if support.distance.sign() == Sign::Negative { 1.0 } else { -1.0 };
        directions[side] = normal.map(|v| sign * v);
    }
    let cosine = (0..3).map(|i| directions[0][i] * directions[1][i]).sum::<f64>();
    let denominator = 1.0 + cosine;
    if !denominator.is_finite() || denominator <= 1e-12 {
        return Err(bad("rolling-ball section has opposite contact directions"));
    }
    let weight = (0.5 * denominator).sqrt();
    let centre = node.centre.each_ref().map(Rational::to_f64_nearest);
    let middle: [f64; 3] = std::array::from_fn(|i| {
        centre[i] + radius * (directions[0][i] + directions[1][i]) / denominator
    });
    let mut result: Packed = std::array::from_fn(|_| Rational::zero());
    for (row, (point, w)) in [(points[0], 1.0), (middle, weight), (points[1], 1.0)].into_iter().enumerate() {
        for i in 0..3 { result[row * 4 + i] = rational(point[i] * w)?; }
        result[row * 4 + 3] = rational(w)?;
    }
    for side in 0..2 {
        for axis in 0..2 { result[12 + side * 2 + axis] = node.parameters[side][axis].clone(); }
    }
    Ok(result)
}

/// Derivatives of the local three-station interpolating quadratic. In
/// particular, end derivatives are not first-order secants.
pub(super) fn station_slopes<'a, const N: usize>(
    parameter: impl Fn(usize) -> &'a Rational,
    values: &[[Rational; N]],
) -> Vec<[Rational; N]> {
    let mut slopes = Vec::with_capacity(values.len());
    for i in 0..values.len() {
        let a = i.saturating_sub(1).min(values.len() - 3);
        let h0 = parameter(a + 1) - parameter(a);
        let h1 = parameter(a + 2) - parameter(a + 1);
        let width = &h0 + &h1;
        slopes.push(std::array::from_fn(|j| {
            let d0 = (&values[a + 1][j] - &values[a][j]) / h0.clone();
            let d1 = (&values[a + 2][j] - &values[a + 1][j]) / h1.clone();
            let numerator = if i == 0 {
                &(&(&h0 + &width) * &d0) - &(&h0 * &d1)
            } else if i + 1 == values.len() {
                &(&(&h1 + &width) * &d1) - &(&h1 * &d0)
            } else {
                &(&h1 * &d0) + &(&h0 * &d1)
            };
            numerator / width.clone()
        }));
    }
    slopes
}

fn hermite(nodes: &[NumericalSpineNode], values: &[Packed]) -> (Vec<Rational>, Vec<Packed>) {
    let slopes = station_slopes(|i| &nodes[i].tau, values);
    let mut knots = vec![nodes[0].tau.clone(); 4];
    let mut controls = Vec::with_capacity(3 * (nodes.len() - 1) + 1);
    controls.push(values[0].clone());
    for i in 0..nodes.len() - 1 {
        let third = (&nodes[i + 1].tau - &nodes[i].tau) / Rational::from_i64(3);
        controls.push(std::array::from_fn(|j| &values[i][j] + &(&third * &slopes[i][j])));
        controls.push(std::array::from_fn(|j| &values[i + 1][j] - &(&third * &slopes[i + 1][j])));
        controls.push(values[i + 1].clone());
        if i + 1 < nodes.len() - 1 { knots.extend(std::iter::repeat_n(nodes[i + 1].tau.clone(), 3)); }
    }
    knots.extend(std::iter::repeat_n(nodes.last().expect("stations").tau.clone(), 4));
    (knots, controls)
}

struct Candidate {
    surface: BSplineSurface,
    rails: [BSplineCurve; 2],
    parameters: [BSplineCurve; 2],
}
fn candidate(nodes: &[NumericalSpineNode], values: &[Packed], dec: &mut Decider) -> Result<Candidate, Refusal> {
    dec.charge_work((nodes.len() as u64) * 128)?;
    let (knots, data) = hermite(nodes, values);
    let mut controls = Vec::with_capacity(3 * data.len());
    let mut weights = Vec::with_capacity(3 * data.len());
    let mut rail_points: [Vec<Point3>; 2] = std::array::from_fn(|_| Vec::with_capacity(data.len()));
    let mut rail_weights: [Vec<Real>; 2] = std::array::from_fn(|_| Vec::with_capacity(data.len()));
    let mut parameters: [Vec<Point3>; 2] = std::array::from_fn(|_| Vec::with_capacity(data.len()));
    for entry in &data {
        for row in 0..3 {
            let weight = &entry[4 * row + 3];
            if weight.sign() != Sign::Positive { return Err(bad("numerical carrier has a non-positive homogeneous weight")); }
            let point = Point3::new(
                Real::from_rational(&entry[4 * row] / weight),
                Real::from_rational(&entry[4 * row + 1] / weight),
                Real::from_rational(&entry[4 * row + 2] / weight),
            );
            if row != 1 {
                let side = usize::from(row == 2);
                rail_points[side].push(point.clone());
                rail_weights[side].push(Real::from_rational(weight.clone()));
            }
            controls.push(point);
            weights.push(Real::from_rational(weight.clone()));
        }
        for side in 0..2 {
            parameters[side].push(Point3::new(
                Real::from_rational(entry[12 + 2 * side].clone()),
                Real::from_rational(entry[13 + 2 * side].clone()), Real::zero(),
            ));
        }
    }
    let count = data.len();
    let surface = BSplineSurface::new(
        [3, 2],
        [knots.clone(), [0, 0, 0, 1, 1, 1].map(Rational::from_i64).to_vec()],
        [count, 3], controls, weights, dec,
    )?;
    let [a, b] = rail_points;
    let [wa, wb] = rail_weights;
    let rails = [
        BSplineCurve::new(3, knots.clone(), a, wa, dec)?,
        BSplineCurve::new(3, knots.clone(), b, wb, dec)?,
    ];
    let [a, b] = parameters;
    let parameters = [
        BSplineCurve::new(3, knots.clone(), a, vec![Real::one(); count], dec)?,
        BSplineCurve::new(3, knots, b, vec![Real::one(); count], dec)?,
    ];
    Ok(Candidate { surface, rails, parameters })
}
fn section(row: &Packed, parameter: &Rational) -> [Rational; 3] {
    let other = Rational::one() - parameter.clone();
    let basis = [other.clone() * other.clone(), Rational::from_i64(2) * other * parameter.clone(), parameter.clone() * parameter.clone()];
    let h: [Rational; 4] = std::array::from_fn(|i| (0..3).fold(Rational::zero(), |sum, j| sum + &basis[j] * &row[j * 4 + i]));
    std::array::from_fn(|i| &h[i] / &h[3])
}
fn distance(point: &[Ball; 3], target: &[Rational; 3]) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| {
        sum + (&point[i].lo() - &target[i]).abs().max((&point[i].hi() - &target[i]).abs())
    })
}
fn rail_error(problem: &Problem, candidate: &Candidate, tau: &Rational, side: usize, dec: &mut Decider) -> Result<Rational, Refusal> {
    let uv = candidate.parameters[side].evaluate_box([tau, tau], 128, dec)?;
    let at = [uv[0].mid(), uv[1].mid()];
    let support = &problem.supports[side];
    let (_, chart) = support.chart(at.each_ref().map(Rational::to_f64_nearest))?;
    let mut meter = Meter::default();
    let contact = chart.contact(&at.each_ref().map(|p| Jet::constant(ball(p), 0)), &mut meter)?;
    meter.charge(dec)?;
    let point = candidate.rails[side].evaluate_box([tau, tau], 128, dec)?;
    Ok((0..3).fold(Rational::zero(), |sum, i| {
        sum + (&point[i].lo() - &contact.point[i].0[0].hi()).abs()
            .max((&point[i].hi() - &contact.point[i].0[0].lo()).abs())
    }))
}

impl NumericalSplineFillet {
    /// Fits equal-radius circular sections to this numerical branch. Quarter
    /// points of every longitudinal span and five points of each section check
    /// the fit. Both native rails are checked on their original supports there.
    /// No complete deviation bound or exact rolling-ball authority is claimed.
    pub fn fit(mut spine: NumericalFilletSpine, dec: &mut Decider) -> Result<Self, Refusal> {
        let radius = spine.problem.supports[0].distance.abs();
        if radius.sign() != Sign::Positive || spine.problem.supports[1].distance.abs() != radius {
            return Err(bad("a rolling-ball carrier requires equal positive offset radii"));
        }
        let mut samples = 0;
        for round in 0..=MAX_DEPTH {
            let values = spine.nodes.iter().map(|node| packed(&spine.problem, node, dec)).collect::<Result<Vec<_>, _>>()?;
            let fitted = candidate(&spine.nodes, &values, dec)?;
            let mut inserted = Vec::new();
            let mut upper = Rational::zero();
            let mut centre_upper = Rational::zero();
            for (index, pair) in spine.nodes.windows(2).enumerate() {
                let mut span_upper = Rational::zero();
                let mut middle = None;
                for i in 1..4 {
                    let share = Rational::from_ratio(i, 4).unwrap();
                    let tau = &pair[0].tau + &(&(&pair[1].tau - &pair[0].tau) * &share);
                    let node = spine.problem.solve(&tau, spine.problem.guess(&pair[0], &pair[1], i as f64 / 4.0), dec)?;
                    let centre_error = (0..3).fold(Rational::zero(), |sum, k| {
                        let chord = &pair[0].centre[k] + &(&(&pair[1].centre[k] - &pair[0].centre[k]) * &share);
                        sum + (&node.centre[k] - &chord).abs()
                    });
                    centre_upper = centre_upper.max(centre_error);
                    let truth = packed(&spine.problem, &node, dec)?;
                    for j in 0..=4 {
                        let section_at = Rational::from_ratio(j, 4).unwrap();
                        let actual = fitted.surface.evaluate_box([&tau, &tau], [&section_at, &section_at], 128, dec)?;
                        let error = distance(&actual, &section(&truth, &section_at));
                        span_upper = span_upper.max(error);
                        samples += 1;
                    }
                    for side in 0..2 { span_upper = span_upper.max(rail_error(&spine.problem, &fitted, &tau, side, dec)?); }
                    if i == 2 { middle = Some(node); }
                }
                if span_upper > spine.tolerance_mm.mul_pow2(-1) {
                    inserted.push((index + 1, middle.expect("checked midpoint")));
                }
                upper = upper.max(span_upper);
            }
            if inserted.is_empty() {
                let [surface] = <[_; 1]>::try_from(RationalAtlas::partition(fitted.surface, [false, false], dec)?)
                    .map_err(|_| bad("numerical fillet carrier does not form one regular G1 logical face"))?;
                spine.sampled_error_mm = centre_upper;
                spine.samples += samples / 5;
                return Ok(Self { spine, surface, rails: fitted.rails, parameter_rails: fitted.parameters, sampled_error_mm: upper, samples });
            }
            if round == MAX_DEPTH || spine.nodes.len() + inserted.len() > MAX_NODES {
                return Err(bad("numerical fillet carrier cannot meet its stated sampled tolerance within the subdivision limit"));
            }
            for (index, node) in inserted.into_iter().rev() { spine.nodes.insert(index, node); }
        }
        unreachable!("fit returns at its subdivision limit")
    }
    /// Exact regular atlas of the fitted carrier, not an exact rolling-ball face.
    pub fn surface(&self) -> &RationalAtlas { &self.surface }
    /// Numerical spine and its explicitly stated tolerance.
    pub fn spine(&self) -> &NumericalFilletSpine { &self.spine }
    /// Three-dimensional boundary curve on the first or second support.
    pub fn rail(&self, side: usize) -> Option<&BSplineCurve> { self.rails.get(side) }
    /// Native support parameters corresponding to a boundary's curve parameter.
    pub fn parameter_rail(&self, side: usize) -> Option<&BSplineCurve> { self.parameter_rails.get(side) }
    /// Largest checked position or rail residual, in L1 millimetres.
    pub fn sampled_error_mm(&self) -> &Rational { &self.sampled_error_mm }
    /// Number of interior section sample positions checked during fitting.
    pub fn samples(&self) -> usize { self.samples }
}
