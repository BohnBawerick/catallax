//! Finite subdivision and Newton refinement without exact incidence authority.

use catallax_numbers::{
    Ball, Decider, NumericalCheck, NumericalEvidence, NumericalSample, Rational, Refusal, Sign,
    operation_budget_active, with_operation_budget,
};
use catallax_shapes::{BSplineCurve, CurveId, GeomStore, Line, Point3};

const BITS: u32 = 128;
const MAX_DEPTH: usize = 64;

/// One numerical meeting, retaining both parameters and a rational point.
#[derive(Clone, Debug)]
pub struct NumericalSplineLinePoint {
    /// Parameter in the spline's closed finite knot domain.
    pub spline_parameter: Rational,
    /// Parameter on the original, unbounded line.
    pub line_parameter: Rational,
    /// Approximate position, checked against both retained definitions.
    pub position: Point3,
}

/// Numerical meetings and the checks actually performed. No exact facts are added.
#[derive(Clone, Debug)]
pub struct NumericalSplineLineMeetings {
    points: Vec<NumericalSplineLinePoint>,
    evidence: NumericalEvidence,
    cells: usize,
}

impl NumericalSplineLineMeetings {
    /// Checked approximate points, ordered by spline parameter.
    pub fn points(&self) -> &[NumericalSplineLinePoint] {
        &self.points
    }
    /// Preset tolerance and residual samples. These are not an exact root count.
    pub fn evidence(&self) -> &NumericalEvidence {
        &self.evidence
    }
    /// Number of parameter intervals examined, including rejected intervals.
    pub fn cells(&self) -> usize {
        self.cells
    }
}

fn refusal(detail: &str) -> Refusal {
    Refusal::new("intersect.numerical-spline-line", detail)
}

fn scalar(value: &Rational) -> Ball {
    Ball::from_rational(value, BITS)
}
fn upper(value: &Ball) -> Rational {
    value.lo().abs().max(value.hi().abs())
}
fn contains_zero(value: &Ball) -> bool {
    value.lo().sign() != Sign::Positive && value.hi().sign() != Sign::Negative
}
fn midpoint(lo: &Rational, hi: &Rational) -> Rational {
    (lo + hi) / Rational::from_i64(2)
}

struct Frame {
    origin: [Ball; 3],
    direction: [Ball; 3],
    axis: usize,
    other: [usize; 2],
    ratio: [Ball; 2],
}

impl Frame {
    fn new(line: &Line, dec: &Decider) -> Result<Self, Refusal> {
        let origin = [&line.origin.x, &line.origin.y, &line.origin.z]
            .map(|x| dec.enclose_at(x, BITS))
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?;
        let direction = [&line.dir.x, &line.dir.y, &line.dir.z]
            .map(|x| dec.enclose_at(x, BITS))
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?;
        let axis = (0..3)
            .filter(|&i| !contains_zero(&direction[i]))
            .max_by_key(|&i| upper(&direction[i]))
            .ok_or_else(|| refusal("line direction is zero or unresolved at 128 bits"))?;
        let other: [usize; 2] = (0..3)
            .filter(|&i| i != axis)
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        let ratio = other.map(|i| direction[i].checked_div(&direction[axis]).unwrap());
        Ok(Self {
            origin: origin.try_into().unwrap(),
            direction: direction.try_into().unwrap(),
            axis,
            other,
            ratio,
        })
    }

    fn transverse(&self, p: &[Ball; 3], derivative: bool) -> [Ball; 2] {
        let offset: [Ball; 3] = std::array::from_fn(|i| {
            if derivative {
                p[i].clone()
            } else {
                &p[i] - &self.origin[i]
            }
        });
        std::array::from_fn(|i| &offset[self.other[i]] - &(&self.ratio[i] * &offset[self.axis]))
    }
}

struct Cell {
    lo: Rational,
    hi: Rational,
    depth: usize,
    transverse: [Vec<Ball>; 2],
}

fn split(controls: &[Ball], dec: &mut Decider) -> Result<[Vec<Ball>; 2], Refusal> {
    let half = scalar(&Rational::from_ratio(1, 2).unwrap());
    let mut row = controls.to_vec();
    let mut left = vec![row[0].clone()];
    let mut right = vec![row[row.len() - 1].clone()];
    while row.len() > 1 {
        dec.charge_work(row.len() as u64)?;
        row = row
            .windows(2)
            .map(|pair| &(&pair[0] + &pair[1]) * &half)
            .collect();
        left.push(row[0].clone());
        right.push(row[row.len() - 1].clone());
    }
    right.reverse();
    Ok([left, right])
}
struct Candidate {
    lo: Rational,
    hi: Rational,
    parameter: Rational,
}

fn refine(
    curve: &BSplineCurve,
    frame: &Frame,
    cell: &Cell,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let mut at = midpoint(&cell.lo, &cell.hi);
    for _ in 0..12 {
        dec.charge_work(32)?;
        let jet = curve.jet_box([&at, &at], BITS, dec)?;
        let values = frame.transverse(&jet.value, false);
        let derivatives = frame.transverse(&jet.first, true);
        let Some(axis) = (0..2)
            .filter(|&i| !contains_zero(&derivatives[i]))
            .max_by_key(|&i| upper(&derivatives[i]))
        else {
            break;
        };
        let step = values[axis].checked_div(&derivatives[axis]).unwrap().mid();
        let next = (&at - &step).max(cell.lo.clone()).min(cell.hi.clone());
        if next == at {
            break;
        }
        at = next;
    }
    Ok(at)
}

/// Numerically intersects a native finite spline with a line in either order.
/// Uses the preset 1e-7 mm tolerance, 128-bit interval evaluations, subdivision
/// and clamped Newton refinement. Overlap intervals and unresolved cells refuse.
/// This query never publishes incidence facts or supplies exact Boolean points.
pub fn numerical_spline_line_points(
    geom: &GeomStore,
    a: CurveId,
    b: CurveId,
    dec: &mut Decider,
) -> Result<Option<NumericalSplineLineMeetings>, Refusal> {
    if operation_budget_active() {
        return query(geom, a, b, dec);
    }
    let (result, _) = with_operation_budget(20_000_000, || query(geom, a, b, dec));
    result
}

fn query(
    geom: &GeomStore,
    a: CurveId,
    b: CurveId,
    dec: &mut Decider,
) -> Result<Option<NumericalSplineLineMeetings>, Refusal> {
    let Some((curve, line)) = super::pair(geom, a, b) else {
        return Ok(None);
    };
    let tolerance = Rational::parse("1/10000000").unwrap();
    let leaf_size = &tolerance / &Rational::from_i64(8);
    let flat_limit = &tolerance / &Rational::from_i64(1_000_000);
    let frame = Frame::new(line, dec)?;
    let mut pending = Vec::new();
    for (span, [lo, hi]) in curve.spans() {
        let hull =
            crate::span_hull::SpanHull::new(curve, span, [lo.clone().into(), hi.clone().into()]);
        let controls = hull.span_controls(BITS, dec)?;
        let transverse = std::array::from_fn(|j| {
            controls
                .iter()
                .map(|p| {
                    let other = &p[frame.other[j]] - &(&frame.origin[frame.other[j]] * &p[3]);
                    let axis = &p[frame.axis] - &(&frame.origin[frame.axis] * &p[3]);
                    &other - &(&frame.ratio[j] * &axis)
                })
                .collect()
        });
        pending.push(Cell {
            lo: lo.clone(),
            hi: hi.clone(),
            depth: 0,
            transverse,
        });
    }
    let mut candidates = Vec::<Candidate>::new();
    let mut cells = 0;
    while let Some(cell) = pending.pop() {
        dec.charge_work(64)?;
        cells += 1;
        if cell.transverse.iter().any(|net| {
            net.iter().all(|p| p.lo().sign() == Sign::Positive)
                || net.iter().all(|p| p.hi().sign() == Sign::Negative)
        }) {
            continue;
        }
        let jet = curve.jet_box([&cell.lo, &cell.hi], BITS, dec)?;
        let residual = frame.transverse(&jet.value, false);
        if residual.iter().any(|r| !contains_zero(r)) {
            continue;
        }
        let diameter = jet
            .value
            .iter()
            .fold(Rational::zero(), |sum, x| sum + x.hi() - x.lo());
        let derivatives = frame.transverse(&jet.first, true);
        if diameter > leaf_size
            && residual.iter().all(|r| upper(r) <= flat_limit)
            && derivatives
                .iter()
                .all(|r| upper(r) * (&cell.hi - &cell.lo) <= flat_limit)
        {
            return Err(Refusal::new(
                "intersect.spline-line-overlap",
                "a finite interval is numerically indistinguishable from line overlap; interval ownership is required",
            ));
        }
        if diameter <= leaf_size {
            let parameter = refine(curve, &frame, &cell, dec)?;
            candidates.push(Candidate {
                lo: cell.lo,
                hi: cell.hi,
                parameter,
            });
        } else if cell.depth >= MAX_DEPTH {
            return Err(refusal(
                "finite interval remains unresolved after 64 subdivisions",
            ));
        } else {
            let middle = midpoint(&cell.lo, &cell.hi);
            let [left_a, right_a] = split(&cell.transverse[0], dec)?;
            let [left_b, right_b] = split(&cell.transverse[1], dec)?;
            pending.push(Cell {
                lo: middle.clone(),
                hi: cell.hi,
                depth: cell.depth + 1,
                transverse: [right_a, right_b],
            });
            pending.push(Cell {
                lo: cell.lo,
                hi: middle,
                depth: cell.depth + 1,
                transverse: [left_a, left_b],
            });
        }
    }
    candidates.sort_by(|a, b| a.lo.cmp(&b.lo));
    let mut merged = Vec::<Candidate>::new();
    for candidate in candidates {
        if let Some(last) = merged.last_mut().filter(|last| last.hi >= candidate.lo) {
            last.hi = last.hi.clone().max(candidate.hi);
        } else {
            merged.push(candidate);
        }
    }
    let mut points = Vec::new();
    let mut curve_samples = Vec::new();
    let mut line_samples = Vec::new();
    for candidate in merged {
        dec.charge_work(64)?;
        let at = candidate.parameter;
        let p = curve.evaluate_box([&at, &at], BITS, dec)?;
        let position = p.each_ref().map(Ball::mid);
        let t = (&p[frame.axis] - &frame.origin[frame.axis])
            .checked_div(&frame.direction[frame.axis])
            .unwrap()
            .mid();
        let q: [Ball; 3] =
            std::array::from_fn(|i| &frame.origin[i] + &(&frame.direction[i] * &scalar(&t)));
        let distance = |target: &[Ball; 3]| {
            (0..3).fold(Rational::zero(), |sum, i| {
                sum + upper(&(&scalar(&position[i]) - &target[i]))
            })
        };
        curve_samples.push(NumericalSample {
            parameters: vec![at.clone()],
            distance_upper_mm: distance(&p),
        });
        line_samples.push(NumericalSample {
            parameters: vec![at.clone(), t.clone()],
            distance_upper_mm: distance(&q),
        });
        points.push(NumericalSplineLinePoint {
            spline_parameter: at,
            line_parameter: t,
            position: Point3::new(
                position[0].clone().into(),
                position[1].clone().into(),
                position[2].clone().into(),
            ),
        });
    }
    let mut checks = vec![NumericalCheck::Passed {
        name: "finite-knot-domain-subdivision".into(),
        scope: format!("spline-line:{a:?}:{b:?};cells={cells}"),
    }];
    if !points.is_empty() {
        checks.push(NumericalCheck::Samples {
            name: "point-on-retained-spline".into(),
            scope: format!("{a:?}:{b:?}"),
            samples: curve_samples,
        });
        checks.push(NumericalCheck::Samples {
            name: "point-on-retained-line".into(),
            scope: format!("{a:?}:{b:?}"),
            samples: line_samples,
        });
    }
    let evidence = NumericalEvidence::new(tolerance, checks, vec![])?;
    Ok(Some(NumericalSplineLineMeetings {
        points,
        evidence,
        cells,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{Real, with_operation_budget};
    use catallax_shapes::{Curve, Vector3};

    fn r(value: &str) -> Real {
        Rational::parse(value).unwrap().into()
    }
    fn point(x: &str, y: &str) -> Point3 {
        Point3::new(r(x), r(y), Real::zero())
    }

    #[test]
    fn cubic_with_algebraic_line_builds_numerically_in_both_orders() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let mut geom = GeomStore::new();
            let curve = BSplineCurve::new(
                3,
                [0, 0, 0, 0, 1, 1, 1, 1].map(Rational::from_i64).to_vec(),
                vec![
                    point("0", "0"),
                    point("1/3", "0"),
                    point("2/3", "0"),
                    point("1", "1"),
                ],
                vec![Real::one(); 4],
                &mut dec,
            )
            .unwrap();
            let curve = geom.add_curve(Curve::BSpline(curve));
            let height = Real::from_i64(2).sqrt() / Real::from_i64(4);
            let line = geom.add_curve(Curve::Line(Line {
                origin: Point3::new(Real::zero(), height.clone(), Real::zero()),
                dir: Vector3::new(Real::from_i64(2).sqrt(), Real::zero(), Real::zero()),
            }));
            assert!(super::super::native_spline_line_points(&geom, curve, line, &mut dec).is_err());
            for (a, b) in [(curve, line), (line, curve)] {
                let result = numerical_spline_line_points(&geom, a, b, &mut dec)
                    .unwrap()
                    .unwrap();
                assert_eq!(result.points().len(), 1);
                assert_eq!(
                    result.evidence().tolerance_mm(),
                    &Rational::parse("1/10000000").unwrap()
                );
                assert_eq!(result.evidence().checks().len(), 3);
                let p = &result.points()[0];
                assert!(
                    p.spline_parameter > Rational::zero() && p.spline_parameter < Rational::one()
                );
                let error = dec
                    .enclose_at(
                        &(&p.position.x - &(height.clone() * Real::from_i64(2))),
                        BITS,
                    )
                    .unwrap();
                assert!(upper(&error) < Rational::parse("1/10000000").unwrap());
            }
            for z in ["1/1000000", "-1/1000000"] {
                let miss = geom.add_curve(Curve::Line(Line {
                    origin: Point3::new(Real::zero(), height.clone(), r(z)),
                    dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
                }));
                assert!(
                    numerical_spline_line_points(&geom, curve, miss, &mut dec)
                        .unwrap()
                        .unwrap()
                        .points()
                        .is_empty()
                );
            }
            let outside = geom.add_curve(Curve::Line(Line {
                origin: point("0", "2"),
                dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
            }));
            assert!(
                numerical_spline_line_points(&geom, curve, outside, &mut dec)
                    .unwrap()
                    .unwrap()
                    .points()
                    .is_empty()
            );
        });
        assert!(!work.exhausted, "{} work units", work.used);
        eprintln!(
            "algebraic cubic numerical controls: {} work units",
            work.used
        );
    }

    #[test]
    fn subdivision_keeps_all_three_components_and_refuses_partial_overlap() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let mut geom = GeomStore::new();
            let curve = BSplineCurve::new(
                3,
                [0, 0, 0, 0, 1, 1, 1, 1].map(Rational::from_i64).to_vec(),
                vec![
                    point("0", "-3/32"),
                    point("1/3", "13/96"),
                    point("2/3", "-13/96"),
                    point("1", "3/32"),
                ],
                vec![Real::one(); 4],
                &mut dec,
            )
            .unwrap();
            let curve = geom.add_curve(Curve::BSpline(curve));
            let line = geom.add_curve(Curve::Line(Line {
                origin: point("0", "0"),
                dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
            }));
            let result = numerical_spline_line_points(&geom, curve, line, &mut dec)
                .unwrap()
                .unwrap();
            assert_eq!(result.points().len(), 3);
            for (p, expected) in result.points().iter().zip(["1/4", "1/2", "3/4"]) {
                assert!(
                    (&p.spline_parameter - &Rational::parse(expected).unwrap()).abs()
                        < Rational::parse("1/10000000").unwrap()
                );
            }
            let overlap = BSplineCurve::new(
                1,
                [0, 0, 1, 2, 2].map(Rational::from_i64).to_vec(),
                vec![point("0", "0"), point("1", "0"), point("2", "1")],
                vec![Real::one(); 3],
                &mut dec,
            )
            .unwrap();
            let overlap = geom.add_curve(Curve::BSpline(overlap));
            assert_eq!(
                numerical_spline_line_points(&geom, overlap, line, &mut dec)
                    .unwrap_err()
                    .code,
                "intersect.spline-line-overlap"
            );
        });
        assert!(!work.exhausted, "{} work units", work.used);
        eprintln!(
            "three-component and partial-overlap numerical controls: {} work units",
            work.used
        );
    }

    #[test]
    fn tangent_and_tiny_crossing_survive_pose_scale_and_finite_domain_checks() {
        use catallax_shapes::{Pose, Rotation};
        for scale in [1, 7] {
            let ((), work) = with_operation_budget(20_000_000, || {
                let mut dec = Decider::default();
                let mut geom = GeomStore::new();
                let rotation = Rotation::from_quaternion(
                    Rational::one(),
                    Rational::one(),
                    Rational::one(),
                    Rational::zero(),
                )
                .unwrap();
                let pose =
                    Pose::new(rotation, Vector3::new(r("-7"), r("3"), r("11")), &mut dec).unwrap();
                let s = Real::from_i64(scale);
                let controls = [(-1, 1), (0, -1), (1, 1)]
                    .map(|(x, y)| {
                        Point3::new(
                            Real::from_i64(x * scale),
                            Real::from_i64(y * scale),
                            Real::zero(),
                        )
                    })
                    .to_vec();
                let curve = BSplineCurve::new(
                    2,
                    [0, 0, 0, 1, 1, 1].map(Rational::from_i64).to_vec(),
                    controls,
                    vec![Real::one(); 3],
                    &mut dec,
                )
                .unwrap()
                .transformed(&pose);
                let curve = geom.add_curve(Curve::BSpline(curve));
                for (height, expected) in
                    [("0", 1), ("1/1000000000000", 2), ("-1/1000000000000", 0)]
                {
                    let line = geom.add_curve(Curve::Line(Line {
                        origin: pose.apply_point(&Point3::new(
                            Real::zero(),
                            r(height) * s.clone(),
                            Real::zero(),
                        )),
                        dir: pose.apply_vector(&Vector3::new(
                            s.clone(),
                            Real::zero(),
                            Real::zero(),
                        )),
                    }));
                    let result = numerical_spline_line_points(&geom, line, curve, &mut dec)
                        .unwrap()
                        .unwrap();
                    assert_eq!(
                        result.points().len(),
                        expected,
                        "scale={scale}, height={height}"
                    );
                }
            });
            assert!(!work.exhausted, "scale={scale}: {} units", work.used);
            eprintln!(
                "posed tangent controls, scale={scale}: {} work units",
                work.used
            );
        }
    }

    #[test]
    fn zero_work_cannot_return_numerical_points() {
        let mut dec = Decider::default();
        let mut geom = GeomStore::new();
        let curve = BSplineCurve::new(
            1,
            [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
            vec![point("0", "-1"), point("1", "1")],
            vec![Real::one(); 2],
            &mut dec,
        )
        .unwrap();
        let curve = geom.add_curve(Curve::BSpline(curve));
        let line = geom.add_curve(Curve::Line(Line {
            origin: point("0", "0"),
            dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
        }));
        let (answer, _) = with_operation_budget(0, || {
            numerical_spline_line_points(&geom, curve, line, &mut dec)
        });
        assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
    }
}
