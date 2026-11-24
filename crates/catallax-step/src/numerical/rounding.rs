//! The rounding error of a written numerical file.
//!
//! Check once (roadmap item 3.0A, step 5; Checkpoint 2 decisions, section
//! 4.1 item 5): "Export writes the numbers and adds the rounding error to
//! the tolerance; it does not re-admit the body." This module is that
//! rounding error: how far the rounding of the written numbers can move a
//! point of a written curve or surface, to first order.
//!
//! The method is the forward error bound of a perturbed input (Higham,
//! "Accuracy and Stability of Numerical Algorithms", 2nd edition, section
//! 1.6): the largest error of each kind of written number, times how far a
//! unit change of such a number can move a point. The second factors are the
//! levers below, read from the values the file holds.
//!
//! OpenCASCADE and SMLib write their numbers and declare the tolerance they
//! hold; neither measures the written values again. This module follows
//! them and adds the one term they leave out.

use super::*;
use catallax_shapes::{BSplineCurve, BSplineSurface};

/// The largest rounding error of each kind of written number.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RoundingErrors {
    /// A coordinate or a radius, in millimetres.
    pub length_mm: Rational,
    /// A line's magnitude: its length per parameter unit.
    pub line_magnitude: Rational,
    /// A direction ratio, or a cone's semi-angle in radians.
    pub direction: Rational,
    /// A knot, a trim parameter, a u/v coordinate or a coefficient of a u/v
    /// line, in that parameter's own unit.
    pub parameter: Rational,
    /// A control weight.
    pub weight: Rational,
}

/// The sizes of a body that turn a rounding error into a distance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RoundingLevers {
    /// `M`: the largest coordinate or radius, in millimetres. Every written
    /// point lies in the cube of half side `M` about the origin.
    pub reach_mm: f64,
    /// `T`: the largest edge parameter.
    pub parameter: f64,
    /// `V`: the largest speed of a curve or surface, in millimetres per unit
    /// of its parameter.
    pub speed_mm: f64,
    /// `U`: the largest speed of a written face copy, in u/v units per unit
    /// of the edge parameter.
    pub face_copy_speed: f64,
    /// `Q`: the largest u or v of a written face copy.
    pub face_copy_reach: f64,
    /// `w`: the smallest control weight, 1 where there is none.
    pub least_weight: f64,
}

/// The rounding error of one written file.
#[derive(Clone, Debug, PartialEq)]
pub struct RoundingError {
    /// The largest error of each kind of written number.
    pub errors: RoundingErrors,
    /// The body's levers.
    pub levers: RoundingLevers,
    /// The bound on how far rounding moves a point of the written model, in
    /// millimetres. Export adds it to the largest entity tolerance.
    pub bound_mm: Rational,
}

fn size(value: &Real) -> f64 {
    let enclosure = value.approx();
    enclosure.lo().abs().max(enclosure.hi().abs())
}

fn point_size(point: &Point3) -> f64 {
    size(&point.x).max(size(&point.y)).max(size(&point.z))
}

fn length(vector: [f64; 3]) -> f64 {
    (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt()
}

/// The largest speed of a spline, and its smallest weight.
///
/// The derivative of a B-spline of degree `p` is a B-spline whose controls
/// are `p (P[i+1] - P[i]) / (u[i+p+1] - u[i+1])`, and it lies in their hull
/// (Piegl and Tiller, "The NURBS Book", equations 3.4 to 3.6). A rational
/// Bezier curve is faster than its control steps by at most the square of
/// its largest weight over its smallest (Floater, "Derivatives of rational
/// Bezier curves", Computer Aided Geometric Design 9, 1992). The same
/// factor is taken here for a rational spline: a lever of a first-order
/// bound, not a proof.
fn spline_speed(
    degree: usize,
    knots: &[Rational],
    controls: &[[f64; 3]],
    stride: usize,
    weights: &[f64],
) -> (f64, f64) {
    let mut fastest = 0.0f64;
    for (index, pair) in controls
        .iter()
        .zip(controls.iter().skip(stride))
        .enumerate()
    {
        let row = index / stride;
        let span = (&knots[row + degree + 1] - &knots[row + 1]).to_f64_nearest();
        if span > 0.0 {
            let step = std::array::from_fn(|axis| pair.1[axis] - pair.0[axis]);
            fastest = fastest.max(degree as f64 * length(step) / span);
        }
    }
    let least = weights.iter().copied().fold(f64::INFINITY, f64::min);
    let most = weights.iter().copied().fold(0.0f64, f64::max);
    let ratio = if least > 0.0 && least.is_finite() {
        most / least
    } else {
        1.0
    };
    (fastest * ratio * ratio, least)
}

fn curve_controls(curve: &BSplineCurve) -> (Vec<[f64; 3]>, Vec<f64>) {
    (
        curve
            .controls()
            .map(|p| [p.x.to_f64(), p.y.to_f64(), p.z.to_f64()])
            .collect(),
        curve.weights().map(size).collect(),
    )
}

impl RoundingLevers {
    fn take_curve(&mut self, curve: &BSplineCurve, dec: &mut Decider) -> Result<f64, Refusal> {
        dec.charge_work(curve.control_count() as u64)?;
        let (controls, weights) = curve_controls(curve);
        let (speed, least) = spline_speed(curve.degree(), curve.knots(), &controls, 1, &weights);
        self.least_weight = self.least_weight.min(least);
        Ok(speed)
    }

    fn take_surface(&mut self, surface: &BSplineSurface, dec: &mut Decider) -> Result<(), Refusal> {
        dec.charge_work(surface.controls().len() as u64)?;
        let controls: Vec<[f64; 3]> = surface
            .controls()
            .iter()
            .map(|p| {
                self.reach_mm = self.reach_mm.max(point_size(p));
                [p.x.to_f64(), p.y.to_f64(), p.z.to_f64()]
            })
            .collect();
        let weights: Vec<f64> = surface.weights().iter().map(size).collect();
        let [degree_u, degree_v] = surface.degrees();
        let [knots_u, knots_v] = surface.knots();
        let columns = surface.control_count()[1];
        // Controls are u-major: a step of one row is a step in u.
        let (along_u, least) = spline_speed(degree_u, knots_u, &controls, columns, &weights);
        let mut along_v = 0.0f64;
        for row in controls.chunks(columns) {
            let row_weights = vec![1.0; row.len()];
            along_v = along_v.max(spline_speed(degree_v, knots_v, row, 1, &row_weights).0);
        }
        let most = weights.iter().copied().fold(0.0f64, f64::max);
        let ratio = if least > 0.0 && least.is_finite() {
            most / least
        } else {
            1.0
        };
        self.speed_mm = self.speed_mm.max(along_u).max(along_v * ratio * ratio);
        self.least_weight = self.least_weight.min(least);
        Ok(())
    }

    /// The levers of the values a file holds: its tables, and the u/v
    /// curves it retains. Two bodies that write the same file have the same
    /// levers, whatever parameters they were built with.
    pub(super) fn of(
        definition: &NumericalDefinition,
        native_trims: &[NumericalNativeTrim],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut levers = Self {
            reach_mm: 0.0,
            parameter: 0.0,
            speed_mm: 0.0,
            face_copy_speed: 0.0,
            face_copy_reach: 0.0,
            least_weight: 1.0,
        };
        dec.charge_work((definition.vertices.len() + definition.edges.len()) as u64)?;
        for vertex in &definition.vertices {
            for coordinate in &vertex.position {
                levers.reach_mm = levers.reach_mm.max(coordinate.to_f64_nearest().abs());
            }
        }
        for edge in &definition.edges {
            for parameter in &edge.parameters {
                levers.parameter = levers.parameter.max(parameter.to_f64_nearest().abs());
            }
            match &edge.curve {
                Curve::Line(line) => {
                    levers.reach_mm = levers.reach_mm.max(point_size(&line.origin));
                    levers.speed_mm = levers.speed_mm.max(length([
                        size(&line.dir.x),
                        size(&line.dir.y),
                        size(&line.dir.z),
                    ]));
                }
                Curve::Circle(circle) => {
                    let radius = size(&circle.radius_sq).sqrt();
                    levers.reach_mm = levers.reach_mm.max(point_size(&circle.center)).max(radius);
                    // A circle's parameter is its angle in radians.
                    levers.speed_mm = levers.speed_mm.max(radius);
                }
                Curve::BSpline(curve) => {
                    for control in curve.controls() {
                        levers.reach_mm = levers.reach_mm.max(point_size(control));
                    }
                    let speed = levers.take_curve(curve, dec)?;
                    levers.speed_mm = levers.speed_mm.max(speed);
                }
                // The writers refuse every other curve before they finish.
                _ => {}
            }
        }
        for face in &definition.faces {
            dec.charge_work(1)?;
            match &face.surface {
                Surface::Plane(plane) => {
                    levers.reach_mm = levers.reach_mm.max(point_size(&plane.origin));
                    // A STEP plane's u and v are millimetres along unit axes.
                    levers.speed_mm = levers.speed_mm.max(1.0);
                }
                Surface::Cylinder(cylinder) => {
                    let radius = size(&cylinder.radius_sq).sqrt();
                    levers.reach_mm = levers
                        .reach_mm
                        .max(point_size(&cylinder.origin))
                        .max(radius);
                    // u is the angle; v is millimetres along the axis.
                    levers.speed_mm = levers.speed_mm.max(radius).max(1.0);
                }
                Surface::Cone(cone) => {
                    levers.reach_mm = levers.reach_mm.max(point_size(&cone.apex));
                    // v is the height along the axis, so a point moves
                    // sqrt(1 + slope^2) per unit of v. Its speed in u is its
                    // distance from the axis, taken below from the reach.
                    let slope = size(&cone.slope);
                    levers.speed_mm = levers.speed_mm.max((1.0 + slope * slope).sqrt());
                }
                Surface::Sphere(sphere) => {
                    let radius = size(&sphere.radius);
                    levers.reach_mm = levers.reach_mm.max(point_size(&sphere.center)).max(radius);
                    levers.speed_mm = levers.speed_mm.max(radius);
                }
                Surface::Torus(torus) => {
                    let radii = size(&torus.major) + size(&torus.minor);
                    levers.reach_mm = levers.reach_mm.max(point_size(&torus.center)).max(radii);
                    levers.speed_mm = levers.speed_mm.max(radii);
                }
                Surface::BSpline(support) => levers.take_surface(support.carrier(), dec)?,
                // The writers refuse every other surface before they finish.
                _ => {}
            }
        }
        // No point of the body is farther from a cone's axis than two points
        // of the cube of half side M are from each other.
        if definition
            .faces
            .iter()
            .any(|face| matches!(face.surface, Surface::Cone(_)))
        {
            levers.speed_mm = levers.speed_mm.max(2.0 * 3f64.sqrt() * levers.reach_mm);
        }
        // Face copies: the u/v curves a file carries. A plane has one only
        // where the body retains a u/v curve on it.
        let mut native = std::collections::BTreeMap::new();
        for trim in native_trims {
            dec.charge_work(1)?;
            native.insert((trim.face.as_str(), trim.boundary, trim.use_index), trim);
        }
        for face in &definition.faces {
            for (boundary, uses) in face.loops.iter().enumerate() {
                for (use_index, trim) in uses.iter().enumerate() {
                    dec.charge_work(1)?;
                    let retained: Option<&&NumericalNativeTrim> =
                        native.get(&(face.key.as_str(), boundary, use_index));
                    if matches!(face.surface, Surface::Plane(_)) && retained.is_none() {
                        continue;
                    }
                    for end in &trim.uv {
                        for value in end {
                            levers.face_copy_reach =
                                levers.face_copy_reach.max(value.to_f64_nearest().abs());
                        }
                    }
                    if let Some(retained) = retained {
                        for control in retained.curve.controls() {
                            levers.face_copy_reach = levers
                                .face_copy_reach
                                .max(size(&control.x))
                                .max(size(&control.y));
                        }
                        let speed = levers.take_curve(&retained.curve, dec)?;
                        levers.face_copy_speed = levers.face_copy_speed.max(speed);
                    } else if let Some(edge) = definition.edges.get(trim.edge) {
                        let width = (&edge.parameters[1] - &edge.parameters[0])
                            .to_f64_nearest()
                            .abs();
                        if width > 0.0 {
                            let step = [
                                (&trim.uv[1][0] - &trim.uv[0][0]).to_f64_nearest(),
                                (&trim.uv[1][1] - &trim.uv[0][1]).to_f64_nearest(),
                                0.0,
                            ];
                            levers.face_copy_speed =
                                levers.face_copy_speed.max(length(step) / width);
                        }
                    }
                }
            }
        }
        Ok(levers)
    }
}

impl RoundingError {
    /// The bound for the largest errors of a written file and the levers of
    /// its body. Each term is one kind of written number.
    ///
    /// - **Lengths.** A point has three coordinates, so it moves by at most
    ///   `sqrt(3) e`; a torus adds its two radii. `(sqrt(3) + 2) e`.
    /// - **Line magnitudes.** A line's point at parameter `t` is its origin
    ///   plus `t` times the magnitude along the direction. `T e`.
    /// - **Directions.** A direction of three ratios, the largest of size
    ///   one, turns by at most `sqrt(3) e`; a placement has two, so its frame
    ///   turns by `2 sqrt(3) e`. No point is farther from a placement's
    ///   origin than `4 sqrt(3) M`: both lie within `M`, a radius and a
    ///   second radius of the cube of half side `M`. The product is `24 M e`.
    ///   A cone's semi-angle is one more turn of the same size or less.
    /// - **Parameters.** A knot moves a spline's point by at most the speed
    ///   over the degree, and a point depends on `2p + 2` knots (Piegl and
    ///   Tiller, "Computing the derivative of NURBS with respect to a knot",
    ///   Computer Aided Geometric Design 15, 1998): `4 V e` for a curve, and
    ///   `8 V e` for the two directions of a surface. A trim parameter moves
    ///   an edge's end by `V e`. A face copy's u/v point moves by its own
    ///   controls (`e`), its knots (`4 U e`), its trim parameters (`U e`),
    ///   or, for a u/v line, by its origin (`e`), its magnitude (`T e`) and
    ///   its direction (`sqrt(2) T U e`): `(1 + T + (5 + 2T) U) e` covers
    ///   both. The surface carries that to space at `2 V`, one `V` for each
    ///   of u and v. The sum is `(15 + 2T + (10 + 4T) U) V e`.
    /// - **Weights.** A weight moves a rational point by its share of the
    ///   distance to a control, at most `2 sqrt(3) M`, over the smallest
    ///   weight: once for the curve and once for the surface. A face copy's
    ///   weights move its u/v point by `2 Q e / w`, carried to space at
    ///   `2 V`. `(4 sqrt(3) M + 4 V Q) e / w`.
    pub(super) fn bound(errors: RoundingErrors, levers: RoundingLevers) -> Result<Self, Refusal> {
        let RoundingLevers {
            reach_mm: m,
            parameter: t,
            speed_mm: v,
            face_copy_speed: u,
            face_copy_reach: q,
            least_weight: w,
        } = levers;
        let root3 = 3f64.sqrt();
        let terms = [
            (root3 + 2.0, &errors.length_mm),
            (t, &errors.line_magnitude),
            (24.0 * m, &errors.direction),
            (
                (15.0 + 2.0 * t + (10.0 + 4.0 * t) * u) * v,
                &errors.parameter,
            ),
            ((4.0 * root3 * m + 4.0 * v * q) / w, &errors.weight),
        ];
        let mut bound_mm = Rational::zero();
        for (lever, error) in terms {
            if error.is_zero() {
                continue;
            }
            let lever = Rational::from_f64_exact(lever).ok_or_else(|| {
                bad(
                    "unsupported",
                    "the body's size has no finite binary64 value for its rounding error",
                )
            })?;
            bound_mm = bound_mm + lever * error.clone();
        }
        Ok(Self {
            errors,
            levers,
            bound_mm,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(text: &str) -> Rational {
        Rational::parse(text).unwrap()
    }

    #[test]
    fn a_file_with_no_rounding_has_no_rounding_error() {
        let levers = RoundingLevers {
            reach_mm: 100.0,
            parameter: 6.5,
            speed_mm: 40.0,
            face_copy_speed: 2.0,
            face_copy_reach: 7.0,
            least_weight: 0.5,
        };
        let exact = RoundingError::bound(RoundingErrors::default(), levers).unwrap();
        assert_eq!(exact.bound_mm, Rational::zero());
    }

    #[test]
    fn each_kind_of_number_adds_its_error_times_its_lever() {
        let levers = RoundingLevers {
            reach_mm: 100.0,
            parameter: 2.0,
            speed_mm: 8.0,
            face_copy_speed: 0.5,
            face_copy_reach: 4.0,
            least_weight: 0.5,
        };
        let only = |errors: RoundingErrors| {
            RoundingError::bound(errors, levers)
                .unwrap()
                .bound_mm
                .to_f64_nearest()
        };
        let e = q("0.001");
        let length = only(RoundingErrors {
            length_mm: e.clone(),
            ..Default::default()
        });
        assert!((length - (3f64.sqrt() + 2.0) * 1e-3).abs() < 1e-15);
        let magnitude = only(RoundingErrors {
            line_magnitude: e.clone(),
            ..Default::default()
        });
        assert!((magnitude - 2.0e-3).abs() < 1e-15);
        let direction = only(RoundingErrors {
            direction: e.clone(),
            ..Default::default()
        });
        assert!((direction - 2.4).abs() < 1e-12);
        // (15 + 2 T + (10 + 4 T) U) V = (15 + 4 + 18 * 0.5) * 8 = 224.
        let parameter = only(RoundingErrors {
            parameter: e.clone(),
            ..Default::default()
        });
        assert!((parameter - 0.224).abs() < 1e-12);
        // (4 sqrt(3) M + 4 V Q) / w = (400 sqrt(3) + 128) / 0.5.
        let weight = only(RoundingErrors {
            weight: e,
            ..Default::default()
        });
        assert!((weight - (800.0 * 3f64.sqrt() + 256.0) * 1e-3).abs() < 1e-12);
    }

    #[test]
    fn a_spline_is_no_faster_than_its_control_steps_over_their_knot_spans() {
        // A cubic Bezier from 0 to 9 on [0, 1]: steps of 3, speed 3 * 3 = 9.
        let knots: Vec<Rational> = [0, 0, 0, 0, 1, 1, 1, 1].map(Rational::from_i64).to_vec();
        let controls = [0.0, 3.0, 6.0, 9.0].map(|x| [x, 0.0, 0.0]);
        let (speed, least) = spline_speed(3, &knots, &controls, 1, &[1.0; 4]);
        assert_eq!(speed.to_bits(), 9.0_f64.to_bits());
        assert_eq!(least.to_bits(), 1.0_f64.to_bits());
        // Weights of 1 and 2 allow four times that.
        let (rational, least) = spline_speed(3, &knots, &controls, 1, &[1.0, 2.0, 2.0, 1.0]);
        assert_eq!(rational.to_bits(), 36.0_f64.to_bits());
        assert_eq!(least.to_bits(), 1.0_f64.to_bits());
    }
}
