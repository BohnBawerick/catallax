//! Sample the authored law before binary64 proposal conversion.
//! Coordinates are relative to the authored spine start, including comparisons.

use super::{Frame, Path, SweepPath, SweepRecipe, V, affine_row, exact_controls, q, rational};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal};
use catallax_recipe::Point3 as RecipePoint;
use catallax_shapes::{BSplineCurve, Point3, Vector3};

const BITS: u32 = 128;
type Bounds = [Ball; 3];

fn scalar(value: &Rational) -> Ball {
    Ball::from_rational(value, BITS)
}
fn integer(value: i64) -> Ball {
    scalar(&q(value))
}
fn sum(a: &Bounds, b: &Bounds) -> Bounds {
    std::array::from_fn(|i| &a[i] + &b[i])
}
fn difference(a: &Bounds, b: &Bounds) -> Bounds {
    std::array::from_fn(|i| &a[i] - &b[i])
}
fn scale(a: &Bounds, b: &Ball) -> Bounds {
    std::array::from_fn(|i| &a[i] * b)
}
fn product(a: &Bounds, b: &Bounds) -> Ball {
    let x = &a[0] * &b[0];
    let y = &a[1] * &b[1];
    let z = &a[2] * &b[2];
    &(&x + &y) + &z
}
fn crossed(a: &Bounds, b: &Bounds) -> Bounds {
    [
        &(&a[1] * &b[2]) - &(&a[2] * &b[1]),
        &(&a[2] * &b[0]) - &(&a[0] * &b[2]),
        &(&a[0] * &b[1]) - &(&a[1] * &b[0]),
    ]
}
fn squared(a: &Bounds) -> Ball {
    a.iter()
        .fold(integer(0), |sum, value| &sum + &value.square())
}
fn precision() -> Refusal {
    Refusal::new(
        "sweep.motion-precision",
        "the authored motion cannot be resolved at the working precision",
    )
}
fn unit_bounds(a: &Bounds) -> Result<Bounds, Refusal> {
    let length = squared(a).sqrt().ok_or_else(precision)?;
    let inverse = length.recip().ok_or_else(precision)?;
    Ok(scale(a, &inverse))
}
fn enclose(value: &Vector3, dec: &mut Decider) -> Result<Bounds, Refusal> {
    Ok([
        dec.enclose_budgeted_at(&value.x, BITS)?,
        dec.enclose_budgeted_at(&value.y, BITS)?,
        dec.enclose_budgeted_at(&value.z, BITS)?,
    ])
}
fn normalized(value: &Vector3, dec: &mut Decider) -> Result<Bounds, Refusal> {
    let length = value.norm_squared().sqrt();
    let length = dec.enclose_budgeted_at(&length, BITS)?;
    let inverse = length.recip().ok_or_else(precision)?;
    Ok(scale(&enclose(value, dec)?, &inverse))
}
fn authored_point(value: &RecipePoint) -> Point3 {
    let [x, y, z] = value
        .0
        .each_ref()
        .map(|value| Real::from(value.value().clone()));
    Point3::new(x, y, z)
}
fn start(path: &SweepPath) -> Point3 {
    match path {
        SweepPath::Line { start, .. } | SweepPath::Helix { start, .. } => authored_point(start),
        SweepPath::Spline { control_points, .. } => authored_point(&control_points[0]),
    }
}
fn sincos(angle: &Ball, dec: &mut Decider) -> Result<(Ball, Ball), Refusal> {
    dec.charge_work(1)?;
    let result = (angle.sin(), angle.cos());
    dec.check_deadline()?;
    Ok(result)
}

struct Line {
    start: Vector3,
    delta: Vector3,
    tangent: Bounds,
}
struct Helix {
    start: Bounds,
    radial: Bounds,
    around: Bounds,
    axis: Bounds,
    unit_radial: Bounds,
    unit_around: Bounds,
    omega: Ball,
    rise: Ball,
    inverse_speed: Ball,
    spin: Ball,
}
struct Location {
    origin: Bounds,
    tangent: Option<Bounds>,
    radial: Option<Bounds>,
}
// Only the matching cache is populated. Keeping caches inline avoids a heap
// allocation for each source path solely to balance enum variant sizes.
struct SourcePath<'a> {
    line: Option<Line>,
    helix: Option<Helix>,
    spline: Option<&'a BSplineCurve>,
    anchor: Point3,
}
impl<'a> SourcePath<'a> {
    fn new(
        raw: &SweepPath,
        path: &'a Path,
        anchor: &Point3,
        pi: Option<&Ball>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut result = Self {
            line: None,
            helix: None,
            spline: None,
            anchor: anchor.clone(),
        };
        match (raw, path) {
            (SweepPath::Line { start, end }, _) => {
                let start = authored_point(start);
                let delta = authored_point(end).minus(&start);
                result.line = Some(Line {
                    start: start.minus(anchor),
                    tangent: normalized(&delta, dec)?,
                    delta,
                });
            }
            (
                SweepPath::Helix {
                    start,
                    centre,
                    axis,
                    angle,
                    rise,
                },
                _,
            ) => {
                let start = authored_point(start);
                let radial = start.minus(&authored_point(centre));
                let axis = Vector3::new(
                    axis.0[0].value().clone().into(),
                    axis.0[1].value().clone().into(),
                    axis.0[2].value().clone().into(),
                );
                let axis = axis.scaled(&(&Real::one() / &axis.norm_squared().sqrt()));
                let around = axis.cross(&radial);
                let omega = pi.expect("a helix requests pi") * &scalar(&(angle.value() / &q(180)));
                let rise = scalar(rise.value());
                let radius_squared = dec.enclose_budgeted_at(&radial.norm_squared(), BITS)?;
                let speed_squared = &(&radius_squared * &omega.square()) + &rise.square();
                let inverse_speed = speed_squared
                    .sqrt()
                    .and_then(|speed| speed.recip())
                    .ok_or_else(precision)?;
                let spin = -&(&(&omega * &rise) * &inverse_speed);
                result.helix = Some(Helix {
                    start: enclose(&start.minus(anchor), dec)?,
                    unit_radial: normalized(&radial, dec)?,
                    unit_around: normalized(&around, dec)?,
                    radial: enclose(&radial, dec)?,
                    around: enclose(&around, dec)?,
                    axis: enclose(&axis, dec)?,
                    omega,
                    rise,
                    inverse_speed,
                    spin,
                });
            }
            (SweepPath::Spline { .. }, Path::Spline(curve)) => result.spline = Some(curve),
            _ => unreachable!("validated path and source kind agree"),
        }
        Ok(result)
    }
    fn at(&self, t: &Rational, with_frame: bool, dec: &mut Decider) -> Result<Location, Refusal> {
        dec.charge_work(1)?;
        if let Some(line) = &self.line {
            let value = line.start.plus(&line.delta.scaled(&Real::from(t.clone())));
            return Ok(Location {
                origin: enclose(&value, dec)?,
                tangent: with_frame.then(|| line.tangent.clone()),
                radial: None,
            });
        }
        if let Some(helix) = &self.helix {
            let parameter = scalar(t);
            let (s, c) = sincos(&(&helix.omega * &parameter), dec)?;
            let origin = sum(
                &helix.start,
                &sum(
                    &sum(
                        &scale(&helix.radial, &(&c - &integer(1))),
                        &scale(&helix.around, &s),
                    ),
                    &scale(&helix.axis, &(&helix.rise * &parameter)),
                ),
            );
            let tangent = with_frame.then(|| {
                scale(
                    &sum(
                        &sum(
                            &scale(&helix.radial, &-&(&helix.omega * &s)),
                            &scale(&helix.around, &(&helix.omega * &c)),
                        ),
                        &scale(&helix.axis, &helix.rise),
                    ),
                    &helix.inverse_speed,
                )
            });
            let radial = with_frame.then(|| {
                sum(
                    &scale(&helix.unit_radial, &c),
                    &scale(&helix.unit_around, &s),
                )
            });
            return Ok(Location {
                origin,
                tangent,
                radial,
            });
        }
        let curve = self.spline.expect("one validated source path cache");
        let [lo, hi] = curve.domain();
        let parameter = Real::from(lo + &(&(hi - lo) * t));
        let point = curve.point_at(&parameter, dec)?;
        let tangent = if with_frame {
            Some(normalized(&curve.first_at(&parameter, dec)?, dec)?)
        } else {
            None
        };
        Ok(Location {
            origin: enclose(&point.minus(&self.anchor), dec)?,
            tangent,
            radial: None,
        })
    }
}

pub(super) struct FrameSample {
    origin: Bounds,
    x: Bounds,
    y: Bounds,
    scale: Ball,
}
pub(super) struct SourceLaw<'a> {
    path: SourcePath<'a>,
    guide: Option<SourcePath<'a>>,
    anchor: Point3,
    section_axes: Option<[Vector3; 2]>,
    reference: Bounds,
    first: [Bounds; 3],
    initial: [Bounds; 2],
    normal_offset: Ball,
    twist: Option<Ball>,
    end_scale: &'a Rational,
    fixed: bool,
}
impl<'a> SourceLaw<'a> {
    pub(super) fn new(
        recipe: &'a SweepRecipe,
        path: &'a Path,
        guide: Option<&'a Path>,
        profile: &crate::Profile,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let needs_pi = matches!(recipe.path, SweepPath::Helix { .. })
            || matches!(recipe.guide, Some(SweepPath::Helix { .. }))
            || !recipe.twist_degrees.value().is_zero();
        let pi = needs_pi.then(|| Ball::pi_checked(BITS, dec)).transpose()?;
        let anchor = start(&recipe.path);
        let path = SourcePath::new(&recipe.path, path, &anchor, pi.as_ref(), dec)?;
        let guide = recipe
            .guide
            .as_ref()
            .zip(guide)
            .map(|(raw, path)| SourcePath::new(raw, path, &anchor, pi.as_ref(), dec))
            .transpose()?;
        let section_axes = recipe.plane.as_ref().map(|_| {
            [
                profile
                    .pose
                    .rotation
                    .apply(&Vector3::new(1.into(), 0.into(), 0.into())),
                profile
                    .pose
                    .rotation
                    .apply(&Vector3::new(0.into(), 1.into(), 0.into())),
            ]
        });
        let normal_offset = if let Some([x, y]) = &section_axes {
            let origin = Point3::new(
                profile.pose.translation.x.clone(),
                profile.pose.translation.y.clone(),
                profile.pose.translation.z.clone(),
            );
            dec.enclose_budgeted_at(&origin.minus(&anchor).dot(&x.cross(y)), BITS)?
        } else {
            integer(0)
        };
        let twist = if recipe.twist_degrees.value().is_zero() {
            None
        } else {
            Some(
                pi.as_ref().expect("a twist requests pi")
                    * &scalar(&(recipe.twist_degrees.value() / &q(180))),
            )
        };
        let reference = if let Some([x, _]) = &section_axes {
            enclose(x, dec)?
        } else {
            enclose(
                &Vector3::new(
                    recipe.reference[0].value().clone().into(),
                    recipe.reference[1].value().clone().into(),
                    recipe.reference[2].value().clone().into(),
                ),
                dec,
            )?
        };
        let first = Self::base(&path, guide.as_ref(), &q(0), &reference, dec)?.1;
        let initial = if let Some([x, y]) = &section_axes {
            [enclose(x, dec)?, enclose(y, dec)?]
        } else if guide.is_some() {
            [first[0].clone(), first[1].clone()]
        } else {
            let x = unit_bounds(&difference(
                &reference,
                &scale(&first[2], &product(&reference, &first[2])),
            ))?;
            let y = crossed(&first[2], &x);
            [x, y]
        };
        Ok(Self {
            path,
            guide,
            anchor,
            section_axes,
            reference,
            first,
            initial,
            normal_offset,
            twist,
            end_scale: recipe.end_scale.value(),
            fixed: recipe.orientation == catallax_recipe::SweepOrientation::Fixed
                && recipe.guide.is_none(),
        })
    }
    fn base(
        path: &SourcePath<'_>,
        guide: Option<&SourcePath<'_>>,
        t: &Rational,
        reference: &Bounds,
        dec: &mut Decider,
    ) -> Result<(Bounds, [Bounds; 3]), Refusal> {
        let location = path.at(t, true, dec)?;
        let tangent = location.tangent.expect("a source frame requests a tangent");
        if guide.is_none()
            && let Some(helix) = &path.helix
        {
            let radial = location
                .radial
                .expect("a helix frame retains its radial direction");
            let binormal = crossed(&tangent, &radial);
            let (s, c) = sincos(&(&helix.spin * &scalar(t)), dec)?;
            let x = sum(&scale(&radial, &c), &scale(&binormal, &s));
            let y = crossed(&tangent, &x);
            return Ok((location.origin, [x, y, tangent]));
        }
        let direction = if let Some(guide) = guide {
            difference(&guide.at(t, false, dec)?.origin, &location.origin)
        } else {
            reference.clone()
        };
        let x = unit_bounds(&difference(
            &direction,
            &scale(&tangent, &product(&direction, &tangent)),
        ))?;
        let y = crossed(&tangent, &x);
        Ok((location.origin, [x, y, tangent]))
    }
    pub(super) fn control(&self, point: &Point3, dec: &mut Decider) -> Result<[Ball; 2], Refusal> {
        let local = if let Some([x, y]) = &self.section_axes {
            let delta = point.minus(&self.anchor);
            [delta.dot(x), delta.dot(y)]
        } else {
            [point.x.clone(), point.y.clone()]
        };
        Ok([
            dec.enclose_budgeted_at(&local[0], BITS)?,
            dec.enclose_budgeted_at(&local[1], BITS)?,
        ])
    }
    pub(super) fn at(&self, t: &Rational, dec: &mut Decider) -> Result<FrameSample, Refusal> {
        let (origin, current) = if self.fixed && self.twist.is_none() {
            (self.path.at(t, false, dec)?.origin, None)
        } else {
            let (origin, current) =
                Self::base(&self.path, self.guide.as_ref(), t, &self.reference, dec)?;
            (origin, Some(current))
        };
        let carry = |axis: &Bounds| {
            if self.fixed || t.is_zero() {
                return axis.clone();
            }
            let current = current.as_ref().expect("moving source frame");
            sum(
                &sum(
                    &scale(&current[0], &product(axis, &self.first[0])),
                    &scale(&current[1], &product(axis, &self.first[1])),
                ),
                &scale(&current[2], &product(axis, &self.first[2])),
            )
        };
        let mut x = carry(&self.initial[0]);
        let mut y = carry(&self.initial[1]);
        if let Some(twist) = &self.twist {
            let tangent = &current.as_ref().expect("twist source frame")[2];
            let (s, c) = sincos(&(twist * &scalar(t)), dec)?;
            let turn = |axis: &Bounds| {
                sum(
                    &sum(&scale(axis, &c), &scale(&crossed(tangent, axis), &s)),
                    &scale(tangent, &(&product(axis, tangent) * &(&integer(1) - &c))),
                )
            };
            x = turn(&x);
            y = turn(&y);
        }
        let scale_value = &q(1) + &(&(self.end_scale - &q(1)) * t);
        let scale_value = scalar(&scale_value);
        let offset = &self.normal_offset * &scale_value;
        let origin = if offset.is_exact_zero() {
            origin
        } else {
            sum(&origin, &scale(&unit_bounds(&crossed(&x, &y))?, &offset))
        };
        Ok(FrameSample {
            origin,
            x,
            y,
            scale: scale_value,
        })
    }
    pub(super) fn error(
        &self,
        frame: &FrameSample,
        control: &[Ball; 2],
        fitted: &Point3,
        dec: &mut Decider,
    ) -> Result<Rational, Refusal> {
        let expected = sum(
            &frame.origin,
            &scale(
                &sum(&scale(&frame.x, &control[0]), &scale(&frame.y, &control[1])),
                &frame.scale,
            ),
        );
        let fitted = enclose(&fitted.minus(&self.anchor), dec)?;
        Ok(squared(&difference(&fitted, &expected)).hi())
    }
}

pub(super) struct Control {
    pub fitted: V,
    pub authored: [Ball; 2],
}
pub(super) enum Fit {
    Line([Point3; 2]),
    Cubic([Point3; 4]),
}
impl Fit {
    pub(super) fn new(frames: [Frame; 4], point: V, linear: bool) -> Result<Self, Refusal> {
        if linear {
            let make = |frame| {
                affine_row(frame, &[rational(point[0])?, rational(point[1])?])
                    .map(|[x, y, z]| Point3::new(x.into(), y.into(), z.into()))
            };
            Ok(Self::Line([make(frames[0])?, make(frames[3])?]))
        } else {
            Ok(Self::Cubic(exact_controls(frames, point)?))
        }
    }
    pub(super) fn at(&self, t: &Rational) -> Point3 {
        let s = &q(1) - t;
        let linear;
        let cubic;
        let (points, weights): (&[Point3], &[Rational]) = match self {
            Self::Line(points) => {
                linear = [s, t.clone()];
                (points, &linear)
            }
            Self::Cubic(points) => {
                let ss = &s * &s;
                let tt = t * t;
                cubic = [&ss * &s, &q(3) * &(&ss * t), &q(3) * &(&s * &tt), &tt * t];
                (points, &cubic)
            }
        };
        let coordinate = |get: fn(&Point3) -> &Real| {
            points
                .iter()
                .zip(weights)
                .fold(q(0), |sum, (point, weight)| {
                    sum + get(point).rational_ref().expect("fitted affine row") * weight
                })
        };
        Point3::new(
            coordinate(|p| &p.x).into(),
            coordinate(|p| &p.y).into(),
            coordinate(|p| &p.z).into(),
        )
    }
}

pub(super) fn error_bound(upper: &Rational) -> Result<f64, Refusal> {
    let error = if upper.is_zero() {
        0.0
    } else {
        libm::sqrt(upper.to_f64_nearest()).next_up().next_up()
    };
    let rounded = rational(error)?;
    if !error.is_finite() || &rounded * &rounded < *upper {
        return Err(precision());
    }
    Ok(error)
}
