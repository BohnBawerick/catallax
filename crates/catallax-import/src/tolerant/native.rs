//! Sampled source-carrier fits. These types confer no exact-model authority.

mod family_seeds;

use family_seeds::FamilySeeds;

use catallax_numbers::{Decider, Rational, Refusal};
use libm::hypot;
use serde::Serialize;

use crate::geometry_source::{
    AngleUnit, Binary64Chart, CoordinateSpace, EntityId, MappingError, ParameterChart, ParsedStep,
};

use super::Error;

const SAMPLE_BITS: u32 = 96;
const NEWTON_STEPS: usize = 24;
/// The error a display sample may carry, from rounding and from a reduced
/// periodic parameter together: a tenth of the width projection stops at. A
/// sample binary64 cannot place this well is enclosed.
const DISPLAY_ROUNDING_MM: f64 = 1e-8;

/// A location in one retained native chart, in original source parameters.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct NativeLocation {
    /// Index in the retained carrier's chart cover.
    pub chart: usize,
    /// The second coordinate is zero for a physical curve.
    pub parameters: [f64; 2],
}

/// One evaluated physical sample, not a whole-domain deviation bound.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct NativeSample {
    /// Cartesian position in millimetres.
    pub point_mm: [f64; 3],
    /// Evaluation enclosure and floating conversion uncertainty at this sample.
    pub uncertainty_mm: f64,
}

/// A sampled fit. The residual is not asserted to be the nearest distance.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct NativeProjection {
    /// Chosen source parameters and chart.
    pub location: NativeLocation,
    /// Actual evaluation at the chosen parameters.
    pub sample: NativeSample,
    /// Observed distance plus this evaluation's uncertainty, in millimetres.
    pub residual_mm: f64,
    /// Refinement iterations on the selected path.
    pub iterations: usize,
    /// Whether the sampled residual reached the numerical stopping width.
    pub converged: bool,
}

#[derive(Clone, Debug)]
struct Carrier {
    source: EntityId,
    charts: Vec<ParameterChart>,
    dimension: usize,
    period: f64,
    /// The charts compiled for binary64, in a carrier made for display.
    display: Option<Vec<Binary64Chart>>,
}

fn broken(detail: &str) -> Error {
    Refusal::new("import.tolerant-broken-geometry", detail).into()
}

pub(super) fn mapped(error: MappingError) -> Error {
    match error {
        MappingError::Source(error) => error.into(),
        MappingError::Geometry(error) => error.into(),
    }
}

pub(super) fn control(error: &Error) -> bool {
    match error {
        Error::Source(crate::ImportRefusal::ResourceLimit) => true,
        Error::Numerical(error) => {
            error.code.contains("budget")
                || error.code.contains("deadline")
                || error.code.contains("precision-limit")
        }
        _ => false,
    }
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    hypot(hypot(a[0] - b[0], a[1] - b[1]), a[2] - b[2])
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn difference(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|axis| a[axis] - b[axis])
}

fn same_location(a: &NativeLocation, b: &NativeLocation) -> bool {
    a.chart == b.chart && a.parameters.map(f64::to_bits) == b.parameters.map(f64::to_bits)
}

impl Carrier {
    fn new(
        source: EntityId,
        charts: Vec<ParameterChart>,
        dimension: usize,
        space: CoordinateSpace,
        angles: &AngleUnit,
    ) -> Result<Self, Error> {
        if charts.is_empty()
            || charts
                .iter()
                .any(|chart| chart.dimension() != dimension || chart.coordinate_space() != space)
        {
            return Err(broken(
                "source carrier has inconsistent dimension or coordinate units",
            ));
        }
        let period = match angles {
            AngleUnit::Radians => std::f64::consts::TAU,
            AngleUnit::Degrees => 360.0,
            AngleUnit::ScaledRadians(factor) => {
                std::f64::consts::TAU / factor.to_f64_nearest().abs()
            }
        };
        if !period.is_finite() || period <= 0.0 {
            return Err(broken("source angle unit is not finite and positive"));
        }
        Ok(Self {
            source,
            charts,
            dimension,
            period,
            display: None,
        })
    }

    /// This carrier with its samples evaluated in binary64 where that
    /// resolves them. Only [`DisplaySurface`] and [`DisplayCurve`] hold one.
    fn for_display(&self, dec: &mut Decider) -> Result<Self, Error> {
        let display = self
            .charts
            .iter()
            .map(|chart| chart.binary64(dec).map_err(mapped))
            .collect::<Result<_, _>>()?;
        Ok(Self {
            display: Some(display),
            ..self.clone()
        })
    }

    fn location(&self, parameters: [f64; 2]) -> Result<NativeLocation, Error> {
        if parameters[..self.dimension]
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(broken("native parameters are not finite"));
        }
        let chart = self
            .charts
            .iter()
            .position(|chart| {
                (0..self.dimension).all(|axis| {
                    let bounds = &chart.domain()[axis];
                    bounds[0]
                        .as_ref()
                        .is_none_or(|low| parameters[axis] >= low.to_f64_nearest())
                        && bounds[1]
                            .as_ref()
                            .is_none_or(|high| parameters[axis] <= high.to_f64_nearest())
                })
            })
            .ok_or_else(|| broken("parameters are outside every native source chart"))?;
        Ok(NativeLocation { chart, parameters })
    }

    fn periodic_location(
        &self,
        parameters: [f64; 2],
        periods: [Option<f64>; 2],
    ) -> Result<NativeLocation, Error> {
        if let Ok(location) = self.location(parameters) {
            return Ok(location);
        }
        let canonical = std::array::from_fn(|axis| {
            periods[axis].map_or(parameters[axis], |period| {
                parameters[axis].rem_euclid(period)
            })
        });
        let mut location = self.location(canonical)?;
        location.parameters = parameters;
        Ok(location)
    }

    fn periodic_sample(
        &self,
        location: NativeLocation,
        periods: [Option<f64>; 2],
        dec: &mut Decider,
    ) -> Result<NativeSample, Error> {
        if let Ok(canonical) = self.location(location.parameters) {
            return self.sample(canonical, dec);
        }
        let mut ranges = [[0.0; 2]; 2];
        for axis in 0..self.dimension {
            let parameter = location.parameters[axis];
            if !parameter.is_finite() {
                return Err(broken("periodic parameter is not finite"));
            }
            ranges[axis] = if let Some(period) = periods[axis] {
                let radius = 16.0 * f64::EPSILON * (parameter.abs() + period);
                if !radius.is_finite() || radius >= period / 8.0 {
                    return Err(broken("periodic parameter reduction cannot be resolved"));
                }
                let canonical = parameter.rem_euclid(period);
                [canonical - radius, canonical + radius]
            } else {
                [parameter; 2]
            };
        }
        let mut result: Option<NativeSample> = None;
        for (index, chart) in self.charts.iter().enumerate() {
            for x in -1..=1 {
                if periods[0].is_none() && x != 0 {
                    continue;
                }
                for y in -1..=1 {
                    if (self.dimension == 1 || periods[1].is_none()) && y != 0 {
                        continue;
                    }
                    let offsets = [x, y];
                    let mut bounds: [[Rational; 2]; 2] =
                        std::array::from_fn(|_| [Rational::zero(), Rational::zero()]);
                    let mut intersects = true;
                    for axis in 0..self.dimension {
                        let offset = periods[axis].unwrap_or(0.0) * f64::from(offsets[axis]);
                        let low = Rational::from_f64_exact(ranges[axis][0] + offset)
                            .ok_or_else(|| broken("periodic bound is not finite"))?;
                        let high = Rational::from_f64_exact(ranges[axis][1] + offset)
                            .ok_or_else(|| broken("periodic bound is not finite"))?;
                        let domain = &chart.domain()[axis];
                        bounds[axis] = [
                            domain[0]
                                .as_ref()
                                .map_or(low.clone(), |limit| low.max(limit.clone())),
                            domain[1]
                                .as_ref()
                                .map_or(high.clone(), |limit| high.min(limit.clone())),
                        ];
                        intersects &= bounds[axis][0] <= bounds[axis][1];
                    }
                    if !intersects {
                        continue;
                    }
                    dec.charge_work(1)?;
                    let sample = self.sample_region(index, &bounds[..self.dimension], dec)?;
                    match &mut result {
                        Some(total) => {
                            total.uncertainty_mm = total.uncertainty_mm.max(
                                distance(total.point_mm, sample.point_mm) + sample.uncertainty_mm,
                            )
                        }
                        None => result = Some(sample),
                    }
                }
            }
        }
        result.ok_or_else(|| broken("periodic parameters are outside the source chart cover"))
    }

    fn sample(&self, location: NativeLocation, dec: &mut Decider) -> Result<NativeSample, Error> {
        dec.check_deadline()?;
        dec.charge_work(1)?;
        let chart = self
            .charts
            .get(location.chart)
            .ok_or_else(|| broken("native chart index is out of range"))?;
        let mut bounds = std::array::from_fn::<_, 2, _>(|_| [Rational::zero(), Rational::zero()]);
        for (axis, bound) in bounds.iter_mut().enumerate().take(self.dimension) {
            let mut parameter = Rational::from_f64_exact(location.parameters[axis])
                .ok_or_else(|| broken("native parameter is not finite"))?;
            if let Some(low) = &chart.domain()[axis][0]
                && parameter < *low
            {
                parameter = low.clone();
            }
            if let Some(high) = &chart.domain()[axis][1]
                && parameter > *high
            {
                parameter = high.clone();
            }
            *bound = [parameter.clone(), parameter];
        }
        self.sample_region(location.chart, &bounds[..self.dimension], dec)
    }

    fn sample_region(
        &self,
        chart: usize,
        bounds: &[[Rational; 2]],
        dec: &mut Decider,
    ) -> Result<NativeSample, Error> {
        if let Some(display) = &self.display {
            // The middle of the box, which the callers have already kept
            // inside the chart's domain. A box wider than a point is a
            // periodic parameter's reduction, known only to that width: its
            // radius goes into the evaluation, so the bound says how far the
            // reduction can move the point.
            let mut parameters = [0.0; 2];
            let mut radii = [0.0; 2];
            for (axis, range) in bounds.iter().enumerate() {
                let middle = (&range[0] + &range[1]).mul_pow2(-1);
                let parameter = middle.to_f64_nearest();
                let half = (&range[1] - &range[0]).mul_pow2(-1).to_f64_nearest();
                let rounded = Rational::from_f64_exact(parameter) != Some(middle);
                parameters[axis] = parameter;
                radii[axis] = half * (1.0 + 4.0 * f64::EPSILON)
                    + if rounded {
                        f64::EPSILON * parameter.abs()
                    } else {
                        0.0
                    };
            }
            dec.charge_work(display[chart].work())?;
            let (point_mm, bound) = display[chart].evaluate(&parameters, &radii);
            // Past this bound binary64 has lost the carrier to cancellation,
            // or the reduction leaves the point too loosely known, or the
            // point is singular: the enclosing evaluation below answers, at
            // its price.
            if bound <= DISPLAY_ROUNDING_MM {
                let size = point_mm
                    .iter()
                    .fold(0.0f64, |size, value| size.max(value.abs()));
                return Ok(NativeSample {
                    point_mm,
                    uncertainty_mm: bound + 8.0 * f64::EPSILON * size,
                });
            }
        }
        let chart = &self.charts[chart];
        let point = chart
            .evaluate_box(
                bounds,
                SAMPLE_BITS.min(dec.policy().max_precision_bits),
                dec,
            )
            .map_err(mapped)?;
        let point_mm = point.each_ref().map(|value| value.mid().to_f64_nearest());
        let mut uncertainty_mm: f64 = 0.0;
        for (axis, value) in point.iter().enumerate() {
            let width = (&value.hi() - &value.lo()).to_f64_nearest().abs();
            uncertainty_mm = hypot(
                uncertainty_mm,
                width + 8.0 * f64::EPSILON * point_mm[axis].abs(),
            );
        }
        if point_mm.iter().any(|value| !value.is_finite()) || !uncertainty_mm.is_finite() {
            return Err(broken("source carrier sample has non-finite coordinates"));
        }
        Ok(NativeSample {
            point_mm,
            uncertainty_mm,
        })
    }

    fn try_sample(
        &self,
        location: NativeLocation,
        dec: &mut Decider,
    ) -> Result<Option<NativeSample>, Error> {
        match self.sample(location, dec) {
            Ok(sample) => Ok(Some(sample)),
            Err(error) if control(&error) => Err(error),
            // A chart seed can meet a coordinate singularity outside the finite trim.
            Err(_) => Ok(None),
        }
    }

    fn clamp(&self, mut location: NativeLocation) -> NativeLocation {
        for axis in 0..self.dimension {
            let limits = &self.charts[location.chart].domain()[axis];
            if let Some(low) = &limits[0] {
                location.parameters[axis] = location.parameters[axis].max(low.to_f64_nearest());
            }
            if let Some(high) = &limits[1] {
                location.parameters[axis] = location.parameters[axis].min(high.to_f64_nearest());
            }
        }
        location
    }

    fn jacobian(
        &self,
        location: NativeLocation,
        dec: &mut Decider,
    ) -> Result<Option<[[f64; 3]; 2]>, Error> {
        let mut derivatives = [[0.0; 3]; 2];
        for (axis, derivative) in derivatives.iter_mut().enumerate().take(self.dimension) {
            let limits = &self.charts[location.chart].domain()[axis];
            let width = match (&limits[0], &limits[1]) {
                (Some(low), Some(high)) => (high - low).to_f64_nearest().abs(),
                _ => 1.0 + location.parameters[axis].abs(),
            };
            let step = (width * 1e-5).max(1e-8 * (1.0 + location.parameters[axis].abs()));
            let mut left = location;
            let mut right = location;
            left.parameters[axis] -= step;
            right.parameters[axis] += step;
            left = self.clamp(left);
            right = self.clamp(right);
            let span = right.parameters[axis] - left.parameters[axis];
            if span == 0.0 || !span.is_finite() {
                return Ok(None);
            }
            let (Some(a), Some(b)) = (self.try_sample(left, dec)?, self.try_sample(right, dec)?)
            else {
                return Ok(None);
            };
            *derivative = std::array::from_fn(|coordinate| {
                (b.point_mm[coordinate] - a.point_mm[coordinate]) / span
            });
        }
        Ok(Some(derivatives))
    }

    fn refine(
        &self,
        target: [f64; 3],
        location: NativeLocation,
        dec: &mut Decider,
    ) -> Result<Option<NativeProjection>, Error> {
        let mut location = self.clamp(location);
        let Some(sample) = self.try_sample(location, dec)? else {
            return Ok(None);
        };
        let mut best = NativeProjection {
            location,
            residual_mm: distance(sample.point_mm, target) + sample.uncertainty_mm,
            sample,
            iterations: 0,
            converged: false,
        };
        let stopping =
            1e-7 + 32.0 * f64::EPSILON * target.iter().map(|x| x.abs()).fold(1.0, f64::max);
        for iteration in 0..NEWTON_STEPS {
            if best.residual_mm <= stopping {
                best.converged = true;
                break;
            }
            let Some(jacobian) = self.jacobian(best.location, dec)? else {
                break;
            };
            let residual = difference(target, best.sample.point_mm);
            let a = dot(jacobian[0], jacobian[0]);
            let mut delta = [0.0; 2];
            if self.dimension == 1 {
                if a <= 0.0 || !a.is_finite() {
                    break;
                }
                delta[0] = dot(jacobian[0], residual) / a;
            } else {
                let b = dot(jacobian[0], jacobian[1]);
                let d = dot(jacobian[1], jacobian[1]);
                let determinant = a * d - b * b;
                if !determinant.is_finite() {
                    break;
                }
                if determinant <= 1e-24 * (a + d).powi(2) {
                    let axis = usize::from(d > a);
                    let squared_speed = [a, d][axis];
                    if squared_speed <= 0.0 || !squared_speed.is_finite() {
                        break;
                    }
                    delta[axis] = dot(jacobian[axis], residual) / squared_speed;
                } else {
                    let first = dot(jacobian[0], residual);
                    let second = dot(jacobian[1], residual);
                    delta = [
                        (d * first - b * second) / determinant,
                        (a * second - b * first) / determinant,
                    ];
                }
            }
            let mut improved = false;
            for reduction in 0..8 {
                let scale = 2.0f64.powi(-reduction);
                location = best.location;
                for axis in 0..self.dimension {
                    location.parameters[axis] += scale * delta[axis];
                }
                location = self.clamp(location);
                // The retained best sample already evaluates these exact bits;
                // evaluating it again cannot satisfy the strict improvement test.
                if same_location(&location, &best.location) {
                    dec.check_deadline()?;
                    continue;
                }
                if let Some(sample) = self.try_sample(location, dec)? {
                    let residual_mm = distance(sample.point_mm, target) + sample.uncertainty_mm;
                    if residual_mm < best.residual_mm {
                        best = NativeProjection {
                            location,
                            sample,
                            residual_mm,
                            iterations: iteration + 1,
                            converged: residual_mm <= stopping,
                        };
                        improved = true;
                        break;
                    }
                }
            }
            if !improved {
                break;
            }
        }
        Ok(Some(best))
    }

    fn project(
        &self,
        target: [f64; 3],
        seed: Option<NativeLocation>,
        family: Option<&FamilySeeds>,
        dec: &mut Decider,
    ) -> Result<NativeProjection, Error> {
        if target.iter().any(|value| !value.is_finite()) {
            return Err(broken("projection target is not finite"));
        }
        let mut best = None;
        if let Some(seed) = seed.filter(|seed| seed.chart < self.charts.len())
            && let Some(candidate) = self.refine(target, seed, dec)?
        {
            if candidate.converged {
                return Ok(candidate);
            }
            best = Some(candidate);
        }
        if let Some(family) = family {
            let point = target.map(Rational::from_f64_exact);
            if let [Some(x), Some(y), Some(z)] = point {
                for parameters in family.seeds(&[x, y, z], dec).map_err(mapped)? {
                    let parameters = parameters.each_ref().map(Rational::to_f64_nearest);
                    let Ok(location) = self.location(parameters) else {
                        continue;
                    };
                    if let Some(candidate) = self.refine(target, location, dec)? {
                        if candidate.converged {
                            return Ok(candidate);
                        }
                        if best
                            .as_ref()
                            .is_none_or(|old| candidate.residual_mm < old.residual_mm)
                        {
                            best = Some(candidate);
                        }
                    }
                }
            }
        }
        for (chart, data) in self.charts.iter().enumerate() {
            let center = std::array::from_fn(|axis| {
                match (&data.domain()[axis][0], &data.domain()[axis][1]) {
                    (Some(low), Some(high)) => ((low + high).mul_pow2(-1)).to_f64_nearest(),
                    _ => 0.0,
                }
            });
            // Finite limits already determine the scan ranges. Origin samples
            // and derivatives contribute only to an unbounded search envelope.
            let (base, derivatives) = if data
                .domain()
                .iter()
                .take(self.dimension)
                .any(|bounds| bounds.iter().any(Option::is_none))
            {
                let origin = NativeLocation {
                    chart,
                    parameters: center,
                };
                (self.try_sample(origin, dec)?, self.jacobian(origin, dec)?)
            } else {
                (None, None)
            };
            let mut ranges = [[0.0; 2]; 2];
            for axis in 0..self.dimension {
                ranges[axis] = match (&data.domain()[axis][0], &data.domain()[axis][1]) {
                    (Some(low), Some(high)) => [low.to_f64_nearest(), high.to_f64_nearest()],
                    _ => {
                        let physical = base.map_or(1.0, |sample| distance(sample.point_mm, target));
                        let speed = derivatives
                            .map_or(1.0, |values| dot(values[axis], values[axis]).sqrt());
                        let half = ((physical + 1.0) / speed.max(1e-12)).max(self.period / 2.0);
                        [center[axis] - half, center[axis] + half]
                    }
                };
            }
            let mut seed_best: Option<(NativeLocation, f64)> = None;
            for first in 0..5 {
                for second in 0..if self.dimension == 1 { 1 } else { 5 } {
                    let mut parameters = center;
                    for (axis, station) in
                        [first, second].into_iter().enumerate().take(self.dimension)
                    {
                        let fraction = if data.domain()[axis].iter().all(Option::is_some) {
                            station as f64 / 4.0
                        } else {
                            0.1 + station as f64 / 5.0
                        };
                        parameters[axis] =
                            ranges[axis][0] + fraction * (ranges[axis][1] - ranges[axis][0]);
                    }
                    let location = NativeLocation { chart, parameters };
                    if let Some(sample) = self.try_sample(location, dec)? {
                        let residual = distance(sample.point_mm, target) + sample.uncertainty_mm;
                        if seed_best.as_ref().is_none_or(|(_, old)| residual < *old) {
                            seed_best = Some((location, residual));
                        }
                    }
                }
            }
            if let Some((seed, _)) = seed_best
                && let Some(candidate) = self.refine(target, seed, dec)?
                && best
                    .as_ref()
                    .is_none_or(|old: &NativeProjection| candidate.residual_mm < old.residual_mm)
            {
                best = Some(candidate);
            }
        }
        best.ok_or_else(|| broken("no finite sample on any retained native chart"))
    }
}

/// A physical source surface with sampled UV fitting, not an embedding proof.
#[derive(Debug)]
pub struct NativeSurface(Carrier, [Option<f64>; 2], Option<FamilySeeds>);
impl NativeSurface {
    /// Compile source charts without any whole-domain topology admission.
    pub fn from_source(
        source: &ParsedStep,
        entity: EntityId,
        length_mm: &Rational,
        angles: AngleUnit,
        dec: &mut Decider,
    ) -> Result<Self, Error> {
        let charts = source
            .surface_charts(entity, length_mm, angles.clone(), dec)
            .map_err(mapped)?;
        let carrier = Carrier::new(entity, charts, 2, CoordinateSpace::Millimetres, &angles)?;
        let periodic = surface_periodic_axes(source, entity, dec)?;
        let periods = periodic.map(|periodic| periodic.then_some(carrier.period));
        let family =
            FamilySeeds::from_source(source, entity, length_mm, &angles, dec).map_err(mapped)?;
        Ok(Self(carrier, periods, family))
    }
    /// This surface for display. See [`DisplaySurface`].
    pub fn for_display(&self, dec: &mut Decider) -> Result<DisplaySurface, Error> {
        Ok(DisplaySurface(Self(
            self.0.for_display(dec)?,
            self.1,
            self.2.clone(),
        )))
    }
    /// Original source surface identity.
    pub fn source(&self) -> EntityId {
        self.0.source
    }
    /// Actual native source chart count.
    pub fn chart_count(&self) -> usize {
        self.0.charts.len()
    }
    /// Original exact limits; an absent bound is not an inferred trim.
    pub fn domain(&self, chart: usize) -> &[[Option<Rational>; 2]; 2] {
        self.0.charts[chart].domain()
    }
    /// Angular period in the source's plane-angle unit.
    pub fn angular_period(&self) -> f64 {
        self.0.period
    }
    /// Known angular periods in original source UV units. A closed spline
    /// declaration alone does not authorize a periodic parameter extension.
    pub fn periods(&self) -> [Option<f64>; 2] {
        self.1
    }
    /// Keep equivalent periodic UV on one sheet. The first location uses the
    /// native carrier's nonnegative angular branch.
    pub fn align_parameters(
        &self,
        mut parameters: [f64; 2],
        reference: Option<[f64; 2]>,
    ) -> [f64; 2] {
        for axis in 0..2 {
            if let Some(period) = self.1[axis] {
                parameters[axis] = match reference {
                    Some(reference) => {
                        parameters[axis]
                            + ((reference[axis] - parameters[axis]) / period).round() * period
                    }
                    None => parameters[axis].rem_euclid(period),
                };
            }
        }
        parameters
    }
    /// Resolve original source UV to its native knot span.
    pub fn location(&self, parameters: [f64; 2]) -> Result<NativeLocation, Error> {
        self.0.periodic_location(parameters, self.1)
    }
    /// One physical sample in original source UV.
    pub fn sample(
        &self,
        location: NativeLocation,
        dec: &mut Decider,
    ) -> Result<NativeSample, Error> {
        self.0.periodic_sample(location, self.1, dec)
    }
    /// Finite-difference normal direction at one sample. `None` means the
    /// sampled frame collapsed; it is not a whole-face regularity theorem.
    pub fn sampled_normal(
        &self,
        location: NativeLocation,
        dec: &mut Decider,
    ) -> Result<Option<[f64; 3]>, Error> {
        let parameters = std::array::from_fn(|axis| {
            self.1[axis].map_or(location.parameters[axis], |period| {
                location.parameters[axis].rem_euclid(period)
            })
        });
        let location = self.0.location(parameters)?;
        let Some([mut a, mut b]) = self.0.jacobian(location, dec)? else {
            return Ok(None);
        };
        let scale = |vector: &[f64; 3]| vector.iter().map(|value| value.abs()).fold(0.0, f64::max);
        let first = scale(&a);
        let second = scale(&b);
        if first == 0.0 || second == 0.0 || !first.is_finite() || !second.is_finite() {
            return Ok(None);
        }
        a = a.map(|value| value / first);
        b = b.map(|value| value / second);
        let normal = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        Ok(normal.iter().any(|value| *value != 0.0).then_some(normal))
    }
    /// Fit source UV and report the actual sampled residual and convergence.
    pub fn project(
        &self,
        point_mm: [f64; 3],
        seed: Option<NativeLocation>,
        dec: &mut Decider,
    ) -> Result<NativeProjection, Error> {
        let canonical = seed
            .map(|location| {
                let parameters = std::array::from_fn(|axis| {
                    self.1[axis].map_or(location.parameters[axis], |period| {
                        location.parameters[axis].rem_euclid(period)
                    })
                });
                self.0.location(parameters)
            })
            .transpose()?;
        let mut projection = self.0.project(point_mm, canonical, self.2.as_ref(), dec)?;
        if let Some(seed) = seed {
            let parameters =
                self.align_parameters(projection.location.parameters, Some(seed.parameters));
            projection.location = self.location(parameters)?;
            projection.sample = self.sample(projection.location, dec)?;
            projection.residual_mm =
                distance(projection.sample.point_mm, point_mm) + projection.sample.uncertainty_mm;
        }
        Ok(projection)
    }
}

/// A [`NativeSurface`] for display: the same charts, periods and seeds, with
/// every sample, normal and projection evaluated in binary64.
///
/// A display takes thousands of samples per face and states no bound on any
/// of them. An enclosing evaluation of each costs a body's whole work
/// allowance, so these samples are rounded evaluations that carry a running
/// rounding bound, and only a sample binary64 cannot place within 1e-8 mm is
/// enclosed. It is a separate type so that nothing which fits, admits or
/// certifies geometry can be handed one. The edge-gap checker has a binary64
/// reading of its own, [`Binary64Carrier`], whose bounds go into what it
/// reports.
#[derive(Debug)]
pub struct DisplaySurface(NativeSurface);
impl DisplaySurface {
    /// Known angular periods in original source UV units.
    pub fn periods(&self) -> [Option<f64>; 2] {
        self.0.periods()
    }
    /// Resolve original source UV to its native knot span.
    pub fn location(&self, parameters: [f64; 2]) -> Result<NativeLocation, Error> {
        self.0.location(parameters)
    }
    /// One physical sample in original source UV.
    pub fn sample(
        &self,
        location: NativeLocation,
        dec: &mut Decider,
    ) -> Result<NativeSample, Error> {
        self.0.sample(location, dec)
    }
    /// Finite-difference normal direction at one sample, or `None` where the
    /// sampled frame collapsed.
    pub fn sampled_normal(
        &self,
        location: NativeLocation,
        dec: &mut Decider,
    ) -> Result<Option<[f64; 3]>, Error> {
        self.0.sampled_normal(location, dec)
    }
    /// Fit source UV to a point and report the sampled residual.
    pub fn project(
        &self,
        point_mm: [f64; 3],
        seed: Option<NativeLocation>,
        dec: &mut Decider,
    ) -> Result<NativeProjection, Error> {
        self.0.project(point_mm, seed, dec)
    }
}

/// A [`NativeCurve`] for display, sampled as a [`DisplaySurface`] is.
#[derive(Debug)]
pub struct DisplayCurve(NativeCurve);
impl DisplayCurve {
    /// Resolve an original source parameter to its native knot span.
    pub fn location(&self, parameter: f64) -> Result<NativeLocation, Error> {
        self.0.location(parameter)
    }
    /// One physical curve sample.
    pub fn sample(
        &self,
        location: NativeLocation,
        dec: &mut Decider,
    ) -> Result<NativeSample, Error> {
        self.0.sample(location, dec)
    }
}

/// A physical source curve, independent of any finite-edge certificate.
#[derive(Debug)]
pub struct NativeCurve(Carrier, bool);
impl NativeCurve {
    /// Compile physical source curves without monotonicity or embedding gates.
    pub fn from_source(
        source: &ParsedStep,
        entity: EntityId,
        length_mm: &Rational,
        angles: AngleUnit,
        dec: &mut Decider,
    ) -> Result<Self, Error> {
        let charts = source
            .curve_charts(entity, length_mm, angles.clone(), dec)
            .map_err(mapped)?;
        let angular = angular_curve(source, entity, dec)?;
        Ok(Self(
            Carrier::new(entity, charts, 1, CoordinateSpace::Millimetres, &angles)?,
            angular,
        ))
    }
    /// This curve for display. See [`DisplayCurve`].
    pub fn for_display(&self, dec: &mut Decider) -> Result<DisplayCurve, Error> {
        Ok(DisplayCurve(Self(self.0.for_display(dec)?, self.1)))
    }
    /// Original curve identity.
    pub fn source(&self) -> EntityId {
        self.0.source
    }
    /// Native source chart count.
    pub fn chart_count(&self) -> usize {
        self.0.charts.len()
    }
    /// Original exact limits, not reconstructed finite-edge parameters.
    pub fn domain(&self, chart: usize) -> &[[Option<Rational>; 2]; 2] {
        self.0.charts[chart].domain()
    }
    /// Resolve an original source parameter to its native knot span.
    pub fn location(&self, parameter: f64) -> Result<NativeLocation, Error> {
        self.0
            .periodic_location([parameter, 0.0], [self.1.then_some(self.0.period), None])
    }
    /// Source angular period for periodic analytic curves.
    pub fn angular_period(&self) -> f64 {
        self.0.period
    }
    /// Resolve a finite source EDGE_CURVE traversal from fitted endpoint
    /// parameters. Closed analytic curves use the source angular period;
    /// closed splines use their declared finite native domain.
    pub fn finite_range(
        &self,
        mut endpoints: [f64; 2],
        same_sense: bool,
        closed: bool,
    ) -> Result<[f64; 2], Error> {
        if self.1 {
            endpoints = endpoints.map(|parameter| parameter.rem_euclid(self.0.period));
        }
        if closed {
            if self.1 {
                return Ok(if same_sense {
                    [endpoints[0], endpoints[0] + self.0.period]
                } else {
                    [endpoints[0] + self.0.period, endpoints[0]]
                });
            }
            let mut domain: Option<[f64; 2]> = None;
            for chart in &self.0.charts {
                if let [Some(low), Some(high)] = &chart.domain()[0] {
                    let low = low.to_f64_nearest();
                    let high = high.to_f64_nearest();
                    domain = Some(match domain {
                        None => [low, high],
                        Some([a, b]) => [a.min(low), b.max(high)],
                    });
                }
            }
            return if let Some([low, high]) = domain {
                Ok(if same_sense { [low, high] } else { [high, low] })
            } else {
                Err(broken(
                    "closed edge has no finite or periodic source curve domain",
                ))
            };
        }
        if self.1 {
            let period = self.0.period;
            if same_sense {
                endpoints[1] += ((endpoints[0] - endpoints[1]) / period).floor().max(0.0) * period;
                if endpoints[1] <= endpoints[0] {
                    endpoints[1] += period;
                }
            } else {
                endpoints[1] -= ((endpoints[1] - endpoints[0]) / period).floor().max(0.0) * period;
                if endpoints[1] >= endpoints[0] {
                    endpoints[1] -= period;
                }
                if endpoints[1] < 0.0 {
                    endpoints[0] += period;
                    endpoints[1] += period;
                }
            }
        }
        Ok(endpoints)
    }
    /// One physical curve sample.
    pub fn sample(
        &self,
        location: NativeLocation,
        dec: &mut Decider,
    ) -> Result<NativeSample, Error> {
        self.0
            .periodic_sample(location, [self.1.then_some(self.0.period), None], dec)
    }
    /// Fit a source parameter and retain the observed residual.
    pub fn project(
        &self,
        point_mm: [f64; 3],
        seed: Option<NativeLocation>,
        dec: &mut Decider,
    ) -> Result<NativeProjection, Error> {
        self.0.project(point_mm, seed, None, dec)
    }
}

fn surface_periodic_axes(
    source: &ParsedStep,
    mut entity: EntityId,
    dec: &mut Decider,
) -> Result<[bool; 2], Error> {
    use crate::geometry_source::{Geometry, Surface};
    for _ in 0..32 {
        match source.geometry(entity, dec)? {
            Geometry::Surface(Surface::Cylinder(..) | Surface::Cone(..) | Surface::Sphere(..)) => {
                return Ok([true, false]);
            }
            Geometry::Surface(Surface::Torus(..)) => return Ok([true; 2]),
            Geometry::Surface(Surface::Extrusion(curve, _)) => {
                return Ok([angular_curve(source, curve, dec)?, false]);
            }
            Geometry::Surface(Surface::Revolution(curve, _)) => {
                return Ok([true, angular_curve(source, curve, dec)?]);
            }
            Geometry::Surface(
                Surface::Replica(parent, _)
                | Surface::Offset(parent, ..)
                | Surface::CurveBounded(parent, ..),
            ) => entity = parent,
            _ => return Ok([false; 2]),
        }
    }
    Err(broken("cyclic or excessive surface wrappers"))
}

fn angular_curve(
    source: &ParsedStep,
    mut entity: EntityId,
    dec: &mut Decider,
) -> Result<bool, Error> {
    use crate::geometry_source::{Curve, Geometry};
    for _ in 0..32 {
        match source.geometry(entity, dec)? {
            Geometry::Curve(Curve::Circle(..) | Curve::Ellipse(..)) => return Ok(true),
            Geometry::Curve(
                Curve::Replica(parent, _)
                | Curve::Offset3d(parent, ..)
                | Curve::SurfaceCurve(parent, ..),
            ) => entity = parent,
            Geometry::Curve(Curve::Pcurve(_, representation)) => {
                let representation = source.representation(representation)?;
                if representation.items.len() != 1 {
                    return Ok(false);
                }
                entity = representation.items[0];
            }
            _ => return Ok(false),
        }
    }
    Err(broken("cyclic or excessive curve wrappers"))
}

/// A source PCURVE. Its values are UV coordinates, never Cartesian millimetres.
#[derive(Debug)]
pub struct NativeUvCurve(Carrier, bool);
impl NativeUvCurve {
    /// Compile a retained PCURVE in its original mixed surface-parameter units.
    pub fn from_source(
        source: &ParsedStep,
        entity: EntityId,
        length_mm: &Rational,
        angles: AngleUnit,
        dec: &mut Decider,
    ) -> Result<Self, Error> {
        let charts = source
            .curve_charts(entity, length_mm, angles.clone(), dec)
            .map_err(mapped)?;
        let angular = angular_curve(source, entity, dec)?;
        Ok(Self(
            Carrier::new(
                entity,
                charts,
                1,
                CoordinateSpace::SurfaceParameters,
                &angles,
            )?,
            angular,
        ))
    }
    /// Evaluate UV; the uncertainty is in those same parameter units.
    pub fn sample(
        &self,
        location: NativeLocation,
        dec: &mut Decider,
    ) -> Result<([f64; 2], f64), Error> {
        let value =
            self.0
                .periodic_sample(location, [self.1.then_some(self.0.period), None], dec)?;
        Ok(([value.point_mm[0], value.point_mm[1]], value.uncertainty_mm))
    }

    /// The u/v at the parameter of this curve nearest to `parameter`, and
    /// how far that parameter is from it: zero inside the curve's domain.
    /// A STEP edge has no parameter range of its own, so an end placed on
    /// its 3D curve from its vertex can fall just outside the interval the
    /// file gives the u/v curve.
    pub(crate) fn sample_nearest(
        &self,
        parameter: f64,
        dec: &mut Decider,
    ) -> Result<([f64; 2], f64), Error> {
        let at = |chart: usize, parameter: f64| NativeLocation {
            chart,
            parameters: [parameter, 0.0],
        };
        match self.sample(at(0, parameter), dec) {
            Ok((uv, _)) => return Ok((uv, 0.0)),
            Err(error) if control(&error) => return Err(error),
            Err(_) => {}
        }
        let mut nearest: Option<(usize, f64, f64)> = None;
        for (chart, cover) in self.0.charts.iter().enumerate() {
            let limits = &cover.domain()[0];
            let mut inside = parameter;
            if let Some(low) = &limits[0] {
                inside = inside.max(low.to_f64_nearest());
            }
            if let Some(high) = &limits[1] {
                inside = inside.min(high.to_f64_nearest());
            }
            let moved = (inside - parameter).abs();
            if nearest.is_none_or(|(_, _, least)| moved < least) {
                nearest = Some((chart, inside, moved));
            }
        }
        let (chart, inside, moved) =
            nearest.ok_or_else(|| broken("a u/v curve has no parameter chart"))?;
        let (uv, _) = self.sample(at(chart, inside), dec)?;
        Ok((uv, moved))
    }

    pub(crate) fn hint(
        &self,
        parameter: f64,
        dec: &mut Decider,
    ) -> Result<Option<[f64; 2]>, Error> {
        let Ok(location) = self.0.location([parameter, 0.0]) else {
            return Ok(None);
        };
        Ok(self
            .0
            .try_sample(location, dec)?
            .map(|sample| [sample.point_mm[0], sample.point_mm[1]]))
    }
}

/// A carrier's charts in binary64, for the edge-gap checker.
///
/// The checker of roadmap item 3.0A reads both copies of an edge in binary64
/// and puts the rounding bound of its samples into the measurement error it
/// reports (Checkpoint 2 decisions, section 4.2 step 7: "Record a noise floor
/// from the f64 error estimate instead of escalating precision"). This is
/// that reading of a native carrier: every point comes with the bound its
/// chart's binary64 program states. A point binary64 cannot place has an
/// infinite bound, and the checker's reading of its edge use is then not
/// finite, so the caller measures that use the enclosing way.
///
/// It is a type of its own, like [`DisplaySurface`]: the fits of a body's
/// own samples still enclose, and cannot be handed one.
#[derive(Debug)]
pub(crate) struct Binary64Carrier {
    charts: Vec<crate::geometry_source::Binary64Chart>,
    /// Each chart's limits, rounded to binary64.
    domains: Vec<[[Option<f64>; 2]; 2]>,
    dimension: usize,
    periods: [Option<f64>; 2],
}

impl Binary64Carrier {
    fn new(carrier: &Carrier, periods: [Option<f64>; 2], dec: &mut Decider) -> Result<Self, Error> {
        let charts = carrier
            .charts
            .iter()
            .map(|chart| chart.binary64(dec).map_err(mapped))
            .collect::<Result<_, _>>()?;
        let domains = carrier
            .charts
            .iter()
            .map(|chart| {
                std::array::from_fn(|axis| {
                    chart.domain()[axis]
                        .each_ref()
                        .map(|limit| limit.as_ref().map(Rational::to_f64_nearest))
                })
            })
            .collect();
        Ok(Self {
            charts,
            domains,
            dimension: carrier.dimension,
            periods,
        })
    }

    /// The chart whose domain holds `parameters`.
    fn holding(&self, parameters: [f64; 2]) -> Option<usize> {
        self.domains.iter().position(|domain| {
            (0..self.dimension).all(|axis| {
                domain[axis][0].is_none_or(|low| parameters[axis] >= low)
                    && domain[axis][1].is_none_or(|high| parameters[axis] <= high)
            })
        })
    }

    /// `parameters` moved inside one chart, and how far they moved.
    fn clamped(&self, chart: usize, mut parameters: [f64; 2]) -> ([f64; 2], f64) {
        let mut moved: f64 = 0.0;
        for axis in 0..self.dimension {
            let [low, high] = self.domains[chart][axis];
            let inside = parameters[axis]
                .max(low.unwrap_or(f64::NEG_INFINITY))
                .min(high.unwrap_or(f64::INFINITY));
            moved = moved.max((inside - parameters[axis]).abs());
            parameters[axis] = inside;
        }
        (parameters, moved)
    }

    /// The point at `parameters` and a bound on its rounding error. A
    /// periodic parameter outside the charts is reduced by its period, and
    /// the bound carries what the rounding of that reduction moves the point
    /// by, as [`Carrier::periodic_sample`] states it for an enclosure. Past
    /// the charts of a parameter with no period the nearest border answers.
    /// The bound is infinite where binary64 cannot place the point.
    pub(crate) fn evaluate(&self, parameters: [f64; 2]) -> ([f64; 3], f64) {
        let unplaced = ([0.0; 3], f64::INFINITY);
        if parameters[..self.dimension]
            .iter()
            .any(|value| !value.is_finite())
        {
            return unplaced;
        }
        let mut at = parameters;
        let mut radii = [0.0; 2];
        let mut chart = self.holding(at);
        if chart.is_none() {
            for axis in 0..self.dimension {
                if let Some(period) = self.periods[axis] {
                    radii[axis] = 16.0 * f64::EPSILON * (parameters[axis].abs() + period);
                    at[axis] = parameters[axis].rem_euclid(period);
                }
            }
            chart = self.holding(at);
        }
        let chart = chart.or_else(|| {
            (0..self.charts.len())
                .min_by(|a, b| self.clamped(*a, at).1.total_cmp(&self.clamped(*b, at).1))
        });
        let Some(chart) = chart else {
            return unplaced;
        };
        let (at, _) = self.clamped(chart, at);
        let (point, bound) =
            self.charts[chart].evaluate(&at[..self.dimension], &radii[..self.dimension]);
        if !bound.is_finite() || point.iter().any(|value| !value.is_finite()) {
            return unplaced;
        }
        let size = point
            .iter()
            .fold(0.0f64, |size, value| size.max(value.abs()));
        (point, bound + 8.0 * f64::EPSILON * size)
    }

    /// The work one evaluation counts: what a display sample charges, a
    /// sixty-fourth of a unit per binary64 step, for the costliest chart. The
    /// checker is told this before it evaluates, so it takes the work from
    /// its edge's allowance first, whichever chart a point falls in.
    pub(crate) fn work(&self) -> u64 {
        self.charts
            .iter()
            .map(crate::geometry_source::Binary64Chart::work)
            .max()
            .unwrap_or(0)
    }

    /// Every finite border of a chart along one parameter: the knots of a
    /// spline, the ends of a trimmed carrier.
    pub(crate) fn borders(&self, axis: usize) -> Vec<f64> {
        let mut borders: Vec<f64> = self
            .domains
            .iter()
            .flat_map(|domain| domain[axis])
            .flatten()
            .collect();
        borders.sort_by(f64::total_cmp);
        borders.dedup();
        borders
    }

    /// The least and greatest value the charts reach along one parameter.
    /// `None` on a side where a chart has no limit, and on both sides of a
    /// periodic parameter.
    pub(crate) fn extent(&self, axis: usize) -> [Option<f64>; 2] {
        if self.periods[axis].is_some() {
            return [None, None];
        }
        let side = |which: usize, pick: fn(f64, f64) -> f64| {
            self.domains
                .iter()
                .map(|domain| domain[axis][which])
                .try_fold(None, |best: Option<f64>, limit| {
                    limit.map(|limit| Some(best.map_or(limit, |best| pick(best, limit))))
                })
                .flatten()
        };
        [side(0, f64::min), side(1, f64::max)]
    }

    /// The period of each parameter, where it has one.
    pub(crate) fn periods(&self) -> [Option<f64>; 2] {
        self.periods
    }
}

impl NativeCurve {
    /// This curve for the edge-gap checker. See [`Binary64Carrier`].
    pub(crate) fn binary64(&self, dec: &mut Decider) -> Result<Binary64Carrier, Error> {
        Binary64Carrier::new(&self.0, [self.1.then_some(self.0.period), None], dec)
    }
}

impl NativeSurface {
    /// This surface for the edge-gap checker. See [`Binary64Carrier`].
    pub(crate) fn binary64(&self, dec: &mut Decider) -> Result<Binary64Carrier, Error> {
        Binary64Carrier::new(&self.0, self.1, dec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn projection_residual_has_platform_independent_bits() {
        // System hypot on Linux and macOS round these pairs differently.
        for (x, y, expected) in [
            (
                0x3ede_0504_6f1c_e3c4,
                0x4036_b87a_4d48_a1ba,
                0x4036_b87a_4d48_a2f8,
            ),
            (
                0x3fe0_9d39_31a5_8d35,
                0x4027_5b8f_67a1_6f11,
                0x4027_6177_4dbc_bb3e,
            ),
        ] {
            let a = [f64::from_bits(x), f64::from_bits(y), 0.0];
            assert_eq!(distance(a, [0.0; 3]).to_bits(), expected);
        }
    }

    // An open cubic about a metre long on the knot span [0, 1]. Sampled with
    // a period of 1 it is a carrier whose periodic axis has a finite domain,
    // which is where a parameter outside the domain is reduced to a box.
    fn carriers(closed: bool) -> (Carrier, Carrier) {
        let last = if closed { "0.,0.,0." } else { "1000.,0.,0." };
        let text = format!(
            "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('display'),'2;1');\nFILE_NAME('display.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n#1=CARTESIAN_POINT('',(0.,0.,0.));\n#2=CARTESIAN_POINT('',(400.,300.,0.));\n#3=CARTESIAN_POINT('',(800.,-300.,50.));\n#4=CARTESIAN_POINT('',({last}));\n#5=B_SPLINE_CURVE_WITH_KNOTS('',3,(#1,#2,#3,#4),.UNSPECIFIED.,.F.,.F.,(4,4),(0.,1.),.UNSPECIFIED.);\nENDSEC;\nEND-ISO-10303-21;\n"
        );
        let mut dec = Decider::default();
        let parsed = ParsedStep::read(text.as_bytes(), &mut dec).unwrap();
        let charts = parsed
            .curve_charts(
                5,
                &Rational::from_ratio(1, 1).unwrap(),
                AngleUnit::Radians,
                &mut dec,
            )
            .unwrap();
        let enclosing = Carrier::new(
            5,
            charts,
            1,
            CoordinateSpace::Millimetres,
            &AngleUnit::Radians,
        )
        .unwrap();
        let display = enclosing.for_display(&mut dec).unwrap();
        (enclosing, display)
    }

    fn wrapped(carrier: &Carrier, parameter: f64) -> NativeSample {
        let periods = [Some(1.0), None];
        let location = carrier
            .periodic_location([parameter, 0.0], periods)
            .unwrap();
        carrier
            .periodic_sample(location, periods, &mut Decider::default())
            .unwrap()
    }

    // The reduction of a wrapped parameter is known to 16 units of rounding
    // of the parameter and the period, and this curve never moves slower
    // than 600 mm per unit of parameter: the display sample has to say so.
    #[test]
    fn a_wrapped_display_sample_carries_the_reduction_of_its_parameter() {
        let (enclosing, display) = carriers(false);
        for parameter in [1000.25, -640.5, 77.125] {
            let enclosed = wrapped(&enclosing, parameter);
            let shown = wrapped(&display, parameter);
            let miss = distance(enclosed.point_mm, shown.point_mm);
            assert!(miss <= shown.uncertainty_mm, "{parameter}: {miss}");
            let reduction = 16.0 * f64::EPSILON * (parameter.abs() + 1.0);
            assert!(
                shown.uncertainty_mm >= 600.0 * reduction,
                "{parameter}: {} under {}",
                shown.uncertainty_mm,
                600.0 * reduction
            );
            assert!(shown.uncertainty_mm <= 4.0 * enclosed.uncertainty_mm);
        }
    }

    // A million periods out the reduction is known to 3.6e-9 of a parameter
    // unit, micrometres on this curve and past what a display sample may
    // carry, though binary64 rounds the point itself to 1e-12 mm.
    #[test]
    fn a_reduction_past_the_display_bound_is_enclosed() {
        let (enclosing, display) = carriers(false);
        for parameter in [1_000_000.25, -3_000_000.5] {
            let enclosed = wrapped(&enclosing, parameter);
            let shown = wrapped(&display, parameter);
            assert!(enclosed.uncertainty_mm > DISPLAY_ROUNDING_MM);
            assert_eq!(
                enclosed.point_mm.map(f64::to_bits),
                shown.point_mm.map(f64::to_bits)
            );
            assert_eq!(
                enclosed.uncertainty_mm.to_bits(),
                shown.uncertainty_mm.to_bits()
            );
        }
    }

    // A whole number of periods lands on the seam, where the reduced box is
    // cut into one piece at each end of the domain and the sample spans both.
    // On an open curve it is uncertain by the gap between the ends. On a
    // closed one the two middles are 1,084 reductions apart, and the piece at
    // the far end is 1,283 more wide: a sample that kept only the middles
    // would stop at the first figure.
    #[test]
    fn a_display_sample_on_the_seam_spans_both_ends_like_the_enclosure() {
        for closed in [true, false] {
            let (enclosing, display) = carriers(closed);
            for parameter in [512.0, -100.0, 4096.0] {
                let enclosed = wrapped(&enclosing, parameter);
                let shown = wrapped(&display, parameter);
                let miss = distance(enclosed.point_mm, shown.point_mm);
                assert!(miss <= shown.uncertainty_mm, "{parameter}: {miss}");
                assert!(shown.uncertainty_mm <= 4.0 * enclosed.uncertainty_mm);
                assert!(enclosed.uncertainty_mm <= 4.0 * shown.uncertainty_mm);
                if closed {
                    let reduction = 16.0 * f64::EPSILON * (parameter.abs() + 1.0);
                    assert!(
                        shown.uncertainty_mm >= 1700.0 * reduction,
                        "{parameter}: {} under {}",
                        shown.uncertainty_mm,
                        1700.0 * reduction
                    );
                } else {
                    assert!(shown.uncertainty_mm > 999.999, "{parameter}");
                }
            }
        }
    }
}
