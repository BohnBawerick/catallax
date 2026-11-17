//! Sampled source-carrier fits. These types confer no exact-model authority.

mod family_seeds;

use family_seeds::FamilySeeds;

use catallax_numbers::{Decider, Rational, Refusal};
use serde::Serialize;

use crate::geometry_source::{
    AngleUnit, CoordinateSpace, EntityId, MappingError, ParameterChart, ParsedStep,
};

use super::Error;

const SAMPLE_BITS: u32 = 96;
const NEWTON_STEPS: usize = 24;

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

#[derive(Debug)]
struct Carrier {
    source: EntityId,
    charts: Vec<ParameterChart>,
    dimension: usize,
    period: f64,
}

fn broken(detail: &str) -> Error {
    Refusal::new("import.tolerant-broken-geometry", detail).into()
}

fn mapped(error: MappingError) -> Error {
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
    (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2])
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
        let point = chart
            .evaluate_box(
                &bounds[..self.dimension],
                SAMPLE_BITS.min(dec.policy().max_precision_bits),
                dec,
            )
            .map_err(mapped)?;
        let point_mm = point.each_ref().map(|value| value.mid().to_f64_nearest());
        let mut uncertainty_mm: f64 = 0.0;
        for (axis, value) in point.iter().enumerate() {
            let width = (&value.hi() - &value.lo()).to_f64_nearest().abs();
            uncertainty_mm =
                uncertainty_mm.hypot(width + 8.0 * f64::EPSILON * point_mm[axis].abs());
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
                if determinant <= 1e-24 * (a + d).powi(2) || !determinant.is_finite() {
                    break;
                }
                let first = dot(jacobian[0], residual);
                let second = dot(jacobian[1], residual);
                delta = [
                    (d * first - b * second) / determinant,
                    (a * second - b * first) / determinant,
                ];
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
    /// Finite-difference normal direction at one sample. `None` means the
    /// sampled frame collapsed; it is not a whole-face regularity theorem.
    pub fn sampled_normal(
        &self,
        location: NativeLocation,
        dec: &mut Decider,
    ) -> Result<Option<[f64; 3]>, Error> {
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
        self.0.project(point_mm, seed, self.2.as_ref(), dec)
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
        self.0.location([parameter, 0.0])
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
            } else if self.1 {
                Ok(if same_sense {
                    [endpoints[0], endpoints[0] + self.0.period]
                } else {
                    [endpoints[0] + self.0.period, endpoints[0]]
                })
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
        self.0.sample(location, dec)
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
            _ => return Ok(false),
        }
    }
    Err(broken("cyclic or excessive curve wrappers"))
}

/// A source PCURVE. Its values are UV coordinates, never Cartesian millimetres.
#[derive(Debug)]
pub struct NativeUvCurve(Carrier);
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
        Ok(Self(Carrier::new(
            entity,
            charts,
            1,
            CoordinateSpace::SurfaceParameters,
            &angles,
        )?))
    }
    /// Evaluate UV; the uncertainty is in those same parameter units.
    pub fn sample(
        &self,
        location: NativeLocation,
        dec: &mut Decider,
    ) -> Result<([f64; 2], f64), Error> {
        let value = self.0.sample(location, dec)?;
        Ok(([value.point_mm[0], value.point_mm[1]], value.uncertainty_mm))
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
