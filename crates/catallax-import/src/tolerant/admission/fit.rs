use std::collections::BTreeMap;
use std::sync::Arc;

use catallax_numbers::{Decider, Rational};

use crate::geometry_source::{AngleUnit, Curve as SourceCurve, Geometry, ParsedStep, Topology};
use crate::{Curve, Placement, Surface};

use super::super::native::{
    NativeCurve, NativeLocation, NativeProjection, NativeSample, NativeSurface, NativeUvCurve,
};
use super::{
    EntityTolerance, Error, NativeTrim, NativeWarning, SampledCheck, TolerantBody, TolerantCounts,
    TolerantEvidence, UvSample, invalid,
};

mod crossing;
mod poles;
mod trace;

const STATIONS: usize = 16;

pub(super) struct Context<'a> {
    parsed: &'a ParsedStep,
    scale: &'a Rational,
    angles: AngleUnit,
    pub(super) curves: BTreeMap<u64, Arc<NativeCurve>>,
    surfaces: BTreeMap<u64, Arc<NativeSurface>>,
    pcurves: BTreeMap<u64, Arc<NativeUvCurve>>,
    repair_cap: f64,
    pub(super) stage: super::NativeStage,
}
#[derive(Debug)]
pub(super) struct Carriers {
    pub(super) curves: BTreeMap<u64, Arc<NativeCurve>>,
    pub(super) surfaces: BTreeMap<u64, Arc<NativeSurface>>,
    pub(super) angles: AngleUnit,
}
impl<'a> Context<'a> {
    pub fn new(
        parsed: &'a ParsedStep,
        scale: &'a Rational,
        angles: AngleUnit,
        curves: BTreeMap<u64, Arc<NativeCurve>>,
        repair_cap: f64,
    ) -> Self {
        Self {
            parsed,
            scale,
            angles,
            curves,
            surfaces: BTreeMap::new(),
            pcurves: BTreeMap::new(),
            repair_cap,
            stage: super::NativeStage::BodyBuild,
        }
    }
    pub fn into_carriers(self) -> Carriers {
        Carriers {
            curves: self.curves,
            surfaces: self.surfaces,
            angles: self.angles,
        }
    }
    fn curve(&mut self, entity: u64, dec: &mut Decider) -> Result<Arc<NativeCurve>, Error> {
        if let Some(curve) = self.curves.get(&entity) {
            return Ok(Arc::clone(curve));
        }
        let curve = Arc::new(NativeCurve::from_source(
            self.parsed,
            entity,
            self.scale,
            self.angles.clone(),
            dec,
        )?);
        self.curves.insert(entity, Arc::clone(&curve));
        Ok(curve)
    }
    fn surface(&mut self, entity: u64, dec: &mut Decider) -> Result<Arc<NativeSurface>, Error> {
        if let Some(surface) = self.surfaces.get(&entity) {
            return Ok(Arc::clone(surface));
        }
        let surface = Arc::new(NativeSurface::from_source(
            self.parsed,
            entity,
            self.scale,
            self.angles.clone(),
            dec,
        )?);
        self.surfaces.insert(entity, Arc::clone(&surface));
        Ok(surface)
    }
    fn pcurve(&mut self, entity: u64, dec: &mut Decider) -> Result<Arc<NativeUvCurve>, Error> {
        if let Some(curve) = self.pcurves.get(&entity) {
            return Ok(Arc::clone(curve));
        }
        let curve = Arc::new(NativeUvCurve::from_source(
            self.parsed,
            entity,
            self.scale,
            self.angles.clone(),
            dec,
        )?);
        self.pcurves.insert(entity, Arc::clone(&curve));
        Ok(curve)
    }
}

#[derive(Clone, Copy)]
struct Accumulator {
    tolerance: f64,
    residual: f64,
    samples: usize,
}
impl Accumulator {
    fn new(floor: f64) -> Self {
        Self {
            tolerance: floor,
            residual: 0.0,
            samples: 0,
        }
    }
    fn measure(&mut self, residual: f64) -> Result<(), Error> {
        if !residual.is_finite() || residual < 0.0 {
            return Err(invalid(
                "import.tolerant-broken-geometry",
                "non-finite local fit measurement",
            ));
        }
        self.residual = self.residual.max(residual);
        self.tolerance = self.tolerance.max(residual * (1.0 + 16.0 * f64::EPSILON));
        self.samples += 1;
        Ok(())
    }
    fn record(self, body: &str, kind: &'static str, entity: &str) -> EntityTolerance {
        EntityTolerance {
            body: body.into(),
            kind,
            entity: entity.into(),
            tolerance_mm: self.tolerance,
            maximum_sample_residual_mm: self.residual,
            samples: self.samples,
        }
    }
}

fn id(name: &str) -> Option<u64> {
    name.strip_prefix('#')?.parse().ok()
}
fn point(point: &crate::Vec3) -> Result<NativeSample, Error> {
    let coordinates = point.each_ref().map(Rational::to_f64_nearest);
    if coordinates.iter().any(|value| !value.is_finite()) {
        return Err(invalid(
            "import.tolerant-broken-geometry",
            "non-finite source vertex",
        ));
    }
    Ok(NativeSample {
        point_mm: coordinates,
        uncertainty_mm: coordinates
            .iter()
            .map(|x| 8.0 * f64::EPSILON * x.abs())
            .sum(),
    })
}
fn subtract(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|axis| a[axis] - b[axis])
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: [f64; 3]) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}
fn unit(a: [f64; 3]) -> Result<[f64; 3], Error> {
    let length = norm(a);
    if length <= 0.0 || !length.is_finite() {
        return Err(invalid(
            "import.tolerant-broken-geometry",
            "zero or non-finite native frame axis",
        ));
    }
    Ok(a.map(|x| x / length))
}

struct PlaneFrame {
    origin: [f64; 3],
    x: [f64; 3],
    y: [f64; 3],
    parameter_mm: f64,
}
impl PlaneFrame {
    fn new(placement: &Placement, parameter_mm: f64) -> Result<Self, Error> {
        if !parameter_mm.is_finite() || parameter_mm <= 0.0 {
            return Err(invalid(
                "import.tolerant-broken-geometry",
                "invalid plane parameter scale",
            ));
        }
        let z = unit(point(&placement.axis)?.point_mm)?;
        let reference = point(&placement.ref_dir)?.point_mm;
        let projection = dot(reference, z);
        let x = unit(std::array::from_fn(|i| reference[i] - projection * z[i]))?;
        let y = [
            z[1] * x[2] - z[2] * x[1],
            z[2] * x[0] - z[0] * x[2],
            z[0] * x[1] - z[1] * x[0],
        ];
        Ok(Self {
            origin: point(&placement.origin)?.point_mm,
            x,
            y,
            parameter_mm,
        })
    }
    fn project(&self, target: NativeSample) -> NativeProjection {
        let delta = subtract(target.point_mm, self.origin);
        let parameters = [dot(delta, self.x), dot(delta, self.y)];
        let position = std::array::from_fn(|i| {
            self.origin[i] + parameters[0] * self.x[i] + parameters[1] * self.y[i]
        });
        let uncertainty = 64.0
            * f64::EPSILON
            * (norm(self.origin)
                + norm(target.point_mm)
                + norm(position)
                + parameters[0].abs()
                + parameters[1].abs());
        NativeProjection {
            location: NativeLocation {
                chart: 0,
                parameters: parameters.map(|value| value / self.parameter_mm),
            },
            sample: NativeSample {
                point_mm: position,
                uncertainty_mm: uncertainty,
            },
            residual_mm: norm(subtract(position, target.point_mm))
                + target.uncertainty_mm
                + uncertainty,
            iterations: 1,
            converged: true,
        }
    }
}

enum PhysicalCurve {
    Native(Arc<NativeCurve>),
    Line {
        origin: [f64; 3],
        direction: [f64; 3],
    },
    Point([f64; 3]),
}
impl PhysicalCurve {
    fn sample_at(
        &self,
        parameter: f64,
        preferred: usize,
        dec: &mut Decider,
    ) -> Result<(NativeLocation, NativeSample), Error> {
        match self {
            Self::Native(curve) => {
                let location = NativeLocation {
                    chart: chart_at(curve, parameter, preferred),
                    parameters: [parameter, 0.0],
                };
                Ok((location, curve.sample(location, dec)?))
            }
            Self::Point(point) => Ok((
                NativeLocation {
                    chart: 0,
                    parameters: [parameter, 0.0],
                },
                NativeSample {
                    point_mm: *point,
                    uncertainty_mm: 8.0 * f64::EPSILON * norm(*point),
                },
            )),
            Self::Line { origin, direction } => {
                let position =
                    std::array::from_fn(|axis| origin[axis] + parameter * direction[axis]);
                Ok((
                    NativeLocation {
                        chart: 0,
                        parameters: [parameter, 0.0],
                    },
                    NativeSample {
                        point_mm: position,
                        uncertainty_mm: 8.0 * f64::EPSILON * norm(position),
                    },
                ))
            }
        }
    }
}
struct EdgeSamples {
    points: Vec<(NativeLocation, NativeSample)>,
    source_curve: Option<u64>,
    carrier: PhysicalCurve,
}
fn chart_at(curve: &NativeCurve, parameter: f64, preferred: usize) -> usize {
    let contains = |index| {
        let limits = &curve.domain(index)[0];
        limits[0]
            .as_ref()
            .is_none_or(|low| parameter >= low.to_f64_nearest())
            && limits[1]
                .as_ref()
                .is_none_or(|high| parameter <= high.to_f64_nearest())
    };
    if contains(preferred) {
        preferred
    } else {
        (0..curve.chart_count())
            .find(|&index| contains(index))
            .unwrap_or(preferred)
    }
}

fn finite_edge(
    context: &mut Context<'_>,
    body: &TolerantBody,
    edge_index: usize,
    vertices: &[NativeSample],
    vertex_tolerances: &mut [Accumulator],
    evidence: &mut TolerantEvidence,
    dec: &mut Decider,
) -> Result<(EdgeSamples, Accumulator), Error> {
    let edge = &body.edges[edge_index];
    let mut edge_tolerance = Accumulator::new(edge.tolerance_mm().to_f64_nearest());
    let source_curve = match id(&edge.name).filter(|_| !edge.degenerate) {
        Some(entity) => match context.parsed.topology(entity, dec)? {
            Topology::Edge { curve, .. } => curve,
            _ => None,
        },
        None => None,
    };
    let carrier = if edge.degenerate {
        PhysicalCurve::Point(vertices[edge.vertices[0]].point_mm)
    } else if let Curve::Line(origin, direction) = &edge.curve {
        if let Some(entity) = source_curve {
            context.curve(entity, dec)?;
        }
        PhysicalCurve::Line {
            origin: point(origin)?.point_mm,
            direction: point(direction)?.point_mm,
        }
    } else if let Some(entity) = source_curve {
        PhysicalCurve::Native(context.curve(entity, dec)?)
    } else {
        return Err(invalid(
            "import.tolerant-broken-reference",
            format!("{} has no source curve declaration", edge.name),
        ));
    };
    let mut endpoint_fits = Vec::with_capacity(2);
    for (endpoint, &vertex) in edge.vertices.iter().enumerate() {
        let target = vertices[vertex];
        let fit = match &carrier {
            PhysicalCurve::Native(curve) => curve.project(target.point_mm, None, dec)?,
            PhysicalCurve::Point(point) => NativeProjection {
                location: NativeLocation {
                    chart: 0,
                    parameters: [endpoint as f64, 0.0],
                },
                sample: NativeSample {
                    point_mm: *point,
                    uncertainty_mm: 8.0 * f64::EPSILON * norm(*point),
                },
                residual_mm: norm(subtract(*point, target.point_mm)),
                iterations: 0,
                converged: true,
            },
            PhysicalCurve::Line { origin, direction } => {
                let length = dot(*direction, *direction);
                if length <= 0.0 || !length.is_finite() {
                    return Err(invalid(
                        "import.tolerant-broken-geometry",
                        "zero source line direction",
                    ));
                }
                let parameter = dot(subtract(target.point_mm, *origin), *direction) / length;
                let evaluated = std::array::from_fn(|i| origin[i] + parameter * direction[i]);
                NativeProjection {
                    location: NativeLocation {
                        chart: 0,
                        parameters: [parameter, 0.0],
                    },
                    sample: NativeSample {
                        point_mm: evaluated,
                        uncertainty_mm: 8.0 * f64::EPSILON * norm(evaluated),
                    },
                    residual_mm: norm(subtract(evaluated, target.point_mm)),
                    iterations: 1,
                    converged: true,
                }
            }
        };
        let residual = fit.residual_mm + target.uncertainty_mm;
        vertex_tolerances[vertex].measure(residual)?;
        edge_tolerance.measure(residual)?;
        if !fit.converged {
            evidence.warnings.push(NativeWarning {
                body: Some(body.name.clone()),
                scope: edge.name.clone(),
                code: "endpoint-fit-not-converged",
                detail: format!(
                    "observed endpoint residual {residual:.17e} mm; not a nearest-distance proof"
                ),
            });
        }
        endpoint_fits.push(fit);
    }
    let mut parameters = [
        endpoint_fits[0].location.parameters[0],
        endpoint_fits[1].location.parameters[0],
    ];
    if let PhysicalCurve::Native(curve) = &carrier {
        parameters = curve.finite_range(
            parameters,
            edge.curve_same_sense,
            edge.vertices[0] == edge.vertices[1],
        )?;
    }
    if let Some(range) = &edge.source_range {
        parameters = range.each_ref().map(Rational::to_f64_nearest);
    }
    let mut points = Vec::with_capacity(STATIONS + 1);
    let mut preferred = endpoint_fits[0].location.chart;
    for station in 0..=STATIONS {
        dec.charge_work(1)?;
        let t = parameters[0] + (parameters[1] - parameters[0]) * station as f64 / STATIONS as f64;
        let (location, sample) = carrier.sample_at(t, preferred, dec)?;
        preferred = location.chart;
        edge_tolerance.measure(sample.uncertainty_mm)?;
        points.push((location, sample));
        if station == 0 || station == STATIONS {
            let endpoint = usize::from(station == STATIONS);
            let vertex = edge.vertices[endpoint];
            let target = vertices[vertex];
            let residual = norm(subtract(sample.point_mm, target.point_mm))
                + sample.uncertainty_mm
                + target.uncertainty_mm;
            vertex_tolerances[vertex].measure(residual)?;
            edge_tolerance.measure(residual)?;
        }
    }
    if edge.vertices[0] == edge.vertices[1] {
        let first = points.first().expect("finite samples").1;
        let last = points.last().expect("finite samples").1;
        let gap = norm(subtract(first.point_mm, last.point_mm))
            + first.uncertainty_mm
            + last.uncertainty_mm;
        if gap > body.options.initial_mm.to_f64_nearest() {
            return Err(invalid(
                "import.tolerant-broken-topology",
                format!(
                    "{}: finite self-loop curve endpoints are not closed",
                    edge.name
                ),
            ));
        }
        evidence.sampled_checks.push(SampledCheck {
            body: body.name.clone(),
            scope: edge.name.clone(),
            check: "edge.finite-self-loop-endpoints-closed",
            samples: 2,
            maximum_residual_mm: gap,
            converged: 2,
        });
    }
    evidence.sampled_checks.push(SampledCheck {
        body: body.name.clone(),
        scope: edge.name.clone(),
        check: "edge.native-endpoint-fit-and-finite-samples",
        samples: STATIONS + 3,
        maximum_residual_mm: edge_tolerance.residual,
        converged: endpoint_fits.iter().filter(|fit| fit.converged).count() + STATIONS + 1,
    });
    Ok((
        EdgeSamples {
            points,
            source_curve,
            carrier,
        },
        edge_tolerance,
    ))
}

fn closed_surface_samples(
    context: &Context<'_>,
    mut entity: u64,
    surface: &NativeSurface,
    tolerance: &mut Accumulator,
    dec: &mut Decider,
) -> Result<(usize, [[f64; 2]; 2]), Error> {
    use crate::geometry_source::Surface as SourceSurface;
    let mut angular = [false; 2];
    let mut sphere = false;
    for depth in 0..32 {
        match context.parsed.geometry(entity, dec)? {
            Geometry::Surface(
                SourceSurface::Replica(parent, _)
                | SourceSurface::Offset(parent, ..)
                | SourceSurface::CurveBounded(parent, ..),
            ) => entity = parent,
            Geometry::Surface(SourceSurface::RectangularTrim(parent, ..)) => entity = parent,
            Geometry::Surface(SourceSurface::Sphere(..)) => {
                angular[0] = true;
                sphere = true;
                break;
            }
            Geometry::Surface(SourceSurface::Torus(..)) => {
                angular = [true; 2];
                break;
            }
            Geometry::Surface(
                SourceSurface::Cylinder(..)
                | SourceSurface::Cone(..)
                | SourceSurface::Revolution(..),
            ) => {
                angular[0] = true;
                break;
            }
            Geometry::Surface(SourceSurface::BSpline(_)) => break,
            Geometry::Surface(_) => break,
            _ => {
                return Err(invalid(
                    "import.tolerant-broken-reference",
                    "face carrier is not a surface",
                ));
            }
        }
        if depth == 31 {
            return Err(invalid(
                "import.tolerant-broken-reference",
                "cyclic or excessive surface wrappers",
            ));
        }
    }
    let mut limits = [[f64::INFINITY, f64::NEG_INFINITY]; 2];
    for chart in 0..surface.chart_count() {
        for (axis, range) in limits.iter_mut().enumerate() {
            let fallback = if sphere && axis == 1 {
                [
                    -surface.angular_period() / 4.0,
                    surface.angular_period() / 4.0,
                ]
            } else if angular[axis] {
                [0.0, surface.angular_period()]
            } else {
                return_finite_limits(surface.domain(chart)[axis].each_ref())?
            };
            let low = surface.domain(chart)[axis][0]
                .as_ref()
                .map_or(fallback[0], Rational::to_f64_nearest);
            let high = surface.domain(chart)[axis][1]
                .as_ref()
                .map_or(fallback[1], Rational::to_f64_nearest);
            range[0] = range[0].min(low);
            range[1] = range[1].max(high);
        }
    }
    let mut count = 0;
    let closure_width = tolerance.tolerance;
    for axis in 0..2 {
        let mut rims = [
            Vec::with_capacity(STATIONS + 1),
            Vec::with_capacity(STATIONS + 1),
        ];
        for station in 0..=STATIONS {
            for side in 0..2 {
                let mut parameters = [0.0; 2];
                parameters[axis] = limits[axis][side];
                parameters[1 - axis] = limits[1 - axis][0]
                    + (limits[1 - axis][1] - limits[1 - axis][0]) * station as f64
                        / STATIONS as f64;
                let sample = surface.sample(surface.location(parameters)?, dec)?;
                tolerance.measure(sample.uncertainty_mm)?;
                rims[side].push(sample);
                count += 1;
            }
        }
        let seam = rims[0]
            .iter()
            .zip(&rims[1])
            .map(|(a, b)| {
                norm(subtract(a.point_mm, b.point_mm)) + a.uncertainty_mm + b.uncertainty_mm
            })
            .fold(0.0, f64::max);
        let collapsed = rims.iter().all(|rim| {
            rim.iter().all(|sample| {
                norm(subtract(sample.point_mm, rim[0].point_mm))
                    + sample.uncertainty_mm
                    + rim[0].uncertainty_mm
                    <= closure_width
            })
        });
        if seam <= closure_width {
            tolerance.measure(seam)?;
        } else if !collapsed {
            return Err(invalid(
                "import.tolerant-broken-topology",
                "edgeless face has an unmatched, noncollapsed native surface boundary",
            ));
        }
    }
    Ok((count, limits))
}

fn return_finite_limits(limits: [&Option<Rational>; 2]) -> Result<[f64; 2], Error> {
    match limits {
        [Some(low), Some(high)] => Ok([low.to_f64_nearest(), high.to_f64_nearest()]),
        _ => Err(invalid(
            "import.tolerant-broken-topology",
            "edgeless face has an unbounded nonangular parameter",
        )),
    }
}

fn pcurves(
    context: &Context<'_>,
    mut curve: u64,
    surface: u64,
    dec: &mut Decider,
) -> Result<Vec<u64>, Error> {
    for _ in 0..32 {
        match context.parsed.geometry(curve, dec)? {
            Geometry::Curve(SourceCurve::SurfaceCurve(parent, associated, _, _)) => {
                let mut result = Vec::new();
                for entity in associated {
                    if let Geometry::Curve(SourceCurve::Pcurve(basis, _)) =
                        context.parsed.geometry(entity, dec)?
                        && basis == surface
                    {
                        result.push(entity);
                    }
                }
                if !result.is_empty() {
                    return Ok(result);
                }
                curve = parent;
            }
            _ => return Ok(Vec::new()),
        }
    }
    Err(invalid(
        "import.tolerant-broken-reference",
        "cyclic source surface-curve association",
    ))
}

fn manifold(body: &TolerantBody, dec: &mut Decider) -> Result<TolerantCounts, Error> {
    let mut active_edges = vec![false; body.edges.len()];
    let mut active_vertices = vec![false; body.vertices.len()];
    let mut face_count = 0;
    for shell in &body.shells {
        let mut uses = vec![Vec::new(); body.edges.len()];
        for (face_index, face) in shell.faces.iter().enumerate() {
            face_count += 1;
            for boundary in &face.loops {
                if let Some(vertex) = boundary.vertex {
                    active_vertices[vertex] = true;
                }
                for usage in &boundary.edges {
                    dec.charge_work(1)?;
                    active_edges[usage.edge] = true;
                    uses[usage.edge].push((face_index, usage.forward == boundary.orientation));
                    for vertex in body.edges[usage.edge].vertices {
                        active_vertices[vertex] = true;
                    }
                }
            }
        }
        if shell.faces.is_empty() {
            return Err(invalid(
                "import.tolerant-non-manifold",
                "source shell has no faces",
            ));
        }
        for (edge, usage) in uses.iter().enumerate() {
            if usage.is_empty() {
                continue;
            }
            if usage.len() != 2 || usage[0].1 == usage[1].1 {
                return Err(invalid(
                    "import.tolerant-non-manifold",
                    format!(
                        "{} has {} non-opposite shell uses",
                        body.edges[edge].name,
                        usage.len()
                    ),
                ));
            }
        }
        let mut vertex_edges = vec![Vec::new(); body.vertices.len()];
        for (edge_index, usage) in uses.iter().enumerate() {
            if !usage.is_empty() {
                for vertex in body.edges[edge_index].vertices {
                    vertex_edges[vertex].push(edge_index);
                }
            }
        }
        let (mut adjacent, mut visited, mut queue) = (Vec::new(), Vec::new(), Vec::new());
        for incident in vertex_edges {
            if incident.is_empty() {
                continue;
            }
            adjacent.clear();
            visited.clear();
            queue.clear();
            for &edge in &incident {
                for &(face, _) in &uses[edge] {
                    if !adjacent.contains(&face) {
                        adjacent.push(face);
                    }
                }
            }
            queue.push(adjacent[0]);
            while let Some(face) = queue.pop() {
                if visited.contains(&face) {
                    continue;
                }
                visited.push(face);
                for &edge in &incident {
                    if uses[edge].iter().any(|&(value, _)| value == face) {
                        for &(other, _) in &uses[edge] {
                            if !visited.contains(&other) {
                                queue.push(other);
                            }
                        }
                    }
                }
            }
            if visited.len() != adjacent.len() {
                return Err(invalid(
                    "import.tolerant-non-manifold",
                    "source vertex has disconnected face fans",
                ));
            }
        }
    }
    Ok(TolerantCounts {
        faces: face_count,
        edges: active_edges.into_iter().filter(|active| *active).count(),
        vertices: active_vertices.into_iter().filter(|active| *active).count(),
        shells: body.shells.len(),
    })
}

pub(super) fn check_body(
    context: &mut Context<'_>,
    body: &mut TolerantBody,
    evidence: &mut TolerantEvidence,
    dec: &mut Decider,
) -> Result<TolerantCounts, Error> {
    context.stage = super::NativeStage::ManifoldIncidence;
    let mut counts = manifold(body, dec)?;
    context.stage = super::NativeStage::CurveFit;
    let mut vertices = body
        .vertices
        .iter()
        .map(|vertex| point(&vertex.source.point))
        .collect::<Result<Vec<_>, _>>()?;
    let mut vertex_tolerances: Vec<_> = body
        .vertices
        .iter()
        .map(|vertex| Accumulator::new(vertex.tolerance_mm().to_f64_nearest()))
        .collect();
    for (accumulator, sample) in vertex_tolerances.iter_mut().zip(&vertices) {
        accumulator.measure(sample.uncertainty_mm)?;
    }
    let mut edge_tolerances: Vec<_> = body
        .edges
        .iter()
        .map(|edge| Accumulator::new(edge.tolerance_mm().to_f64_nearest()))
        .collect();
    let mut active = vec![false; body.edges.len()];
    for shell in &body.shells {
        for face in &shell.faces {
            for boundary in &face.loops {
                for usage in &boundary.edges {
                    active[usage.edge] = true;
                }
            }
        }
    }
    let mut edges = Vec::with_capacity(body.edges.len());
    for (index, is_active) in active.iter().enumerate() {
        edges.push(if *is_active {
            let (samples, tolerance) = finite_edge(
                context,
                body,
                index,
                &vertices,
                &mut vertex_tolerances,
                evidence,
                dec,
            )?;
            edge_tolerances[index] = tolerance;
            Some(samples)
        } else {
            None
        });
    }
    context.stage = super::NativeStage::CrossingCut;
    crossing::cut(
        context,
        body,
        crossing::State {
            vertices: &mut vertices,
            edges: &mut edges,
            vertex_tolerances: &mut vertex_tolerances,
            edge_tolerances: &mut edge_tolerances,
        },
        evidence,
        dec,
    )?;
    for shell in &body.shells {
        for face in &shell.faces {
            context.stage = super::NativeStage::SurfaceFit;
            dec.check_deadline()?;
            let surface_entity = match id(&face.name) {
                Some(entity) => match context.parsed.topology(entity, dec)? {
                    Topology::Face { surface, .. } => surface,
                    _ => None,
                },
                None => None,
            };
            let surface = surface_entity
                .map(|entity| context.surface(entity, dec))
                .transpose()?;
            let plane = match &face.surface {
                Surface::Plane(placement) => Some(PlaneFrame::new(
                    placement,
                    if surface.is_some() {
                        context.scale.to_f64_nearest()
                    } else {
                        1.0
                    },
                )?),
                _ if surface.is_none() => {
                    return Err(invalid(
                        "import.tolerant-broken-reference",
                        "face has no retained surface declaration",
                    ));
                }
                _ => None,
            };
            let mut tolerance = Accumulator::new(body.options.initial_mm.to_f64_nearest());
            let mut sample_count = 0;
            let mut converged_count = 0;
            let mut fit_residual: f64 = 0.0;
            let mut repeated = BTreeMap::<usize, usize>::new();
            for (boundary_index, boundary) in face.loops.iter().enumerate() {
                let mut previous = None;
                if let Some(vertex) = boundary.vertex {
                    let target = vertices[vertex];
                    let projection = if let Some(plane) = &plane {
                        plane.project(target)
                    } else {
                        surface.as_ref().expect("native surface").project(
                            target.point_mm,
                            None,
                            dec,
                        )?
                    };
                    let residual = projection.residual_mm + target.uncertainty_mm;
                    vertex_tolerances[vertex].measure(residual)?;
                    tolerance.measure(projection.sample.uncertainty_mm)?;
                    fit_residual = fit_residual.max(residual);
                    sample_count += 1;
                    converged_count += usize::from(projection.converged);
                }
                for (use_index, usage) in boundary.edges.iter().enumerate() {
                    let edge = edges[usage.edge].as_ref().expect("active source edge");
                    let associated = match (edge.source_curve, surface_entity) {
                        (Some(curve), Some(surface)) => match pcurves(context, curve, surface, dec)
                        {
                            Ok(curves) => curves,
                            Err(error) if super::super::native::control(&error) => {
                                return Err(error);
                            }
                            Err(_) => Vec::new(),
                        },
                        _ => Vec::new(),
                    };
                    let occurrence = repeated.entry(usage.edge).or_default();
                    let source_pcurve = associated
                        .get(*occurrence % associated.len().max(1))
                        .copied();
                    *occurrence += 1;
                    let pcurve = match source_pcurve
                        .map(|entity| context.pcurve(entity, dec))
                        .transpose()
                    {
                        Ok(curve) => curve,
                        Err(error) if super::super::native::control(&error) => return Err(error),
                        Err(_) => None,
                    };
                    let mut preserve_pcurve = pcurve.is_some();
                    let mut samples = Vec::with_capacity(STATIONS + 1);
                    for station in 0..=STATIONS {
                        dec.charge_work(1)?;
                        let index = if usage.forward {
                            station
                        } else {
                            STATIONS - station
                        };
                        let (curve_location, target) = edge.points[index];
                        let vertex = match index {
                            0 => Some(body.edges[usage.edge].vertices[0]),
                            STATIONS => Some(body.edges[usage.edge].vertices[1]),
                            _ => None,
                        };
                        let fit_target = target;
                        let seed = if station == 0 {
                            match &pcurve {
                                Some(curve) => {
                                    let sampled = curve.sample(
                                        NativeLocation {
                                            chart: 0,
                                            parameters: [curve_location.parameters[0], 0.0],
                                        },
                                        dec,
                                    );
                                    let uv = match sampled {
                                        Ok((uv, _)) => uv,
                                        Err(error) if super::super::native::control(&error) => {
                                            return Err(error);
                                        }
                                        Err(_) => {
                                            preserve_pcurve = false;
                                            previous.map_or([0.0; 2], |value: NativeLocation| {
                                                value.parameters
                                            })
                                        }
                                    };
                                    Some(NativeLocation {
                                        chart: previous
                                            .map_or(0, |value: NativeLocation| value.chart),
                                        parameters: uv,
                                    })
                                }
                                None => previous,
                            }
                        } else {
                            previous
                        };
                        let mut pole_seed = false;
                        let seed = match (&surface, seed) {
                            (Some(surface), Some(seed)) => {
                                pole_seed = station == 0
                                    && previous.is_some()
                                    && poles::collapsed(
                                        surface,
                                        [
                                            previous.expect("previous trace").parameters,
                                            seed.parameters,
                                        ],
                                        fit_target,
                                        body.options.initial_mm.to_f64_nearest(),
                                        dec,
                                    )?
                                    .is_some();
                                let parameters = if pole_seed {
                                    seed.parameters
                                } else {
                                    surface.align_parameters(
                                        seed.parameters,
                                        previous.map(|value| value.parameters),
                                    )
                                };
                                Some(match surface.location(parameters) {
                                    Ok(location) => location,
                                    Err(error) if super::super::native::control(&error) => {
                                        return Err(error);
                                    }
                                    Err(_) => surface.location(
                                        surface.align_parameters(seed.parameters, None),
                                    )?,
                                })
                            }
                            _ => seed,
                        };
                        let mut projection = match (&surface, seed) {
                            (Some(surface), Some(seed)) if pole_seed => {
                                let sample = surface.sample(seed, dec)?;
                                NativeProjection {
                                    location: seed,
                                    sample,
                                    residual_mm: norm(subtract(
                                        sample.point_mm,
                                        fit_target.point_mm,
                                    )),
                                    iterations: 0,
                                    converged: true,
                                }
                            }
                            _ if plane.is_some() => {
                                plane.as_ref().expect("plane").project(fit_target)
                            }
                            (Some(surface), _) => {
                                surface.project(fit_target.point_mm, seed, dec)?
                            }
                            _ => unreachable!("face carrier established"),
                        };
                        if let Some(surface) = &surface
                            && !pole_seed
                            && plane.is_none()
                        {
                            let parameters = surface.align_parameters(
                                projection.location.parameters,
                                previous.map(|value: NativeLocation| value.parameters),
                            );
                            if let Ok(location) = surface.location(parameters) {
                                projection.location = location;
                                projection.sample = surface.sample(location, dec)?;
                                projection.residual_mm =
                                    norm(subtract(projection.sample.point_mm, fit_target.point_mm))
                                        + projection.sample.uncertainty_mm;
                            }
                        }
                        if preserve_pcurve {
                            let retained =
                                (|| {
                                    let (uv, _) = pcurve.as_ref().expect("retained trace").sample(
                                        NativeLocation {
                                            chart: 0,
                                            parameters: [curve_location.parameters[0], 0.0],
                                        },
                                        dec,
                                    )?;
                                    let surface = surface.as_ref().expect("source pcurve surface");
                                    if surface.periods().into_iter().enumerate().any(
                                        |(axis, period)| {
                                            period.is_some_and(|period| {
                                                ((projection.location.parameters[axis] - uv[axis])
                                                    / period)
                                                    .round()
                                                    .abs()
                                                    >= 1.0
                                            })
                                        },
                                    ) {
                                        return Ok::<_, Error>(false);
                                    }
                                    let point = surface.sample(surface.location(uv)?, dec)?;
                                    Ok(norm(subtract(point.point_mm, projection.sample.point_mm))
                                        + point.uncertainty_mm
                                        + projection.sample.uncertainty_mm
                                        <= body.options.initial_mm.to_f64_nearest().max(
                                            64.0 * f64::EPSILON * norm(point.point_mm).max(1.0),
                                        ))
                                })();
                            preserve_pcurve = match retained {
                                Ok(valid) => valid,
                                Err(error) if super::super::native::control(&error) => {
                                    return Err(error);
                                }
                                Err(_) => false,
                            };
                        }
                        let residual = norm(subtract(projection.sample.point_mm, target.point_mm))
                            + projection.sample.uncertainty_mm
                            + target.uncertainty_mm;
                        tolerance.measure(projection.sample.uncertainty_mm)?;
                        edge_tolerances[usage.edge].measure(residual + tolerance.tolerance)?;
                        if let Some(vertex) = vertex {
                            vertex_tolerances[vertex].measure(
                                norm(subtract(
                                    projection.sample.point_mm,
                                    vertices[vertex].point_mm,
                                )) + projection.sample.uncertainty_mm
                                    + vertices[vertex].uncertainty_mm,
                            )?;
                        }
                        fit_residual = fit_residual
                            .max(residual)
                            .max(projection.residual_mm + fit_target.uncertainty_mm);
                        sample_count += 1;
                        converged_count += usize::from(projection.converged);
                        previous = Some(projection.location);
                        samples.push(UvSample {
                            fraction: station as f64 / STATIONS as f64,
                            curve: curve_location,
                            surface: projection.location,
                            residual_mm: residual,
                            converged: projection.converged,
                        });
                    }
                    if let Some(surface) = &surface {
                        let endpoints = body.edges[usage.edge].vertices;
                        let endpoints = if usage.forward {
                            endpoints
                        } else {
                            [endpoints[1], endpoints[0]]
                        };
                        if poles::align_endpoints(
                            surface,
                            &mut samples,
                            endpoints.map(|vertex| vertices[vertex]),
                            body.options.initial_mm.to_f64_nearest(),
                            dec,
                        )? {
                            preserve_pcurve = false;
                        }
                    }
                    let reconstructed = !preserve_pcurve && surface.is_some();
                    let projected = if reconstructed {
                        context.stage = super::NativeStage::TraceRebuild;
                        trace::planar::project(
                            context,
                            surface.as_ref().expect("native surface"),
                            plane.as_ref(),
                            edge,
                            &samples,
                            body.options.initial_mm.to_f64_nearest(),
                            dec,
                        )?
                    } else {
                        None
                    };
                    if reconstructed {
                        let (trace_samples, maximum_residual, trace_converged) = if projected
                            .is_some()
                            || (plane.is_some()
                                && matches!(edge.carrier, PhysicalCurve::Line { .. }))
                        {
                            (
                                samples.len(),
                                samples
                                    .iter()
                                    .map(|sample| sample.residual_mm)
                                    .fold(0.0, f64::max),
                                samples.len(),
                            )
                        } else {
                            let mut fit = trace::Fit {
                                curve: &edge.carrier,
                                surface: surface.as_ref().expect("native trace surface"),
                                precision: body.options.initial_mm.to_f64_nearest(),
                                maximum_residual: 0.0,
                                samples: 0,
                                converged: 0,
                            };
                            fit.rebuild(&mut samples, dec)?;
                            (fit.samples, fit.maximum_residual, fit.converged)
                        };
                        edge_tolerances[usage.edge].measure(maximum_residual)?;
                        sample_count += trace_samples;
                        converged_count += trace_converged;
                        evidence
                            .repairs
                            .push(super::super::healing::HealingChange::TraceRebuilt {
                                body: body.name.clone(),
                                face: face.name.clone(),
                                edge: body.edges[usage.edge].name.clone(),
                                source_pcurve,
                                samples: trace_samples,
                                maximum_residual_mm: format!("{maximum_residual:.17e}"),
                            });
                    }
                    evidence.native_trims.push(NativeTrim {
                        body: body.name.clone(),
                        face: face.name.clone(),
                        boundary: boundary_index,
                        use_index,
                        edge: body.edges[usage.edge].name.clone(),
                        forward: usage.forward,
                        source_pcurve: source_pcurve.filter(|_| preserve_pcurve),
                        samples,
                        reconstructed,
                        projected,
                    });
                }
            }
            if face.loops.iter().all(|boundary| boundary.edges.is_empty()) {
                let surface = surface.as_ref().ok_or_else(|| {
                    invalid(
                        "import.tolerant-broken-reference",
                        "closed face has no native surface",
                    )
                })?;
                let entity = surface_entity.expect("native surface source identity");
                let (samples, bounds) =
                    closed_surface_samples(context, entity, surface, &mut tolerance, dec)?;
                evidence.native_domains.push(super::NativeDomain {
                    body: body.name.clone(),
                    face: face.name.clone(),
                    bounds,
                });
                sample_count += samples;
                converged_count += samples;
            }
            if converged_count < sample_count {
                evidence.warnings.push(NativeWarning { body: Some(body.name.clone()), scope: face.name.clone(), code: "native-uv-fits-not-all-converged", detail: format!("{converged_count}/{sample_count} fit samples reached the stopping width; actual residuals are retained in local tolerances") });
            }
            evidence.sampled_checks.push(SampledCheck {
                body: body.name.clone(),
                scope: face.name.clone(),
                check: "face.native-uv-sampled-fit",
                samples: sample_count,
                maximum_residual_mm: fit_residual.max(tolerance.residual),
                converged: converged_count,
            });
            evidence
                .tolerances
                .push(tolerance.record(&body.name, "face", &face.name));
            for boundary in &face.loops {
                if let Some(vertex) = boundary.vertex {
                    vertex_tolerances[vertex].tolerance =
                        vertex_tolerances[vertex].tolerance.max(tolerance.tolerance);
                }
                for usage in &boundary.edges {
                    edge_tolerances[usage.edge].tolerance = edge_tolerances[usage.edge]
                        .tolerance
                        .max(tolerance.tolerance);
                }
            }
        }
    }
    context.stage = super::NativeStage::ToleranceGrowth;
    for (index, accumulator) in edge_tolerances.iter().enumerate() {
        if !active[index] {
            continue;
        }
        for vertex in body.edges[index].vertices {
            vertex_tolerances[vertex].tolerance = vertex_tolerances[vertex]
                .tolerance
                .max(accumulator.tolerance);
        }
        evidence
            .tolerances
            .push(accumulator.record(&body.name, "edge", &body.edges[index].name));
    }
    for (index, accumulator) in vertex_tolerances.iter().enumerate() {
        evidence.tolerances.push(accumulator.record(
            &body.name,
            "vertex",
            &body.vertices[index].source.name,
        ));
    }
    for (index, accumulator) in edge_tolerances.iter().enumerate() {
        if active[index] {
            body.edges[index].tolerance_mm = Rational::from_f64_exact(accumulator.tolerance)
                .ok_or_else(|| {
                    invalid(
                        "import.tolerant-broken-geometry",
                        "non-finite edge tolerance",
                    )
                })?;
        }
    }
    for (vertex, accumulator) in body.vertices.iter_mut().zip(&vertex_tolerances) {
        vertex.tolerance_mm = Rational::from_f64_exact(accumulator.tolerance).ok_or_else(|| {
            invalid(
                "import.tolerant-broken-geometry",
                "non-finite vertex tolerance",
            )
        })?;
    }
    context.stage = super::NativeStage::PoleBoundary;
    counts.edges += poles::insert(context, body, evidence, dec)?;
    let attained = evidence
        .tolerances
        .iter()
        .filter(|entry| entry.body == body.name)
        .map(|entry| entry.tolerance_mm)
        .fold(0.0, f64::max);
    body.maximum_attained_mm = Rational::from_f64_exact(attained).ok_or_else(|| {
        invalid(
            "import.tolerant-broken-geometry",
            "non-finite body tolerance",
        )
    })?;
    body.options.maximum_mm = body
        .options
        .maximum_mm
        .clone()
        .max(body.maximum_attained_mm.clone());
    Ok(counts)
}
