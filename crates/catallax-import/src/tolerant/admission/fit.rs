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

const STATIONS: usize = 16;

pub(super) struct Context<'a> {
    parsed: &'a ParsedStep,
    scale: &'a Rational,
    angles: AngleUnit,
    curves: BTreeMap<u64, Arc<NativeCurve>>,
    surfaces: BTreeMap<u64, Arc<NativeSurface>>,
    pcurves: BTreeMap<u64, Arc<NativeUvCurve>>,
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
    ) -> Self {
        Self {
            parsed,
            scale,
            angles,
            curves,
            surfaces: BTreeMap::new(),
            pcurves: BTreeMap::new(),
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
    fn pcurve(
        &mut self,
        entity: u64,
        dec: &mut Decider,
    ) -> Result<Option<Arc<NativeUvCurve>>, Error> {
        if let Some(curve) = self.pcurves.get(&entity) {
            return Ok(Some(Arc::clone(curve)));
        }
        let curve = match NativeUvCurve::from_source(
            self.parsed,
            entity,
            self.scale,
            self.angles.clone(),
            dec,
        ) {
            Ok(curve) => Arc::new(curve),
            Err(error) if crate::tolerant::native::control(&error) => return Err(error),
            // This carrier is only an optional hint. Required 3D curve and
            // surface compilation above still propagates its own refusals.
            Err(_) => return Ok(None),
        };
        self.pcurves.insert(entity, Arc::clone(&curve));
        Ok(Some(curve))
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
}
impl PlaneFrame {
    fn new(placement: &Placement) -> Result<Self, Error> {
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
        })
    }
    fn project(&self, target: NativeSample) -> NativeProjection {
        let delta = subtract(target.point_mm, self.origin);
        let parameters = [dot(delta, self.x), dot(delta, self.y)];
        let position = std::array::from_fn(|i| {
            self.origin[i] + parameters[0] * self.x[i] + parameters[1] * self.y[i]
        });
        NativeProjection {
            location: NativeLocation {
                chart: 0,
                parameters,
            },
            sample: NativeSample {
                point_mm: position,
                uncertainty_mm: 8.0 * f64::EPSILON * norm(position),
            },
            residual_mm: norm(subtract(position, target.point_mm))
                + target.uncertainty_mm
                + 8.0 * f64::EPSILON * norm(position),
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
}
struct EdgeSamples {
    points: Vec<(NativeLocation, NativeSample)>,
    source_curve: Option<u64>,
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
    let source_curve = match id(&edge.name) {
        Some(entity) => match context.parsed.topology(entity, dec)? {
            Topology::Edge { curve, .. } => curve,
            _ => None,
        },
        None => None,
    };
    let carrier = match source_curve {
        Some(entity) => PhysicalCurve::Native(context.curve(entity, dec)?),
        None => match &edge.curve {
            Curve::Line(origin, direction) => PhysicalCurve::Line {
                origin: point(origin)?.point_mm,
                direction: point(direction)?.point_mm,
            },
            _ => {
                return Err(invalid(
                    "import.tolerant-broken-reference",
                    format!("{} has no source curve declaration", edge.name),
                ));
            }
        },
    };
    let mut endpoint_fits = Vec::with_capacity(2);
    for &vertex in &edge.vertices {
        let target = vertices[vertex];
        let fit = match &carrier {
            PhysicalCurve::Native(curve) => curve.project(target.point_mm, None, dec)?,
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
    let mut points = Vec::with_capacity(STATIONS + 1);
    let mut preferred = endpoint_fits[0].location.chart;
    for station in 0..=STATIONS {
        dec.charge_work(1)?;
        let t = parameters[0] + (parameters[1] - parameters[0]) * station as f64 / STATIONS as f64;
        let (location, sample) = match &carrier {
            PhysicalCurve::Native(curve) => {
                preferred = chart_at(curve, t, preferred);
                let location = NativeLocation {
                    chart: preferred,
                    parameters: [t, 0.0],
                };
                (location, curve.sample(location, dec)?)
            }
            PhysicalCurve::Line { origin, direction } => {
                let position = std::array::from_fn(|i| origin[i] + t * direction[i]);
                (
                    NativeLocation {
                        chart: 0,
                        parameters: [t, 0.0],
                    },
                    NativeSample {
                        point_mm: position,
                        uncertainty_mm: 8.0 * f64::EPSILON * norm(position),
                    },
                )
            }
        };
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

fn pcurve_seed(
    surface: &NativeSurface,
    pcurve: Option<&NativeUvCurve>,
    parameter: f64,
    previous: Option<NativeLocation>,
    dec: &mut Decider,
) -> Result<Option<NativeLocation>, Error> {
    let Some(uv) = pcurve
        .map(|curve| curve.hint(parameter, dec))
        .transpose()?
        .flatten()
    else {
        return Ok(previous);
    };
    let uv = surface.align_parameters(uv, previous.map(|location| location.parameters));
    Ok(surface.location(uv).ok().or(previous))
}

#[cfg(test)]
mod tests;

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
    let counts = manifold(body, dec)?;
    let vertices = body
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
    for shell in &body.shells {
        for face in &shell.faces {
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
            let plane = if surface.is_none() {
                match &face.surface {
                    Surface::Plane(placement) => Some(PlaneFrame::new(placement)?),
                    _ => {
                        return Err(invalid(
                            "import.tolerant-broken-reference",
                            "face has no retained surface declaration",
                        ));
                    }
                }
            } else {
                None
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
                    let projection = match &surface {
                        Some(surface) => surface.project(target.point_mm, None, dec)?,
                        None => plane.as_ref().expect("inferred plane").project(target),
                    };
                    let residual = projection.residual_mm + target.uncertainty_mm;
                    vertex_tolerances[vertex].measure(residual)?;
                    tolerance.measure(projection.sample.uncertainty_mm)?;
                    fit_residual = fit_residual.max(residual);
                    sample_count += 1;
                    converged_count += usize::from(projection.converged);
                }
                for usage in &boundary.edges {
                    let edge = edges[usage.edge].as_ref().expect("active source edge");
                    let associated = match (edge.source_curve, surface_entity) {
                        (Some(curve), Some(surface)) => pcurves(context, curve, surface, dec)?,
                        _ => Vec::new(),
                    };
                    let occurrence = repeated.entry(usage.edge).or_default();
                    let source_pcurve = associated
                        .get(*occurrence % associated.len().max(1))
                        .copied();
                    *occurrence += 1;
                    let pcurve = source_pcurve
                        .map(|entity| context.pcurve(entity, dec))
                        .transpose()?
                        .flatten();
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
                        let fit_target = vertex.map_or(target, |vertex| vertices[vertex]);
                        let seed = match &surface {
                            Some(surface) => pcurve_seed(
                                surface,
                                pcurve.as_deref(),
                                curve_location.parameters[0],
                                previous,
                                dec,
                            )?,
                            None => previous,
                        };
                        let mut projection = match &surface {
                            Some(surface) => surface.project(fit_target.point_mm, seed, dec)?,
                            None => plane.as_ref().expect("inferred plane").project(fit_target),
                        };
                        if let Some(surface) = &surface {
                            let parameters = surface.align_parameters(
                                projection.location.parameters,
                                previous.map(|location: NativeLocation| location.parameters),
                            );
                            if parameters.map(f64::to_bits)
                                != projection.location.parameters.map(f64::to_bits)
                            {
                                projection.location = surface.location(parameters)?;
                                projection.sample = surface.sample(projection.location, dec)?;
                                projection.residual_mm =
                                    norm(subtract(projection.sample.point_mm, fit_target.point_mm))
                                        + projection.sample.uncertainty_mm;
                            }
                        }
                        let residual = norm(subtract(projection.sample.point_mm, target.point_mm))
                            + projection.sample.uncertainty_mm
                            + target.uncertainty_mm;
                        tolerance.measure(projection.sample.uncertainty_mm)?;
                        edge_tolerances[usage.edge].measure(residual + tolerance.tolerance)?;
                        if let Some(vertex) = vertex {
                            vertex_tolerances[vertex]
                                .measure(projection.residual_mm + fit_target.uncertainty_mm)?;
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
                    evidence.native_trims.push(NativeTrim {
                        body: body.name.clone(),
                        face: face.name.clone(),
                        boundary: boundary_index,
                        edge: body.edges[usage.edge].name.clone(),
                        forward: usage.forward,
                        source_pcurve,
                        samples,
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
