//! Display sampling of admitted native carriers. No exact Model is fabricated.
use std::collections::BTreeMap;

use catallax_import::{
    Curve,
    geometry_source::Topology,
    tolerant::{
        admission::{NativeTrim, TolerantImport},
        native::{NativeCurve, NativeLocation, NativeSurface},
    },
};
use catallax_numbers::{Decider, Rational, Refusal};
use spade::Triangulation;

use super::imported::{
    Cdt, FaceOutput, Node, SourceChart, SurfaceSamples, finish_face, insert, push, seed_rectangles,
};
use crate::{PreviewEdge, PreviewMesh, PreviewOptions};

type P = [f64; 3];
type Uv = [f64; 2];
fn failure(detail: impl Into<String>) -> Refusal {
    Refusal::new("mesh.tolerant-preview-unresolved", detail)
}
fn source_error(error: catallax_import::tolerant::Error) -> Refusal {
    match error {
        catallax_import::tolerant::Error::Numerical(error) => error,
        other => failure(format!("{other:?}")),
    }
}
fn distance(a: P, b: P) -> f64 {
    (0..3).fold(0.0, |norm, axis| norm.hypot(a[axis] - b[axis]))
}
fn point(value: &[Rational; 3]) -> Result<P, Refusal> {
    let point = value.each_ref().map(Rational::to_f64_nearest);
    if point.iter().all(|coordinate| coordinate.is_finite()) {
        Ok(point)
    } else {
        Err(failure("vertex exceeds binary64 range"))
    }
}

enum Carrier<'a> {
    Native(&'a NativeSurface),
    Inferred(SourceChart),
}
impl SurfaceSamples for Carrier<'_> {
    fn position(&self, uv: Uv, dec: &mut Decider) -> Result<P, Refusal> {
        match self {
            Self::Native(surface) => surface
                .sample(surface.location(uv).map_err(source_error)?, dec)
                .map(|sample| sample.point_mm)
                .map_err(source_error),
            Self::Inferred(surface) => SurfaceSamples::position(surface, uv, dec),
        }
    }
    fn normal(&self, uv: Uv, dec: &mut Decider) -> Result<P, Refusal> {
        let Self::Native(surface) = self else {
            let Self::Inferred(surface) = self else {
                unreachable!()
            };
            return SurfaceSamples::normal(surface, uv, dec);
        };
        let location = surface.location(uv).map_err(source_error)?;
        surface
            .sampled_normal(location, dec)
            .map_err(source_error)?
            .ok_or_else(|| failure("native sampled normal collapsed"))
    }
}

enum EdgeCarrier<'a> {
    Native(&'a NativeCurve, [f64; 2]),
    Line([P; 2]),
}
impl EdgeCarrier<'_> {
    fn sample(&self, fraction: f64, dec: &mut Decider) -> Result<P, Refusal> {
        match self {
            Self::Native(curve, range) => {
                let parameter = range[0] + fraction * (range[1] - range[0]);
                curve
                    .sample(curve.location(parameter).map_err(source_error)?, dec)
                    .map(|sample| sample.point_mm)
                    .map_err(source_error)
            }
            Self::Line(ends) => Ok(std::array::from_fn(|axis| {
                ends[0][axis] + fraction * (ends[1][axis] - ends[0][axis])
            })),
        }
    }
}

fn edge_samples(
    carrier: &EdgeCarrier<'_>,
    ends: [usize; 2],
    target: f64,
    mesh: &mut PreviewMesh<String, String>,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<Vec<(f64, usize)>, Refusal> {
    let mut output = vec![(0.0, ends[0])];
    let mut stack = vec![(0.0, 1.0, ends, 0)];
    while let Some((a, b, indices, depth)) = stack.pop() {
        dec.charge_work(1)?;
        let mid = (a + b) * 0.5;
        let actual = carrier.sample(mid, dec)?;
        let linear = std::array::from_fn(|axis| {
            (mesh.positions[indices[0]][axis] + mesh.positions[indices[1]][axis]) * 0.5
        });
        if distance(actual, linear) <= target * 0.5 {
            output.push((b, indices[1]));
        } else {
            if depth == 24 || mid <= a || mid >= b {
                return Err(failure("edge subdivision reached its sampling limit"));
            }
            let index = push(mesh, actual, options)?;
            stack.push((mid, b, [index, indices[1]], depth + 1));
            stack.push((a, mid, [indices[0], index], depth + 1));
        }
    }
    Ok(output)
}

fn trim_seed(trim: &NativeTrim, fraction: f64) -> NativeLocation {
    let next = trim
        .samples
        .partition_point(|sample| sample.fraction < fraction)
        .min(trim.samples.len() - 1)
        .max(1);
    let first = next - 1;
    let width = trim.samples[next].fraction - trim.samples[first].fraction;
    let weight = ((fraction - trim.samples[first].fraction) / width).clamp(0.0, 1.0);
    let a = trim.samples[first].surface;
    let b = trim.samples[next].surface;
    NativeLocation {
        chart: a.chart,
        parameters: std::array::from_fn(|axis| {
            a.parameters[axis] + weight * (b.parameters[axis] - a.parameters[axis])
        }),
    }
}

fn domain_loop(
    carrier: &Carrier<'_>,
    bounds: [[f64; 2]; 2],
    mesh: &mut PreviewMesh<String, String>,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<Vec<Node>, Refusal> {
    let corners = [
        [bounds[0][0], bounds[1][0]],
        [bounds[0][1], bounds[1][0]],
        [bounds[0][1], bounds[1][1]],
        [bounds[0][0], bounds[1][1]],
    ];
    let mut nodes = Vec::new();
    for side in 0..4 {
        let a = corners[side];
        let b = corners[(side + 1) % 4];
        let start = carrier.position(a, dec)?;
        let end = carrier.position(b, dec)?;
        let mut stack = vec![(a, b, start, end, 0)];
        while let Some((a, b, start, end, depth)) = stack.pop() {
            dec.charge_work(1)?;
            let mid = std::array::from_fn(|axis| (a[axis] + b[axis]) * 0.5);
            let actual = carrier.position(mid, dec)?;
            let linear = std::array::from_fn(|axis| (start[axis] + end[axis]) * 0.5);
            if distance(actual, linear) <= options.chord_tolerance.to_f64_nearest() * 0.5 {
                nodes.push(Node {
                    uv: a,
                    mesh: push(mesh, start, options)?,
                });
            } else {
                if depth == 24
                    || mid.map(f64::to_bits) == a.map(f64::to_bits)
                    || mid.map(f64::to_bits) == b.map(f64::to_bits)
                {
                    return Err(failure(
                        "natural seam subdivision reached its sampling limit",
                    ));
                }
                stack.push((mid, b, actual, end, depth + 1));
                stack.push((a, mid, start, actual, depth + 1));
            }
        }
    }
    Ok(nodes)
}

// Cut two noncontractible rims along one common source-parameter seam. The
// additional cut nodes are display samples, not new pickable B-rep edges.
fn cut_band(
    mut rims: Vec<Vec<Node>>,
    carrier: &Carrier<'_>,
    mesh: &mut PreviewMesh<String, String>,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<Vec<Node>, Refusal> {
    if rims.len() != 2 {
        return Err(failure("periodic trim needs a paired-rim seam"));
    }
    let delta = |rim: &[Node]| {
        [
            rim.last().unwrap().uv[0] - rim[0].uv[0],
            rim.last().unwrap().uv[1] - rim[0].uv[1],
        ]
    };
    let first = delta(&rims[0]);
    let axis = usize::from(first[1].abs() > first[0].abs());
    let period = first[axis].abs();
    if !period.is_finite() || period == 0.0 {
        return Err(failure("periodic rim has no finite winding"));
    }
    let cut = rims[0][0].uv[axis];
    for rim in &mut rims {
        let winding = delta(rim);
        if (winding[axis].abs() - period).abs() > period * 1e-6
            || winding[1 - axis].abs() > period * 1e-6
        {
            return Err(failure("paired rims have incompatible native winding"));
        }
        if winding[axis] < 0.0 {
            rim.reverse();
        }
        let shift = ((cut - rim[0].uv[axis]) / period).floor() * period;
        for node in rim.iter_mut() {
            node.uv[axis] += shift;
        }
        let target = if cut < rim[0].uv[axis] {
            cut + period
        } else {
            cut
        };
        let split = rim
            .windows(2)
            .position(|pair| pair[0].uv[axis] <= target && target <= pair[1].uv[axis])
            .ok_or_else(|| failure("native periodic rim is not monotone at the chosen seam"))?;
        let a = rim[split];
        let b = rim[split + 1];
        let fraction = if (b.uv[axis] - a.uv[axis]).abs() == 0.0 {
            0.0
        } else {
            (target - a.uv[axis]) / (b.uv[axis] - a.uv[axis])
        };
        let uv = std::array::from_fn(|i| a.uv[i] + fraction * (b.uv[i] - a.uv[i]));
        let cut_node = if fraction == 0.0 {
            a
        } else if fraction.to_bits() == 1.0_f64.to_bits() {
            b
        } else {
            Node {
                uv,
                mesh: push(mesh, carrier.position(uv, dec)?, options)?,
            }
        };
        let mut rotated = Vec::with_capacity(rim.len() + 1);
        rotated.push(cut_node);
        rotated.extend(rim[split + 1..rim.len() - 1].iter().copied());
        for node in &rim[..=split] {
            let mut node = *node;
            node.uv[axis] += period;
            rotated.push(node);
        }
        let mut end = cut_node;
        end.uv[axis] += period;
        rotated.push(end);
        rotated.dedup_by(|a, b| a.uv.iter().zip(b.uv).all(|(a, b)| (*a - b).abs() == 0.0));
        let shift = cut - rotated[0].uv[axis];
        if (shift / period - (shift / period).round()).abs() > 1e-6 {
            return Err(failure(
                "periodic seam alignment changed source parameterization",
            ));
        }
        for node in &mut rotated {
            node.uv[axis] += shift;
        }
        *rim = rotated;
    }
    let mut first = rims.remove(0);
    first.extend(rims.remove(0).into_iter().rev());
    Ok(first)
}

/// Preview one admitted native definition in its source-body coordinates.
/// Assembly placement remains in the source document, not in this body-local
/// preview. The returned triangles have no watertightness or deviation claim.
pub fn preview_tolerant_body(
    import: &TolerantImport,
    index: usize,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<PreviewMesh<String, String>, Refusal> {
    super::validate(options)?;
    let body = import
        .bodies()
        .get(index)
        .ok_or_else(|| failure("body index is out of range"))?
        .definition();
    let trims: Vec<_> = import
        .evidence()
        .native_trims
        .iter()
        .filter(|trim| trim.body == body.name())
        .collect();
    let mut first_trim = BTreeMap::new();
    let mut face_trims = BTreeMap::<(&str, usize), Vec<&NativeTrim>>::new();
    for trim in &trims {
        first_trim.entry(trim.edge.as_str()).or_insert(*trim);
        face_trims
            .entry((&trim.face, trim.boundary))
            .or_default()
            .push(*trim);
    }
    for trims in face_trims.values_mut() {
        trims.sort_by_key(|trim| trim.use_index);
    }
    let mut mesh = PreviewMesh::default();
    let vertices = body
        .vertices()
        .iter()
        .map(|vertex| push(&mut mesh, point(&vertex.source.point)?, options))
        .collect::<Result<Vec<_>, _>>()?;
    let mut edges = BTreeMap::new();
    for edge in body.edges() {
        let Some(trim) = first_trim.get(edge.name.as_str()) else {
            continue;
        };
        let ends = edge.vertices.map(|vertex| vertices[vertex]);
        let source_id = edge
            .name
            .strip_prefix('#')
            .and_then(|name| name.parse().ok());
        let curve = match source_id
            .map(|entity| import.parsed_source().topology(entity, dec))
            .transpose()
            .map_err(|error| failure(format!("{error:?}")))?
        {
            Some(Topology::Edge {
                curve: Some(curve), ..
            }) => import.native_curve(curve),
            _ => None,
        };
        let carrier = if let Some(curve) = curve {
            let first = trim.samples.first().unwrap().curve.parameters[0];
            let last = trim.samples.last().unwrap().curve.parameters[0];
            EdgeCarrier::Native(
                curve,
                if trim.forward {
                    [first, last]
                } else {
                    [last, first]
                },
            )
        } else if edge.degenerate || matches!(edge.curve, Curve::Line(..)) {
            EdgeCarrier::Line(ends.map(|index| mesh.positions[index]))
        } else {
            return Err(failure("active edge has no retained native carrier"));
        };
        let samples = edge_samples(
            &carrier,
            ends,
            options
                .chord_tolerance
                .to_f64_nearest()
                .max(edge.tolerance_mm().to_f64_nearest()),
            &mut mesh,
            options,
            dec,
        )?;
        mesh.edges.push(PreviewEdge {
            edge: edge.name.clone(),
            vertices: samples.iter().map(|(_, index)| *index).collect(),
        });
        edges.insert(edge.name.as_str(), samples);
    }
    for face in body.shells().iter().flat_map(|shell| &shell.faces) {
        let source_id = face
            .name
            .strip_prefix('#')
            .and_then(|name| name.parse().ok());
        let surface = match source_id
            .map(|entity| import.parsed_source().topology(entity, dec))
            .transpose()
            .map_err(|error| failure(format!("{error:?}")))?
        {
            Some(Topology::Face {
                surface: Some(surface),
                ..
            }) => import.native_surface(surface),
            _ => None,
        };
        let carrier = match surface {
            Some(surface) => Carrier::Native(surface),
            None => Carrier::Inferred(SourceChart::new(&face.surface, dec)?),
        };
        let mut face_options = options.clone();
        if let Some(tolerance) = import.evidence().tolerances.iter().find(|entry| {
            entry.body == body.name() && entry.kind == "face" && entry.entity == face.name
        }) {
            face_options.chord_tolerance = face_options.chord_tolerance.max(
                Rational::from_f64_exact(tolerance.tolerance_mm)
                    .ok_or_else(|| failure("nonfinite face tolerance"))?,
            );
        }
        let mut loops = Vec::new();
        let mut winding = Vec::new();
        for (boundary_index, boundary) in face.loops.iter().enumerate() {
            if boundary.edges.is_empty() {
                continue;
            }
            let trim_rows = face_trims
                .get(&(face.name.as_str(), boundary_index))
                .ok_or_else(|| failure("native trim evidence is missing"))?;
            let mut nodes: Vec<Node> = Vec::new();
            for offset in 0..boundary.edges.len() {
                let usage_index = if boundary.orientation {
                    offset
                } else {
                    boundary.edges.len() - 1 - offset
                };
                let usage = &boundary.edges[usage_index];
                let trim = trim_rows[usage_index];
                let edge = &body.edges()[usage.edge];
                let row = &edges[edge.name.as_str()];
                for offset in 0..row.len() {
                    let (fraction, index) = row[if usage.forward == boundary.orientation {
                        offset
                    } else {
                        row.len() - 1 - offset
                    }];
                    let seed = trim_seed(
                        trim,
                        if usage.forward {
                            fraction
                        } else {
                            1.0 - fraction
                        },
                    );
                    let uv = match &carrier {
                        Carrier::Native(surface) => {
                            surface
                                .project(mesh.positions[index], Some(seed), dec)
                                .map_err(source_error)?
                                .location
                                .parameters
                        }
                        Carrier::Inferred(_) => {
                            trim_seed(
                                trim,
                                if usage.forward {
                                    fraction
                                } else {
                                    1.0 - fraction
                                },
                            )
                            .parameters
                        }
                    };
                    if nodes.last().is_some_and(|node| {
                        node.mesh == index && (node.uv[0] - uv[0]).hypot(node.uv[1] - uv[1]) <= 1e-8
                    }) {
                        continue;
                    }
                    nodes.push(Node { uv, mesh: index });
                }
            }
            if nodes.len() < 4 || nodes[0].mesh != nodes.last().unwrap().mesh {
                return Err(failure("native trim loop is not closed"));
            }
            let gap = (nodes[0].uv[0] - nodes.last().unwrap().uv[0])
                .hypot(nodes[0].uv[1] - nodes.last().unwrap().uv[1]);
            if gap > 1e-8 {
                winding.push(nodes);
            } else {
                nodes.pop();
                loops.push(nodes);
            }
        }
        let natural = import
            .evidence()
            .native_domains
            .iter()
            .find(|domain| domain.body == body.name() && domain.face == face.name);
        if let Some(domain) = natural {
            loops.push(domain_loop(
                &carrier,
                domain.bounds,
                &mut mesh,
                &face_options,
                dec,
            )?);
        }
        if !winding.is_empty() {
            loops.push(cut_band(winding, &carrier, &mut mesh, &face_options, dec)?);
        }
        let first = loops
            .first()
            .and_then(|nodes| nodes.first())
            .ok_or_else(|| failure("face has no sampled parameter region"))?;
        let mut bounds = [first.uv; 2];
        let mut cdt = Cdt::new();
        for boundary in &loops {
            for i in 0..boundary.len() {
                dec.charge_work(1)?;
                for axis in 0..2 {
                    bounds[0][axis] = bounds[0][axis].min(boundary[i].uv[axis]);
                    bounds[1][axis] = bounds[1][axis].max(boundary[i].uv[axis]);
                }
                let a = insert(&mut cdt, boundary[i])?;
                let b = insert(&mut cdt, boundary[(i + 1) % boundary.len()])?;
                if a == b || !cdt.can_add_constraint(a, b) {
                    return Err(failure("native preview trim constraints cross or collapse"));
                }
                cdt.add_constraint(a, b);
            }
        }
        seed_rectangles(
            &carrier,
            bounds,
            &loops,
            &mut cdt,
            &mut mesh,
            &face_options,
            dec,
        )?;
        finish_face(
            FaceOutput {
                name: &face.name,
                same_sense: face.same_sense,
                skip_collapsed: natural.is_some()
                    || face
                        .loops
                        .iter()
                        .flat_map(|boundary| &boundary.edges)
                        .any(|usage| body.edges()[usage.edge].degenerate),
            },
            &carrier,
            cdt,
            &loops,
            &mut mesh,
            &face_options,
            dec,
        )?;
    }
    dec.check_deadline()?;
    Ok(mesh)
}
