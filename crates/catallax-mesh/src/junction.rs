//! Structured display grids for checked finite spline faces.
//!
//! Each supported face supplies a grid and a chord bound. Shared boundary rows
//! seed the edge polylines, and analytic neighbours triangulate up to those rows.

use std::collections::BTreeMap;

mod bezier;
mod rational_patch;

use catallax_numbers::{Ball, Decider, Rational, Real, Refusal};
use catallax_shapes::{
    BoundaryRow, Curve, EnvelopeGrid, LinearFillet, RollingMitre, RollingPipe, RollingSaddle,
    SplineSupport, Surface,
};
use catallax_topology::{
    BodyId, EdgeId, FaceId, LinearFilletBand, Model, RollingMitreBand, Sense, envelope_rows,
    pipe_rows, saddle_rows,
};

use crate::{
    FaceMesh, MeshControl, MeshOptions, Triangle, UnavailableControl,
    chart::{self, BITS},
    edges::Boundaries,
    failure,
    preview::{PreviewEdge, PreviewFace, PreviewMesh, PreviewOptions},
};

/// One boundary row's grid points, with the row's distance to its grid polyline.
pub(crate) struct Seed {
    pub(crate) points: Vec<[Ball; 3]>,
    pub(crate) error: Rational,
    /// Whether the points run from the edge's start vertex.
    pub(crate) forward: bool,
    /// An exact boundary row's points. A planar face needs them to triangulate, so
    /// a rail with no planar neighbour may leave this empty.
    pub(crate) exact: Vec<[Real; 3]>,
}

struct Envelope {
    grid: EnvelopeGrid,
    rows: Vec<(EdgeId, BoundaryRow, bool)>,
    // Whether the face boundary walks the rail along increasing chart u.
    along_u: bool,
    ends: Ends,
    rational: Option<rational_patch::Proof>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Ends {
    Meridians,
    Seam,
    Poles,
    ZeroStart,
    ZeroEnd,
}

impl Ends {
    fn collapsed(self, column: usize, columns: usize) -> bool {
        (column == 0 && matches!(self, Self::Poles | Self::ZeroStart))
            || (column + 1 == columns && matches!(self, Self::Poles | Self::ZeroEnd))
    }
}

pub(crate) struct Junction {
    envelopes: BTreeMap<FaceId, Envelope>,
    seeds: BTreeMap<EdgeId, Seed>,
}

fn row_points(grid: &EnvelopeGrid, row: BoundaryRow) -> Vec<[Ball; 3]> {
    let at = |column: usize, line: usize| grid.points[column * grid.rows + line].clone();
    match row {
        BoundaryRow::Rail => (0..grid.columns).map(|c| at(c, 0)).collect(),
        BoundaryRow::Mouth => (0..grid.columns).map(|c| at(c, grid.rows - 1)).collect(),
        BoundaryRow::Negative => (0..grid.rows).map(|r| at(0, r)).collect(),
        BoundaryRow::Positive => (0..grid.rows).map(|r| at(grid.columns - 1, r)).collect(),
    }
}

/// Whether every spline face has a checked structured display route.
pub(crate) fn present(model: &Model, body: BodyId) -> bool {
    let mut found = false;
    for face in model.faces_of_body(body) {
        match model.geom().surface(model.face(face).surface) {
            Surface::BSpline(
                SplineSupport::RollingEnvelope(_)
                | SplineSupport::BezierExtrusion(_)
                | SplineSupport::RollingPipe(_)
                | SplineSupport::RollingSaddle(_)
                | SplineSupport::LinearFillet(_)
                | SplineSupport::RollingMitre(_)
                | SplineSupport::RationalPatch(_),
            ) => found = true,
            Surface::BSpline(_) => return false,
            _ => {}
        }
    }
    found
}

/// The display mesh at the preview's target, without its bounds.
pub(crate) fn preview(
    model: &Model,
    body: BodyId,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<PreviewMesh, Refusal> {
    let mesh = crate::triangulate::build(
        model,
        body,
        &MeshOptions {
            chord_tolerance: options.chord_tolerance.clone(),
            angular_deviation: None,
            max_vertices: options.max_vertices,
            max_edge_length: None,
            max_aspect_ratio: None,
        },
        dec,
    )?
    .value;
    Ok(PreviewMesh {
        positions: mesh.positions,
        triangles: mesh.triangles,
        faces: mesh
            .faces
            .into_iter()
            .map(|face| PreviewFace {
                face: face.face,
                triangles: face.triangles,
            })
            .collect(),
        edges: mesh
            .edges
            .into_iter()
            .map(|edge| PreviewEdge {
                edge: edge.edge,
                vertices: edge.vertices,
            })
            .collect(),
    })
}

/// The grids of a body's envelope faces. `None` for a body without one, or with
/// any other spline face, which keeps its own route.
pub(crate) fn prepare(
    model: &Model,
    faces: &[FaceId],
    options: &MeshOptions,
    dec: &mut Decider,
) -> Result<Option<Junction>, Refusal> {
    let mut found = Vec::new();
    let mut bands = Vec::new();
    let mut saddles = Vec::new();
    let mut variable = Vec::new();
    let mut mitres = Vec::new();
    let mut beziers = Vec::new();
    let mut rational = Vec::new();
    for &face in faces {
        match model.geom().surface(model.face(face).surface) {
            Surface::BSpline(SplineSupport::RationalPatch(_)) => rational.push(face),
            Surface::BSpline(SplineSupport::BezierExtrusion(_)) => beziers.push(face),
            Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) => {
                found.push((face, envelope))
            }
            Surface::BSpline(SplineSupport::RollingPipe(pipe)) => bands.push((face, pipe)),
            Surface::BSpline(SplineSupport::RollingSaddle(saddle)) => saddles.push((face, saddle)),
            Surface::BSpline(SplineSupport::LinearFillet(band)) => variable.push((face, band)),
            Surface::BSpline(SplineSupport::RollingMitre(band)) => mitres.push((face, band)),
            Surface::BSpline(_) => {
                if !rational.is_empty()
                    || faces.iter().any(|face| {
                        matches!(
                            model.geom().surface(model.face(*face).surface),
                            Surface::BSpline(SplineSupport::RationalPatch(_))
                        )
                    })
                {
                    return Err(Refusal::new(
                        "mesh.rational-patch-junction",
                        "mixed generic and law-specific spline grids need a common boundary certificate",
                    ));
                }
                return Ok(None);
            }
            _ => {}
        }
    }
    if !rational.is_empty() {
        if !found.is_empty()
            || !bands.is_empty()
            || !saddles.is_empty()
            || !variable.is_empty()
            || !mitres.is_empty()
            || !beziers.is_empty()
        {
            return Err(Refusal::new(
                "mesh.rational-patch-junction",
                "mixed generic and law-specific spline grids need a common boundary certificate",
            ));
        }
        return rational_patch::prepare(model, &rational, options, dec).map(Some);
    }
    if found.is_empty()
        && bands.is_empty()
        && saddles.is_empty()
        && variable.is_empty()
        && mitres.is_empty()
        && beziers.is_empty()
    {
        return Ok(None);
    }
    let tolerance = options.chord_tolerance.mul_pow2(-1);
    let mut envelopes = BTreeMap::new();
    let mut seeds = BTreeMap::new();
    for face in beziers {
        envelopes.insert(
            face,
            bezier::band(model, face, &tolerance, options, &mut seeds, dec)?,
        );
    }
    for (face, pipe) in bands {
        envelopes.insert(
            face,
            band(model, face, pipe, &tolerance, options, &mut seeds, dec)?,
        );
    }
    for (face, saddle) in saddles {
        envelopes.insert(
            face,
            saddle_band(model, face, saddle, &tolerance, options, &mut seeds, dec)?,
        );
    }
    for (face, support) in variable {
        envelopes.insert(
            face,
            linear_band(model, face, support, &tolerance, options, &mut seeds, dec)?,
        );
    }
    for (face, support) in mitres {
        envelopes.insert(
            face,
            mitre_band(model, face, support, &tolerance, options, &mut seeds, dec)?,
        );
    }
    for (face, envelope) in found {
        let grid = envelope.display_grid(&tolerance, options.max_vertices, dec)?;
        let mut rows = Vec::with_capacity(4);
        let mut along_u = None;
        for (edge, row, sense) in envelope_rows(model, face, dec)? {
            let points = row_points(&grid, row);
            let definition = model.edge(edge);
            let (Some(start), Some(end)) = (definition.start, definition.end) else {
                return Err(failure("a rolling-envelope row has two ends"));
            };
            let [start, end] = [start, end].map(|v| model.geom().position(model.vertex(v).point));
            let distance = |point: &catallax_shapes::Point3| -> Result<Rational, Refusal> {
                let xyz = [&point.x, &point.y, &point.z].map(|x| dec.enclose_at(x, BITS));
                let mut sum = Rational::zero();
                for (i, value) in xyz.into_iter().enumerate() {
                    sum = sum + (&value?.mid() - &points[0][i].mid()).abs();
                }
                Ok(sum)
            };
            let forward = distance(&start)? <= distance(&end)?;
            if row == BoundaryRow::Rail {
                along_u = Some((sense == Sense::Same) == forward);
            }
            rows.push((edge, row, forward));
            seeds.insert(
                edge,
                Seed {
                    points,
                    error: grid.chord.clone(),
                    forward,
                    exact: Vec::new(),
                },
            );
        }
        let along_u = along_u.ok_or_else(|| failure("a rolling-envelope face has no rail"))?;
        envelopes.insert(
            face,
            Envelope {
                grid,
                rows,
                along_u,
                ends: Ends::Meridians,
                rational: None,
            },
        );
    }
    Ok(Some(Junction { envelopes, seeds }))
}

// A band over a mitre: its rails are the grid's first and last rows, its
// end section the first column, and its last column collapses to the apex.
fn mitre_band(
    model: &Model,
    face: FaceId,
    support: &RollingMitre,
    tolerance: &Rational,
    options: &MeshOptions,
    seeds: &mut BTreeMap<EdgeId, Seed>,
    dec: &mut Decider,
) -> Result<Envelope, Refusal> {
    let checked = RollingMitreBand::of_face(model, face, dec)?
        .ok_or_else(|| failure("a rolling mitre band has no proved complete domain"))?;
    let grid = support.display_grid(tolerance, options.max_vertices, dec)?;
    let mut rows = Vec::with_capacity(3);
    let mut along_u = None;
    for (edge, row, sense, forward) in checked.rows() {
        let points = row_points(&grid, row);
        let count = points.len();
        let mut exact = Vec::with_capacity(count);
        for index in 0..count {
            let t = Rational::from_ratio(index as i64, count as i64 - 1).unwrap();
            let (u, v) = match row {
                BoundaryRow::Rail => (t, Rational::zero()),
                BoundaryRow::Mouth => (t, Rational::one()),
                BoundaryRow::Negative => (Rational::zero(), t),
                BoundaryRow::Positive => (Rational::one(), t),
            };
            let point = support.point(&u, &v);
            exact.push([point.x, point.y, point.z]);
        }
        if row == BoundaryRow::Rail {
            along_u = Some((sense == Sense::Same) == forward);
        }
        rows.push((edge, row, forward));
        seeds.insert(
            edge,
            Seed {
                points,
                error: grid.chord.clone(),
                forward,
                exact,
            },
        );
    }
    Ok(Envelope {
        grid,
        rows,
        along_u: along_u.ok_or_else(|| failure("a rolling mitre band has no first rail"))?,
        ends: Ends::ZeroEnd,
        rational: None,
    })
}

fn linear_band(
    model: &Model,
    face: FaceId,
    support: &LinearFillet,
    tolerance: &Rational,
    options: &MeshOptions,
    seeds: &mut BTreeMap<EdgeId, Seed>,
    dec: &mut Decider,
) -> Result<Envelope, Refusal> {
    let checked = LinearFilletBand::of_face(model, face, dec)?
        .ok_or_else(|| failure("variable band has no proved complete domain"))?;
    let grid = support.display_grid(tolerance, options.max_vertices, dec)?;
    let mut rows = Vec::with_capacity(4);
    let mut along_u = None;
    for (edge, row, sense, forward) in checked.rows() {
        let points = row_points(&grid, row);
        let count = points.len();
        let mut exact = Vec::with_capacity(count);
        for index in 0..count {
            let t = Rational::from_ratio(index as i64, count as i64 - 1).unwrap();
            let (arc, station) = match row {
                BoundaryRow::Rail => (Rational::zero(), t),
                BoundaryRow::Mouth => (Rational::one(), t),
                BoundaryRow::Negative => (t, Rational::zero()),
                BoundaryRow::Positive => (t, Rational::one()),
            };
            let point = support.point(&arc, &station)?;
            exact.push([point.x, point.y, point.z]);
        }
        if row == BoundaryRow::Rail {
            along_u = Some((sense == Sense::Same) == forward);
        }
        rows.push((edge, row, forward));
        seeds.insert(
            edge,
            Seed {
                points,
                error: grid.chord.clone(),
                forward,
                exact,
            },
        );
    }
    let ends = if support.zero_end(0) {
        Ends::ZeroStart
    } else if support.zero_end(1) {
        Ends::ZeroEnd
    } else {
        Ends::Meridians
    };
    Ok(Envelope {
        grid,
        rows,
        along_u: along_u.ok_or_else(|| failure("variable band has no first contact rail"))?,
        ends,
        rational: None,
    })
}

// A closed band's grid and its two rings. A ring has no end vertices, so its
// direction comes from the sense of its ellipse against the chart's turn.
fn band(
    model: &Model,
    face: FaceId,
    pipe: &RollingPipe,
    tolerance: &Rational,
    options: &MeshOptions,
    seeds: &mut BTreeMap<EdgeId, Seed>,
    dec: &mut Decider,
) -> Result<Envelope, Refusal> {
    let grid = pipe.display_grid(tolerance, options.max_vertices, dec)?;
    let rows = pipe_rows(model, face, dec)?
        .ok_or_else(|| failure("a pipe face has no complete band of two rails"))?;
    let turn = pipe.chart_normal();
    let last = Rational::from_i64(grid.columns as i64 - 1);
    let mut rings = Vec::with_capacity(2);
    for (index, (edge, sense)) in rows.into_iter().enumerate() {
        let definition = model.edge(edge);
        // A cone's wall rail is the carrier's last row and runs with the chart.
        let same = match model.geom().curve(definition.curve) {
            Curve::Ellipse(rail) => {
                let turns = dec.enclose_at(&rail.normal.dot(&turn), BITS)?;
                turns
                    .sign()
                    .ok_or_else(|| failure("a pipe rail's turn is not separated from zero"))?
                    == catallax_numbers::StrictSign::Positive
            }
            Curve::BSpline(_) => true,
            _ => return Err(failure("a pipe rail is not its ellipse or row")),
        };
        let forward = same == (definition.sense == Sense::Same);
        // Each span's parameter is cut in equal halves, so a column is an exact
        // fraction of the chart.
        let line = Rational::from_i64(index as i64);
        let exact = (0..grid.columns)
            .map(|column| {
                let u = Rational::from_i64(4 * column as i64) / last.clone();
                let p = pipe.point(&u, &line);
                [p.x, p.y, p.z]
            })
            .collect();
        rings.push(Ring {
            edge,
            sense,
            forward,
            exact,
        });
    }
    Ok(seed_band(grid, rings, seeds, Ends::Seam))
}

// A closed saddle's rings are the carrier's boundary rows, so each runs up the
// chart in its edge's sense.
fn saddle_band(
    model: &Model,
    face: FaceId,
    saddle: &RollingSaddle,
    tolerance: &Rational,
    options: &MeshOptions,
    seeds: &mut BTreeMap<EdgeId, Seed>,
    dec: &mut Decider,
) -> Result<Envelope, Refusal> {
    let grid = saddle.display_grid(tolerance, options.max_vertices, dec)?;
    let columns = saddle.display_columns(tolerance, options.max_vertices, dec)?;
    let rows = saddle_rows(model, face)
        .ok_or_else(|| failure("a saddle face has no complete band of two rails"))?;
    let mut rings = Vec::with_capacity(2);
    for (edge, sense) in rows {
        let definition = model.edge(edge);
        let Curve::BSpline(rail) = model.geom().curve(definition.curve) else {
            return Err(failure("a saddle rail is not its row"));
        };
        // Each exact point costs a rational evaluation of the rail, and only a
        // planar neighbour reads them.
        let mut exact = Vec::new();
        if meets_plane(model, edge) {
            exact.reserve(grid.columns);
            for u in &columns {
                let p = rail.point_at(&Real::from(u.clone()), dec)?;
                exact.push([p.x, p.y, p.z]);
            }
        }
        rings.push(Ring {
            edge,
            sense,
            forward: definition.sense == Sense::Same,
            exact,
        });
    }
    let ends = if matches!(
        saddle.definition().law().ring,
        catallax_shapes::SaddleRing::Crossing { .. }
    ) {
        Ends::Poles
    } else {
        Ends::Seam
    };
    Ok(seed_band(grid, rings, seeds, ends))
}

// Whether a planar face uses the edge.
fn meets_plane(model: &Model, edge: EdgeId) -> bool {
    model.coedges_of_edge(edge).into_iter().any(|coedge| {
        model.owner_of_coedge(coedge).is_some_and(|(_, face)| {
            matches!(
                model.geom().surface(model.face(face).surface),
                Surface::Plane(_)
            )
        })
    })
}

// One closed rail of a band, in its grid row's column order.
struct Ring {
    edge: EdgeId,
    sense: Sense,
    forward: bool,
    exact: Vec<[Real; 3]>,
}

// Seeds both rails, first on line zero and second on the last line.
fn seed_band(
    grid: EnvelopeGrid,
    rings: Vec<Ring>,
    seeds: &mut BTreeMap<EdgeId, Seed>,
    ends: Ends,
) -> Envelope {
    let mut rows = Vec::with_capacity(2);
    let mut along_u = false;
    for (index, ring) in rings.into_iter().enumerate() {
        let row = if index == 0 {
            BoundaryRow::Rail
        } else {
            BoundaryRow::Mouth
        };
        if index == 0 {
            along_u = (ring.sense == Sense::Same) == ring.forward;
        }
        rows.push((ring.edge, row, ring.forward));
        seeds.insert(
            ring.edge,
            Seed {
                points: row_points(&grid, row),
                error: grid.chord.clone(),
                forward: ring.forward,
                exact: ring.exact,
            },
        );
    }
    Envelope {
        grid,
        rows,
        along_u,
        ends,
        rational: None,
    }
}

impl Junction {
    pub(crate) fn contains(&self, face: FaceId) -> bool {
        self.envelopes.contains_key(&face)
    }

    pub(crate) fn seed(&self, edge: EdgeId) -> Option<&Seed> {
        self.seeds.get(&edge)
    }

    pub(crate) fn triangulate(
        &self,
        face: FaceId,
        boundaries: &mut Boundaries,
        options: &MeshOptions,
    ) -> Result<(), Refusal> {
        let envelope = &self.envelopes[&face];
        let grid = &envelope.grid;
        let (columns, rows) = (grid.columns, grid.rows);
        let at = |column: usize, line: usize| column * rows + line;
        let mut index = vec![usize::MAX; columns * rows];
        let mut error = grid.chord.clone();
        for (edge, row, forward) in &envelope.rows {
            let mesh = &boundaries.mesh.edges[boundaries.edges[edge]];
            error = error.max(mesh.max_chord_error.clone());
            let mut vertices = mesh.vertices.clone();
            if !forward {
                vertices.reverse();
            }
            for (i, vertex) in vertices.into_iter().enumerate() {
                let slot = match row {
                    BoundaryRow::Rail => at(i, 0),
                    BoundaryRow::Mouth => at(i, rows - 1),
                    BoundaryRow::Negative => at(0, i),
                    BoundaryRow::Positive => at(columns - 1, i),
                };
                index[slot] = vertex;
            }
        }
        for column in [0, columns - 1] {
            if envelope.ends.collapsed(column, columns) {
                let pole = index[at(column, 0)];
                if pole == usize::MAX || index[at(column, rows - 1)] != pole {
                    return Err(failure("a collapsed band does not share its rail-end pole"));
                }
                for line in 1..rows - 1 {
                    index[at(column, line)] = pole;
                }
            }
        }
        let first_column = usize::from(envelope.ends != Ends::Seam);
        for column in first_column..columns - 1 {
            for line in 1..rows - 1 {
                if boundaries.mesh.positions.len() >= options.max_vertices {
                    return Err(Refusal::new(
                        "mesh.resource-limit",
                        "envelope grid reached max_vertices",
                    ));
                }
                let point = &grid.points[at(column, line)];
                let position = chart::position(point)?;
                if envelope.rational.is_none() {
                    error = error.max(&grid.chord + &chart::residual(&position, point));
                }
                index[at(column, line)] = boundaries.mesh.positions.len();
                boundaries.mesh.positions.push(position);
            }
        }
        if envelope.ends == Ends::Seam {
            for line in 0..rows {
                index[at(columns - 1, line)] = index[at(0, line)];
            }
        }
        if let Some(proof) = &envelope.rational {
            let mut residual = Rational::zero();
            for (slot, &vertex) in index.iter().enumerate() {
                residual = residual.max(chart::residual(
                    &boundaries.mesh.positions[vertex],
                    &grid.points[slot],
                ));
            }
            if !proof.embedding.accepts_rounding(&residual) {
                return Err(Refusal::new(
                    "mesh.rational-patch-unproven",
                    "rounded grid has no positive whole-chart embedding margin",
                ));
            }
            error = error.max(&grid.chord + &residual);
        }
        if error > options.chord_tolerance {
            return Err(failure(
                "proved envelope mesh error exceeds the requested output accuracy",
            ));
        }
        let start = boundaries.mesh.triangles.len();
        let mut longest = Rational::zero();
        let mut aspect = Rational::zero();
        let mut angular = None::<Rational>;
        for column in 0..columns - 1 {
            for line in 0..rows - 1 {
                let a = index[at(column, line)];
                let c = index[at(column + 1, line)];
                let d = index[at(column, line + 1)];
                let e = index[at(column + 1, line + 1)];
                let pair = if envelope.along_u {
                    [[a, c, e], [a, e, d]]
                } else {
                    [[a, e, c], [a, d, e]]
                };
                for vertices in pair {
                    // One half of each pole cell has zero-area linear image.
                    // Its image is contained in the other half's fan triangle;
                    // the same whole-cell chord bound therefore still applies.
                    if (envelope.ends.collapsed(column, columns)
                        || envelope.ends.collapsed(column + 1, columns))
                        && (vertices[0] == vertices[1]
                            || vertices[1] == vertices[2]
                            || vertices[0] == vertices[2])
                    {
                        continue;
                    }
                    let (ratio, length, normal) = quality(&boundaries.mesh.positions, vertices)?;
                    if let Some(boxes) = envelope
                        .rational
                        .as_ref()
                        .and_then(|proof| proof.normal_boxes.as_ref())
                    {
                        let bound = rational_patch::angular(
                            &normal,
                            &boxes[column * (rows - 1) + line],
                            envelope.along_u,
                        )?;
                        if options
                            .angular_deviation
                            .as_ref()
                            .is_some_and(|limit| bound > *limit)
                        {
                            return Err(Refusal::new(
                                "mesh.rational-patch-unproven",
                                "rounded facet angular bound exceeds the requested limit",
                            ));
                        }
                        angular =
                            Some(angular.map_or(bound.clone(), |previous| previous.max(bound)));
                    }
                    aspect = aspect.max(ratio);
                    longest = longest.max(length);
                    boundaries.mesh.triangles.push(Triangle { vertices, face });
                }
            }
        }
        let mut unavailable = |control, reason| {
            boundaries
                .mesh
                .unavailable_controls
                .push(UnavailableControl {
                    face,
                    control,
                    reason,
                })
        };
        if options.angular_deviation.is_some() && angular.is_none() {
            unavailable(
                MeshControl::AngularDeviation,
                "the envelope grid certifies position, not normal variation",
            );
        }
        if options
            .max_edge_length
            .as_ref()
            .is_some_and(|limit| longest > limit * limit)
        {
            unavailable(
                MeshControl::EdgeLength,
                "the envelope grid is refined for chord error only",
            );
        }
        if options
            .max_aspect_ratio
            .as_ref()
            .is_some_and(|limit| aspect > *limit)
        {
            unavailable(
                MeshControl::AspectRatio,
                "the certified aspect-ratio bound exceeds the requested limit",
            );
        }
        boundaries.mesh.faces.push(FaceMesh {
            face,
            triangles: start..boundaries.mesh.triangles.len(),
            max_chord_error: error,
            max_aspect_ratio: Some(aspect),
            max_angular_deviation: angular,
        });
        Ok(())
    }
}

// Longest edge over shortest altitude, and the longest squared edge.
fn quality(
    positions: &[[f64; 3]],
    triangle: [usize; 3],
) -> Result<(Rational, Rational, [Ball; 3]), Refusal> {
    let points =
        triangle.map(|index| positions[index].map(|x| Ball::from_rational(&chart::q(x), BITS)));
    let a: [Ball; 3] = std::array::from_fn(|i| &points[1][i] - &points[0][i]);
    let c: [Ball; 3] = std::array::from_fn(|i| &points[2][i] - &points[0][i]);
    let normal = [
        &a[1] * &c[2] - &a[2] * &c[1],
        &a[2] * &c[0] - &a[0] * &c[2],
        &a[0] * &c[1] - &a[1] * &c[0],
    ];
    let area = (&(&normal[0].square() + &normal[1].square()) + &normal[2].square())
        .sqrt()
        .ok_or_else(|| failure("rounded triangle area is not enclosed"))?;
    if area.lo() <= Rational::zero() {
        return Err(failure(
            "rounded triangle quality exceeds certificate precision",
        ));
    }
    let longest = (0..3)
        .map(|i| {
            crate::controls::length_squared(
                positions[triangle[i]],
                positions[triangle[(i + 1) % 3]],
            )
        })
        .max()
        .unwrap();
    Ok((&longest / &area.lo(), longest, normal))
}
