//! Display of a body with rolling-envelope faces.
//!
//! No chart covers a fitted carrier, so each envelope face takes its structured
//! display grid, whose chord bound comes from subdivided controls. The grid's
//! boundary rows seed the shared edge polylines, which keeps one mesh across
//! every boundary: the analytic neighbours triangulate up to those rows.

use std::collections::BTreeMap;

use catallax_numbers::{Ball, Decider, Rational, Real, Refusal};
use catallax_shapes::{
    BoundaryRow, Curve, EnvelopeGrid, RollingPipe, RollingSaddle, SplineSupport, Surface,
};
use catallax_topology::{
    BodyId, EdgeId, FaceId, Model, Sense, envelope_rows, pipe_rows, saddle_rows,
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
    /// A closed rail's exact points. A planar face needs them to triangulate.
    pub(crate) exact: Vec<[Real; 3]>,
}

struct Envelope {
    grid: EnvelopeGrid,
    rows: Vec<(EdgeId, BoundaryRow, bool)>,
    // Whether the face boundary walks the rail along increasing chart u.
    along_u: bool,
    // A closed band: the last column is the seam again, and there are no
    // meridian rows to supply its interior points.
    closed: bool,
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

/// Whether every spline face of the body is a rolling envelope.
pub(crate) fn present(model: &Model, body: BodyId) -> bool {
    let mut found = false;
    for face in model.faces_of_body(body) {
        match model.geom().surface(model.face(face).surface) {
            Surface::BSpline(
                SplineSupport::RollingEnvelope(_)
                | SplineSupport::RollingPipe(_)
                | SplineSupport::RollingSaddle(_),
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
    for &face in faces {
        match model.geom().surface(model.face(face).surface) {
            Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) => {
                found.push((face, envelope))
            }
            Surface::BSpline(SplineSupport::RollingPipe(pipe)) => bands.push((face, pipe)),
            Surface::BSpline(SplineSupport::RollingSaddle(saddle)) => saddles.push((face, saddle)),
            Surface::BSpline(_) => return Ok(None),
            _ => {}
        }
    }
    if found.is_empty() && bands.is_empty() && saddles.is_empty() {
        return Ok(None);
    }
    let tolerance = options.chord_tolerance.mul_pow2(-1);
    let mut envelopes = BTreeMap::new();
    let mut seeds = BTreeMap::new();
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
                closed: false,
            },
        );
    }
    Ok(Some(Junction { envelopes, seeds }))
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
        let Curve::Ellipse(rail) = model.geom().curve(definition.curve) else {
            return Err(failure("a pipe rail is not its ellipse"));
        };
        let turns = dec.enclose_at(&rail.normal.dot(&turn), BITS)?;
        let same = turns
            .sign()
            .ok_or_else(|| failure("a pipe rail's turn is not separated from zero"))?
            == catallax_numbers::StrictSign::Positive;
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
    Ok(closed(grid, rings, seeds))
}

// A closed saddle's rings are the carrier's boundary rows, so each runs up the
// chart in its edge's sense. Columns are equal cuts of `[-1, 3]`.
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
    let rows = saddle_rows(model, face)
        .ok_or_else(|| failure("a saddle face has no complete band of two rails"))?;
    let last = Rational::from_i64(grid.columns as i64 - 1);
    let mut rings = Vec::with_capacity(2);
    for (edge, sense) in rows {
        let definition = model.edge(edge);
        let Curve::BSpline(rail) = model.geom().curve(definition.curve) else {
            return Err(failure("a saddle rail is not its row"));
        };
        let mut exact = Vec::with_capacity(grid.columns);
        for column in 0..grid.columns {
            let u = Rational::from_i64(4 * column as i64) / last.clone() - Rational::one();
            let p = rail.point_at(&Real::from(u), dec)?;
            exact.push([p.x, p.y, p.z]);
        }
        rings.push(Ring {
            edge,
            sense,
            forward: definition.sense == Sense::Same,
            exact,
        });
    }
    Ok(closed(grid, rings, seeds))
}

// One closed rail of a band, in its grid row's column order.
struct Ring {
    edge: EdgeId,
    sense: Sense,
    forward: bool,
    exact: Vec<[Real; 3]>,
}

// Seeds both rails, first on line zero and second on the last line.
fn closed(grid: EnvelopeGrid, rings: Vec<Ring>, seeds: &mut BTreeMap<EdgeId, Seed>) -> Envelope {
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
        closed: true,
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
        let first_column = usize::from(!envelope.closed);
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
                error = error.max(&grid.chord + &chart::residual(&position, point));
                index[at(column, line)] = boundaries.mesh.positions.len();
                boundaries.mesh.positions.push(position);
            }
        }
        if envelope.closed {
            for line in 0..rows {
                index[at(columns - 1, line)] = index[at(0, line)];
            }
        }
        if error > options.chord_tolerance {
            return Err(failure(
                "proved envelope mesh error exceeds the requested output accuracy",
            ));
        }
        let start = boundaries.mesh.triangles.len();
        let mut longest = Rational::zero();
        let mut aspect = Rational::zero();
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
                    let (ratio, length) = quality(&boundaries.mesh.positions, vertices)?;
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
        if options.angular_deviation.is_some() {
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
            max_angular_deviation: None,
        });
        Ok(())
    }
}

// Longest edge over shortest altitude, and the longest squared edge.
fn quality(positions: &[[f64; 3]], triangle: [usize; 3]) -> Result<(Rational, Rational), Refusal> {
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
    Ok((&longest / &area.lo(), longest))
}
