//! One shared certified grid for an exact ruled side and its two cap disks.
use crate::{EdgeMesh, FaceMesh, Mesh, MeshOptions, Triangle, chart, failure};
use catallax_numbers::{Answer, Answered, Decider, Rational, Refusal};
use catallax_topology::{BodyId, DraftedCylinder, Model};

pub(crate) fn build(
    model: &Model,
    body: BodyId,
    options: &MeshOptions,
    dec: &mut Decider,
) -> Answered<Mesh> {
    if options.angular_deviation.is_some() {
        return Err(Refusal::new(
            "mesh.oblique-draft-angular-coverage",
            "whole-ruled-face normal variation is not certified for the requested angular control",
        ));
    }
    let checked = DraftedCylinder::check(model, body, dec)?;
    let support = checked.support();
    let grid = support.display_grid(
        &options.chord_tolerance.mul_pow2(-2),
        options.max_vertices.saturating_sub(2),
        dec,
    )?;
    let columns = grid.columns - 1;
    let rows = grid.rows;
    if columns.saturating_mul(rows).saturating_add(2) > options.max_vertices {
        return Err(failure("oblique draft display exceeds the vertex limit"));
    }
    let mut mesh = Mesh::default();
    let mut rounding = Rational::zero();
    for p in grid.points.iter().take(columns * rows) {
        let position = chart::position(p)?;
        rounding = rounding.max(chart::residual(&position, p));
        mesh.positions.push(position);
    }
    // The last chart column is the first physical column. Its independently
    // enclosed seam must fit the same rounded coordinates, not a second wire.
    for j in 0..rows {
        rounding = rounding.max(chart::residual(
            &mesh.positions[j],
            &grid.points[columns * rows + j],
        ));
    }
    let mut centers = [0, 0];
    for (center, height) in centers.iter_mut().zip(support.heights()) {
        let p = support.origin().plus(&support.axis().scaled(height));
        let enclosure = [
            dec.enclose_at(&p.x, chart::BITS)?,
            dec.enclose_at(&p.y, chart::BITS)?,
            dec.enclose_at(&p.z, chart::BITS)?,
        ];
        let position = chart::position(&enclosure)?;
        rounding = rounding.max(chart::residual(&position, &enclosure));
        *center = mesh.positions.len();
        mesh.positions.push(position);
    }
    let error = &grid.chord + &rounding;
    if error > options.chord_tolerance {
        return Err(failure(
            "rounded oblique draft grid exceeds the requested chord tolerance",
        ));
    }
    let index = |i: usize, j: usize| (i % columns) * rows + j;
    let faces = checked.faces();
    let mut grouped: [Vec<[usize; 3]>; 3] = std::array::from_fn(|_| Vec::new());
    for i in 0..columns {
        grouped[0].push([centers[0], index(i + 1, 0), index(i, 0)]);
        grouped[1].push([centers[1], index(i, rows - 1), index(i + 1, rows - 1)]);
        for j in 0..rows - 1 {
            grouped[2].push([index(i, j), index(i + 1, j), index(i + 1, j + 1)]);
            grouped[2].push([index(i, j), index(i + 1, j + 1), index(i, j + 1)]);
        }
    }
    for face in model.faces_of_body(body) {
        let slot = faces.iter().position(|f| *f == face).ok_or_else(|| {
            failure("ruled display face inventory differs from its complete cell")
        })?;
        let start = mesh.triangles.len();
        let mut aspect = Rational::zero();
        for (t, vertices) in grouped[slot].iter().copied().enumerate() {
            let normal = if slot < 2 {
                support
                    .axis()
                    .scaled(&catallax_numbers::Real::from_i64(if slot == 0 {
                        -1
                    } else {
                        1
                    }))
            } else {
                let cell = t / 2;
                let i = cell / (rows - 1);
                let j = cell % (rows - 1);
                let u = Rational::from_ratio((8 * i + 4) as i64, columns as i64)
                    .unwrap()
                    .mul_pow2(-1);
                let v = Rational::from_ratio((2 * j + 1) as i64, (2 * (rows - 1)) as i64).unwrap();
                support.normal_at(&support.chart_point(&u, &v)?, dec)?
            };
            let n = [
                dec.enclose_at(&normal.x, chart::BITS)?,
                dec.enclose_at(&normal.y, chart::BITS)?,
                dec.enclose_at(&normal.z, chart::BITS)?,
            ];
            let (ratio, longest) = crate::spline_prism::oriented(&mesh.positions, vertices, &n)?;
            aspect = aspect.max(ratio);
            if options
                .max_edge_length
                .as_ref()
                .is_some_and(|limit| longest > limit * limit)
            {
                return Err(Refusal::new(
                    "mesh.oblique-draft-edge-length",
                    "ruled display needs interior subdivision to meet the requested edge length",
                ));
            }
            mesh.triangles.push(Triangle { vertices, face });
        }
        if options
            .max_aspect_ratio
            .as_ref()
            .is_some_and(|limit| aspect > *limit)
        {
            mesh.unavailable_controls.push(crate::UnavailableControl {
                face,
                control: crate::MeshControl::AspectRatio,
                reason: "the certified ruled-grid aspect bound exceeds the requested limit",
            });
        }
        mesh.faces.push(FaceMesh {
            face,
            triangles: start..mesh.triangles.len(),
            max_chord_error: error.clone(),
            max_aspect_ratio: Some(aspect),
            max_angular_deviation: None,
        });
    }
    for (i, edge) in checked.rings().into_iter().enumerate() {
        let mut vertices: Vec<_> = (0..=columns)
            .map(|j| index(j, if i == 0 { 0 } else { rows - 1 }))
            .collect();
        if model.edge(edge).sense == catallax_topology::Sense::Reversed {
            vertices.reverse();
        }
        mesh.edges.push(EdgeMesh {
            edge,
            vertices,
            max_chord_error: error.clone(),
        });
    }
    crate::triangulate::check_closed(&mesh)?;
    Ok(Answer::with_marks(mesh, model.body(body).marks.clone()))
}

pub(crate) fn preview(
    model: &Model,
    body: BodyId,
    options: &crate::PreviewOptions,
    dec: &mut Decider,
) -> Result<crate::PreviewMesh, Refusal> {
    // Preview angular targets are not certificates. Bound the chord grid and
    // discard its certificates rather than promote an unproved normal target.
    let mesh = build(
        model,
        body,
        &MeshOptions {
            chord_tolerance: options.chord_tolerance.clone(),
            max_vertices: options.max_vertices,
            ..MeshOptions::default()
        },
        dec,
    )?
    .value;
    Ok(crate::PreviewMesh {
        positions: mesh.positions,
        triangles: mesh.triangles,
        faces: mesh
            .faces
            .into_iter()
            .map(|f| crate::PreviewFace {
                face: f.face,
                triangles: f.triangles,
            })
            .collect(),
        edges: mesh
            .edges
            .into_iter()
            .map(|e| crate::PreviewEdge {
                edge: e.edge,
                vertices: e.vertices,
            })
            .collect(),
    })
}
