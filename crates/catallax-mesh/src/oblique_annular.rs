//! A common tensor grid for both walls and the annular cap strips.
use super::*;
use catallax_numbers::Real;

pub(super) fn build(
    model: &Model,
    body: BodyId,
    checked: DraftedCylinder<'_>,
    options: &MeshOptions,
    dec: &mut Decider,
) -> Answered<Mesh> {
    let (inner, inner_face, inner_rings) = checked.hole().expect("checked annular body");
    let supports = [checked.support(), inner];
    let clearance = checked.section_clearance(dec)?.unwrap();
    let clearance = dec.enclose_at(&clearance, chart::BITS)?.lo();
    let tolerance = options
        .chord_tolerance
        .mul_pow2(-2)
        .min(clearance.mul_pow2(-3));
    let grids =
        supports.map(|support| support.display_grid(&tolerance, options.max_vertices / 2, dec));
    let [outer, inner] = grids;
    let (outer, inner) = (outer?, inner?);
    let columns = outer.columns.max(inner.columns) - 1;
    let rows = outer.rows.max(inner.rows);
    let count = columns
        .checked_mul(rows)
        .filter(|n| n.saturating_mul(2) <= options.max_vertices)
        .ok_or_else(|| failure("annular ruled grid exceeds the vertex limit"))?;
    let mut mesh = Mesh::default();
    let mut rounding = Rational::zero();
    for support in supports {
        let start = mesh.positions.len();
        for i in 0..=columns {
            dec.check_deadline()?;
            let u = Rational::from_ratio((4 * i) as i64, columns as i64).unwrap();
            for j in 0..rows {
                let v = Rational::from_ratio(j as i64, (rows - 1) as i64).unwrap();
                let p = support.chart_point(&u, &v)?;
                let enclosure = [
                    dec.enclose_at(&p.x, chart::BITS)?,
                    dec.enclose_at(&p.y, chart::BITS)?,
                    dec.enclose_at(&p.z, chart::BITS)?,
                ];
                if i == columns {
                    rounding =
                        rounding.max(chart::residual(&mesh.positions[start + j], &enclosure));
                } else {
                    let position = chart::position(&enclosure)?;
                    rounding = rounding.max(chart::residual(&position, &enclosure));
                    mesh.positions.push(position);
                }
            }
        }
    }
    let error = outer.chord.max(inner.chord) + rounding;
    if error > options.chord_tolerance || error.mul_pow2(1) >= clearance {
        return Err(failure(
            "rounded annular grid exceeds its chord or wall clearance bound",
        ));
    }
    let index = |wall: usize, i: usize, j: usize| wall * count + (i % columns) * rows + j;
    let mut grouped: [Vec<[usize; 3]>; 4] = std::array::from_fn(|_| Vec::new());
    for i in 0..columns {
        for cap in 0..2 {
            let j = if cap == 0 { 0 } else { rows - 1 };
            let a = index(0, i, j);
            let b = index(0, i + 1, j);
            let c = index(1, i, j);
            let d = index(1, i + 1, j);
            let triangles = if cap == 0 {
                [[a, c, d], [a, d, b]]
            } else {
                [[a, b, d], [a, d, c]]
            };
            grouped[cap].extend(triangles);
        }
        for wall in 0..2 {
            for j in 0..rows - 1 {
                let a = index(wall, i, j);
                let b = index(wall, i + 1, j);
                let c = index(wall, i, j + 1);
                let d = index(wall, i + 1, j + 1);
                let triangles = if wall == 0 {
                    [[a, b, d], [a, d, c]]
                } else {
                    [[a, d, b], [a, c, d]]
                };
                grouped[2 + wall].extend(triangles);
            }
        }
    }
    let [lower, upper, outer_face] = checked.faces();
    let faces = [lower, upper, outer_face, inner_face];
    for face in model.faces_of_body(body) {
        let slot = faces.iter().position(|f| *f == face).ok_or_else(|| {
            failure("annular display face inventory differs from its checked body")
        })?;
        let start = mesh.triangles.len();
        let mut aspect = Rational::zero();
        for (t, vertices) in grouped[slot].iter().copied().enumerate() {
            let normal = if slot < 2 {
                supports[0]
                    .axis()
                    .scaled(&Real::from_i64(if slot == 0 { -1 } else { 1 }))
            } else {
                let cell = t / 2;
                let i = cell / (rows - 1);
                let j = cell % (rows - 1);
                let u = Rational::from_ratio((4 * i + 2) as i64, columns as i64).unwrap();
                let v = Rational::from_ratio((2 * j + 1) as i64, (2 * (rows - 1)) as i64).unwrap();
                let support = supports[slot - 2];
                support
                    .normal_at(&support.chart_point(&u, &v)?, dec)?
                    .scaled(&Real::from_i64(if slot == 2 { 1 } else { -1 }))
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
                    "annular display needs interior subdivision to meet the requested edge length",
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
                reason: "the certified annular-grid aspect bound exceeds the requested limit",
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
    for (wall, rings) in [checked.rings(), inner_rings].into_iter().enumerate() {
        for (cap, edge) in rings.into_iter().enumerate() {
            let mut vertices: Vec<_> = (0..=columns)
                .map(|i| index(wall, i, if cap == 0 { 0 } else { rows - 1 }))
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
    }
    crate::triangulate::check_closed(&mesh)?;
    Ok(Answer::with_marks(mesh, model.body(body).marks.clone()))
}
