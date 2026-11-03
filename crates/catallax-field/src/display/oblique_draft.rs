//! Closed display directly from the checked oblique draft definition.
use super::*;
use crate::{draft::Cell as DraftCell, geometry::Leaf};

fn rounded(p: &[Ball; 3]) -> Result<([f64; 3], Rational), Refusal> {
    let position = p.each_ref().map(|x| x.mid().to_f64_nearest());
    if !position.iter().all(|x| x.is_finite()) {
        return Err(failure(
            "field.display-coordinate",
            "oblique draft coordinate is not finite binary64",
        ));
    }
    let error = (0..3).fold(Rational::zero(), |sum, i| {
        let q = Rational::from_f64_exact(position[i]).unwrap();
        sum + (&q - &p[i].lo()).abs().max((&q - &p[i].hi()).abs())
    });
    Ok((position, error))
}

pub(super) fn mesh(
    plan: &Plan,
    options: &Options,
    dec: &mut Decider,
) -> Result<Option<Mesh>, Refusal> {
    let Instruction::Leaf(leaf) = plan.instructions[plan.root] else {
        return Ok(None);
    };
    let Leaf::Draft(region) = &plan.leaves[leaf] else {
        return Ok(None);
    };
    let DraftCell::Oblique(cell) = &region.cell else {
        return Ok(None);
    };
    let support = &cell.support;
    let mut faces = [None, None, None];
    for (index, face) in plan.faces.iter().enumerate() {
        if face.leaf != leaf || !face.in_result {
            continue;
        }
        let Some(origin) = region.source_origin(&face.origin) else {
            continue;
        };
        let slot = if cell.origins.contains(origin) {
            Some(2)
        } else {
            cell.caps
                .iter()
                .position(|(_, names)| names.contains(origin))
        };
        if let Some(slot) = slot {
            faces[slot].get_or_insert(FaceId(index));
        }
    }
    let [Some(lower), Some(upper), Some(side)] = faces else {
        return Err(failure(
            "field.display-origin",
            "oblique draft is missing a complete named boundary",
        ));
    };
    let capacity = options
        .max_vertices
        .saturating_sub(2)
        .min(options.max_cells / 2);
    let grid = support.display_grid(&options.tolerance.mul_pow2(-2), capacity, dec)?;
    let columns = grid.columns - 1;
    let rows = grid.rows;
    let rectangles = columns.saturating_mul(rows - 1);
    let visited = rectangles;
    if visited > options.max_cells
        || columns.saturating_mul(rows).saturating_add(2) > options.max_vertices
    {
        return Err(failure(
            "field.display-resource-limit",
            "oblique draft grid exceeds the fixed display allowance",
        ));
    }
    let mut positions = Vec::with_capacity(columns * rows + 2);
    let mut rounding = Rational::zero();
    for p in grid.points.iter().take(columns * rows) {
        let (position, error) = rounded(p)?;
        positions.push(position);
        rounding = rounding.max(error);
    }
    for (position, p) in positions
        .iter()
        .take(rows)
        .zip(&grid.points[columns * rows..])
    {
        let error = (0..3).fold(Rational::zero(), |sum, i| {
            let q = Rational::from_f64_exact(position[i]).unwrap();
            sum + (&q - &p[i].lo()).abs().max((&q - &p[i].hi()).abs())
        });
        rounding = rounding.max(error);
    }
    let mut centers = [0, 0];
    for (center, height) in centers.iter_mut().zip(support.heights()) {
        let p = support.origin().plus(&support.axis().scaled(height));
        let enclosed = [
            dec.enclose_at(&p.x, 128)?,
            dec.enclose_at(&p.y, 128)?,
            dec.enclose_at(&p.z, 128)?,
        ];
        let (position, error) = rounded(&enclosed)?;
        rounding = rounding.max(error);
        *center = positions.len();
        positions.push(position);
    }
    let max_deviation = &grid.chord + &rounding;
    if max_deviation > options.tolerance {
        return Err(failure(
            "field.display-accuracy",
            "rounded oblique draft display exceeds the requested tolerance",
        ));
    }
    let at = |i: usize, j: usize| (i % columns) * rows + j;
    let mut triangles = Vec::with_capacity(2 * rectangles + 2 * columns);
    for i in 0..columns {
        triangles.push(Triangle {
            vertices: [centers[0], at(i + 1, 0), at(i, 0)],
            face: lower,
        });
        triangles.push(Triangle {
            vertices: [centers[1], at(i, rows - 1), at(i + 1, rows - 1)],
            face: upper,
        });
        for j in 0..rows - 1 {
            triangles.push(Triangle {
                vertices: [at(i, j), at(i + 1, j), at(i + 1, j + 1)],
                face: side,
            });
            triangles.push(Triangle {
                vertices: [at(i, j), at(i + 1, j + 1), at(i, j + 1)],
                face: side,
            });
        }
    }
    // The complete tensor cover and convex cap sections supply the two-way
    // bound. Verify that binary64 rounding did not collapse or reverse a facet.
    // A shared index map, including the periodic column, closes every wire.
    for (triangle_index, triangle) in triangles.iter().enumerate() {
        dec.check_deadline()?;
        let p = triangle
            .vertices
            .map(|i| positions[i].map(|x| Rational::from_f64_exact(x).unwrap()));
        let a: [Rational; 3] = std::array::from_fn(|i| &p[1][i] - &p[0][i]);
        let b: [Rational; 3] = std::array::from_fn(|i| &p[2][i] - &p[0][i]);
        let n: [Rational; 3] = std::array::from_fn(|i| {
            &(&a[(i + 1) % 3] * &b[(i + 2) % 3]) - &(&a[(i + 2) % 3] * &b[(i + 1) % 3])
        });
        let normal = if triangle.face == lower {
            support.axis().scaled(&catallax_numbers::Real::from_i64(-1))
        } else if triangle.face == upper {
            support.axis().clone()
        } else {
            let i = triangle_index / (2 * rows);
            let j = (triangle_index % (2 * rows) - 2) / 2;
            let u = Rational::from_ratio((4 * i + 2) as i64, columns as i64).unwrap();
            let v = Rational::from_ratio((2 * j + 1) as i64, (2 * (rows - 1)) as i64).unwrap();
            support.normal_at(&support.chart_point(&u, &v)?, dec)?
        };
        let actual = catallax_shapes::Vector3::new(
            n[0].clone().into(),
            n[1].clone().into(),
            n[2].clone().into(),
        );
        if crate::geometry::sign(&actual.dot(&normal), dec)? != catallax_numbers::Sign::Positive {
            return Err(failure(
                "field.display-degenerate",
                "rounded oblique draft facet has no proved outward orientation",
            ));
        }
    }
    Ok(Some(Mesh {
        positions,
        triangles,
        max_deviation,
        cells: visited,
    }))
}
