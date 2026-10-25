//! Certified affine images of an upright field prism. This covers finite signed
//! algebraic extrusions without widening rational charts or consulting B-rep.
use super::{Diagnostics, Mesh, Options, failure, mesh_plan_inner, q};
use crate::{Instruction, Plan, geometry};
use catallax_numbers::{Ball, Decider, Rational, Refusal, Sign};
use catallax_recipe::Role;
use catallax_shapes::RPoint3;

fn source_role(mut role: &Role) -> &Role {
    while let Role::Copy { of, .. } = role {
        role = &of.role;
    }
    role
}

pub(super) fn mesh(
    plan: &Plan,
    options: &Options,
    dec: &mut Decider,
    report: Option<&mut Diagnostics>,
) -> Result<Option<Mesh>, Refusal> {
    let Instruction::Leaf(leaf) = plan.instructions[plan.root] else {
        return Ok(None);
    };
    let geometry::Leaf::Feature(feature) = &plan.leaves[leaf] else {
        return Ok(None);
    };
    // Historical Boolean ownership gates may refer to other, untransported
    // definitions. Restrict this chart transport to complete source domains.
    if plan
        .faces
        .iter()
        .any(|f| f.in_result && f.leaf == leaf && f.domain != 0)
    {
        return Ok(None);
    }
    let Some(height) = feature.algebraic_display_height(dec)? else {
        return Ok(None);
    };
    let Some(first) = plan.faces.iter().find(|f| f.in_result && f.leaf == leaf) else {
        return Ok(None);
    };
    let Some((definition, faces)) =
        feature.algebraic_display_source(&first.origin.step, &height, dec)?
    else {
        return Ok(None);
    };
    // Preserve registry handles, origins, domain ownership and the result tape.
    // Only the reference geometry changes; the returned mesh is then mapped by
    // the exact affine definition back into the actual result.
    let mut source = plan.clone();
    source.leaves[leaf] = geometry::Leaf::Feature(definition);
    for (index, face) in source
        .faces
        .iter_mut()
        .enumerate()
        .filter(|(_, f)| f.leaf == leaf)
    {
        let role = source_role(plan.definition_role(index));
        let reference = faces
            .iter()
            .find(|f| &f.origin.role == role)
            .ok_or_else(|| {
                failure(
                    "field.display-unsupported",
                    "extrusion reference profile has no matching source face",
                )
            })?;
        geometry::face_geometry(face, reference.surface.clone(), dec)?;
    }
    let lipschitz = feature.algebraic_display_lipschitz(&height, dec)?;
    let mut controls = options.clone();
    controls.tolerance = &options.tolerance / &(&q(2) * &lipschitz);
    let mut mesh = mesh_plan_inner(&source, &controls, dec, report)?;
    // Enclose the affine coefficients once. Every output vertex and facet
    // reuses these exact coefficient enclosures; no per-vertex root isolation.
    let (origin, matrix) = feature.algebraic_display_affine(&height);
    let origin = [&origin.x, &origin.y, &origin.z]
        .map(|x| geometry::enclose(x, 192, dec))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;
    let mut columns = Vec::new();
    for column in matrix {
        columns.push(
            [&column.x, &column.y, &column.z]
                .map(|x| geometry::enclose(x, 192, dec))
                .into_iter()
                .collect::<Result<Vec<_>, _>>()?,
        );
    }
    let mut rounding = q(0);
    let mut enclosed = Vec::with_capacity(mesh.positions.len());
    for position in &mut mesh.positions {
        dec.check_deadline()?;
        let source = position.map(|x| {
            Ball::from_rational(
                &Rational::from_f64_exact(x).expect("certified finite source coordinate"),
                192,
            )
        });
        let mapped: [Ball; 3] = std::array::from_fn(|i| {
            (0..3).fold(origin[i].clone(), |sum, j| {
                sum + &columns[j][i] * &source[j]
            })
        });
        let mut error = q(0);
        for axis in 0..3 {
            let b = &mapped[axis];
            let rounded = b.mid().to_f64_nearest();
            let value = Rational::from_f64_exact(rounded).ok_or_else(|| {
                failure(
                    "field.display-precision",
                    "non-finite affine extrusion display coordinate",
                )
            })?;
            error = error + (&b.lo() - &value).abs().max((&b.hi() - &value).abs());
            position[axis] = if rounded == 0.0 { 0.0 } else { rounded };
        }
        rounding = rounding.max(error);
        enclosed.push(mapped);
    }
    let deviation = &lipschitz * &mesh.max_deviation + rounding;
    if deviation > options.tolerance {
        return Err(failure(
            "field.display-precision",
            "affine extrusion rounding exceeds requested display tolerance",
        ));
    }
    let reversing = feature.algebraic_display_reversing();
    for triangle in &mut mesh.triangles {
        dec.check_deadline()?;
        let [a, b, c] = triangle.vertices;
        // Affine maps commute with triangle interpolation. Check that final
        // coordinate rounding also preserves each nondegenerate facet's sense.
        let u: [Ball; 3] = std::array::from_fn(|i| &enclosed[b][i] - &enclosed[a][i]);
        let v: [Ball; 3] = std::array::from_fn(|i| &enclosed[c][i] - &enclosed[a][i]);
        let normal: [Ball; 3] = std::array::from_fn(|i| {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            &u[j] * &v[k] - &u[k] * &v[j]
        });
        let rounded = [a, b, c].map(|i| {
            let [x, y, z] = mesh.positions[i].map(|x| Rational::from_f64_exact(x).unwrap());
            RPoint3::new(x, y, z)
        });
        let rounded_normal = rounded[1]
            .minus(&rounded[0])
            .cross(&rounded[2].minus(&rounded[0]));
        let mut dot = Ball::from_rational(&q(0), 192);
        for (n, c) in normal
            .iter()
            .zip([rounded_normal.x, rounded_normal.y, rounded_normal.z])
        {
            dot = dot + n * &Ball::from_rational(&c, 192);
        }
        if dot.lo().sign() != Sign::Positive {
            return Err(failure(
                "field.display-precision",
                "affine extrusion facet orientation is unresolved after rounding",
            ));
        }
        if reversing {
            triangle.vertices.swap(1, 2);
        }
    }
    mesh.max_deviation = deviation;
    Ok(Some(mesh))
}
