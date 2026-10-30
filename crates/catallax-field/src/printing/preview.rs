//! Sampled previews never participate in an export certificate.
use super::*;

pub(crate) fn mesh(
    plan: &Plan,
    options: &Options,
    dec: &mut Decider,
) -> Result<Option<Mesh>, Refusal> {
    // This sampling pitch resolves the supported constant cubic beams and one
    // eighth of a gyroid period. Other fields keep their existing interval-driven
    // preview route.
    let mut feature: Option<Rational> = None;
    for leaf in &plan.leaves {
        let size = match leaf {
            geometry::Leaf::Lattice(lattice) => {
                let Some((_, _, r)) = lattice.axis_cubic_grid() else {
                    return Ok(None);
                };
                r
            }
            geometry::Leaf::Gyroid(gyroid) => gyroid.period().mul_pow2(-3),
            _ => continue,
        };
        feature = Some(feature.map_or(size.clone(), |previous| previous.min(size)));
    }
    let Some(feature) = feature else {
        return Ok(None);
    };
    if options.tolerance <= q(0) || options.max_cells == 0 || options.max_vertices < 4 {
        return Err(refusal(
            "field.print-options",
            "positive tolerance and limits required",
        ));
    }
    let target = options.tolerance.clone().min(feature);
    let mut edge = q(1);
    while edge > target {
        edge = edge.mul_pow2(-1);
    }
    while edge.mul_pow2(1) <= target {
        edge = edge.mul_pow2(1);
    }
    let bounds = bounds(plan, dec)?;
    let origin: [Rational; 3] = std::array::from_fn(|i| {
        &(&edge * &Rational::from_integer((&bounds[i][0] / &edge).floor())) - &edge
    });
    let mut counts = [0u32; 3];
    let mut cells = 1usize;
    for i in 0..3 {
        let count = (&(&bounds[i][1] - &origin[i]) / &edge)
            .to_f64_nearest()
            .ceil()
            + 1.0;
        if !count.is_finite() || count > (u32::MAX - 1) as f64 {
            return Err(refusal(
                "field.print-grid",
                "preview grid index range exhausted",
            ));
        }
        counts[i] = count as u32;
        cells = cells
            .checked_mul(counts[i] as usize)
            .ok_or_else(|| refusal("field.print-budget", "preview cell limit reached"))?;
    }
    if cells > options.max_cells {
        return Err(refusal("field.print-budget", "preview cell limit reached"));
    }
    let mut builder = Builder {
        classifier: Classifier::for_print(plan, dec)?,
        workspace: Workspace::default(),
        origin,
        edge,
        corner_signs: BTreeMap::new(),
        edge_vertices: BTreeMap::new(),
        grid_vertices: Vec::new(),
        mesh: Mesh {
            positions: vec![],
            triangles: vec![],
            max_deviation: q(0),
            cells,
        },
        options,
        root_steps: 0,
    };
    let mut workspace = Workspace::default();
    for x in 0..=counts[0] {
        for y in 0..=counts[1] {
            for z in 0..=counts[2] {
                charge(dec)?;
                let index = [x, y, z];
                let point = builder.point(index).map(|v| v.to_f64_nearest());
                let inside = plan.fast(point, &mut workspace)?.value < 0.0;
                builder.corner_signs.insert(index, inside);
            }
        }
    }
    for x in 0..counts[0] {
        for y in 0..counts[1] {
            for z in 0..counts[2] {
                builder.polygonise([x, y, z], dec)?;
            }
        }
    }
    builder.mesh.max_deviation = &builder.edge * &q(4);
    Ok(Some(builder.mesh))
}
