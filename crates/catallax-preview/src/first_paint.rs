//! Coarse first paint of the field display from batch-sampled values.
//!
//! The field is sampled on a power-of-two grid by an [`Evaluator`] and a
//! surface-nets mesh is extracted on the processor. The mesh is an unproved
//! display picture. It has no writer, no face ownership, no deviation bound and
//! no conversion to a certified mesh.
use crate::batch::{Evaluator, Grid, Hardware, Samples};
use catallax_field::Plan;
use catallax_numbers::{Decider, Refusal};

/// Whole pitches of padding around the field's bounds, so every sign change
/// lies inside the grid and the mesh closes.
pub const MARGIN: u32 = 2;

/// Default largest grid a first paint samples: 2^28 points, about 2 GB of
/// binary64 values. Larger grids refuse with `preview.grid-size`.
pub const MAX_POINTS: usize = 1 << 28;

/// Unproved display triangles sampled from the field.
///
/// ```compile_fail
/// use catallax_preview::first_paint::SampledMesh;
/// fn export(mesh: SampledMesh) { mesh.write_obj(&mut Vec::new()).unwrap(); }
/// ```
///
/// ```compile_fail
/// use catallax_preview::first_paint::SampledMesh;
/// fn certify(mesh: SampledMesh) -> catallax_field::display::Mesh { mesh.into() }
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct SampledMesh {
    /// Binary64 millimetre coordinates, one per surface grid cube.
    pub positions: Vec<[f64; 3]>,
    /// Outward-wound triangles, two per sign-changing grid edge.
    pub triangles: Vec<[u32; 3]>,
    /// Sign changes on the grid's outer layer, which cannot close. Zero when
    /// the margin holds; a positive count means the picture has holes.
    pub open_edges: usize,
}

impl SampledMesh {
    /// Proof status of the picture: always unproved.
    pub const STATUS: &'static str = "unproved";

    /// Signed volume enclosed by the triangles, in cubic millimetres.
    pub fn volume(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|v| self.positions[v as usize]);
                a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                    + a[2] * (b[0] * c[1] - b[1] * c[0])
            })
            .sum::<f64>()
            / 6.0
    }
}

/// A first paint: the grid, its sampled values and the mesh built from them.
#[derive(Clone, Debug)]
pub struct FirstPaint {
    /// The sampled grid.
    pub grid: Grid,
    /// Unproved values at every grid point.
    pub samples: Samples,
    /// The extracted picture.
    pub mesh: SampledMesh,
}

impl FirstPaint {
    /// Where the samples were computed.
    pub fn hardware(&self) -> &Hardware {
        &self.samples.hardware
    }
}

/// The grid at `2^pitch_exponent` mm covering the plan's exact material bounds
/// ([`catallax_field::printing::bounds`]) with [`MARGIN`] pitches on each side.
pub fn grid(plan: &Plan, pitch_exponent: i32, dec: &mut Decider) -> Result<Grid, Refusal> {
    let bounds = catallax_field::printing::bounds(plan, dec)?;
    let lo = std::array::from_fn(|i| bounds[i][0].to_f64_nearest());
    let hi = std::array::from_fn(|i| bounds[i][1].to_f64_nearest());
    // Nearest rounding moves a bound by far less than the margin.
    Grid::covering(lo, hi, pitch_exponent, MARGIN)
}

/// Refuses with `preview.grid-size` when `grid` has more than `max_points`
/// points, naming the counts, before any evaluation.
pub fn check_size(grid: &Grid, max_points: usize) -> Result<(), Refusal> {
    if grid.len() > max_points {
        return Err(Refusal::new(
            "preview.grid-size",
            format!(
                "{} grid points ({} x {} x {}) at pitch {} mm exceed the {max_points}-point limit",
                grid.len(),
                grid.counts[0],
                grid.counts[1],
                grid.counts[2],
                grid.pitch()
            ),
        ));
    }
    Ok(())
}

/// Samples the plan on [`grid`] with `evaluator` and meshes the result.
/// A grid above `max_points` refuses with `preview.grid-size`.
pub fn first_paint(
    plan: &Plan,
    evaluator: &Evaluator<'_, '_>,
    pitch_exponent: i32,
    max_points: usize,
    dec: &mut Decider,
) -> Result<FirstPaint, Refusal> {
    let grid = grid(plan, pitch_exponent, dec)?;
    check_size(&grid, max_points)?;
    let samples = evaluator.grid(&grid)?;
    let mesh = surface_nets(&grid, &samples.values)?;
    Ok(FirstPaint {
        grid,
        samples,
        mesh,
    })
}

/// Surface nets over grid values. A point is inside when its value is
/// negative, matching the ray marcher; zero counts as outside. Each cube with
/// mixed corners gets the mean of its edge crossings, found by linear
/// interpolation; each sign-changing edge gets a quad of its four cubes. The
/// output depends only on the values and the grid, in a fixed order.
pub fn surface_nets(grid: &Grid, values: &[f64]) -> Result<SampledMesh, Refusal> {
    if values.len() != grid.len() {
        return Err(Refusal::new(
            "preview.grid-values",
            format!("{} values for {} grid points", values.len(), grid.len()),
        ));
    }
    let [nx, ny, nz] = grid.counts.map(|n| n as usize);
    if nx < 2 || ny < 2 || nz < 2 {
        return Ok(SampledMesh {
            positions: Vec::new(),
            triangles: Vec::new(),
            open_edges: 0,
        });
    }
    let index = |i: usize, j: usize, k: usize| i + nx * (j + ny * k);
    let (cx, cy) = (nx - 1, ny - 1);
    let cube = |i: usize, j: usize, k: usize| i + cx * (j + cy * k);
    let inside = |v: f64| v < 0.0;
    const NONE: u32 = u32::MAX;
    let mut vertex = vec![NONE; cx * cy * (nz - 1)];
    let mut positions = Vec::new();
    const CORNERS: [[usize; 3]; 8] = [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [1, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [0, 1, 1],
        [1, 1, 1],
    ];
    const EDGES: [[usize; 2]; 12] = [
        [0, 1],
        [2, 3],
        [4, 5],
        [6, 7],
        [0, 2],
        [1, 3],
        [4, 6],
        [5, 7],
        [0, 4],
        [1, 5],
        [2, 6],
        [3, 7],
    ];
    for k in 0..nz - 1 {
        for j in 0..ny - 1 {
            for i in 0..nx - 1 {
                let v = CORNERS.map(|c| values[index(i + c[0], j + c[1], k + c[2])]);
                let count = v.iter().filter(|x| inside(**x)).count();
                if count == 0 || count == 8 {
                    continue;
                }
                let mut sum = [0.0; 3];
                let mut crossings = 0.0;
                for [a, b] in EDGES {
                    if inside(v[a]) == inside(v[b]) {
                        continue;
                    }
                    let t = v[a] / (v[a] - v[b]);
                    for axis in 0..3 {
                        let pa = CORNERS[a][axis] as f64;
                        let pb = CORNERS[b][axis] as f64;
                        sum[axis] += pa + t * (pb - pa);
                    }
                    crossings += 1.0;
                }
                let base = grid.at(i as u32, j as u32, k as u32);
                let pitch = grid.pitch();
                if positions.len() >= NONE as usize {
                    return Err(Refusal::new(
                        "preview.mesh-size",
                        "sampled mesh exceeds u32 vertex indices",
                    ));
                }
                vertex[cube(i, j, k)] = positions.len() as u32;
                positions.push(std::array::from_fn(|axis| {
                    base[axis] + pitch * (sum[axis] / crossings)
                }));
            }
        }
    }
    let mut triangles = Vec::new();
    let mut open_edges = 0;
    let n = [nx, ny, nz];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let p = [i, j, k];
                let here = inside(values[index(i, j, k)]);
                for axis in 0..3 {
                    if p[axis] + 1 >= n[axis] {
                        continue;
                    }
                    let mut q = p;
                    q[axis] += 1;
                    if here == inside(values[index(q[0], q[1], q[2])]) {
                        continue;
                    }
                    let (u, w) = ((axis + 1) % 3, (axis + 2) % 3);
                    if p[u] == 0 || p[w] == 0 || p[u] + 1 >= n[u] || p[w] + 1 >= n[w] {
                        open_edges += 1;
                        continue;
                    }
                    // The four cubes round the edge, counter-clockwise seen from +axis.
                    let quad = [(1, 1), (0, 1), (0, 0), (1, 0)].map(|(du, dw)| {
                        let mut c = p;
                        c[u] -= du;
                        c[w] -= dw;
                        vertex[cube(c[0], c[1], c[2])]
                    });
                    debug_assert!(quad.iter().all(|&v| v != NONE));
                    // Material on the lower end means the outward normal is +axis.
                    let [a, b, c, d] = if here {
                        quad
                    } else {
                        [quad[3], quad[2], quad[1], quad[0]]
                    };
                    triangles.push([a, b, c]);
                    triangles.push([a, c, d]);
                }
            }
        }
    }
    Ok(SampledMesh {
        positions,
        triangles,
        open_edges,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sphere_grid(radius: f64, pitch_exponent: i32) -> (Grid, Vec<f64>) {
        let grid = Grid::covering([-radius; 3], [radius; 3], pitch_exponent, MARGIN).unwrap();
        let values = (0..grid.len())
            .map(|i| {
                let p = grid.point(i);
                (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() - radius
            })
            .collect();
        (grid, values)
    }

    #[test]
    fn sphere_is_closed_outward_and_near_its_volume() {
        let (grid, values) = sphere_grid(5.0, -2);
        let mesh = surface_nets(&grid, &values).unwrap();
        assert_eq!(mesh.open_edges, 0);
        let exact = 4.0 / 3.0 * std::f64::consts::PI * 125.0;
        let volume = mesh.volume();
        assert!((volume - exact).abs() < 0.01 * exact, "{volume} vs {exact}");
        // Every undirected edge is used once in each direction.
        let mut edges = std::collections::HashMap::new();
        for t in &mesh.triangles {
            for e in 0..3 {
                let (a, b) = (t[e], t[(e + 1) % 3]);
                *edges.entry((a.min(b), a.max(b))).or_insert(0i64) += if a < b { 1 } else { -1 };
            }
        }
        assert!(edges.values().all(|&c| c == 0));
        // Vertices lie within one pitch of the sphere.
        for p in &mesh.positions {
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!((r - 5.0).abs() < grid.pitch());
        }
    }

    #[test]
    fn same_values_give_identical_meshes_and_wrong_lengths_refuse() {
        let (grid, values) = sphere_grid(3.0, -1);
        assert_eq!(
            surface_nets(&grid, &values).unwrap(),
            surface_nets(&grid, &values).unwrap()
        );
        assert_eq!(
            surface_nets(&grid, &values[1..]).unwrap_err().code,
            "preview.grid-values"
        );
    }

    #[test]
    fn p4_first_paint_closes_and_refuses_oversized_grids() {
        use crate::batch::Options;
        use catallax_recipe::Recipe;
        let plan = Plan::compile(
            &Recipe::parse(include_str!("../../../parts/P4-shelled-tray.json")).unwrap(),
        )
        .unwrap();
        let options = Options {
            workers: 2,
            ..Options::default()
        };
        let processor = Evaluator::new(&plan, None, options.clone()).unwrap();
        let paint =
            first_paint(&plan, &processor, -1, MAX_POINTS, &mut Decider::default()).unwrap();
        assert!(matches!(paint.hardware(), Hardware::Processor { .. }));
        assert_eq!(paint.mesh.open_edges, 0);
        assert!(!paint.mesh.triangles.is_empty());
        // A coarse picture of the tray's material, not a measurement.
        assert!(paint.mesh.volume() > 0.0);
        assert_eq!(
            first_paint(
                &plan,
                &processor,
                -1,
                paint.grid.len() - 1,
                &mut Decider::default()
            )
            .unwrap_err()
            .code,
            "preview.grid-size"
        );
        if let Ok(gpu) = crate::gpu::Gpu::new(false) {
            let device = Evaluator::new(&plan, Some(&gpu), options).unwrap();
            let gpu_paint =
                first_paint(&plan, &device, -1, MAX_POINTS, &mut Decider::default()).unwrap();
            assert!(matches!(
                gpu_paint.hardware(),
                Hardware::GraphicsCard { .. }
            ));
            assert_eq!(gpu_paint.grid, paint.grid);
            assert_eq!(gpu_paint.mesh.triangles, paint.mesh.triangles);
        }
    }

    #[test]
    fn material_touching_the_grid_edge_is_counted_open() {
        let grid = Grid {
            origin: [0; 3],
            counts: [4; 3],
            pitch_exponent: 0,
        };
        let values: Vec<_> = (0..grid.len())
            .map(|i| if grid.point(i)[0] < 1.0 { -1.0 } else { 1.0 })
            .collect();
        assert!(surface_nets(&grid, &values).unwrap().open_edges > 0);
    }
}
