//! Coordinate-plane propagation keeps locally seeded refinements conforming.
use super::*;
pub(super) type Cell = [[u32; 2]; 3];
pub(super) struct Grid {
    origin: [Rational; 3],
    edge: Rational,
    axes: [Vec<u32>; 3],
}
impl Grid {
    pub(super) fn new(bounds: &[[Rational; 2]; 3], options: &Options) -> Result<Self, Refusal> {
        let accuracy = &options.tolerance;
        let mut coarse = q(1);
        while coarse > accuracy * &q(4) {
            coarse = coarse.mul_pow2(-1);
        }
        while coarse.mul_pow2(1) <= accuracy * &q(4) {
            coarse = coarse.mul_pow2(1);
        }
        let origin = std::array::from_fn(|i| {
            &(&coarse * &Rational::from_integer((&bounds[i][0] / &coarse).floor())) - &coarse
        });
        let mut axes: [Vec<u32>; 3] = std::array::from_fn(|_| Vec::new());
        for i in 0..3 {
            let n = (&(&bounds[i][1] - &origin[i]) / &coarse)
                .to_f64_nearest()
                .ceil()
                + 1.0;
            if !n.is_finite() || n > (1u32 << 25) as f64 {
                return Err(refusal(
                    "field.print-grid",
                    "adaptive grid index range exhausted",
                ));
            }
            if n > options.max_cells as f64 {
                return Err(refusal(
                    "field.print-budget",
                    "initial adaptive grid exceeds the cell limit",
                ));
            }
            axes[i] = (0..=n as u32).map(|j| j * 16).collect();
            // Accuracy sets this initial clip grid independently of the detail
            // floor, just as it sets the initial grid of the uniform mesher.
            // Subsequent feature-chasing refinements still obey the floor.
            for boundary in &bounds[i] {
                let bin = (&(boundary - &origin[i]) / &coarse)
                    .to_f64_nearest()
                    .floor() as u32;
                for j in bin.saturating_sub(1)..=bin.min(n as u32 - 1) {
                    axes[i].extend([j * 16 + 4, j * 16 + 8, j * 16 + 12]);
                }
            }
            axes[i].sort_unstable();
            axes[i].dedup();
        }
        Ok(Self {
            origin,
            edge: coarse.mul_pow2(-4),
            axes,
        })
    }
    pub(super) fn refine(
        &mut self,
        cells: &BTreeSet<Cell>,
        bounds: &[[Rational; 2]; 3],
        floor: &Rational,
    ) -> Result<(), Refusal> {
        let mut additions: [BTreeSet<u32>; 3] = std::array::from_fn(|_| BTreeSet::new());
        for cell in cells {
            // Prefer a clip-boundary normal. Splits propagate along that plane
            // so adjacent rectangles use identical edge vertices and diagonals.
            let axis = (0..3)
                .min_by(|&a, &b| {
                    let gap = |i: usize| {
                        let lo = &self.origin[i] + &(&self.edge * &q(cell[i][0] as i64));
                        let hi = &self.origin[i] + &(&self.edge * &q(cell[i][1] as i64));
                        let left = (&lo - &bounds[i][0]).abs();
                        let right = (&hi - &bounds[i][1]).abs();
                        left.min(right).to_f64_nearest()
                    };
                    gap(a).total_cmp(&gap(b))
                })
                .unwrap();
            let [lo, hi] = cell[axis];
            if hi - lo < 2 || &self.edge * &q(((hi - lo) / 2) as i64) < *floor {
                return Err(refusal(
                    "field.print-unresolved",
                    "local boundary refinement needs a grid finer than the minimum detail allows",
                ));
            }
            additions[axis].insert(lo + (hi - lo) / 2);
        }
        for (axis, extra) in self.axes.iter_mut().zip(additions) {
            axis.extend(extra);
            axis.sort_unstable();
            axis.dedup();
        }
        Ok(())
    }
    pub(super) fn mesh(
        &self,
        plan: &Plan,
        options: &Options,
        visited: usize,
        dec: &mut Decider,
    ) -> Result<(Mesh, Vec<(usize, Cell)>), Refusal> {
        let mut builder = Builder {
            classifier: Classifier::for_print(plan, dec)?,
            workspace: Workspace::default(),
            origin: self.origin.clone(),
            edge: self.edge.clone(),
            corner_signs: BTreeMap::new(),
            edge_vertices: BTreeMap::new(),
            grid_vertices: Vec::new(),
            mesh: Mesh {
                positions: vec![],
                triangles: vec![],
                max_deviation: q(0),
                cells: visited,
            },
            options,
            root_steps: 8,
        };
        let root: [[usize; 2]; 3] = std::array::from_fn(|i| [0, self.axes[i].len() - 1]);
        let mut pending = vec![root];
        let mut owners = Vec::new();
        while let Some(range) = pending.pop() {
            spend(&mut builder.mesh.cells, options, dec)?;
            let cell: Cell = std::array::from_fn(|i| range[i].map(|j| self.axes[i][j]));
            let lo = builder.point(cell.map(|c| c[0]));
            let hi = builder.point(cell.map(|c| c[1]));
            let box_ = std::array::from_fn(|i| {
                Interval::new(
                    Interval::from_rational(&lo[i]).lo(),
                    Interval::from_rational(&hi[i]).hi(),
                )
                .unwrap()
            });
            if builder.classifier.classify(&box_, dec)? != BoxState::Unknown {
                continue;
            }
            if range.iter().all(|r| r[1] - r[0] == 1) {
                let corners =
                    std::array::from_fn(|bits| std::array::from_fn(|i| cell[i][(bits >> i) & 1]));
                let before = builder.mesh.triangles.len();
                builder.polygonise_corners(corners, dec)?;
                if builder.mesh.triangles.len() > before {
                    owners.push((builder.mesh.triangles.len(), cell));
                }
            } else {
                let axis = (0..3)
                    .filter(|&i| range[i][1] - range[i][0] > 1)
                    .max_by_key(|&i| cell[i][1] - cell[i][0])
                    .unwrap();
                let mid = (range[axis][0] + range[axis][1]) / 2;
                let mut a = range;
                let mut b = range;
                a[axis][1] = mid;
                b[axis][0] = mid;
                pending.push(a);
                pending.push(b);
            }
        }
        check_topology(&builder.mesh)?;
        Ok((builder.mesh, owners))
    }
}
