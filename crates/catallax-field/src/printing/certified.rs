//! A coarse conforming surface with independently refined distance certificates.
//!
//! Surface-to-mesh coverage excludes only interval-certified interior/exterior
//! boxes. Every remaining box is wholly within the bound of one mesh triangle.
//! Mesh-to-surface coverage places every triangle patch between an inside box
//! and an outside box, translated by a vector whose Euclidean length is at most the bound.
use super::*;
mod refined;

type Point = [f64; 3];
type Triangle = [Point; 3];
fn point(v: f64) -> Interval {
    Interval::point(v).unwrap()
}
fn sub(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}
fn dot(a: Point, b: Point) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn triangle(mesh: &Mesh, id: usize) -> Triangle {
    mesh.triangles[id].map(|v| mesh.positions[v].map(f64::from))
}

// These weights select a point, not an authoritative nearest-point answer.
// Quantisation makes nonnegativity and u+v<=1 exact in the interval calculation.
fn weights(p: Point, t: Triangle) -> [f64; 2] {
    let a = sub(t[1], t[0]);
    let b = sub(t[2], t[0]);
    let d = sub(p, t[0]);
    let aa = dot(a, a);
    let ab = dot(a, b);
    let bb = dot(b, b);
    let ad = dot(a, d);
    let bd = dot(b, d);
    let determinant = aa * bb - ab * ab;
    let (u, v) = (
        (bb * ad - ab * bd) / determinant,
        (aa * bd - ab * ad) / determinant,
    );
    let mut candidates = [[0.0, 0.0]; 4];
    candidates[0] = if u >= 0.0 && v >= 0.0 && u + v <= 1.0 {
        [u, v]
    } else {
        [0.0, 0.0]
    };
    candidates[1] = [(ad / aa).clamp(0.0, 1.0), 0.0];
    candidates[2] = [0.0, (bd / bb).clamp(0.0, 1.0)];
    let along = ((bd - ad - ab + aa) / (aa - 2.0 * ab + bb)).clamp(0.0, 1.0);
    candidates[3] = [1.0 - along, along];
    let [u, v] = candidates
        .into_iter()
        .min_by(|x, y| {
            let distance = |w: &[f64; 2]| {
                (0..3)
                    .map(|i| (t[0][i] + w[0] * a[i] + w[1] * b[i] - p[i]).powi(2))
                    .sum::<f64>()
            };
            distance(x).total_cmp(&distance(y))
        })
        .unwrap();
    let scale = 16777216.0;
    let u = (u * scale).floor().clamp(0.0, scale);
    let v = (v * scale).floor().clamp(0.0, scale - u);
    [u / scale, v / scale]
}
fn distance_squared(p: Point, t: Triangle) -> f64 {
    let [u, v] = weights(p, t);
    let mut sum = point(0.0);
    for i in 0..3 {
        let nearest = point(t[0][i])
            + point(u) * (point(t[1][i]) - point(t[0][i]))
            + point(v) * (point(t[2][i]) - point(t[0][i]));
        sum = sum + (nearest - point(p[i])).square();
    }
    sum.hi()
}
struct Node {
    bounds: [[f64; 2]; 3],
    range: std::ops::Range<usize>,
    children: Option<[usize; 2]>,
}
struct Tree {
    nodes: Vec<Node>,
    ids: Vec<usize>,
    pending: Vec<usize>,
}
impl Tree {
    fn new(mesh: &Mesh) -> Self {
        let mut tree = Self {
            nodes: Vec::new(),
            ids: (0..mesh.triangles.len()).collect(),
            pending: Vec::new(),
        };
        tree.build(mesh, 0..tree.ids.len());
        tree
    }
    fn build(&mut self, mesh: &Mesh, range: std::ops::Range<usize>) -> usize {
        let bounds = std::array::from_fn(|axis| {
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for &id in &self.ids[range.clone()] {
                for p in triangle(mesh, id) {
                    lo = lo.min(p[axis]);
                    hi = hi.max(p[axis]);
                }
            }
            [lo, hi]
        });
        let id = self.nodes.len();
        self.nodes.push(Node {
            bounds,
            range: range.clone(),
            children: None,
        });
        if range.len() > 8 {
            let axis = (0..3)
                .max_by(|&a, &b| {
                    (bounds[a][1] - bounds[a][0]).total_cmp(&(bounds[b][1] - bounds[b][0]))
                })
                .unwrap();
            let middle = range.len() / 2;
            self.ids[range.clone()].select_nth_unstable_by(middle, |&a, &b| {
                let centre = |id| triangle(mesh, id).iter().map(|p| p[axis]).sum::<f64>();
                centre(a).total_cmp(&centre(b)).then_with(|| a.cmp(&b))
            });
            let middle = range.start + middle;
            let a = self.build(mesh, range.start..middle);
            let b = self.build(mesh, middle..range.end);
            self.nodes[id].children = Some([a, b]);
        }
        id
    }
    fn nearest(&mut self, mesh: &Mesh, p: Point) -> usize {
        self.pending.clear();
        self.pending.push(0);
        let mut best = f64::INFINITY;
        let mut result = 0;
        while let Some(id) = self.pending.pop() {
            let node = &self.nodes[id];
            let lower = (0..3)
                .map(|i| (p[i] - p[i].clamp(node.bounds[i][0], node.bounds[i][1])).powi(2))
                .sum::<f64>();
            // Pruning only chooses the witness triangle. A wrong choice can make
            // certification slower or refuse, but cannot certify a false bound.
            if lower > best {
                continue;
            }
            if let Some([a, b]) = node.children {
                let distance = |id: usize| {
                    (0..3)
                        .map(|i| {
                            let bounds = self.nodes[id].bounds[i];
                            (p[i] - p[i].clamp(bounds[0], bounds[1])).powi(2)
                        })
                        .sum::<f64>()
                };
                let children = if distance(a) < distance(b) {
                    [b, a]
                } else {
                    [a, b]
                };
                self.pending.extend(children);
            } else {
                for &id in &self.ids[node.range.clone()] {
                    let t = triangle(mesh, id);
                    let [u, v] = weights(p, t);
                    let distance = (0..3)
                        .map(|i| {
                            (t[0][i] + u * (t[1][i] - t[0][i]) + v * (t[2][i] - t[0][i]) - p[i])
                                .powi(2)
                        })
                        .sum::<f64>();
                    if distance < best {
                        best = distance;
                        result = id;
                    }
                }
            }
        }
        result
    }
}
fn spend(cells: &mut usize, options: &Options, dec: &mut Decider) -> Result<(), Refusal> {
    charge(dec)?;
    *cells += 1;
    if *cells > options.max_cells {
        return Err(refusal(
            "field.print-budget",
            "adaptive print certificate cell limit reached",
        ));
    }
    Ok(())
}
type FacetKey = [[u32; 3]; 3];
fn bracket(
    classifier: &mut Classifier<'_>,
    box_: [Interval; 3],
    bound: Interval,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for direction in [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 1, 0],
        [1, -1, 0],
        [1, 0, 1],
        [1, 0, -1],
        [0, 1, 1],
        [0, 1, -1],
        [1, 1, 1],
        [1, 1, -1],
        [1, -1, 1],
        [1, -1, -1],
    ] {
        let nonzero = direction.iter().filter(|&&v| v != 0).count();
        // Exact dyadic scales with k*scale^2 <= 1 keep every translation
        // inside the Euclidean accuracy ball, including diagonal directions.
        let scale = match nonzero {
            1 => 1.0,
            2 => 181.0 / 256.0,
            3 => 591.0 / 1024.0,
            _ => unreachable!(),
        };
        let step = bound * point(scale);
        let shift: [Interval; 3] = direction.map(|v| step * point(v as f64));
        let lo = std::array::from_fn(|i| box_[i] - shift[i]);
        let hi = std::array::from_fn(|i| box_[i] + shift[i]);
        let a = classifier.classify(&lo, dec)?;
        let b = classifier.classify(&hi, dec)?;
        if matches!(
            (a, b),
            (BoxState::Inside, BoxState::Outside) | (BoxState::Outside, BoxState::Inside)
        ) {
            return Ok(true);
        }
    }
    Ok(false)
}
fn cover_mesh(
    mesh: &mut Mesh,
    classifier: &mut Classifier<'_>,
    bound: Interval,
    cache: &mut BTreeSet<FacetKey>,
    options: &Options,
    dec: &mut Decider,
) -> Result<BTreeSet<usize>, Refusal> {
    let mut patches = Vec::new();
    let mut failed = BTreeSet::new();
    for id in 0..mesh.triangles.len() {
        let key = mesh.triangles[id].map(|v| mesh.positions[v].map(f32::to_bits));
        if cache.contains(&key) {
            continue;
        }
        patches.push((triangle(mesh, id), 0u32));
        while let Some((patch, depth)) = patches.pop() {
            spend(&mut mesh.cells, options, dec)?;
            let box_: [Interval; 3] = std::array::from_fn(|i| {
                Interval::new(
                    patch.iter().map(|p| p[i]).fold(f64::INFINITY, f64::min),
                    patch.iter().map(|p| p[i]).fold(f64::NEG_INFINITY, f64::max),
                )
                .unwrap()
            });
            if bracket(classifier, box_, bound, dec)? {
                continue;
            }
            let centre: [Interval; 3] = std::array::from_fn(|i| {
                (point(patch[0][i]) + point(patch[1][i]) + point(patch[2][i])) / point(3.0)
            });
            if depth == 12 || (depth >= 2 && !bracket(classifier, centre, bound, dec)?) {
                failed.insert(id);
                patches.clear();
                break;
            }
            let mid: [Point; 3] = std::array::from_fn(|i| {
                std::array::from_fn(|k| (patch[i][k] + patch[(i + 1) % 3][k]) * 0.5)
            });
            patches.extend([
                ([patch[0], mid[0], mid[2]], depth + 1),
                ([mid[0], patch[1], mid[1]], depth + 1),
                ([mid[2], mid[1], patch[2]], depth + 1),
                ([mid[0], mid[1], mid[2]], depth + 1),
            ]);
        }
        if !failed.contains(&id) {
            cache.insert(key);
        }
    }
    Ok(failed)
}
fn cover_surface(
    mesh: &mut Mesh,
    classifier: &mut Classifier<'_>,
    bounds: &[[Rational; 2]; 3],
    limit_squared: f64,
    options: &Options,
    dec: &mut Decider,
) -> Result<Option<usize>, Refusal> {
    let mut tree = Tree::new(mesh);
    let root: [Interval; 3] = std::array::from_fn(|i| {
        Interval::new(
            Interval::from_rational(&bounds[i][0]).lo(),
            Interval::from_rational(&bounds[i][1]).hi(),
        )
        .unwrap()
    });
    let mut pending = vec![(root, 0u32)];
    while let Some((box_, depth)) = pending.pop() {
        spend(&mut mesh.cells, options, dec)?;
        if classifier.classify(&box_, dec)? != BoxState::Unknown {
            continue;
        }
        let centre = box_.map(|b| b.lo() * 0.5 + b.hi() * 0.5);
        let id = tree.nearest(mesh, centre);
        let t = triangle(mesh, id);
        if (0..8).all(|bits| {
            distance_squared(
                std::array::from_fn(|i| {
                    if bits & (1 << i) == 0 {
                        box_[i].lo()
                    } else {
                        box_[i].hi()
                    }
                }),
                t,
            ) <= limit_squared
        }) {
            continue;
        }
        if depth == 60 {
            return Ok(Some(id));
        }
        let axis = (0..3)
            .max_by(|&a, &b| {
                (box_[a].hi() - box_[a].lo()).total_cmp(&(box_[b].hi() - box_[b].lo()))
            })
            .unwrap();
        let mid = centre[axis];
        let mut a = box_;
        let mut b = box_;
        a[axis] = Interval::new(box_[axis].lo(), mid).unwrap();
        b[axis] = Interval::new(mid, box_[axis].hi()).unwrap();
        pending.push((a, depth + 1));
        pending.push((b, depth + 1));
    }
    Ok(None)
}

pub(crate) fn mesh(
    plan: &Plan,
    options: &Options,
    floor: &Rational,
    dec: &mut Decider,
) -> Result<Mesh, Refusal> {
    if plan.is_proven_empty() {
        return Ok(Mesh {
            positions: vec![],
            triangles: vec![],
            max_deviation: q(0),
            cells: 0,
        });
    }
    let body_bounds = bounds(plan, dec)?;
    let mut grid = refined::Grid::new(&body_bounds, options)?;
    let rounded = options.tolerance.to_f64_nearest();
    let bound = Rational::from_f64_exact(rounded.next_down()).ok_or_else(|| {
        refusal(
            "field.print-grid",
            "certificate bound is outside binary64 range",
        )
    })?;
    if bound <= q(0) {
        return Err(refusal("field.print-grid", "certificate bound underflows"));
    }
    let interval = Interval::from_rational(&bound);
    let mut visited = 0;
    let mut cache = BTreeSet::new();
    for refinement in 0..=3 {
        let (mut mesh, owners) = grid.mesh(plan, options, visited, dec)?;
        if mesh.triangles.is_empty() {
            return Err(refusal(
                "field.print-unresolved",
                "the candidate grid did not resolve a boundary; minimum detail was retained",
            ));
        }
        let mut classifier = Classifier::for_print(plan, dec)?;
        let mut failed = cover_mesh(
            &mut mesh,
            &mut classifier,
            interval,
            &mut cache,
            options,
            dec,
        )
        .map_err(|e| {
            refusal(
                e.code,
                format!(
                    "{}; cells={} vertices={}",
                    e.detail,
                    mesh.cells,
                    mesh.positions.len()
                ),
            )
        })?;
        if failed.is_empty() {
            let coverage = cover_surface(
                &mut mesh,
                &mut classifier,
                &body_bounds,
                interval.square().lo(),
                options,
                dec,
            )
            .map_err(|e| {
                refusal(
                    e.code,
                    format!(
                        "{}; cells={} vertices={}",
                        e.detail,
                        mesh.cells,
                        mesh.positions.len()
                    ),
                )
            })?;
            if let Some(id) = coverage {
                failed.insert(id);
            } else {
                mesh.max_deviation = bound;
                return Ok(mesh);
            }
        }
        let seeds = failed
            .iter()
            .map(|&id| owners[owners.partition_point(|(end, _)| *end <= id)].1)
            .collect();
        visited = mesh.cells;
        if refinement == 3 {
            return Err(refusal(
                "field.print-unresolved",
                format!(
                    "local refinement depth exhausted; cells={} vertices={} failing_facets={}",
                    mesh.cells,
                    mesh.positions.len(),
                    failed.len()
                ),
            ));
        }
        grid.refine(&seeds, &body_bounds, floor).map_err(|e| {
            refusal(
                e.code,
                format!(
                    "{}; cells={} vertices={} failing_facets={}",
                    e.detail,
                    mesh.cells,
                    mesh.positions.len(),
                    failed.len()
                ),
            )
        })?;
    }
    unreachable!()
}
