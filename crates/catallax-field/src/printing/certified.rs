//! A coarse conforming surface with independently refined distance certificates.
//!
//! Surface-to-mesh coverage excludes only interval-certified interior/exterior
//! boxes. Every remaining box is wholly within the bound of one mesh triangle.
//! Mesh-to-surface coverage places every triangle patch between an inside box
//! and an outside box, translated by a vector whose Euclidean length is at most the bound.
use super::*;
mod refined;
mod separation;

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
    clip: Option<&separation::Clip<'_>>,
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
            if let Some(clip) = clip
                && separation::bracket_along_gradient(clip, &patch, classifier, bound, dec)?
            {
                continue;
            }
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
    let mut classifier = Classifier::for_print(plan, dec)?;
    let charts = parametric::Charts::new(&classifier, &body_bounds);
    let clip = separation::Clip::new(&classifier, &body_bounds);
    // A graded gyroid keeps one thirty-second of its period as the coarse pitch at
    // every accuracy, so a coarse request still meets the curvature bound.
    let pitch = plan
        .leaves
        .iter()
        .filter_map(|leaf| match leaf {
            geometry::Leaf::Gyroid(gyroid) if gyroid.graded() => {
                Some(gyroid.period() * &Rational::from_ratio(1, 32).unwrap())
            }
            _ => None,
        })
        .min();
    let mut grid = refined::Grid::new(
        &body_bounds,
        options,
        charts.is_some(),
        clip.is_some(),
        pitch.as_ref(),
    )?;
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
        let mut failed = cover_mesh(
            &mut mesh,
            &mut classifier,
            clip.as_ref(),
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
            let coverage = if let Some(charts) = &charts {
                charts.cover(
                    &mut mesh,
                    &mut classifier,
                    interval.square().lo(),
                    options,
                    dec,
                )
            } else if let Some(clip) = &clip {
                match separation::Separation::new(clip, &grid, interval) {
                    Some(test) => {
                        test.cover(&mut mesh, &mut classifier, &body_bounds, options, dec)
                    }
                    None => Ok(Some(0)),
                }
            } else {
                cover_surface(
                    &mut mesh,
                    &mut classifier,
                    &body_bounds,
                    interval.square().lo(),
                    options,
                    dec,
                )
            }
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

/// An axis plane of a convex support: its axis, coefficient and enclosed position.
fn axis_plane(support: &Terms) -> Option<(usize, Interval, Interval)> {
    let mut constant = point(0.0);
    let mut linear = None;
    for (powers, coefficient) in support {
        if *powers == [0, 0, 0] {
            constant = constant + *coefficient;
        } else {
            let axis = (0..3).find(|&i| *powers == std::array::from_fn(|j| u32::from(i == j)))?;
            if linear.is_some() || coefficient.lo() <= 0.0 && coefficient.hi() >= 0.0 {
                return None;
            }
            linear = Some((axis, *coefficient));
        }
    }
    let (normal, coefficient) = linear?;
    let position = (point(0.0) - constant) / coefficient;
    (position.lo().is_finite() && position.hi().is_finite()).then_some((
        normal,
        coefficient,
        position,
    ))
}

/// Complete boundary charts for an intersection of axis planes and one constant
/// origin-aligned cubic lattice. Unsupported definitions use volume coverage.
mod parametric {
    use super::*;

    #[derive(Clone, Copy)]
    enum Patch {
        Plane {
            normal: usize,
            position: Interval,
            uv: [Interval; 2],
        },
        Cylinder {
            axis: usize,
            centre: [Interval; 2],
            radius: Interval,
            quadrant: usize,
            uv: [Interval; 2],
        },
    }
    impl Patch {
        fn uv(self) -> [Interval; 2] {
            match self {
                Self::Plane { uv, .. } | Self::Cylinder { uv, .. } => uv,
            }
        }
        fn with_uv(self, uv: [Interval; 2]) -> Self {
            match self {
                Self::Plane {
                    normal, position, ..
                } => Self::Plane {
                    normal,
                    position,
                    uv,
                },
                Self::Cylinder {
                    axis,
                    centre,
                    radius,
                    quadrant,
                    ..
                } => Self::Cylinder {
                    axis,
                    centre,
                    radius,
                    quadrant,
                    uv,
                },
            }
        }
        fn bounds(self) -> [Interval; 3] {
            match self {
                Self::Plane {
                    normal,
                    position,
                    uv,
                } => {
                    let mut bounds = [position; 3];
                    bounds[(normal + 1) % 3] = uv[0];
                    bounds[(normal + 2) % 3] = uv[1];
                    bounds
                }
                Self::Cylinder {
                    axis,
                    centre,
                    radius,
                    quadrant,
                    uv,
                } => {
                    // The rational quarter-circle chart has decreasing x and
                    // increasing y on [0,1]. Outward endpoint evaluations enclose
                    // the whole arc, not just the evaluated endpoints.
                    let xy = |t: f64| {
                        let t = point(t);
                        let denominator = point(1.0) + t.square();
                        [
                            (point(1.0) - t.square()) / denominator,
                            point(2.0) * t / denominator,
                        ]
                    };
                    let a = xy(uv[1].lo());
                    let b = xy(uv[1].hi());
                    let x = Interval::new(b[0].lo(), a[0].hi()).unwrap() * radius;
                    let y = Interval::new(a[1].lo(), b[1].hi()).unwrap() * radius;
                    let [x, y] = match quadrant {
                        0 => [x, y],
                        1 => [point(0.0) - y, x],
                        2 => [point(0.0) - x, point(0.0) - y],
                        _ => [y, point(0.0) - x],
                    };
                    let mut bounds = [uv[0]; 3];
                    bounds[(axis + 1) % 3] = centre[0] + x;
                    bounds[(axis + 2) % 3] = centre[1] + y;
                    bounds
                }
            }
        }
    }

    pub(super) struct Charts {
        planes: Vec<Patch>,
        clip: [Interval; 3],
        lines: [[i64; 2]; 3],
        size: Rational,
        radius: Interval,
    }
    impl Charts {
        pub(super) fn new(
            classifier: &Classifier<'_>,
            bounds: &[[Rational; 2]; 3],
        ) -> Option<Self> {
            if classifier.plan.instructions.iter().any(|i| {
                matches!(i, Instruction::Boolean { kind, .. } if *kind != BooleanKind::Intersect)
            }) {
                return None;
            }
            let clip = std::array::from_fn(|i| {
                Interval::new(
                    Interval::from_rational(&bounds[i][0]).lo(),
                    Interval::from_rational(&bounds[i][1]).hi(),
                )
                .unwrap()
            });
            let mut planes = Vec::new();
            let mut parameters = None;
            for (id, definition) in classifier.definitions.iter().enumerate() {
                match definition {
                    None => (),
                    Some(Definition::Lattice) => {
                        if parameters.is_some() {
                            return None;
                        }
                        let geometry::Leaf::Lattice(lattice) = &classifier.plan.leaves[id] else {
                            unreachable!("compiled lattice");
                        };
                        let (translation, size, radius) = lattice.axis_cubic_grid()?;
                        // Signed axis permutations describe the same cubic beam
                        // set. Nonzero translations retain the general route.
                        if translation.iter().any(|v| !v.is_zero()) {
                            return None;
                        }
                        parameters = Some((size, radius));
                    }
                    Some(Definition::Cell(supports)) => {
                        for support in supports {
                            let (normal, _, position) = axis_plane(support)?;
                            planes.push(Patch::Plane {
                                normal,
                                position,
                                uv: [clip[(normal + 1) % 3], clip[(normal + 2) % 3]],
                            });
                        }
                    }
                    _ => return None,
                }
            }
            let (size, radius) = parameters?;
            let mut lines = [[0; 2]; 3];
            for i in 0..3 {
                // Rational floor/ceiling includes every beam that can meet the
                // clip box. Extra outside lines are excluded by the classifier.
                lines[i] = [
                    (&(&bounds[i][0] - &radius) / &size).floor().to_i64()?,
                    (&(&bounds[i][1] + &radius) / &size).ceil().to_i64()?,
                ];
            }
            Some(Self {
                planes,
                clip,
                lines,
                size,
                radius: Interval::from_rational(&radius),
            })
        }

        pub(super) fn cover(
            &self,
            mesh: &mut Mesh,
            classifier: &mut Classifier<'_>,
            limit_squared: f64,
            options: &Options,
            dec: &mut Decider,
        ) -> Result<Option<usize>, Refusal> {
            let mut tree = Tree::new(mesh);
            let mut pending = Vec::new();
            // Generate roots lazily. Even an enormous beam range hits the
            // existing cell/work guards without allocating all its charts.
            let cylinders = (0..3).flat_map(|axis| {
                let j = (axis + 1) % 3;
                let k = (axis + 2) % 3;
                (self.lines[j][0]..=self.lines[j][1]).flat_map(move |a| {
                    (self.lines[k][0]..=self.lines[k][1]).flat_map(move |b| {
                        let centre = [
                            Interval::from_rational(&(&self.size * &q(a))),
                            Interval::from_rational(&(&self.size * &q(b))),
                        ];
                        (0..4).map(move |quadrant| Patch::Cylinder {
                            axis,
                            centre,
                            radius: self.radius,
                            quadrant,
                            uv: [self.clip[axis], Interval::new(0.0, 1.0).unwrap()],
                        })
                    })
                })
            });
            for root in self.planes.iter().copied().chain(cylinders) {
                pending.push((root, 0u32));
                while let Some((patch, depth)) = pending.pop() {
                    spend(&mut mesh.cells, options, dec)?;
                    let bounds = patch.bounds();
                    if classifier.classify(&bounds, dec)? != BoxState::Unknown {
                        continue;
                    }
                    let centre = bounds.map(|v| v.lo() * 0.5 + v.hi() * 0.5);
                    let id = tree.nearest(mesh, centre);
                    let t = triangle(mesh, id);
                    // Distance to one convex triangle is convex. Certified
                    // upper bounds at every corner therefore cover this whole
                    // box and every boundary point in the enclosed chart.
                    if (0..8).all(|bits| {
                        distance_squared(
                            std::array::from_fn(|i| {
                                if bits & (1 << i) == 0 {
                                    bounds[i].lo()
                                } else {
                                    bounds[i].hi()
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
                    let uv = patch.uv();
                    let axis = match patch {
                        Patch::Plane { .. } => {
                            usize::from(uv[1].hi() - uv[1].lo() > uv[0].hi() - uv[0].lo())
                        }
                        Patch::Cylinder { axis, .. } => {
                            let width = |i: usize| bounds[i].hi() - bounds[i].lo();
                            usize::from(
                                width((axis + 1) % 3).max(width((axis + 2) % 3))
                                    > uv[0].hi() - uv[0].lo(),
                            )
                        }
                    };
                    let mid = uv[axis].lo() * 0.5 + uv[axis].hi() * 0.5;
                    let mut a = uv;
                    let mut b = uv;
                    a[axis] = Interval::new(uv[axis].lo(), mid).unwrap();
                    b[axis] = Interval::new(mid, uv[axis].hi()).unwrap();
                    pending.push((patch.with_uv(a), depth + 1));
                    pending.push((patch.with_uv(b), depth + 1));
                }
            }
            Ok(None)
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn cylinder_chart_boxes_enclose_exact_interior_points_and_quadrant_seams() {
            let radius = Rational::from_ratio(3, 2).unwrap();
            let centre = [
                Rational::from_ratio(-5, 2).unwrap(),
                Rational::from_ratio(5, 3).unwrap(),
            ];
            for axis in 0..3 {
                for quadrant in 0..4 {
                    for part in 0..4 {
                        let bounds = Patch::Cylinder {
                            axis,
                            centre: centre.each_ref().map(Interval::from_rational),
                            radius: Interval::from_rational(&radius),
                            quadrant,
                            uv: [
                                Interval::new(-1.0, 2.0).unwrap(),
                                Interval::new(part as f64 / 4.0, (part + 1) as f64 / 4.0).unwrap(),
                            ],
                        }
                        .bounds();
                        for sample in 0..=16 {
                            let t = Rational::from_ratio(part * 16 + sample, 64).unwrap();
                            let t2 = &t * &t;
                            let denominator = &q(1) + &t2;
                            let x = &radius * &(&(&q(1) - &t2) / &denominator);
                            let y = &radius * &(&t.mul_pow2(1) / &denominator);
                            let [x, y] = match quadrant {
                                0 => [x, y],
                                1 => [-y, x],
                                2 => [-x, -y],
                                _ => [y, -x],
                            };
                            for height in [q(-1), Rational::from_ratio(1, 3).unwrap(), q(2)] {
                                let mut p: [Rational; 3] = std::array::from_fn(|_| height.clone());
                                p[(axis + 1) % 3] = &centre[0] + &x;
                                p[(axis + 2) % 3] = &centre[1] + &y;
                                for i in 0..3 {
                                    assert!(
                                        Rational::from_f64_exact(bounds[i].lo()).unwrap() <= p[i]
                                    );
                                    assert!(
                                        Rational::from_f64_exact(bounds[i].hi()).unwrap() >= p[i]
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
