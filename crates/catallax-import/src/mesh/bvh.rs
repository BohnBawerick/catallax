//! A bounding-volume tree over triangle boxes in binary64.
//!
//! Boxes are rounded outward from the nearest binary64 copies of exact
//! coordinates, so every box contains its exact triangle. The tree only
//! chooses candidates; every decision on a candidate is exact.

const LEAF: usize = 4;

#[derive(Clone, Debug)]
pub(crate) struct Node {
    pub lo: [f64; 3],
    pub hi: [f64; 3],
    /// Leaf: first index into `order`. Inner: index of the left child; the
    /// right child follows the left subtree.
    first: u32,
    /// Leaf: number of triangles. Inner: zero.
    count: u32,
    right: u32,
}

#[derive(Clone, Debug, Default)]
pub(crate) struct Bvh {
    nodes: Vec<Node>,
    order: Vec<u32>,
    boxes: Vec<([f64; 3], [f64; 3])>,
}

/// Whether two closed boxes meet.
pub(crate) fn overlap(a: &([f64; 3], [f64; 3]), b: &([f64; 3], [f64; 3])) -> bool {
    (0..3).all(|i| a.0[i] <= b.1[i] && b.0[i] <= a.1[i])
}

impl Bvh {
    /// Builds the tree over outward triangle boxes. Deterministic for equal input.
    pub(crate) fn new(boxes: Vec<([f64; 3], [f64; 3])>) -> Self {
        let mut tree = Self {
            nodes: Vec::with_capacity(boxes.len() / LEAF * 2 + 1),
            order: (0..boxes.len() as u32).collect(),
            boxes,
        };
        if !tree.boxes.is_empty() {
            tree.build(0, tree.boxes.len());
        }
        tree
    }

    fn build(&mut self, start: usize, end: usize) -> usize {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &t in &self.order[start..end] {
            let b = &self.boxes[t as usize];
            for i in 0..3 {
                lo[i] = lo[i].min(b.0[i]);
                hi[i] = hi[i].max(b.1[i]);
            }
        }
        let index = self.nodes.len();
        self.nodes.push(Node {
            lo,
            hi,
            first: start as u32,
            count: (end - start) as u32,
            right: 0,
        });
        if end - start <= LEAF {
            return index;
        }
        let axis = (0..3)
            .max_by(|&a, &b| (hi[a] - lo[a]).total_cmp(&(hi[b] - lo[b])))
            .unwrap_or(0);
        let middle = start + (end - start) / 2;
        let boxes = &self.boxes;
        let centre = |t: &u32| {
            let b = &boxes[*t as usize];
            b.0[axis] + b.1[axis]
        };
        self.order[start..end].select_nth_unstable_by(middle - start, |a, b| {
            centre(a).total_cmp(&centre(b)).then(a.cmp(b))
        });
        let left = self.build(start, middle);
        let right = self.build(middle, end);
        let node = &mut self.nodes[index];
        node.first = left as u32;
        node.count = 0;
        node.right = right as u32;
        index
    }

    fn children(&self, node: usize) -> Option<(usize, usize)> {
        let n = &self.nodes[node];
        (n.count == 0).then_some((n.first as usize, n.right as usize))
    }

    fn leaf(&self, node: usize) -> &[u32] {
        let n = &self.nodes[node];
        &self.order[n.first as usize..(n.first + n.count) as usize]
    }

    fn node_box(&self, node: usize) -> ([f64; 3], [f64; 3]) {
        (self.nodes[node].lo, self.nodes[node].hi)
    }

    /// Visits every unordered pair of distinct triangles whose boxes meet, in a
    /// fixed order. Stops when `visit` returns `Some`.
    pub(crate) fn pairs<T>(&self, mut visit: impl FnMut(u32, u32) -> Option<T>) -> Option<T> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut stack = vec![(0usize, 0usize)];
        while let Some((a, b)) = stack.pop() {
            if a != b && !overlap(&self.node_box(a), &self.node_box(b)) {
                continue;
            }
            match (self.children(a), self.children(b)) {
                (None, None) => {
                    let (la, lb) = (self.leaf(a), self.leaf(b));
                    for (i, &s) in la.iter().enumerate() {
                        let rest = if a == b { &lb[i + 1..] } else { lb };
                        for &t in rest {
                            if overlap(&self.boxes[s as usize], &self.boxes[t as usize])
                                && let Some(found) = visit(s.min(t), s.max(t))
                            {
                                return Some(found);
                            }
                        }
                    }
                }
                (Some((l, r)), _) if a == b => {
                    stack.push((r, r));
                    stack.push((l, r));
                    stack.push((l, l));
                }
                (Some((l, r)), None) => {
                    stack.push((r, b));
                    stack.push((l, b));
                }
                (None, Some((l, r))) => {
                    stack.push((a, r));
                    stack.push((a, l));
                }
                (Some((al, ar)), Some((bl, br))) => {
                    stack.push((ar, br));
                    stack.push((ar, bl));
                    stack.push((al, br));
                    stack.push((al, bl));
                }
            }
        }
        None
    }

    /// Visits every triangle whose box meets the closed query box.
    pub(crate) fn query(&self, query: &([f64; 3], [f64; 3]), mut visit: impl FnMut(u32)) {
        if self.nodes.is_empty() {
            return;
        }
        let mut stack = vec![0usize];
        while let Some(node) = stack.pop() {
            if !overlap(&self.node_box(node), query) {
                continue;
            }
            match self.children(node) {
                Some((l, r)) => {
                    stack.push(r);
                    stack.push(l);
                }
                None => {
                    for &t in self.leaf(node) {
                        if overlap(&self.boxes[t as usize], query) {
                            visit(t);
                        }
                    }
                }
            }
        }
    }

    /// Best-first search by a rigorous lower bound on squared distance from a
    /// box to the query. `exact` returns an upper bound for one triangle and
    /// may lower `best`; nodes whose lower bound exceeds `best` are skipped.
    pub(crate) fn nearest(
        &self,
        lower: impl Fn(&([f64; 3], [f64; 3])) -> f64,
        mut exact: impl FnMut(u32) -> f64,
    ) {
        if self.nodes.is_empty() {
            return;
        }
        let mut best = f64::INFINITY;
        let mut heap = std::collections::BinaryHeap::new();
        heap.push(Pending(lower(&self.node_box(0)), 0));
        while let Some(Pending(bound, node)) = heap.pop() {
            if bound > best {
                break;
            }
            match self.children(node) {
                Some((l, r)) => {
                    for child in [l, r] {
                        let b = lower(&self.node_box(child));
                        if b <= best {
                            heap.push(Pending(b, child));
                        }
                    }
                }
                None => {
                    let mut leaf: Vec<(f64, u32)> = self
                        .leaf(node)
                        .iter()
                        .map(|&t| (lower(&self.boxes[t as usize]), t))
                        .collect();
                    leaf.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
                    for (b, t) in leaf {
                        if b <= best {
                            best = best.min(exact(t));
                        }
                    }
                }
            }
        }
    }
}

/// Min-heap entry on the lower bound, ties broken by node index.
struct Pending(f64, usize);
impl PartialEq for Pending {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}
impl Eq for Pending {}
impl PartialOrd for Pending {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Pending {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .0
            .total_cmp(&self.0)
            .then_with(|| other.1.cmp(&self.1))
    }
}
