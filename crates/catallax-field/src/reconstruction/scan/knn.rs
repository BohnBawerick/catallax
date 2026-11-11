//! A k-d tree for nearest-neighbour graphs over scan samples. Splits and ties
//! are broken by sample index, so the neighbours of a point depend only on the
//! input, never on thread count or memory layout.

use super::linalg::V3;

const LEAF: usize = 8;

pub(super) struct KdTree<'a> {
    points: &'a [V3],
    order: Vec<u32>,
    /// Inner nodes: (axis, split value, left end, right start); leaves: ranges.
    nodes: Vec<Node>,
}

enum Node {
    Leaf(usize, usize),
    Inner {
        axis: usize,
        split: f64,
        left: usize,
        right: usize,
    },
}

impl<'a> KdTree<'a> {
    pub(super) fn new(points: &'a [V3]) -> Self {
        let mut tree = Self {
            points,
            order: (0..points.len() as u32).collect(),
            nodes: Vec::new(),
        };
        if !points.is_empty() {
            tree.build(0, points.len());
        }
        tree
    }

    fn build(&mut self, start: usize, end: usize) -> usize {
        let index = self.nodes.len();
        if end - start <= LEAF {
            self.nodes.push(Node::Leaf(start, end));
            return index;
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &i in &self.order[start..end] {
            let p = self.points[i as usize];
            for a in 0..3 {
                lo[a] = lo[a].min(p[a]);
                hi[a] = hi[a].max(p[a]);
            }
        }
        let mut axis = 0;
        for a in 1..3 {
            if hi[a] - lo[a] > hi[axis] - lo[axis] {
                axis = a;
            }
        }
        let middle = start + (end - start) / 2;
        let points = self.points;
        self.order[start..end].select_nth_unstable_by(middle - start, |&i, &j| {
            points[i as usize][axis]
                .total_cmp(&points[j as usize][axis])
                .then(i.cmp(&j))
        });
        let split = points[self.order[middle] as usize][axis];
        self.nodes.push(Node::Leaf(0, 0));
        let left = self.build(start, middle);
        let right = self.build(middle, end);
        self.nodes[index] = Node::Inner {
            axis,
            split,
            left,
            right,
        };
        index
    }

    /// The `k` nearest samples to sample `i`, itself first, ordered by squared
    /// distance and then by index.
    pub(super) fn nearest(&self, i: usize, k: usize) -> Vec<u32> {
        let q = self.points[i];
        let mut best: Vec<(f64, u32)> = Vec::with_capacity(k + 1);
        if !self.nodes.is_empty() {
            self.search(0, q, k, &mut best);
        }
        best.into_iter().map(|(_, j)| j).collect()
    }

    fn search(&self, node: usize, q: V3, k: usize, best: &mut Vec<(f64, u32)>) {
        match self.nodes[node] {
            Node::Leaf(start, end) => {
                for &j in &self.order[start..end] {
                    let p = self.points[j as usize];
                    let d = (0..3).map(|a| (p[a] - q[a]).powi(2)).sum::<f64>();
                    let before = |e: &(f64, u32)| d.total_cmp(&e.0).then(j.cmp(&e.1)).is_lt();
                    if best.len() < k || before(&best[k - 1]) {
                        let at = best.partition_point(|e| !before(e));
                        best.insert(at, (d, j));
                        best.truncate(k);
                    }
                }
            }
            Node::Inner {
                axis,
                split,
                left,
                right,
            } => {
                let gap = q[axis] - split;
                let (near, far) = if gap < 0.0 {
                    (left, right)
                } else {
                    (right, left)
                };
                self.search(near, q, k, best);
                // Equal coordinates may sit on either side of the split.
                if best.len() < k || gap * gap <= best[k - 1].0 {
                    self.search(far, q, k, best);
                }
            }
        }
    }
}
