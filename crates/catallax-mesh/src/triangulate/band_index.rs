use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Decider, Interval, Refusal};
use spade::Point2;

type Point = Point2<f64>;
pub(super) type Segment = [Point; 2];

#[derive(Clone, Copy)]
pub(super) struct Bounds([Interval; 2]);

impl Bounds {
    fn point(at: Point) -> Self {
        Self([at.x, at.y].map(|v| Interval::point(v).unwrap()))
    }

    fn union(self, other: Self) -> Self {
        Self(std::array::from_fn(|i| self.0[i].hull(&other.0[i])))
    }

    pub(super) fn segment(segment: Segment) -> Self {
        Self::point(segment[0]).union(Self::point(segment[1]))
    }
}

fn disk_dot(at: Bounds, a: Bounds, b: Bounds) -> Interval {
    let [x, y] = at.0;
    (x - a.0[0]) * (x - b.0[0]) + (y - a.0[1]) * (y - b.0[1])
}

fn clears_segment(at: Point, [a, b]: Segment) -> bool {
    if (at.x <= a.x.min(b.x) || at.x >= a.x.max(b.x))
        && (at.y <= a.y.min(b.y) || at.y >= a.y.max(b.y))
    {
        return true;
    }
    disk_dot(Bounds::point(at), Bounds::point(a), Bounds::point(b)).lo() >= 0.0
}

struct Entry<T> {
    bounds: Bounds,
    item: T,
}

enum Children<T> {
    Leaf(T),
    Branch(Box<[Node<T>; 2]>),
}

struct Node<T> {
    bounds: Bounds,
    children: Children<T>,
}

impl<T: Copy> Node<T> {
    fn build(entries: &mut [Entry<T>], dec: &mut Decider) -> Result<Self, Refusal> {
        dec.charge_traversal(1)?;
        if let [entry] = entries {
            return Ok(Self {
                bounds: entry.bounds,
                children: Children::Leaf(entry.item),
            });
        }
        dec.charge_traversal(entries.len() as u64)?;
        let mut lo = [f64::INFINITY; 2];
        let mut hi = [f64::NEG_INFINITY; 2];
        for entry in entries.iter() {
            for i in 0..2 {
                let coordinate = entry.bounds.0[i].mid();
                lo[i] = lo[i].min(coordinate);
                hi[i] = hi[i].max(coordinate);
            }
        }
        let axis = usize::from(hi[1] - lo[1] > hi[0] - lo[0]);
        let mid = entries.len() / 2;
        let mut comparisons = 0;
        entries.select_nth_unstable_by(mid, |a, b| {
            comparisons += 1;
            let key = |item: &Entry<T>, i: usize| item.bounds.0[i].mid();
            key(a, axis)
                .total_cmp(&key(b, axis))
                .then_with(|| key(a, 1 - axis).total_cmp(&key(b, 1 - axis)))
        });
        dec.charge_traversal(comparisons)?;
        let (left, right) = entries.split_at_mut(mid);
        let left = Self::build(left, dec)?;
        let right = Self::build(right, dec)?;
        Ok(Self {
            bounds: left.bounds.union(right.bounds),
            children: Children::Branch(Box::new([left, right])),
        })
    }

    fn any(
        &self,
        possible: &impl Fn(Bounds) -> bool,
        matches: &impl Fn(T) -> bool,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        dec.charge_traversal(1)?;
        if !possible(self.bounds) {
            return Ok(false);
        }
        match &self.children {
            Children::Leaf(item) => {
                dec.charge_traversal(1)?;
                Ok(matches(*item))
            }
            Children::Branch(children) => Ok(children[0].any(possible, matches, dec)?
                || children[1].any(possible, matches, dec)?),
        }
    }
}

pub(super) struct Index<T>(Option<Node<T>>);

impl<T: Copy> Index<T> {
    fn build(
        items: impl IntoIterator<Item = T>,
        bounds: impl Fn(T) -> Bounds,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut entries = Vec::new();
        for item in items {
            dec.charge_traversal(1)?;
            entries.push(Entry {
                bounds: bounds(item),
                item,
            });
        }
        let root = if entries.is_empty() {
            None
        } else {
            Some(Node::build(&mut entries, dec)?)
        };
        Ok(Self(root))
    }

    fn any(
        &self,
        possible: impl Fn(Bounds) -> bool,
        matches: impl Fn(T) -> bool,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        self.0
            .as_ref()
            .map_or(Ok(false), |root| root.any(&possible, &matches, dec))
    }
}

impl Index<Segment> {
    pub(super) fn segments(items: Vec<Segment>, dec: &mut Decider) -> Result<Self, Refusal> {
        Self::build(items, Bounds::segment, dec)
    }

    pub(super) fn clears(&self, at: Point, dec: &mut Decider) -> Result<bool, Refusal> {
        // Each endpoint lies in the box. A nonnegative lower bound covers
        // every segment in the node, including the rounding in the leaf test.
        Ok(!self.any(
            |bounds| disk_dot(Bounds::point(at), bounds, bounds).lo() < 0.0,
            |segment| !clears_segment(at, segment),
            dec,
        )?)
    }
}

impl Index<Point> {
    pub(super) fn anchors(
        items: impl IntoIterator<Item = Point>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        Self::build(items, Bounds::point, dec)
    }

    pub(super) fn clears(&self, segment: Segment, dec: &mut Decider) -> Result<bool, Refusal> {
        let [a, b] = segment.map(Bounds::point);
        Ok(!self.any(
            |bounds| disk_dot(bounds, a, b).lo() < 0.0,
            |at| !clears_segment(at, segment),
            dec,
        )?)
    }
}

// Ordinary comparison treats the two zeros as equal. The integer key has
// that same equality and orders negative coordinates before positive ones.
fn ordered(value: f64) -> u64 {
    let sign = 1 << 63;
    let bits = value.to_bits();
    if value == 0.0 {
        sign
    } else if bits & sign == 0 {
        bits | sign
    } else {
        !bits
    }
}

#[derive(Default)]
pub(super) struct Points(BTreeMap<u64, BTreeSet<u64>>);

impl Points {
    pub(super) fn insert(&mut self, at: Point, dec: &mut Decider) -> Result<(), Refusal> {
        dec.charge_traversal(1)?;
        self.0
            .entry(ordered(at.x))
            .or_default()
            .insert(ordered(at.y));
        Ok(())
    }

    pub(super) fn separated(
        &self,
        at: Point,
        radius: f64,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        dec.charge_traversal(1)?;
        let range = |value: f64| {
            ordered((value - radius).next_down())..=ordered((value + radius).next_up())
        };
        let rows = range(at.y);
        for (_, column) in self.0.range(range(at.x)) {
            dec.charge_traversal(1)?;
            if column.range(rows.clone()).next().is_some() {
                dec.charge_traversal(1)?;
                return Ok(false);
            }
        }
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    #[test]
    fn indexed_trim_matches_a_linear_check_and_charges_visits() {
        let segments: Vec<_> = (0..1024)
            .map(|i| {
                let x = f64::from(i) * 4.0;
                [
                    Point::new(x, 0.0),
                    Point::new(x + 1.0, f64::from(i % 3) / 2.0),
                ]
            })
            .collect();
        let mut dec = Decider::default();
        let index = Index::segments(segments.clone(), &mut dec).unwrap();
        for x in [0.5, 19.0, 1800.0, 4095.5] {
            for y in [-0.25, 0.0, 0.25, 1.125] {
                let at = Point::new(x, y);
                let (clear, work) = with_operation_budget(60, || index.clears(at, &mut dec));
                assert_eq!(
                    clear.unwrap(),
                    segments.iter().all(|&segment| clears_segment(at, segment))
                );
                assert!(!work.exhausted);
                assert!(work.used > 0);
            }
        }
        let (clear, work) =
            with_operation_budget(1, || index.clears(Point::new(-3.0, -3.0), &mut dec));
        assert!(clear.unwrap());
        assert_eq!(work.used, 1);
        assert!(!work.exhausted);
        let (refusal, _) =
            with_operation_budget(1, || index.clears(Point::new(0.5, 0.0), &mut dec));
        assert_eq!(refusal.unwrap_err().code, "numbers.operation-budget");
    }

    #[test]
    fn point_index_checks_insertions_at_the_outward_rounded_limit() {
        let mut dec = Decider::default();
        let mut points = Points::default();
        let at = Point::new(-0.0, 1.0);
        let radius = 0.125;
        assert!(points.separated(at, radius, &mut dec).unwrap());
        points
            .insert(Point::new(0.0, 1.125_f64.next_up()), &mut dec)
            .unwrap();
        assert!(!points.separated(at, radius, &mut dec).unwrap());
        assert!(
            points
                .separated(Point::new(0.0, 0.0), radius, &mut dec)
                .unwrap()
        );
    }

    #[test]
    fn tilted_segment_and_anchor_queries_include_overhang() {
        let mut dec = Decider::default();
        let segment = [Point::new(0.0, 0.0), Point::new(2.0, 2.0)];
        let at = Point::new(2.125, 0.5);
        let segments = Index::segments(vec![segment], &mut dec).unwrap();
        let anchors = Index::anchors([at], &mut dec).unwrap();
        assert!(!segments.clears(at, &mut dec).unwrap());
        assert!(!anchors.clears(segment, &mut dec).unwrap());
    }
}
