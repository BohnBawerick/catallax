//! Certify the affine binary64-to-binary32 vertex motion, not just its endpoints.
//! Separating planes keep every triangle pair disjoint except its shared simplex.
//! Bernstein coefficients bound the complete motion parameter in [0, 1].
use catallax_numbers::{Decider, Rational, Refusal};
use std::collections::{BTreeMap, BTreeSet};

type Vector = [Polynomial; 3];
#[derive(Clone)]
struct Polynomial([Rational; 5]);
impl Polynomial {
    fn constant(x: f64) -> Self {
        let mut p = Self::zero();
        p.0[0] = Rational::from_f64_exact(x).expect("finite mesh coordinate");
        p
    }
    fn zero() -> Self {
        Self(std::array::from_fn(|_| Rational::zero()))
    }
    fn add(&self, b: &Self) -> Self {
        Self(std::array::from_fn(|i| &self.0[i] + &b.0[i]))
    }
    fn sub(&self, b: &Self) -> Self {
        Self(std::array::from_fn(|i| &self.0[i] - &b.0[i]))
    }
    fn mul(&self, b: &Self) -> Self {
        let mut p = Self::zero();
        for i in 0..5 {
            if self.0[i].is_zero() {
                continue;
            }
            for j in 0..5 - i {
                if !b.0[j].is_zero() {
                    p.0[i + j] = &p.0[i + j] + &(&self.0[i] * &b.0[j]);
                }
            }
        }
        p
    }
    // 2 = strictly positive, 1 = nonnegative, 0 = identically zero,
    // -1 = nonpositive, -2 = strictly negative, 3 = unresolved.
    fn sign(&self) -> i8 {
        const CHOOSE: [[i64; 5]; 5] = [
            [1, 0, 0, 0, 0],
            [1, 1, 0, 0, 0],
            [1, 2, 1, 0, 0],
            [1, 3, 3, 1, 0],
            [1, 4, 6, 4, 1],
        ];
        let n = self.0.iter().rposition(|v| !v.is_zero()).unwrap_or(0);
        let b: [Rational; 5] = std::array::from_fn(|k| {
            if k > n {
                return Rational::zero();
            }
            (0..=k).fold(Rational::zero(), |v, i| {
                &v + &(&self.0[i] * &Rational::from_ratio(CHOOSE[k][i], CHOOSE[n][i]).unwrap())
            })
        });
        let zero = Rational::zero();
        let b = &b[..=n];
        if b.iter().all(|x| x > &zero) {
            2
        } else if b.iter().all(|x| x < &zero) {
            -2
        } else if b.iter().all(|x| x == &zero) {
            0
        } else if b.iter().all(|x| x >= &zero) {
            1
        } else if b.iter().all(|x| x <= &zero) {
            -1
        } else {
            3
        }
    }
}
fn sub(a: &Vector, b: &Vector) -> Vector {
    std::array::from_fn(|i| a[i].sub(&b[i]))
}
fn add(a: &Vector, b: &Vector) -> Vector {
    std::array::from_fn(|i| a[i].add(&b[i]))
}
fn cross(a: &Vector, b: &Vector) -> Vector {
    std::array::from_fn(|i| {
        a[(i + 1) % 3]
            .mul(&b[(i + 2) % 3])
            .sub(&a[(i + 2) % 3].mul(&b[(i + 1) % 3]))
    })
}
fn dot(a: &Vector, b: &Vector) -> Polynomial {
    (0..3).fold(Polynomial::zero(), |v, i| v.add(&a[i].mul(&b[i])))
}
fn failure(detail: impl Into<String>) -> Refusal {
    Refusal::new("print.brep-encoding-topology", detail)
}

pub(super) fn topology(
    triangles: &[[usize; 3]],
    count: usize,
    dec: &Decider,
) -> Result<(), Refusal> {
    let mut edges = BTreeMap::<(usize, usize), (usize, i32)>::new();
    let mut links = vec![Vec::new(); count];
    for &[a, b, c] in triangles {
        dec.check_deadline()?;
        for (v, x, y) in [(a, b, c), (b, c, a), (c, a, b)] {
            if v == x || x == y || y == v || v >= count || x >= count || y >= count {
                return Err(failure("invalid facet indices"));
            }
            links[v].push((x, y));
            let entry = edges.entry((v.min(x), v.max(x))).or_default();
            entry.0 += 1;
            entry.1 += if v < x { 1 } else { -1 };
        }
    }
    if edges.values().any(|&e| e != (2, 0)) {
        return Err(failure("open or non-oppositely-paired facet edge"));
    }
    for link in links {
        dec.check_deadline()?;
        let mut adjacency = BTreeMap::<usize, Vec<usize>>::new();
        for (a, b) in link {
            adjacency.entry(a).or_default().push(b);
            adjacency.entry(b).or_default().push(a);
        }
        if adjacency.is_empty() || adjacency.values().any(|v| v.len() != 2) {
            return Err(failure("vertex link is not a cycle"));
        }
        let mut pending = vec![*adjacency.first_key_value().unwrap().0];
        let mut seen = BTreeSet::new();
        while let Some(v) = pending.pop() {
            if seen.insert(v) {
                pending.extend(&adjacency[&v]);
            }
        }
        if seen.len() != adjacency.len() {
            return Err(failure("vertex link joins disconnected sheets"));
        }
    }
    Ok(())
}

fn separated(a: &[usize; 3], b: &[usize; 3], p: &[Vector], dec: &Decider) -> Result<bool, Refusal> {
    let ea: [Vector; 3] = std::array::from_fn(|i| sub(&p[a[(i + 1) % 3]], &p[a[i]]));
    let eb: [Vector; 3] = std::array::from_fn(|i| sub(&p[b[(i + 1) % 3]], &p[b[i]]));
    let na = cross(&ea[0], &ea[1]);
    let nb = cross(&eb[0], &eb[1]);
    let mut shared = [0; 3];
    let mut count = 0;
    for &vertex in a {
        if b.contains(&vertex) {
            shared[count] = vertex;
            count += 1;
        }
    }
    let common = &shared[..count];
    if common.len() == 3 {
        return Err(failure("duplicate facets"));
    }
    let separates = |axis: &Vector| -> Result<bool, Refusal> {
        dec.check_deadline()?;
        if common.is_empty() {
            let mut positive = true;
            let mut negative = true;
            for &i in a {
                for &j in b {
                    let s = dot(axis, &sub(&p[j], &p[i])).sign();
                    positive &= s == 2;
                    negative &= s == -2;
                    if !positive && !negative {
                        return Ok(false);
                    }
                }
            }
            Ok(positive || negative)
        } else {
            let anchor = &p[common[0]];
            if common
                .iter()
                .skip(1)
                .any(|&i| dot(axis, &sub(&p[i], anchor)).sign() != 0)
            {
                return Ok(false);
            }
            for direction in [1, -1] {
                let mut strict_a = true;
                let mut strict_b = true;
                let mut valid = true;
                for &i in a.iter().filter(|i| !common.contains(i)) {
                    let s = dot(axis, &sub(&p[i], anchor)).sign() * direction;
                    valid &= matches!(s, 0..=2);
                    strict_a &= s == 2;
                }
                for &i in b.iter().filter(|i| !common.contains(i)) {
                    let s = dot(axis, &sub(&p[i], anchor)).sign() * direction;
                    valid &= matches!(s, -2..=0);
                    strict_b &= s == -2;
                }
                if valid && (strict_a || strict_b) {
                    return Ok(true);
                }
            }
            Ok(false)
        }
    };
    if common.len() == 2 {
        let edge = sub(&p[common[1]], &p[common[0]]);
        if separates(&cross(&edge, &add(&na, &nb)))? || separates(&cross(&edge, &sub(&na, &nb)))? {
            return Ok(true);
        }
    }
    if separates(&na)? || separates(&nb)? {
        return Ok(true);
    }
    for x in &ea {
        for y in &eb {
            if separates(&cross(x, y))? {
                return Ok(true);
            }
        }
    }
    // Coplanar triangles need in-plane separating axes too.
    for (normal, edges) in [(&na, &ea), (&nb, &eb)] {
        for edge in edges {
            if separates(&cross(normal, edge))? {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

#[derive(Clone, Copy)]
struct Bounds {
    lo: [f64; 3],
    hi: [f64; 3],
}
impl Bounds {
    fn overlaps(&self, b: &Self) -> bool {
        (0..3).all(|i| self.lo[i] <= b.hi[i] && b.lo[i] <= self.hi[i])
    }
    fn union(self, b: Self) -> Self {
        Self {
            lo: std::array::from_fn(|i| self.lo[i].min(b.lo[i])),
            hi: std::array::from_fn(|i| self.hi[i].max(b.hi[i])),
        }
    }
}
struct Node {
    bounds: Bounds,
    children: Option<(usize, usize)>,
    triangle: usize,
}
fn tree(
    ids: &mut [usize],
    bounds: &[Bounds],
    nodes: &mut Vec<Node>,
    dec: &Decider,
) -> Result<usize, Refusal> {
    dec.check_deadline()?;
    let total = ids
        .iter()
        .map(|&i| bounds[i])
        .reduce(Bounds::union)
        .unwrap();
    let id = nodes.len();
    nodes.push(Node {
        bounds: total,
        children: None,
        triangle: ids[0],
    });
    if ids.len() > 1 {
        let axis = (0..3)
            .max_by(|&i, &j| (total.hi[i] - total.lo[i]).total_cmp(&(total.hi[j] - total.lo[j])))
            .unwrap();
        ids.sort_unstable_by(|&a, &b| {
            (bounds[a].lo[axis] + bounds[a].hi[axis])
                .total_cmp(&(bounds[b].lo[axis] + bounds[b].hi[axis]))
                .then(a.cmp(&b))
        });
        let (a, b) = ids.split_at_mut(ids.len() / 2);
        nodes[id].children = Some((tree(a, bounds, nodes, dec)?, tree(b, bounds, nodes, dec)?));
    }
    Ok(id)
}

pub(super) fn motion(
    original: &[[f64; 3]],
    encoded: &[[f32; 3]],
    triangles: &[[usize; 3]],
    dec: &Decider,
) -> Result<(), Refusal> {
    let mut positions = BTreeSet::new();
    for point in encoded {
        dec.check_deadline()?;
        let key = point.map(|x| if x == 0.0 { 0 } else { x.to_bits() });
        if !positions.insert(key) {
            return Err(failure("binary32 rounding joins distinct vertices"));
        }
    }
    let p: Vec<Vector> = original
        .iter()
        .zip(encoded)
        .map(|(a, b)| {
            std::array::from_fn(|i| {
                let mut v = Polynomial::constant(a[i]);
                v.0[1] = &Rational::from_f64_exact(f64::from(b[i])).unwrap() - &v.0[0];
                v
            })
        })
        .collect();
    for t in triangles {
        dec.check_deadline()?;
        let n = cross(&sub(&p[t[1]], &p[t[0]]), &sub(&p[t[2]], &p[t[0]]));
        let initial: Vector = std::array::from_fn(|i| {
            let mut v = Polynomial::zero();
            v.0[0] = n[i].0[0].clone();
            v
        });
        if dot(&n, &initial).sign() != 2 {
            return Err(failure(
                "facet collapses or reverses during binary32 encoding",
            ));
        }
    }
    if triangles.is_empty() {
        return Ok(());
    }
    let bounds: Vec<_> = triangles
        .iter()
        .map(|t| Bounds {
            lo: std::array::from_fn(|axis| {
                t.iter()
                    .map(|&i| original[i][axis].min(f64::from(encoded[i][axis])))
                    .fold(f64::INFINITY, f64::min)
            }),
            hi: std::array::from_fn(|axis| {
                t.iter()
                    .map(|&i| original[i][axis].max(f64::from(encoded[i][axis])))
                    .fold(f64::NEG_INFINITY, f64::max)
            }),
        })
        .collect();
    let mut nodes = Vec::with_capacity(triangles.len() * 2);
    let mut ids: Vec<_> = (0..triangles.len()).collect();
    let root = tree(&mut ids, &bounds, &mut nodes, dec)?;
    let mut pending = vec![(root, root)];
    while let Some((a, b)) = pending.pop() {
        dec.check_deadline()?;
        let (a_node, b_node) = (&nodes[a], &nodes[b]);
        if !a_node.bounds.overlaps(&b_node.bounds) {
            continue;
        }
        if a == b {
            if let Some((x, y)) = a_node.children {
                pending.extend([(x, x), (x, y), (y, y)]);
            }
        } else if let Some((x, y)) = a_node.children {
            pending.extend([(x, b), (y, b)]);
        } else if let Some((x, y)) = b_node.children {
            pending.extend([(a, x), (a, y)]);
        } else {
            let (x, y) = (a_node.triangle, b_node.triangle);
            if !separated(&triangles[x], &triangles[y], &p, dec)? {
                return Err(failure(format!(
                    "no complete separation proof for facets {x} and {y} during encoding"
                )));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rounding_must_not_merge_separate_sheets_or_collapse_facets() {
        let original = [[1.0, 0.0, 0.0], [1.0 + 1e-9, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let encoded = original.map(|p| p.map(|v| v as f32));
        assert_eq!(
            motion(&original, &encoded, &[[0, 1, 2]], &Decider::default())
                .unwrap_err()
                .code,
            "print.brep-encoding-topology"
        );
    }
    #[test]
    fn disjoint_indices_do_not_authorize_crossing_triangles() {
        let original = [
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.25, -1.0],
            [0.0, 0.25, 1.0],
            [0.0, 0.75, 0.0],
        ];
        let encoded = original.map(|p| p.map(|v| v as f32));
        assert!(
            motion(
                &original,
                &encoded,
                &[[0, 1, 2], [3, 4, 5]],
                &Decider::default()
            )
            .is_err()
        );
    }
    #[test]
    fn vertex_contact_between_closed_sheets_is_not_a_manifold() {
        let tetrahedron = [[0, 1, 2], [0, 3, 1], [1, 3, 2], [0, 2, 3]];
        let mut triangles = tetrahedron.to_vec();
        triangles.extend(tetrahedron.map(|t| t.map(|i| if i == 0 { 0 } else { i + 3 })));
        assert!(topology(&triangles, 7, &Decider::default()).is_err());
    }
    #[test]
    fn rounding_may_join_facet_interiors_without_merging_vertices() {
        let original = [
            [0.0, 0.0, 1.0 - 1e-9],
            [2.0, 0.0, 1.0 - 1e-9],
            [0.0, 2.0, 1.0 - 1e-9],
            [0.25, 0.25, 1.0 + 1e-9],
            [2.25, 0.25, 1.0 + 1e-9],
            [0.25, 2.25, 1.0 + 1e-9],
        ];
        let encoded = original.map(|p| p.map(|v| v as f32));
        assert!(
            motion(
                &original,
                &encoded,
                &[[0, 1, 2], [3, 4, 5]],
                &Decider::default()
            )
            .is_err()
        );
    }
}
