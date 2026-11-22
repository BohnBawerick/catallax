//! Exact planar cells of the affine field. As in OpenCASCADE and SMLib's
//! trimmed-face tessellators, triangulate only after sharing boundary edges.
//! https://dev.opencascade.org/doc/overview/html/occt_user_guides__mesh.html
//! https://github.com/NVIDIA-Omniverse/usd-brep/blob/main/.agents/operations/tessellation.md

use super::{Box3, Budget, Model, Point, Support, q, refusal};
use crate::Instruction;
use catallax_numbers::{Decider, Interval, Rational, Refusal, Sign};
use catallax_recipe::BooleanKind;
use std::collections::{BTreeMap, BTreeSet};

mod regions;

type Polygon = Vec<Point>;
type Triangles = (Vec<Point>, Vec<[usize; 3]>);

fn difference(a: &Point, b: &Point) -> Point {
    std::array::from_fn(|i| &a[i] - &b[i])
}
fn dot(a: &Point, b: &Point) -> Rational {
    (0..3).fold(q(0), |sum, i| &sum + &(&a[i] * &b[i]))
}
fn cross(a: &Point, b: &Point) -> Point {
    std::array::from_fn(|i| {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        &(&a[j] * &b[k]) - &(&a[k] * &b[j])
    })
}
fn nonzero(p: &Point) -> bool {
    p.iter().any(|v| !v.is_zero())
}
fn box_of<'a>(points: impl Iterator<Item = &'a Point>) -> Box3 {
    let mut points = points;
    let first = points.next().expect("nonempty polygon");
    let mut bounds = first.clone().map(|x| [x.clone(), x]);
    for p in points {
        for i in 0..3 {
            if p[i] < bounds[i][0] {
                bounds[i][0] = p[i].clone();
            }
            if p[i] > bounds[i][1] {
                bounds[i][1] = p[i].clone();
            }
        }
    }
    bounds
}
fn overlaps(a: &Box3, b: &Box3) -> bool {
    (0..3).all(|i| a[i][0] <= b[i][1] && b[i][0] <= a[i][1])
}

struct Plane {
    normal: Point,
    offset: Rational,
    enclosed_normal: [Interval; 3],
    enclosed_offset: Interval,
}
impl Plane {
    fn new(poly: &[Point]) -> Self {
        let normal = cross(
            &difference(&poly[1], &poly[0]),
            &difference(&poly[2], &poly[0]),
        );
        assert!(nonzero(&normal));
        let offset = dot(&normal, &poly[0]);
        Self {
            enclosed_offset: Interval::from_rational(&offset),
            enclosed_normal: normal.each_ref().map(Interval::from_rational),
            offset,
            normal,
        }
    }
    fn value(&self, p: &Point) -> Rational {
        &dot(&self.normal, p) - &self.offset
    }
    fn sign(&self, p: &Point, enclosure: &[Interval; 3]) -> Sign {
        let value = (0..3).fold(super::negative(self.enclosed_offset), |v, i| {
            v + self.enclosed_normal[i] * enclosure[i]
        });
        if value.lo() > 0.0 {
            Sign::Positive
        } else if value.hi() < 0.0 {
            Sign::Negative
        } else {
            self.value(p).sign()
        }
    }
}

struct Solid {
    leaf: usize,
    bounds: Box3,
    faces: Vec<Polygon>,
    planes: Vec<Plane>,
}
impl Solid {
    fn from_faces(leaf: usize, faces: Vec<Polygon>) -> Self {
        let bounds = box_of(faces.iter().flatten());
        let planes = faces.iter().map(|face| Plane::new(face)).collect();
        Self {
            leaf,
            bounds,
            faces,
            planes,
        }
    }
    fn block(leaf: usize, bounds: &Box3) -> Self {
        let points: [Point; 8] =
            std::array::from_fn(|bits| std::array::from_fn(|i| bounds[i][(bits >> i) & 1].clone()));
        let faces = (0..3)
            .flat_map(|axis| {
                let a = 1 << axis;
                let j = 1 << ((axis + 1) % 3);
                let k = 1 << ((axis + 2) % 3);
                [[0, k, j + k, j], [a, a + j, a + j + k, a + k]]
            })
            .map(|idx| idx.map(|i| points[i].clone()).to_vec())
            .collect();
        Self::from_faces(leaf, faces)
    }
    fn prism(
        leaf: usize,
        axis: usize,
        local: [&Rational; 2],
        axial: [&Rational; 2],
        ring: &[[Rational; 2]],
    ) -> Self {
        let point = |p: &[Rational; 2], z: &Rational| {
            let mut value = [q(0), q(0), q(0)];
            value[axis] = z.clone();
            value[(axis + 1) % 3] = local[0] + &p[0];
            value[(axis + 2) % 3] = local[1] + &p[1];
            value
        };
        let low: Polygon = ring.iter().map(|p| point(p, axial[0])).collect();
        let high: Polygon = ring.iter().map(|p| point(p, axial[1])).collect();
        let mut faces = vec![low.iter().rev().cloned().collect(), high.clone()];
        for i in 0..ring.len() {
            let j = (i + 1) % ring.len();
            faces.push(vec![
                low[i].clone(),
                low[j].clone(),
                high[j].clone(),
                high[i].clone(),
            ]);
        }
        Self::from_faces(leaf, faces)
    }
    // Membership at a point infinitesimally before/after the candidate face
    // along its normal. Exact ties use the linear directional coefficient.
    fn sides(&self, p: &Point, normal: &Point, dec: &mut Decider) -> Result<[bool; 2], Refusal> {
        let mut inside = [true; 2];
        let enclosure = p.each_ref().map(Interval::from_rational);
        for plane in &self.planes {
            dec.charge_work(1)?;
            match plane.sign(p, &enclosure) {
                Sign::Positive => return Ok([false; 2]),
                Sign::Zero => {
                    let slope = dot(&plane.normal, normal);
                    inside[0] &= slope >= q(0);
                    inside[1] &= slope <= q(0);
                }
                Sign::Negative => {}
            }
        }
        Ok(inside)
    }
}

// Subdividing at both transverse knot sets traces the complete level polygon
// of L(x²)+L(y²)=r². Between knots the residual is affine.
fn ring(
    radius: &Rational,
    pitch: &Rational,
    budget: &mut Budget<'_>,
    dec: &mut Decider,
) -> Result<Vec<[Rational; 2]>, Refusal> {
    let r2 = radius * radius;
    let mut knots = vec![q(0)];
    while knots.last().expect("initial knot") < radius {
        // Each transverse knot gives at most eight reflected endpoints.
        budget.visit_many(8, dec)?;
        knots.push(knots.last().unwrap() + pitch);
    }
    let squares: Vec<_> = knots.iter().map(|v| v * v).collect();
    let mut radial = BTreeSet::new();
    for (x, x2) in knots.iter().zip(&squares) {
        if x2 > &r2 {
            break;
        }
        let y2 = &r2 - x2;
        let i = squares.partition_point(|value| value < &y2);
        let y = if squares[i] == y2 {
            knots[i].clone()
        } else {
            &knots[i - 1] + &(&(&y2 - &squares[i - 1]) / &(&knots[i - 1] + &knots[i]))
        };
        for a in [x.clone(), -x] {
            for b in [y.clone(), -&y] {
                radial.insert([a.clone(), b.clone()]);
                radial.insert([b, a.clone()]);
            }
        }
    }
    let mut polygon: Vec<_> = radial.into_iter().collect();
    let upper = |p: &[Rational; 2]| p[1] > q(0) || p[1].is_zero() && p[0] >= q(0);
    polygon.sort_by(|a, b| {
        upper(b)
            .cmp(&upper(a))
            .then_with(|| (&(&b[0] * &a[1]) - &(&b[1] * &a[0])).cmp(&q(0)))
    });
    // At e.g. x = 0 the quadratic interpolant can have a flat polygon side
    // broken by a knot of the other coordinate. Merge straight consecutive edges.
    let mut corners = Vec::new();
    for (i, p) in polygon.iter().enumerate() {
        let a = &polygon[(i + polygon.len() - 1) % polygon.len()];
        let b = &polygon[(i + 1) % polygon.len()];
        let first = &(&p[0] - &a[0]) * &(&b[1] - &a[1]);
        let second = &(&p[1] - &a[1]) * &(&b[0] - &a[0]);
        if first != second {
            corners.push(p.clone());
        }
    }
    Ok(corners)
}
fn supports(
    model: &Model<'_>,
    regions: &[Vec<Box3>],
    budget: &mut Budget<'_>,
    dec: &mut Decider,
) -> Result<Vec<Solid>, Refusal> {
    let mut result = Vec::new();
    for (i, support) in model.supports.iter().enumerate() {
        match support {
            None => {}
            Some(Support::Block(b)) => {
                budget.visit_many(6, dec)?;
                result.push(Solid::block(i, b));
            }
            Some(Support::Cylinder {
                centre,
                radius,
                height,
            }) => {
                let ring = ring(radius, &model.pitch, budget, dec)?;
                budget.visit_many(ring.len() + 2, dec)?;
                result.push(Solid::prism(
                    i,
                    2,
                    [&centre[0], &centre[1]],
                    [&centre[2], &(&centre[2] + height)],
                    &ring,
                ));
            }
            Some(Support::Cubic {
                centre,
                radius,
                period,
            }) => {
                let ring = ring(radius, &model.pitch, budget, dec)?;
                let mut strands = BTreeMap::<(usize, i64, i64), Vec<[Rational; 2]>>::new();
                for zone in &regions[i] {
                    let mut indices = [[0; 2]; 3];
                    for j in 0..3 {
                        let lower = &(&(&zone[j][0] - radius) - &centre[j]) / period;
                        let upper = &(&(&zone[j][1] + radius) - &centre[j]) / period;
                        let outside = || {
                            refusal(
                                "field.print-budget",
                                "lattice index exceeds the print range",
                            )
                        };
                        indices[j] = [
                            lower.ceil().to_i64().ok_or_else(outside)?,
                            upper.floor().to_i64().ok_or_else(outside)?,
                        ];
                    }
                    for axis in 0..3 {
                        let [j, k] = [(axis + 1) % 3, (axis + 2) % 3];
                        for a in indices[j][0]..=indices[j][1] {
                            for b in indices[k][0]..=indices[k][1] {
                                budget.visit(dec)?;
                                strands.entry((axis, a, b)).or_default().push([
                                    &zone[axis][0] - &model.pitch,
                                    &zone[axis][1] + &model.pitch,
                                ]);
                            }
                        }
                    }
                }
                for ((axis, a, b), mut spans) in strands {
                    spans.sort_unstable();
                    let mut merged: Vec<[Rational; 2]> = Vec::new();
                    for [low, high] in spans {
                        if let Some(last) = merged.last_mut().filter(|span| span[1] >= low) {
                            last[1] = last[1].clone().max(high);
                        } else {
                            merged.push([low, high]);
                        }
                    }
                    for [low, high] in merged {
                        budget.visit_many(ring.len() + 2, dec)?;
                        let at = |j, n| &centre[j] + &(period * &q(n));
                        let [j, k] = [(axis + 1) % 3, (axis + 2) % 3];
                        result.push(Solid::prism(
                            i,
                            axis,
                            [&at(j, a), &at(k, b)],
                            [&low, &high],
                            &ring,
                        ));
                    }
                }
            }
        }
    }
    Ok(result)
}

fn append(poly: &mut Polygon, p: Point) {
    if poly.last() != Some(&p) {
        poly.push(p);
    }
}
fn split(poly: &[Point], plane: &Plane, dec: &mut Decider) -> Result<(Polygon, Polygon), Refusal> {
    dec.charge_work(poly.len() as u64)?;
    let signs: Vec<_> = poly
        .iter()
        .map(|p| plane.sign(p, &p.each_ref().map(Interval::from_rational)))
        .collect();
    if !signs.contains(&Sign::Positive) {
        return Ok((poly.to_vec(), vec![]));
    }
    if !signs.contains(&Sign::Negative) {
        return Ok((vec![], poly.to_vec()));
    }
    let mut sides = [Vec::new(), Vec::new()];
    for i in 0..poly.len() {
        let j = (i + 1) % poly.len();
        let a = signs[i];
        let b = signs[j];
        if a != Sign::Positive {
            append(&mut sides[0], poly[i].clone());
        }
        if a != Sign::Negative {
            append(&mut sides[1], poly[i].clone());
        }
        if a != Sign::Zero && b != Sign::Zero && a != b {
            let a = plane.value(&poly[i]);
            let b = plane.value(&poly[j]);
            let t = &a / &(&a - &b);
            let p: Point = std::array::from_fn(|axis| {
                &poly[i][axis] + &(&(&poly[j][axis] - &poly[i][axis]) * &t)
            });
            append(&mut sides[0], p.clone());
            append(&mut sides[1], p);
        }
    }
    for side in &mut sides {
        if side.first() == side.last() {
            side.pop();
        }
    }
    let [a, b] = sides;
    Ok((a, b))
}

// A convex clip produces convex interior and exterior cells. Where the tool
// only touches an edge or a point the original polygon can be retained.
fn divide(
    poly: Polygon,
    solid: &Solid,
    budget: &mut Budget<'_>,
    dec: &mut Decider,
) -> Result<Vec<Polygon>, Refusal> {
    let mut inner = poly.clone();
    let mut outside = Vec::new();
    for plane in &solid.planes {
        let (new_inner, outer) = split(&inner, plane, dec)?;
        if new_inner.is_empty() {
            return Ok(vec![poly]);
        }
        if !outer.is_empty() {
            budget.visit(dec)?;
            outside.push(outer);
        }
        inner = new_inner;
    }
    budget.visit(dec)?;
    outside.push(inner);
    Ok(outside)
}

fn centroid(poly: &[Point]) -> Point {
    let n = q(poly.len() as i64);
    std::array::from_fn(|i| &poly.iter().fold(q(0), |s, p| &s + &p[i]) / &n)
}

fn emit(
    model: &Model<'_>,
    solids: &[Solid],
    near: &[usize],
    source: usize,
    normal: &Point,
    poly: &mut Polygon,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let p = centroid(poly);
    let mut leaves = vec![[false; 2]; model.supports.len()];
    leaves[solids[source].leaf][0] = true;
    for &i in near {
        if i == source {
            continue;
        }
        let solid = &solids[i];
        if !(0..3).all(|i| solid.bounds[i][0] <= p[i] && p[i] <= solid.bounds[i][1]) {
            continue;
        }
        let sides = solid.sides(&p, normal, dec)?;
        // Coincident finite faces use their first owner. Their vertex joins
        // are otherwise split by the same geometric supports.
        if i < source && sides[0] != sides[1] {
            return Ok(false);
        }
        for j in 0..2 {
            leaves[solid.leaf][j] |= sides[j];
        }
    }
    let mut tape: Vec<[bool; 2]> = Vec::with_capacity(model.plan.instructions.len());
    for instruction in &model.plan.instructions {
        let sides = match *instruction {
            Instruction::Empty => [false; 2],
            Instruction::Leaf(i) => leaves[i],
            Instruction::Boolean { kind, a, b } => std::array::from_fn(|i| match kind {
                BooleanKind::Union => tape[a][i] || tape[b][i],
                BooleanKind::Intersect => tape[a][i] && tape[b][i],
                BooleanKind::Subtract => tape[a][i] && !tape[b][i],
            }),
        };
        tape.push(sides);
    }
    let [inward, outward] = tape[model.plan.root];
    if !inward && outward {
        poly.reverse();
    }
    Ok(inward != outward)
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct Line {
    axis: usize,
    direction: Point,
    intercept: Point,
}
impl Line {
    fn through(a: &Point, b: &Point) -> Self {
        let d = difference(b, a);
        let axis = d
            .iter()
            .position(|v| !v.is_zero())
            .expect("nonzero polygon edge");
        let direction = d.each_ref().map(|v| v / &d[axis]);
        let intercept = std::array::from_fn(|i| &a[i] - &(&direction[i] * &a[axis]));
        Self {
            axis,
            direction,
            intercept,
        }
    }
    fn point(&self, t: &Rational) -> Point {
        std::array::from_fn(|i| &self.intercept[i] + &(&self.direction[i] * t))
    }
}

fn triangulate(
    polygons: Vec<Polygon>,
    budget: &mut Budget<'_>,
    dec: &mut Decider,
) -> Result<Triangles, Refusal> {
    // Extended clip planes can introduce a collinear partition on only one
    // side of an edge. Make all uses share the union of exact partition points.
    let mut lines: BTreeMap<Line, BTreeSet<Rational>> = BTreeMap::new();
    for poly in &polygons {
        for (a, b) in poly.iter().zip(poly.iter().cycle().skip(1)) {
            budget.visit(dec)?;
            let line = Line::through(a, b);
            let i = line.axis;
            lines
                .entry(line)
                .or_default()
                .extend([a[i].clone(), b[i].clone()]);
        }
    }
    let mut points = Vec::new();
    let mut triangles = Vec::new();
    let mut ids = BTreeMap::new();
    let mut insert = |p: Point| {
        if let Some(i) = ids.get(&p) {
            return Ok(*i);
        }
        if points.len() == budget.options.max_vertices {
            return Err(refusal(
                "field.print-budget",
                "lattice print vertex limit reached",
            ));
        }
        let i = points.len();
        points.push(p.clone());
        ids.insert(p, i);
        Ok(i)
    };
    for poly in polygons {
        let centre = insert(centroid(&poly))?;
        let mut boundary = Vec::new();
        for (a, b) in poly.iter().zip(poly.iter().cycle().skip(1)) {
            let line = Line::through(a, b);
            let i = line.axis;
            let knots = &lines[&line];
            let start = boundary.len();
            if a[i] < b[i] {
                for t in knots.range(&a[i]..&b[i]) {
                    boundary.push(insert(line.point(t))?);
                }
            } else {
                for t in knots.range(&b[i]..=&a[i]).rev().take_while(|&t| t > &b[i]) {
                    boundary.push(insert(line.point(t))?);
                }
            }
            dec.charge_work((boundary.len() - start) as u64)?;
        }
        triangles.extend(
            boundary
                .iter()
                .zip(boundary.iter().cycle().skip(1))
                .map(|(&a, &b)| [centre, a, b]),
        );
    }
    Ok((points, triangles))
}

pub(super) fn mesh(
    model: &Model<'_>,
    budget: &mut Budget<'_>,
    dec: &mut Decider,
) -> Result<Triangles, Refusal> {
    let regions = regions::leaf_regions(model, budget, dec)?;
    let solids = supports(model, &regions, budget, dec)?;
    let mut boundary = Vec::new();
    for (source, solid) in solids.iter().enumerate() {
        for (seed, plane) in solid.faces.iter().zip(&solid.planes) {
            budget.visit(dec)?;
            let bounds = box_of(seed.iter());
            if !regions[solid.leaf]
                .iter()
                .any(|region| overlaps(&bounds, region))
            {
                continue;
            }
            let mut cells = vec![seed.clone()];
            dec.charge_work(solids.len() as u64)?;
            let near: Vec<_> = solids
                .iter()
                .enumerate()
                .filter_map(|(i, s)| overlaps(&bounds, &s.bounds).then_some(i))
                .collect();
            for &i in &near {
                if i == source {
                    continue;
                }
                let mut divided = Vec::new();
                for poly in cells {
                    if overlaps(&box_of(poly.iter()), &solids[i].bounds) {
                        divided.extend(divide(poly, &solids[i], budget, dec)?);
                    } else {
                        divided.push(poly);
                    }
                }
                cells = divided;
            }
            for mut poly in cells {
                if emit(model, &solids, &near, source, &plane.normal, &mut poly, dec)? {
                    boundary.push(poly);
                }
            }
        }
    }
    triangulate(boundary, budget, dec)
}
