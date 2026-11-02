//! Closed finite planar regions bounded by line segments and circular arcs.

use super::{Bounds, Options};
use crate::{Circle, Plane, Point3, Vector3, math};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};

/// One complete edge of a closed planar trim chain. Arc endpoints follow the
/// circle's natural counter-clockwise direction, independently of loop winding.
#[derive(Clone, Debug)]
pub enum PlanarEdge {
    /// Closed finite straight segment, indexing `PlanarBoundary::vertices`.
    Segment([usize; 2]),
    /// A circular arc, or the complete circle when `ends` is absent.
    Arc {
        /// Exact circle carrier.
        circle: Circle,
        /// Ordered endpoints on the circle, indexing the boundary vertices.
        ends: Option<[usize; 2]>,
    },
}

/// Complete closed trim chains on one plane, including holes and disconnected
/// regions. The producer certifies edge incidence, closure and planar support.
/// Even-odd membership does not depend on loop traversal orientation.
#[derive(Clone, Debug)]
pub struct PlanarBoundary {
    /// Exact planar carrier.
    pub plane: Plane,
    /// Exact finite endpoints, shared across edges.
    pub vertices: Vec<Point3>,
    /// All physical boundary edges, never triangulation or construction seams.
    pub edges: Vec<PlanarEdge>,
}

fn min_candidate(
    best: &mut Option<(Real, Point3)>,
    point: &Point3,
    witness: Point3,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let distance = point.minus(&witness).norm_squared();
    if best
        .as_ref()
        .map(|(old, _)| math::sign(&(&distance - old), dec).map(|s| s == Sign::Negative))
        .transpose()?
        .unwrap_or(true)
    {
        *best = Some((distance, witness));
    }
    Ok(())
}

fn certificate(
    squared: &Real,
    witnesses: [Point3; 2],
    options: &Options,
    dec: &mut Decider,
) -> Result<Bounds, Refusal> {
    let mut bits = 128u32;
    loop {
        let enclosure = dec
            .enclose_budgeted_at(squared, bits)?
            .sqrt()
            .ok_or_else(|| {
                Refusal::new(
                    "shapes.distance-domain",
                    "negative squared planar distance enclosure",
                )
            })?;
        let lower = enclosure.lo().max(Rational::zero());
        let upper = enclosure.hi();
        if &upper - &lower <= options.width {
            return Ok(Bounds {
                lower,
                upper: Some(upper),
                witnesses: Some(witnesses),
                converged: true,
                cells: 0,
            });
        }
        bits = bits.checked_mul(2).ok_or_else(|| {
            Refusal::new(
                "numbers.refinement-budget",
                "planar distance precision overflow",
            )
        })?;
    }
}

impl PlanarBoundary {
    /// Owns the used vertices of a complete straight-edge polygonal trim. Input
    /// loops index `vertices`; their outer/hole orientation may be reversed.
    pub fn polygon(plane: Plane, vertices: &[Point3], loops: &[Vec<usize>]) -> Self {
        use std::collections::{BTreeMap, BTreeSet};
        let used: BTreeSet<_> = loops.iter().flatten().copied().collect();
        let indices: BTreeMap<_, _> = used
            .iter()
            .enumerate()
            .map(|(i, source)| (*source, i))
            .collect();
        let points = used.into_iter().map(|i| vertices[i].clone()).collect();
        let edges = loops
            .iter()
            .flat_map(|cycle| {
                (0..cycle.len()).map(|i| {
                    PlanarEdge::Segment([
                        indices[&cycle[i]],
                        indices[&cycle[(i + 1) % cycle.len()]],
                    ])
                })
            })
            .collect();
        Self {
            plane,
            vertices: points,
            edges,
        }
    }

    fn edge_contains(
        &self,
        edge: &PlanarEdge,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        match edge {
            PlanarEdge::Segment([a, b]) => {
                let a = point.minus(&self.vertices[*a]);
                let b = point.minus(&self.vertices[*b]);
                Ok(math::sign(&a.cross(&b).norm_squared(), dec)? == Sign::Zero
                    && math::sign(&a.dot(&b), dec)? != Sign::Positive)
            }
            PlanarEdge::Arc { circle, ends } => {
                let delta = point.minus(&circle.center);
                Ok(math::sign(&delta.dot(&circle.normal), dec)? == Sign::Zero
                    && math::sign(&(delta.norm_squared() - circle.radius_sq.clone()), dec)?
                        == Sign::Zero
                    && self.arc_contains(circle, *ends, point, dec)?)
            }
        }
    }

    fn arc_contains(
        &self,
        circle: &Circle,
        ends: Option<[usize; 2]>,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let Some([a, b]) = ends else {
            return Ok(true);
        };
        let a = self.vertices[a].minus(&circle.center);
        let b = self.vertices[b].minus(&circle.center);
        let p = point.minus(&circle.center);
        let side = |v: &Vector3, dec: &mut Decider| -> Result<bool, Refusal> {
            let y = math::sign(&a.cross(v).dot(&circle.normal), dec)?;
            Ok(y == Sign::Negative
                || (y == Sign::Zero && math::sign(&a.dot(v), dec)? == Sign::Negative))
        };
        let bp = side(&b, dec)?;
        let pp = side(&p, dec)?;
        Ok(if bp != pp {
            !pp
        } else {
            math::sign(&p.cross(&b).dot(&circle.normal), dec)? != Sign::Negative
        })
    }

    fn edge_closest(
        &self,
        edge: &PlanarEdge,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Point3, Refusal> {
        match edge {
            PlanarEdge::Segment([a, b]) => {
                let a = &self.vertices[*a];
                let direction = self.vertices[*b].minus(a);
                let norm = direction.norm_squared();
                if math::sign(&norm, dec)? == Sign::Zero {
                    return Ok(a.clone());
                }
                let parameter = point.minus(a).dot(&direction) / norm;
                let parameter = if math::sign(&parameter, dec)? == Sign::Negative {
                    Real::zero()
                } else if math::sign(&(&parameter - &Real::one()), dec)? == Sign::Positive {
                    Real::one()
                } else {
                    parameter
                };
                Ok(a.plus(&direction.scaled(&parameter)))
            }
            PlanarEdge::Arc { circle, ends } => {
                let delta = point.minus(&circle.center);
                let radial = delta.minus(
                    &circle
                        .normal
                        .scaled(&(delta.dot(&circle.normal) / circle.normal.norm_squared())),
                );
                let norm = radial.norm_squared();
                if math::sign(&norm, dec)? == Sign::Zero {
                    return Ok(match ends {
                        Some([a, _]) => self.vertices[*a].clone(),
                        None => {
                            circle.center.plus(&circle.ref_dir.scaled(
                                &(&circle.radius_sq / &circle.ref_dir.norm_squared()).sqrt(),
                            ))
                        }
                    });
                }
                let projected = circle
                    .center
                    .plus(&radial.scaled(&(&circle.radius_sq / &norm).sqrt()));
                if self.arc_contains(circle, *ends, &projected, dec)? {
                    return Ok(projected);
                }
                let [a, b] = ends.expect("a whole circle admits every radial projection");
                let a = &self.vertices[a];
                let b = &self.vertices[b];
                Ok(
                    if math::sign(
                        &(point.minus(a).norm_squared() - point.minus(b).norm_squared()),
                        dec,
                    )? == Sign::Positive
                    {
                        b.clone()
                    } else {
                        a.clone()
                    },
                )
            }
        }
    }

    /// Exact closed trim membership, including line/arc junctions and hole rims.
    pub fn contains(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        if math::sign(
            &point.minus(&self.plane.origin).dot(&self.plane.normal),
            dec,
        )? != Sign::Zero
        {
            return Ok(false);
        }
        for edge in &self.edges {
            if self.edge_contains(edge, point, dec)? {
                return Ok(true);
            }
        }
        let u = &self.plane.ref_dir;
        let v = self.plane.normal.cross(u);
        let mut slope = 0_i64;
        let direction = loop {
            let direction = u.plus(&v.scaled(&Real::from_i64(slope)));
            let mut separated = true;
            for vertex in &self.vertices {
                // Coplanar named regions may share unused vertices. The ray
                // origin and vertices behind it cannot be forward crossings.
                let delta = vertex.minus(point);
                if math::sign(&delta.cross(&direction).dot(&self.plane.normal), dec)? == Sign::Zero
                    && math::sign(&delta.dot(&direction), dec)? == Sign::Positive
                {
                    separated = false;
                    break;
                }
            }
            if separated {
                for edge in &self.edges {
                    if let PlanarEdge::Arc { circle, .. } = edge {
                        let delta = circle.center.minus(point);
                        if math::sign(
                            &(delta.cross(&direction).norm_squared()
                                - &circle.radius_sq * &direction.norm_squared()),
                            dec,
                        )? == Sign::Zero
                        {
                            separated = false;
                            break;
                        }
                    }
                }
            }
            if separated {
                break direction;
            }
            slope = slope.checked_add(1).ok_or_else(|| {
                Refusal::new("shapes.distance-domain", "planar trim ray index overflow")
            })?;
        };
        let mut crossings = 0usize;
        for edge in &self.edges {
            match edge {
                PlanarEdge::Segment([a, b]) => {
                    let start = &self.vertices[*a];
                    let side = self.vertices[*b].minus(start);
                    let determinant = direction.cross(&side).dot(&self.plane.normal);
                    if math::sign(&determinant, dec)? == Sign::Zero {
                        continue;
                    }
                    let delta = start.minus(point);
                    let t = delta.cross(&side).dot(&self.plane.normal) / determinant.clone();
                    let s = delta.cross(&direction).dot(&self.plane.normal) / determinant;
                    if math::sign(&t, dec)? == Sign::Positive
                        && math::sign(&s, dec)? == Sign::Positive
                        && math::sign(&(&s - &Real::one()), dec)? == Sign::Negative
                    {
                        crossings += 1;
                    }
                }
                PlanarEdge::Arc { circle, ends } => {
                    let delta = point.minus(&circle.center);
                    let a = direction.norm_squared();
                    let b = delta.dot(&direction);
                    let discriminant =
                        b.square() - &a * &(&delta.norm_squared() - &circle.radius_sq);
                    if math::sign(&discriminant, dec)? != Sign::Positive {
                        continue;
                    }
                    let root = discriminant.sqrt();
                    for sign in [-1, 1] {
                        let t = (-&b + &root * &Real::from_i64(sign)) / a.clone();
                        if math::sign(&t, dec)? != Sign::Positive {
                            continue;
                        }
                        let hit = point.plus(&direction.scaled(&t));
                        if self.arc_contains(circle, *ends, &hit, dec)? {
                            crossings += 1;
                        }
                    }
                }
            }
        }
        Ok(crossings % 2 == 1)
    }

    /// Exact minimum squared distance and one closest point on the closed face.
    pub fn closest(&self, point: &Point3, dec: &mut Decider) -> Result<(Real, Point3), Refusal> {
        let height = point.minus(&self.plane.origin).dot(&self.plane.normal);
        let projection = point.plus(
            &self
                .plane
                .normal
                .scaled(&(-height / self.plane.normal.norm_squared())),
        );
        if self.contains(&projection, dec)? {
            return Ok((point.minus(&projection).norm_squared(), projection));
        }
        let mut best = None;
        for edge in &self.edges {
            min_candidate(&mut best, point, self.edge_closest(edge, point, dec)?, dec)?;
        }
        best.ok_or_else(|| {
            Refusal::new(
                "shapes.distance-domain",
                "a finite planar face needs a complete trim boundary",
            )
        })
    }

    /// Certified finite-face distance. Analytic trim candidates, not carrier
    /// distance or a spatial grid, supply both the lower bound and its witness.
    pub fn point_distance(
        &self,
        point: &Point3,
        options: &Options,
        dec: &mut Decider,
    ) -> Result<Bounds, Refusal> {
        super::validate(options)?;
        let (squared, witness) = self.closest(point, dec)?;
        certificate(&squared, [point.clone(), witness], options, dec)
    }

    /// Complete finite polygon-to-polygon minimum distance, including crossings,
    /// holes and skew planes. Circular edge pairs keep the general certified
    /// analytic search rather than using a polygonal substitute.
    pub fn face_distance(
        &self,
        other: &Self,
        options: &Options,
        dec: &mut Decider,
    ) -> Result<Option<Bounds>, Refusal> {
        super::validate(options)?;
        if self
            .edges
            .iter()
            .chain(&other.edges)
            .any(|edge| !matches!(edge, PlanarEdge::Segment(_)))
        {
            return Ok(None);
        }
        let mut best: Option<(Real, [Point3; 2])> = None;
        fn admit(
            best: &mut Option<(Real, [Point3; 2])>,
            pair: [Point3; 2],
            dec: &mut Decider,
        ) -> Result<(), Refusal> {
            let squared = pair[0].minus(&pair[1]).norm_squared();
            if best
                .as_ref()
                .map(|(old, _)| math::sign(&(&squared - old), dec).map(|s| s == Sign::Negative))
                .transpose()?
                .unwrap_or(true)
            {
                *best = Some((squared, pair));
            }
            Ok(())
        }
        for (first, second, reverse) in [(self, other, false), (other, self, true)] {
            let mut visited = vec![false; first.vertices.len()];
            for edge in &first.edges {
                let PlanarEdge::Segment([a, b]) = edge else {
                    unreachable!()
                };
                for vertex in [*a, *b] {
                    if visited[vertex] {
                        continue;
                    }
                    visited[vertex] = true;
                    let point = &first.vertices[vertex];
                    let (_, witness) = second.closest(point, dec)?;
                    admit(
                        &mut best,
                        if reverse {
                            [witness, point.clone()]
                        } else {
                            [point.clone(), witness]
                        },
                        dec,
                    )?;
                }
                let a = &first.vertices[*a];
                let b = &first.vertices[*b];
                let direction = b.minus(a);
                let denominator = direction.dot(&second.plane.normal);
                if math::sign(&denominator, dec)? != Sign::Zero {
                    let parameter =
                        -a.minus(&second.plane.origin).dot(&second.plane.normal) / denominator;
                    if math::sign(&parameter, dec)? != Sign::Negative
                        && math::sign(&(&parameter - &Real::one()), dec)? != Sign::Positive
                    {
                        let hit = a.plus(&direction.scaled(&parameter));
                        if second.contains(&hit, dec)? {
                            return certificate(&Real::zero(), [hit.clone(), hit], options, dec)
                                .map(Some);
                        }
                    }
                }
            }
        }
        for first in &self.edges {
            let PlanarEdge::Segment([a, b]) = first else {
                unreachable!()
            };
            let a = &self.vertices[*a];
            let u = self.vertices[*b].minus(a);
            for second in &other.edges {
                let PlanarEdge::Segment([c, d]) = second else {
                    unreachable!()
                };
                let c = &other.vertices[*c];
                let v = other.vertices[*d].minus(c);
                let w = a.minus(c);
                let aa = u.norm_squared();
                let bb = u.dot(&v);
                let cc = v.norm_squared();
                let dd = u.dot(&w);
                let ee = v.dot(&w);
                let determinant = &aa * &cc - bb.square();
                if math::sign(&determinant, dec)? != Sign::Positive {
                    continue;
                }
                let s = (&bb * &ee - &cc * &dd) / determinant.clone();
                let t = (&aa * &ee - &bb * &dd) / determinant;
                if math::sign(&s, dec)? != Sign::Negative
                    && math::sign(&(&s - &Real::one()), dec)? != Sign::Positive
                    && math::sign(&t, dec)? != Sign::Negative
                    && math::sign(&(&t - &Real::one()), dec)? != Sign::Positive
                {
                    admit(
                        &mut best,
                        [a.plus(&u.scaled(&s)), c.plus(&v.scaled(&t))],
                        dec,
                    )?;
                }
            }
        }
        let (squared, pair) = best.ok_or_else(|| {
            Refusal::new(
                "shapes.distance-domain",
                "a planar polygon pair needs finite boundaries",
            )
        })?;
        certificate(&squared, pair, options, dec).map(Some)
    }
}

#[cfg(test)]
mod tests;
