//! Exact convex profile halfplanes and their finite swept column.

use super::*;
use catallax_recipe::Origin;
use catallax_shapes::{Curve, Line, PointDef};
use catallax_topology::{Body, BodyId, BodyKind, Edge, Lump, Model, Sense, Shell, Vertex};
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
pub(super) struct Polygon {
    // Each authored segment runs from the previous endpoint to this endpoint.
    points: Vec<[Rational; 2]>,
    sides: Vec<[Rational; 3]>,
    winding: Sign,
}

impl Polygon {
    pub(super) fn new(points: Vec<[Rational; 2]>, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut seen = BTreeSet::new();
        let mut area = Rational::zero();
        for (i, b) in points.iter().enumerate() {
            dec.check_deadline()?;
            if !seen.insert(b.clone()) {
                return Err(coverage("convex cap profile repeats a vertex"));
            }
            let a = &points[(i + points.len() - 1) % points.len()];
            area = area + &a[0] * &b[1] - &a[1] * &b[0];
        }
        let winding = area.sign();
        if winding == Sign::Zero {
            return Err(coverage(
                "convex cap profile has no two-dimensional interior",
            ));
        }
        let mut sides = Vec::with_capacity(points.len());
        for (i, b) in points.iter().enumerate() {
            let a = &points[(i + points.len() - 1) % points.len()];
            let [du, dv] = std::array::from_fn(|j| &b[j] - &a[j]);
            // For CCW winding, the right side of each directed edge is outside.
            let mut side = [dv.clone(), -&du, &du * &a[1] - &dv * &a[0]];
            if winding == Sign::Negative {
                side = side.map(|c| -c);
            }
            for p in &points {
                dec.check_deadline()?;
                if (&side[0] * &p[0] + &side[1] * &p[1] + side[2].clone()).sign() == Sign::Positive
                {
                    return Err(coverage("cap profile is not a convex line polygon"));
                }
            }
            sides.push(side);
        }
        Ok(Self {
            points,
            sides,
            winding,
        })
    }

    pub(super) fn sides(&self) -> &[[Rational; 3]] {
        &self.sides
    }

    fn forward(&self, frame: &Frame, dec: &mut Decider) -> Result<bool, Refusal> {
        let orientation = sign(&frame.u.cross(&frame.v).dot(&frame.direction), dec)?;
        if orientation == Sign::Zero {
            return Err(coverage("construction polygon column is degenerate"));
        }
        Ok(orientation == self.winding)
    }

    /// Outward supports in authored side order, then start and artificial end.
    pub(super) fn surfaces(
        &self,
        frame: &Frame,
        extent: &Rational,
        dec: &mut Decider,
    ) -> Result<Vec<Surface>, Refusal> {
        let forward = self.forward(frame, dec)?;
        let offset = frame.direction.scaled(&Real::from(extent.clone()));
        let mut surfaces = Vec::with_capacity(self.points.len() + 2);
        for (i, b) in self.points.iter().enumerate() {
            let a = &self.points[(i + self.points.len() - 1) % self.points.len()];
            let a = frame.point(&a.clone().map(Real::from));
            let edge = frame.point(&b.clone().map(Real::from)).minus(&a);
            let normal = edge.cross(&offset);
            surfaces.push(Surface::Plane(Plane::new(
                a,
                if forward {
                    normal
                } else {
                    normal.scaled(&Real::from_i64(-1))
                },
                edge,
                dec,
            )?));
        }
        let mut normal = frame.u.cross(&frame.v);
        if sign(&normal.dot(&frame.direction), dec)? == Sign::Negative {
            normal = normal.scaled(&Real::from_i64(-1));
        }
        for end in [false, true] {
            surfaces.push(Surface::Plane(Plane::new(
                if end {
                    frame.origin.plus(&offset)
                } else {
                    frame.origin.clone()
                },
                if end {
                    normal.clone()
                } else {
                    normal.scaled(&Real::from_i64(-1))
                },
                frame.u.clone(),
                dec,
            )?));
        }
        Ok(surfaces)
    }

    pub(super) fn build(
        &self,
        model: &mut Model,
        frame: &Frame,
        extent: &Rational,
        names: &[Vec<Origin>],
        dec: &mut Decider,
    ) -> Result<BodyId, Refusal> {
        // Maximal straight sides avoid coplanar face seams in the Boolean
        // input. Their exact support and every authored contribution survive.
        // Analytic cells may retain the equivalent duplicate halfplanes.
        let authored_n = self.points.len();
        let mut corners = Vec::new();
        for (i, b) in self.points.iter().enumerate() {
            dec.check_deadline()?;
            let a = &self.points[(i + authored_n - 1) % authored_n];
            let c = &self.points[(i + 1) % authored_n];
            let ab: [Rational; 2] = std::array::from_fn(|j| &b[j] - &a[j]);
            let bc: [Rational; 2] = std::array::from_fn(|j| &c[j] - &b[j]);
            if &ab[0] * &bc[1] != &ab[1] * &bc[0] {
                corners.push(i);
            }
        }
        let n = corners.len();
        let authored_surfaces = self.surfaces(frame, extent, dec)?;
        let mut surfaces = Vec::with_capacity(n + 2);
        let mut merged_names = Vec::with_capacity(n + 2);
        for (j, &i) in corners.iter().enumerate() {
            let mut k = (corners[(j + n - 1) % n] + 1) % authored_n;
            let mut contributions = Vec::new();
            loop {
                contributions.extend(names[k].iter().cloned());
                if k == i {
                    break;
                }
                k = (k + 1) % authored_n;
            }
            contributions.sort();
            contributions.dedup();
            merged_names.push(contributions);
            surfaces.push(authored_surfaces[i].clone());
        }
        for i in authored_n..authored_n + 2 {
            surfaces.push(authored_surfaces[i].clone());
            merged_names.push(names[i].clone());
        }
        let forward = self.forward(frame, dec)?;
        let offset = frame.direction.scaled(&Real::from(extent.clone()));
        let mut points: Vec<_> = self
            .points
            .iter()
            .enumerate()
            .filter(|(i, _)| corners.contains(i))
            .map(|(_, p)| frame.point(&p.clone().map(Real::from)))
            .collect();
        points.extend((0..n).map(|i| points[i].plus(&offset)).collect::<Vec<_>>());
        let vertices: Vec<_> = points
            .iter()
            .map(|p| {
                let point = model.geom_mut().add_point(PointDef::Explicit(p.clone()));
                model.add_vertex(Vertex { point })
            })
            .collect();
        let mut cycles = Vec::with_capacity(n + 2);
        for i in 0..n {
            let a = (i + n - 1) % n;
            let mut cycle = vec![a, i, i + n, a + n];
            if !forward {
                cycle.reverse();
            }
            cycles.push(cycle);
        }
        let mut start: Vec<_> = (0..n).collect();
        let mut end: Vec<_> = (n..2 * n).collect();
        if forward {
            start.reverse();
        } else {
            end.reverse();
        }
        cycles.extend([start, end]);
        let mut edges = BTreeMap::new();
        let mut faces = Vec::with_capacity(n + 2);
        for ((cycle, surface), names) in cycles.into_iter().zip(surfaces).zip(merged_names) {
            dec.check_deadline()?;
            let mut path = Vec::with_capacity(cycle.len());
            for (j, &a) in cycle.iter().enumerate() {
                let b = cycle[(j + 1) % cycle.len()];
                let key = (a.min(b), a.max(b));
                let edge = *edges.entry(key).or_insert_with(|| {
                    let curve = model.geom_mut().add_curve(Curve::Line(Line {
                        origin: points[key.0].clone(),
                        dir: points[key.1].minus(&points[key.0]),
                    }));
                    model.add_edge(Edge {
                        curve,
                        sense: Sense::Same,
                        start: Some(vertices[key.0]),
                        end: Some(vertices[key.1]),
                    })
                });
                path.push((edge, if a < b { Sense::Same } else { Sense::Reversed }));
            }
            // The shared sweep helper records exact vertex/edge/support incidences.
            faces.push(crate::extrude::face(
                model,
                surface,
                Sense::Same,
                vec![path],
                None,
                names.clone(),
            ));
        }
        let shell = model.add_shell(Shell { faces });
        let lump = model.add_lump(Lump {
            shells: vec![shell],
        });
        Ok(model.add_body(Body {
            kind: BodyKind::Solid,
            lumps: vec![lump],
            marks: Default::default(),
        }))
    }
}
