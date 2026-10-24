//! Bounded convex plane arrangements with exact vertices and immutable face origins.

use std::collections::BTreeMap;

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{Curve, Line, Plane, Point3, PointDef, Surface, Vector3};
use catallax_topology::{Body, BodyId, BodyKind, Edge, Lump, Model, Sense, Shell, Vertex};

use crate::profile::sign;

/// One outward-oriented halfspace, with material on its nonpositive side.
#[derive(Clone, Debug)]
pub struct PlaneFace {
    /// Exact plane law. Its normal need not have unit length.
    pub plane: Plane,
    /// Generators of this boundary, in canonical order.
    pub origins: Vec<Origin>,
}

/// The complete vertex arrangement of a bounded intersection of halfspaces.
/// The caller must supply a bounded set, for example a subset of an input solid.
#[derive(Clone, Debug)]
pub struct ConvexCell {
    /// Defining halfspaces, including those that become redundant.
    pub faces: Vec<PlaneFace>,
    points: Vec<Point3>,
    pub(crate) polygons: Vec<Vec<usize>>,
}

fn same(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(&a.minus(b).norm_squared(), dec)? == Sign::Zero)
}
fn coordinates(p: &Point3) -> Vector3 {
    Vector3::new(p.x.clone(), p.y.clone(), p.z.clone())
}
fn intersection(a: &Plane, b: &Plane, c: &Plane, dec: &mut Decider) -> Result<Option<Point3>, Refusal> {
    let bc = b.normal.cross(&c.normal);
    let determinant = a.normal.dot(&bc);
    if sign(&determinant, dec)? == Sign::Zero {
        return Ok(None);
    }
    let h = |p: &Plane| coordinates(&p.origin).dot(&p.normal);
    let v = bc.scaled(&h(a))
        .plus(&c.normal.cross(&a.normal).scaled(&h(b)))
        .plus(&a.normal.cross(&b.normal).scaled(&h(c)))
        .scaled(&(&Real::one() / &determinant));
    Ok(Some(Point3::new(v.x, v.y, v.z)))
}

impl ConvexCell {
    /// Enumerates all feasible triple meetings, without a numerical acceptance band.
    pub fn new(faces: Vec<PlaneFace>, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut unique: Vec<PlaneFace> = Vec::with_capacity(faces.len());
        for mut face in faces {
            let mut duplicate = None;
            for (i, other) in unique.iter().enumerate() {
                if sign(&face.plane.normal.cross(&other.plane.normal).norm_squared(), dec)? == Sign::Zero
                    && sign(&face.plane.normal.dot(&other.plane.normal), dec)? == Sign::Positive
                    && sign(&face.plane.origin.minus(&other.plane.origin).dot(&other.plane.normal), dec)? == Sign::Zero {
                    duplicate = Some(i);
                    break;
                }
            }
            if let Some(i) = duplicate {
                unique[i].origins.append(&mut face.origins);
                unique[i].origins.sort();
                unique[i].origins.dedup();
            } else {
                face.origins.sort();
                face.origins.dedup();
                unique.push(face);
            }
        }
        let faces = unique;
        let mut points = Vec::new();
        for i in 0..faces.len() {
            dec.check_deadline()?;
            for j in i + 1..faces.len() {
                for k in j + 1..faces.len() {
                    let Some(p) = intersection(&faces[i].plane, &faces[j].plane, &faces[k].plane, dec)? else {
                        continue;
                    };
                    let mut inside = true;
                    for face in &faces {
                        if sign(&p.minus(&face.plane.origin).dot(&face.plane.normal), dec)? == Sign::Positive {
                            inside = false;
                            break;
                        }
                    }
                    if !inside {
                        continue;
                    }
                    let mut duplicate = false;
                    for q in &points {
                        if same(&p, q, dec)? {
                            duplicate = true;
                            break;
                        }
                    }
                    if !duplicate {
                        points.push(p);
                    }
                }
            }
        }
        // A bounded cell of dimension below three has no regularised material.
        let mut full_dimension = false;
        if let Some(a) = points.first() {
            'dimension: for b in &points[1..] {
                for c in &points[1..] {
                    let cross = b.minus(a).cross(&c.minus(a));
                    for d in &points[1..] {
                        if sign(&cross.dot(&d.minus(a)), dec)? != Sign::Zero {
                            full_dimension = true;
                            break 'dimension;
                        }
                    }
                }
            }
        }
        if !full_dimension {
            return Ok(Self { faces, points: Vec::new(), polygons: Vec::new() });
        }
        let mut polygons = Vec::with_capacity(faces.len());
        for face in &faces {
            let mut ids = Vec::new();
            for (i, p) in points.iter().enumerate() {
                if sign(&p.minus(&face.plane.origin).dot(&face.plane.normal), dec)? == Sign::Zero {
                    ids.push(i);
                }
            }
            if ids.len() < 3 {
                polygons.push(Vec::new());
                continue;
            }
            let mut sum = Vector3::zero();
            for i in &ids {
                sum = sum.plus(&coordinates(&points[*i]));
            }
            let centre = sum.scaled(&(&Real::one() / &Real::from_i64(ids.len() as i64)));
            let centre = Point3::new(centre.x, centre.y, centre.z);
            let u = &face.plane.ref_dir;
            let v = face.plane.normal.cross(u);
            let half = |p: &Point3, dec: &mut Decider| -> Result<bool, Refusal> {
                let q = p.minus(&centre);
                let y = sign(&q.dot(&v), dec)?;
                Ok(y == Sign::Positive || (y == Sign::Zero && sign(&q.dot(u), dec)? != Sign::Negative))
            };
            // Fallible insertion ordering keeps every angular comparison exact.
            for i in 1..ids.len() {
                let mut j = i;
                while j > 0 {
                    let a = &points[ids[j - 1]];
                    let b = &points[ids[j]];
                    let ha = half(a, dec)?;
                    let hb = half(b, dec)?;
                    let before = if ha != hb {
                        hb
                    } else {
                        sign(&a.minus(&centre).cross(&b.minus(&centre)).dot(&face.plane.normal), dec)? == Sign::Negative
                    };
                    if !before { break; }
                    ids.swap(j - 1, j);
                    j -= 1;
                }
            }
            polygons.push(ids);
        }
        Ok(Self { faces, points, polygons })
    }

    /// Whether regularisation removes the entire cell.
    pub fn is_empty(&self) -> bool { self.points.is_empty() }

    /// Exact vertices, used for containment and projection proofs.
    pub fn vertices(&self) -> &[Point3] { &self.points }

    /// Builds one shared edge per arrangement segment. Redundant planes vanish.
    pub fn build(&self, model: &mut Model) -> BodyId {
        if self.is_empty() { return catallax_topology::empty(model); }
        let vertices: Vec<_> = self.points.iter().map(|p| {
            let point = model.geom_mut().add_point(PointDef::Explicit(p.clone()));
            model.add_vertex(Vertex { point })
        }).collect();
        let mut edges = BTreeMap::new();
        let mut faces = Vec::new();
        for (face, polygon) in self.faces.iter().zip(&self.polygons) {
            if polygon.is_empty() { continue; }
            let mut path = Vec::with_capacity(polygon.len());
            for i in 0..polygon.len() {
                let a = polygon[i];
                let b = polygon[(i + 1) % polygon.len()];
                let key = (a.min(b), a.max(b));
                let edge = *edges.entry(key).or_insert_with(|| {
                    let curve = model.geom_mut().add_curve(Curve::Line(Line {
                        origin: self.points[key.0].clone(),
                        dir: self.points[key.1].minus(&self.points[key.0]),
                    }));
                    model.add_edge(Edge { curve, sense: Sense::Same,
                        start: Some(vertices[key.0]), end: Some(vertices[key.1]) })
                });
                path.push((edge, if a < b { Sense::Same } else { Sense::Reversed }));
            }
            faces.push(crate::extrude::face(model, Surface::Plane(face.plane.clone()), Sense::Same,
                vec![path], None, face.origins.clone()));
        }
        let shell = model.add_shell(Shell { faces });
        let lump = model.add_lump(Lump { shells: vec![shell] });
        model.add_body(Body { kind: BodyKind::Solid, lumps: vec![lump], marks: Default::default() })
    }
}
