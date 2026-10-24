//! Exact reductions of distance features on a proved convex planar skin.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{EdgeName, Origin, Role, StepId, ThickenSide};
use catallax_shapes::{Plane, Point3, Surface};

use crate::{polyhedron::{ConvexCell, PlaneFace}, profile::sign};

fn bad(detail: &str) -> Refusal { Refusal::new("distance.planar-domain", detail) }
fn shifted(face: &PlaneFace, d: &Real, dec: &mut Decider) -> Result<PlaneFace, Refusal> {
    let Some(Surface::Plane(plane)) = Surface::Plane(face.plane.clone()).parallel_support(d, dec)? else {
        unreachable!("a plane parallel never collapses")
    };
    Ok(PlaneFace { plane, origins: face.origins.clone() })
}

impl ConvexCell {
    /// Erodes a convex polyhedron. Halfspace clearances are exact distances to its exterior.
    pub fn eroded(&self, distance: &Real, dec: &mut Decider) -> Result<Self, Refusal> {
        let faces = self.faces.iter().map(|f| shifted(f, &(-distance), dec)).collect::<Result<_, _>>()?;
        Self::new(faces, dec)
    }

    /// Constructs a cavity when a minimum kept-plane clearance has a trimmed witness.
    /// A minimum-clearance projection cannot cross another kept halfspace: that
    /// neighbour would have smaller clearance. Only open halfspaces need testing.
    pub fn cavity(&self, open: &[bool], thickness: &Real, step: &StepId, dec: &mut Decider) -> Result<Self, Refusal> {
        if open.len() != self.faces.len() { return Err(bad("selection and face inventory differ")); }
        let witnesses = if open.iter().any(|v| *v) { self.vertices() } else { &[] };
        let mut faces = self.faces.clone();
        for (i, face) in self.faces.iter().enumerate() {
            if open[i] { continue; }
            let plane = &face.plane;
            for p in witnesses {
                let projection = p.plus(&plane.normal.scaled(
                    &(-p.minus(&plane.origin).dot(&plane.normal) / plane.normal.norm_squared())));
                for (j, trim) in self.faces.iter().enumerate() {
                    if !open[j] { continue; }
                    if sign(&projection.minus(&trim.plane.origin).dot(&trim.plane.normal), dec)? == Sign::Positive {
                        return Err(bad("a kept-face projection leaves its closed trim; edge and vertex distance envelopes are required"));
                    }
                }
            }
            faces[i] = shifted(face, &(-thickness), dec)?;
            faces[i].origins = face.origins.iter().map(|of| Origin {
                step: step.clone(), role: Role::Inner { of: Box::new(of.clone()) },
            }).collect();
        }
        Self::new(faces, dec)
    }

    /// A selected planar face's complete normal sweep, including its straight rims.
    pub fn face_sweep(&self, selected: usize, thickness: &Real, side: ThickenSide, step: &StepId, dec: &mut Decider) -> Result<Self, Refusal> {
        let face = &self.faces[selected];
        let base = &face.plane;
        let positive = matches!(side, ThickenSide::Front);
        let distance = if positive { thickness.clone() } else { -thickness };
        let mut far = shifted(face, &distance, dec)?;
        far.origins = face.origins.iter().map(|of| Origin {
            step: step.clone(), role: Role::Offset { of: Box::new(of.clone()) },
        }).collect();
        let mut start = face.clone();
        if positive { start.plane.normal = start.plane.normal.scaled(&Real::from_i64(-1)); }
        else { far.plane.normal = far.plane.normal.scaled(&Real::from_i64(-1)); }
        let mut faces = vec![start, far];
        for (i, neighbour) in self.faces.iter().enumerate() {
            if i == selected { continue; }
            let normal = neighbour.plane.normal.minus(&base.normal.scaled(
                &(neighbour.plane.normal.dot(&base.normal) / base.normal.norm_squared())));
            let length = normal.norm_squared();
            if sign(&length, dec)? == Sign::Zero { continue; }
            let offset = neighbour.plane.origin.minus(&base.origin).dot(&neighbour.plane.normal);
            let origin = base.origin.plus(&normal.scaled(&(offset / length)));
            let origins = match (face.origins.first(), neighbour.origins.first()) {
                (Some(a), Some(b)) => vec![Origin {
                    step: step.clone(), role: Role::Rim { edge: EdgeName::new(a.clone(), b.clone()).map_err(|e| bad(&e.to_string()))? },
                }],
                _ => return Err(bad("a normal-sweep rim requires named input faces")),
            };
            faces.push(PlaneFace { plane: Plane { origin, normal, ref_dir: base.normal.clone() }, origins });
        }
        Self::new(faces, dec)
    }

    /// Proves a complete input vertex inventory lies in every oriented halfspace.
    pub fn contains_vertices(&self, vertices: &[Point3], dec: &mut Decider) -> Result<bool, Refusal> {
        for p in vertices {
            for f in &self.faces {
                if sign(&p.minus(&f.plane.origin).dot(&f.plane.normal), dec)? == Sign::Positive {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    /// Resolves final trimmed faces, including adjacency and unique-nearest filters.
    pub fn select(&self, references: &[catallax_recipe::FaceRef], dec: &mut Decider) -> Result<Vec<bool>, Refusal> {
        let mut selected = vec![false; self.faces.len()];
        for reference in references {
            let mut candidates = Vec::new();
            for (i, face) in self.faces.iter().enumerate() {
                if self.polygons[i].is_empty() || !face.origins.iter().any(|o| reference.face.matches(o)) {
                    continue;
                }
                let adjacent = reference.touching.iter().all(|pattern| {
                    self.faces.iter().enumerate().any(|(j, other)| {
                        j != i && other.origins.iter().any(|o| pattern.matches(o))
                            && self.polygons[i].iter().filter(|v| self.polygons[j].contains(v)).count() >= 2
                    })
                });
                if adjacent { candidates.push(i); }
            }
            if candidates.is_empty() { return Err(Refusal::new("reference.lost-face", format!("no input face resolves {reference:?}"))); }
            if let Some(near) = &reference.near {
                let [x, y, z] = near.0.each_ref().map(|n| n.value().clone().into());
                let point = Point3::new(x, y, z);
                let mut best = candidates[0];
                let mut distance = self.face_distance_squared(best, &point, dec)?;
                let mut tie = false;
                for candidate in &candidates[1..] {
                    let d = self.face_distance_squared(*candidate, &point, dec)?;
                    match sign(&(&d - &distance), dec)? {
                        Sign::Negative => { best = *candidate; distance = d; tie = false; }
                        Sign::Zero => tie = true,
                        Sign::Positive => (),
                    }
                }
                if tie { return Err(Refusal::new("reference.ambiguous-face", "nearest trimmed faces are tied")); }
                candidates = vec![best];
            }
            for i in candidates { selected[i] = true; }
        }
        Ok(selected)
    }

    fn face_distance_squared(&self, face: usize, point: &Point3, dec: &mut Decider) -> Result<Real, Refusal> {
        let plane = &self.faces[face].plane;
        let height = point.minus(&plane.origin).dot(&plane.normal);
        let projection = point.plus(&plane.normal.scaled(&(-&height / plane.normal.norm_squared())));
        let mut inside = true;
        for face in &self.faces {
            if sign(&projection.minus(&face.plane.origin).dot(&face.plane.normal), dec)? == Sign::Positive {
                inside = false;
                break;
            }
        }
        if inside { return Ok(height.square() / plane.normal.norm_squared()); }
        let polygon = &self.polygons[face];
        let mut best: Option<Real> = None;
        for i in 0..polygon.len() {
            let a = &self.vertices()[polygon[i]];
            let b = &self.vertices()[polygon[(i + 1) % polygon.len()]];
            let direction = b.minus(a);
            let t = point.minus(a).dot(&direction) / direction.norm_squared();
            let t = if sign(&t, dec)? == Sign::Negative { Real::zero() }
                else if sign(&(&t - &Real::one()), dec)? == Sign::Positive { Real::one() } else { t };
            let d = point.minus(&a.plus(&direction.scaled(&t))).norm_squared();
            if best.as_ref().map(|old| sign(&(&d - old), dec).map(|s| s == Sign::Negative)).transpose()?.unwrap_or(true) {
                best = Some(d);
            }
        }
        best.ok_or_else(|| bad("a selected face has no regular polygon"))
    }
}
