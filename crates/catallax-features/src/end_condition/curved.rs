//! Certified single-branch curved caps over rectangular sketch domains.
//!
//! This reduction keeps analytic supports exact. It proves the same positive
//! quadratic hit branch over the complete closed rectangle, including every
//! edge and corner. Silhouettes, changing start states, projected trim events
//! and nonquadratic targets remain named coverage gaps. No sampled witness
//! licenses a branch or a missing-hit refusal.

use catallax_numbers::{Decider, Outcome, Rational, Real, Refusal, Sign, StrictSign};
use catallax_recipe::{ElementId, FaceRef, Operation, Recipe, Segment, StepId};
use catallax_shapes::{Plane, Point3, Surface, Vector3};
use std::collections::BTreeSet;

#[path = "curved/build.rs"]
mod build;
#[path = "curved/projected.rs"]
mod projected;
#[path = "curved/scene.rs"]
mod scene;

pub use build::{CapConstruction, ConvexCell};
pub use projected::Frame;
pub use scene::Scene;

/// Which positive quadratic boundary terminates the sweep.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Branch {
    /// Start outside material and stop at the first entry.
    Entry,
    /// Start inside material and stop at the first exit.
    Exit,
}

/// An exact curved first-hit definition with a whole-domain branch certificate.
#[derive(Clone, Debug)]
pub struct CurvedCapped {
    scene: Scene,
    frame: Frame,
    bounds: [[Rational; 2]; 2],
    branch: Branch,
    midpoint: Plane,
    extent: Rational,
}

/// A compiled authored rectangle, retaining its segment ids and finite bound.
#[derive(Clone, Debug)]
pub struct AuthoredCap {
    definition: CurvedCapped,
    side_ids: [ElementId; 4],
    radius: Rational,
}

impl AuthoredCap {
    /// The immutable whole-domain first-hit definition.
    pub fn definition(&self) -> &CurvedCapped {
        &self.definition
    }

    /// Authored side ids in u-min, u-max, v-min, v-max order.
    pub fn side_ids(&self) -> &[ElementId; 4] {
        &self.side_ids
    }

    /// Certified result sphere about world zero.
    pub fn radius(&self) -> &Rational {
        &self.radius
    }

    /// Compiles a recipe's rectangular sketch and stock quadric target. The
    /// same certified crossing is the first hit for both face and next modes.
    pub fn compile(
        recipe: &Recipe,
        target: &StepId,
        profile: &StepId,
        direction: Option<&catallax_recipe::Vector3>,
        limit: Option<&FaceRef>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        recipe
            .check()
            .map_err(|e| Refusal::new(e.code(), e.to_string()))?;
        dec.check_deadline()?;
        let source = Scene::compile(recipe, target, dec)?;
        if let Some(reference) = limit {
            source.resolve_wall(reference)?;
        }
        let Some(Operation::Sketch { plane, loops }) = recipe.step(profile).map(|s| &s.op) else {
            return Err(Refusal::new(
                "recipe.input-failed",
                "curved extrusion input is not a sketch",
            ));
        };
        if loops.len() != 1 || loops[0].0.len() != 4 {
            return Err(coverage(
                "curved cap profile currently requires one four-line rectangle",
            ));
        }
        let mut points = Vec::new();
        let mut ids = Vec::new();
        for segment in &loops[0].0 {
            let Segment::Line { id, to } = segment else {
                return Err(coverage("curved cap profile has a non-line segment"));
            };
            ids.push(id.clone());
            points.push(to.0.each_ref().map(|n| n.value().clone()));
        }
        let bounds: [[Rational; 2]; 2] = std::array::from_fn(|i| {
            [
                points
                    .iter()
                    .map(|p| p[i].clone())
                    .min()
                    .expect("four endpoints"),
                points
                    .iter()
                    .map(|p| p[i].clone())
                    .max()
                    .expect("four endpoints"),
            ]
        });
        let mut seen = BTreeSet::new();
        let mut sides: [Option<ElementId>; 4] = std::array::from_fn(|_| None);
        for i in 0..4 {
            let a = &points[(i + 3) % 4];
            let b = &points[i];
            if !seen.insert(b.clone())
                || (0..2).any(|j| b[j] != bounds[j][0] && b[j] != bounds[j][1])
            {
                return Err(coverage("the sketch is not a nondegenerate rectangle"));
            }
            let wall = if a[0] == b[0] && a[1] != b[1] {
                usize::from(b[0] == bounds[0][1])
            } else if a[1] == b[1] && a[0] != b[0] {
                2 + usize::from(b[1] == bounds[1][1])
            } else {
                return Err(coverage("the sketch rectangle has a diagonal or zero edge"));
            };
            if sides[wall].replace(ids[i].clone()).is_some() {
                return Err(coverage("the sketch repeats a rectangle side"));
            }
        }
        let pose = scene::placement(plane)?;
        let direction = direction
            .map(|v| {
                let [x, y, z] = v.0.each_ref().map(|n| Real::from(n.value().clone()));
                Vector3::new(x, y, z)
            })
            .unwrap_or_else(|| pose.apply_vector(&axis(2)));
        let frame = Frame {
            origin: pose.apply_point(&Point3::origin()),
            u: pose.apply_vector(&axis(0)),
            v: pose.apply_vector(&axis(1)),
            direction,
        };
        let radius = source
            .radius
            .clone()
            .max(build::frame_radius(&frame, &bounds, dec)?);
        let definition = CurvedCapped::rectangle(source, frame, bounds, dec)?;
        Ok(Self {
            definition,
            side_ids: sides.map(|s| s.expect("each rectangle side occurs once")),
            radius,
        })
    }

    /// Builds the complete analytic result with the authored segment origins.
    pub fn build(
        &self,
        model: &mut catallax_topology::Model,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<catallax_topology::BodyId, Refusal> {
        self.definition
            .build(model, step, self.side_ids.clone(), dec)
    }
}

/// Exact classification relative to the closed capped rectangle.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Location {
    /// Three-dimensional interior.
    Inside,
    /// Outside the closed set.
    Outside,
    /// On the actual result boundary.
    Boundary,
}

fn coverage(detail: &str) -> Refusal {
    Refusal::new("features.unsupported-curved-limit", detail)
}

fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    // Even rational shortcuts must honour the whole-request allowance.
    dec.check_deadline()?;
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::Refused(r) => Err(r),
        Outcome::NotExact(_) => Err(coverage(
            "the exact curved-limit predicate was not certified",
        )),
    }
}

fn axis(i: usize) -> Vector3 {
    let [x, y, z] = std::array::from_fn(|j| Real::from_i64(i64::from(i == j)));
    Vector3::new(x, y, z)
}

impl CurvedCapped {
    /// Certifies a capped rectangle. The bounds are the profile definition,
    /// rather than a sample-derived approximation to a different domain.
    pub fn rectangle(
        scene: Scene,
        frame: Frame,
        bounds: [[Rational; 2]; 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let certificate = projected::certify(&scene, &frame, &bounds, dec)?;
        Ok(Self {
            scene,
            frame,
            bounds,
            branch: certificate.branch,
            midpoint: certificate.midpoint,
            extent: certificate.extent,
        })
    }

    /// The branch proved constant over the entire profile.
    pub fn branch(&self) -> Branch {
        self.branch
    }

    /// The exact source definition; no B-rep is needed for queries.
    pub fn scene(&self) -> &Scene {
        &self.scene
    }

    /// The immutable authored profile frame and raw ray direction.
    pub fn frame(&self) -> &Frame {
        &self.frame
    }

    /// The closed rectangular profile domain.
    pub fn bounds(&self) -> &[[Rational; 2]; 2] {
        &self.bounds
    }

    /// A strict rational upper bound on every selected line parameter.
    pub fn extent(&self) -> &Rational {
        &self.extent
    }

    /// The plane separating the entry and exit roots; its negative halfspace
    /// contains the entry branch. It is construction geometry, not a real cap.
    pub fn midpoint(&self) -> &Plane {
        &self.midpoint
    }

    /// The exact positive first hit at a point of the certified rectangle.
    pub fn height(&self, q: &[Real; 2], dec: &mut Decider) -> Result<Real, Refusal> {
        for (x, [lo, hi]) in q.iter().zip(&self.bounds) {
            if sign(&(x - &Real::from(lo.clone())), dec)? == Sign::Negative
                || sign(&(x - &Real::from(hi.clone())), dec)? == Sign::Positive
            {
                return Err(coverage("height query is outside its certified domain"));
            }
        }
        let [c, b, a] =
            projected::coefficients(&self.scene.surface, &self.frame).map(|p| p.eval_real(q));
        let d = &b * &b - (&a * &c) * Real::from_i64(4);
        // Construction already proves d>0 and a>0 over the complete domain.
        let root = d.sqrt();
        let signed = if self.branch == Branch::Entry {
            -root
        } else {
            root
        };
        Ok((-b + signed) / (Real::from_i64(2) * a))
    }

    /// Classifies profile/ray coordinates without constructing topology.
    /// The field caller transforms world points through the retained frame.
    pub fn location(
        &self,
        q: &[Real; 2],
        t: &Real,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        let mut boundary = false;
        for (x, [lo, hi]) in q.iter().zip(&self.bounds) {
            let a = sign(&(x - &Real::from(lo.clone())), dec)?;
            let b = sign(&(x - &Real::from(hi.clone())), dec)?;
            if a == Sign::Negative || b == Sign::Positive {
                return Ok(Location::Outside);
            }
            boundary |= a == Sign::Zero || b == Sign::Zero;
        }
        let start = sign(t, dec)?;
        if start == Sign::Negative {
            return Ok(Location::Outside);
        }
        let height = self.height(q, dec)?;
        let end = sign(&(t - &height), dec)?;
        if end == Sign::Positive {
            return Ok(Location::Outside);
        }
        Ok(if boundary || start == Sign::Zero || end == Sign::Zero {
            Location::Boundary
        } else {
            Location::Inside
        })
    }
}
