//! Certified single-branch curved caps over convex polygons and coaxial disks.
//!
//! This reduction keeps analytic supports exact. It proves the same positive
//! quadratic hit branch over a closed rectangular enclosure, including every
//! edge and corner. Exact polygon sides define the swept domain. Coaxial
//! circular profiles on spheres instead use a whole-disk radial certificate
//! and an exact cylindrical wall. Silhouettes,
//! changing start states, projected trim events
//! and nonquadratic targets remain named coverage gaps. No sampled witness
//! licenses a branch or a missing-hit refusal.

use catallax_numbers::{Decider, Outcome, Rational, Real, Refusal, Sign, StrictSign};
use catallax_recipe::{ElementId, FaceRef, Operation, Recipe, Segment, StepId};
use catallax_shapes::{Plane, Point3, Surface, Vector3};
use std::collections::BTreeSet;

#[path = "curved/build.rs"]
mod build;
#[path = "curved/circular.rs"]
mod circular;
#[path = "curved/circular_build.rs"]
mod circular_build;
#[path = "curved/polygon.rs"]
mod polygon;
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
    sides: Vec<[Rational; 3]>,
    polygon: Option<polygon::Polygon>,
    circle: Option<circular::Circle>,
    branch: Branch,
    midpoint: Plane,
    extent: Rational,
}

/// A compiled authored convex polygon or coaxial disk, with segment ids and bound.
#[derive(Clone, Debug)]
pub struct AuthoredCap {
    definition: CurvedCapped,
    side_ids: Vec<ElementId>,
    radius: Rational,
}

impl AuthoredCap {
    /// The immutable whole-domain first-hit definition.
    pub fn definition(&self) -> &CurvedCapped {
        &self.definition
    }

    /// Authored side ids in the definition's side-constraint order. Axis-aligned
    /// rectangles retain their u-min, u-max, v-min, v-max order; other polygons
    /// retain authored segment order, in either winding. A disk has one circle id.
    pub fn side_ids(&self) -> &[ElementId] {
        &self.side_ids
    }

    /// Certified result sphere about world zero.
    pub fn radius(&self) -> &Rational {
        &self.radius
    }

    /// Compiles a convex polygon against a stock quadric, or a coaxial disk
    /// against a stock sphere. The same certified crossing is the first hit
    /// for both face and next modes.
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
        if loops.len() == 1
            && let [Segment::Circle { id, centre, radius }] = loops[0].0.as_slice()
        {
            return circular::compile(source, plane, id, centre, radius, direction, dec);
        }
        if loops.len() != 1 || loops[0].0.len() < 3 {
            return Err(coverage(
                "curved cap profile requires one convex rational line polygon",
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
                    .expect("at least three endpoints"),
                points
                    .iter()
                    .map(|p| p[i].clone())
                    .max()
                    .expect("at least three endpoints"),
            ]
        });
        // Keep the established rectangle construction and physical side order.
        let rectangle_ids = rectangle_ids(&points, &ids, &bounds);
        let polygon = if rectangle_ids.is_none() {
            Some(polygon::Polygon::new(points, dec)?)
        } else {
            None
        };
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
        let mut definition = CurvedCapped::rectangle(source, frame, bounds, dec)?;
        if let Some(polygon) = polygon {
            definition.sides = polygon.sides().to_vec();
            definition.polygon = Some(polygon);
        }
        Ok(Self {
            definition,
            side_ids: rectangle_ids.map(Vec::from).unwrap_or(ids),
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
        self.definition.build(model, step, &self.side_ids, dec)
    }
}

/// Exact classification relative to the closed capped profile.
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

fn rectangle_ids(
    points: &[[Rational; 2]],
    ids: &[ElementId],
    bounds: &[[Rational; 2]; 2],
) -> Option<[ElementId; 4]> {
    if points.len() != 4 {
        return None;
    }
    let mut seen = BTreeSet::new();
    let mut sides: [Option<ElementId>; 4] = std::array::from_fn(|_| None);
    for i in 0..4 {
        let a = &points[(i + 3) % 4];
        let b = &points[i];
        if !seen.insert(b.clone()) || (0..2).any(|j| b[j] != bounds[j][0] && b[j] != bounds[j][1]) {
            return None;
        }
        let wall = if a[0] == b[0] && a[1] != b[1] {
            usize::from(b[0] == bounds[0][1])
        } else if a[1] == b[1] && a[0] != b[0] {
            2 + usize::from(b[1] == bounds[1][1])
        } else {
            return None;
        };
        if sides[wall].replace(ids[i].clone()).is_some() {
            return None;
        }
    }
    Some(sides.map(|s| s.expect("each rectangle side occurs once")))
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
        let sides = vec![
            [(-Rational::one()), Rational::zero(), bounds[0][0].clone()],
            [Rational::one(), Rational::zero(), -&bounds[0][1]],
            [Rational::zero(), -Rational::one(), bounds[1][0].clone()],
            [Rational::zero(), Rational::one(), -&bounds[1][1]],
        ];
        Ok(Self {
            scene,
            frame,
            bounds,
            sides,
            polygon: None,
            circle: None,
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

    /// A closed rectangular enclosure. Polygon branch proofs use this bound;
    /// disk proofs use the radial domain returned by `circular_radius`.
    pub fn bounds(&self) -> &[[Rational; 2]; 2] {
        &self.bounds
    }

    /// Exact chart halfplanes `[a, b, c]`, with `a*u + b*v + c <= 0`
    /// inside the polygon, in authored side-id order. Empty for a disk, whose
    /// complete domain is defined by `circular_radius` instead.
    pub fn side_constraints(&self) -> &[[Rational; 3]] {
        &self.sides
    }

    /// The radius of a complete circular profile, centred at the retained
    /// frame origin. Such profiles use `u² + v² <= radius²` instead of linear
    /// side constraints, and have one authored side id.
    pub fn circular_radius(&self) -> Option<&Rational> {
        self.circle.as_ref().map(|c| &c.radius)
    }

    fn profile_location(&self, q: &[Real; 2], dec: &mut Decider) -> Result<Location, Refusal> {
        if let Some(radius) = self.circular_radius() {
            return Ok(
                match sign(
                    &(&q[0] * &q[0] + &q[1] * &q[1] - Real::from(radius * radius)),
                    dec,
                )? {
                    Sign::Positive => Location::Outside,
                    Sign::Zero => Location::Boundary,
                    Sign::Negative => Location::Inside,
                },
            );
        }
        let mut boundary = false;
        for [a, b, c] in &self.sides {
            match sign(
                &(&Real::from(a.clone()) * &q[0]
                    + &Real::from(b.clone()) * &q[1]
                    + Real::from(c.clone())),
                dec,
            )? {
                Sign::Positive => return Ok(Location::Outside),
                Sign::Zero => boundary = true,
                Sign::Negative => (),
            }
        }
        Ok(if boundary {
            Location::Boundary
        } else {
            Location::Inside
        })
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

    /// The exact positive first hit at a point of the certified profile.
    pub fn height(&self, q: &[Real; 2], dec: &mut Decider) -> Result<Real, Refusal> {
        if self.profile_location(q, dec)? == Location::Outside {
            return Err(coverage("height query is outside its certified domain"));
        }
        let [c, b, a] =
            projected::coefficients(&self.scene.surface, &self.frame)?.map(|p| p.eval_real(q));
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
        let profile = self.profile_location(q, dec)?;
        if profile == Location::Outside {
            return Ok(Location::Outside);
        }
        let boundary = profile == Location::Boundary;
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
