//! Exact analytic draft supports from RECIPE.md section 4.10.
//!
//! A support alone is not a drafted solid. The arrangement builder must retain
//! its nappe, material orientation, trims and neutral-connected components.

use catallax_numbers::{
    Decider, Outcome, Rational, Real, Refusal, Sign, StrictSign, sin_cos_degrees,
};
use catallax_shapes::{Cone, Plane, Point3, Surface, Vector3};

/// Exact hollow axial cells with strictly separated radial walls.
#[path = "draft_annular.rs"]
pub mod annular;
/// Oriented multi-cell planar support arrangements and patch provenance.
#[path = "draft_arrangement.rs"]
pub mod arrangement;
/// One unchanged bore strictly inside a convex planar draft envelope.
#[path = "draft_bored.rs"]
pub mod bored;
/// Exact cylindrical cells with independently extended planar caps.
#[path = "draft_capped.rs"]
pub mod capped;
/// Exact bounded cells used by analytic support arrangements.
#[path = "draft_cell.rs"]
pub mod cell;
/// Exact nearest-face distances for plane-capped axial cylinders.
#[path = "draft_near.rs"]
pub mod near;
/// Exact axial cells bounded by drafted cylindrical supports.
#[path = "draft_round.rs"]
pub mod round;

pub(crate) fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    if let Some(value) = value.rational_ref() {
        return Ok(value.sign());
    }
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        Outcome::Refused(r) => Err(r),
    }
}

fn unit(vector: &Vector3, dec: &mut Decider) -> Result<Vector3, Refusal> {
    let length = reduced(vector.norm_squared().sqrt(), dec)?;
    let scaled = vector.scaled(&(&Real::one() / &length));
    Ok(Vector3::new(
        reduced(scaled.x, dec)?,
        reduced(scaled.y, dec)?,
        reduced(scaled.z, dec)?,
    ))
}

fn reduced(value: Real, dec: &mut Decider) -> Result<Real, Refusal> {
    if value.rational_ref().is_some() {
        return Ok(value);
    }
    Ok(dec.rational_value(&value)?.map(Real::from).unwrap_or(value))
}

fn reduced_point(point: Point3, dec: &mut Decider) -> Result<Point3, Refusal> {
    Ok(Point3::new(
        reduced(point.x, dec)?,
        reduced(point.y, dec)?,
        reduced(point.z, dec)?,
    ))
}

fn homogeneous_direction(vector: Vector3, dec: &mut Decider) -> Result<Vector3, Refusal> {
    let mut magnitude = None;
    for component in [&vector.x, &vector.y, &vector.z] {
        match sign(component, dec)? {
            Sign::Zero => (),
            Sign::Positive => {
                magnitude = Some(component.clone());
                break;
            }
            Sign::Negative => {
                magnitude = Some(-component);
                break;
            }
        }
    }
    let magnitude = magnitude.expect("validated nonzero draft direction");
    let scaled = vector.scaled(&(&Real::one() / &magnitude));
    Ok(Vector3::new(
        reduced(scaled.x, dec)?,
        reduced(scaled.y, dec)?,
        reduced(scaled.z, dec)?,
    ))
}

/// The saved pull, neutral plane and exact algebraic angle of a draft.
#[derive(Clone, Debug)]
pub struct DraftDefinition {
    /// Unit pull direction, retaining algebraic normalisation exactly.
    pub pull: Vector3,
    /// The authored neutral plane, never replaced with an axial section.
    pub neutral: Plane,
    /// Sine of the authored rational angle in degrees.
    pub sine: Real,
    /// Cosine of the authored rational angle in degrees.
    pub cosine: Real,
}

/// One oriented replacement support. It retains the source face's origins.
#[derive(Clone, Debug)]
pub struct DraftSupport {
    /// Exact replacement plane, cylinder or single cone nappe.
    pub surface: Surface,
    /// Whether the material is on the negative side of the support.
    pub same_sense: bool,
}

impl DraftDefinition {
    /// Validates the authored domains before any empty-body simplification.
    pub fn new(
        pull: Vector3,
        neutral: Plane,
        angle: &Rational,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if angle.abs() >= Rational::from_i64(90) {
            return Err(Refusal::new(
                "draft.angle",
                "draft angle must lie strictly between -90 and 90 degrees",
            ));
        }
        if sign(&pull.norm_squared(), dec)? != Sign::Positive {
            return Err(Refusal::new(
                "draft.zero-pull",
                "draft pull direction must be nonzero",
            ));
        }
        if sign(&neutral.normal.norm_squared(), dec)? != Sign::Positive {
            return Err(Refusal::new(
                "draft.zero-neutral-normal",
                "neutral plane normal must be nonzero",
            ));
        }
        let [sine, cosine] = sin_cos_degrees(angle, dec)?;
        Ok(Self {
            pull: unit(&pull, dec)?,
            neutral,
            sine,
            cosine,
        })
    }

    /// Constructs only the support specified by the neutral-curve generators.
    /// This does not extend neighbours, select material or publish topology.
    pub fn support(
        &self,
        source: &Surface,
        same_sense: bool,
        dec: &mut Decider,
    ) -> Result<DraftSupport, Refusal> {
        match source {
            Surface::Plane(plane) => self.plane(plane, same_sense, dec),
            Surface::Cylinder(cylinder) => {
                if sign(&cylinder.axis.cross(&self.pull).norm_squared(), dec)? != Sign::Zero {
                    return Err(Refusal::new(
                        "draft.unsupported-cylinder-axis",
                        "non-axial cylindrical draft needs a derived ruled support",
                    ));
                }
                let normal_rate = self.pull.dot(&self.neutral.normal);
                let zero_angle = sign(&self.sine, dec)? == Sign::Zero;
                if sign(&normal_rate, dec)? == Sign::Zero {
                    let offset = self
                        .neutral
                        .origin
                        .minus(&cylinder.origin)
                        .dot(&self.neutral.normal);
                    if sign(
                        &(offset.square()
                            - &cylinder.radius_sq * &self.neutral.normal.norm_squared()),
                        dec,
                    )? == Sign::Positive
                    {
                        return Err(Refusal::new(
                            "draft.missing-neutral",
                            "selected cylindrical support misses the neutral plane",
                        ));
                    }
                    if zero_angle {
                        return Err(Refusal::new(
                            "draft.degenerate-support",
                            "axial neutral generators at zero angle sweep only lines",
                        ));
                    }
                }
                // At zero angle every axial generator meeting N reproduces
                // the original cylinder, even when its neutral curve is oblique.
                if zero_angle {
                    return Ok(DraftSupport {
                        surface: source.clone(),
                        same_sense,
                    });
                }
                if sign(&self.neutral.normal.cross(&self.pull).norm_squared(), dec)? != Sign::Zero {
                    return Err(Refusal::new(
                        "draft.unsupported-oblique-neutral",
                        "an oblique cylindrical neutral curve does not reduce to an exact cone",
                    ));
                }
                let radius = reduced(cylinder.radius_sq.sqrt(), dec)?;
                let centre = cylinder.origin.plus(
                    &self.pull.scaled(
                        &(&self
                            .neutral
                            .origin
                            .minus(&cylinder.origin)
                            .dot(&self.neutral.normal)
                            / &normal_rate),
                    ),
                );
                let signed_sine = if same_sense {
                    self.sine.clone()
                } else {
                    -&self.sine
                };
                let axis = self.pull.scaled(&Real::from_i64(
                    if sign(&signed_sine, dec)? == Sign::Positive {
                        -1
                    } else {
                        1
                    },
                ));
                let slope = &self.sine / &self.cosine;
                let slope = if sign(&slope, dec)? == Sign::Negative {
                    -slope
                } else {
                    slope
                };
                let cotangent = reduced(&self.cosine / &signed_sine, dec)?;
                let apex =
                    reduced_point(centre.plus(&self.pull.scaled(&(&radius * &cotangent))), dec)?;
                Ok(DraftSupport {
                    surface: Surface::Cone(Cone::new(
                        apex,
                        axis,
                        reduced(slope, dec)?,
                        cylinder.ref_dir.clone(),
                        dec,
                    )?),
                    same_sense,
                })
            }
            _ => Err(Refusal::new(
                "draft.unsupported-surface",
                "analytic draft supports are implemented only for planes and pull-parallel cylinders",
            )),
        }
    }

    fn plane(
        &self,
        plane: &Plane,
        same_sense: bool,
        dec: &mut Decider,
    ) -> Result<DraftSupport, Refusal> {
        let normal = plane
            .normal
            .scaled(&Real::from_i64(if same_sense { 1 } else { -1 }));
        if sign(&normal.norm_squared(), dec)? != Sign::Positive {
            return Err(Refusal::new(
                "draft.zero-face-normal",
                "selected face normal must be nonzero",
            ));
        }
        let radial = normal.minus(&self.pull.scaled(&normal.dot(&self.pull)));
        if sign(&radial.norm_squared(), dec)? == Sign::Zero {
            return Err(Refusal::new(
                "draft.parallel-face-normal",
                "selected face normal is parallel to pull on its neutral curve",
            ));
        }
        let generator = self
            .pull
            .scaled(&self.cosine)
            .minus(&unit(&radial, dec)?.scaled(&self.sine));
        let line = normal.cross(&self.neutral.normal);
        if sign(&line.norm_squared(), dec)? == Sign::Zero {
            let separation = self.neutral.origin.minus(&plane.origin).dot(&normal);
            if sign(&separation, dec)? == Sign::Zero
                && sign(&normal.dot(&generator), dec)? == Sign::Zero
            {
                return Ok(DraftSupport {
                    surface: Surface::Plane(plane.clone()),
                    same_sense,
                });
            }
            return Err(if sign(&separation, dec)? == Sign::Zero {
                Refusal::new(
                    "draft.degenerate-neutral",
                    "coincident support and neutral plane do not define a neutral line",
                )
            } else {
                Refusal::new(
                    "draft.missing-neutral",
                    "selected support misses the neutral plane",
                )
            });
        }
        // Move within the source plane, perpendicular to the neutral line.
        let across = self
            .neutral
            .normal
            .minus(&normal.scaled(&(&self.neutral.normal.dot(&normal) / &normal.norm_squared())));
        let origin = plane.origin.plus(
            &across.scaled(
                &(&self
                    .neutral
                    .origin
                    .minus(&plane.origin)
                    .dot(&self.neutral.normal)
                    / &across.norm_squared()),
            ),
        );
        let mut outward = generator.cross(&line);
        if sign(&outward.norm_squared(), dec)? == Sign::Zero {
            return Err(Refusal::new(
                "draft.degenerate-support",
                "draft generators collapse onto the neutral line",
            ));
        }
        // Preserve the source's material-outward radial direction. An oblique
        // neutral can cross the generator direction as the angle changes; using
        // its trace normal here would silently reverse the solid at that event.
        match sign(&outward.dot(&radial), dec)? {
            Sign::Negative => outward = outward.scaled(&Real::from_i64(-1)),
            Sign::Zero => {
                return Err(Refusal::new(
                    "draft.unsupported-neutral-orientation",
                    "neutral curve is parallel to the source's pull-tangent direction",
                ));
            }
            Sign::Positive => (),
        }
        Ok(DraftSupport {
            surface: Surface::Plane(Plane {
                origin: reduced_point(origin, dec)?,
                // Positive homogeneous rescaling removes common algebraic
                // factors before Boolean predicates see the support equation.
                normal: homogeneous_direction(outward, dec)?,
                ref_dir: homogeneous_direction(line, dec)?,
            }),
            same_sense: true,
        })
    }
}

/// Converts the recipe's exact rational neutral plane without fitting a frame.
pub fn neutral_plane(plane: &catallax_recipe::Plane) -> Plane {
    let [x, y, z] = plane.point.0.each_ref().map(|n| n.value().clone().into());
    let [nx, ny, nz] = plane.normal.0.each_ref().map(|n| n.value().clone().into());
    let normal = Vector3::new(nx, ny, nz);
    let basis = if plane.normal.0[0].value().is_zero() {
        Vector3::new(Real::one(), Real::zero(), Real::zero())
    } else {
        Vector3::new(Real::zero(), Real::one(), Real::zero())
    };
    let ref_dir = normal.cross(&basis);
    Plane {
        origin: Point3::new(x, y, z),
        normal,
        ref_dir,
    }
}
