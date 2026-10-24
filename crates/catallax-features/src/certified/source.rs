//! Proves the exact recipe target before fitting: a circular sweep clipped by
//! an axial box, including its inactive lateral bounds and original oblique caps.

use super::{geometry, sign, unproven, unsupported};
use crate::Profile;
use catallax_numbers::{Decider, MPoly, Rational, Real, Refusal, Sign};
use catallax_recipe::{BooleanKind, ElementId, Operation, Origin, Recipe, Role, StepId};
use catallax_shapes::{EllipticCylinder, RPoint3, RVector3, Surface, Vector3};
use catallax_topology::{BodyId, Model};
use std::{collections::BTreeMap, sync::Arc};

#[derive(Clone, Debug)]
pub(super) struct PrismSource {
    pub original: Surface,
    pub axis: Vector3,
    pub axial: [Real; 2],
    /// Lower cap, upper cap, side; offsets preserve these selectable origins.
    pub origins: Arc<[Vec<Origin>; 3]>,
}

fn operation<'a>(recipe: &'a Recipe, step: &StepId) -> Result<&'a Operation, Refusal> {
    recipe
        .steps
        .iter()
        .find(|s| &s.id == step)
        .map(|s| &s.op)
        .ok_or_else(|| unproven("retained recipe is missing a source step"))
}

fn normalize(polynomial: &MPoly) -> Result<MPoly, Refusal> {
    let terms = polynomial.terms();
    let coefficient = &terms
        .first()
        .ok_or_else(|| unproven("source support polynomial vanished"))?
        .1;
    Ok(polynomial.scale(&(&Rational::one() / coefficient)))
}

fn box_role(axis: usize, high: bool) -> Role {
    match (axis, high) {
        (2, false) => Role::Start,
        (2, true) => Role::End,
        _ => Role::Side {
            segment: ElementId::new(match (axis, high) {
                (0, false) => "x_min",
                (0, true) => "x_max",
                (1, false) => "y_min",
                (1, true) => "y_max",
                _ => unreachable!("three box axes"),
            })
            .expect("fixed stock role"),
            path: None,
        },
    }
}

/// Support radius of the complete analytic ellipse in a supplied direction.
/// Raw perpendicular basis vectors avoid assuming rational unit principal axes.
pub(super) fn support_radius(original: &EllipticCylinder, direction: &Vector3) -> Real {
    let major = &original.major_dir;
    let minor = original.axis.cross(major);
    (original.major_sq.clone() * (direction.dot(major).square() / major.norm_squared())
        + original.minor_sq.clone() * (direction.dot(&minor).square() / minor.norm_squared()))
    .sqrt()
}

/// Proves a box only clips the axial interval. No mark, sample point, fitted
/// boundary or bounding sphere substitutes for its four complete lateral bounds.
pub(super) fn axial_box(
    original: &EllipticCylinder,
    axis: &Vector3,
    op: &Operation,
    step: &StepId,
    dec: &mut Decider,
) -> Result<([Real; 2], [Vec<Origin>; 2]), Refusal> {
    let Operation::Block { size, pose } = op else {
        return Err(unsupported(
            "the supported planar trim is an axial stock box",
        ));
    };
    let placement = crate::profile::placement(pose, dec)?;
    let corner = crate::profile::rational_point(&placement, &RPoint3::origin()).to_real();
    let delta = original.origin.minus(&corner);
    let mut axial = None;
    for (index, extent) in size.iter().enumerate() {
        dec.check_deadline()?;
        if extent.value() <= &Rational::zero() {
            return Err(unproven("clipping box has a nonpositive extent"));
        }
        let basis = std::array::from_fn(|i| {
            if i == index {
                Rational::one()
            } else {
                Rational::zero()
            }
        });
        let [x, y, z] = basis;
        let direction =
            crate::profile::rational_vector(&placement, &RVector3::new(x, y, z)).to_real();
        let along = direction.dot(axis);
        let center = direction.dot(&delta);
        if sign(&direction.cross(axis).norm_squared(), dec)? == Sign::Zero {
            if axial.is_some() {
                return Err(unproven("box has two parallel axial directions"));
            }
            let positive = match sign(&along, dec)? {
                Sign::Positive => true,
                Sign::Negative => false,
                Sign::Zero => return Err(unproven("axial clipping direction vanished")),
            };
            let a = -center.clone() / along.clone();
            let b = (Real::from_rational(extent.value().clone()) - center) / along;
            let bounds = if positive { [a, b] } else { [b, a] };
            let origins = [!positive, positive].map(|high| {
                vec![Origin {
                    step: step.clone(),
                    role: box_role(index, high),
                }]
            });
            axial = Some((bounds, origins));
        } else {
            if sign(&along, dec)? != Sign::Zero {
                return Err(unsupported(
                    "clipping box is not aligned with the exact sweep axis",
                ));
            }
            let radius = support_radius(original, &direction);
            if sign(&(&center - &radius), dec)? == Sign::Negative
                || sign(
                    &(&center + &radius - Real::from_rational(extent.value().clone())),
                    dec,
                )? == Sign::Positive
            {
                return Err(unsupported("a lateral box plane cuts the retained ellipse"));
            }
        }
    }
    axial.ok_or_else(|| unsupported("clipping box supplies no exact axial cap pair"))
}

pub(super) fn recognize(
    recipe: &Recipe,
    input: &StepId,
    profiles: &BTreeMap<StepId, Profile>,
    model: &Model,
    body: BodyId,
    dec: &mut Decider,
) -> Result<PrismSource, Refusal> {
    if !model.body(body).marks.is_empty() {
        return Err(unproven(
            "an unscoped marked source cannot be treated as an exact prism",
        ));
    }
    let Operation::Boolean {
        kind: BooleanKind::Intersect,
        target,
        tool,
    } = operation(recipe, input)?
    else {
        return Err(unsupported(
            "first offset family needs an exact circular sweep clipped by an axial box",
        ));
    };
    let (sweep_step, block_step) = match (operation(recipe, target)?, operation(recipe, tool)?) {
        (Operation::Extrude { .. }, Operation::Block { .. }) => (target, tool),
        (Operation::Block { .. }, Operation::Extrude { .. }) => (tool, target),
        _ => {
            return Err(unsupported(
                "source intersection is not a circular extrusion and stock box",
            ));
        }
    };
    let Operation::Extrude {
        profile,
        distance,
        direction,
    } = operation(recipe, sweep_step)?
    else {
        unreachable!("matched extrusion source");
    };
    let profile = profiles
        .get(profile)
        .ok_or_else(|| unproven("missing validated source sketch"))?;
    let [boundary] = profile.loops.as_slice() else {
        return Err(unsupported(
            "source sketch must be one solid circular region without holes",
        ));
    };
    let [segment] = boundary.segments.as_slice() else {
        return Err(unsupported(
            "source sketch is not one complete circular carrier",
        ));
    };
    let Some((center, radius_sq, _)) = &segment.circle else {
        return Err(unsupported("source sketch does not define a circle"));
    };
    if !segment.full || boundary.depth != 0 {
        return Err(unsupported(
            "source circle does not bound a complete material disk",
        ));
    }
    let normal = crate::profile::rational_vector(
        &profile.pose,
        &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
    );
    let raw_direction = direction.as_ref().map(|v| {
        RVector3::new(
            v.0[0].value().clone(),
            v.0[1].value().clone(),
            v.0[2].value().clone(),
        )
    });
    let raw_direction = raw_direction.as_ref().unwrap_or(&normal);
    let length = crate::extrude::sqrt(&raw_direction.norm_squared()).ok_or_else(|| {
        unsupported("source direction needs the existing algebraic-placement extension")
    })?;
    if length.is_zero() || distance.value().is_zero() {
        return Err(unproven("source sweep has zero length"));
    }
    let vector = raw_direction.scaled(&(distance.value() / &length));
    let height = normal.dot(&vector);
    if height.is_zero() {
        return Err(unproven("source sweep lies in its sketch plane"));
    }
    let mut original = None;
    for face in model.faces_of_body(body) {
        if let Surface::EllipticCylinder(surface) = model.geom().surface(model.face(face).surface)
            && original.replace(surface).is_some()
        {
            return Err(unproven("source body has more than one elliptic side"));
        }
    }
    let original =
        original.ok_or_else(|| unsupported("source body has no strict elliptic sweep side"))?;
    let original = EllipticCylinder::new(
        original.origin.clone(),
        original.axis.clone(),
        original.major_dir.clone(),
        original.major_sq.clone(),
        original.minor_sq.clone(),
        dec,
    )?;
    let original_surface = Surface::EllipticCylinder(original.clone());
    let circle_center = profile.world(center, &RVector3::zero());
    let coordinates = [&circle_center.x, &circle_center.y, &circle_center.z];
    let delta: [MPoly; 3] =
        std::array::from_fn(|i| MPoly::var(3, i) - MPoly::constant(3, coordinates[i].clone()));
    let n = [&normal.x, &normal.y, &normal.z];
    let v = [&vector.x, &vector.y, &vector.z];
    let along = delta
        .iter()
        .zip(n)
        .fold(MPoly::constant(3, Rational::zero()), |sum, (p, n)| {
            sum + p.scale(n)
        });
    let projected: [MPoly; 3] =
        std::array::from_fn(|i| &delta[i] - &along.scale(&(v[i] / &height)));
    let expected = projected
        .iter()
        .fold(MPoly::constant(3, -radius_sq), |sum, p| sum + p * p);
    let actual = original_surface
        .implicit()
        .and_then(|p| p.into_rational())
        .ok_or_else(|| unproven("exact source support is not the rational sweep law"))?;
    if normalize(&actual)? != normalize(&expected)? {
        return Err(unproven(
            "actual elliptic side does not represent the immutable circular sweep",
        ));
    }
    let axis = original
        .axis
        .scaled(&(Real::one() / original.axis.norm_squared().sqrt()));
    let (axial, caps) = axial_box(
        &original,
        &axis,
        operation(recipe, block_step)?,
        block_step,
        dec,
    )?;
    let n = normal.to_real();
    let center_projection = original.origin.minus(&circle_center.to_real()).dot(&n);
    let travel = axis.dot(&n);
    let endpoints = axial
        .each_ref()
        .map(|z| &center_projection + &(z * &travel));
    let (lower, upper) = if sign(&travel, dec)? == Sign::Negative {
        (&endpoints[1], &endpoints[0])
    } else {
        (&endpoints[0], &endpoints[1])
    };
    let radius = support_radius(&original, &n);
    let source_lower = Rational::zero().min(height.clone());
    let source_upper = Rational::zero().max(height);
    if sign(&(lower - &radius - Real::from_rational(source_lower)), dec)? == Sign::Negative
        || sign(&(upper + &radius - Real::from_rational(source_upper)), dec)? == Sign::Positive
    {
        return Err(unsupported(
            "the axial prism reaches an original oblique source cap",
        ));
    }
    let mut side: Vec<_> = segment
        .ids
        .iter()
        .map(|id| Origin {
            step: sweep_step.clone(),
            role: Role::Side {
                segment: id.clone(),
                path: None,
            },
        })
        .collect();
    side.sort();
    side.dedup();
    let [lower, upper] = caps;
    let origins = Arc::new([lower, upper, side]);
    geometry::check_prism(
        model,
        body,
        &original,
        axial.each_ref(),
        geometry::ProfileSupport::Analytic,
        &origins,
        dec,
    )?;
    Ok(PrismSource {
        original: original_surface,
        axis,
        axial,
        origins,
    })
}
