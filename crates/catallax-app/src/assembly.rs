//! One reported enclosure for independent answers on the same placed definition.
//!
//! Agreement is at the requested numeric resolution, not a proof that arbitrary
//! models have equal geometry. Callers must query the same recipe, placement and
//! density in both kernels. The reported hull encloses both values, even if the
//! underlying definitions differ within their construction bounds. Disjoint or
//! excessively wide enclosures refuse; intersecting intervals alone do not pass.

use catallax_features::Built;
pub use catallax_field::assembly::Mass;
use catallax_field::{assembly::Relation, mass};
use catallax_measure::{Accuracy, Enclosure};
use catallax_numbers::{Answer, Decider, Rational, Refusal};

use crate::queries::{BodyRelation, QueryRefusal, query};

/// Uniform mass of the retained B-rep without requesting inertia. A field-only
/// region envelope cannot supply mass for its lattice or gyroid material.
/// The caller's volume and centroid accuracies keep their respective units.
pub fn uniform_mass(
    built: &Built,
    density: &Rational,
    accuracy: &mass::Accuracy,
    dec: &mut Decider,
) -> Result<Answer<Mass>, QueryRefusal> {
    query(&built.model, &[built.result], || {
        if density <= &Rational::zero() {
            return Err(Refusal::new(
                "assembly.density",
                "density must be positive in g/mm^3",
            ));
        }
        validate(&accuracy.volume)?;
        validate(&accuracy.centroid)?;
        if !built.stand_ins.is_empty() {
            return Err(Refusal::new(
                "assembly.mass-envelope",
                format!(
                    "{} field-only steps have region envelopes, not material mass",
                    built.stand_ins.len()
                ),
            ));
        }
        let volume = catallax_measure::volume(&built.model, built.result, &accuracy.volume, dec)?;
        let centroid =
            catallax_measure::centroid(&built.model, built.result, &accuracy.centroid, dec)?;
        require_width(&volume.value, &accuracy.volume, "volume in mm^3")?;
        if let Some(centre) = &centroid.value {
            for (i, coordinate) in centre.iter().enumerate() {
                require_width(
                    coordinate,
                    &accuracy.centroid,
                    &format!("centre coordinate {i} in mm"),
                )?;
            }
        }
        let mut marks = volume.marks;
        marks.extend(&centroid.marks);
        Ok(Answer::with_marks(
            Mass {
                grams: Enclosure {
                    lo: &volume.value.lo * density,
                    hi: &volume.value.hi * density,
                },
                centre_of_gravity: centroid.value,
            },
            marks,
        ))
    })
}

fn validate(accuracy: &Accuracy) -> Result<(), Refusal> {
    if accuracy.relative < Rational::zero() || accuracy.absolute < Rational::zero() {
        Err(Refusal::new(
            "assembly.agreement-accuracy",
            "accuracy must be nonnegative",
        ))
    } else {
        Ok(())
    }
}

fn allowed_width(value: &Enclosure, accuracy: &Accuracy) -> Rational {
    let magnitude = if value.lo <= Rational::zero() && value.hi >= Rational::zero() {
        Rational::zero()
    } else {
        value.lo.abs().min(value.hi.abs())
    };
    accuracy
        .absolute
        .clone()
        .max(&accuracy.relative * &magnitude)
}

fn bounds(value: &Enclosure) -> String {
    format!(
        "[{}, {}]",
        value.lo.to_decimal_bounds(10).0,
        value.hi.to_decimal_bounds(10).1
    )
}

fn require_width(value: &Enclosure, accuracy: &Accuracy, quantity: &str) -> Result<(), Refusal> {
    let allowed = allowed_width(value, accuracy);
    if value.width() > allowed {
        Err(Refusal::new(
            "assembly.mass-accuracy",
            format!(
                "{quantity}: enclosure {}; width exceeds {}",
                bounds(value),
                allowed.to_decimal_bounds(10).1
            ),
        ))
    } else {
        Ok(())
    }
}

fn agreed(
    a: &Enclosure,
    b: &Enclosure,
    accuracy: &Accuracy,
    quantity: &str,
) -> Result<Enclosure, Refusal> {
    validate(accuracy)?;
    if a.lo > a.hi || b.lo > b.hi {
        return Err(Refusal::new(
            "assembly.agreement-enclosure",
            "inverted enclosure",
        ));
    }
    if a.hi < b.lo || b.hi < a.lo {
        return Err(Refusal::new(
            "assembly.disagreement",
            format!("{quantity}: B-rep {}; field {}", bounds(a), bounds(b)),
        ));
    }
    let hull = Enclosure {
        lo: a.lo.clone().min(b.lo.clone()),
        hi: a.hi.clone().max(b.hi.clone()),
    };
    let allowed = allowed_width(&hull, accuracy);
    if hull.width() > allowed {
        return Err(Refusal::new(
            "assembly.agreement-width",
            format!(
                "{quantity}: B-rep {}; field {}; combined width {} exceeds {}",
                bounds(a),
                bounds(b),
                hull.width().to_decimal_bounds(10).1,
                allowed.to_decimal_bounds(10).1
            ),
        ));
    }
    Ok(hull)
}

/// Reconcile two proven mass answers for the same density and world placement.
/// Returns one enclosure per quantity, preserving B-rep construction marks.
/// The gram accuracy is the requested volume accuracy multiplied by density.
/// A refusal in either query must be handled before calling this function.
pub fn mass_agreement(
    brep: &Answer<Mass>,
    field: &Mass,
    density: &Rational,
    accuracy: &mass::Accuracy,
) -> Result<Answer<Mass>, Refusal> {
    if density <= &Rational::zero() {
        return Err(Refusal::new(
            "assembly.density",
            "density must be positive in g/mm^3",
        ));
    }
    validate(&accuracy.volume)?;
    validate(&accuracy.centroid)?;
    let grams = agreed(
        &brep.value.grams,
        &field.grams,
        &Accuracy {
            relative: accuracy.volume.relative.clone(),
            absolute: &accuracy.volume.absolute * density,
        },
        "mass in grams",
    )?;
    let centre_of_gravity = match (&brep.value.centre_of_gravity, &field.centre_of_gravity) {
        (None, None) if grams.lo.is_zero() && grams.hi.is_zero() => None,
        (Some(a), Some(b)) if grams.lo > Rational::zero() => Some([
            agreed(&a[0], &b[0], &accuracy.centroid, "centre x in mm")?,
            agreed(&a[1], &b[1], &accuracy.centroid, "centre y in mm")?,
            agreed(&a[2], &b[2], &accuracy.centroid, "centre z in mm")?,
        ]),
        _ => {
            return Err(Refusal::new(
                "assembly.mass-centroid",
                "mass and centroid presence disagree or positive mass is unresolved",
            ));
        }
    };
    Ok(Answer::with_marks(
        Mass {
            grams,
            centre_of_gravity,
        },
        brep.marks.clone(),
    ))
}

/// Reconcile independent clash or clear results on the same placed operands.
/// A field refusal propagates unchanged. A clear answer has one bounded gap.
pub fn relation_agreement(
    brep: &Answer<BodyRelation>,
    field: &Relation,
    accuracy: &Accuracy,
) -> Result<Answer<BodyRelation>, Refusal> {
    validate(accuracy)?;
    let value = match (&brep.value, field) {
        (_, Relation::Unresolved(r)) => return Err(r.clone()),
        (BodyRelation::Clash { overlap }, Relation::Clash { .. })
            if overlap.lo > Rational::zero() =>
        {
            brep.value.clone()
        }
        (BodyRelation::Clear { gap: a }, Relation::Clear { gap: b })
            if a.lo > Rational::zero() && b.lo > Rational::zero() =>
        {
            BodyRelation::Clear {
                gap: agreed(a, b, accuracy, "minimum gap in mm")?,
            }
        }
        _ => {
            return Err(Refusal::new(
                "assembly.relation-disagreement",
                format!("B-rep {:?}; field {field:?}", brep.value),
            ));
        }
    };
    Ok(Answer::with_marks(value, brep.marks.clone()))
}
