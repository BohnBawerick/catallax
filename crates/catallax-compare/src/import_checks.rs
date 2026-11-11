//! Requested imported-part membership and volume comparison.

use crate::{ImportAgreement, ImportPoint};
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::{Decider, Outcome, Rational, Refusal};
use catallax_shapes::RPoint3;
use catallax_topology::{BodyId, Model};

pub(crate) fn check(
    model: &mut Model,
    body: BodyId,
    field: &Plan,
    requested: &catallax_measure::Accuracy,
    dec: &mut Decider,
) -> Result<ImportAgreement, Refusal> {
    let tight = catallax_measure::Accuracy {
        relative: requested.relative.mul_pow2(-2),
        absolute: requested.absolute.mul_pow2(-2),
    };
    let field_volume = catallax_field::mass::volume(field, &tight, dec)?;
    let measured = catallax_measure::volume(model, body, &tight, dec)?;
    if !measured.marks.is_empty() {
        return Err(Refusal::new(
            "compare.import-volume",
            "B-rep volume did not reach the requested bound",
        ));
    }
    let brep_volume = measured.value;
    let difference = (&brep_volume.lo - &field_volume.hi)
        .abs()
        .max((&brep_volume.hi - &field_volume.lo).abs());
    let magnitude = brep_volume.lo.abs().min(field_volume.lo.abs());
    let allowed = (&requested.relative * &magnitude).max(requested.absolute.clone());
    let volume_agrees = difference <= allowed;
    let bounds = catallax_measure::bounding_box(model, body, &Default::default(), dec)?
        .value
        .ok_or_else(|| Refusal::new("compare.import-bounds", "import has no bounded solid"))?;
    let mut points = Vec::new();
    let mut storage = Workspace::default();
    for i in [-1, 1, 2, 3, 5] {
        for j in [-1, 1, 2, 3, 5] {
            for k in [-1, 1, 2, 3, 5] {
                let [x, y, z] =
                    [i, j, k].map(|v| Rational::from_ratio(v, 4).expect("fixed denominator"));
                let factors = [x, y, z];
                let [x, y, z] = std::array::from_fn(|axis| {
                    let span = &bounds.max[axis] - &bounds.min[axis];
                    &bounds.min[axis] + &(&span * &factors[axis])
                });
                let point = RPoint3::new(x, y, z);
                let brep = match catallax_boolean::classify_point(model, body, &point, dec) {
                    Outcome::Proven(catallax_boolean::Containment::Inside) => Location::Inside,
                    Outcome::Proven(catallax_boolean::Containment::Outside) => Location::Outside,
                    Outcome::ProvenCoincident(_) => Location::Boundary,
                    Outcome::Refused(error) => return Err(error),
                    other => {
                        return Err(Refusal::new(
                            "compare.import-membership",
                            format!("B-rep classification unresolved: {other:?}"),
                        ));
                    }
                };
                let field = field.proven_checked(&point, &mut storage, dec)?.location;
                if field == Location::Indeterminate {
                    return Err(Refusal::new(
                        "compare.import-membership",
                        "field classification is indeterminate",
                    ));
                }
                points.push(ImportPoint { point, brep, field });
            }
        }
    }
    Ok(ImportAgreement {
        points,
        brep_volume,
        field_volume,
        volume_difference_upper: difference,
        volume_difference_allowed: allowed,
        volume_agrees,
    })
}
