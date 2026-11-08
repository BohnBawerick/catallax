//! Whole-ellipse polynomial integrals, with quadrature retained for arc length.

use super::{dot, offset, orientation, parallel, unit};
use crate::integrate::{Builder, Program};
use crate::planar::real_direction;
use crate::surface::Property;
use catallax_numbers::Real;
use catallax_shapes::{Ellipse, RPoint3, Surface};

pub(crate) fn compile(
    surface: &Surface,
    ellipse: &Ellipse,
    anchor: &RPoint3,
    indices: &[usize],
    property: Property,
) -> Option<Program> {
    let first = real_direction(&ellipse.major_dir);
    let across = ellipse.normal.cross(&ellipse.major_dir);
    let second = real_direction(&across);
    let first_squared = ellipse.major_dir.norm_squared();
    let second_squared = across.norm_squared();
    if [&first_squared, &second_squared]
        .iter()
        .any(|n| n.as_rational().is_some_and(|n| n.is_zero()))
    {
        return None;
    }
    // This is four times the filled ellipse's covariance, without unit-axis radicals.
    let tensor: [[Real; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            (&ellipse.major_sq * &(first[i] * first[j])) / first_squared.clone()
                + (&ellipse.minor_sq * &(second[i] * second[j])) / second_squared.clone()
        })
    });
    let area_factor = (&ellipse.major_sq * &ellipse.minor_sq).sqrt();
    let b = Builder::default();
    let pi = b.input(6);
    let factor = b.exact(&area_factor);
    let anchor = anchor.to_real();
    match surface {
        Surface::Plane(plane) => {
            if !parallel(&plane.normal, &ellipse.normal)? {
                return None;
            }
            let normal = unit(&b, &plane.normal);
            let center = offset(&b, &ellipse.center, &anchor);
            let area = orientation(&b, &ellipse.normal, &plane.normal)? * pi * factor;
            let values = match property {
                Property::SecondMoment(i, j) => {
                    let moment = dot(center, normal)
                        * area
                        * (center[i] * center[j] + b.exact(&tensor[i][j]) / b.k(4))
                        / b.k(5);
                    return Some(b.finish(&[moment]));
                }
                Property::FaceCentroid => [
                    area,
                    area,
                    area * center[0],
                    area * center[1],
                    area * center[2],
                ],
                Property::Area | Property::Volume | Property::Centroid => {
                    let volume = dot(center, normal) * area / b.k(3);
                    let moments = std::array::from_fn::<_, 3, _>(|i| {
                        normal[i] * area * (center[i].square() + b.exact(&tensor[i][i]) / b.k(4))
                            / b.k(2)
                    });
                    [area, volume, moments[0], moments[1], moments[2]]
                }
            };
            Some(b.finish(&indices.iter().map(|i| values[*i]).collect::<Vec<_>>()))
        }
        Surface::EllipticCylinder(cylinder) => {
            if !parallel(&cylinder.axis, &ellipse.normal)?
                || !parallel(&cylinder.major_dir, &ellipse.major_dir)?
                || !ellipse
                    .center
                    .minus(&cylinder.origin)
                    .cross(&cylinder.axis)
                    .norm_squared()
                    .as_rational()?
                    .is_zero()
                || !(&ellipse.major_sq - &cylinder.major_sq)
                    .as_rational()?
                    .is_zero()
                || !(&ellipse.minor_sq - &cylinder.minor_sq)
                    .as_rational()?
                    .is_zero()
            {
                return None;
            }
            if matches!(property, Property::Area | Property::FaceCentroid) {
                return None;
            }
            let axis = unit(&b, &cylinder.axis);
            let z = dot(offset(&b, &ellipse.center, &cylinder.origin), axis);
            let origin = offset(&b, &cylinder.origin, &anchor);
            let signed = orientation(&b, &ellipse.normal, &cylinder.axis)? * pi * factor;
            if let Property::SecondMoment(i, j) = property {
                let axial = dot(origin, axis);
                let perpendicular: [_; 3] = std::array::from_fn(|k| origin[k] - axial * axis[k]);
                let second = b.k(2) * origin[i] * origin[j]
                    + b.exact(&tensor[i][j])
                    + origin[i] * perpendicular[j]
                    + origin[j] * perpendicular[i];
                let first = |k| b.k(2) * origin[k] + perpendicular[k];
                let moment = -signed
                    * (second * z
                        + (axis[j] * first(i) + axis[i] * first(j)) * z.square() / b.k(2)
                        + b.k(2) * axis[i] * axis[j] * z.square() * z / b.k(3))
                    / b.k(5);
                return Some(b.finish(&[moment]));
            }
            let volume = -b.k(2) * signed * z / b.k(3);
            let moments = std::array::from_fn::<_, 3, _>(|i| {
                -signed
                    * (b.k(1) - axis[i].square())
                    * (origin[i] * z + axis[i] * z.square() / b.k(2))
            });
            let values = [b.k(0), volume, moments[0], moments[1], moments[2]];
            Some(b.finish(&indices.iter().map(|i| values[*i]).collect::<Vec<_>>()))
        }
        _ => None,
    }
}
