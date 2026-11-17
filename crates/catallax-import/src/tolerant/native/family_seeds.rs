//! Initial guesses for sampled native projection; no distance or chart authority.

use crate::geometry_source::{
    AngleUnit, Curve, EntityId, Frame, Geometry, MappingError, ParsedStep, Surface,
};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{Point3, Vector3};

#[derive(Clone, Debug)]
enum Kind {
    Torus {
        frame: Frame,
        major_mm: Rational,
        spindle: bool,
    },
    ConicExtrusion {
        frame: Frame,
        radii_mm: [Rational; 2],
        sweep_mm: Vector3,
    },
    LineExtrusion {
        origin: Point3,
        tangent_mm: Vector3,
        sweep_mm: Vector3,
    },
}

#[derive(Clone, Debug)]
pub(crate) struct FamilySeeds {
    kind: Kind,
    angle_unit: AngleUnit,
    offset: bool,
}

fn invalid(detail: &str) -> MappingError {
    Refusal::new("import.native-family-seed", detail).into()
}

impl FamilySeeds {
    pub(crate) fn from_source(
        source: &ParsedStep,
        mut id: EntityId,
        length_mm: &Rational,
        angle_unit: &AngleUnit,
        dec: &mut Decider,
    ) -> Result<Option<Self>, MappingError> {
        let mut seen = Vec::new();
        let mut offset = false;
        let kind = loop {
            dec.charge_work(1)?;
            if seen.contains(&id) {
                return Err(invalid("cyclic offset basis"));
            }
            seen.push(id);
            match source.geometry(id, dec)? {
                Geometry::Surface(Surface::Offset(parent, _, _)) => {
                    id = parent;
                    offset = true;
                }
                Geometry::Surface(Surface::Torus(placement, major, minor)) => {
                    if major.sign() != Sign::Positive || minor.sign() != Sign::Positive {
                        return Ok(None);
                    }
                    break Kind::Torus {
                        frame: source.axis_placement(placement, length_mm, dec)?,
                        major_mm: &major * length_mm,
                        spindle: minor > major,
                    };
                }
                Geometry::Surface(Surface::Extrusion(profile, vector)) => {
                    let Geometry::Vector(direction, magnitude) = source.geometry(vector, dec)?
                    else {
                        return Err(invalid("extrusion has no VECTOR"));
                    };
                    if magnitude.is_zero() {
                        return Err(invalid("zero extrusion VECTOR"));
                    }
                    let sign = if magnitude.sign() == Sign::Negative {
                        -1
                    } else {
                        1
                    };
                    let sweep_mm = source
                        .unit_direction(direction, dec)?
                        .scaled(&Real::from(length_mm * &Rational::from_i64(sign)));
                    match source.geometry(profile, dec)? {
                        Geometry::Curve(Curve::Circle(placement, radius)) => {
                            let radius_mm = &radius * length_mm;
                            break Kind::ConicExtrusion {
                                frame: source.axis_placement(placement, length_mm, dec)?,
                                radii_mm: [radius_mm.clone(), radius_mm],
                                sweep_mm,
                            };
                        }
                        Geometry::Curve(Curve::Ellipse(placement, a, b)) => {
                            break Kind::ConicExtrusion {
                                frame: source.axis_placement(placement, length_mm, dec)?,
                                radii_mm: [&a * length_mm, &b * length_mm],
                                sweep_mm,
                            };
                        }
                        Geometry::Curve(Curve::Line(_, _)) => {
                            let catallax_shapes::Curve::Line(line) =
                                source.native_curve(profile, length_mm, dec)?
                            else {
                                return Err(invalid("line profile did not retain a line"));
                            };
                            break Kind::LineExtrusion {
                                origin: line.origin,
                                tangent_mm: line.dir,
                                sweep_mm,
                            };
                        }
                        _ => return Ok(None),
                    }
                }
                _ => return Ok(None),
            }
        };
        Ok(Some(Self {
            kind,
            angle_unit: angle_unit.clone(),
            offset,
        }))
    }

    pub(crate) fn seeds(
        &self,
        point_mm: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Vec<[Rational; 2]>, MappingError> {
        let bits = dec.policy().max_precision_bits.min(128);
        if bits == 0 {
            return Ok(Vec::new());
        }
        dec.charge_work(32)?;
        let pi = Ball::pi(bits);
        let angular_scale = match &self.angle_unit {
            AngleUnit::Radians => Ball::from_rational(&Rational::one(), bits),
            AngleUnit::Degrees => {
                &pi * &Ball::from_rational(
                    &Rational::from_ratio(1, 180).expect("nonzero denominator"),
                    bits,
                )
            }
            AngleUnit::ScaledRadians(scale) => Ball::from_rational(scale, bits),
        };
        let Some(inverse_angle) =
            Ball::from_rational(&Rational::one(), bits).checked_div(&angular_scale)
        else {
            return Ok(Vec::new());
        };
        let mut result = Vec::with_capacity(4);
        let offset_axis = match &self.kind {
            Kind::Torus {
                frame,
                major_mm,
                spindle,
            } => {
                let local = local(point_mm, frame, bits, dec)?;
                let Some(u) = angle(&local[0], &local[1], bits) else {
                    return Ok(result);
                };
                let Some(radial) = (local[0].square() + local[1].square()).sqrt() else {
                    return Ok(result);
                };
                let major = Ball::from_rational(major_mm, bits);
                if let Some(v) = angle(&(&radial - &major), &local[2], bits) {
                    result.push([(&u * &inverse_angle).mid(), (&v * &inverse_angle).mid()]);
                }
                if *spindle {
                    let inner_radial = -&radial;
                    if let Some(v) = angle(&(&inner_radial - &major), &local[2], bits) {
                        let inner_u = &u + &pi;
                        result.push([
                            (&inner_u * &inverse_angle).mid(),
                            (&v * &inverse_angle).mid(),
                        ]);
                    }
                }
                Some(1)
            }
            Kind::ConicExtrusion {
                frame,
                radii_mm,
                sweep_mm,
            } => {
                let delta = relative(point_mm, &frame.origin, bits, dec)?;
                let z = vector(&frame.axes[2], bits, dec)?;
                let sweep = vector(sweep_mm, bits, dec)?;
                let Some(v) = dot(&delta, &z).checked_div(&dot(&sweep, &z)) else {
                    return Ok(result);
                };
                let projected = std::array::from_fn(|i| &delta[i] - &(&sweep[i] * &v));
                let x_axis = vector(&frame.axes[0], bits, dec)?;
                let y_axis = vector(&frame.axes[1], bits, dec)?;
                let Some(x) =
                    dot(&projected, &x_axis).checked_div(&Ball::from_rational(&radii_mm[0], bits))
                else {
                    return Ok(result);
                };
                let Some(y) =
                    dot(&projected, &y_axis).checked_div(&Ball::from_rational(&radii_mm[1], bits))
                else {
                    return Ok(result);
                };
                if let Some(u) = angle(&x, &y, bits) {
                    result.push([(&u * &inverse_angle).mid(), v.mid()]);
                }
                Some(0)
            }
            Kind::LineExtrusion {
                origin,
                tangent_mm,
                sweep_mm,
            } => {
                let delta = relative(point_mm, origin, bits, dec)?;
                let u = vector(tangent_mm, bits, dec)?;
                let v = vector(sweep_mm, bits, dec)?;
                let a = dot(&u, &u);
                let b = dot(&u, &v);
                let c = dot(&v, &v);
                let determinant = &a * &c - b.square();
                let r = dot(&u, &delta);
                let s = dot(&v, &delta);
                if let (Some(u), Some(v)) = (
                    (&c * &r - &b * &s).checked_div(&determinant),
                    (&a * &s - &b * &r).checked_div(&determinant),
                ) {
                    result.push([u.mid(), v.mid()]);
                }
                None
            }
        };
        if self.offset
            && let Some(axis) = offset_axis
        {
            let half_turn = (&pi * &inverse_angle).mid();
            let count = result.len();
            for i in 0..count {
                let mut alternate = result[i].clone();
                alternate[axis] = &alternate[axis] + &half_turn;
                result.push(alternate);
            }
        }
        Ok(result)
    }
}

fn vector(v: &Vector3, bits: u32, dec: &mut Decider) -> Result<[Ball; 3], MappingError> {
    Ok([
        dec.enclose_budgeted_at(&v.x, bits)?,
        dec.enclose_budgeted_at(&v.y, bits)?,
        dec.enclose_budgeted_at(&v.z, bits)?,
    ])
}

fn relative(
    p: &[Rational; 3],
    origin: &Point3,
    bits: u32,
    dec: &mut Decider,
) -> Result<[Ball; 3], MappingError> {
    let origin = [
        dec.enclose_budgeted_at(&origin.x, bits)?,
        dec.enclose_budgeted_at(&origin.y, bits)?,
        dec.enclose_budgeted_at(&origin.z, bits)?,
    ];
    Ok(std::array::from_fn(|i| {
        let coordinate = Ball::from_rational(&p[i], bits);
        &coordinate - &origin[i]
    }))
}

fn dot(a: &[Ball; 3], b: &[Ball; 3]) -> Ball {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}

fn local(
    p: &[Rational; 3],
    frame: &Frame,
    bits: u32,
    dec: &mut Decider,
) -> Result<[Ball; 3], MappingError> {
    let delta = relative(p, &frame.origin, bits, dec)?;
    Ok([
        dot(&delta, &vector(&frame.axes[0], bits, dec)?),
        dot(&delta, &vector(&frame.axes[1], bits, dec)?),
        dot(&delta, &vector(&frame.axes[2], bits, dec)?),
    ])
}

fn angle(x: &Ball, y: &Ball, bits: u32) -> Option<Ball> {
    // A seed uses a single representative, not an angular enclosure or proof.
    Ball::atan2(
        &Ball::from_rational(&y.mid(), bits),
        &Ball::from_rational(&x.mid(), bits),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    fn q(value: &str) -> Rational {
        Rational::parse(value).unwrap()
    }

    #[test]
    fn source_units_inner_sheets_and_signed_offsets_keep_the_original_geometry() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let source = ParsedStep::read(
                b"DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=AXIS2_PLACEMENT_3D('',#1,$,$);
#3=TOROIDAL_SURFACE('',#2,4.,1.);
#4=TOROIDAL_SURFACE('',#2,1.,2.);
#5=ELLIPSE('',#2,2.,1.);
#6=DIRECTION('',(1.,0.,1.));
#7=VECTOR('',#6,7.);
#8=SURFACE_OF_LINEAR_EXTRUSION('',#5,#7);
#9=OFFSET_SURFACE('',#8,1.,.F.);
#10=OFFSET_SURFACE('',#3,-3.,.F.);
ENDSEC;",
                &mut dec,
            )
            .unwrap();
            for (surface, scale, units, point) in [
                (3, "1", AngleUnit::Degrees, ["0", "5", "0"]),
                (4, "1", AngleUnit::Radians, ["-1", "0", "0"]),
                (8, "10", AngleUnit::Degrees, ["30", "10", "30"]),
                (9, "10", AngleUnit::Degrees, ["30", "20", "30"]),
                (10, "1", AngleUnit::Radians, ["2", "0", "0"]),
            ] {
                let scale = q(scale);
                let point = point.map(q);
                let compiled = FamilySeeds::from_source(&source, surface, &scale, &units, &mut dec)
                    .unwrap()
                    .unwrap();
                let candidates = compiled.seeds(&point, &mut dec).unwrap();
                let charts = source
                    .surface_charts(surface, &scale, units, &mut dec)
                    .unwrap();
                let mut best: Option<Rational> = None;
                for uv in candidates {
                    let parameters = uv.each_ref().map(|value| Ball::from_rational(value, 128));
                    let sampled = charts[0].evaluate(&parameters, &mut dec).unwrap();
                    let error = sampled.iter().zip(&point).fold(
                        Rational::zero(),
                        |sum, (value, expected)| {
                            let difference = value - &Ball::from_rational(expected, 128);
                            sum + difference.lo().abs().max(difference.hi().abs())
                        },
                    );
                    best = Some(match best {
                        Some(old) => old.min(error),
                        None => error,
                    });
                }
                assert!(
                    best.is_some_and(|error| error < q("0.000000000000000001")),
                    "surface #{surface}"
                );
            }
        });
        assert!(!work.exhausted);
    }
}
