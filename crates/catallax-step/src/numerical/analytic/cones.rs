//! STEP cone v is axial height; OCCT converts it to generator distance on import.

use super::*;
use catallax_shapes::Cone;

pub(super) struct Parameters {
    sine: Ball,
    cosine: Ball,
    pub angle_rounding_error: Rational,
}

pub(super) fn rounded(
    w: &mut Writer<'_>,
    name: &str,
    source: &Cone,
) -> Result<(Id, Surface, Parameters), Refusal> {
    let slope = w.dec.enclose_budgeted_at(&source.slope, BITS)?;
    let angle = Ball::atan2(&slope, &Ball::from_rational(&Rational::one(), BITS))
        .ok_or_else(|| bad("rounding-topology", "cone semi-angle is unresolved"))?;
    let original = angle.mid();
    // Angular rounding has its own record. Geometry checks enforce the physical tolerance.
    let number = crate::number::rational(&original, w.options.significant_digits);
    let written_angle = written(&original, w.options)?;
    let angle = Ball::from_rational(&written_angle, BITS);
    let sine = angle.sin();
    let cosine = angle.cos();
    if angle.lo().sign() != Sign::Positive
        || cosine.lo().sign() != Sign::Positive
        || angle.hi()
            >= (Ball::pi(BITS) * Ball::from_rational(&Rational::from_ratio(1, 2).unwrap(), BITS))
                .lo()
    {
        return Err(bad(
            "rounding-topology",
            "written cone semi-angle must remain acute",
        ));
    }
    let slope = sine
        .checked_div(&cosine)
        .ok_or_else(|| bad("rounding-topology", "written cone slope is unresolved"))?
        .mid();
    let WrittenFrame {
        entity,
        origin,
        normal,
        reference,
    } = placement(w, &source.apex, &source.axis, &source.ref_dir)?;
    let id = w.graph.add(
        "CONICAL_SURFACE",
        vec![
            text(string(name)),
            entity.into(),
            text("0."),
            text(number.text),
        ],
    );
    Ok((
        id,
        Surface::Cone(Cone {
            apex: point_from(&origin),
            axis: vector_from(&normal),
            slope: Real::from(slope),
            ref_dir: vector_from(&reference),
        }),
        Parameters {
            sine,
            cosine,
            angle_rounding_error: number.error,
        },
    ))
}

pub(super) fn uv(
    cone: &Cone,
    p: &[Rational; 3],
    previous: &[Rational; 2],
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    let [x, y, z] = local(&cone.apex, &cone.axis, &cone.ref_dir, p)?;
    Ok([angle(&x, &y, &previous[0], dec)?, z.mid()])
}

impl Parameters {
    /// Evaluate the actual decimal STEP cone, including its trigonometric enclosure.
    pub fn evaluate(
        &self,
        surface: &Surface,
        u: &Ball,
        v: &Ball,
        dec: &mut Decider,
    ) -> Result<[Ball; 3], Refusal> {
        let Surface::Cone(cone) = surface else {
            return Err(bad(
                "rounding-topology",
                "cone parameter record has the wrong carrier",
            ));
        };
        let origin = rational_point(&cone.apex)?;
        let axis = rational_vector(&cone.axis)?;
        let reference = rational_vector(&cone.ref_dir)?;
        let second = charged_cross(&axis, &reference, dec)?;
        let mut basis = Vec::new();
        for direction in [&axis, &reference, &second] {
            let norm = Ball::from_rational(&dot(direction, direction, dec)?, BITS)
                .sqrt()
                .ok_or_else(|| bad("rounding-topology", "written cone frame is degenerate"))?;
            let mut unit = Vec::new();
            for component in direction {
                unit.push(
                    Ball::from_rational(component, BITS)
                        .checked_div(&norm)
                        .ok_or_else(|| {
                            bad("rounding-topology", "written cone frame is unresolved")
                        })?,
                );
            }
            basis.push(unit);
        }
        let slope = self
            .sine
            .checked_div(&self.cosine)
            .ok_or_else(|| bad("rounding-topology", "written cone slope is unresolved"))?;
        let radius = v * &slope;
        let (cosine, sine) = (u.cos(), u.sin());
        Ok(std::array::from_fn(|i| {
            Ball::from_rational(&origin[i], BITS)
                + &basis[0][i] * v
                + (&basis[1][i] * &cosine + &basis[2][i] * &sine) * radius.clone()
        }))
    }
}
