//! Representation-local unit conversion, following STEP measure-with-unit chains.
use super::*;
use crate::part21::Param;
use catallax_numbers::Rational;

/// Geometric units from one context, without filling missing declarations.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct ContextUnits {
    /// Millimetres per written Cartesian length unit.
    pub length_mm: Option<Rational>,
    /// Plane-angle conversion, retaining the source's exact conversion measure.
    pub plane_angle: Option<AngleUnit>,
}
impl ParsedStep {
    /// Resolve only this context's length and angle assignments. Different
    /// representations may use different units in the same Part 21 file.
    pub fn context_units(
        &self,
        id: EntityId,
        dec: &mut Decider,
    ) -> Result<ContextUnits, ImportRefusal> {
        let context = self.representation_context(id)?;
        let mut out = ContextUnits {
            length_mm: None,
            plane_angle: None,
        };
        for unit in context.units {
            let e = self.entity(unit)?;
            decode::charge(e, dec)?;
            if e.part("LENGTH_UNIT").is_some() {
                let factor = crate::recognize::unit_millimetres(&self.graph, unit, 0)?;
                assign(&mut out.length_mm, factor, unit)?;
            }
            if e.part("PLANE_ANGLE_UNIT").is_some() {
                let factor = angle(self, unit, &mut Vec::new(), dec)?;
                assign(&mut out.plane_angle, AngleUnit::ScaledRadians(factor), unit)?;
            }
        }
        Ok(out)
    }

    pub(crate) fn model_unit_contexts(
        &self,
        models: &[ModelDeclaration],
        dec: &mut Decider,
    ) -> Result<Option<Vec<EntityId>>, ImportRefusal> {
        let roots: std::collections::BTreeSet<_> =
            models.iter().map(|model| model.entity).collect();
        let mut covered = std::collections::BTreeSet::new();
        let mut contexts = std::collections::BTreeSet::new();
        for (&id, entity) in &self.graph {
            dec.charge_work(1)
                .map_err(|_| ImportRefusal::ResourceLimit)?;
            if !entity
                .parts
                .iter()
                .any(|(kind, _)| kind.ends_with("REPRESENTATION"))
            {
                continue;
            }
            let Ok(representation) = self.representation(id) else {
                continue;
            };
            let mut used = false;
            for item in representation.items {
                if roots.contains(&item) {
                    covered.insert(item);
                    used = true;
                }
            }
            if used {
                if self
                    .representation_context(representation.context)?
                    .dimension
                    != Some(3)
                {
                    return Err(malformed(id, "body prototype requires a 3D unit context"));
                }
                contexts.insert(representation.context);
            }
        }
        // Unbound prototypes still require the former globally unambiguous
        // units; a representation-local choice cannot assign their geometry.
        if covered != roots {
            return Ok(None);
        }
        Ok(Some(contexts.into_iter().collect()))
    }
}
fn assign<T: Eq>(slot: &mut Option<T>, value: T, id: EntityId) -> Result<(), ImportRefusal> {
    if slot.as_ref().is_some_and(|old| *old != value) {
        return Err(malformed(
            id,
            "conflicting units in one representation context",
        ));
    }
    *slot = Some(value);
    Ok(())
}
fn angle(
    source: &ParsedStep,
    id: EntityId,
    active: &mut Vec<EntityId>,
    dec: &mut Decider,
) -> Result<Rational, ImportRefusal> {
    if active.len() >= 32 || active.contains(&id) {
        return Err(malformed(
            id,
            "cyclic or excessively deep plane-angle unit conversion",
        ));
    }
    active.push(id);
    let e = source.entity(id)?;
    decode::charge(e, dec)?;
    let factor = if let Some(a) = e.part("SI_UNIT") {
        if !matches!(a.get(1),Some(Param::Enum(name)) if name=="RADIAN") {
            return Err(malformed(id, "plane-angle SI unit is not RADIAN"));
        }
        let power = match a.first() {
            Some(Param::Unset) => 0,
            Some(Param::Enum(p)) => match p.as_str() {
                "EXA" => 18,
                "PETA" => 15,
                "TERA" => 12,
                "GIGA" => 9,
                "MEGA" => 6,
                "KILO" => 3,
                "HECTO" => 2,
                "DECA" => 1,
                "DECI" => -1,
                "CENTI" => -2,
                "MILLI" => -3,
                "MICRO" => -6,
                "NANO" => -9,
                "PICO" => -12,
                "FEMTO" => -15,
                "ATTO" => -18,
                _ => return Err(malformed(id, "unknown SI prefix")),
            },
            _ => return Err(malformed(id, "invalid SI prefix")),
        };
        Rational::from_i64(10)
            .pow(power)
            .ok_or_else(|| malformed(id, "invalid SI prefix exponent"))?
    } else if let Some(a) = e.part("CONVERSION_BASED_UNIT") {
        let measure = decode::reference(a, 1, id)?;
        let m = source.entity(measure)?;
        decode::charge(m, dec)?;
        let args = m
            .part("MEASURE_WITH_UNIT")
            .or_else(|| m.part("PLANE_ANGLE_MEASURE_WITH_UNIT"))
            .ok_or_else(|| malformed(measure, "expected plane-angle measure with unit"))?;
        let value = match args.first() {
            Some(Param::Typed(kind, values))
                if kind == "PLANE_ANGLE_MEASURE" && values.len() == 1 =>
            {
                decode::number(&values[0], measure)?
            }
            _ => return Err(malformed(measure, "expected PLANE_ANGLE_MEASURE")),
        };
        if value <= Rational::zero() {
            return Err(malformed(
                measure,
                "plane-angle conversion must be positive",
            ));
        }
        &value * &angle(source, decode::reference(args, 1, measure)?, active, dec)?
    } else {
        return Err(malformed(id, "unsupported plane-angle unit declaration"));
    };
    active.pop();
    Ok(factor)
}
