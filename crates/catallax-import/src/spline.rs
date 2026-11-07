//! STEP spline encoding recognition. Geometric algebra lives in shapes;
//! this module only translates Part 21 storage, units and logical declarations.
use super::*;
use catallax_numbers::Refusal;

pub(super) fn source_error(error: Refusal) -> ImportRefusal {
    if error.code.starts_with("numbers.") {
        ImportRefusal::ResourceLimit
    } else {
        malformed(error.detail)
    }
}

fn integer(p: &Param) -> Result<usize, ImportRefusal> {
    match p {
        Param::Number(s) => s
            .parse()
            .map_err(|_| malformed("expected a nonnegative integer")),
        _ => Err(malformed("expected a nonnegative integer")),
    }
}

fn list(p: &Param) -> Result<&[Param], ImportRefusal> {
    match p {
        Param::List(v) => Ok(v),
        _ => Err(malformed("expected a list")),
    }
}

fn logical(p: &Param) -> Result<Option<bool>, ImportRefusal> {
    match p {
        Param::Enum(e) if e == "T" => Ok(Some(true)),
        Param::Enum(e) if e == "F" => Ok(Some(false)),
        Param::Enum(e) if e == "U" => Ok(None),
        _ => Err(malformed("expected a STEP logical")),
    }
}

fn arg(args: &[Param], i: usize) -> Result<&Param, ImportRefusal> {
    args.get(i)
        .ok_or_else(|| malformed("missing spline argument"))
}

fn expanded(
    mult: &Param,
    values: &Param,
    count: usize,
    degree: usize,
    dec: &Decider,
) -> Result<Vec<Rational>, ImportRefusal> {
    let mult = list(mult)?;
    let values = list(values)?;
    let expected = count
        .checked_add(degree)
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| malformed("expanded knot count overflow"))?;
    if mult.len() != values.len() {
        return Err(malformed("knot and multiplicity counts differ"));
    }
    let mut knots = Vec::new();
    let mut previous = None;
    for (m, value) in mult.iter().zip(values) {
        dec.check_deadline().map_err(source_error)?;
        let m = integer(m)?;
        let value = number(value)?;
        if m == 0 || knots.len().checked_add(m).is_none_or(|n| n > expected) {
            return Err(malformed("invalid knot multiplicity or expanded count"));
        }
        if previous.as_ref().is_some_and(|p| p >= &value) {
            return Err(malformed("distinct knots must increase strictly"));
        }
        previous = Some(value.clone());
        for _ in 0..m {
            dec.check_deadline().map_err(source_error)?;
            knots.push(value.clone());
        }
    }
    if knots.len() != expected {
        return Err(malformed(
            "expanded knot count does not match degree and controls",
        ));
    }
    Ok(knots)
}

fn implicit(
    kind: &str,
    degree: usize,
    count: usize,
    dec: &Decider,
) -> Result<Vec<Rational>, ImportRefusal> {
    if degree == 0 || degree >= count {
        return Err(malformed("invalid implicit spline degree"));
    }
    let len = count
        .checked_add(degree)
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| malformed("knot count overflow"))?;
    let bezier = kind.starts_with("BEZIER");
    if bezier && count != degree + 1 {
        return Err(malformed("Bezier control count does not match degree"));
    }
    let uniform = kind.starts_with("UNIFORM");
    (0..len)
        .map(|i| {
            dec.check_deadline().map_err(source_error)?;
            let k = if uniform {
                i
            } else {
                i.saturating_sub(degree).min(count - degree)
            };
            Ok(Rational::from_integer((k as u64).into()))
        })
        .collect()
}

// Complex instances omit inherited attributes in each subtype. Simple
// instances include the representation-item name and all base attributes.
fn base<'a>(e: &'a Entity, root: &str, kinds: &[&str]) -> Option<(&'a [Param], &'a str, usize)> {
    if let [(kind, args)] = e.parts.as_slice()
        && kinds.contains(&kind.as_str())
    {
        return Some((args, kind, 1));
    }
    let args = e.part(root)?;
    let kind = e
        .parts
        .iter()
        .find(|(kind, _)| kinds.contains(&kind.as_str()))?
        .0
        .as_str();
    Some((args, kind, 0))
}

impl Reader<'_> {
    pub(super) fn spline_curve(
        &self,
        id: u64,
        e: &Entity,
    ) -> Result<Option<StepSplineCurve>, ImportRefusal> {
        let kinds = [
            "B_SPLINE_CURVE_WITH_KNOTS",
            "BEZIER_CURVE",
            "UNIFORM_CURVE",
            "QUASI_UNIFORM_CURVE",
        ];
        let Some((a, kind, offset)) = base(e, "B_SPLINE_CURVE", &kinds) else {
            return Ok(None);
        };
        let degree = integer(arg(a, offset)?)?;
        let controls = list(arg(a, offset + 1)?)?
            .iter()
            .map(|p| {
                let Param::Ref(id) = p else {
                    return Err(malformed("expected spline control reference"));
                };
                self.point(*id)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let closed = logical(arg(a, offset + 3)?)?;
        let self_intersect = logical(arg(a, offset + 4)?)?;
        let knots = if kind == kinds[0] {
            let (args, start) = if offset == 1 {
                (a, 6)
            } else {
                (e.part(kind).ok_or_else(|| unsupported(id, e))?, 0)
            };
            expanded(
                arg(args, start)?,
                arg(args, start + 1)?,
                controls.len(),
                degree,
                self.dec,
            )?
        } else {
            implicit(kind, degree, controls.len(), self.dec)?
        };
        let weights = match e.part("RATIONAL_B_SPLINE_CURVE") {
            Some(args) => list(arg(args, 0)?)?
                .iter()
                .map(number)
                .collect::<Result<_, _>>()?,
            None => vec![Rational::one(); controls.len()],
        };
        let geometry = RationalSplineCurve::new(degree, knots, controls, weights, self.dec)
            .map_err(source_error)?;
        let geometry = geometry.into_clamped(self.dec).map_err(source_error)?;
        Ok(Some(StepSplineCurve {
            geometry,
            closed,
            self_intersect,
        }))
    }

    pub(super) fn spline_surface(
        &self,
        id: u64,
        e: &Entity,
    ) -> Result<Option<StepSplineSurface>, ImportRefusal> {
        let kinds = [
            "B_SPLINE_SURFACE_WITH_KNOTS",
            "BEZIER_SURFACE",
            "UNIFORM_SURFACE",
            "QUASI_UNIFORM_SURFACE",
        ];
        let Some((a, kind, offset)) = base(e, "B_SPLINE_SURFACE", &kinds) else {
            return Ok(None);
        };
        let degrees = [integer(arg(a, offset)?)?, integer(arg(a, offset + 1)?)?];
        let rows = list(arg(a, offset + 2)?)?;
        let mut count = [rows.len(), 0];
        let mut controls = Vec::new();
        for (index, row) in rows.iter().enumerate() {
            let row = list(row)?;
            if index == 0 {
                count[1] = row.len();
            }
            if row.len() != count[1] {
                return Err(malformed("ragged spline control net"));
            }
            for p in row {
                let Param::Ref(id) = p else {
                    return Err(malformed("expected surface control reference"));
                };
                controls.push(self.point(*id)?);
            }
        }
        let closed = [logical(arg(a, offset + 4)?)?, logical(arg(a, offset + 5)?)?];
        let self_intersect = logical(arg(a, offset + 6)?)?;
        let knots = if kind == kinds[0] {
            let (args, start) = if offset == 1 {
                (a, 8)
            } else {
                (e.part(kind).ok_or_else(|| unsupported(id, e))?, 0)
            };
            [
                expanded(
                    arg(args, start)?,
                    arg(args, start + 2)?,
                    count[0],
                    degrees[0],
                    self.dec,
                )?,
                expanded(
                    arg(args, start + 1)?,
                    arg(args, start + 3)?,
                    count[1],
                    degrees[1],
                    self.dec,
                )?,
            ]
        } else {
            [
                implicit(kind, degrees[0], count[0], self.dec)?,
                implicit(kind, degrees[1], count[1], self.dec)?,
            ]
        };
        let weights = match e.part("RATIONAL_B_SPLINE_SURFACE") {
            Some(args) => {
                let rows = list(arg(args, 0)?)?;
                if rows.len() != count[0] {
                    return Err(malformed("weight row count does not match surface"));
                }
                let mut weights = Vec::with_capacity(controls.len());
                for row in rows {
                    let row = list(row)?;
                    if row.len() != count[1] {
                        return Err(malformed("weight column count does not match surface"));
                    }
                    for p in row {
                        weights.push(number(p)?);
                    }
                }
                weights
            }
            None => vec![Rational::one(); controls.len()],
        };
        let geometry =
            RationalSplineSurface::new(degrees, knots, count, controls, weights, self.dec)
                .map_err(source_error)?;
        let geometry = geometry.into_clamped(self.dec).map_err(source_error)?;
        Ok(Some(StepSplineSurface {
            geometry,
            closed,
            self_intersect,
        }))
    }
}
