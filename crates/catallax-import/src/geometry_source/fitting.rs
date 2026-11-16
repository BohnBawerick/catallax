//! Finite numerical geometry fits with whole-chart interpolation bounds.
use super::*;
use catallax_numbers::{Ball, Rational, Refusal};
use catallax_shapes::spline_tools::{RationalSplineCurve, RationalSplineSurface};

/// Limits fixed before fitting one finite parameter rectangle or interval.
#[derive(Clone, Debug)]
pub struct FitOptions {
    /// Maximum source-to-fit distance, in millimetres.
    pub tolerance_mm: Rational,
    /// Maximum returned intervals or rectangles. Exhaustion returns no fit.
    pub maximum_pieces: usize,
}
/// A numerical secant fit. Its bound covers the complete source interval.
#[derive(Clone, Debug)]
pub struct FittedCurve {
    source: EntityId,
    wrappers: Vec<EntityId>,
    /// Exact control data of the numerical approximation, not the source curve.
    geometry: RationalSplineCurve,
    /// Source parameter interval covered by this piece.
    interval: [Rational; 2],
    /// Conservative maximum source-to-fit distance in millimetres.
    maximum_deviation_mm: Rational,
}
/// A numerical bilinear fit. Its bound covers the complete source rectangle.
#[derive(Clone, Debug)]
pub struct FittedSurface {
    source: EntityId,
    wrappers: Vec<EntityId>,
    /// Exact control data of the numerical approximation, not the source surface.
    geometry: RationalSplineSurface,
    /// Source parameter rectangle in u/v order.
    rectangle: [[Rational; 2]; 2],
    /// Conservative maximum source-to-fit distance in millimetres.
    maximum_deviation_mm: Rational,
}
impl FittedCurve {
    /// Original source curve instance.
    pub fn source(&self) -> EntityId {
        self.source
    }
    /// Retained source wrappers, including trim and composite-segment instances.
    pub fn wrappers(&self) -> &[EntityId] {
        &self.wrappers
    }
    /// Exact definition of the numerical approximation.
    pub fn geometry(&self) -> &RationalSplineCurve {
        &self.geometry
    }
    /// Covered source interval.
    pub fn interval(&self) -> &[Rational; 2] {
        &self.interval
    }
    /// Whole-interval source-to-fit distance bound in millimetres.
    pub fn maximum_deviation_mm(&self) -> &Rational {
        &self.maximum_deviation_mm
    }
}
impl FittedSurface {
    /// Original source surface instance.
    pub fn source(&self) -> EntityId {
        self.source
    }
    /// Retained source wrappers, including trims and surface patches.
    pub fn wrappers(&self) -> &[EntityId] {
        &self.wrappers
    }
    /// Exact definition of the numerical approximation.
    pub fn geometry(&self) -> &RationalSplineSurface {
        &self.geometry
    }
    /// Covered source rectangle.
    pub fn rectangle(&self) -> &[[Rational; 2]; 2] {
        &self.rectangle
    }
    /// Whole-rectangle source-to-fit distance bound in millimetres.
    pub fn maximum_deviation_mm(&self) -> &Rational {
        &self.maximum_deviation_mm
    }
}
fn invalid(detail: &str) -> MappingError {
    Refusal::new("import.geometry-fit", detail).into()
}
fn limit() -> MappingError {
    Refusal::new(
        "import.geometry-fit-limit",
        "finite geometry did not meet the fixed tolerance within the piece allowance",
    )
    .into()
}
fn validate(options: &FitOptions) -> Result<(), MappingError> {
    if options.tolerance_mm <= Rational::zero() || options.maximum_pieces == 0 {
        return Err(invalid(
            "fitting requires a positive tolerance and piece allowance",
        ));
    }
    Ok(())
}
fn magnitude(b: &Ball) -> Rational {
    b.lo().abs().max(b.hi().abs())
}
fn midpoint(a: &Rational, b: &Rational) -> Rational {
    (a + b) / Rational::from_i64(2)
}
fn sample(
    chart: &ParameterChart,
    at: &[Rational],
    bits: u32,
    dec: &mut Decider,
) -> Result<([Rational; 3], [Rational; 3]), MappingError> {
    let bounds = at
        .iter()
        .map(|p| [p.clone(), p.clone()])
        .collect::<Vec<_>>();
    let point = chart.evaluate_box(&bounds, bits, dec)?;
    let value = point.each_ref().map(|b| midpoint(&b.lo(), &b.hi()));
    let error = std::array::from_fn(|i| {
        (&value[i] - &point[i].lo())
            .abs()
            .max((&point[i].hi() - &value[i]).abs())
    });
    Ok((value, error))
}
fn splittable(error: &MappingError) -> bool {
    matches!(error,MappingError::Geometry(r) if r.code=="import.geometry-parameter-domain")
}
impl ParameterChart {
    /// Fit secants with a second-derivative remainder bound `M h² / 8` plus
    /// endpoint rounding. Zero attained error still returns a numerical fit.
    pub fn fit_curve(
        &self,
        interval: [Rational; 2],
        options: &FitOptions,
        dec: &mut Decider,
    ) -> Result<Vec<FittedCurve>, MappingError> {
        validate(options)?;
        if self.coordinate_space() != CoordinateSpace::Millimetres {
            return Err(invalid(
                "physical-distance fitting requires Cartesian millimetre coordinates",
            ));
        }
        if self.dimension() != 1 || interval[0] >= interval[1] {
            return Err(invalid(
                "curve fit requires an increasing one-dimensional interval",
            ));
        }
        let bits = dec.policy().max_precision_bits.min(128);
        if bits == 0 {
            return Err(invalid("no fitting arithmetic precision"));
        }
        let mut pending = vec![interval];
        let mut out = Vec::new();
        while let Some(cell) = pending.pop() {
            dec.charge_work(1)?;
            let result = (|| {
                let (a, ea) = sample(self, &[cell[0].clone()], bits, dec)?;
                let (b, eb) = sample(self, &[cell[1].clone()], bits, dec)?;
                let second =
                    self.evaluate_box_derivative(std::slice::from_ref(&cell), &[0, 0], bits, dec)?;
                let factor = (&cell[1] - &cell[0]).pow(2).expect("nonnegative exponent")
                    / Rational::from_i64(8);
                let bound = (0..3).fold(Rational::zero(), |sum, i| {
                    sum + (&magnitude(&second[i]) * &factor) + ea[i].clone().max(eb[i].clone())
                });
                Ok::<_, MappingError>((a, b, bound))
            })();
            match result {
                Ok((a, b, bound)) if bound <= options.tolerance_mm => {
                    let geometry = RationalSplineCurve::new(
                        1,
                        vec![
                            cell[0].clone(),
                            cell[0].clone(),
                            cell[1].clone(),
                            cell[1].clone(),
                        ],
                        vec![a, b],
                        vec![Rational::one(); 2],
                        dec,
                    )?;
                    out.push(FittedCurve {
                        source: self.source(),
                        wrappers: self.wrappers().to_vec(),
                        geometry,
                        interval: cell,
                        maximum_deviation_mm: bound,
                    });
                    continue;
                }
                Err(e) if !splittable(&e) => return Err(e),
                _ => {}
            }
            if out.len() + pending.len() + 2 > options.maximum_pieces {
                return Err(limit());
            }
            let middle = midpoint(&cell[0], &cell[1]);
            pending.push([middle.clone(), cell[1].clone()]);
            pending.push([cell[0].clone(), middle]);
        }
        Ok(out)
    }
    /// Fit bilinear rectangles. Tensor interpolation bounds the complete error
    /// by `M_uu hu² / 8 + M_vv hv² / 8` plus corner rounding in each coordinate.
    /// The L1 coordinate bound also bounds Euclidean distance.
    pub fn fit_surface(
        &self,
        rectangle: [[Rational; 2]; 2],
        options: &FitOptions,
        dec: &mut Decider,
    ) -> Result<Vec<FittedSurface>, MappingError> {
        validate(options)?;
        if self.coordinate_space() != CoordinateSpace::Millimetres {
            return Err(invalid(
                "physical-distance fitting requires Cartesian millimetre coordinates",
            ));
        }
        if self.dimension() != 2 || rectangle.iter().any(|r| r[0] >= r[1]) {
            return Err(invalid("surface fit requires increasing u/v intervals"));
        }
        let bits = dec.policy().max_precision_bits.min(128);
        if bits == 0 {
            return Err(invalid("no fitting arithmetic precision"));
        }
        let mut pending = vec![rectangle];
        let mut out = Vec::new();
        while let Some(cell) = pending.pop() {
            dec.charge_work(1)?;
            let widths = cell.each_ref().map(|[a, b]| b - a);
            let mut split = usize::from(widths[1] > widths[0]);
            let result = (|| {
                let mut controls = Vec::new();
                let mut rounding = [Rational::zero(), Rational::zero(), Rational::zero()];
                for u in &cell[0] {
                    for v in &cell[1] {
                        let (p, e) = sample(self, &[u.clone(), v.clone()], bits, dec)?;
                        controls.push(p);
                        for i in 0..3 {
                            rounding[i] = rounding[i].clone().max(e[i].clone());
                        }
                    }
                }
                let du = self.evaluate_box_derivative(&cell, &[0, 0], bits, dec)?;
                let dv = self.evaluate_box_derivative(&cell, &[1, 1], bits, dec)?;
                let factors = widths
                    .each_ref()
                    .map(|h| h.pow(2).expect("nonnegative exponent") / Rational::from_i64(8));
                let error = [du, dv]
                    .each_ref()
                    .iter()
                    .zip(&factors)
                    .map(|(d, f)| {
                        d.iter()
                            .fold(Rational::zero(), |sum, b| sum + (&magnitude(b) * f))
                    })
                    .collect::<Vec<_>>();
                split = usize::from(error[1] > error[0]);
                let bound = error[0].clone()
                    + error[1].clone()
                    + rounding.into_iter().fold(Rational::zero(), |a, b| a + b);
                Ok::<_, MappingError>((controls, bound))
            })();
            match result {
                Ok((controls, bound)) if bound <= options.tolerance_mm => {
                    let knots = cell
                        .each_ref()
                        .map(|[a, b]| vec![a.clone(), a.clone(), b.clone(), b.clone()]);
                    let geometry = RationalSplineSurface::new(
                        [1, 1],
                        knots,
                        [2, 2],
                        controls,
                        vec![Rational::one(); 4],
                        dec,
                    )?;
                    out.push(FittedSurface {
                        source: self.source(),
                        wrappers: self.wrappers().to_vec(),
                        geometry,
                        rectangle: cell,
                        maximum_deviation_mm: bound,
                    });
                    continue;
                }
                Err(e) if !splittable(&e) => return Err(e),
                _ => {}
            }
            if out.len() + pending.len() + 2 > options.maximum_pieces {
                return Err(limit());
            }
            let middle = midpoint(&cell[split][0], &cell[split][1]);
            let mut right = cell.clone();
            right[split][0] = middle.clone();
            let mut left = cell;
            left[split][1] = middle;
            pending.push(right);
            pending.push(left);
        }
        Ok(out)
    }
}
