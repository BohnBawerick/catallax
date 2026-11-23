//! Reconstruct one spherical part from the measured points of an isolated scan.

use catallax_import::scan::ScanReference;
use catallax_numbers::{
    Decider, Rational, Refusal, operation_budget_active, with_operation_budget,
};
use catallax_recipe::{Number, Operation, Pose, Recipe, Step, StepId, Vector3};

use super::linalg::{V3, add, dot, eigen, moments_checked, norm, scale, solve, sub};
use super::{DeviationSide, MAX_SAMPLES, MIN_SAMPLES, stats};
use crate::Plan;

/// A fitted sphere and its measured deviation in both directions.
pub struct SphereFit {
    /// The scan reference and the fitted sphere, as an ordinary recipe.
    pub recipe: Recipe,
    /// The validated B-rep of the fitted sphere.
    pub brep: catallax_features::Built,
    /// The same sphere built independently in the field.
    pub field: Plan,
    /// Distance from every measured point to the fitted sphere.
    pub scan_to_part: DeviationSide,
    /// Distance from 512 directions on the fitted sphere to the scan.
    pub part_to_scan: DeviationSide,
}

fn domain(detail: &'static str) -> Refusal {
    Refusal::new("reconstruction.sphere-domain", detail)
}

/// Fits a single isolated spherical part, including its unobserved surface.
/// All supplied vertices must fit within `max_deviation_mm` of the written
/// sphere. The reverse deviation includes unsampled or occluded regions.
/// Both reported deviations are sampled measurements.
pub fn fit_sphere_scan(
    step: &Operation,
    bytes: &[u8],
    max_deviation_mm: f64,
    dec: &mut Decider,
) -> Result<SphereFit, Refusal> {
    if !max_deviation_mm.is_finite() || max_deviation_mm <= 0.0 {
        return Err(domain(
            "a spherical fit needs a positive finite scan deviation",
        ));
    }
    if operation_budget_active() {
        return fit(step, bytes, max_deviation_mm, dec);
    }
    let (answer, work) = with_operation_budget(dec.policy().max_work, || {
        fit(step, bytes, max_deviation_mm, dec)
    });
    if work.exhausted {
        return Err(Refusal::new(
            "numbers.operation-budget",
            "spherical scan fitting exhausted the caller allowance",
        ));
    }
    answer
}

fn estimate(samples: &[V3], dec: &mut Decider) -> Result<(V3, f64), Refusal> {
    let (mean, covariance, count) = moments_checked(samples.iter().copied(), dec)?;
    if !covariance.iter().flatten().all(|v| v.is_finite()) {
        return Err(domain(
            "the scan does not have a finite three-dimensional extent",
        ));
    }
    let ([smallest, _, largest], _) = eigen(covariance);
    if largest <= 0.0 || smallest <= largest * 0.01 {
        return Err(domain("the scanned directions do not constrain a sphere"));
    }

    // In coordinates relative to the sample mean, solve the linear sphere
    // equations for the centre. Subtracting the mean removes the constant term.
    let mut rhs = [0.0; 3];
    for &sample in samples {
        dec.charge_work(1)?;
        let at = sub(sample, mean);
        rhs = add(rhs, scale(at, dot(at, at) * (0.5 / count as f64)));
    }
    let offset = solve(covariance, rhs).ok_or_else(|| domain("the sphere fit is singular"))?;
    let center = add(mean, offset);
    let mut radius = 0.0;
    for &sample in samples {
        dec.charge_work(1)?;
        radius += norm(sub(sample, center)) / count as f64;
    }
    if !radius.is_finite() || radius <= 0.0 || smallest <= 0.025 * radius * radius {
        return Err(domain("the scan covers too little of the fitted sphere"));
    }
    Ok((center, radius))
}

fn fit(
    step: &Operation,
    bytes: &[u8],
    max_deviation: f64,
    dec: &mut Decider,
) -> Result<SphereFit, Refusal> {
    let scan = ScanReference::from_step(step, bytes, dec)?;
    let samples = scan.approx_points();
    if !(MIN_SAMPLES..=MAX_SAMPLES).contains(&samples.len()) {
        return Err(domain(
            "the scan sample count is outside the fitter's limits",
        ));
    }
    let (center, radius) = estimate(samples, dec)?;
    let number = |value| {
        Rational::from_f64_exact(value)
            .map(Number::new)
            .ok_or_else(|| domain("the fit is outside the finite numerical range"))
    };
    let [x, y, z] = center;
    let scan_id = StepId::new("scan").expect("static step id");
    let sphere_id = StepId::new("sphere").expect("static step id");
    let mut recipe = Recipe::new(
        "Sphere from scan",
        vec![
            Step {
                id: scan_id,
                op: step.clone(),
            },
            Step {
                id: sphere_id.clone(),
                op: Operation::Sphere {
                    radius: number(radius)?,
                    pose: Pose {
                        translation: Vector3([number(x)?, number(y)?, number(z)?]),
                        ..Pose::identity()
                    },
                },
            },
        ],
        sphere_id,
    );
    recipe.format = catallax_recipe::FORMAT_V4.into();

    let mut distances = Vec::with_capacity(samples.len());
    for &sample in samples {
        dec.charge_work(1)?;
        let distance = (norm(sub(sample, center)) - radius).abs();
        if !distance.is_finite() {
            return Err(domain("a measured point has nonfinite deviation"));
        }
        distances.push(distance);
    }
    let scan_to_part = stats(distances, dec)?;
    if scan_to_part.max_mm > max_deviation {
        return Err(Refusal::new(
            "reconstruction.sphere-deviation",
            format!(
                "maximum measured radial deviation {} mm exceeds {max_deviation} mm",
                scan_to_part.max_mm
            ),
        ));
    }

    let mut reverse = Vec::with_capacity(512);
    let angle = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    for i in 0..512 {
        dec.charge_work(1)?;
        let z = 1.0 - 2.0 * (f64::from(i) + 0.5) / 512.0;
        let radial = (1.0 - z * z).sqrt();
        let (sin, cos) = (angle * f64::from(i)).sin_cos();
        let sample = add(center, scale([radial * cos, radial * sin, z], radius));
        let distance = scan.surface().distance_checked(sample, dec)?;
        if !distance.is_finite() {
            return Err(domain("a sphere sample has nonfinite deviation"));
        }
        reverse.push(distance);
    }
    let part_to_scan = stats(reverse, dec)?;

    let brep = catallax_features::evaluate(&recipe, dec).map_err(|error| match error {
        catallax_features::BuildError::Geometry { refusal, .. } => refusal,
        error => Refusal::new("reconstruction.recipe", error.to_string()),
    })?;
    if !brep.model.check(brep.result).is_empty() {
        return Err(Refusal::new(
            "reconstruction.invalid-part",
            "the reconstructed sphere failed structural validation",
        ));
    }
    let field = Plan::compile_checked(&recipe, dec)?;
    dec.check_deadline()?;
    Ok(SphereFit {
        recipe,
        brep,
        field,
        scan_to_part,
        part_to_scan,
    })
}
