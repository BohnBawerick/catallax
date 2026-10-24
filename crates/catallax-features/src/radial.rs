//! Exact distance and normal-sweep laws for a complete spherical or ring-torus skin.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{Origin, Role, StepId, ThickenSide};
use catallax_shapes::Surface;

use crate::profile::sign;

/// A regularised radial band. Its inner support, when present, bounds a cavity.
#[derive(Clone, Debug)]
pub struct RadialBand {
    /// Outward support and all its generating origins.
    pub outer: (Surface, Vec<Origin>),
    /// Cavity support with its outward-from-cavity orientation.
    pub inner: Option<(Surface, Vec<Origin>)>,
}

fn radius(surface: &Surface) -> Result<&Real, Refusal> {
    match surface {
        Surface::Sphere(s) => Ok(&s.radius),
        Surface::Torus(s) => Ok(&s.minor),
        _ => Err(Refusal::new("distance.radial-domain", "a complete sphere or ring torus is required")),
    }
}

impl RadialBand {
    /// Keeps exactly the material no farther than `thickness` from the closed skin.
    pub fn shell(surface: &Surface, origins: &[Origin], thickness: &Real, step: &StepId, dec: &mut Decider) -> Result<Self, Refusal> {
        radius(surface)?;
        let inner = surface.parallel_support(&(-thickness), dec)?.map(|surface| (surface,
            origins.iter().map(|of| Origin { step: step.clone(), role: Role::Inner { of: Box::new(of.clone()) } }).collect()));
        Ok(Self { outer: (surface.clone(), origins.to_vec()), inner })
    }

    /// Takes the full normal image, including a back sweep folding through its core.
    pub fn thicken(surface: &Surface, origins: &[Origin], thickness: &Real, side: ThickenSide, step: &StepId, dec: &mut Decider) -> Result<Self, Refusal> {
        let r = radius(surface)?;
        let far_origins: Vec<_> = origins.iter().map(|of| Origin {
            step: step.clone(), role: Role::Offset { of: Box::new(of.clone()) },
        }).collect();
        if matches!(side, ThickenSide::Front) {
            return Ok(Self {
                outer: (surface.parallel_support(thickness, dec)?.expect("positive radial growth"), far_origins),
                inner: Some((surface.clone(), origins.to_vec())),
            });
        }
        let twice = r * &Real::from_i64(2);
        match sign(&(thickness - &twice), dec)? {
            Sign::Positive => Ok(Self {
                outer: (surface.parallel_support(&(thickness - &twice), dec)?.expect("positive folded growth"), far_origins),
                inner: None,
            }),
            Sign::Zero => {
                let mut all = origins.to_vec();
                all.extend(far_origins);
                all.sort();
                all.dedup();
                Ok(Self { outer: (surface.clone(), all), inner: None })
            }
            Sign::Negative => Ok(Self {
                outer: (surface.clone(), origins.to_vec()),
                inner: surface.parallel_support(&(-thickness), dec)?.map(|surface| (surface, far_origins)),
            }),
        }
    }
}
