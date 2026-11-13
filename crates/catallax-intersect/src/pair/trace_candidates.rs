//! SMLib proposes selectors. Complete Catallax components remain the geometry.

use super::*;
use catallax_shapes::{BSplineSurface, Isometry, PatchPairCurve, PatchSectionBranch};
use catallax_trace::{
    spline::{BSplineSurface as CandidateSurface, ControlPoint, PointForm},
    surface_intersector::Intersector,
    vector::Vector3 as CandidateVector,
};

pub(super) struct Candidates {
    seeds: Vec<Point3>,
}

fn refusal(detail: impl Into<String>) -> Refusal {
    Refusal::new("intersect.trace-candidates", detail)
}

fn approximate(value: &Real) -> Result<f64, Refusal> {
    let value = value.to_f64();
    if !value.is_finite() {
        return Err(refusal(
            "the proposal input needs finite binary64 coordinates",
        ));
    }
    Ok(value)
}

fn surface(
    source: &BSplineSurface,
    map: &Isometry,
    dec: &mut Decider,
) -> Result<CandidateSurface, Refusal> {
    let mut points = Vec::with_capacity(source.controls().len());
    for (point, weight) in source.controls().iter().zip(source.weights()) {
        dec.check_deadline()?;
        let point = map.apply_point(point);
        let weight = approximate(weight)?;
        if weight <= 0.0 {
            return Err(refusal(
                "the unproved proposal adapter lost a positive spline weight",
            ));
        }
        let coordinates = [
            approximate(&point.x)? * weight,
            approximate(&point.y)? * weight,
            approximate(&point.z)? * weight,
            weight,
        ];
        if coordinates.iter().any(|value| !value.is_finite()) {
            return Err(refusal(
                "the homogeneous proposal control overflowed binary64",
            ));
        }
        points.push(ControlPoint { coordinates });
    }
    let knots = source.knots().map(|knots| {
        knots
            .iter()
            .map(Rational::to_f64_nearest)
            .collect::<Vec<_>>()
    });
    if knots.iter().flatten().any(|value| !value.is_finite()) {
        return Err(refusal(
            "the proposal knot domain is outside finite binary64",
        ));
    }
    CandidateSurface::new(source.degrees(), knots, source.control_count(), points)
        .map_err(|error| refusal(format!("unproved surface adaptation failed: {error:?}")))
}

impl Candidates {
    pub(super) fn generate(
        sources: [(&BSplineSurface, &Isometry); 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let first = surface(sources[0].0, sources[0].1, dec)?;
        let second = surface(sources[1].0, sources[1].1, dec)?;
        Self::trace(first, second, dec)
    }

    pub(super) fn plane_section(
        plane: &Plane,
        source: &BSplineSurface,
        swap: bool,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let vector = |value: &Vector3| -> Result<CandidateVector, Refusal> {
            Ok(CandidateVector {
                x: approximate(&value.x)?,
                y: approximate(&value.y)?,
                z: approximate(&value.z)?,
            })
        };
        let origin = CandidateVector {
            x: approximate(&plane.origin.x)?,
            y: approximate(&plane.origin.y)?,
            z: approximate(&plane.origin.z)?,
        };
        let x_axis = plane.normal.cross(&plane.ref_dir);
        let native = CandidateSurface::canonical_plane(
            origin,
            vector(&x_axis)?,
            vector(&plane.normal.cross(&x_axis))?,
        )
        .map_err(|error| refusal(format!("unproved plane adaptation failed: {error:?}")))?;
        let fitted = surface(source, &Isometry::identity(), dec)?;
        if swap {
            Self::trace(fitted, native, dec)
        } else {
            Self::trace(native, fitted, dec)
        }
    }

    fn trace(
        first: CandidateSurface,
        second: CandidateSurface,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        let curves = Intersector::new(
            &first,
            first.natural_domain(),
            &second,
            second.natural_domain(),
        )
        .do_intersection()
        .map_err(|error| refusal(format!("unproved branch tracing failed: {error:?}")))?;
        dec.check_deadline()?;
        let mut seeds = Vec::with_capacity(curves.len());
        for curve in curves {
            dec.check_deadline()?;
            let (point, _) = curve
                .curve
                .control_point(0, PointForm::NonRational)
                .map_err(|error| refusal(format!("proposal selector unavailable: {error:?}")))?;
            let exact = |value| {
                Rational::from_f64_exact(value)
                    .map(Real::from)
                    .ok_or_else(|| refusal("a proposal selector is not finite"))
            };
            seeds.push(Point3::new(
                exact(point.x)?,
                exact(point.y)?,
                exact(point.z)?,
            ));
        }
        Ok(Self { seeds })
    }

    /// Each proposal chooses an exact retained witness, not an accepted fit.
    /// The independent complete cover must be represented by these selectors.
    pub(super) fn prove_selection(
        &self,
        components: &[PatchPairCurve],
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        self.select(components.len(), |index| components[index].witnesses(), dec)
    }

    pub(super) fn prove_sections(
        &self,
        components: &[PatchSectionBranch],
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        self.select(components.len(), |index| components[index].witnesses(), dec)
    }

    fn select<F, I>(&self, count: usize, witnesses: F, dec: &mut Decider) -> Result<(), Refusal>
    where
        F: Fn(usize) -> I,
        I: Iterator<Item = Point3>,
    {
        let mut selected = vec![false; count];
        for seed in &self.seeds {
            dec.check_deadline()?;
            let mut best: Option<(usize, Real)> = None;
            let mut tied = false;
            for index in 0..count {
                for witness in witnesses(index) {
                    let distance = witness.minus(seed).norm_squared();
                    let Some((previous, previous_distance)) = &best else {
                        best = Some((index, distance));
                        continue;
                    };
                    match sign(&(&distance - previous_distance), dec)? {
                        Sign::Negative => {
                            best = Some((index, distance));
                            tied = false;
                        }
                        Sign::Zero if *previous != index => tied = true,
                        _ => (),
                    }
                }
            }
            let Some((index, _)) = best else {
                return Err(refusal(
                    "a traced proposal has no certified component selector",
                ));
            };
            if tied {
                return Err(refusal(
                    "a proposal has no strict exact retained-witness selector",
                ));
            }
            selected[index] = true;
        }
        if selected.iter().any(|selected| !selected) {
            return Err(Refusal::new(
                "intersect.trace-candidate-cover",
                "the proposals do not select every component of the exact complete cover",
            ));
        }
        Ok(())
    }
}
