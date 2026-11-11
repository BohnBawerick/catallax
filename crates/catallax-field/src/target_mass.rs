//! Uniform-density mass and centre enclosures on the complete F7 definition.
//!
//! This opt-in query does not change `mass::properties` or require inertia.
//! Unresolved grid cells retain their entire possible first moment. No midpoint
//! is presented as a measurement. Density is explicit and uniform in world mm.
use crate::solver::{Definition, GridOptions, Materials};
use catallax_measure::Enclosure;
use catallax_numbers::{Decider, Rational, Refusal};
use catallax_recipe::Recipe;

/// Complete certified record, including unresolved numerical width.
#[derive(Clone, Debug, PartialEq)]
pub struct Mass {
    /// True only after independent unmarked B-rep measurements of this recipe
    /// have been intersected with the field enclosures.
    pub brep_agreed: bool,
    /// Identity of the complete recipe and homogeneous material assignment.
    pub revision: String,
    /// Explicit uniform density, grams per cubic millimetre.
    pub density: Rational,
    /// Volume in cubic millimetres.
    pub volume: Enclosure,
    /// Mass in grams.
    pub grams: Enclosure,
    /// Unit-density first moments in world coordinates, mm^4.
    pub first_moments: [Enclosure; 3],
    /// World-coordinate centre, mm; None only for certified zero volume.
    pub centre: Option<[Enclosure; 3]>,
    /// Actual classifier work in the canonical parallel jobs.
    pub classifier_work: u64,
    /// Work reserved for those jobs, independent of worker count.
    pub reserved_work: u64,
}
/// Execution metadata is separate from semantic equality.
pub struct Run {
    /// Complete answer or original typed refusal, never a partial sum.
    pub answer: Result<Mass, Refusal>,
    /// Number of nonempty spawned workers (zero if preparation refused).
    pub spawned_workers: usize,
    /// Actual classifier work, including failed jobs.
    pub classifier_work: u64,
}
fn zero() -> Enclosure {
    Enclosure {
        lo: Rational::zero(),
        hi: Rational::zero(),
    }
}
fn product(a: &Enclosure, b: &Enclosure) -> Enclosure {
    let p = [&a.lo * &b.lo, &a.lo * &b.hi, &a.hi * &b.lo, &a.hi * &b.hi];
    Enclosure {
        lo: p.iter().min().unwrap().clone(),
        hi: p.iter().max().unwrap().clone(),
    }
}
/// Measure a complete native recipe with an explicitly assigned uniform material.
/// Requires the caller's existing operation scope. Grid resolution controls the
/// returned uncertainty, not a modelling tolerance. No B-rep fallback or material
/// density inference occurs. Caller cancellation/deadlines propagate via Decider
/// and the grid's transferred allowance; wall-clock cancellation is not a numeric
/// determinism claim. All available cores are used by default.
pub fn measure(
    recipe: &Recipe,
    density: &Rational,
    options: &GridOptions,
    dec: &mut Decider,
) -> Run {
    let mut spawned_workers = 0;
    let mut classifier_work = 0;
    let answer = (|| {
        if density <= &Rational::zero() {
            return Err(Refusal::new(
                "target-mass.density",
                "uniform density must be positive in g/mm^3",
            ));
        }
        let definition = Definition::compile(
            recipe,
            Materials::Homogeneous {
                region: "body".into(),
                material: "explicit-uniform-density".into(),
            },
            dec,
        )?;
        let queries = definition.prepare(dec)?;
        let run = queries.grid(options, dec);
        spawned_workers = run.spawned_workers;
        classifier_work = run.actual_work;
        let grid = run.answer?;
        let mut moments = std::array::from_fn::<_, 3, _>(|_| zero());
        // Canonical sequential reduction. Each cell is charged before rational
        // arithmetic; this is additive to, not multiplied by, worker allowances.
        for (index, cell) in grid.cells.iter().enumerate() {
            dec.check_deadline()?;
            dec.charge_work(64)?;
            let bounds = grid.cell_bounds(index).expect("complete grid index");
            let volume = bounds
                .iter()
                .fold(Rational::from(1), |v, b| v * (&b[1] - &b[0]));
            let occupied = Enclosure {
                lo: &volume * &cell.fraction.lo,
                hi: &volume * &cell.fraction.hi,
            };
            for axis in 0..3 {
                let contribution = if cell.fraction.lo == Rational::from(1) {
                    // Full measure: the exact box first moment, even if the
                    // closed box includes a measure-zero boundary.
                    let value = &volume * &(&bounds[axis][0] + &bounds[axis][1]).mul_pow2(-1);
                    Enclosure {
                        lo: value.clone(),
                        hi: value,
                    }
                } else {
                    product(
                        &occupied,
                        &Enclosure {
                            lo: bounds[axis][0].clone(),
                            hi: bounds[axis][1].clone(),
                        },
                    )
                };
                moments[axis].lo = &moments[axis].lo + &contribution.lo;
                moments[axis].hi = &moments[axis].hi + &contribution.hi;
            }
        }
        dec.charge_work(32)?;
        let centre = if grid.volume.hi.is_zero() {
            None
        } else {
            if grid.volume.lo <= Rational::zero() {
                return Err(Refusal::new(
                    "target-mass.positive-volume",
                    "grid does not certify positive volume; centre remains unresolved",
                ));
            }
            let inverse = Enclosure {
                lo: Rational::from(1).checked_div(&grid.volume.hi).unwrap(),
                hi: Rational::from(1).checked_div(&grid.volume.lo).unwrap(),
            };
            Some(std::array::from_fn(|axis| {
                let c = product(&moments[axis], &inverse);
                // The barycentre of nonnegative mass lies in its bounding box.
                Enclosure {
                    lo: c.lo.max(grid.bounds[axis][0].clone()),
                    hi: c.hi.min(grid.bounds[axis][1].clone()),
                }
            }))
        };
        let grams = Enclosure {
            lo: &grid.volume.lo * density,
            hi: &grid.volume.hi * density,
        };
        Ok(Mass {
            brep_agreed: false,
            revision: grid.revision,
            density: density.clone(),
            volume: grid.volume,
            grams,
            first_moments: moments,
            centre,
            classifier_work,
            reserved_work: grid.reserved_work,
        })
    })();
    Run {
        answer,
        spawned_workers,
        classifier_work,
    }
}

/// One result from independent field and B-rep builds of the same complete recipe.
/// Both exact-definition routes must succeed. Marked B-rep measurements are not
/// silently treated as proof. Intersection of their certified enclosures is a
/// common enclosure, not evidence that two arbitrary overlapping estimates agree.
/// A disjoint result refuses; neither route's existing answer is modified.
pub fn measure_agreed(
    recipe: &Recipe,
    density: &Rational,
    options: &GridOptions,
    accuracy: &catallax_measure::Accuracy,
    dec: &mut Decider,
) -> Run {
    let mut run = measure(recipe, density, options, dec);
    run.answer = run.answer.and_then(|mut result| {
        let built = catallax_features::evaluate(recipe, dec).map_err(|error| match error {
            catallax_features::BuildError::Geometry { refusal, .. } => refusal,
            catallax_features::BuildError::Recipe(error) => {
                Refusal::new("target-mass.recipe", error.to_string())
            }
        })?;
        if !built.stand_ins.is_empty() || !built.model.body(built.result).marks.is_empty() {
            return Err(Refusal::new(
                "target-mass.brep-definition",
                "a stand-in or marked B-rep is not the same exact solid for mass agreement",
            ));
        }
        let volume = catallax_measure::volume(&built.model, built.result, accuracy, dec)?;
        let centre = catallax_measure::centroid(&built.model, built.result, accuracy, dec)?;
        if !volume.is_exact() || !centre.is_exact() {
            return Err(Refusal::new(
                "target-mass.brep-marked",
                "agreement requires unmarked B-rep measurements of the exact definition",
            ));
        }
        fn common(a: &Enclosure, b: &Enclosure) -> Result<Enclosure, Refusal> {
            let c = Enclosure {
                lo: a.lo.clone().max(b.lo.clone()),
                hi: a.hi.clone().min(b.hi.clone()),
            };
            if c.lo > c.hi {
                return Err(Refusal::new(
                    "target-mass.disagreement",
                    "independent exact-definition mass enclosures are disjoint",
                ));
            }
            Ok(c)
        }
        result.volume = common(&result.volume, &volume.value)?;
        result.centre = match (result.centre, centre.value) {
            (Some(a), Some(b)) => Some([
                common(&a[0], &b[0])?,
                common(&a[1], &b[1])?,
                common(&a[2], &b[2])?,
            ]),
            (None, None) => None,
            _ => {
                return Err(Refusal::new(
                    "target-mass.disagreement",
                    "independent definitions disagree on whether a centre exists",
                ));
            }
        };
        result.grams = Enclosure {
            lo: &result.volume.lo * density,
            hi: &result.volume.hi * density,
        };
        result.brep_agreed = true;
        Ok(result)
    });
    run
}
