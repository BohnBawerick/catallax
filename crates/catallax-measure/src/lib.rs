//! Piece 6 of catallax: measurement.
//!
//! Volume, area and centre of mass of a body, and the area of a face, as
//! proven enclosures computed from exact definitions.
//!
//! The contract is `docs/pieces/06-measurement.md`.
//!
//! Planar rational boundaries use exact polynomial integration. Closed circles
//! and polynomial forms on rationally framed line, circle and ellipse arcs use
//! exact antiderivatives. Other arcs use piece 4's certified traces and degree-19
//! Taylor quadrature with a degree-20 whole-panel remainder.
//! Straight cone generators ending at the apex have identically zero boundary
//! integrands; their exact zero avoids a removable radial-frame singularity.
//!
//! Every face uses the divergence fields `p / 3` for volume and
//! `p_i^2 e_i / 2` for first moments, in one translated frame per body.
//! Spherical harmonics reduce sphere moments to polynomial boundary integrals;
//! an area gauge supplies the pole term. Torus integration lifts the tube angle
//! around each loop and includes the noncontractible-loop correction.
//! Whole-surface lift bounds retain both the empty and full endpoints, including
//! a closed surface represented with coincident seam boundaries.
//!
//! Ball operations share a fresh deterministic budget for each property.
//! Subdivision and precision refine independently. The last certified enclosure
//! survives budget exhaustion, with one width mark for each unmet property.
//! A centroid is returned only after its volume enclosure excludes zero.
//!
//! `cargo bench -p catallax-measure --bench speed` reports wall time, enclosure widths
//! and public Decider statistics. `CATALLAX_MEASURE_BENCH_QUICK=1` selects default
//! accuracy; `CATALLAX_MEASURE_BENCH_CASE=cylinder` selects a case by name.

use catallax_numbers::{
    Answer, Answered, Ball, BoundKind, Decider, Mark, MarkLog, Rational, Refusal, Sign,
};
use catallax_topology::{BodyId, BodyKind, FaceId, Model};

mod engine;
mod integrate;
mod planar;
mod ring;
mod surface;
mod torus;

/// A proven enclosure `[lo, hi]` with exact rational ends.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Enclosure {
    /// Lower end.
    pub lo: Rational,
    /// Upper end.
    pub hi: Rational,
}

impl Enclosure {
    /// `hi - lo`.
    pub fn width(&self) -> Rational {
        &self.hi - &self.lo
    }

    /// Whether `value` lies in the enclosure.
    pub fn contains(&self, value: &Rational) -> bool {
        self.lo <= *value && *value <= self.hi
    }

    /// Whether `lo == hi`.
    pub fn is_exact(&self) -> bool {
        self.lo == self.hi
    }

    fn exact(value: Rational) -> Self {
        Self {
            lo: value.clone(),
            hi: value,
        }
    }

    fn from_ball(value: &Ball) -> Self {
        Self {
            lo: value.lo(),
            hi: value.hi(),
        }
    }
}

/// Requested width: at most `max(relative * |value|, absolute)`.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Accuracy {
    /// Relative part.
    pub relative: Rational,
    /// Absolute part, in model units of the quantity.
    pub absolute: Rational,
}

impl Default for Accuracy {
    fn default() -> Self {
        Self {
            relative: Rational::parse("1e-12").unwrap(),
            absolute: Rational::parse("1e-30").unwrap(),
        }
    }
}

impl Accuracy {
    fn validate(&self) {
        assert!(
            self.relative.sign() != Sign::Negative && self.absolute.sign() != Sign::Negative,
            "negative measurement accuracy"
        );
    }

    fn reached(&self, enclosure: &Enclosure) -> bool {
        let magnitude = if enclosure.lo.sign() == enclosure.hi.sign() {
            enclosure.lo.abs().min(enclosure.hi.abs())
        } else {
            Rational::zero()
        };
        enclosure.width() <= (&self.relative * &magnitude).max(self.absolute.clone())
    }

    fn mark(&self, marks: &mut MarkLog, values: &[&Enclosure], property: &str) {
        if values.iter().any(|value| !self.reached(value)) {
            let width = values.iter().map(|value| value.width()).max().unwrap();
            marks.push(Mark::new(
                "measure.accuracy-not-reached",
                BoundKind::Value,
                width,
                format!("{property} width exceeds the requested accuracy"),
            ));
        }
    }
}

/// Mass properties of a body.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct MassProperties {
    /// Volume.
    pub volume: Enclosure,
    /// Surface area.
    pub area: Enclosure,
    /// Centre of mass; `None` exactly when the body is `Empty`.
    pub centroid: Option<[Enclosure; 3]>,
    /// Unit-density inertia about the centroid, in model axes. `None` for `Empty`.
    pub inertia: Option<InertiaTensor>,
}

/// The three original properties, without requesting second-moment integration.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct VolumeAreaCentroid {
    /// Volume.
    pub volume: Enclosure,
    /// Surface area.
    pub area: Enclosure,
    /// Centre of mass; `None` exactly for an empty body.
    pub centroid: Option<[Enclosure; 3]>,
}

/// Volume, area, centre of mass and centroidal unit-density inertia.
pub fn mass_properties(
    model: &Model,
    body: BodyId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<MassProperties> {
    let mut traces = engine::TraceCache::grouped();
    let base = base_properties(model, body, accuracy, &mut traces, dec)?;
    let inertia = inertia_from_centroid(
        model,
        body,
        accuracy,
        &base.value.volume,
        base.value.centroid.as_ref(),
        &mut traces,
        dec,
    )?;
    let mut marks = base.marks;
    marks.extend(&inertia.marks);
    Ok(Answer::with_marks(
        MassProperties {
            volume: base.value.volume,
            area: base.value.area,
            centroid: base.value.centroid,
            inertia: inertia.value,
        },
        marks,
    ))
}

/// Requests volume, area and centroid without requesting inertia. Each property
/// retains its own budget; edge traces and failed sphere-gauge attempts are shared.
pub fn volume_area_centroid(
    model: &Model,
    body: BodyId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<VolumeAreaCentroid> {
    base_properties(
        model,
        body,
        accuracy,
        &mut engine::TraceCache::grouped(),
        dec,
    )
}

fn base_properties(
    model: &Model,
    body: BodyId,
    accuracy: &Accuracy,
    traces: &mut engine::TraceCache,
    dec: &mut Decider,
) -> Answered<VolumeAreaCentroid> {
    accuracy.validate();
    let volume = scalar(
        model,
        body,
        accuracy,
        surface::Property::Volume,
        traces,
        dec,
    )?;
    let area = scalar(model, body, accuracy, surface::Property::Area, traces, dec)?;
    let mut marks = volume.marks;
    marks.extend(&area.marks);
    let centroid = if model.body(body).kind == BodyKind::Empty {
        None
    } else {
        if volume.value.contains(&Rational::zero()) {
            return Err(Refusal::new(
                "measure.centroid-not-separated",
                "public volume enclosure includes zero after exhausting its budget",
            ));
        }
        let faces = solid_faces(model, body)?;
        let anchor = engine::anchor(model, &faces);
        let answer = engine::measure(
            model,
            &faces,
            &anchor,
            accuracy,
            surface::Property::Centroid,
            traces,
            dec,
        )?;
        let values: [Enclosure; 3] = answer.value.try_into().unwrap();
        accuracy.mark(&mut marks, &values.each_ref(), "centroid");
        Some(values)
    };
    Ok(Answer::with_marks(
        VolumeAreaCentroid {
            volume: volume.value,
            area: area.value,
            centroid,
        },
        marks,
    ))
}

/// Symmetric unit-density inertia tensor, in model axes and model length units
/// to the fifth power. Off-diagonal entries are negative products of inertia.
pub type InertiaTensor = [[Enclosure; 3]; 3];

/// Centre of mass. Empty bodies have no centre.
pub fn centroid(
    model: &Model,
    body: BodyId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<Option<[Enclosure; 3]>> {
    accuracy.validate();
    let faces = solid_faces(model, body)?;
    if model.body(body).kind == BodyKind::Empty {
        return Ok(Answer::exact(None));
    }
    let anchor = engine::anchor(model, &faces);
    let answer = engine::measure(
        model,
        &faces,
        &anchor,
        accuracy,
        surface::Property::Centroid,
        &mut engine::TraceCache::new(),
        dec,
    )?;
    let values: [Enclosure; 3] = answer.value.try_into().unwrap();
    let mut marks = answer.marks;
    accuracy.mark(&mut marks, &values.each_ref(), "centroid");
    Ok(Answer::with_marks(Some(values), marks))
}

/// Proven unit-density inertia about the centroid. Measurement does not mutate
/// the body. Unsupported bodies retain the ordinary measurement refusal.
pub fn inertia(
    model: &Model,
    body: BodyId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<Option<InertiaTensor>> {
    let volume = volume(model, body, accuracy, dec)?;
    let center = centroid(model, body, accuracy, dec)?;
    inertia_from_centroid(
        model,
        body,
        accuracy,
        &volume.value,
        center.value.as_ref(),
        &mut engine::TraceCache::new(),
        dec,
    )
}

#[allow(clippy::too_many_arguments)]
fn inertia_from_centroid(
    model: &Model,
    body: BodyId,
    accuracy: &Accuracy,
    volume: &Enclosure,
    center: Option<&[Enclosure; 3]>,
    traces: &mut engine::TraceCache,
    dec: &mut Decider,
) -> Answered<Option<InertiaTensor>> {
    let faces = solid_faces(model, body)?;
    let Some(center) = center else {
        return Ok(Answer::exact(None));
    };
    // A centroid accurate enough for position need not be accurate enough for
    // inertia: its error is squared and multiplied by volume. Tighten only when
    // that parallel-axis uncertainty could exceed the tensor's absolute floor.
    let volume_magnitude = volume.lo.abs().max(volume.hi.abs());
    let center_width = center.iter().map(Enclosure::width).max().unwrap();
    let correction_bound = &volume_magnitude * &(&center_width * &center_width);
    let refined_center: Option<[Enclosure; 3]> = if !accuracy.absolute.is_zero()
        && correction_bound > &accuracy.absolute / &Rational::from_i64(16)
    {
        let tight = Accuracy {
            relative: Rational::zero(),
            absolute: (&accuracy.absolute
                / &((&volume_magnitude + &Rational::one()) * Rational::from_i64(16)))
                .min(Rational::one()),
        };
        let anchor = engine::anchor(model, &faces);
        match engine::measure(
            model,
            &faces,
            &anchor,
            &tight,
            surface::Property::Centroid,
            traces,
            dec,
        ) {
            Ok(answer) => Some(answer.value.try_into().unwrap()),
            // Keep the already certified centre if refinement cannot start.
            Err(refusal) if refusal.code == "measure.budget-exhausted" => None,
            Err(refusal) => return Err(refusal),
        }
    } else {
        None
    };
    let center = refined_center.as_ref().unwrap_or(center);
    let midpoint = center
        .each_ref()
        .map(|c| (&c.lo + &c.hi) / Rational::from_i64(2));
    let anchor = catallax_shapes::RPoint3::new(
        midpoint[0].clone(),
        midpoint[1].clone(),
        midpoint[2].clone(),
    );
    let offset: [Enclosure; 3] = std::array::from_fn(|i| Enclosure {
        lo: &center[i].lo - &midpoint[i],
        hi: &center[i].hi - &midpoint[i],
    });
    let raw_accuracy = Accuracy {
        relative: &accuracy.relative / &Rational::from_i64(4),
        absolute: &accuracy.absolute / &Rational::from_i64(4),
    };
    let mut moments: InertiaTensor =
        std::array::from_fn(|_| std::array::from_fn(|_| Enclosure::exact(Rational::zero())));
    for i in 0..3 {
        for j in i..3 {
            let answer = engine::measure(
                model,
                &faces,
                &anchor,
                &raw_accuracy,
                surface::Property::SecondMoment(i, j),
                traces,
                dec,
            )?;
            let correction = product(volume, &product(&offset[i], &offset[j]));
            moments[i][j] = difference(&answer.value[0], &correction);
            moments[j][i] = moments[i][j].clone();
        }
    }
    let tensor: InertiaTensor = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            if i == j {
                let a = &moments[(i + 1) % 3][(i + 1) % 3];
                let b = &moments[(i + 2) % 3][(i + 2) % 3];
                Enclosure {
                    lo: &a.lo + &b.lo,
                    hi: &a.hi + &b.hi,
                }
            } else {
                Enclosure {
                    lo: -&moments[i][j].hi,
                    hi: -&moments[i][j].lo,
                }
            }
        })
    });
    let mut marks = MarkLog::new();
    accuracy.mark(
        &mut marks,
        &[
            &tensor[0][0],
            &tensor[1][1],
            &tensor[2][2],
            &tensor[0][1],
            &tensor[0][2],
            &tensor[1][2],
        ],
        "inertia",
    );
    Ok(Answer::with_marks(Some(tensor), marks))
}

fn product(a: &Enclosure, b: &Enclosure) -> Enclosure {
    let values = [&a.lo * &b.lo, &a.lo * &b.hi, &a.hi * &b.lo, &a.hi * &b.hi];
    Enclosure {
        lo: values.iter().min().unwrap().clone(),
        hi: values.iter().max().unwrap().clone(),
    }
}

fn difference(a: &Enclosure, b: &Enclosure) -> Enclosure {
    Enclosure {
        lo: &a.lo - &b.hi,
        hi: &a.hi - &b.lo,
    }
}

/// Volume.
pub fn volume(
    model: &Model,
    body: BodyId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<Enclosure> {
    scalar(
        model,
        body,
        accuracy,
        surface::Property::Volume,
        &mut engine::TraceCache::new(),
        dec,
    )
}

/// Surface area.
pub fn area(
    model: &Model,
    body: BodyId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<Enclosure> {
    scalar(
        model,
        body,
        accuracy,
        surface::Property::Area,
        &mut engine::TraceCache::new(),
        dec,
    )
}

/// Area of one face.
pub fn face_area(
    model: &Model,
    face: FaceId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<Enclosure> {
    accuracy.validate();
    let faces = [face];
    let anchor = engine::anchor(model, &faces);
    let mut answer = engine::measure(
        model,
        &faces,
        &anchor,
        accuracy,
        surface::Property::Area,
        &mut engine::TraceCache::new(),
        dec,
    )?;
    let enclosure = answer.value.remove(0);
    accuracy.mark(&mut answer.marks, &[&enclosure], "face area");
    Ok(Answer::with_marks(enclosure, answer.marks))
}

fn solid_faces(model: &Model, body: BodyId) -> Result<Vec<FaceId>, Refusal> {
    match model.body(body).kind {
        BodyKind::Solid | BodyKind::Empty => Ok(model.faces_of_body(body)),
        _ => Err(Refusal::new(
            "measure.not-solid",
            "measurement requires a solid or empty body",
        )),
    }
}

fn scalar(
    model: &Model,
    body: BodyId,
    accuracy: &Accuracy,
    property: surface::Property,
    traces: &mut engine::TraceCache,
    dec: &mut Decider,
) -> Answered<Enclosure> {
    accuracy.validate();
    let faces = solid_faces(model, body)?;
    if model.body(body).kind == BodyKind::Empty {
        return Ok(Answer::exact(Enclosure::exact(Rational::zero())));
    }
    let anchor = engine::anchor(model, &faces);
    let mut answer = engine::measure(model, &faces, &anchor, accuracy, property, traces, dec)?;
    let enclosure = answer.value.remove(0);
    accuracy.mark(
        &mut answer.marks,
        &[&enclosure],
        if property == surface::Property::Area {
            "area"
        } else {
            "volume"
        },
    );
    Ok(Answer::with_marks(enclosure, answer.marks))
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod fixtures;
