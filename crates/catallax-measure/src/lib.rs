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
}

/// Volume, area and centre of mass.
pub fn mass_properties(
    model: &Model,
    body: BodyId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<MassProperties> {
    accuracy.validate();
    let mut traces = engine::TraceCache::new();
    let volume = scalar(
        model,
        body,
        accuracy,
        surface::Property::Volume,
        &mut traces,
        dec,
    )?;
    let area = scalar(
        model,
        body,
        accuracy,
        surface::Property::Area,
        &mut traces,
        dec,
    )?;
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
            &mut traces,
            dec,
        )?;
        let values: [Enclosure; 3] = answer.value.try_into().unwrap();
        accuracy.mark(&mut marks, &values.each_ref(), "centroid");
        Some(values)
    };
    Ok(Answer::with_marks(
        MassProperties {
            volume: volume.value,
            area: area.value,
            centroid,
        },
        marks,
    ))
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
