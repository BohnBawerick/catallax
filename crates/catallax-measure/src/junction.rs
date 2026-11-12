//! Volume, area and centroid of a blend with rolling-envelope, rolling-saddle
//! or sphere-band faces.
//!
//! The engine has no chart for a fitted carrier or for a cylinder trimmed by its
//! rail. The blend changes its unfilleted base only inside one ring about the
//! cylinder axis, so each property is the base's plus what the ring adds: the
//! exact section over the ring's untouched angle and one exact sector per mouth.
//! A saddle's change is one set of integrals round its ring, in `crate::saddle`,
//! and a sphere band's in `crate::saddle::sphere`.
//!
//! Volume and centroid enclose both the exact definition and the fitted carrier.
//! A homotopy from envelope to carrier moves each point at most U and keeps the
//! boundary rows on their neighbouring surfaces, so any point between the two
//! bodies lies within U of the exact envelope. A point of that tube reaches the
//! envelope along a normal or at its boundary. With principal curvatures 1/r and
//! k, the normal part is at most `2 U (1 + U / r) (A + U K)` for envelope area A
//! and K, the integral of |k| over it. The boundary part is at most `pi U^2 P`
//! for boundary length P: half-tubes around its curves and four corner balls,
//! given U far below every curve radius.
//!
//! A position bound says nothing about a carrier's area, so the area needs its
//! own argument. Where every fitted face is a rolling envelope that keeps its
//! whole chart, the represented area encloses both the definition's and the
//! carriers': each carrier's area is enclosed directly, and the faces it
//! borders change only by the strip its rail moves along the cylinder.
//! Elsewhere the area is the exact definition's, and every represented answer
//! that states it says so with an [`AREA_OF_DEFINITION`] mark, so no area
//! overstates what its enclosure covers.

use catallax_numbers::{
    Answer, Answered, Ball, BoundKind, Decider, Mark, MarkLog, Outcome, Rational, Real, Refusal,
};
use catallax_shapes::{
    BSplineSurface, ContactMoments, Curve, RollingEnvelope, RollingSaddle, RollingSphere,
    SplineSupport, Surface, Vector3,
};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, saddle_added};
use std::sync::{Condvar, Mutex, atomic::AtomicUsize};

use crate::{Accuracy, Enclosure, VolumeAreaCentroid, engine, product, surface};

const BITS: u32 = 128;
// Most cells for one carrier's area enclosure. A rolling envelope's area
// element is not separated from zero along its mouth row, so cells there are
// enclosed coarsely between zero and their largest element. That bound and the
// rail strips dominate the represented area on P5 whatever the cell grid: at
// 32 bits, one cell per carrier gives the same seven represented-area widths
// as sixteen (186, 226, 361, 44, 174, 6446 and 185 mm²; the requested accuracy
// is not reached at either, so the answer carries the accuracy mark with the
// width it did reach). One cell keeps the carrier stage at its floor: a
// hosted 2-core runner gives the writer one helper, so carriers are enclosed
// one at a time, and every extra cell there only adds work to the test's
// ten-minute limit without narrowing the answer.
const CARRIER_CELLS: usize = 1;
// The carrier areas are enclosed at 32 bits, not at the answer's [`BITS`]:
// the ball arithmetic's rounding there is relative 2^-32, about 2e-6 on a
// carrier of some thousand square millimetres, far below the enclosure's
// geometric width, and small mantissas make a stage that otherwise dominates
// the answer's time four times cheaper. [`crate::entity::gap`] bounds its
// boxes below the answer's precision for the same reason.
const CARRIER_BITS: u32 = 32;

/// What a measured enclosure has to contain.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Scope {
    /// The exact definition alone.
    Definition,
    /// The definition and the fitted carrier.
    Represented,
}

/// Code of the mark on a fitted blend's area: the enclosure contains the exact
/// definition's area only. Its bound is the worst fitted carrier's distance to
/// the definition; the carriers' own area is not bounded.
pub const AREA_OF_DEFINITION: &str = "measure.area-of-definition";

pub(crate) struct Junction<'a> {
    base: BodyId,
    envelopes: Vec<&'a RollingEnvelope>,
    // The face of each envelope, in the same order.
    envelope_faces: Vec<FaceId>,
    // Each with whether its pocket adds material.
    saddles: Vec<(&'a RollingSaddle, bool)>,
    spheres: Vec<(&'a RollingSphere, bool)>,
}

// What the ring adds, in world coordinates.
struct Ring {
    volume: Ball,
    area: Ball,
    moment: [Ball; 3],
    // Largest volume between the exact envelopes and their carriers, and the
    // first moment it can carry along each axis.
    tube: Rational,
    moment_slack: [Rational; 3],
    // Each envelope's exact definition area, in the junction's order.
    envelope_areas: Vec<Ball>,
}

fn exact(value: &Rational) -> Ball {
    Ball::from_rational(value, BITS)
}

fn zero_ball() -> Ball {
    exact(&Rational::zero())
}

fn ball(value: &Enclosure) -> Ball {
    exact(&value.lo).hull(&exact(&value.hi))
}

fn coordinates(v: &Vector3, dec: &Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        dec.enclose_at(&v.x, BITS)?,
        dec.enclose_at(&v.y, BITS)?,
        dec.enclose_at(&v.z, BITS)?,
    ])
}

/// Refuses unless `a` and `c` are proven equal. An undecided difference keeps
/// its decision.
fn same_ring(a: &Real, c: &Real, dec: &mut Decider) -> Result<(), Refusal> {
    match dec.sign(&(a - c)) {
        Outcome::ProvenCoincident(_) => Ok(()),
        Outcome::Refused(refusal) => Err(refusal),
        other => {
            Err(refused("rolling envelopes do not share one fillet ring").keeping_undecided(other))
        }
    }
}

fn refused(detail: &'static str) -> Refusal {
    Refusal::new("measure.unsupported-surface", detail)
}

fn fitted_rail(model: &Model, edge: EdgeId) -> bool {
    matches!(
        model.geom().curve(model.edge(edge).curve),
        Curve::BSpline(_)
    ) && model.coedges_of_edge(edge).into_iter().any(|coedge| {
        model.owner_of_coedge(coedge).is_some_and(|(_, face)| {
            matches!(
                model.geom().surface(model.face(face).surface),
                Surface::BSpline(
                    SplineSupport::RollingEnvelope(_)
                        | SplineSupport::RollingSaddle(_)
                        | SplineSupport::RollingSphere(_)
                )
            )
        })
    })
}

/// A fitted rail has no certified arc integral, so a face it trims refuses
/// before the engine spends its allowance on the trace.
pub(crate) fn refuse_fitted_face(model: &Model, face: FaceId) -> Result<(), Refusal> {
    let railed = model.face(face).loops.iter().any(|boundary| {
        matches!(model.loop_(*boundary), Loop::Coedges(coedges) if coedges
            .iter()
            .any(|coedge| fitted_rail(model, model.coedge(*coedge).edge)))
    });
    if railed {
        return Err(refused(
            "a face trimmed by a fitted rail has no single-face measure",
        ));
    }
    Ok(())
}

/// The same refusal for the rail itself.
pub(crate) fn refuse_fitted_edge(model: &Model, edge: EdgeId) -> Result<(), Refusal> {
    if fitted_rail(model, edge) {
        return Err(refused("a fitted rail has no certified arc integral"));
    }
    Ok(())
}

/// The body's junction, when any face is a rolling envelope. Such a body with
/// no recorded base keeps the typed refusal.
pub(crate) fn of(model: &Model, body: BodyId) -> Result<Option<Junction<'_>>, Refusal> {
    let (envelopes, envelope_faces): (Vec<&RollingEnvelope>, Vec<FaceId>) = model
        .faces_of_body(body)
        .into_iter()
        .filter_map(
            |face| match model.geom().surface(model.face(face).surface) {
                Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) => {
                    Some((envelope, face))
                }
                _ => None,
            },
        )
        .unzip();
    let mut saddles = Vec::new();
    let mut spheres = Vec::new();
    for face in model.faces_of_body(body) {
        match model.geom().surface(model.face(face).surface) {
            Surface::BSpline(SplineSupport::RollingSaddle(saddle)) => {
                saddles.push((saddle, saddle_added(model, face)?));
            }
            Surface::BSpline(SplineSupport::RollingSphere(band)) => {
                spheres.push((band, saddle_added(model, face)?));
            }
            _ => {}
        }
    }
    if envelopes.is_empty() && saddles.is_empty() && spheres.is_empty() {
        return Ok(None);
    }
    let base = model
        .blend_base(body)
        .ok_or_else(|| refused("this fitted blend has no recorded unfilleted base"))?;
    if model.blend_base(base).is_some() {
        return Err(refused("a fitted blend's base is itself a fitted blend"));
    }
    Ok(Some(Junction {
        base,
        envelopes,
        envelope_faces,
        saddles,
        spheres,
    }))
}

impl Junction<'_> {
    fn ring(&self, width: &Rational, moments: bool, dec: &mut Decider) -> Result<Ring, Refusal> {
        let mut ring = if self.envelopes.is_empty() {
            Ring {
                volume: zero_ball(),
                area: zero_ball(),
                moment: std::array::from_fn(|_| zero_ball()),
                tube: Rational::zero(),
                moment_slack: std::array::from_fn(|_| Rational::zero()),
                envelope_areas: Vec::new(),
            }
        } else {
            self.envelope_ring(width, dec)?
        };
        let mut changes: Vec<(catallax_shapes::SaddleLaw, crate::saddle::Change)> = Vec::new();
        for (saddle, added) in &self.saddles {
            let mut reused = None;
            for (source, change) in &changes {
                if let Some(change) = change.reflected(source, saddle, dec)? {
                    reused = Some(change);
                    break;
                }
            }
            let change = match reused {
                Some(change) => change,
                None => {
                    let change = crate::saddle::change(saddle, width, moments, dec)?;
                    if changes.len() < 32
                        && matches!(
                            saddle.definition().law().ring,
                            catallax_shapes::SaddleRing::Crossing { .. }
                        )
                    {
                        changes.push((saddle.definition().law(), change.clone()));
                    }
                    change
                }
            };
            let sign = exact(&Rational::from_i64(if *added { 1 } else { -1 }));
            ring.volume = &ring.volume + &(&change.pocket * &sign);
            ring.area = &ring.area + &change.area;
            let tube = saddle_tube(saddle, &change, dec)?;
            let bound = saddle.bound().upper.clone();
            if let Some(moment) = &change.moment {
                for (i, part) in moment.iter().enumerate() {
                    ring.moment[i] = &ring.moment[i] + &(part * &sign);
                    if !tube.is_zero() {
                        let [[low, _], [_, high]] = saddle.range(i, &bound, dec)?;
                        let far = low.abs().max(high.abs()) + bound.clone();
                        ring.moment_slack[i] = &ring.moment_slack[i] + &(far * tube.clone());
                    }
                }
            }
            ring.tube = ring.tube + tube;
        }
        let mut measured: Vec<(
            &catallax_shapes::RollingSphereDefinition,
            crate::saddle::Change,
        )> = Vec::new();
        for (band, added) in &self.spheres {
            let mut reused = None;
            for (source, change) in &measured {
                if let Some(change) = crate::saddle::sphere::reflected(change, source, band, dec)? {
                    reused = Some(change);
                    break;
                }
            }
            let change = match reused {
                Some(change) => change,
                None => {
                    let change = crate::saddle::sphere::change(band, width, moments, dec)?;
                    measured.push((band.definition(), change.clone()));
                    change
                }
            };
            let sign = exact(&Rational::from_i64(if *added { 1 } else { -1 }));
            ring.volume = &ring.volume + &(&change.pocket * &sign);
            ring.area = &ring.area + &change.area;
            let bound = band.bound().upper.clone();
            let tube = tube(&bound, band.definition().radius(), &change, dec)?;
            if let Some(moment) = &change.moment {
                for (i, part) in moment.iter().enumerate() {
                    ring.moment[i] = &ring.moment[i] + &(part * &sign);
                    if !tube.is_zero() {
                        let [[low, _], [_, high]] = band.range(i, &bound, dec)?;
                        let far = low.abs().max(high.abs()) + bound.clone();
                        ring.moment_slack[i] = &ring.moment_slack[i] + &(far * tube.clone());
                    }
                }
            }
            ring.tube = ring.tube + tube;
        }
        Ok(ring)
    }

    fn envelope_ring(&self, width: &Rational, dec: &mut Decider) -> Result<Ring, Refusal> {
        let first = self.envelopes[0].definition();
        let section = first.section(dec)?;
        let pi = Ball::pi_checked(BITS, dec)?;
        let two = exact(&Rational::from_i64(2));
        let origin = coordinates(
            &first.origin().minus(&catallax_shapes::Point3::origin()),
            dec,
        )?;
        let axis = coordinates(first.axis(), dec)?;
        let mut angle = &two * &pi;
        let mut volume = zero_ball();
        let mut area = zero_ball();
        let mut axial = zero_ball();
        let mut radial: [Ball; 3] = std::array::from_fn(|_| zero_ball());
        let mut tube = Rational::zero();
        let mut envelope_areas = Vec::with_capacity(self.envelopes.len());
        for envelope in &self.envelopes {
            let definition = envelope.definition();
            for (a, c) in [
                (&definition.origin().x, &first.origin().x),
                (&definition.origin().y, &first.origin().y),
                (&definition.origin().z, &first.origin().z),
                (&definition.axis().x, &first.axis().x),
                (&definition.axis().y, &first.axis().y),
                (&definition.axis().z, &first.axis().z),
                (definition.body(), first.body()),
                (definition.radius(), first.radius()),
            ] {
                same_ring(a, c, dec)?;
            }
            let moments = envelope.contact_moments(width, dec)?;
            let half = &moments.half_angle;
            angle = &angle - &(&two * half);
            volume = &volume + &moments.volume;
            envelope_areas.push(moments.envelope_area.clone());
            area =
                &(&area + &moments.envelope_area) - &(&moments.plane_area + &moments.cylinder_area);
            axial = &axial + &moments.axial_moment;
            let toward = coordinates(definition.toward(), dec)?;
            let sine = half.sin();
            let arm = &moments.radial_moment - &(&(&two * &section.radial_moment) * &sine);
            for i in 0..3 {
                radial[i] = &radial[i] + &(&toward[i] * &arm);
            }
            tube = tube + transfer(envelope, &moments, dec)?;
        }
        let volume = &volume + &(&section.volume * &angle);
        let area = &area + &(&section.area_change * &angle);
        let axial = &axial + &(&section.axial_moment * &angle);
        let moment =
            std::array::from_fn(|i| &(&(&origin[i] * &volume) + &(&axis[i] * &axial)) + &radial[i]);
        let reach = {
            let radius = dec.enclose_at(first.radius(), BITS)?.hi();
            let mouth = dec.enclose_at(first.mouth(), BITS)?.hi();
            let offset = dec.enclose_at(first.offset(), BITS)?.hi();
            offset + mouth + radius.mul_pow2(1)
        };
        let moment_slack = std::array::from_fn(|i| {
            let far = origin[i].lo().abs().max(origin[i].hi().abs());
            (far + reach.clone()) * tube.clone()
        });
        Ok(Ring {
            volume,
            area,
            moment,
            tube,
            moment_slack,
            envelope_areas,
        })
    }

    pub(crate) fn scalar(
        &self,
        model: &Model,
        accuracy: &Accuracy,
        property: surface::Property,
        traces: &mut engine::TraceCache,
        dec: &mut Decider,
    ) -> Answered<Enclosure> {
        let base = traces.with_native_conics(!self.saddles.is_empty(), |traces| {
            crate::scalar(model, self.base, accuracy, property, traces, dec)
        })?;
        let width = self.sector_width(&base.value, accuracy, Scope::Represented);
        let ring = self.ring(&width, false, dec)?;
        let mut represented = None;
        let value = if property == surface::Property::Area {
            let definition = add(&base.value, &ring.area, &Rational::zero());
            represented = self.represented_area(model, &ring, accuracy, &definition, dec)?;
            represented.clone().unwrap_or(definition)
        } else {
            add(&base.value, &ring.volume, &ring.tube)
        };
        let mut marks = base.marks;
        accuracy.mark(
            &mut marks,
            &[&value],
            if property == surface::Property::Area {
                "area"
            } else {
                "volume"
            },
        );
        if property == surface::Property::Area && represented.is_none() {
            self.mark_area_scope(&mut marks);
        }
        Ok(Answer::with_marks(value, marks))
    }

    pub(crate) fn properties(
        &self,
        model: &Model,
        accuracy: &Accuracy,
        scope: Scope,
        traces: &mut engine::TraceCache,
        dec: &mut Decider,
    ) -> Answered<VolumeAreaCentroid> {
        let base = traces.with_native_conics(!self.saddles.is_empty(), |traces| {
            crate::base_properties(model, self.base, accuracy, traces, dec)
        })?;
        let width = self.sector_width(&base.value.volume, accuracy, scope);
        let mut ring = self.ring(&width, true, dec)?;
        if scope == Scope::Definition {
            ring.tube = Rational::zero();
            ring.moment_slack = std::array::from_fn(|_| Rational::zero());
        }
        let volume = add(&base.value.volume, &ring.volume, &ring.tube);
        let area = add(&base.value.area, &ring.area, &Rational::zero());
        let represented = if scope == Scope::Represented {
            self.represented_area(model, &ring, accuracy, &area, dec)?
        } else {
            None
        };
        let area = represented.clone().unwrap_or(area);
        let center = base
            .value
            .centroid
            .as_ref()
            .ok_or_else(|| refused("a blend base is a nonempty solid"))?;
        let mut centroid = Vec::with_capacity(3);
        for (i, center) in center.iter().enumerate() {
            let moment = add(
                &product(center, &base.value.volume),
                &ring.moment[i],
                &ring.moment_slack[i],
            );
            let quotient = ball(&moment).checked_div(&ball(&volume)).ok_or_else(|| {
                Refusal::new(
                    "measure.centroid-not-separated",
                    "public volume enclosure includes zero after exhausting its budget",
                )
            })?;
            centroid.push(Enclosure::from_ball(&quotient));
        }
        let centroid: [Enclosure; 3] = centroid.try_into().expect("three coordinates");
        let mut marks = base.marks;
        accuracy.mark(&mut marks, &[&volume], "volume");
        accuracy.mark(&mut marks, &[&area], "area");
        if scope == Scope::Represented && represented.is_none() {
            self.mark_area_scope(&mut marks);
        }
        accuracy.mark(&mut marks, &centroid.each_ref(), "centroid");
        Ok(Answer::with_marks(
            VolumeAreaCentroid {
                volume,
                area,
                centroid: Some(centroid),
            },
            marks,
        ))
    }
}

fn add(base: &Enclosure, ring: &Ball, slack: &Rational) -> Enclosure {
    Enclosure {
        lo: &(&base.lo + &ring.lo()) - slack,
        hi: &(&base.hi + &ring.hi()) + slack,
    }
}

impl Junction<'_> {
    // The area enclosure that holds both the exact definition and the fitted
    // carriers, or `None` where the carriers' area is not covered yet: a
    // saddle or sphere band, or an envelope face that is not its whole patch.
    // Each envelope face is its whole carrier patch, whose area is enclosed
    // directly. Its mouth row stays on the mouth circle, so the neighbouring
    // mouth face is unchanged; its rail row moves only along the cylinder
    // axis, by at most the fit bound U, so the cylinder face changes by at
    // most U times the rail's angular length, below 2 pi body.
    fn represented_area(
        &self,
        model: &Model,
        ring: &Ring,
        accuracy: &Accuracy,
        definition: &Enclosure,
        dec: &mut Decider,
    ) -> Result<Option<Enclosure>, Refusal> {
        if !self.saddles.is_empty()
            || !self.spheres.is_empty()
            || ring.envelope_areas.len() != self.envelopes.len()
        {
            return Ok(None);
        }
        for &face in &self.envelope_faces {
            if !catallax_topology::envelope_retains_chart(model, face, dec)? {
                return Ok(None);
            }
        }
        let count = Rational::from_i64(2 * self.envelopes.len().max(1) as i64);
        let target =
            (&accuracy.relative * &definition.lo.abs()).max(accuracy.absolute.clone()) / count;
        let pi = Ball::pi_checked(BITS, dec)?;
        let mut change = zero_ball();
        let mut strips = Rational::zero();
        let carriers = self.carrier_areas(&target, dec)?;
        for ((envelope, exact_area), carrier) in self
            .envelopes
            .iter()
            .zip(&ring.envelope_areas)
            .zip(carriers)
        {
            change = &change + &(&carrier - exact_area);
            let body = dec
                .enclose_at(envelope.definition().body(), BITS)?
                .hi()
                .abs();
            let length = (&exact(&Rational::from_i64(2)) * &pi).hi() * body;
            strips = strips + length * envelope.bound().upper.clone();
        }
        let represented = add(definition, &change, &strips);
        Ok(Some(Enclosure {
            lo: definition.lo.clone().min(represented.lo),
            hi: definition.hi.clone().max(represented.hi),
        }))
    }

    // One carrier's area enclosure, one per envelope, in the envelopes' order,
    // with the request's own work accounting.
    fn carrier_areas(&self, target: &Rational, dec: &mut Decider) -> Result<Vec<Ball>, Refusal> {
        carrier_areas_of(
            self.envelopes.iter().map(|envelope| envelope.carrier()),
            target,
            dec,
        )
    }

    // The worst measured distance of a fitted carrier to its definition.
    fn worst_fit(&self) -> Rational {
        self.envelopes
            .iter()
            .map(|envelope| envelope.bound())
            .chain(self.saddles.iter().map(|(saddle, _)| saddle.bound()))
            .chain(self.spheres.iter().map(|(band, _)| band.bound()))
            .fold(Rational::zero(), |worst, bound| {
                worst.max(bound.upper.clone())
            })
    }

    // States that a represented area encloses the exact definition only.
    fn mark_area_scope(&self, marks: &mut MarkLog) {
        marks.push(Mark::new(
            AREA_OF_DEFINITION,
            BoundKind::Distance,
            self.worst_fit(),
            "the area encloses the exact definition only; the fitted carriers lie within this \
             distance of it and their own area is not bounded",
        ));
    }

    // A sector shares the requested width with its neighbours. A represented
    // answer needs no more than the carrier's own displacement allows.
    fn sector_width(&self, magnitude: &Enclosure, accuracy: &Accuracy, scope: Scope) -> Rational {
        let least = Rational::from_ratio(1, 1 << 40).expect("nonzero");
        let requested = (&accuracy.relative * &magnitude.lo.abs())
            .max(accuracy.absolute.clone())
            .mul_pow2(-6);
        let width = if scope == Scope::Definition {
            requested
        } else {
            requested.max(self.worst_fit())
        };
        width.max(least)
    }
}

// One area enclosure per carrier, in the carriers' order. A carrier's
// enclosure is a pure function of the carrier, the width and the precision, so
// every value is the same on any number of threads and the sum in
// `represented_area` keeps the carriers' order. Helpers take the carriers as
// they come and each records its whole computation; the writer adopts their
// records in order, so the request's allowance, refusals and work count stay
// those of one thread. A carrier whose record cannot be adopted, whose
// helper unwound, and every carrier when no helper spawned, is computed on
// the writer.
fn carrier_areas_of<'a>(
    carriers: impl IntoIterator<Item = &'a BSplineSurface>,
    target: &Rational,
    dec: &mut Decider,
) -> Result<Vec<Ball>, Refusal> {
    let carriers: Vec<&BSplineSurface> = carriers.into_iter().collect();
    let evaluate = |carrier: &BSplineSurface, worker: &Decider| {
        carrier.certified_area(target, CARRIER_CELLS, CARRIER_BITS, worker)
    };
    let Some(helpers) = (carriers.len() > 1)
        .then(|| catallax_numbers::Helpers::new(dec))
        .flatten()
    else {
        return carriers
            .iter()
            .map(|carrier| evaluate(carrier, dec))
            .collect();
    };
    let board = Board {
        slots: (0..carriers.len())
            .map(|_| Mutex::new(Slot::Open))
            .collect(),
        changed: Condvar::new(),
        next: AtomicUsize::new(0),
    };
    std::thread::scope(|scope| {
        struct Stop<'h>(&'h catallax_numbers::Helpers);
        impl Drop for Stop<'_> {
            fn drop(&mut self) {
                self.0.stop();
            }
        }
        let _stop = Stop(&helpers);
        let mut started = 0;
        for _ in 0..helpers.count().min(carriers.len()) {
            let spawned = std::thread::Builder::new()
                .name("catallax-carrier-area".into())
                .stack_size(helpers.stack_size())
                .spawn_scoped(scope, || {
                    helpers.run(|worker| explore(&carriers, &board, &helpers, &evaluate, worker));
                });
            if spawned.is_err() {
                break;
            }
            started += 1;
        }
        let mut areas = Vec::with_capacity(carriers.len());
        for (index, carrier) in carriers.iter().enumerate() {
            let area = if started == 0 {
                evaluate(carrier, dec)
            } else {
                take(&board, index, dec, |dec| evaluate(carrier, dec))
            };
            areas.push(area?);
        }
        Ok(areas)
    })
}

// A carrier's place in the in-order board: free, taken by one thread, done
// with its area and the charges the helper recorded for the writer to adopt,
// or failed before it could record anything, its evaluation having unwound.
enum Slot {
    Open,
    Taken,
    Done(
        Box<Result<Ball, Refusal>>,
        Option<catallax_numbers::Charges>,
    ),
    Failed,
}

struct Board {
    slots: Vec<Mutex<Slot>>,
    changed: Condvar,
    next: AtomicUsize,
}

fn locked<T>(mutex: &Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    mutex
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

// Settles a taken slot when its helper leaves, whether it returns or unwinds:
// a slot left taken would be waited for forever.
struct Settle<'a> {
    board: &'a Board,
    index: usize,
}

impl Drop for Settle<'_> {
    fn drop(&mut self) {
        let mut slot = locked(&self.board.slots[self.index]);
        if matches!(*slot, Slot::Taken) {
            *slot = Slot::Failed;
        }
        drop(slot);
        self.board.changed.notify_all();
    }
}

// Takes an open carrier, or the next one, and records its whole computation.
fn explore<F>(
    carriers: &[&BSplineSurface],
    board: &Board,
    helpers: &catallax_numbers::Helpers,
    evaluate: &F,
    worker: &mut Decider,
) where
    F: Fn(&BSplineSurface, &Decider) -> Result<Ball, Refusal>,
{
    loop {
        if helpers.stopped() {
            return;
        }
        let index = board
            .next
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let Some(carrier) = carriers.get(index) else {
            return;
        };
        {
            let mut slot = locked(&board.slots[index]);
            if !matches!(*slot, Slot::Open) {
                continue;
            }
            *slot = Slot::Taken;
        }
        // The slot is settled on every path out, a panic included. A helper's
        // unwind stops at the board: the writer computes that carrier itself,
        // as it does for a record it cannot adopt.
        let settle = Settle { board, index };
        let recorded = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            catallax_numbers::record(|| evaluate(carrier, worker))
        }));
        if let Ok((value, charges)) = recorded {
            *locked(&board.slots[index]) = Slot::Done(Box::new(value), charges);
        }
        drop(settle);
    }
}

// The writer's area for one carrier: the helper's record adopted in order, or
// the carrier computed here when no record can be adopted or the helper
// failed before recording one.
fn take<F>(board: &Board, index: usize, dec: &mut Decider, evaluate: F) -> Result<Ball, Refusal>
where
    F: FnOnce(&mut Decider) -> Result<Ball, Refusal>,
{
    let mut slot = locked(&board.slots[index]);
    loop {
        if matches!(*slot, Slot::Failed) {
            drop(slot);
            return evaluate(dec);
        }
        if let Slot::Done(..) = *slot {
            break;
        }
        slot = board
            .changed
            .wait(slot)
            .unwrap_or_else(std::sync::PoisonError::into_inner);
    }
    let Slot::Done(value, charges) = std::mem::replace(&mut *slot, Slot::Taken) else {
        unreachable!("the wait stops only at a done slot")
    };
    drop(slot);
    if let Some(charges) = charges
        && charges.adopt(dec)
    {
        return *value;
    }
    evaluate(dec)
}

// The tube bound of the module note. The boundary length P is at most
// 2 pi (body + mouth) + 8 r: the rail turns less than once about the cylinder
// and climbs r on each side, the mouth row is part of one circle, and each
// meridian is a quarter of the tube circle.
fn transfer(
    envelope: &RollingEnvelope,
    moments: &ContactMoments,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let definition = envelope.definition();
    let upper = envelope.bound().upper.clone();
    let radius = dec.enclose_at(definition.radius(), BITS)?;
    if upper.mul_pow2(8) > radius.lo() {
        return Err(refused(
            "a carrier's bound is too wide for the tube transfer",
        ));
    }
    let body = dec.enclose_at(definition.body(), BITS)?.hi();
    let mouth = dec.enclose_at(definition.mouth(), BITS)?.hi();
    let pi = Ball::pi_checked(BITS, dec)?.hi();
    let length = pi.mul_pow2(1) * (body + mouth) + radius.hi().mul_pow2(3);
    let normal = (&upper * &(&moments.envelope_area.hi() + &(&upper * &moments.turning)))
        .mul_pow2(1)
        * (Rational::one() + &upper / &radius.lo());
    Ok(normal + pi * (&upper * &upper) * length)
}

// The same bound for a saddle, whose rails are all its boundary.
fn saddle_tube(
    saddle: &RollingSaddle,
    change: &crate::saddle::Change,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    tube(
        &saddle.bound().upper,
        saddle.definition().radius(),
        change,
        dec,
    )
}

// The tube bound of a ring band with carrier bound `upper` and ball `radius`.
fn tube(
    upper: &Rational,
    radius: &Real,
    change: &crate::saddle::Change,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let upper = upper.clone();
    if upper.is_zero() {
        dec.check_deadline()?;
        return Ok(Rational::zero());
    }
    let radius = dec.enclose_at(radius, BITS)?;
    if upper.mul_pow2(8) > radius.lo() {
        return Err(refused(
            "a carrier's bound is too wide for the tube transfer",
        ));
    }
    let pi = Ball::pi_checked(BITS, dec)?.hi();
    let normal = (&upper * &(&change.band.hi() + &(&upper * &change.turning.hi()))).mul_pow2(1)
        * (Rational::one() + &upper / &radius.lo());
    Ok(normal + pi * (&upper * &upper) * change.length.hi())
}

/// A starved decision must reach the caller with its code and bound intact.
/// See docs/NOT-EXACT-INVENTORY.md.
#[cfg(test)]
mod undecided_evidence {
    use super::*;
    use catallax_numbers::{BudgetPolicy, NotExact};

    fn starved() -> Decider {
        Decider::new(BudgetPolicy {
            max_precision_bits: 64,
            max_work: 0,
        })
    }

    fn undecided(value: &Real) -> NotExact {
        match starved().sign(value) {
            Outcome::NotExact(decision) => decision,
            other => panic!("an empty allowance must stop undecided: {other:?}"),
        }
    }

    #[test]
    fn an_undecided_ring_coordinate_keeps_its_decision() {
        let [a, c] = [Real::from(2).sqrt(), Real::from(1)];
        let decision = undecided(&(&a - &c));
        let refusal = same_ring(&a, &c, &mut starved()).unwrap_err();
        assert_eq!(refusal.code, "measure.unsupported-surface");
        assert_eq!(refusal.undecided(), Some(&decision));
    }
}

/// The tests here set the process's thread count, so they share one lock.
#[cfg(test)]
static THREAD_COUNT_LOCK: Mutex<()> = Mutex::new(());

/// Carrier areas are request work, so the request's allowance and refusals do
/// not depend on how many threads enclose them.
#[cfg(test)]
mod carrier_area_accounting {
    use super::*;
    use catallax_numbers::{set_threads, with_operation_budget};
    use catallax_shapes::Point3;
    use std::num::NonZeroUsize;

    // A quarter cylinder of radius 2 over height 3: a cheap carrier whose net
    // construction and cells still charge the enclosing request.
    pub(super) fn carrier(dec: &mut Decider) -> BSplineSurface {
        let zero = Real::zero();
        let radius = Real::from_i64(2);
        let length = Real::from_i64(3);
        let half = Real::from_i64(2).sqrt() / Real::from_i64(2);
        let arc = |z: &Real| {
            (
                [
                    Point3::new(radius.clone(), zero.clone(), z.clone()),
                    Point3::new(radius.clone(), radius.clone(), z.clone()),
                    Point3::new(zero.clone(), radius.clone(), z.clone()),
                ],
                [Real::one(), half.clone(), Real::one()],
            )
        };
        let (bottom, weights) = arc(&zero);
        let (top, _) = arc(&length);
        let mut controls = Vec::new();
        let mut all_weights = Vec::new();
        for i in 0..3 {
            controls.extend([bottom[i].clone(), top[i].clone()]);
            all_weights.extend([weights[i].clone(), weights[i].clone()]);
        }
        let knots = |degree: usize| {
            let mut knots = vec![Rational::zero(); degree + 1];
            knots.extend(vec![Rational::one(); degree + 1]);
            knots
        };
        BSplineSurface::new(
            [2, 1],
            [knots(2), knots(1)],
            [3, 2],
            controls,
            all_weights,
            dec,
        )
        .unwrap()
    }

    #[test]
    fn carrier_work_is_charged_whatever_the_thread_count() {
        let _threads = super::THREAD_COUNT_LOCK
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        struct Restore;
        impl Drop for Restore {
            fn drop(&mut self) {
                set_threads(None);
            }
        }
        let _restore = Restore;
        let mut setup = Decider::default();
        let first = carrier(&mut setup);
        let second = carrier(&mut setup);
        let carriers = [&first, &second];
        let target = Rational::from_ratio(1, 100).unwrap();

        // A successful request charges the carriers' work, whatever threads
        // the kernel may run their cells on.
        set_threads(Some(NonZeroUsize::new(1).unwrap()));
        let (answer, serial) = with_operation_budget(1_000_000_000, || {
            carrier_areas_of(carriers.iter().copied(), &target, &mut Decider::default())
        });
        assert!(answer.is_ok());
        assert!(serial.used > 0);
        assert!(!serial.exhausted);
        for threads in [2, 4] {
            set_threads(Some(NonZeroUsize::new(threads).unwrap()));
            let (answer, work) = with_operation_budget(1_000_000_000, || {
                carrier_areas_of(carriers.iter().copied(), &target, &mut Decider::default())
            });
            assert!(answer.is_ok(), "{threads} threads");
            assert!(work.used > 0, "{threads} threads charged no carrier work");
            assert_eq!(work.used, serial.used, "{threads} threads charge");
            assert_eq!(work.exhausted, serial.exhausted, "{threads} threads");
        }

        // An allowance too small for the carriers' net construction refuses
        // the same way, after the same charge, at every thread count.
        let limit = 128;
        assert!(limit < serial.used);
        set_threads(Some(NonZeroUsize::new(1).unwrap()));
        let (serial_answer, serial_work) = with_operation_budget(limit, || {
            carrier_areas_of(carriers.iter().copied(), &target, &mut Decider::default())
        });
        assert!(serial_work.exhausted, "{serial_work:?}");
        let serial_refusal = match serial_answer {
            Ok(areas) => panic!("a net-construction allowance enclosed {areas:?}"),
            Err(refusal) => refusal,
        };
        for threads in [2, 4] {
            set_threads(Some(NonZeroUsize::new(threads).unwrap()));
            let (answer, work) = with_operation_budget(limit, || {
                carrier_areas_of(carriers.iter().copied(), &target, &mut Decider::default())
            });
            assert_eq!(work.used, serial_work.used, "{threads} threads");
            assert_eq!(work.exhausted, serial_work.exhausted, "{threads} threads");
            match answer {
                Ok(areas) => panic!("{threads} threads enclosed {areas:?} uncharged"),
                Err(refusal) => assert_eq!(refusal.code, serial_refusal.code, "{threads} threads"),
            }
        }
    }
}

/// A helper's evaluation may unwind. The board must then settle that carrier
/// so the writer neither waits for it forever nor loses its work accounting:
/// the writer computes the carrier itself, as it does for a record it cannot
/// adopt.
#[cfg(test)]
mod carrier_unwind {
    use super::carrier_area_accounting::carrier;
    use super::*;
    use catallax_numbers::{Helpers, set_threads};
    use std::num::NonZeroUsize;
    use std::sync::mpsc;
    use std::time::Duration;

    fn scene() -> (Vec<Ball>, bool) {
        let mut setup = Decider::default();
        let first = carrier(&mut setup);
        let second = carrier(&mut setup);
        let carriers = [&first, &second];
        let target = Rational::from_ratio(1, 100).unwrap();
        let evaluate = |carrier: &BSplineSurface, worker: &Decider| {
            if std::thread::current().name() == Some("carrier-unwind-helper") {
                panic!("injected carrier evaluation panic");
            }
            carrier.certified_area(&target, CARRIER_CELLS, CARRIER_BITS, worker)
        };
        let mut dec = Decider::default();
        let helpers = Helpers::new(&dec).expect("two threads run helpers");
        let board = Board {
            slots: (0..carriers.len())
                .map(|_| Mutex::new(Slot::Open))
                .collect(),
            changed: Condvar::new(),
            next: AtomicUsize::new(0),
        };
        std::thread::scope(|scope| {
            let helper = std::thread::Builder::new()
                .name("carrier-unwind-helper".into())
                .stack_size(helpers.stack_size())
                .spawn_scoped(scope, || {
                    helpers.run(|worker| explore(&carriers, &board, &helpers, &evaluate, worker));
                })
                .expect("helper thread");
            let areas: Vec<_> = (0..carriers.len())
                .map(|index| {
                    take(&board, index, &mut dec, |dec| {
                        evaluate(carriers[index], dec)
                    })
                })
                .collect();
            let joined = helper.join().is_ok();
            let areas: Vec<Ball> = areas
                .into_iter()
                .map(|area| area.expect("the writer's own evaluation encloses the carrier"))
                .collect();
            // The writer computes each carrier exactly as a serial request
            // would: the helpers' work is not lost and no ball changes.
            let direct: Vec<Ball> = carriers
                .iter()
                .map(|carrier| {
                    carrier
                        .certified_area(&target, CARRIER_CELLS, CARRIER_BITS, &dec)
                        .unwrap()
                })
                .collect();
            for (index, (area, direct)) in areas.iter().zip(&direct).enumerate() {
                assert_eq!(
                    (area.lo(), area.hi()),
                    (direct.lo(), direct.hi()),
                    "carrier {index}"
                );
            }
            (areas, joined)
        })
    }

    #[test]
    fn a_panicking_carrier_helper_leaves_the_writer_a_terminal_outcome() {
        let _threads = super::THREAD_COUNT_LOCK
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        struct Restore;
        impl Drop for Restore {
            fn drop(&mut self) {
                set_threads(None);
            }
        }
        let _restore = Restore;
        set_threads(Some(NonZeroUsize::new(2).unwrap()));

        let (done, finished) = mpsc::channel();
        let scene = std::thread::spawn(move || {
            let outcome = scene();
            let _ = done.send((outcome.0.clone(), outcome.1));
            outcome
        });
        let (areas, joined) = match finished.recv_timeout(Duration::from_secs(20)) {
            Ok(outcome) => outcome,
            Err(mpsc::RecvTimeoutError::Timeout) => {
                panic!("a panicking carrier helper left the writer waiting")
            }
            Err(mpsc::RecvTimeoutError::Disconnected) => match scene.join() {
                Ok(_) => unreachable!("the scene sends its outcome before it returns"),
                Err(panic) => std::panic::resume_unwind(panic),
            },
        };
        scene.join().unwrap();

        // The helper's unwind stopped at the board: every carrier is settled,
        // the helper joined cleanly, and the writer computed each area itself.
        assert!(joined, "the helper panicked out of its board");
        assert_eq!(areas.len(), 2);
    }
}
