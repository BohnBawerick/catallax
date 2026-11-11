//! Volume, area and centroid of a blend with rolling-envelope or rolling-saddle
//! faces.
//!
//! The engine has no chart for a fitted carrier or for a cylinder trimmed by its
//! rail. The blend changes its unfilleted base only inside one ring about the
//! cylinder axis, so each property is the base's plus what the ring adds: the
//! exact section over the ring's untouched angle and one exact sector per mouth.
//! A saddle's change is one set of integrals round its ring, in `crate::saddle`.
//!
//! Volume and centroid enclose both the exact definition and the fitted carrier.
//! A homotopy from envelope to carrier moves each point at most U and keeps the
//! boundary rows on their neighbouring surfaces, so any point between the two
//! bodies lies within U of the exact envelope. A point of that tube reaches the
//! envelope along a normal or at its boundary. With principal curvatures 1/r and
//! k, the normal part is at most `2 U (1 + U / r) (A + U K)` for envelope area A
//! and K, the integral of |k| over it. The boundary part is at most `pi U^2 P`
//! for boundary length P: half-tubes around its curves and four corner balls,
//! given U far below every curve radius. The area is the exact definition's,
//! because a position bound says nothing about the carrier's area.

use catallax_numbers::{Answer, Answered, Ball, Decider, Outcome, Rational, Real, Refusal};
use catallax_shapes::{
    ContactMoments, Curve, RollingEnvelope, RollingSaddle, SplineSupport, Surface, Vector3,
};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, saddle_added};

use crate::{Accuracy, Enclosure, VolumeAreaCentroid, engine, product, surface};

const BITS: u32 = 128;

/// What a measured enclosure has to contain.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Scope {
    /// The exact definition alone.
    Definition,
    /// The definition and the fitted carrier.
    Represented,
}

pub(crate) struct Junction<'a> {
    base: BodyId,
    envelopes: Vec<&'a RollingEnvelope>,
    // Each with whether its pocket adds material.
    saddles: Vec<(&'a RollingSaddle, bool)>,
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
    let envelopes: Vec<&RollingEnvelope> = model
        .faces_of_body(body)
        .into_iter()
        .filter_map(
            |face| match model.geom().surface(model.face(face).surface) {
                Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) => Some(envelope),
                _ => None,
            },
        )
        .collect();
    let mut saddles = Vec::new();
    for face in model.faces_of_body(body) {
        match model.geom().surface(model.face(face).surface) {
            Surface::BSpline(SplineSupport::RollingSaddle(saddle)) => {
                saddles.push((saddle, saddle_added(model, face)?));
            }
            Surface::BSpline(SplineSupport::RollingSphere(_)) => {
                return Err(refused(
                    "a band between a sphere and a cylinder has no certified change integral",
                ));
            }
            _ => {}
        }
    }
    if envelopes.is_empty() && saddles.is_empty() {
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
        saddles,
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
        let value = if property == surface::Property::Area {
            add(&base.value, &ring.area, &Rational::zero())
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
            let bounds = self
                .envelopes
                .iter()
                .map(|envelope| envelope.bound())
                .chain(self.saddles.iter().map(|(saddle, _)| saddle.bound()));
            bounds.fold(requested, |worst, bound| worst.max(bound.upper.clone()))
        };
        width.max(least)
    }
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
    let upper = saddle.bound().upper.clone();
    if upper.is_zero() {
        dec.check_deadline()?;
        return Ok(Rational::zero());
    }
    let radius = dec.enclose_at(saddle.definition().radius(), BITS)?;
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
