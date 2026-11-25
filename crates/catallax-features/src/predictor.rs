//! Exact or fitted for one operation: the quick look and the cap tiers.
//!
//! Checkpoint 2 decisions, section 4.5, the predictor. An exact attempt that
//! has a fitted construction behind it runs under a work cap
//! (`catallax_numbers::capped_attempt`). The cap is one of two tiers, chosen
//! by a quick look at the attempt's inputs. A wrong tier costs at most the
//! cap, because the construction that takes over is the one that follows any
//! refusal of the attempt.
//!
//! The look is taken when the attempt reaches the smaller of the two caps
//! (`catallax_numbers::cap_attempt_when_reached`). An attempt that ends
//! before that finishes the same under either cap, so it is never looked at
//! and pays nothing. One that gets there has charged 2,000,000 units, and the
//! look reads at most a hundredth of that ([`LOOK_READS`]): "It costs at most
//! 1% of the operation." The charge at which it is taken is a count of work
//! units, the same on every machine and at every thread count.
//!
//! The design is `docs/archive/research/cp2-grill/predictor-research.md`, design (e):
//! a prediction sets the cap and a deterministic limit enforces it (Purohit,
//! Svitkina and Kumar, "Improving online algorithms via ML predictions",
//! NeurIPS 2018). No kernel read for that note predicts the cost of a whole
//! operation before running it: OpenCASCADE and SMLib run one construction
//! and report its failure, and Parasolid and ACIS fall back to a tolerant
//! edge after the exact intersection fails.
//!
//! What the look reads:
//!
//! - D, the product of the two participants' implicit degrees. A plane has
//!   degree 1, a quadric 2 and a torus 4. A tensor-product patch of degree
//!   (p, q) has implicit degree 2pq, so a bicubic has 18 (Sederberg, Anderson
//!   and Goldman, "Implicit representation of parametric curves and
//!   surfaces", Computer Vision, Graphics, and Image Processing 28, 1984).
//!   Elimination between two surfaces works in degree D (Bezout).
//! - L, the largest limb count of the rationals that define the
//!   participants, numerator and denominator added, at 64 bits a limb. A
//!   product of two such values costs L^2 limb operations (Knuth, The Art of
//!   Computer Programming, volume 2, section 4.3.1, algorithm M), which is
//!   how the decider prices it.
//! - Whether the two participants are proven apart. A pair that is cannot
//!   meet, so it adds nothing to the cost (Sederberg and Meyers, "Loop
//!   detection in surface patch intersections", Computer Aided Geometric
//!   Design 5, 1988, for the same use of a cheap separation test).
//!
//! Everything the look does is a read and is counted: each number of a
//! surface, each surface kind, and each step through the topology that
//! finds the pairs, which is one for a selected edge, one for each use of
//! it and one for each pair. A look that would pass [`LOOK_READS`] stops
//! there and takes the high tier, the smaller cap, so its cost is bounded
//! whatever the body.
//!
//! The score of a pair is D x L^2 and the score of an operation is that of
//! its costliest pair. The look reads no normal angle: the decisions allow
//! one only to demote spline pairs, and D already places those. It keeps no
//! memo between steps and has no learned part.

use catallax_numbers::Real;
use catallax_shapes::Surface;
use catallax_shapes::SurfaceId;
use catallax_topology::{EdgeId, Model};
use std::{cell::RefCell, collections::BTreeMap};

/// The cap of an exact attempt on low-complexity inputs, in work units: the
/// work target of one build. Decisions 4.5, predictor step 3: "low-complexity
/// inputs get the 20M target for the exact attempt".
pub const LOW_COMPLEXITY_CAP: u64 = 20_000_000;

/// The cap of an exact attempt on high-complexity inputs. Decisions 4.5,
/// predictor step 3: "fitted splines stored as 53-bit fractions get about
/// 2M", enough for a case that turns out to be trivial.
pub const HIGH_COMPLEXITY_CAP: u64 = 2_000_000;

/// The largest degree product of the low tier. The research note keeps the
/// pairs of a plane, a quadric and a torus there ("Plane 1, quadric 2, torus
/// 4 keep D <= 8"), after ESOLID's limit of algebraic degree 4 for exact
/// boundary evaluation (Keyser, Culver, Foskey, Krishnan and Manocha,
/// "ESOLID: a system for exact boundary evaluation", Computer-Aided Design
/// 36, 2004).
pub const LOW_COMPLEXITY_DEGREE: u64 = 8;

/// The limb count of a small rational: one limb for the numerator and one
/// for the denominator, which is the least a rational has and what a
/// machine integer or a binary64 value takes.
pub const SMALL_LIMBS: u64 = 2;

/// The largest score, D x L^2, of the low tier: the largest low-complexity
/// degree product on small rationals.
pub const LOW_COMPLEXITY_SCORE: u64 = LOW_COMPLEXITY_DEGREE * SMALL_LIMBS * SMALL_LIMBS;

/// The most values one look reads: a hundredth of the work an attempt has
/// charged when it is looked at, which is [`HIGH_COMPLEXITY_CAP`]. Decisions
/// 4.5, predictor step 2: "It costs at most 1% of the operation." A surface
/// with more numbers than are left is not read to its end, and is not taken
/// for a low-complexity one.
pub const LOOK_READS: u64 = HIGH_COMPLEXITY_CAP / 100;

/// The cap tier of an exact attempt.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tier {
    /// Low-complexity inputs: [`LOW_COMPLEXITY_CAP`].
    Low,
    /// High-complexity inputs: [`HIGH_COMPLEXITY_CAP`].
    High,
}

impl Tier {
    /// The work cap of an exact attempt of this tier.
    pub fn cap(self) -> u64 {
        match self {
            Self::Low => LOW_COMPLEXITY_CAP,
            Self::High => HIGH_COMPLEXITY_CAP,
        }
    }

    /// The tier's name in a record.
    pub fn name(self) -> &'static str {
        match self {
            Self::Low => "low",
            Self::High => "high",
        }
    }
}

/// What the quick look read of one operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Look {
    /// D of the costliest pair: the product of its two implicit degrees.
    pub degree: u64,
    /// L of the costliest pair: the largest limb count among its rationals.
    pub limbs: u64,
    /// Everything read: one for each surface kind and each number, and one
    /// for each selected edge, each use of it and each pair. This is the
    /// look's own work, counted in `OperationWork::uncharged` and charged to
    /// nothing.
    pub read: u64,
}

impl Look {
    /// A look that has read nothing: no pair, so the least score.
    pub const NONE: Self = Self {
        degree: 1,
        limbs: SMALL_LIMBS,
        read: 0,
    };

    /// D x L^2 of the costliest pair.
    pub fn score(&self) -> u64 {
        score(self.degree, self.limbs)
    }

    /// The cap tier these inputs get.
    pub fn tier(&self) -> Tier {
        if self.score() <= LOW_COMPLEXITY_SCORE {
            Tier::Low
        } else {
            Tier::High
        }
    }

    /// Takes one more pair into the look. A pair proven apart adds nothing.
    pub fn pair(&mut self, a: Participant, b: Participant, proven_apart: bool) {
        if proven_apart {
            return;
        }
        let degree = a.degree.saturating_mul(b.degree);
        let limbs = a.limbs.max(b.limbs);
        if score(degree, limbs) > self.score() {
            self.degree = degree;
            self.limbs = limbs;
        }
    }

    /// Takes one more number into every pair, such as a radius both
    /// participants of each pair are combined with.
    pub fn number(&mut self, value: &Real) {
        self.read += 1;
        self.limbs = self.limbs.max(limbs(value));
    }

    /// Takes one read for a step through the topology: a selected edge, a
    /// use of it or a pair. `false` once the look has read all it may.
    fn step(&mut self) -> bool {
        if self.read >= LOOK_READS {
            return false;
        }
        self.read += 1;
        true
    }

    /// A look that stopped before it had read every pair. What it did not
    /// read may be anything, so it has no degree to go by and is not taken
    /// for a low-complexity one.
    fn unfinished(mut self) -> Self {
        self.degree = u64::MAX;
        self
    }
}

fn score(degree: u64, limbs: u64) -> u64 {
    degree.saturating_mul(limbs).saturating_mul(limbs)
}

/// One input of an operation as the look reads it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Participant {
    /// The implicit degree of its surface.
    pub degree: u64,
    /// The largest limb count among the rationals that define it.
    pub limbs: u64,
    /// The values read to learn the two.
    pub read: u64,
}

impl Participant {
    /// Reads a surface: its kind and every number of its definition.
    pub fn of(surface: &Surface) -> Self {
        Self::within(surface, u64::MAX)
    }

    /// Reads a surface as [`Self::of`] does, and at most `reads` values. A
    /// surface it could not read to its end has no degree to go by.
    pub fn within(surface: &Surface, reads: u64) -> Self {
        let mut read = Read {
            // One is for the kind.
            left: reads.saturating_sub(1),
            ..Read::default()
        };
        let degree = match surface {
            Surface::Plane(plane) => {
                read.point(&plane.origin);
                read.vector(&plane.normal);
                read.vector(&plane.ref_dir);
                1
            }
            Surface::Cylinder(cylinder) => {
                read.point(&cylinder.origin);
                read.vector(&cylinder.axis);
                read.number(&cylinder.radius_sq);
                read.vector(&cylinder.ref_dir);
                2
            }
            Surface::Cone(cone) => {
                read.point(&cone.apex);
                read.vector(&cone.axis);
                read.number(&cone.slope);
                read.vector(&cone.ref_dir);
                2
            }
            Surface::Sphere(sphere) => {
                read.point(&sphere.center);
                read.number(&sphere.radius);
                read.vector(&sphere.axis);
                read.vector(&sphere.ref_dir);
                2
            }
            Surface::EllipticCylinder(cylinder) => {
                read.point(&cylinder.origin);
                read.vector(&cylinder.axis);
                read.vector(&cylinder.major_dir);
                read.number(&cylinder.major_sq);
                read.number(&cylinder.minor_sq);
                2
            }
            Surface::Torus(torus) => {
                read.point(&torus.center);
                read.vector(&torus.axis);
                read.number(&torus.major);
                read.number(&torus.minor);
                read.vector(&torus.ref_dir);
                4
            }
            Surface::TorusSheet(sheet) => {
                read.point(&sheet.center);
                read.vector(&sheet.axis);
                read.number(&sheet.major);
                read.number(&sheet.minor);
                read.vector(&sheet.ref_dir);
                4
            }
            Surface::BSpline(support) => {
                let carrier = support.carrier();
                for control in carrier.controls() {
                    if read.short {
                        break;
                    }
                    read.point(control);
                }
                for weight in carrier.weights() {
                    if read.short {
                        break;
                    }
                    read.number(weight);
                }
                for knots in carrier.knots() {
                    for knot in knots {
                        if read.short {
                            break;
                        }
                        read.rational(knot);
                    }
                }
                let [p, q] = carrier.degrees();
                2u64.saturating_mul(p as u64).saturating_mul(q as u64)
            }
            // A kind added after this was written has no degree here to
            // read. It is not taken for a cheap one.
            _ => u64::MAX,
        };
        Self {
            // Neither is a surface whose numbers were not all read.
            degree: if read.short { u64::MAX } else { degree },
            limbs: read.limbs,
            // One more for the kind.
            read: (read.count + 1).min(reads),
        }
    }
}

// The numbers of one definition as they are read.
struct Read {
    limbs: u64,
    count: u64,
    // The values that may still be read, and whether one was left unread.
    left: u64,
    short: bool,
}

impl Default for Read {
    fn default() -> Self {
        Self {
            limbs: SMALL_LIMBS,
            count: 0,
            left: u64::MAX,
            short: false,
        }
    }
}

impl Read {
    // Whether one more value may be read, taking it off what is left.
    fn one_more(&mut self) -> bool {
        if self.left == 0 {
            self.short = true;
            return false;
        }
        self.left -= 1;
        self.count += 1;
        true
    }
    fn rational(&mut self, value: &catallax_numbers::Rational) {
        if self.one_more() {
            self.limbs = self.limbs.max(value.limbs());
        }
    }
    fn number(&mut self, value: &Real) {
        if self.one_more() {
            self.limbs = self.limbs.max(limbs(value));
        }
    }
    fn point(&mut self, point: &catallax_shapes::Point3) {
        self.number(&point.x);
        self.number(&point.y);
        self.number(&point.z);
    }
    fn vector(&mut self, vector: &catallax_shapes::Vector3) {
        self.number(&vector.x);
        self.number(&vector.y);
        self.number(&vector.z);
    }
}

// The limb count of a number that is a rational. One that is not, a root or
// an angle function, is decided through enclosures at a working precision
// and has no limbs of its own to read, so it counts as a small rational.
fn limbs(value: &Real) -> u64 {
    value
        .rational_ref()
        .map_or(SMALL_LIMBS, catallax_numbers::Rational::limbs)
}

/// Reads the edges a fillet rounds: for each, the surfaces of the faces that
/// use it. They meet along the edge, so no pair is proven apart.
///
/// The faces of an edge are found through the edge's own uses, so the look
/// costs by the selection and not by the size of the body. Each selected
/// edge, each use and each pair is one read, beside the numbers of the
/// surfaces. At [`LOOK_READS`] the look stops and is
/// [unfinished](Look::unfinished): the high tier.
///
/// A use in another body of the model is read as well. The model keeps
/// earlier bodies, and a face of one on the same edge lies on a surface
/// through the same curve. A pair can only raise the score.
pub fn of_fillet(model: &Model, edges: &[EdgeId], radius: &Real) -> Look {
    let mut look = Look::NONE;
    // One read is kept for the radius.
    look.number(radius);
    let mut read: BTreeMap<SurfaceId, Participant> = BTreeMap::new();
    for edge in edges {
        if !look.step() {
            return look.unfinished();
        }
        let mut surfaces: Vec<SurfaceId> = Vec::new();
        for coedge in model.coedges_of_edge(*edge) {
            if !look.step() {
                return look.unfinished();
            }
            let Some((_, face)) = model.owner_of_coedge(coedge) else {
                continue;
            };
            let surface = model.face(face).surface;
            if !surfaces.contains(&surface) {
                surfaces.push(surface);
            }
        }
        let mut participants = Vec::with_capacity(surfaces.len());
        for surface in surfaces {
            let participant = match read.get(&surface) {
                Some(participant) => *participant,
                None => {
                    let participant =
                        Participant::within(model.geom().surface(surface), LOOK_READS - look.read);
                    look.read += participant.read;
                    read.insert(surface, participant);
                    participant
                }
            };
            participants.push(participant);
        }
        for (index, a) in participants.iter().enumerate() {
            for b in &participants[index + 1..] {
                if !look.step() {
                    return look.unfinished();
                }
                look.pair(*a, *b, false);
            }
        }
    }
    look
}

/// Caps the exact fillet attempt that is open on this thread, if one is.
/// The attempt may charge the high tier's cap. If it gets there, the look
/// reads the selected edges and the attempt goes on under the cap of their
/// tier, which for the high tier stops it at that charge. Outside a capped
/// attempt this does nothing. The exact fillet builders call it once they
/// have resolved their selection, so the cap covers the whole attempt and
/// the selection is resolved once.
pub(crate) fn cap_fillet_attempt(
    model: &Model,
    edges: &[EdgeId],
    radius: &Real,
    step: &impl ToString,
) {
    if !catallax_numbers::capped_attempt_open() {
        return;
    }
    let step = step.to_string();
    LOG.with(|log| {
        if let Some(choices) = log.borrow_mut().as_mut() {
            choices.push(Choice {
                step: step.clone(),
                look: None,
                path: Path::Open,
                attempt_work: 0,
                operation_work: 0,
            });
        }
    });
    // The model is a persistent value, so the copy the look keeps shares its
    // tables with this one.
    let (model, edges, radius) = (model.clone(), edges.to_vec(), radius.clone());
    catallax_numbers::cap_attempt_when_reached(HIGH_COMPLEXITY_CAP, move || {
        let look = of_fillet(&model, &edges, &radius);
        catallax_numbers::count_uncharged_work(look.read);
        LOG.with(|log| {
            if let Some(choice) = log.borrow_mut().as_mut().and_then(|choices| {
                choices
                    .iter_mut()
                    .rev()
                    .find(|choice| choice.step == step && choice.path == Path::Open)
            }) {
                choice.look = Some(look);
            }
        });
        look.tier().cap()
    });
}

/// What became of an exact attempt that has a fitted construction behind it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Path {
    /// The attempt has not ended.
    Open,
    /// The attempt finished: the result is exact.
    Exact,
    /// The attempt refused or reached its cap, and the fitted construction
    /// built the result.
    Fitted,
    /// The attempt refused or reached its cap, and the fitted construction
    /// refused as well.
    Refused,
}

impl Path {
    /// The path's name in a record.
    pub fn name(self) -> &'static str {
        match self {
            Self::Open => "open",
            Self::Exact => "exact",
            Self::Fitted => "fitted",
            Self::Refused => "refused",
        }
    }
}

/// One exact attempt with a fitted construction behind it, as recorded by
/// [`with_log`].
#[derive(Clone, Debug)]
pub struct Choice {
    /// The recipe step.
    pub step: String,
    /// What the quick look read, or `None` when the attempt ended before it
    /// reached the high tier's cap and so was never looked at.
    pub look: Option<Look>,
    /// What became of the attempt.
    pub path: Path,
    /// The work the attempt charged: kept in `used` where the path is exact,
    /// and returned to the allowance and counted in `refunded` otherwise.
    pub attempt_work: u64,
    /// The whole operation's work, the attempt and the fitted construction
    /// after it: what it added to `used`, `refunded` and `uncharged`
    /// together. The quick look's `read` is part of it.
    pub operation_work: u64,
}

thread_local! {
    static LOG: RefCell<Option<Vec<Choice>>> = const { RefCell::new(None) };
}

/// Records every capped exact attempt `run` makes on this thread. Recording
/// reads the counters and changes no answer and no work.
pub fn with_log<T>(run: impl FnOnce() -> T) -> (T, Vec<Choice>) {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            LOG.with(|log| *log.borrow_mut() = None);
        }
    }
    LOG.with(|log| {
        assert!(log.borrow().is_none(), "predictor logs cannot nest");
        *log.borrow_mut() = Some(Vec::new());
    });
    let _reset = Reset;
    let result = run();
    let choices = LOG.with(|log| log.borrow_mut().take().unwrap_or_default());
    (result, choices)
}

/// Completes the record of the attempt of `step` that just ended.
pub(crate) fn record(step: &impl ToString, path: Path, attempt_work: u64, operation_work: u64) {
    LOG.with(|log| {
        if let Some(choices) = log.borrow_mut().as_mut() {
            let step = step.to_string();
            match choices
                .iter_mut()
                .rev()
                .find(|choice| choice.step == step && choice.path == Path::Open)
            {
                Some(choice) => {
                    choice.path = path;
                    choice.attempt_work = attempt_work;
                    choice.operation_work = operation_work;
                }
                None => choices.push(Choice {
                    step,
                    look: None,
                    path,
                    attempt_work,
                    operation_work,
                }),
            }
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::Rational;
    use catallax_shapes::{Plane, Point3, Vector3};

    fn plane(height: &str) -> Surface {
        let q = |text: &str| Real::from_rational(Rational::parse(text).unwrap());
        Surface::Plane(Plane {
            origin: Point3 {
                x: q("0"),
                y: q("0"),
                z: q(height),
            },
            normal: Vector3 {
                x: q("0"),
                y: q("0"),
                z: q("1"),
            },
            ref_dir: Vector3 {
                x: q("1"),
                y: q("0"),
                z: q("0"),
            },
        })
    }

    #[test]
    fn the_tier_is_set_by_degree_and_limbs_together() {
        let participant = |degree, limbs| Participant {
            degree,
            limbs,
            read: 0,
        };
        let tier = |a, b| {
            let mut look = Look::NONE;
            look.pair(a, b, false);
            look.tier()
        };
        // Planes, quadrics and tori on small rationals: D <= 8.
        assert_eq!(tier(participant(1, 2), participant(1, 2)), Tier::Low);
        assert_eq!(tier(participant(2, 2), participant(2, 2)), Tier::Low);
        assert_eq!(tier(participant(2, 2), participant(4, 2)), Tier::Low);
        // Two tori are past ESOLID's pairs, and so are two bicubic patches.
        assert_eq!(tier(participant(4, 2), participant(4, 2)), Tier::High);
        assert_eq!(tier(participant(18, 2), participant(18, 2)), Tier::High);
        // A patch against a plane is its own degree.
        assert_eq!(tier(participant(18, 2), participant(1, 2)), Tier::High);
        // Long rationals raise the cost of a low degree: 1 x 6^2 is past 32.
        assert_eq!(tier(participant(1, 5), participant(1, 5)), Tier::Low);
        assert_eq!(tier(participant(1, 6), participant(1, 2)), Tier::High);
        assert_eq!(tier(participant(2, 4), participant(1, 2)), Tier::Low);
        assert_eq!(tier(participant(2, 5), participant(1, 2)), Tier::High);
    }

    #[test]
    fn a_pair_proven_apart_adds_nothing() {
        let patch = Participant {
            degree: 18,
            limbs: 2,
            read: 0,
        };
        let mut look = Look::NONE;
        look.pair(patch, patch, true);
        assert_eq!(look, Look::NONE);
        assert_eq!(look.tier(), Tier::Low);
        look.pair(patch, patch, false);
        assert_eq!((look.degree, look.tier()), (324, Tier::High));
    }

    #[test]
    fn a_surface_is_read_by_kind_and_by_its_longest_rational() {
        let small = Participant::of(&plane("5"));
        // Nine numbers and the kind.
        assert_eq!((small.degree, small.limbs, small.read), (1, 2, 10));
        // A numerator of 130 bits takes three limbs, and the denominator one.
        let long = Participant::of(&plane("1361129467683753853853498429727072845824"));
        assert_eq!((long.degree, long.limbs), (1, 4));
    }

    #[test]
    fn a_surface_with_more_numbers_than_may_be_read_is_not_a_cheap_one() {
        // A plane is its kind and nine numbers. With five reads left the kind
        // and four numbers are read, and the plane has no degree to go by.
        let short = Participant::within(&plane("5"), 5);
        assert_eq!((short.degree, short.read), (u64::MAX, 5));
        let mut look = Look::NONE;
        look.pair(short, Participant::of(&plane("5")), false);
        assert_eq!(look.tier(), Tier::High);
        // Exactly enough reads is the whole plane.
        assert_eq!(
            Participant::within(&plane("5"), 10),
            Participant::of(&plane("5"))
        );
        // A hundredth of the charge at which a look is taken.
        assert_eq!(LOOK_READS * 100, HIGH_COMPLEXITY_CAP);
    }

    #[test]
    fn the_caps_are_the_decisions_tiers() {
        assert_eq!(Tier::Low.cap(), 20_000_000);
        assert_eq!(Tier::High.cap(), 2_000_000);
        assert_eq!(LOW_COMPLEXITY_SCORE, 32);
    }
}
