//! Distance to the short circular sections of a certified rolling-ball path.
//! This uses the interval panel cover of `rolling_region::path`, specialized
//! to the spine's stored binary64 enclosures. SMLib's
//! `SmCircularCrossSectionFSG` also uses the short circular sections between
//! the contact rails. OpenCASCADE searches for point-to-curve extrema by
//! refining sampled brackets in `Extrema_GExtPC`. Here every panel bounds
//! the distance and only a certified station can supply an upper witness.

use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;

use catallax_shapes::{RPoint3, RVector3};

use super::*;
use crate::rolling_region::Range;

pub(crate) mod loft;

// The finite path cover's existing panel and depth limits.
const PANELS: usize = 4096;
const DEPTH: usize = 64;

type V = [Interval; 3];

fn norm(v: &V) -> Interval {
    (v[0].square() + v[1].square() + v[2].square())
        .sqrt()
        .expect("sum of squares")
}

fn reserve_panel(dec: &mut Decider, products: u64) -> Result<(), Refusal> {
    dec.charge_work(products)
}

fn gap(detail: &str) -> Refusal {
    Refusal::new("field.spline-sweep", detail)
}

impl Support {
    pub(super) fn normal_products(&self, multiply: u64, scale: u64) -> u64 {
        let [m, n] = self.degrees.map(|degree| degree as u64);
        let powers = m.saturating_sub(1) + n.saturating_sub(1) + (m + 1) * n;
        let coefficients = 12 * (m + 1) * (n + 1);
        (powers + 18) * multiply + coefficients * scale
    }

    fn normal_numerator<N: Num>(&self, u: N, v: N) -> [N; 3] {
        let f = self.values(&u, &v, 3);
        let xyz = |row: usize| std::array::from_fn(|i| f[row][i].clone());
        let w = &f[0][3];
        let x = xyz(0);
        let tu = tangent(w, &xyz(1), &f[1][3], &x);
        let tv = tangent(w, &xyz(2), &f[2][3], &x);
        cross(&tu, &tv)
    }

    fn radial(&self, u: Interval, v: Interval, length: Interval) -> V {
        // P - C = -signed_radius * M / |M|.
        mul3(&(-self.contact_side / length), &self.normal_numerator(u, v))
    }

    fn radial_slope(&self, u: Directional, v: Directional, length: Directional) -> V {
        let [l, dl, _] = length.0;
        let scale = -self.contact_side / l;
        let ds = -(scale * dl / l);
        self.normal_numerator(u, v)
            .map(|d| product(scale, d.0[1]) + product(ds, d.0[0]))
    }
}

/// An oriented cap plane through `point`, with `inward` into the retained
/// half space. Its normal need not be a unit vector.
pub type Cap = (RPoint3, RVector3);

#[derive(Clone, Debug)]
struct HalfSpace {
    origin: V,
    inward: V,
    length: Interval,
}

struct Arc {
    centre: V,
    feet: [V; 2],
    normal: V,
}

impl Arc {
    fn sector(&self, v: &V) -> [Interval; 2] {
        let [a, b] = &self.feet;
        [
            dot(&cross(a, v), &self.normal),
            dot(&cross(v, b), &self.normal),
        ]
    }
}

impl HalfSpace {
    fn new((p, n): &Cap) -> Result<Self, Refusal> {
        let origin = [&p.x, &p.y, &p.z].map(Interval::from_rational);
        let inward = [&n.x, &n.y, &n.z].map(Interval::from_rational);
        let length = norm(&inward);
        if length.lo() <= 0.0 || !length.is_finite() || origin.iter().any(|x| !x.is_finite()) {
            return Err(cap_gap(
                "the cap plane has no finite nonzero normal or origin",
            ));
        }
        Ok(Self {
            origin,
            inward,
            length,
        })
    }

    fn at(&self, centre: &V) -> Interval {
        dot(&self.inward, &sub3(centre, &self.origin))
    }
}

struct ClippedArc<'a> {
    arc: &'a Arc,
    caps: &'a [HalfSpace; 2],
    query: &'a V,
    lower: f64,
    upper: f64,
}

impl ClippedArc<'_> {
    fn candidate(&mut self, radial: &V, exists: bool, foot: bool, on_cap: Option<usize>) {
        let q = add3(&self.arc.centre, radial);
        let mut certain = exists;
        if !foot {
            for side in self.arc.sector(radial) {
                if side.hi() < 0.0 {
                    return;
                }
                certain &= side.lo() >= 0.0;
            }
        }
        for (i, cap) in self.caps.iter().enumerate() {
            if on_cap == Some(i) {
                continue;
            }
            let side = cap.at(&q);
            if side.hi() < 0.0 {
                return;
            }
            certain &= side.lo() >= 0.0;
        }
        let distance = norm(&sub3(self.query, &q));
        self.lower = self.lower.min(distance.lo());
        if certain {
            self.upper = self.upper.min(distance.hi());
        }
    }
}

fn cap_gap(detail: &str) -> Refusal {
    Refusal::new("field.spline-sweep-cap", detail)
}

/// Proves the whole-support spine, then a continuation on extensions of
/// the same supports along their common edge parameter axis.
pub fn cap_continuation(
    definition: &SpineDefinition,
    dec: &mut Decider,
) -> Result<[Spine; 2], Refusal> {
    let mut plane = None;
    for (i, support) in definition.supports.iter().enumerate() {
        if flat_normal(support, dec)?.is_some() {
            plane = Some(i);
            break;
        }
    }
    let original = Spine::prove_offsets(definition, plane, dec)?;
    let axis = definition.axis;
    // Use the support extension of M11 (docs/M11-OFFSET-PROOF.md), one
    // quarter of the edge length on either end of its parameter interval.
    let quarter = Rational::one().mul_pow2(-2);
    let scale = &Rational::one() + &quarter.mul_pow2(1);
    let mut extended = definition.clone();
    for support in &mut extended.supports {
        let mut domain = support.carrier().domain().map(|d| d.map(Clone::clone));
        let margin = exact::mul(
            &exact::sub(&domain[axis][1], &domain[axis][0], dec)?,
            &quarter,
            dec,
        )?;
        domain[axis][0] = exact::sub(&domain[axis][0], &margin, dec)?;
        domain[axis][1] = exact::add(&domain[axis][1], &margin, dec)?;
        *support = support.extended(domain, dec)?;
    }
    extended.seed = exact::div(&exact::add(&definition.seed, &quarter, dec)?, &scale, dec)?;
    let extension = Spine::prove_offsets(&extended, plane, dec)?;

    // The old proved seed root, expressed in the new normalized charts,
    // must fall inside a new span's uniqueness box. The normal numerator
    // and its length each pick up one factor of the chart scale.
    let mut root = original.seed_root;
    let t = extension.system.tau();
    let old_tau = Rational::from_f64_exact(root[t].lo()).expect("finite proved seed");
    let new_tau = exact::div(&exact::add(&old_tau, &quarter, dec)?, &scale, dec)?;
    let scale = Interval::from_rational(&scale);
    let quarter = Interval::from_rational(&quarter);
    reserve_panel(dec, 2)?;
    for i in 0..2 {
        let slot = PARAMETERS + 2 * i + axis;
        root[slot] = (root[slot] + quarter) / scale;
        root[LENGTHS + i] = product(root[LENGTHS + i], scale);
    }
    root[t] = Interval::from_rational(&new_tau);
    let mut joined = false;
    for span in &extension.spans {
        if root[t].lo() < span.tau[0] || root[t].hi() > span.tau[1] {
            continue;
        }
        reserve_panel(dec, Sweep::enclosure_products(span))?;
        let s = span.variable(root[t]);
        if free(t).into_iter().enumerate().all(|(i, k)| {
            let residual = root[k] - span.model_at(i, s);
            residual.lo() > -span.radius[i] && residual.hi() < span.radius[i]
        }) {
            joined = true;
            break;
        }
    }

    if !joined {
        return Err(cap_gap(
            "the extended spine is not tied to the original seed",
        ));
    }
    Ok([original, extension])
}

/// A finite spine whose normal sections reach at least as far as the ball.
#[derive(Clone, Debug)]
pub struct Sweep<'a> {
    spine: &'a Spine,
    caps: Option<[HalfSpace; 2]>,
    tail_bound: f64,
}

impl<'a> Sweep<'a> {
    /// Checks reach and the contact sector on every span.
    pub fn new(spine: &'a Spine, dec: &mut Decider) -> Result<Self, Refusal> {
        Self::admit(spine, None, dec)
    }

    fn admit(
        spine: &'a Spine,
        caps: Option<&[Cap; 2]>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let checked = (|| {
            let caps = caps
                .map(|[a, b]| -> Result<_, Refusal> {
                    reserve_panel(dec, 6)?;
                    Ok([HalfSpace::new(a)?, HalfSpace::new(b)?])
                })
                .transpose()?;
            let mut sweep = Self {
                spine,
                caps,
                tail_bound: 0.0,
            };
            sweep.tail_bound = sweep.check(dec)?;
            Ok(sweep)
        })();

        checked
    }

    /// Maximum distance from a near-tangent section to either contact rail
    /// in panels admitted by the chord bound, in millimetres.
    pub fn tail_bound(&self) -> Option<Rational> {
        (self.tail_bound > 0.0)
            .then(|| Rational::from_f64_exact(self.tail_bound).expect("proved finite bound"))
    }

    /// Checks and clips a continuation against the two cap planes, in spine
    /// order. Each cap must cut every arc parameter once, and the two cuts
    /// must be ordered along the whole face.
    pub fn between_caps(
        spine: &'a Spine,
        caps: &[Cap; 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let sweep = Self::admit(spine, Some(caps), dec)?;
        let checked = sweep.check_caps(dec);

        checked?;
        Ok(sweep)
    }

    fn check_caps(&self, dec: &mut Decider) -> Result<(), Refusal> {
        let caps = self.caps.as_ref().expect("has caps");
        let extent = |p: &[Interval; 9], cap: &HalfSpace| {
            let reach = (self.spine.radius * cap.length).hi();
            cap.at(&[p[0], p[1], p[2]]) + spread(reach)
        };
        let [a, b] = &self.spine.ends;
        reserve_panel(dec, 8)?;
        if extent(&a.point, &caps[0]).hi() >= 0.0 || extent(&b.point, &caps[1]).hi() >= 0.0 {
            return Err(cap_gap("an end sphere is not proved outside its cap"));
        }
        let mut has_interior = false;
        for span in &self.spine.spans {
            let mid = span.tau[0].midpoint(span.tau[1]);
            let t = self.spine.system.tau();
            if mid <= a.point[t].hi() || mid >= b.point[t].lo() {
                continue;
            }
            reserve_panel(dec, Self::enclosure_products(span) + 8)?;
            let station = self.enclose(span, point(mid));
            if caps.iter().all(|cap| extent(&station, cap).lo() > 0.0) {
                has_interior = true;
                break;
            }
        }
        if !has_interior {
            return Err(cap_gap("no section sphere is proved between the caps"));
        }
        for span in &self.spine.spans {
            self.check_cap_slope(span, span.tau, 0, dec)?;
        }
        Ok(())
    }

    fn check_cap_slope(
        &self,
        span: &SpineSpan,
        [lo, hi]: [f64; 2],
        depth: usize,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let slope_products = self
            .spine
            .system
            .supports
            .iter()
            .map(|support| support.normal_products(6, 3) + 7)
            .sum::<u64>();
        reserve_panel(
            dec,
            Self::enclosure_products(span) + self.radial_products() + slope_products + 20,
        )?;
        let pinned = self.spine.system.tau();
        let d = span.slope.expect("reach proved the slope");
        let z = self.enclose(span, Interval::new(lo, hi).unwrap());
        let mut jet = z.map(|x| Directional([x, point(1.0), zero()]));
        for (i, index) in free(pinned).into_iter().enumerate() {
            jet[index].0[1] = d[i];
        }
        let [a, b] = self.radials(&z);
        let radius = self.spine.radius;
        // Parameterize the short arc by normalizing (1-u) A + u B.
        // Its shortest unnormalized vector is at u=1/2, since |A|=|B|.
        let shortest = ((radius.square() + dot(&a, &b)) / point(2.0))
            .sqrt()
            .map_or(0.0, |x| x.lo());
        let drift = (0..2)
            .map(|i| {
                let s = &self.spine.system.supports[i];
                norm(&s.radial_slope(
                    jet[3 + 2 * i].clone(),
                    jet[4 + 2 * i].clone(),
                    jet[7 + i].clone(),
                ))
                .hi()
            })
            .fold(0.0, f64::max);
        if shortest > 0.0
            && self
                .caps
                .as_ref()
                .unwrap()
                .iter()
                .enumerate()
                .all(|(i, cap)| {
                    let rate = dot(&cap.inward, &[d[0], d[1], d[2]]);
                    let variation = (cap.length * radius * point(drift) / point(shortest)).hi();
                    if i == 0 {
                        rate.lo() > variation
                    } else {
                        rate.hi() < -variation
                    }
                })
        {
            return Ok(());
        }
        let mid = lo.midpoint(hi);
        if depth >= (LAST_LEVEL - FIRST_LEVEL) as usize || mid <= lo || mid >= hi {
            return Err(cap_gap(
                "a cap height is not proved monotone along the sweep",
            ));
        }
        self.check_cap_slope(span, [lo, mid], depth + 1, dec)?;
        self.check_cap_slope(span, [mid, hi], depth + 1, dec)
    }

    fn check(&self, dec: &mut Decider) -> Result<f64, Refusal> {
        let reach = || {
            Refusal::new(
                "field.spline-sweep-reach",
                "the normal sections are not proved clear of the other balls",
            )
        };
        let mut speed: Option<V> = None;
        let mut acceleration: f64 = 0.0;
        for span in &self.spine.spans {
            reserve_panel(dec, 3)?;
            let d = span.slope.ok_or_else(reach)?;
            let d = [d[0], d[1], d[2]];
            speed = Some(speed.map_or(d, |a| std::array::from_fn(|i| a[i].hull(&d[i]))));
            let second = span.acceleration.ok_or_else(reach)?;
            acceleration = acceleration.max(norm(&second).hi());
        }
        // Federer's normal-chord reach bound (Curvature Measures, 1959,
        // Theorem 4.18). A fixed unit coordinate direction gives a global
        // chord bound m |dt|. Taylor's theorem bounds its projection onto
        // any unit normal at the station by M dt^2 / 2. Thus every such
        // normal segment of length r is clear when r M < m^2.
        let m = speed
            .ok_or_else(reach)?
            .map(|v| v.lo().max(-v.hi()).max(0.0))
            .into_iter()
            .fold(0.0, f64::max);
        reserve_panel(dec, 2)?;
        if (self.spine.radius * point(acceleration)).hi() >= point(m).square().lo() {
            return Err(reach());
        }
        self.spine.spans.iter().try_fold(0.0_f64, |bound, span| {
            Ok(bound.max(self.check_sector(span, span.tau, 0, dec)?))
        })
    }

    fn check_sector(
        &self,
        span: &SpineSpan,
        [lo, hi]: [f64; 2],
        depth: usize,
        dec: &mut Decider,
    ) -> Result<f64, Refusal> {
        let sector_products = if self.caps.is_some() { 9 } else { 15 };
        reserve_panel(
            dec,
            Self::enclosure_products(span) + self.radial_products() + sector_products,
        )?;
        let z = self.enclose(span, Interval::new(lo, hi).unwrap());
        let [a, b] = self.radials(&z);
        // The chord from either rail increases with central angle on the
        // short arc. Its maximum is the distance to the other rail.
        // This also covers a point section when the two feet coincide.
        if self.caps.is_none() && dot(&a, &b).lo() > 0.0 {
            let chord = norm(&sub3(&a, &b)).hi();
            if chord <= self.spine.target {
                return Ok(chord);
            }
        }
        if norm(&cross(&a, &b)).lo() > 0.0 {
            return Ok(0.0);
        }
        let mid = lo.midpoint(hi);
        if depth >= (LAST_LEVEL - FIRST_LEVEL) as usize || mid <= lo || mid >= hi {
            return Err(Refusal::new(
                "field.spline-sweep-sector",
                "the two contact directions are not proved transverse",
            ));
        }
        let first = self.check_sector(span, [lo, mid], depth + 1, dec)?;
        Ok(first.max(self.check_sector(span, [mid, hi], depth + 1, dec)?))
    }

    fn enclosure_products(span: &SpineSpan) -> u64 {
        (span.model.len() * span.model[0].len()) as u64
    }

    fn radial_products(&self) -> u64 {
        self.spine
            .system
            .supports
            .iter()
            .map(|support| support.normal_products(1, 1) + 3)
            .sum()
    }

    fn enclose(&self, span: &SpineSpan, tau: Interval) -> [Interval; 9] {
        let pinned = self.spine.system.tau();
        let raw = span.variable(tau);
        let s = Interval::new(raw.lo().max(-1.0), raw.hi().min(1.0)).expect("inside span");
        let mut z = [tau; 9];
        for (i, index) in free(pinned).into_iter().enumerate() {
            z[index] = span.model_at(i, s) + span.offset[i];
        }
        z
    }

    fn radials(&self, z: &[Interval; 9]) -> [V; 2] {
        std::array::from_fn(|i| {
            self.spine.system.supports[i].radial(z[3 + 2 * i], z[4 + 2 * i], z[7 + i])
        })
    }

    fn evaluate_arc(
        &self,
        p: &V,
        end: Result<&[Interval; 9], (&SpineSpan, Interval)>,
        dec: &mut Decider,
    ) -> Result<(f64, f64), Refusal> {
        let arc_products = if self.caps.is_some() { 257 } else { 54 };
        let enclosure = end
            .as_ref()
            .err()
            .map_or(0, |(span, _)| Self::enclosure_products(span));
        reserve_panel(dec, self.radial_products() + arc_products + enclosure)?;
        Ok(match end {
            Ok(point) => self.arc(p, point),
            Err((span, domain)) => self.arc(p, &self.enclose(span, domain)),
        })
    }

    // At a station these are lower/upper bounds on distance to the arc.
    // Over a panel the lower bound applies to every station in the panel.
    fn arc(&self, p: &V, z: &[Interval; 9]) -> (f64, f64) {
        let radius = self.spine.radius;
        let centre = [z[0], z[1], z[2]];
        let v = sub3(p, &centre);
        let active_caps = if let Some(caps) = &self.caps {
            let side: [Interval; 2] = std::array::from_fn(|i| {
                caps[i].at(&centre) + spread((radius * caps[i].length).hi())
            });
            if side.iter().any(|s| s.hi() < 0.0) {
                return (f64::INFINITY, f64::INFINITY);
            }
            side.iter().any(|s| s.lo() < 0.0).then_some(caps)
        } else {
            None
        };
        let [a, b] = self.radials(z);
        let feet = [norm(&sub3(&v, &a)), norm(&sub3(&v, &b))];
        let foot_lo = feet[0].lo().min(feet[1].lo());
        let foot_hi = feet[0].hi().min(feet[1].hi());
        let sphere = (norm(&v) - radius).square().sqrt().expect("square").lo();
        let chord = norm(&sub3(&a, &b));
        let rail_bound = || (point(foot_lo) - chord).lo().max(sphere);
        if self.caps.is_none() && dot(&a, &b).lo() > 0.0 && chord.hi() <= self.spine.target {
            return (rail_bound(), foot_hi);
        }

        // Both radii are perpendicular to C'. Their cross product is a
        // section normal with the orientation from rail A toward rail B.
        let n = cross(&a, &b);
        let length = norm(&n);
        if length.lo() <= 0.0 {
            return (rail_bound(), active_caps.map_or(foot_hi, |_| f64::INFINITY));
        }
        let normal = n.map(|x| x / length);
        let height = dot(&v, &normal);
        let planar = sub3(&v, &mul3(&height, &normal));
        let planar_length = norm(&planar);
        let circle = ((planar_length - radius).square() + height.square())
            .sqrt()
            .expect("sum of squares");
        let arc = Arc {
            centre,
            feet: [a, b],
            normal,
        };
        let [first, second] = arc.sector(&v);
        let (lo, hi) = if first.lo() > 0.0 && second.lo() > 0.0 {
            (circle.lo().max(sphere), circle.hi())
        } else if first.hi() < 0.0 || second.hi() < 0.0 {
            (foot_lo.max(sphere), foot_hi)
        } else {
            (circle.lo().min(foot_lo).max(sphere), foot_hi)
        };
        if let Some(caps) = active_caps {
            self.clipped_arc(p, &arc, caps, &planar, planar_length, lo)
        } else {
            (lo, hi)
        }
    }

    fn clipped_arc(
        &self,
        p: &V,
        arc: &Arc,
        caps: &[HalfSpace; 2],
        planar: &V,
        planar_length: Interval,
        lower: f64,
    ) -> (f64, f64) {
        let radius = self.spine.radius;
        let mut clipped = ClippedArc {
            arc,
            caps,
            query: p,
            lower: f64::INFINITY,
            upper: f64::INFINITY,
        };
        for foot in &arc.feet {
            clipped.candidate(foot, true, true, None);
        }
        // The minimum of distance on a clipped circle lies at the usual
        // projection, on a rail, or at a cap intersection. If the query
        // projects onto the circle's centre, the distance is constant.
        let mut complete = planar_length.lo() > 0.0;
        if complete {
            clipped.candidate(&mul3(&(radius / planar_length), planar), true, false, None);
        }
        for (i, cap) in caps.iter().enumerate() {
            let height = cap.at(&arc.centre);
            if (height + spread((radius * cap.length).hi())).lo() >= 0.0 {
                continue;
            }
            // The cap cuts a line in the section plane. Its foot is along
            // the projected cap normal, and its two circle intersections
            // lie one half-chord to either side of that foot.
            let t = &arc.normal;
            let axis = sub3(&cap.inward, &mul3(&dot(&cap.inward, t), t));
            let length = norm(&axis);
            if length.lo() <= 0.0 {
                complete = false;
                continue;
            }
            let unit = axis.map(|x| x / length);
            let offset = -height / length;
            let radicand = radius.square() - offset.square();
            if let Some(half_chord) = radicand.sqrt() {
                let base = mul3(&offset, &unit);
                let shift = mul3(&half_chord, &cross(t, &unit));
                for radial in [add3(&base, &shift), sub3(&base, &shift)] {
                    clipped.candidate(&radial, radicand.lo() >= 0.0, false, Some(i));
                }
            }
        }
        (
            if complete {
                lower.max(clipped.lower)
            } else {
                lower
            },
            clipped.upper,
        )
    }

    /// Propose the perpendicular centre station using the stored model.
    /// This is the stationary equation used by OCCT's Extrema_GExtPC.
    /// Its answer only contributes if the interval arc supplies a witness.
    fn guess(&self, p: &V, dec: &mut Decider) -> Result<f64, Refusal> {
        reserve_panel(dec, 7)?;
        let [a, b] = self.spine.ends.each_ref().map(|e| &e.point);
        let t = self.spine.system.tau();
        let direction: [f64; 3] = std::array::from_fn(|i| b[i].mid() - a[i].mid());
        let length = direction.iter().map(|x| x * x).sum::<f64>();
        let offset = (0..3)
            .map(|i| (p[i].mid() - a[i].mid()) * direction[i])
            .sum::<f64>();

        let start = a[t].hi();
        let end = b[t].lo();
        let mut tau = (start + (end - start) * (offset / length)).clamp(start, end);

        let eval = |coefficients: &[f64], s: f64| {
            let mut value = *coefficients.last().unwrap();
            let mut first = 0.0;
            let mut second = 0.0;
            for c in coefficients[..coefficients.len() - 1].iter().rev() {
                second = second * s + 2.0 * first;
                first = first * s + value;
                value = value * s + c;
            }
            [value, first, second]
        };
        for _ in 0..NEWTON_STEPS {
            let Some(span) = self.span_at(tau) else {
                break;
            };
            reserve_panel(dec, (12 * (TERMS - 1) + 10) as u64)?;
            let half = (span.tau[1] - span.tau[0]) / 2.0;
            let s = (tau - span.tau[0].midpoint(span.tau[1])) / half;
            let mut f = 0.0;
            let mut df = 0.0;
            for (i, query) in p.iter().enumerate() {
                let [x, first, second] = eval(&span.model[i], s);
                let gap = x - query.mid();
                f += gap * first;
                df += first * first + gap * second;
            }

            let next = (tau - half * f / df).clamp(start, end);
            if !next.is_finite() || (next - tau).abs() < NEWTON_STOP {
                break;
            }
            tau = next;
        }
        Ok(tau)
    }

    fn span_at(&self, tau: f64) -> Option<&SpineSpan> {
        let spans = &self.spine.spans;
        spans
            .get(spans.partition_point(|span| span.tau[1] < tau))
            .filter(|span| span.tau[0] <= tau && tau <= span.tau[1])
    }

    /// A guaranteed unsigned distance range in millimetres to this sweep's
    /// short arcs, including the two certified end arcs.
    pub fn distance(
        &self,
        query: &RPoint3,
        width: &Rational,
        dec: &mut Decider,
    ) -> Result<Range, Refusal> {
        self.distance_range(query, width, false, dec)
    }

    /// Refines the distance range until it proves or disproves the given bound.
    /// A range still straddles the bound if the refinement limits are reached.
    pub fn distance_at_most(
        &self,
        query: &RPoint3,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<Range, Refusal> {
        self.distance_range(query, bound, true, dec)
    }

    fn distance_range(
        &self,
        query: &RPoint3,
        limit: &Rational,
        threshold: bool,
        dec: &mut Decider,
    ) -> Result<Range, Refusal> {
        if limit.sign() == Sign::Negative {
            return Err(gap("requested distance limit must be non-negative"));
        }
        let p = [&query.x, &query.y, &query.z].map(Interval::from_rational);
        if p.iter().any(|v| !v.is_finite()) {
            return Err(gap("the query is outside the interval coordinate range"));
        }
        let tau = self.spine.system.tau();
        let [a, b] = &self.spine.ends;
        let mut upper = self
            .evaluate_arc(&p, Ok(&a.point), dec)?
            .1
            .min(self.evaluate_arc(&p, Ok(&b.point), dec)?.1);
        let limit = Interval::from_rational(limit);
        let guess = self.guess(&p, dec)?;
        if a.point[tau].hi() <= guess
            && guess <= b.point[tau].lo()
            && let Some(span) = self.span_at(guess)
        {
            upper = upper.min(self.evaluate_arc(&p, Err((span, point(guess))), dec)?.1);
        }
        if upper.is_finite() && upper <= limit.lo() {
            return Ok(Range {
                lower: Rational::zero(),
                upper: Rational::from_f64_exact(upper),
                converged: true,
            });
        }
        let mut pending = BinaryHeap::new();
        for (index, span) in self.spine.spans.iter().enumerate() {
            let lo = span.tau[0].max(a.point[tau].lo());
            let hi = span.tau[1].min(b.point[tau].hi());
            if let Some(domain) = Interval::new(lo, hi) {
                let (lower, _) = self.evaluate_arc(&p, Err((span, domain)), dec)?;
                if lower.is_infinite() {
                    continue;
                }
                pending.push(Reverse(Panel {
                    lower,
                    lo,
                    hi,
                    index,
                    depth: 0,
                }));
            }
        }
        let mut panels = 0;
        while let Some(Reverse(Panel {
            lower,
            lo,
            hi,
            index,
            depth,
        })) = pending.pop()
        {
            let lower = lower.min(upper);
            let converged = if threshold {
                (upper.is_finite() && upper <= limit.lo()) || lower > limit.hi()
            } else {
                upper.is_finite() && (point(upper) - point(lower)).hi() <= limit.lo()
            };
            let mid = lo.midpoint(hi);
            if converged || panels >= PANELS || depth >= DEPTH || mid <= lo || mid >= hi {
                return Ok(Range {
                    lower: Rational::from_f64_exact(lower)
                        .ok_or_else(|| gap("nonfinite lower bound"))?,
                    upper: Rational::from_f64_exact(upper),
                    converged,
                });
            }
            panels += 1;
            let span = &self.spine.spans[index];
            // The outer interval bounds every actual station. A witness
            // needs a parameter in the inner interval or a certified end.
            if a.point[tau].hi() <= mid && mid <= b.point[tau].lo() {
                upper = upper.min(self.evaluate_arc(&p, Err((span, point(mid))), dec)?.1);
            }
            for (lo, hi) in [(lo, mid), (mid, hi)] {
                let domain = Interval::new(lo, hi).expect("ordered panel");
                let (lower, _) = self.evaluate_arc(&p, Err((span, domain)), dec)?;
                if lower.is_infinite() {
                    continue;
                }
                pending.push(Reverse(Panel {
                    lower,
                    lo,
                    hi,
                    index,
                    depth: depth + 1,
                }));
            }
        }
        Err(gap("the sweep has no stations"))
    }
}

struct Panel {
    lower: f64,
    lo: f64,
    hi: f64,
    index: usize,
    depth: usize,
}

impl Ord for Panel {
    fn cmp(&self, other: &Self) -> Ordering {
        self.lower
            .total_cmp(&other.lower)
            .then(self.index.cmp(&other.index))
            .then(self.lo.total_cmp(&other.lo))
            .then(self.hi.total_cmp(&other.hi))
            .then(self.depth.cmp(&other.depth))
    }
}

impl PartialOrd for Panel {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Panel {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}

impl Eq for Panel {}
