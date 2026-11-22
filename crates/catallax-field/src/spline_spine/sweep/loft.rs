//! A convex ruled loft with one rounded edge. The cap cuts give a
//! horizontal arc in each slice. The arc replaces the corner between its
//! two rails, and the remaining sides stay on the original loft.

use catallax_recipe::Origin;

use crate::Location;

use super::*;

#[cfg(test)]
thread_local! {
    static PREVIEW_SECTIONS: Cell<usize> = const { Cell::new(0) };
}

fn vector(v: &RVector3) -> V {
    [&v.x, &v.y, &v.z].map(Interval::from_rational)
}

fn xyz(p: &RPoint3) -> V {
    [&p.x, &p.y, &p.z].map(Interval::from_rational)
}

fn clearance() -> Refusal {
    Refusal::new(
        "field.spline-sweep-clearance",
        "the rolling ball is not proved to clear the loft's side faces",
    )
}

/// The side inequality of a convex loft in a section of normalized height
/// h. The horizontal outward normal and the point on its side are affine
/// functions of h. In world coordinates the inequality is quadratic.
#[derive(Clone, Debug)]
struct Side {
    start: RPoint3,
    shift: RVector3,
    normal: RVector3,
    turn: RVector3,
}

impl Side {
    fn value(&self, p: &RPoint3, h: &Rational) -> Sign {
        let start = self.start.plus(&self.shift.scaled(h));
        self.normal
            .plus(&self.turn.scaled(h))
            .dot(&p.minus(&start))
            .sign()
    }

    fn project(&self, p: &RPoint3, h: &Rational) -> RPoint3 {
        let start = self.start.plus(&self.shift.scaled(h));
        let normal = self.normal.plus(&self.turn.scaled(h));
        p.plus(&normal.scaled(&(-normal.dot(&p.minus(&start)) / normal.norm_squared())))
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Round {
    pub(crate) origin: Origin,
    pub(crate) selected: [usize; 2],
    spine: Spine,
    caps: [HalfSpace; 2],
    sides: Vec<Side>,
    bottom: RPoint3,
    height: RVector3,
    /// A horizontal ray towards the original loft corner, and the remaining
    /// horizontal coordinate that completes the projection.
    outward: V,
    across: V,
}

impl Round {
    /// The sides are in loop order, and their two endpoints lie on the two
    /// parallel cap planes. The producer has proved all intermediate
    /// sections are strictly convex.
    pub(crate) fn prove(
        origin: Origin,
        selected: [usize; 2],
        patches: &[(&RationalPatch, bool)],
        caps: &[Cap; 2],
        radius: Rational,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let definition = SpineDefinition {
            supports: selected.map(|i| patches[i].0.clone()),
            radius,
            sides: selected.map(|i| {
                if patches[i].1 {
                    Sign::Negative
                } else {
                    Sign::Positive
                }
            }),
            axis: 1,
            seed: Rational::one().mul_pow2(-1),
            // The M11 fitted face's 1 um maximum deviation.
            tolerance: Rational::from_ratio(1, 1000).expect("nonzero"),
        };
        let [_, spine] = cap_continuation(&definition, dec)?;
        let sweep = Sweep::between_caps(&spine, caps, dec)?;
        let (bottom, up) = &caps[0];
        let height = up.scaled(&(&Rational::one() / &up.dot(&caps[1].0.minus(bottom))));
        let sides = patches
            .iter()
            .map(|(patch, forward)| {
                let controls: Vec<_> = patch
                    .carrier()
                    .controls()
                    .iter()
                    .map(catallax_shapes::Point3::as_rational)
                    .collect::<Option<_>>()
                    .ok_or_else(clearance)?;
                let [a, b, c, d] = controls.as_slice() else {
                    return Err(clearance());
                };
                let normal = |edge: RVector3| {
                    if *forward {
                        edge.cross(up)
                    } else {
                        up.cross(&edge)
                    }
                };
                let n = normal(c.minus(a));
                Ok(Side {
                    start: a.clone(),
                    shift: b.minus(a),
                    turn: normal(d.minus(b)).minus(&n),
                    normal: n,
                })
            })
            .collect::<Result<Vec<_>, _>>()?;

        reserve_panel(dec, sweep.radial_products())?;
        let [a, b] = sweep.radials(&spine.seed_root);
        let proposal =
            add3(&a, &b).map(|x| Rational::from_f64_exact(x.mid()).expect("finite seed"));
        let [x, y, z] = proposal;
        let outward = up.cross(&RVector3::new(x, y, z).cross(up));
        let across = vector(&up.cross(&outward));
        let outward = vector(&outward);
        let prepared_caps = sweep.caps.clone().expect("checked caps");
        let round = Self {
            origin,
            selected,
            spine,
            caps: prepared_caps,
            sides,
            bottom: bottom.clone(),
            height,
            outward,
            across,
        };
        let checked = round.check(dec);
        PRODUCTS.with(|c| c.set(0));
        checked?;
        Ok(round)
    }

    fn sweep(&self) -> Sweep<'_> {
        Sweep {
            spine: &self.spine,
            caps: Some(self.caps.clone()),
            tail_bound: 0.0,
        }
    }

    fn check(&self, dec: &mut Decider) -> Result<(), Refusal> {
        let sweep = self.sweep();
        let height = vector(&self.height);
        let bottom = xyz(&self.bottom);
        let r = self.spine.radius;
        for span in &self.spine.spans {
            let mut pending = vec![(span.tau, 0)];
            while let Some(([lo, hi], depth)) = pending.pop() {
                reserve_panel(
                    dec,
                    Sweep::enclosure_products(span)
                        + sweep.radial_products()
                        + 12
                        + 33 * self.sides.len() as u64,
                )?;
                let z = sweep.enclose(span, Interval::new(lo, hi).expect("ordered panel"));
                let c = [z[0], z[1], z[2]];
                let radial = sweep.radials(&z);
                let h = dot(&height, &sub3(&c, &bottom));
                let h_radius = (norm(&height) * r).hi();
                if h.hi() < -h_radius || h.lo() > (point(1.0) + point(h_radius)).hi() {
                    continue;
                }
                let projected = radial.iter().all(|v| dot(&self.outward, v).lo() > 0.0);
                let rails_inside = [z[3], z[5]].iter().all(|u| u.lo() > 0.0 && u.hi() < 1.0);
                let clear = self.sides.iter().enumerate().all(|(i, side)| {
                    let turn = vector(&side.turn);
                    let normal = |h| add3(&vector(&side.normal), &mul3(&h, &turn));
                    let shift = vector(&side.shift);
                    let start = xyz(&side.start);
                    // Exact quadratic Taylor remainder of F(x) =
                    // n(h) . (x - start(h)), with h affine in x.
                    let remainder = norm(&turn) * norm(&height)
                        + dot(&turn, &shift).square().sqrt().unwrap() * norm(&height).square();
                    if let Some(rail) = self.selected.iter().position(|&j| i == j) {
                        let foot = add3(&c, &radial[rail]);
                        let foot_h = dot(&height, &sub3(&foot, &bottom));
                        // At the contact, grad F is outward and parallel
                        // to the radius. On the ball sphere its linear
                        // term is -|grad F| |delta|^2 / (2 r).
                        (point(2.0) * r * remainder).hi() < norm(&normal(foot_h)).lo()
                    } else {
                        let x = sub3(&sub3(&c, &start), &mul3(&h, &shift));
                        let n = normal(h);
                        let grad = add3(&n, &mul3(&(dot(&turn, &x) - dot(&n, &shift)), &height));
                        (dot(&n, &x) + norm(&grad) * r + remainder * r.square()).hi() < 0.0
                    }
                });
                if projected && rails_inside && clear {
                    continue;
                }
                let mid = lo.midpoint(hi);
                if depth >= (LAST_LEVEL - FIRST_LEVEL) as usize || mid <= lo || mid >= hi {
                    return Err(clearance());
                }
                pending.extend([([mid, hi], depth + 1), ([lo, mid], depth + 1)]);
            }
        }
        Ok(())
    }

    pub(crate) fn distance(
        &self,
        p: &RPoint3,
        width: &Rational,
        dec: &mut Decider,
    ) -> Result<Range, Refusal> {
        self.sweep().distance(p, width, dec)
    }

    pub(crate) fn distance_at_most(
        &self,
        p: &RPoint3,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<Range, Refusal> {
        self.sweep().distance_at_most(p, bound, dec)
    }

    pub(crate) fn preview(&self, query: [f64; 3]) -> (bool, f64) {
        PRODUCTS.with(|products| {
            let saved = products.get();
            let p = query.map(point);
            let sweep = self.sweep();
            let max_height = self.caps[0].at(&self.caps[1].origin).mid();
            let height = self.caps[0].at(&p).mid().clamp(0.0, max_height);
            let [a, b] = &self.spine.ends;
            let t = self.spine.system.tau();
            let (start, end) = (a.point[t].lo(), b.point[t].hi());
            let sample = |station: f64, angle: f64| {
                #[cfg(test)]
                PREVIEW_SECTIONS.with(|count| count.set(count.get() + 1));
                let tau = (start + (end - start) * station).clamp(start, end);
                self.section_at(&sweep, tau, point(angle))
            };
            let coordinates = |q: &V| [self.caps[0].at(q).mid(), dot(&self.across, q).mid()];
            let x = coordinates(&p)[1];
            let rail = |angle| {
                let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
                for _ in 0..24 {
                    let mid = lo.midpoint(hi);
                    if coordinates(&sample(mid, angle))[0] < height {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let station = lo.midpoint(hi);
                let q = sample(station, angle);
                (coordinates(&q)[1], [station, angle], q)
            };
            let (mut left, mut right) = (rail(0.0), rail(1.0));
            if left.0 > right.0 {
                std::mem::swap(&mut left, &mut right);
            }
            let (in_fan, nearest) = if x <= left.0 {
                (false, left.2)
            } else if x >= right.0 {
                (false, right.2)
            } else {
                let fraction = (x - left.0) / (right.0 - left.0);
                let mut uv =
                    std::array::from_fn(|i| left.1[i] + fraction * (right.1[i] - left.1[i]));
                let delta = |v| {
                    if v <= 0.5 {
                        1.0 / 65_536.0
                    } else {
                        -1.0 / 65_536.0
                    }
                };
                for _ in 0..12 {
                    let [u, v] = uv;
                    let [h, across] = coordinates(&sample(u, v));
                    let (du, dv) = (delta(u), delta(v));
                    let [hu, xu] = coordinates(&sample(u + du, v));
                    let [hv, xv] = coordinates(&sample(u, v + dv));
                    let (hu, xu) = ((hu - h) / du, (xu - across) / du);
                    let (hv, xv) = ((hv - h) / dv, (xv - across) / dv);
                    let det = hu * xv - hv * xu;
                    let (dh, dx) = (height - h, x - across);
                    let step = [(dh * xv - hv * dx) / det, (hu * dx - dh * xu) / det];
                    if step.iter().any(|d| !d.is_finite()) {
                        break;
                    }
                    uv = std::array::from_fn(|i| (uv[i] + step[i]).clamp(0.0, 1.0));
                    if step.iter().all(|d| d.abs() <= 1.0 / 16_777_216.0) {
                        break;
                    }
                }
                (true, sample(uv[0], uv[1]))
            };
            let cut = in_fan && dot(&self.outward, &sub3(&p, &nearest)).mid() > 0.0;
            let h = dot(&vector(&self.height), &sub3(&p, &xyz(&self.bottom)));
            let outside = cut
                || self.caps.iter().any(|cap| cap.at(&p).mid() < 0.0)
                || self.sides.iter().any(|side| {
                    let start = add3(&xyz(&side.start), &mul3(&h, &vector(&side.shift)));
                    let normal = add3(&vector(&side.normal), &mul3(&h, &vector(&side.turn)));
                    dot(&normal, &sub3(&p, &start)).mid() > 0.0
                });
            let magnitude = query
                .into_iter()
                .zip(nearest)
                .map(|(coordinate, sample)| (coordinate - sample.mid()).powi(2))
                .sum::<f64>()
                .sqrt();
            products.set(saved);
            (cut, if outside { magnitude } else { -magnitude })
        })
    }

    pub(crate) fn boxes(&self) -> impl Iterator<Item = V> + '_ {
        let radius = spread(self.spine.radius.hi());
        self.spine.spans.iter().map(move |span| {
            let s = span.ranges();
            [s[0] + radius, s[1] + radius, s[2] + radius]
        })
    }

    pub(crate) fn side_witnesses(&self, p: &RPoint3) -> impl Iterator<Item = RPoint3> + '_ {
        let height = self.height.dot(&p.minus(&self.bottom));
        let p = p.clone();
        let inside = Rational::zero() <= height && height <= Rational::one();
        self.sides
            .iter()
            .take(if inside { self.sides.len() } else { 0 })
            .map(move |side| side.project(&p, &height))
    }

    /// Sections of the sphere envelope are graphs in the outward ray
    /// direction. The projected Jacobian has the sign of the normal's
    /// outward component (the normal is radial). The two rails have strict
    /// interior parameters on adjacent polygon sides, and every other side
    /// clears the ball. Hence the projection runs monotonically between
    /// the two rail projections, and cuts only the selected corner.
    pub(crate) fn location(&self, q: &RPoint3, dec: &mut Decider) -> Result<Location, Refusal> {
        let d = q.minus(&self.bottom);
        let height = self.height.dot(&d);
        if height < Rational::zero() || height > Rational::one() {
            return Ok(Location::Outside);
        }
        let mut boundary = height.is_zero() || height == Rational::one();
        for side in &self.sides {
            match side.value(q, &height) {
                Sign::Positive => return Ok(Location::Outside),
                Sign::Zero => boundary = true,
                Sign::Negative => (),
            }
        }
        let kept = if boundary {
            Location::Boundary
        } else {
            Location::Inside
        };
        PRODUCTS.with(|c| c.set(0));
        let p = xyz(q);
        let result = self.clip_location(&p, kept, dec);
        charge(dec)?;
        result
    }

    fn clip_location(&self, p: &V, kept: Location, dec: &mut Decider) -> Result<Location, Refusal> {
        let sweep = self.sweep();
        let height = self.caps[0].at(p);
        let side = dot(&self.across, p);
        let first = dot(
            &self.across,
            &self.at_height(&sweep, height, point(0.0), dec)?,
        );
        let second = dot(
            &self.across,
            &self.at_height(&sweep, height, point(1.0), dec)?,
        );
        let (low, high, increasing) = if first.mid() < second.mid() {
            (first, second, true)
        } else {
            (second, first, false)
        };
        if side.hi() < low.lo() || side.lo() > high.hi() {
            return Ok(kept);
        }
        let projection_inside = side.lo() > low.hi() && side.hi() < high.lo();
        let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
        let depth = |v: &V| dot(&self.outward, &sub3(v, p));
        for _ in 0..DEPTH {
            let mid = lo.midpoint(hi);
            let across = Interval::new(lo, hi).expect("ordered section");
            let gap = depth(&self.at_height(&sweep, height, across, dec)?);
            if gap.lo() > 0.0 {
                return Ok(kept);
            }
            if projection_inside && gap.hi() < 0.0 {
                return Ok(Location::Outside);
            }
            let old = (lo, hi);
            for s in [mid, lo.midpoint(mid), mid.midpoint(hi)] {
                if s <= lo || s >= hi {
                    continue;
                }
                let g = dot(
                    &self.across,
                    &self.at_height(&sweep, height, point(s), dec)?,
                );
                if (increasing && g.hi() < side.lo()) || (!increasing && g.lo() > side.hi()) {
                    lo = s;
                } else if (increasing && g.lo() > side.hi()) || (!increasing && g.hi() < side.lo())
                {
                    hi = s;
                }
            }
            if lo <= old.0 && hi >= old.1 {
                break;
            }
        }
        Ok(Location::Indeterminate)
    }

    /// Enclose the section point at the given height for every angle in
    /// the panel. Lower and upper parameter brackets move only after
    /// interval evaluations show they are on the right side of every root.
    fn at_height(
        &self,
        sweep: &Sweep<'_>,
        height: Interval,
        angle: Interval,
        dec: &mut Decider,
    ) -> Result<V, Refusal> {
        let t = self.spine.system.tau();
        let [a, b] = &self.spine.ends;
        let (mut lo, mut hi) = (a.point[t].lo(), b.point[t].hi());
        let point_at = |tau| self.section_at(sweep, tau, angle);
        for _ in 0..DEPTH {
            let mid = lo.midpoint(hi);
            if mid <= lo || mid >= hi {
                break;
            }
            let value = self.caps[0].at(&point_at(mid));
            charge(dec)?;
            if value.hi() < height.lo() {
                lo = mid;
            } else if value.lo() > height.hi() {
                hi = mid;
            } else {
                let left = lo.midpoint(mid);
                let right = mid.midpoint(hi);
                let next_lo = if self.caps[0].at(&point_at(left)).hi() < height.lo() {
                    left
                } else {
                    lo
                };
                let next_hi = if self.caps[0].at(&point_at(right)).lo() > height.hi() {
                    right
                } else {
                    hi
                };
                charge(dec)?;
                if next_lo <= lo && next_hi >= hi {
                    break;
                }
                (lo, hi) = (next_lo, next_hi);
            }
        }
        let mut bound: Option<V> = None;
        for span in &self.spine.spans {
            if let Some(panel) = Interval::new(lo.max(span.tau[0]), hi.min(span.tau[1])) {
                let section = self.section(sweep, &sweep.enclose(span, panel), angle);
                bound = Some(match bound {
                    Some(old) => std::array::from_fn(|i| old[i].hull(&section[i])),
                    None => section,
                });
            }
        }
        Ok(bound.expect("the spans cover every station"))
    }

    fn section_at(&self, sweep: &Sweep<'_>, tau: f64, angle: Interval) -> V {
        let span = sweep.span_at(tau).expect("proved continuation");
        self.section(sweep, &sweep.enclose(span, point(tau)), angle)
    }

    fn section(&self, sweep: &Sweep<'_>, z: &[Interval; 9], angle: Interval) -> V {
        let [a, b] = sweep.radials(z);
        let c = [z[0], z[1], z[2]];
        // The endpoint vectors already have length r. For an open angular
        // panel the linear combination is enclosed by the outward operations.
        if is_zero(angle) {
            return add3(&c, &a);
        }
        if angle == point(1.0) {
            return add3(&c, &b);
        }
        let w = add3(&a, &mul3(&angle, &sub3(&b, &a)));
        let radius = self.spine.radius;
        let shortest = ((radius.square() + dot(&a, &b)) / point(2.0))
            .sqrt()
            .map_or(0.0, |x| x.lo());
        let length = norm(&w);
        let length = Interval::new(length.lo().max(shortest), length.hi().min(radius.hi()))
            .unwrap_or_else(Interval::entire);
        add3(&c, &mul3(&(radius / length), &w))
    }
}

#[cfg(test)]
mod tests {
    use super::PREVIEW_SECTIONS;
    use crate::{Plan, Workspace};
    use catallax_features::ImportFiles;
    use catallax_numbers::Decider;
    use catallax_recipe::Recipe;
    use serde_json::json;

    #[test]
    fn preview_between_the_m27_rails_bounds_section_evaluations() {
        let section = |id, prefix, height, points: [[i64; 2]; 4]| {
            let edges: Vec<_> = points.into_iter().enumerate().map(|(i, [x, y])|
                json!({"kind":"line","id":format!("{prefix}_{i}"),"to":[x.to_string(),y.to_string()]}))
                .collect();
            json!({"id":id,"op":"sketch","loops":[edges],
                "plane":{"rotation":["1","0","0","0"],"translation":["0","0",height]}})
        };
        let recipe = Recipe::parse(
            &json!({
                "format":"catallax-recipe/5","name":"M27 preview","result":"round","steps":[
                    section("bottom", "a", "0", [[1,-1],[1,1],[-1,1],[-1,-1]]),
                    section("top", "b", "2", [[1,0],[0,2],[-2,0],[0,-1]]),
                    {"id":"loft","op":"loft","shape":"ruled","sections":[
                        {"profile":"bottom","start":["a_0"]},{"profile":"top","start":["b_0"]}]},
                    {"id":"round","op":"fillet","body":"loft","radius":"1/10","edges":[{"between":[
                        {"step":"loft","role":{"side":{"segment":"a_0"}}},
                        {"step":"loft","role":{"side":{"segment":"a_1"}}}
                    ]}]}
                ]
            })
            .to_string(),
        )
        .unwrap();
        let plan =
            Plan::compile_for_comparison(&recipe, &ImportFiles::new(), &mut Decider::default())
                .unwrap();
        PREVIEW_SECTIONS.with(|count| count.set(0));
        let value = plan
            .fast([0.999, -0.499, 1.0], &mut Workspace::default())
            .unwrap()
            .value;
        assert!(value > 0.01, "{value}");
        let evaluations = PREVIEW_SECTIONS.with(|count| count.get());
        assert!(
            (51..=435).contains(&evaluations),
            "{evaluations} section evaluations"
        );
    }
}
