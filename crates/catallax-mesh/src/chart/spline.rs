//! A piecewise rational face chart for the constrained triangulation.

use std::{cell::RefCell, collections::BTreeMap};

use catallax_numbers::{Decider, Interval, Rational, Refusal};
use catallax_shapes::{
    BSplineSurface,
    tessellation::{
        RationalSpan, chord_bound_interval, rational_spans, triangle_chord_bound_interval,
    },
};

use super::{Uv, failure};

pub(crate) struct SplineChart {
    spans: Vec<Span>,
    seams: Vec<Seam>,
    // A chart lives for one mesh request. The edge plan uses a projection
    // goal at least as strict as the face that consumes it.
    edge_parameters: RefCell<BTreeMap<[u64; 3], Uv>>,
}

struct Span {
    rational: RationalSpan,
    bounds: [[f64; 2]; 2],
    width: [Interval; 2],
    affine_hessian: [Interval; 3],
    local_affine_hessian: Option<[Interval; 3]>,
    normal: [Interval; 3],
    normal_cells: RefCell<BTreeMap<[u64; 4], [Interval; 3]>>,
}

struct Seam {
    axis: usize,
    at: f64,
    transverse: [f64; 2],
    jump: Interval,
}

fn point(value: f64) -> Interval {
    Interval::point(value).expect("finite chart coordinate")
}

impl Span {
    fn new(rational: RationalSpan, dec: &mut Decider) -> Result<Self, Refusal> {
        let bounds = rational
            .domain()
            .each_ref()
            .map(|range| range.each_ref().map(Rational::to_f64_nearest));
        if bounds
            .iter()
            .any(|[lo, hi]| !lo.is_finite() || !hi.is_finite() || lo >= hi)
        {
            return Err(failure(
                "spline knot spans need distinct finite parameter endpoints",
            ));
        }
        let width = bounds.map(|[lo, hi]| point(hi) - point(lo));
        let [u, v] = width;
        let local = rational.affine_hessian();
        let local_affine_hessian = local
            .iter()
            .zip(rational.chord_coefficients())
            .any(|(affine, homogeneous)| affine.hi() < homogeneous.hi())
            .then_some(local);
        let affine_hessian = [
            local[0] / u.square(),
            local[1] / (u * v),
            local[2] / v.square(),
        ];
        let normal = rational.normal_box([Interval::new(0.0, 1.0).unwrap(); 2], dec)?;
        Ok(Self {
            rational,
            bounds,
            width,
            affine_hessian,
            local_affine_hessian,
            normal,
            normal_cells: RefCell::default(),
        })
    }

    fn contains(&self, uv: Uv) -> bool {
        self.bounds
            .iter()
            .zip(uv)
            .all(|([lo, hi], x)| (*lo..=*hi).contains(&x))
    }

    fn local_seed(&self, uv: Uv) -> Uv {
        std::array::from_fn(|i| ((uv[i] - self.bounds[i][0]) / self.width[i].mid()).clamp(0.0, 1.0))
    }

    fn global(&self, uv: Uv) -> Uv {
        std::array::from_fn(|i| {
            let [lo, hi] = self.bounds[i];
            if uv[i] <= 0.0 {
                lo
            } else if uv[i] >= 1.0 {
                hi
            } else {
                (lo + uv[i] * self.width[i].mid()).clamp(lo, hi)
            }
        })
    }

    fn local_box(&self, bounds: [[f64; 2]; 2]) -> [Interval; 2] {
        std::array::from_fn(|i| {
            let [lo, hi] = bounds[i];
            let value = (Interval::new(lo, hi).unwrap() - point(self.bounds[i][0])) / self.width[i];
            Interval::new(value.lo().max(0.0), value.hi().min(1.0)).unwrap()
        })
    }

    fn distance_score(&self, position: [f64; 3]) -> f64 {
        self.rational
            .coordinate_bounds()
            .into_iter()
            .zip(position)
            .fold(0.0_f64, |score, (bounds, x)| {
                score.max((bounds.lo() - x).max(x - bounds.hi()).next_down())
            })
    }

    fn chord_bound(
        &self,
        local: [Interval; 2],
        edges: Option<[[Interval; 2]; 3]>,
        dec: &mut Decider,
    ) -> Result<Interval, Refusal> {
        let mut bound = |hessian| match edges {
            Some(edges) => {
                // Box formula, diagonal metric and three squared edge lengths.
                dec.charge_work(6 + 2 + 3 * 4)?;
                Ok::<_, Refusal>(triangle_chord_bound_interval(edges, local, &hessian))
            }
            None => {
                dec.charge_work(6)?;
                Ok(chord_bound_interval(local, &hessian))
            }
        };
        let homogeneous = bound(self.rational.chord_coefficients())?;
        if let Some(hessian) = self.local_affine_hessian {
            let affine = bound(hessian)?;
            if affine.hi() < homogeneous.hi() {
                return Ok(affine);
            }
        }
        Ok(homogeneous)
    }

    fn normal_cell(
        &self,
        local: [Interval; 2],
        dec: &mut Decider,
    ) -> Result<[Interval; 3], Refusal> {
        if local == [Interval::new(0.0, 1.0).unwrap(); 2] {
            return Ok(self.normal);
        }
        let key = [
            local[0].lo().to_bits(),
            local[0].hi().to_bits(),
            local[1].lo().to_bits(),
            local[1].hi().to_bits(),
        ];
        if let Some(&normal) = self.normal_cells.borrow().get(&key) {
            return Ok(normal);
        }
        let normal = self.rational.normal_box(local, dec)?;
        self.normal_cells.borrow_mut().insert(key, normal);
        Ok(normal)
    }
}

/// A cover of at most two adjacent cells at the range's dyadic scale.
fn dyadic_cell(range: Interval) -> Interval {
    let width = range.width();
    // Below the spacing at the unit endpoint, there is no useful common grid.
    if width < f64::EPSILON {
        return range;
    }
    let (_, exponent) = libm::frexp(width);
    let step = libm::ldexp(1.0, exponent).min(1.0);
    let lo = (range.lo() / step).floor() * step;
    let hi = (range.hi() / step).ceil() * step;
    Interval::new(lo.min(range.lo()), hi.max(range.hi())).unwrap()
}

impl SplineChart {
    pub(super) fn new(surface: &BSplineSurface, dec: &mut Decider) -> Result<Self, Refusal> {
        // Internal knot multiplicity at most the degree gives a continuous
        // homogeneous boundary. Higher multiplicity needs separate charts.
        for axis in 0..2 {
            let domain = surface.domain()[axis];
            for repeated in surface.knots()[axis].chunk_by(|a, b| a == b) {
                if repeated.len() > surface.degrees()[axis]
                    && &repeated[0] > domain[0]
                    && &repeated[0] < domain[1]
                {
                    return Err(failure("a discontinuous spline knot needs a chart split"));
                }
            }
        }
        let spans: Vec<_> = rational_spans(surface, dec)?
            .into_iter()
            .map(|span| Span::new(span, dec))
            .collect::<Result<_, _>>()?;
        let mut seams = Vec::new();
        for span in &spans {
            let domain = span.rational.domain();
            for axis in 0..2 {
                let other = 1 - axis;
                let Some(next) = spans.iter().find(|next| {
                    next.rational.domain()[axis][0] == domain[axis][1]
                        && next.rational.domain()[other] == domain[other]
                }) else {
                    continue;
                };
                let jump = span.rational.first_jump(
                    &next.rational,
                    axis,
                    [span.width[axis], next.width[axis]],
                    dec,
                )?;
                seams.push(Seam {
                    axis,
                    at: span.bounds[axis][1],
                    transverse: span.bounds[other],
                    jump,
                });
            }
        }
        Ok(Self {
            spans,
            seams,
            edge_parameters: RefCell::default(),
        })
    }

    pub(crate) fn project(
        &self,
        position: [f64; 3],
        seed: Option<Uv>,
        goal: f64,
        dec: &mut Decider,
    ) -> Result<Uv, Refusal> {
        let key = position.map(f64::to_bits);
        if let Some(&uv) = self.edge_parameters.borrow().get(&key) {
            return Ok(uv);
        }
        let uv = self.project_uncached(position, seed, goal, dec)?;
        self.edge_parameters.borrow_mut().insert(key, uv);
        Ok(uv)
    }

    fn project_uncached(
        &self,
        position: [f64; 3],
        seed: Option<Uv>,
        goal: f64,
        dec: &mut Decider,
    ) -> Result<Uv, Refusal> {
        let neighbours = seed.and_then(|seed| self.spans.iter().find(|span| span.contains(seed)));
        if let Some(span) = neighbours {
            let result = span.rational.project(
                position,
                seed.map(|seed| span.local_seed(seed)),
                goal,
                dec,
            )?;
            if result.residual.hi() <= goal {
                return Ok(span.global(result.uv));
            }
        }
        let mut ordered: Vec<_> = self
            .spans
            .iter()
            .filter(|span| neighbours.is_none_or(|near| !std::ptr::eq(*span, near)))
            .map(|span| (span.distance_score(position), span))
            .filter(|(gap, _)| *gap <= goal)
            .collect();
        ordered.sort_by(|a, b| a.0.total_cmp(&b.0));
        for (_, span) in ordered {
            let result = span.rational.project(
                position,
                seed.map(|seed| span.local_seed(seed)),
                goal,
                dec,
            )?;
            if result.residual.hi() <= goal {
                return Ok(span.global(result.uv));
            }
        }
        Err(failure(
            "projected boundary point exceeds the spline chart residual",
        ))
    }

    pub(crate) fn evaluate(&self, uv: Uv, dec: &mut Decider) -> Result<[Interval; 3], Refusal> {
        let span = self
            .spans
            .iter()
            .find(|span| span.contains(uv))
            .ok_or_else(|| failure("spline sample lies outside its chart"))?;
        span.rational
            .evaluate_box(span.local_box(uv.map(|value| [value; 2])), dec)
    }

    pub(crate) fn triangulation_stretch(&self) -> f64 {
        let [du, dv] = self.spans.iter().fold([0.0_f64; 2], |speed, span| {
            let local = span.rational.first_derivative_bounds();
            std::array::from_fn(|i| speed[i].max(local[i].hi() / span.width[i].lo()))
        });
        let stretch = dv / du;
        if stretch.is_finite() && stretch > 0.0 {
            stretch
        } else {
            1.0
        }
    }

    pub(crate) fn chord(&self, uv: &[Uv; 3], dec: &mut Decider) -> Result<Interval, Refusal> {
        let b = parameter_box(uv);
        let width = b.map(|[lo, hi]| point(hi) - point(lo));
        let edges: [[Interval; 2]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|axis| point(uv[(i + 1) % 3][axis]) - point(uv[i][axis]))
        });
        let mut cover = self.spans.iter().filter(|span| overlaps(span.bounds, b));
        let first = cover
            .next()
            .ok_or_else(|| failure("triangle has no spline span"))?;
        let Some(second) = cover.next() else {
            let local = [width[0] / first.width[0], width[1] / first.width[1]];
            let edges = edges.map(|e| [e[0] / first.width[0], e[1] / first.width[1]]);
            return first.chord_bound(local, Some(edges), dec);
        };
        let mut hessian = first.affine_hessian;
        for span in std::iter::once(second).chain(cover) {
            hessian = std::array::from_fn(|i| hessian[i].hull(&span.affine_hessian[i]));
        }
        dec.charge_work(6 + 2 + 3 * 4)?;
        let mut bound = triangle_chord_bound_interval(edges, width, &hessian);
        for seam in &self.seams {
            let [lo, hi] = b[seam.axis];
            let [a, z] = b[1 - seam.axis];
            if lo < seam.at && seam.at < hi && a <= seam.transverse[1] && z >= seam.transverse[0] {
                // Piecewise Taylor expansion adds the kink's hinge error.
                // The interpolation error of max(0, x - knot) is at most
                // one quarter of the coordinate range.
                dec.charge_work(2)?;
                bound = bound + point(0.25) * width[seam.axis] * seam.jump;
            }
        }
        Ok(bound)
    }

    pub(crate) fn edge_within(
        &self,
        [a, middle, b]: [[f64; 3]; 3],
        [first_hint, middle_hint, last_hint]: [Option<Uv>; 3],
        chord: Interval,
        cosine: Option<Interval>,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let goal = (chord * point(0.125)).lo();
        let a = self.project_edge(a, first_hint, None, goal, dec)?;
        let b = self.project_edge(b, last_hint, Some(a), goal, dec)?;
        let seed = std::array::from_fn(|i| a[i].midpoint(b[i]));
        let m = self.project_edge(middle, middle_hint, Some(seed), goal, dec)?;
        let uv = [a, m, b];
        if self.chord(&uv, dec)?.hi() > chord.lo() {
            return Ok(false);
        }
        let bounds = parameter_box(&uv);
        let span = self.spans.iter().find(|span| span.contains(m)).unwrap();
        if bounds
            .iter()
            .zip(span.bounds)
            .all(|([a, b], [lo, hi])| lo <= *a && *b <= hi)
        {
            let local = span.local_box(bounds);
            let cell = local.map(dyadic_cell);
            let normal = span.normal_cell(cell, dec)?;
            if normals_agree(normal.map(|n| point(n.mid())), normal, cosine, dec)? {
                return Ok(true);
            }
        }
        let at = span
            .rational
            .normal_box(span.local_box(m.map(|x| [x; 2])), dec)?;
        let direction = at.map(|n| point(n.mid()));
        self.check_normal_box(bounds, direction, cosine, dec)
    }

    fn project_edge(
        &self,
        position: [f64; 3],
        hint: Option<Uv>,
        seed: Option<Uv>,
        goal: f64,
        dec: &mut Decider,
    ) -> Result<Uv, Refusal> {
        let key = position.map(f64::to_bits);
        if let Some(&uv) = self.edge_parameters.borrow().get(&key) {
            return Ok(uv);
        }
        if let Some(uv) = hint
            && self.spans.iter().any(|span| span.contains(uv))
        {
            let expected = self.evaluate(uv, dec)?;
            let residual = position
                .into_iter()
                .zip(expected)
                .fold(point(0.0), |sum, (x, y)| {
                    let delta = point(x) - y;
                    sum + point(delta.lo().abs().max(delta.hi().abs()))
                });
            if residual.hi() <= goal {
                self.edge_parameters.borrow_mut().insert(key, uv);
                return Ok(uv);
            }
        }
        self.project(position, hint.or(seed), goal, dec)
    }

    pub(crate) fn seed_grid(
        &self,
        limit: &Rational,
        max_vertices: usize,
        dec: &mut Decider,
    ) -> Result<Vec<Uv>, Refusal> {
        let limit = Interval::from_rational(limit).lo();
        let mut total = 0usize;
        let mut grids = Vec::with_capacity(self.spans.len());
        for span in &self.spans {
            let hessian = span
                .local_affine_hessian
                .unwrap_or_else(|| span.rational.chord_coefficients());
            let mut counts = [1usize; 2];
            loop {
                let width = counts.map(|count| point(1.0) / point(count as f64));
                if span.chord_bound(width, None, dec)?.hi() <= limit {
                    break;
                }
                // The larger pure quadratic term gains most by halving its
                // width. The mixed term gains the same amount on either axis.
                let pure = [
                    hessian[0] * width[0].square(),
                    hessian[2] * width[1].square(),
                ];
                dec.charge_work(4)?;
                let axis = if pure[0].hi() > pure[1].hi() {
                    0
                } else if pure[0].hi() < pure[1].hi() {
                    1
                } else {
                    usize::from(counts[1] < counts[0])
                };
                counts[axis] = counts[axis].saturating_mul(2);
                if counts[0].saturating_mul(counts[1]) >= max_vertices {
                    return Err(failure(format!(
                        "spline grid needs more than {max_vertices} points"
                    )));
                }
            }
            // A dyadic grid supplies a feasible upper bound. Search for the
            // fewest vertices under the same chord formula, as in the 3.5
            // display prototype. Each candidate column count needs only a
            // monotone search of the other axis.
            let short = usize::from(counts[1] < counts[0]);
            let vertices =
                |[u, v]: [usize; 2]| u.saturating_add(1).saturating_mul(v.saturating_add(1));
            let mut best = vertices(counts);
            let mut columns = 1;
            while columns + 1 < best.div_ceil(2) {
                let mut candidate = [1; 2];
                candidate[short] = columns;
                let axis = 1 - short;
                let max_rows = (best - 1) / (columns + 1) - 1;
                let mut fits = |rows| {
                    candidate[axis] = rows;
                    let width = candidate.map(|count| point(1.0) / point(count as f64));
                    Ok::<_, Refusal>(span.chord_bound(width, None, dec)?.hi() <= limit)
                };
                if fits(max_rows)? {
                    let (mut low, mut high) = (1, max_rows);
                    while low < high {
                        let mid = low + (high - low) / 2;
                        if fits(mid)? {
                            high = mid;
                        } else {
                            low = mid + 1;
                        }
                    }
                    candidate[axis] = low;
                    counts = candidate;
                    best = vertices(counts);
                }
                columns += 1;
            }
            total = total.saturating_add((counts[0] + 1).saturating_mul(counts[1] + 1));
            if total > max_vertices {
                return Err(failure(format!(
                    "spline grid needs at least {total} points, cap {max_vertices}"
                )));
            }
            grids.push((span, counts));
        }
        let mut coordinates = Vec::with_capacity(total);
        for (span, [u, v]) in grids {
            dec.charge_work(2 * ((u + 1) * (v + 1)) as u64)?;
            for i in 0..=u {
                for j in 0..=v {
                    coordinates.push(span.global([i as f64 / u as f64, j as f64 / v as f64]));
                }
            }
        }
        Ok(coordinates)
    }

    pub(crate) fn normals_within(
        &self,
        uv: &[Uv; 3],
        positions: [[f64; 3]; 3],
        cosine: Option<Interval>,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let bounds = parameter_box(uv);
        let p = positions.map(|p| p.map(point));
        let a: [Interval; 3] = std::array::from_fn(|i| p[1][i] - p[0][i]);
        let b: [Interval; 3] = std::array::from_fn(|i| p[2][i] - p[0][i]);
        dec.charge_work(2 * 3)?;
        let facet = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        self.check_normal_box(bounds, facet, cosine, dec)
    }

    fn check_normal_box(
        &self,
        bounds: [[f64; 2]; 2],
        direction: [Interval; 3],
        cosine: Option<Interval>,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        for span in self
            .spans
            .iter()
            .filter(|span| overlaps(span.bounds, bounds))
        {
            if normals_agree(direction, span.normal, cosine, dec)? {
                continue;
            }
            let intersection = std::array::from_fn(|i| {
                [
                    bounds[i][0].max(span.bounds[i][0]),
                    bounds[i][1].min(span.bounds[i][1]),
                ]
            });
            let local = span.local_box(intersection);
            // Whole-cell certificates can serve adjacent facets with different
            // trim triangles. Fall back to the triangle's tighter rectangle.
            let cell = local.map(dyadic_cell);
            let normal = span.normal_cell(cell, dec)?;
            if !normals_agree(direction, normal, cosine, dec)?
                && (cell == local
                    || !normals_agree(
                        direction,
                        span.rational.normal_box(local, dec)?,
                        cosine,
                        dec,
                    )?)
            {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

fn normals_agree(
    facet: [Interval; 3],
    normal: [Interval; 3],
    cosine: Option<Interval>,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    dec.charge_work(3)?;
    let dot = facet[0] * normal[0] + facet[1] * normal[1] + facet[2] * normal[2];
    if dot.lo() <= 0.0 {
        return Ok(false);
    }
    if let Some(cosine) = cosine.filter(|cosine| cosine.hi() > 0.0) {
        // With a positive dot product, the angle is bounded by the ratio of
        // the cross product to the dot product. Lagrange's identity gives
        // the squared inequality without independently bounding two lengths.
        dec.charge_work(2 * 3 + 7)?;
        let cross = [
            facet[1] * normal[2] - facet[2] * normal[1],
            facet[2] * normal[0] - facet[0] * normal[2],
            facet[0] * normal[1] - facet[1] * normal[0],
        ];
        let cross_sq = cross.into_iter().fold(point(0.0), |s, x| s + x.square());
        let cosine_sq = cosine.square();
        let required = (point(1.0) - cosine_sq) * dot.square();
        Ok((cross_sq * cosine_sq).hi() <= required.lo())
    } else {
        Ok(true)
    }
}

pub(super) fn parameter_box(uv: &[Uv; 3]) -> [[f64; 2]; 2] {
    std::array::from_fn(|axis| {
        [
            uv.iter().map(|p| p[axis]).fold(f64::INFINITY, f64::min),
            uv.iter().map(|p| p[axis]).fold(f64::NEG_INFINITY, f64::max),
        ]
    })
}

fn overlaps(a: [[f64; 2]; 2], b: [[f64; 2]; 2]) -> bool {
    a.into_iter().zip(b).all(|([alo, ahi], [blo, bhi])| {
        if blo < bhi {
            alo < bhi && blo < ahi
        } else {
            alo <= blo && blo <= ahi
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::Real;
    use catallax_shapes::Point3;

    #[test]
    fn straight_trim_splits_when_its_incident_surface_normals_turn() {
        // The edge v = 0 is straight, but S(u, v) = (u, v, uv) has
        // normals turning by pi/4 from one end of the edge to the other.
        let zero = Rational::zero();
        let one = Rational::one();
        let knots = vec![zero.clone(), zero, one.clone(), one];
        let controls = [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 1]]
            .map(|[x, y, z]| Point3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z)));
        let mut dec = Decider::default();
        let surface = BSplineSurface::new(
            [1, 1],
            [knots.clone(), knots],
            [2, 2],
            controls.to_vec(),
            vec![Real::one(); 4],
            &mut dec,
        )
        .unwrap();
        let chart = SplineChart::new(&surface, &mut dec).unwrap();
        let chord = point(0.00125);
        let half_angle = Interval::from_rational(
            &crate::controls::cosine_limit(&Rational::from_ratio(13, 100).unwrap(), &dec).unwrap(),
        );
        let edge = |a: f64, b: f64| [a, (a + b) * 0.5, b].map(|u| [u, 0.0, 0.0]);
        assert!(
            chart
                .edge_within(edge(0.0, 1.0), [None; 3], chord, None, &mut dec)
                .unwrap()
        );
        assert!(
            !chart
                .edge_within(edge(0.0, 1.0), [None; 3], chord, Some(half_angle), &mut dec)
                .unwrap()
        );
        assert!(
            chart
                .edge_within(
                    edge(0.0, 0.0625),
                    [None; 3],
                    chord,
                    Some(half_angle),
                    &mut dec
                )
                .unwrap()
        );
    }

    #[test]
    fn a_triangle_crossing_a_knot_includes_the_derivative_jump() {
        // Each span is planar. The height peaks at 1/2 along its common
        // boundary, above a triangle whose three vertices are at height zero.
        let zero = Rational::zero();
        let one = Rational::one();
        let half = one.mul_pow2(-1);
        let mut controls = Vec::new();
        for (u, z) in [(&zero, &zero), (&half, &half), (&one, &zero)] {
            for v in [&zero, &one] {
                controls.push(Point3::new(
                    u.clone().into(),
                    v.clone().into(),
                    z.clone().into(),
                ));
            }
        }
        let mut dec = Decider::default();
        let surface = BSplineSurface::new(
            [1, 1],
            [
                vec![zero.clone(), zero.clone(), half, one.clone(), one.clone()],
                vec![zero.clone(), zero, one.clone(), one],
            ],
            [3, 2],
            controls,
            vec![Real::one(); 6],
            &mut dec,
        )
        .unwrap();
        let chart = SplineChart::new(&surface, &mut dec).unwrap();
        let bound = chart
            .chord(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], &mut dec)
            .unwrap();
        assert!(bound.hi() >= 0.5);
        assert!(bound.hi() < 0.5 + f64::EPSILON.sqrt());
    }
}
