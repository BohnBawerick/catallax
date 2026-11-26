//! A piecewise rational face chart for the constrained triangulation.

use std::{cell::RefCell, collections::BTreeMap};

use catallax_numbers::{Decider, Interval, Rational, Refusal};
use catallax_shapes::{
    BSplineCurve, BSplineSurface, RollingPipe,
    tessellation::{
        RationalSpan, chord_bound_interval, rational_spans, triangle_chord_bound_interval,
    },
};

use super::{Uv, align, failure, quantize_turn};

mod pipe_normal;
use pipe_normal::PipeNormal;

pub(crate) enum SplineNormals<'a> {
    Carrier,
    // Preserve the proof order and seeds of the existing finite spline charts.
    Patch,
    Saddle,
    Pipe(&'a RollingPipe),
}

pub(crate) struct CollarBridge {
    pub(crate) rail: [Uv; 3],
    pub(crate) apices: [Uv; 2],
    pub(crate) positions: [[f64; 3]; 3],
    pub(crate) residual: f64,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum NormalProof {
    Proven,
    Violation,
    Unresolved,
}

pub(crate) struct SplineChart {
    spans: Vec<Span>,
    rail_spans: Vec<RailSpan>,
    seams: Vec<Seam>,
    periodic: [bool; 2],
    source_domain: [[f64; 2]; 2],
    normal_station: Option<f64>,
    subdivide_normals: bool,
    patch_sampling: bool,
    stretch: f64,
    // A chart lives for one mesh request. The edge plan uses a projection
    // goal at least as strict as the face that consumes it.
    edge_projections: RefCell<BTreeMap<[u64; 3], EdgeProjection>>,
}

pub(crate) struct RailSpan {
    domain: [Rational; 2],
    knots: [f64; 3],
    curves: [RationalSpan; 2],
}

pub(crate) struct RailSample {
    pub(crate) xyz: [Interval; 3],
    uv: Uv,
    chart_offset: Interval,
}

impl RailSpan {
    pub(crate) fn curve(&self, row: usize) -> &RationalSpan {
        &self.curves[row]
    }

    pub(crate) fn sample(
        &self,
        row: usize,
        parameter: Interval,
        dec: &mut Decider,
    ) -> Result<RailSample, Refusal> {
        let curve = self.curve(row);
        let xyz = curve.evaluate_box([parameter, point(0.0)], dec)?;
        let [lo, mid, hi] = self.knots.map(point);
        dec.charge_work(4)?;
        let u = quantize_turn((lo + parameter * (hi - lo)).mid()).clamp(lo.lo(), hi.hi());
        let u = point(u);
        // The chart rounds each knot independently, including the bisection.
        // Bound the displacement between the rail parameter and the chart UV.
        let on_curve = if u.hi() <= mid.lo() {
            ((u - lo) / (mid - lo)) * point(0.5)
        } else {
            point(0.5) + ((u - mid) / (hi - mid)) * point(0.5)
        };
        let delta = on_curve - parameter;
        let chart_offset =
            point(delta.lo().abs().max(delta.hi().abs())) * curve.first_derivative_bounds()[0];
        Ok(RailSample {
            xyz,
            uv: [u.mid(), row as f64],
            chart_offset,
        })
    }
}

struct EdgeProjection {
    uv: Uv,
    residual: Option<Interval>,
}

type NormalRow = [[Interval; 3]; 5];

struct Span {
    rational: RationalSpan,
    bounds: [[f64; 2]; 2],
    width: [Interval; 2],
    periodic: [bool; 2],
    affine_hessian: [Interval; 3],
    local_affine_hessian: Option<[Interval; 3]>,
    normal: [Interval; 3],
    normal_station: Option<f64>,
    pipe_normal: Option<PipeNormal>,
    normal_cells: RefCell<BTreeMap<[u64; 4], [Interval; 3]>>,
    normal_rows: RefCell<BTreeMap<[u64; 3], NormalRow>>,
    shared_normals: RefCell<BTreeMap<[u64; 3], [Interval; 3]>>,
    angular_cones: RefCell<BTreeMap<[u64; 3], [Interval; 3]>>,
}

struct Seam {
    axis: usize,
    at: f64,
    transverse: [f64; 2],
    jump: Interval,
}

struct NormalFootprint {
    rectangle: [Interval; 2],
    // Each inequality is nonnegative on the triangle, in local span UV.
    sides: Option<[[Interval; 3]; 3]>,
    clip_first: bool,
}

#[derive(Clone, Copy)]
struct NormalTriangle<'a> {
    uv: &'a [Uv; 3],
    clip_first: bool,
}

fn point(value: f64) -> Interval {
    Interval::point(value).expect("finite chart coordinate")
}

impl Span {
    fn new(
        rational: RationalSpan,
        domain: [[&Rational; 2]; 2],
        periodic: [bool; 2],
        normal_station: Option<f64>,
        pipe_normal: Option<PipeNormal>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let bounds = std::array::from_fn(|i| {
            rational.domain()[i].each_ref().map(|end| {
                if periodic[i] {
                    let [lo, hi] = domain[i];
                    quantize_turn(((end - lo) / (hi - lo)).to_f64_nearest())
                } else {
                    end.to_f64_nearest()
                }
            })
        });
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
        let full = Interval::new(0.0, 1.0).unwrap();
        let mut span = Self {
            rational,
            bounds,
            width,
            periodic,
            affine_hessian,
            local_affine_hessian,
            normal: [full; 3],
            normal_station,
            pipe_normal,
            normal_cells: RefCell::default(),
            normal_rows: RefCell::default(),
            shared_normals: RefCell::default(),
            angular_cones: RefCell::default(),
        };
        span.normal = span.normal_box([full; 2], dec)?;
        Ok(span)
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
                let at = lo + uv[i] * self.width[i].mid();
                let at = if self.periodic[i] {
                    quantize_turn(at)
                } else {
                    at
                };
                at.clamp(lo, hi)
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

    fn normal_box(
        &self,
        mut local: [Interval; 2],
        dec: &mut Decider,
    ) -> Result<[Interval; 3], Refusal> {
        if let Some(pipe) = &self.pipe_normal {
            let source = self.rational.domain();
            let [lo, hi] = source[0].each_ref().map(Rational::to_f64_nearest);
            let quarter = lo.floor();
            let u = point(lo - quarter) + local[0] * (point(hi) - point(lo));
            let [lo, hi] = source[1].each_ref().map(Rational::to_f64_nearest);
            let v = point(lo) + local[1] * (point(hi) - point(lo));
            let u = Interval::new(u.lo().max(0.0), u.hi().min(1.0)).unwrap();
            let v = Interval::new(v.lo().max(0.0), v.hi().min(1.0)).unwrap();
            return pipe.bounds(quarter as usize, [u, v], dec);
        }
        if let Some(station) = self.normal_station {
            local[1] = point(station);
        }
        self.rational.normal_box(local, dec)
    }

    fn facet_vertex_normal(
        &self,
        triangle: &[Uv; 3],
        i: usize,
        shift: Uv,
        dec: &mut Decider,
    ) -> Result<Option<[Interval; 3]>, Refusal> {
        let vertex = triangle[i];
        let mut local = [point(0.0); 2];
        for axis in 0..2 {
            // Periodic span bounds are quantized turns, so the integer
            // lift is exact. A non-periodic span has zero shift.
            let [lo, hi] = self.bounds[axis].map(|bound| bound + shift[axis]);
            let u = vertex[axis];
            if u < lo || u > hi {
                return Ok(None);
            }
            let [lower, upper] = [lo, hi].map(|bound| (u - bound) == 0.0);
            if lower || upper {
                // A knot can have a different normal on each side. The
                // direction toward the centroid must enter this span.
                let a = point(vertex[axis]);
                let inward = (point(triangle[(i + 1) % 3][axis]) - a)
                    + (point(triangle[(i + 2) % 3][axis]) - a);
                if (lower && inward.lo() <= 0.0) || (upper && inward.hi() >= 0.0) {
                    return Ok(None);
                }
            }
            let fraction = (point(u) - point(lo)) / self.width[axis];
            local[axis] = Interval::new(fraction.lo().max(0.0), fraction.hi().min(1.0)).unwrap();
        }
        self.normal_cell(local, dec).map(Some)
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
        let normal = self.normal_box(local, dec)?;
        self.normal_cells.borrow_mut().insert(key, normal);
        Ok(normal)
    }

    fn normal_row(
        &self,
        u: Interval,
        v: [f64; 2],
        dec: &mut Decider,
    ) -> Result<Option<NormalRow>, Refusal> {
        if self.pipe_normal.is_some() {
            return Ok(None);
        }
        let key = |lo: f64| [u.lo().to_bits(), u.hi().to_bits(), lo.to_bits()];
        let [lo, hi] = v;
        // The midpoint uniquely identifies a dyadic cell within the unit interval.
        let middle = lo.midpoint(hi);
        if let Some(&row) = self.normal_rows.borrow().get(&key(middle)) {
            return Ok(Some(row));
        }
        if lo == 0.0 && hi.to_bits() == 1.0_f64.to_bits() {
            let row = self.rational.quadratic_normal_row(u, dec)?;
            if let Some(row) = row {
                self.normal_rows.borrow_mut().insert(key(middle), row);
            }
            return Ok(row);
        }
        let width = (hi - lo) * 2.0;
        let start = (lo / width).floor() * width;
        let end = start + width;
        let center = start.midpoint(end);
        let Some(mut work) = self.normal_row(u, [start, end], dec)? else {
            return Ok(None);
        };
        let mut children = [work; 2];
        dec.charge_work(3 * 10)?;
        for level in 1..5 {
            let last = 4 - level;
            for i in 0..=last {
                work[i] = std::array::from_fn(|j| (work[i][j] + work[i + 1][j]) * point(0.5));
            }
            children[0][level] = work[0];
            children[1][last] = work[last];
        }
        let mut rows = self.normal_rows.borrow_mut();
        rows.insert(key(start.midpoint(center)), children[0]);
        rows.insert(key(center.midpoint(end)), children[1]);
        Ok(Some(children[usize::from(lo >= center)]))
    }

    fn shared_normal(
        &self,
        [u, v]: [Interval; 2],
        dec: &mut Decider,
    ) -> Result<[Interval; 3], Refusal> {
        let key = [u.lo().to_bits(), u.hi().to_bits(), v.mid().to_bits()];
        if let Some(&normal) = self.shared_normals.borrow().get(&key) {
            return Ok(normal);
        }
        let Some(row) = self.normal_row(u, [v.lo(), v.hi()], dec)? else {
            return self.normal_cell([u, v], dec);
        };
        // Each coefficient may be rescaled by a positive factor before
        // bounding its direction. The new convex hull covers the same rays.
        dec.charge_work(5 * 3)?;
        let scaled = row.map(|normal| {
            let length = normal.iter().map(|c| c.mid().abs()).fold(0.0_f64, f64::max);
            if !length.is_finite() || length == 0.0 {
                return normal;
            }
            let scale = length.recip();
            if !scale.is_finite() || scale == 0.0 {
                return normal;
            }
            normal.map(|x| x * point(scale))
        });
        let normal = std::array::from_fn(|i| {
            scaled
                .iter()
                .skip(1)
                .fold(scaled[0][i], |hull, control| hull.hull(&control[i]))
        });
        self.shared_normals.borrow_mut().insert(key, normal);
        Ok(normal)
    }

    fn angular_cone(
        &self,
        [u, v]: [Interval; 2],
        dec: &mut Decider,
    ) -> Result<[Interval; 3], Refusal> {
        let key = [u.lo().to_bits(), u.hi().to_bits(), v.mid().to_bits()];
        if let Some(&cone) = self.angular_cones.borrow().get(&key) {
            return Ok(cone);
        }
        let Some(row) = self.normal_row(u, [v.lo(), v.hi()], dec)? else {
            return self.normal_cell([u, v], dec);
        };
        let cone = if let Some((axis, sign)) = normal_cone_axis(&row) {
            // Divide each control by its positive projection onto the same
            // axis. The resulting box lies in a plane at distance one, and
            // the original Bernstein normal is in the cone over that box.
            dec.charge_work(5 * 3)?;
            let project = |control: [Interval; 3]| {
                let scale = point(sign) * control[axis];
                std::array::from_fn(|i| {
                    if i == axis {
                        point(sign)
                    } else {
                        control[i] / scale
                    }
                })
            };
            row.into_iter()
                .map(project)
                .reduce(|a, b| std::array::from_fn(|i| a[i].hull(&b[i])))
                .unwrap()
        } else {
            self.shared_normal([u, v], dec)?
        };
        self.angular_cones.borrow_mut().insert(key, cone);
        Ok(cone)
    }
}

fn normal_cone_axis(row: &NormalRow) -> Option<(usize, f64)> {
    let mut best = None;
    let mut margin = 0.0;
    for axis in 0..3 {
        let hull = row
            .iter()
            .skip(1)
            .fold(row[0][axis], |a, b| a.hull(&b[axis]));
        for (component, sign) in [(hull.lo(), 1.0), (-hull.hi(), -1.0)] {
            if component > margin {
                margin = component;
                best = Some((axis, sign));
            }
        }
    }
    best
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
    pub(super) fn new(
        surface: &BSplineSurface,
        periodic: [bool; 2],
        join_scale: [Option<Interval>; 2],
        subdivide: bool,
        normal_station: Option<f64>,
        normals: SplineNormals<'_>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
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
        let source_domain = surface.domain();
        for (i, wrapped) in periodic.into_iter().enumerate() {
            if wrapped && join_scale[i].is_none() && !closed_boundaries(surface, i) {
                return Err(failure(
                    "periodic spline boundaries need identical homogeneous controls",
                ));
            }
        }
        let saddle = matches!(normals, SplineNormals::Saddle);
        let pipe_normal = match normals {
            SplineNormals::Pipe(pipe) => Some(PipeNormal::new(pipe, dec)?),
            _ => None,
        };
        let mut rail_spans = Vec::new();
        let pieces = if saddle {
            let mut pieces = Vec::new();
            for span in rational_spans(surface, dec)? {
                let halves = span.bisect(0, dec)?;
                if span.domain()[1].each_ref() == source_domain[1] {
                    let [lo, hi] = source_domain[0];
                    let knots = [
                        &span.domain()[0][0],
                        &halves[0].domain()[0][1],
                        &span.domain()[0][1],
                    ]
                    .map(|end| quantize_turn(((end - lo) / (hi - lo)).to_f64_nearest()));
                    dec.charge_traversal(1)?;
                    rail_spans.push(RailSpan {
                        domain: span.domain()[0].clone(),
                        knots,
                        curves: [
                            span.boundary_curve(false, dec)?,
                            span.boundary_curve(true, dec)?,
                        ],
                    });
                }
                for half in halves {
                    pieces.extend(half.bisect(1, dec)?);
                }
            }
            pieces
        } else if subdivide {
            split_spans(surface, dec)?
        } else {
            rational_spans(surface, dec)?
        };
        let spans: Vec<_> = pieces
            .into_iter()
            .map(|span| {
                Span::new(
                    span,
                    source_domain,
                    periodic,
                    normal_station,
                    pipe_normal,
                    dec,
                )
            })
            .collect::<Result<_, _>>()?;
        let [du, dv] = spans.iter().fold([0.0_f64; 2], |speed, span| {
            let local = span.rational.first_derivative_bounds();
            std::array::from_fn(|i| speed[i].max(local[i].hi() / span.width[i].lo()))
        });
        let stretch = dv / du;
        let mut seams = Vec::new();
        for span in &spans {
            let domain = span.rational.domain();
            for axis in 0..2 {
                let other = 1 - axis;
                let end = &domain[axis][1];
                let wraps = periodic[axis] && end == source_domain[axis][1];
                let next_start = if wraps { source_domain[axis][0] } else { end };
                let Some(next) = spans.iter().find(|next| {
                    &next.rational.domain()[axis][0] == next_start
                        && next.rational.domain()[other] == domain[other]
                }) else {
                    continue;
                };
                let widths = [span.width[axis], next.width[axis]];
                let jump = if wraps {
                    span.rational.periodic_first_jump(
                        &next.rational,
                        axis,
                        widths,
                        join_scale[axis].unwrap_or_else(|| point(1.0)),
                        dec,
                    )?
                } else {
                    span.rational
                        .first_jump(&next.rational, axis, widths, dec)?
                };
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
            rail_spans,
            seams,
            periodic,
            source_domain: source_domain.map(|bounds| bounds.map(Rational::to_f64_nearest)),
            normal_station,
            subdivide_normals: saddle,
            patch_sampling: matches!(normals, SplineNormals::Patch),
            stretch: if stretch.is_finite() && stretch > 0.0 {
                stretch
            } else {
                1.0
            },
            edge_projections: RefCell::default(),
        })
    }

    pub(crate) fn checked_rail_spans(
        &self,
        curve: &BSplineCurve,
        dec: &mut Decider,
    ) -> Result<Option<&[RailSpan]>, Refusal> {
        if self.rail_spans.is_empty() {
            return Ok(None);
        }
        let mut spans = self.rail_spans.iter();
        for (_, knots) in curve.spans() {
            dec.charge_traversal(1)?;
            if spans
                .next()
                .is_none_or(|span| span.domain.each_ref() != knots)
            {
                return Ok(None);
            }
        }
        Ok(spans.next().is_none().then_some(&self.rail_spans))
    }

    pub(crate) fn admit_rail_sample(
        &self,
        sample: &RailSample,
        position: [f64; 3],
        goal: f64,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        dec.charge_work(3)?;
        let residual = vertex_residual(position, sample.xyz) + sample.chart_offset;
        if residual.hi() <= goal {
            dec.charge_traversal(1)?;
            self.edge_projections.borrow_mut().insert(
                position.map(f64::to_bits),
                EdgeProjection {
                    uv: [sample.uv[0].rem_euclid(1.0), sample.uv[1]],
                    residual: Some(residual),
                },
            );
        }
        Ok(())
    }

    fn lifts(&self, bounds: [[f64; 2]; 2]) -> impl Iterator<Item = Uv> + use<> {
        let [u, v] = std::array::from_fn(|i| {
            let [lo, hi] = bounds[i];
            if self.periodic[i] {
                (lo.ceil() as i32 - 1)..=(hi.floor() as i32)
            } else {
                0..=0
            }
        });
        u.flat_map(move |i| v.clone().map(move |j| [f64::from(i), f64::from(j)]))
    }

    fn cover(&self, bounds: [[f64; 2]; 2]) -> impl Iterator<Item = (&Span, [[f64; 2]; 2], Uv)> {
        self.lifts(bounds).flat_map(move |shift| {
            let local = subtract_shift(bounds, shift);
            self.spans
                .iter()
                .filter(move |span| overlaps(span.bounds, local))
                .map(move |span| (span, local, shift))
        })
    }

    fn normalize_hint(&self, hint: Uv) -> Uv {
        std::array::from_fn(|i| {
            if self.periodic[i] {
                let [lo, hi] = self.source_domain[i];
                quantize_turn((hint[i] - lo) / (hi - lo))
            } else {
                hint[i]
            }
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
        if let Some(projection) = self.edge_projections.borrow().get(&key) {
            return Ok(projection.uv);
        }
        let uv = self.project_uncached(position, seed, None, goal, dec)?;
        self.edge_projections
            .borrow_mut()
            .insert(key, EdgeProjection { uv, residual: None });
        Ok(uv)
    }

    fn project_uncached(
        &self,
        position: [f64; 3],
        seed: Option<Uv>,
        fixed: Option<(usize, f64)>,
        goal: f64,
        dec: &mut Decider,
    ) -> Result<Uv, Refusal> {
        let seed = seed.map(|seed| {
            std::array::from_fn(|i| {
                if self.periodic[i] {
                    seed[i].rem_euclid(1.0)
                } else {
                    seed[i]
                }
            })
        });
        let try_span = |span: &Span, dec: &mut Decider| {
            let result = span.rational.project(
                position,
                seed.map(|seed| span.local_seed(seed)),
                goal,
                dec,
            )?;
            if result.residual.hi() > goal {
                return Ok::<_, Refusal>(None);
            }
            let mut uv = span.global(result.uv);
            if let Some((axis, value)) = fixed {
                uv[axis] = value;
                // A span endpoint can satisfy the projection's Euclidean
                // goal and still miss the rail certificate. Try the other
                // spans before settling the shared point's parameter.
                if self.residual(uv, position, dec)?.hi() > goal {
                    return Ok(None);
                }
            }
            Ok(Some(uv))
        };
        let neighbours = seed.and_then(|seed| self.spans.iter().find(|span| span.contains(seed)));
        if let Some(span) = neighbours
            && let Some(uv) = try_span(span, dec)?
        {
            return Ok(uv);
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
            if let Some(uv) = try_span(span, dec)? {
                return Ok(uv);
            }
        }
        Err(failure(if fixed.is_some() {
            "projected rail point exceeds the spline chart residual"
        } else {
            "projected boundary point exceeds the spline chart residual"
        }))
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
        self.stretch
    }

    pub(crate) fn has_normal_cover(&self) -> bool {
        self.subdivide_normals || self.spans[0].pipe_normal.is_some()
    }

    pub(crate) fn longitudinal_span(&self, u: f64) -> u64 {
        let u = u.rem_euclid(1.0);
        self.spans
            .iter()
            .map(|span| span.bounds[0])
            .find(|[a, b]| *a <= u && u < *b)
            .map_or(u.to_bits(), |[a, _]| a.to_bits())
    }

    pub(crate) fn exact_metric_stretch(&self) -> f64 {
        let (_, exponent) = libm::frexp(self.triangulation_stretch());
        libm::ldexp(1.0, exponent - 1)
    }

    pub(crate) fn rail_anchor(a: Uv, b: Uv, stretch: f64) -> Option<Uv> {
        // The altitude is strictly outside the diametral disk. Keep opposite
        // collars apart, even when the rails have very different spacings.
        let offset = (b[0] - a[0]).abs() * (0.625 / stretch);
        (offset > 0.0 && offset <= 0.25).then(|| {
            [
                quantize_turn(a[0].midpoint(b[0])),
                if a[1] == 0.0 { offset } else { 1.0 - offset },
            ]
        })
    }

    pub(crate) fn collar_bridge_fits(
        &self,
        bridge: CollarBridge,
        options: &crate::MeshOptions,
        cosine: Interval,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let CollarBridge {
            rail: [left, base, right],
            apices: [a, b],
            mut positions,
            residual,
        } = bridge;
        // Both disk centers must be on the rail side of the bridge by at
        // least their radii, using the rounded metric positions seen by CDT.
        dec.charge_work(31)?;
        let stretch = self.exact_metric_stretch();
        let width = point(b[0]) - point(a[0]);
        let rise = point(b[1] * stretch) - point(a[1] * stretch);
        let side = point(if base[1] == 0.0 { -1.0 } else { 1.0 });
        let length_sq = width.square() + rise.square();
        let clear = [left, right].into_iter().all(|end| {
            let center = (point(end[0]) + point(base[0])) * point(0.5);
            let radius = (point(end[0]) - point(base[0])) * point(0.5);
            let distance = side
                * (width * (point(base[1] * stretch) - point(a[1] * stretch))
                    - rise * (center - point(a[0])));
            distance.lo() > 0.0 && distance.square().lo() > (radius.square() * length_sq).hi()
        });
        if !clear {
            return Ok(false);
        }
        let mut uv = [base, a, b];
        if base[1] == 0.0 {
            uv.swap(1, 2);
            positions.swap(1, 2);
        }
        self.collar_facet_fits(&uv, positions, residual, options, cosine, dec)
    }

    pub(crate) fn collar_facet_fits(
        &self,
        uv: &[Uv; 3],
        positions: [[f64; 3]; 3],
        residual: f64,
        options: &crate::MeshOptions,
        cosine: Interval,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        if let Some(limit) = &options.max_edge_length
            && (0..3).any(|i| {
                crate::controls::length_squared(positions[i], positions[(i + 1) % 3])
                    > limit * limit
            })
        {
            return Ok(false);
        }
        let tolerance = Interval::from_rational(&options.chord_tolerance).lo();
        if (self.chord(uv, dec)? + point(residual)).hi() > tolerance
            && self.split_chord(uv, positions, residual, dec)?.hi() > tolerance
        {
            return Ok(false);
        }
        Ok(self.normal_proof(uv, positions, Some(cosine), dec)? == NormalProof::Proven)
    }

    pub(crate) fn rail_pair_fits(
        &self,
        samples: [[f64; 3]; 3],
        row: f64,
        options: &crate::MeshOptions,
        cosine: Interval,
        seam_cosine: Option<Interval>,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let goal = Interval::from_rational(&options.chord_tolerance.mul_pow2(-6)).lo();
        let mut rail = [[0.0; 2]; 3];
        for (uv, sample) in rail.iter_mut().zip(samples) {
            *uv = self.project_edge(sample, None, None, Some((1, row)), goal, dec)?;
        }
        let base = rail[1];
        rail[0] = align(rail[0], base, self.periodic);
        rail[2] = align(rail[2], base, self.periodic);
        if rail[0][0] > rail[2][0] {
            rail.swap(0, 2);
        }
        let stretch = self.exact_metric_stretch();
        let [Some(left), Some(right)] = [
            Self::rail_anchor(rail[0], base, stretch),
            Self::rail_anchor(base, rail[2], stretch),
        ] else {
            return Ok(false);
        };
        let canonical = |uv: Uv| [uv[0].rem_euclid(1.0), uv[1]];
        let left_position = self
            .evaluate(canonical(left), dec)?
            .map(|value| value.mid());
        let right_position = self
            .evaluate(canonical(right), dec)?
            .map(|value| value.mid());
        let positions = [samples[1], left_position, right_position];
        let mut residual = 0.0_f64;
        for (uv, position) in [base, left, right].into_iter().zip(positions) {
            residual = residual.max(self.residual(canonical(uv), position, dec)?.hi());
        }
        if !self.collar_bridge_fits(
            CollarBridge {
                rail,
                apices: [left, right],
                positions,
                residual,
            },
            options,
            cosine,
            dec,
        )? {
            return Ok(false);
        }
        if base[0] != 0.0 {
            return Ok(true);
        }
        // At the chart join, both sides of the wrapped bridge meet a
        // single seam sample. Qualify that seam segment and its two
        // facets before either neighboring rail leaf becomes fixed.
        let fraction = -left[0] / (right[0] - left[0]);
        let tip = [0.0, left[1] + (right[1] - left[1]) * fraction];
        let middle = [0.0, base[1].midpoint(tip[1])];
        let chord = Interval::from_rational(&options.chord_tolerance.mul_pow2(-3));
        if !self.seam_within([base, middle, tip], chord, seam_cosine, dec)? {
            return Ok(false);
        }
        let tip_position = self.evaluate(tip, dec)?.map(|value| value.mid());
        residual = residual.max(self.residual(tip, tip_position, dec)?.hi());
        for (apex, position) in [(left, left_position), (right, right_position)] {
            let mut uv = [base, tip, apex];
            let mut points = [samples[1], tip_position, position];
            if !super::positive_turn(uv[0], uv[1], uv[2]) {
                uv.swap(0, 1);
                points.swap(0, 1);
            }
            if !self.collar_facet_fits(&uv, points, residual, options, cosine, dec)? {
                return Ok(false);
            }
        }
        Ok(true)
    }

    pub(crate) fn chord(&self, uv: &[Uv; 3], dec: &mut Decider) -> Result<Interval, Refusal> {
        self.chord_for(uv, None, dec)
    }

    pub(crate) fn chord_with_limit(
        &self,
        uv: &[Uv; 3],
        limit: f64,
        dec: &mut Decider,
    ) -> Result<Interval, Refusal> {
        self.chord_for(uv, self.patch_sampling.then_some(limit), dec)
    }

    fn chord_for(
        &self,
        uv: &[Uv; 3],
        box_limit: Option<f64>,
        dec: &mut Decider,
    ) -> Result<Interval, Refusal> {
        let b = parameter_box(uv);
        let width = b.map(|[lo, hi]| point(hi) - point(lo));
        let edges: [[Interval; 2]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|axis| point(uv[(i + 1) % 3][axis]) - point(uv[i][axis]))
        });
        let mut cover = self.cover(b);
        let (first, _, _) = cover
            .next()
            .ok_or_else(|| failure("triangle has no spline span"))?;
        let Some(second) = cover.next() else {
            let local = [width[0] / first.width[0], width[1] / first.width[1]];
            if let Some(limit) = box_limit {
                let bound = first.chord_bound(local, None, dec)?;
                if bound.hi() <= limit {
                    return Ok(bound);
                }
            }
            let edges = edges.map(|e| [e[0] / first.width[0], e[1] / first.width[1]]);
            return first.chord_bound(local, Some(edges), dec);
        };
        let mut hessian = first.affine_hessian;
        for (span, _, _) in std::iter::once(second).chain(cover) {
            hessian = std::array::from_fn(|i| hessian[i].hull(&span.affine_hessian[i]));
        }
        dec.charge_work(6 + 2 + 3 * 4)?;
        let mut bound = triangle_chord_bound_interval(edges, width, &hessian);
        for shift in self.lifts(b) {
            let local = subtract_shift(b, shift);
            for seam in &self.seams {
                let [lo, hi] = local[seam.axis];
                let [a, z] = local[1 - seam.axis];
                if lo < seam.at
                    && seam.at < hi
                    && a <= seam.transverse[1]
                    && z >= seam.transverse[0]
                {
                    // Piecewise Taylor expansion adds the kink's hinge error.
                    // The interpolation error of max(0, x - knot) is at most
                    // one quarter of the coordinate range.
                    dec.charge_work(2)?;
                    bound = bound + point(0.25) * width[seam.axis] * seam.jump;
                }
            }
        }
        Ok(bound)
    }

    pub(crate) fn split_chord(
        &self,
        uv: &[Uv; 3],
        positions: [[f64; 3]; 3],
        vertex_error: f64,
        dec: &mut Decider,
    ) -> Result<Interval, Refusal> {
        let bounds = parameter_box(uv);
        let half = point(0.5);
        let widths = bounds.map(|[lo, hi]| (point(hi) - point(lo)) * half);
        let hessian = self
            .cover(bounds)
            .map(|(span, _, _)| span.affine_hessian)
            .reduce(|a, b| std::array::from_fn(|i| a[i].hull(&b[i])))
            .ok_or_else(|| failure("triangle has no spline span"))?;
        dec.charge_work(6)?;
        // The four virtual subtriangles have half the original chart widths.
        // Their common Hessian bound applies to S minus the affine facet, too.
        let mut error = chord_bound_interval(widths, &hessian);
        for shift in self.lifts(bounds) {
            let local = subtract_shift(bounds, shift);
            for seam in &self.seams {
                let [lo, hi] = local[seam.axis];
                let [a, z] = local[1 - seam.axis];
                if lo < seam.at
                    && seam.at < hi
                    && a <= seam.transverse[1]
                    && z >= seam.transverse[0]
                {
                    dec.charge_work(2)?;
                    error = error + point(0.25) * widths[seam.axis] * seam.jump;
                }
            }
        }
        let mut nodal_error = vertex_error;
        for i in 0..3 {
            let j = (i + 1) % 3;
            let midpoint: [Interval; 2] =
                std::array::from_fn(|axis| (point(uv[i][axis]) + point(uv[j][axis])) * half);
            let midpoint = std::array::from_fn(|axis| {
                let [lo, hi] = bounds[axis];
                [midpoint[axis].lo().max(lo), midpoint[axis].hi().min(hi)]
            });
            let mut found = false;
            for shift in self.lifts(midpoint) {
                let local: [Interval; 2] = std::array::from_fn(|axis| {
                    let [lo, hi] = midpoint[axis];
                    Interval::new(lo, hi).unwrap() - point(shift[axis])
                });
                for span in &self.spans {
                    let local_bounds = local.map(|range| [range.lo(), range.hi()]);
                    if !overlaps(span.bounds, local_bounds) {
                        continue;
                    }
                    let clipped = std::array::from_fn(|axis| {
                        [
                            local[axis].lo().max(span.bounds[axis][0]),
                            local[axis].hi().min(span.bounds[axis][1]),
                        ]
                    });
                    let actual = span.rational.evaluate_box(span.local_box(clipped), dec)?;
                    dec.charge_work(12)?;
                    let squared = (0..3).fold(point(0.0), |sum, axis| {
                        let affine = (point(positions[i][axis]) + point(positions[j][axis])) * half;
                        sum + (actual[axis] - affine).square()
                    });
                    let distance = squared.sqrt().expect("sum of squared residuals");
                    nodal_error = nodal_error.max(distance.hi());
                    found = true;
                }
            }
            if !found {
                return Err(failure("triangle midpoint has no spline span"));
            }
        }
        Ok(error + point(nodal_error))
    }

    // The singular station uses x = u*w and y = ±w, where w is the
    // distance in parameter space from the collapsed end.
    fn collapsed_arc(x: f64, w: f64) -> Interval {
        if w == 0.0 || x == 0.0 {
            return point(0.0);
        }
        if x.to_bits() == w.to_bits() {
            return point(1.0);
        }
        let u = point(x) / point(w);
        Interval::new(u.lo().max(0.0), u.hi().min(1.0)).unwrap()
    }

    pub(crate) fn evaluate_collapsed(
        &self,
        [x, y]: Uv,
        dec: &mut Decider,
    ) -> Result<[Interval; 3], Refusal> {
        let w = point(y.abs());
        let v = if self.normal_station == Some(1.0) {
            w
        } else {
            let v = point(1.0) - w;
            Interval::new(v.lo().max(0.0), v.hi().min(1.0)).unwrap()
        };
        self.spans[0]
            .rational
            .evaluate_box([Self::collapsed_arc(x, y.abs()), v], dec)
    }

    pub(crate) fn collapsed_facet(
        &self,
        at: [Uv; 3],
        dec: &mut Decider,
    ) -> Result<(Interval, [[f64; 2]; 2]), Refusal> {
        let mut arc = at
            .into_iter()
            .filter(|uv| uv[1] != 0.0)
            .map(|[x, y]| Self::collapsed_arc(x, y.abs()));
        let first = arc
            .next()
            .ok_or_else(|| failure("linear fillet facet has no regular station"))?;
        let arc = arc.fold(first, |a, b| a.hull(&b));
        let w = at.iter().map(|uv| uv[1].abs()).fold(0.0, f64::max);
        // At barycentric weights λ, S = A + w(B(u) - A),
        // w = Σλ_i*w_i and u = Σ(λ_i*w_i/w)*u_i. The apex has zero
        // weight in the arc interpolation, even for a triangle incident to it.
        let width = point(arc.hi()) - point(arc.lo());
        dec.charge_work(6)?;
        let error = point(w)
            * chord_bound_interval(
                [width, point(0.0)],
                &self.spans[0].rational.affine_hessian(),
            );
        Ok((error, [[arc.lo(), arc.hi()], [0.0, 1.0]]))
    }

    pub(crate) fn edge_within(
        &self,
        positions: [[f64; 3]; 3],
        [first_hint, middle_hint, last_hint]: [Option<Uv>; 3],
        fixed: Option<(usize, f64)>,
        options: &crate::MeshOptions,
        cosine: Option<Interval>,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let [a, middle, b] = positions;
        let chord = Interval::from_rational(&options.chord_tolerance.mul_pow2(-3));
        let goal = (chord * point(0.125)).lo();
        let a = self.project_edge(a, first_hint, None, fixed, goal, dec)?;
        let b = align(
            self.project_edge(b, last_hint, Some(a), fixed, goal, dec)?,
            a,
            self.periodic,
        );
        let seed = std::array::from_fn(|i| a[i].midpoint(b[i]));
        let middle = self.project_edge(middle, middle_hint, Some(seed), fixed, goal, dec)?;
        let m = align(middle, seed, self.periodic);
        let saddle_rail = self.subdivide_normals && matches!(fixed, Some((1, 0.0 | 1.0)));
        // Both adjacent surfaces interpolate the same sampled polyline.
        // Reserve half the edge allowance for each side of the collar.
        let rail_chord = if saddle_rail {
            chord * point(0.5)
        } else {
            chord
        };
        if !self.edge_parameters_within([a, m, b], middle, rail_chord, cosine, dec)? {
            return Ok(false);
        }
        if saddle_rail {
            let Some(inward) = Self::rail_anchor(a, b, self.exact_metric_stretch()) else {
                return Ok(false);
            };
            let uv = if (a[0] < b[0]) == (a[1] == 0.0) {
                [a, b, inward]
            } else {
                [b, a, inward]
            };
            if self.chord(&uv, dec)?.hi() > (point(4.0) * chord).lo() {
                return Ok(false);
            }
            let inside = self.evaluate([inward[0].rem_euclid(1.0), inward[1]], dec)?;
            let inside = inside.map(|value| value.mid());
            let positions = if (a[0] < b[0]) == (a[1] == 0.0) {
                [positions[0], positions[2], inside]
            } else {
                [positions[2], positions[0], inside]
            };
            if let Some(length) = &options.max_edge_length {
                let squared = length * length;
                if positions[..2]
                    .iter()
                    .any(|&endpoint| crate::controls::length_squared(endpoint, inside) > squared)
                {
                    return Ok(false);
                }
            }
            let facet_cosine = cosine
                .filter(|cosine| cosine.hi() > 0.0)
                .map(|cosine| point(2.0) * cosine.square() - point(1.0));
            return self.normals_within(&uv, positions, facet_cosine, dec);
        }
        Ok(true)
    }

    pub(crate) fn seam_within(
        &self,
        uv: [Uv; 3],
        chord: Interval,
        cosine: Option<Interval>,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let middle = std::array::from_fn(|i| {
            if self.periodic[i] {
                uv[1][i].rem_euclid(1.0)
            } else {
                uv[1][i]
            }
        });
        self.edge_parameters_within(uv, middle, chord, cosine, dec)
    }

    fn edge_parameters_within(
        &self,
        uv: [Uv; 3],
        middle: Uv,
        chord: Interval,
        cosine: Option<Interval>,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let m = uv[1];
        let bounds = parameter_box(&uv);
        if (0..2).any(|i| self.periodic[i] && bounds[i][1] - bounds[i][0] > 0.125) {
            // Keep the same local lift limit as the analytic chart sampler.
            return Ok(false);
        }
        if self.chord_with_limit(&uv, chord.lo(), dec)?.hi() > chord.lo() {
            return Ok(false);
        }
        let span = self
            .spans
            .iter()
            .find(|span| span.contains(middle))
            .unwrap();
        let local_bounds = subtract_shift(bounds, std::array::from_fn(|i| m[i] - middle[i]));
        if local_bounds
            .iter()
            .zip(span.bounds)
            .all(|([a, b], [lo, hi])| lo <= *a && *b <= hi)
        {
            let local = span.local_box(local_bounds);
            let cell = local.map(dyadic_cell);
            let normal = span.normal_cell(cell, dec)?;
            if normals_agree(normal.map(|n| point(n.mid())), normal, cosine, dec)? {
                return Ok(true);
            }
        }
        let at = span.normal_box(span.local_box(middle.map(|x| [x; 2])), dec)?;
        let direction = at.map(|n| point(n.mid()));
        Ok(self.check_normal_box(bounds, None, direction, cosine, dec)? == NormalProof::Proven)
    }

    fn project_edge(
        &self,
        position: [f64; 3],
        hint: Option<Uv>,
        seed: Option<Uv>,
        fixed: Option<(usize, f64)>,
        goal: f64,
        dec: &mut Decider,
    ) -> Result<Uv, Refusal> {
        let key = position.map(f64::to_bits);
        let keep_on_rail = |mut uv: Uv| {
            if let Some((axis, value)) = fixed {
                uv[axis] = value;
            }
            uv
        };
        if (self.normal_station.is_none() || hint.is_none())
            && let Some(projection) = self.edge_projections.borrow().get(&key)
        {
            return Ok(keep_on_rail(projection.uv));
        }
        let hint = hint.map(|uv| keep_on_rail(self.normalize_hint(uv)));
        if let Some(uv) = hint
            && self.spans.iter().any(|span| span.contains(uv))
        {
            let residual = self.residual(uv, position, dec)?;
            if residual.hi() <= goal {
                self.edge_projections.borrow_mut().insert(
                    key,
                    EdgeProjection {
                        uv,
                        residual: Some(residual),
                    },
                );
                return Ok(uv);
            }
        }
        if fixed.is_some() {
            // The pipe law puts each complete rail at a fixed v. Retain that
            // parameter line before the CDT classifies the periodic loops.
            let uv = self.project_uncached(position, hint.or(seed), fixed, goal, dec)?;
            self.edge_projections
                .borrow_mut()
                .insert(key, EdgeProjection { uv, residual: None });
            return Ok(uv);
        }
        self.project(position, hint.or(seed), goal, dec)
    }

    fn residual(&self, uv: Uv, position: [f64; 3], dec: &mut Decider) -> Result<Interval, Refusal> {
        Ok(vertex_residual(position, self.evaluate(uv, dec)?))
    }

    pub(crate) fn patch_sampling(&self) -> bool {
        self.patch_sampling
    }

    pub(crate) fn vertex(
        &self,
        uv: Uv,
        replacement: Option<[f64; 3]>,
        dec: &mut Decider,
    ) -> Result<([f64; 3], Interval), Refusal> {
        if self.patch_sampling
            && let Some(position) = replacement
            && let Some(projection) = self
                .edge_projections
                .borrow()
                .get(&position.map(f64::to_bits))
            && super::key(projection.uv) == super::key(uv)
            && let Some(residual) = projection.residual
        {
            return Ok((position, residual));
        }
        let enclosed = self.evaluate(uv, dec)?;
        Ok(Self::rounded_vertex(enclosed, replacement))
    }

    pub(crate) fn rounded_vertex(
        enclosed: [Interval; 3],
        replacement: Option<[f64; 3]>,
    ) -> ([f64; 3], Interval) {
        let position = replacement.unwrap_or_else(|| enclosed.map(|value| value.mid()));
        (position, vertex_residual(position, enclosed))
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
        'spans: for span in &self.spans {
            let remaining = max_vertices - total;
            if self.patch_sampling && remaining < 4 {
                break;
            }
            let vertices =
                |[u, v]: [usize; 2]| u.saturating_add(1).saturating_mul(v.saturating_add(1));
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
                if self.patch_sampling {
                    if vertices(counts) > remaining.saturating_mul(4) {
                        // Leave large seed grids to adaptive facet refinement.
                        continue 'spans;
                    }
                } else if counts[0].saturating_mul(counts[1]) >= max_vertices {
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
            let mut best = vertices(counts);
            if self.patch_sampling {
                best = best.min(remaining.saturating_add(1));
            }
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
            if self.patch_sampling && vertices(counts) > remaining {
                continue;
            }
            total = total.saturating_add(vertices(counts));
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
        Ok(self.triangle_normals(uv, positions, cosine, false, dec)? == NormalProof::Proven)
    }

    pub(crate) fn normal_proof(
        &self,
        uv: &[Uv; 3],
        positions: [[f64; 3]; 3],
        cosine: Option<Interval>,
        dec: &mut Decider,
    ) -> Result<NormalProof, Refusal> {
        self.triangle_normals(uv, positions, cosine, self.subdivide_normals, dec)
    }

    fn triangle_normals(
        &self,
        uv: &[Uv; 3],
        positions: [[f64; 3]; 3],
        cosine: Option<Interval>,
        clip_first: bool,
        dec: &mut Decider,
    ) -> Result<NormalProof, Refusal> {
        if self.has_normal_cover() && !super::positive_turn(uv[0], uv[1], uv[2]) {
            return Err(failure("spline facet has no positive parameter area"));
        }
        let triangle = NormalTriangle {
            uv,
            clip_first: clip_first && cosine.is_some(),
        };
        self.facet_normals(parameter_box(uv), Some(triangle), positions, cosine, dec)
    }

    pub(crate) fn facet_normals_within(
        &self,
        bounds: [[f64; 2]; 2],
        positions: [[f64; 3]; 3],
        cosine: Option<Interval>,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        Ok(self.facet_normals(bounds, None, positions, cosine, dec)? == NormalProof::Proven)
    }

    fn facet_normals(
        &self,
        bounds: [[f64; 2]; 2],
        triangle: Option<NormalTriangle<'_>>,
        positions: [[f64; 3]; 3],
        cosine: Option<Interval>,
        dec: &mut Decider,
    ) -> Result<NormalProof, Refusal> {
        let p = positions.map(|p| p.map(point));
        let a: [Interval; 3] = std::array::from_fn(|i| p[1][i] - p[0][i]);
        let b: [Interval; 3] = std::array::from_fn(|i| p[2][i] - p[0][i]);
        dec.charge_work(2 * 3)?;
        let facet = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        self.check_normal_box(bounds, triangle, facet, cosine, dec)
    }

    fn check_normal_box(
        &self,
        bounds: [[f64; 2]; 2],
        triangle: Option<NormalTriangle<'_>>,
        direction: [Interval; 3],
        cosine: Option<Interval>,
        dec: &mut Decider,
    ) -> Result<NormalProof, Refusal> {
        for (span, bounds, shift) in self.cover(bounds) {
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
            if self.patch_sampling {
                let coarse = local.map(|axis| {
                    Interval::new(
                        (axis.lo() * 4.0).floor() * 0.25,
                        (axis.hi() * 4.0).ceil() * 0.25,
                    )
                    .unwrap()
                });
                if coarse != [Interval::new(0.0, 1.0).unwrap(); 2]
                    && normals_agree(direction, span.normal_cell(coarse, dec)?, cosine, dec)?
                {
                    continue;
                }
            }
            if self.has_normal_cover() {
                let unit = [Interval::new(0.0, 1.0).unwrap(); 2];
                if triangle.is_some() {
                    dec.charge_work(12)?;
                }
                let sides = triangle.map(|triangle| {
                    let uv = triangle.uv;
                    let origin: [Interval; 2] =
                        std::array::from_fn(|i| point(span.bounds[i][0]) + point(shift[i]));
                    std::array::from_fn(|i| {
                        let a = uv[i].map(point);
                        let b = uv[(i + 1) % 3].map(point);
                        let u = b[0] - a[0];
                        let v = b[1] - a[1];
                        [
                            -v * span.width[0],
                            u * span.width[1],
                            u * (origin[1] - a[1]) - v * (origin[0] - a[0]),
                        ]
                    })
                });
                let target = NormalFootprint {
                    rectangle: local,
                    sides,
                    clip_first: triangle.is_some_and(|triangle| triangle.clip_first),
                };
                if cosine.is_some_and(|value| value.lo() > 0.0)
                    && let Some(triangle) = triangle
                {
                    for i in 0..3 {
                        if let Some(normal) =
                            span.facet_vertex_normal(triangle.uv, i, shift, dec)?
                            && normal_relation(direction, normal, cosine, true, dec)?.1
                        {
                            return Ok(NormalProof::Violation);
                        }
                    }
                }
                match normal_subcells(span, unit, &target, direction, cosine, 0, dec)? {
                    NormalProof::Proven => continue,
                    result => return Ok(result),
                }
            }
            // Whole-cell certificates can serve adjacent facets with different
            // trim triangles. Fall back to the triangle's tighter rectangle.
            let cell = local.map(dyadic_cell);
            let normal = span.normal_cell(cell, dec)?;
            if normals_agree(direction, normal, cosine, dec)?
                || (cell != local
                    && normals_agree(direction, span.normal_box(local, dec)?, cosine, dec)?)
            {
                continue;
            }
            return Ok(NormalProof::Unresolved);
        }
        Ok(NormalProof::Proven)
    }
}

fn vertex_residual(position: [f64; 3], enclosed: [Interval; 3]) -> Interval {
    position
        .into_iter()
        .zip(enclosed)
        .fold(point(0.0), |sum, (x, y)| {
            let delta = point(x) - y;
            sum + point(delta.lo().abs().max(delta.hi().abs()))
        })
}

fn normal_subcells(
    span: &Span,
    uv: [Interval; 2],
    target: &NormalFootprint,
    direction: [Interval; 3],
    cosine: Option<Interval>,
    depth: u8,
    dec: &mut Decider,
) -> Result<NormalProof, Refusal> {
    const MAX_DEPTH: u8 = 14;
    if depth == MAX_DEPTH {
        let clipped = std::array::from_fn(|axis| {
            let range = target.rectangle[axis];
            Interval::new(uv[axis].lo().max(range.lo()), uv[axis].hi().min(range.hi()))
                .expect("the shared cell meets the facet rectangle")
        });
        return Ok(
            if normals_agree(direction, span.normal_box(clipped, dec)?, cosine, dec)? {
                NormalProof::Proven
            } else {
                NormalProof::Unresolved
            },
        );
    }
    let axis = usize::from(uv[1].width() > uv[0].width());
    let range = uv[axis];
    let middle = range.mid();
    if middle <= range.lo() || middle >= range.hi() {
        return Ok(NormalProof::Unresolved);
    }
    for (lo, hi) in [(range.lo(), middle), (middle, range.hi())] {
        if hi < target.rectangle[axis].lo() || lo > target.rectangle[axis].hi() {
            continue;
        }
        let mut child = uv;
        child[axis] = Interval::new(lo, hi).unwrap();
        let outside_triangle = |dec: &mut Decider| {
            let Some(sides) = &target.sides else {
                return Ok(false);
            };
            dec.charge_work(6)?;
            Ok::<_, Refusal>(
                sides
                    .iter()
                    .any(|[x, y, constant]| (*x * child[0] + *y * child[1] + *constant).hi() < 0.0),
            )
        };
        if target.clip_first && outside_triangle(dec)? {
            continue;
        }
        let angular = cosine.is_some_and(|value| value.hi() > 0.0);
        let normal = if angular {
            span.angular_cone(child, dec)?
        } else {
            span.shared_normal(child, dec)?
        };
        let (accepts, rejects) = normal_relation(direction, normal, cosine, angular, dec)?;
        if !accepts {
            if !target.clip_first && outside_triangle(dec)? {
                continue;
            }
            if let Some(sides) = &target.sides
                && rejects
            {
                let middle = child.map(|range| point(range.mid()));
                dec.charge_work(6)?;
                if sides.iter().all(|[x, y, constant]| {
                    (*x * middle[0] + *y * middle[1] + *constant).lo() >= 0.0
                }) {
                    return Ok(NormalProof::Violation);
                }
            }
            match normal_subcells(span, child, target, direction, cosine, depth + 1, dec)? {
                NormalProof::Proven => {}
                result => return Ok(result),
            }
        }
    }
    Ok(NormalProof::Proven)
}

// The rolling band has large weight ratios within a quarter turn. Extract
// smaller exact patches so the same derivative proof uses local weights.
fn split_spans(surface: &BSplineSurface, dec: &mut Decider) -> Result<Vec<RationalSpan>, Refusal> {
    let mut pieces = Vec::new();
    for indices in surface.span_indices() {
        let ranges = std::array::from_fn::<_, 2, _>(|axis| {
            let lo = surface.knots()[axis][indices[axis]].clone();
            let hi = surface.knots()[axis][indices[axis] + 1].clone();
            let step = (&hi - &lo).mul_pow2(-2);
            std::array::from_fn::<_, 5, _>(|i| &lo + &(&step * &Rational::from_i64(i as i64)))
        });
        for u in ranges[0].windows(2) {
            for v in ranges[1].windows(2) {
                let patch = surface.bezier_patch(indices, [[&u[0], &u[1]], [&v[0], &v[1]]], dec)?;
                pieces.extend(rational_spans(&patch, dec)?);
            }
        }
    }
    Ok(pieces)
}

fn subtract_shift(bounds: [[f64; 2]; 2], shift: Uv) -> [[f64; 2]; 2] {
    std::array::from_fn(|i| bounds[i].map(|value| value - shift[i]))
}

fn closed_boundaries(surface: &BSplineSurface, axis: usize) -> bool {
    use catallax_shapes::{BSplineCurve, SplineBoundary::*};
    let sides = if axis == 0 {
        [UMin, UMax]
    } else {
        [VMin, VMax]
    };
    let [first, last] = sides.map(|side| BSplineCurve::boundary(surface, side));
    first
        .weights()
        .zip(last.weights())
        .all(|(a, b)| a.same_definition(b))
        && first.controls().zip(last.controls()).all(|(a, b)| {
            [&a.x, &a.y, &a.z]
                .into_iter()
                .zip([&b.x, &b.y, &b.z])
                .all(|(a, b)| a.same_definition(b))
        })
}

pub(super) fn checked_join_scale(
    surface: &BSplineSurface,
    dec: &mut Decider,
) -> Result<Interval, Refusal> {
    use catallax_numbers::{Outcome, Real};
    use catallax_shapes::{BSplineCurve, SplineBoundary};

    let [first, last] = [SplineBoundary::UMin, SplineBoundary::UMax]
        .map(|side| BSplineCurve::boundary(surface, side));
    let scale = last.weight(0) / first.weight(0);
    let mut same = |a: &Real, b: &Real| -> Result<(), Refusal> {
        if a.same_definition(b) {
            return Ok(());
        }
        match dec.sign(&(a - b)) {
            Outcome::ProvenCoincident(_) => Ok(()),
            Outcome::Proven(_) => Err(failure("periodic spline join has different controls")),
            Outcome::NotExact(reason) => Err(reason.into()),
            Outcome::Refused(reason) => Err(reason),
        }
    };
    for (a, b) in first.weights().zip(last.weights()) {
        same(&(a * &scale), b)?;
    }
    for (a, b) in first.controls().zip(last.controls()) {
        for (a, b) in [(&a.x, &b.x), (&a.y, &b.y), (&a.z, &b.z)] {
            same(a, b)?;
        }
    }
    Ok(dec.enclose_at(&scale, super::BITS)?.to_interval())
}

fn normals_agree(
    facet: [Interval; 3],
    normal: [Interval; 3],
    cosine: Option<Interval>,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    normal_relation(facet, normal, cosine, false, dec).map(|(accepts, _)| accepts)
}

fn normal_relation(
    facet: [Interval; 3],
    normal: [Interval; 3],
    cosine: Option<Interval>,
    counterexample: bool,
    dec: &mut Decider,
) -> Result<(bool, bool), Refusal> {
    dec.charge_work(3)?;
    let dot = facet[0] * normal[0] + facet[1] * normal[1] + facet[2] * normal[2];
    if dot.hi() <= 0.0 {
        return Ok((false, counterexample));
    }
    if dot.lo() <= 0.0 && !counterexample {
        return Ok((false, false));
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
        let accepts = dot.lo() > 0.0 && (cross_sq * cosine_sq).hi() <= required.lo();
        let rejects = if counterexample && !accepts && cosine.lo() > 0.0 {
            dec.charge_work(7)?;
            let cosine_sq = point(cosine.lo()).square();
            let required = (point(1.0) - cosine_sq) * point(dot.hi()).square();
            (cross_sq * cosine_sq).lo() > required.hi()
        } else {
            false
        };
        Ok((accepts, rejects))
    } else {
        Ok((dot.lo() > 0.0, false))
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
    fn periodic_lift_checks_both_spans_at_the_join() {
        // A square section extruded along z, with a non-unit source period.
        // The facet has x = 0 and the periodic join has x = 1.
        let mut dec = Decider::default();
        let controls = [[1, 0], [0, 1], [-1, 0], [0, -1], [1, 0]]
            .into_iter()
            .flat_map(|[x, y]| [0, 1].map(|z| Point3::new(x.into(), y.into(), z.into())))
            .collect();
        let surface = BSplineSurface::new(
            [1, 1],
            [
                [2, 2, 3, 4, 5, 6, 6].map(Rational::from_i64).to_vec(),
                [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
            ],
            [5, 2],
            controls,
            vec![Real::one(); 10],
            &mut dec,
        )
        .unwrap();
        let chart = SplineChart::new(
            &surface,
            [true, false],
            [None; 2],
            false,
            None,
            SplineNormals::Carrier,
            &mut dec,
        )
        .unwrap();
        for u in [0.0, 1.0] {
            assert!(chart.evaluate([u, 0.5], &mut dec).unwrap()[0].contains(1.0));
        }
        for shift in [0.0, 1.0] {
            let triangle = [
                [-0.25 + shift, 0.0],
                [0.25 + shift, 0.0],
                [-0.25 + shift, 1.0],
            ];
            assert!(chart.chord(&triangle, &mut dec).unwrap().hi() >= 1.0);
        }
    }

    #[test]
    fn bisected_bezier_span_encloses_the_same_weighted_point() {
        let knots = [0, 0, 1, 1].map(Rational::from_i64).to_vec();
        let controls = [[0, 0], [0, 2], [2, 0], [2, 2]]
            .map(|[x, y]| Point3::new(Real::from_i64(x), Real::from_i64(y), Real::zero()));
        let mut dec = Decider::default();
        let surface = BSplineSurface::new(
            [1, 1],
            [knots.clone(), knots],
            [2, 2],
            controls.to_vec(),
            [1, 1, 1, 2].map(Real::from_i64).to_vec(),
            &mut dec,
        )
        .unwrap();
        let span = rational_spans(&surface, &mut dec).unwrap().pop().unwrap();
        let [left, _] = span.bisect(0, &mut dec).unwrap();
        let [_, upper] = left.bisect(1, &mut dec).unwrap();
        let value = upper.evaluate([0.5, 0.5], &mut dec).unwrap();
        for (bound, numerator) in value.into_iter().zip([14, 30, 0]) {
            let expected = Interval::from_rational(&Rational::from_ratio(numerator, 19).unwrap());
            assert!(bound.lo() <= expected.lo() && expected.hi() <= bound.hi());
        }
    }

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
        let chart = SplineChart::new(
            &surface,
            [false; 2],
            [None; 2],
            false,
            None,
            SplineNormals::Carrier,
            &mut dec,
        )
        .unwrap();
        let options = crate::MeshOptions::default();
        let half_angle = Interval::from_rational(
            &crate::controls::cosine_limit(&Rational::from_ratio(13, 100).unwrap(), &dec).unwrap(),
        );
        let edge = |a: f64, b: f64| [a, (a + b) * 0.5, b].map(|u| [u, 0.0, 0.0]);
        assert!(
            chart
                .edge_within(edge(0.0, 1.0), [None; 3], None, &options, None, &mut dec)
                .unwrap()
        );
        assert!(
            !chart
                .edge_within(
                    edge(0.0, 1.0),
                    [None; 3],
                    None,
                    &options,
                    Some(half_angle),
                    &mut dec
                )
                .unwrap()
        );
        assert!(
            chart
                .edge_within(
                    edge(0.0, 0.0625),
                    [None; 3],
                    None,
                    &options,
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
        let chart = SplineChart::new(
            &surface,
            [false; 2],
            [None; 2],
            false,
            None,
            SplineNormals::Carrier,
            &mut dec,
        )
        .unwrap();
        let bound = chart
            .chord(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], &mut dec)
            .unwrap();
        assert!(bound.hi() >= 0.5);
        assert!(bound.hi() < 0.5 + f64::EPSILON.sqrt());
    }
}
