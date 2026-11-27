use std::{collections::BTreeMap, f64::consts::TAU};

use catallax_shapes::{Curve, SurfaceId};
use catallax_topology::{EdgeId, Model, Sense, VertexId};

use super::{DisplayFrame, DisplayOptions, DisplaySurface, Xyz, evaluate::*, spline};

pub(super) struct Reader<'a> {
    pub model: &'a Model,
    pub origin: Xyz,
    pub options: &'a DisplayOptions,
    pub vertices: &'a BTreeMap<VertexId, Xyz>,
    pub surfaces: &'a BTreeMap<SurfaceId, DisplaySurface>,
}

impl Reader<'_> {
    pub(super) fn edge(&self, id: EdgeId, incident: &[SurfaceId], max: usize) -> Option<Vec<Xyz>> {
        if max < 2 {
            return None;
        }
        let edge = self.model.edge(id);
        let reverse = edge.sense == Sense::Reversed;
        let mut ends = [edge.start, edge.end].map(|v| v.map(|v| self.vertices[&v]));
        if reverse {
            ends.reverse();
        }
        let origin = self.origin;
        let tolerance = self.options.chord_tolerance;
        let geom = self.model.geom();
        let mut points = match geom.curve(edge.curve) {
            Curve::Line(_) => vec![ends[0]?, ends[1]?],
            Curve::Circle(c) => {
                let f = frame(&c.center, &c.normal, &c.ref_dir, origin)?;
                let r = positive(libm::sqrt(c.radius_sq.to_f64()))?;
                round_edge(f, [r, r], ends, tolerance, max)?
            }
            Curve::Ellipse(c) => {
                let f = frame(&c.center, &c.normal, &c.major_dir, origin)?;
                let r = [&c.major_sq, &c.minor_sq].map(|r| libm::sqrt(r.to_f64()));
                round_edge(f, [positive(r[0])?, positive(r[1])?], ends, tolerance, max)?
            }
            Curve::BSpline(c) => {
                spline::SplineEdge::new(c, origin)?.samples(ends, tolerance, max)?
            }
            Curve::Intersection(c) => {
                let [a, b] = c.surfaces.map(|i| self.surfaces.get(&i));
                let [s, t] = c.surfaces.map(|i| geom.surface(i));
                let witness = sub(point(&geom.position(c.witness))?, origin);
                Intersection::new(Projector::new(a?, Some(s))?, Projector::new(b?, Some(t))?)
                    .trace(ends, witness, 1.0, tolerance, max)?
            }
            Curve::PatchSection(c) => {
                let (carrier, [start, end]) = c.first_cell();
                let seed_net = spline::net(carrier, origin)?;
                let mut seed = spline::Patch::new(&seed_net)?;
                let uv = |a: &[catallax_numbers::Real; 2]| a.each_ref().map(|u| u.to_f64());
                let witness = seed.at(uv(start));
                let next = seed.at(uv(end));
                let source = c.surface().surface();
                let cutting = describe(&source, origin)?;
                let patch = incident.iter().find_map(|id| match self.surfaces.get(id) {
                    Some(s @ DisplaySurface::BSpline(_)) => Some((s, geom.surface(*id))),
                    _ => None,
                });
                let owned;
                let (patch, original) = if let Some((patch, source)) = patch {
                    (patch, Some(source))
                } else {
                    owned = DisplaySurface::BSpline(seed_net);
                    (&owned, None)
                };
                let mut trace = Intersection::new(
                    Projector::new(&cutting, Some(&source))?,
                    Projector::new(patch, original)?,
                );
                let (_, tangent) = trace.correct(witness, tolerance * 0.01)?;
                let sense = if dot(tangent, sub(next, witness)) < 0.0 {
                    -1.0
                } else {
                    1.0
                };
                trace.trace(ends, witness, sense, tolerance, max)?
            }
            Curve::PatchPair(c) => {
                let [a, b] = c.carriers();
                let a = DisplaySurface::BSpline(spline::net(a, origin)?);
                let b = DisplaySurface::BSpline(spline::net(b, origin)?);
                let [u, v, _, _] = c.first_parameters().each_ref().map(|u| u.to_f64());
                let mut a = Projector::new(&a, None)?;
                let witness = a.jet([u, v]).0;
                Intersection::new(a, Projector::new(&b, None)?)
                    .trace(ends, witness, 1.0, tolerance, max)?
            }
            _ => return None,
        };
        if reverse {
            points.reverse();
        }
        for p in &mut points {
            let [x, y, z] = p.map(finite);
            *p = [x?, y?, z?];
        }
        Some(points)
    }
}

fn round_edge(
    f: DisplayFrame,
    radii: [f64; 2],
    ends: [Option<Xyz>; 2],
    tolerance: f64,
    max: usize,
) -> Option<Vec<Xyz>> {
    let uv = |p| {
        let d = sub(p, f.origin);
        angle(dot(d, f.y) / radii[1], dot(d, f.x) / radii[0])
    };
    let start = ends[0].map_or(0.0, uv);
    let sweep = ends[1].map_or(1.0, |p| {
        let d = uv(p) - start;
        if d > 0.0 { d } else { d + 1.0 }
    });
    let at = |u| {
        let (s, c) = libm::sincos(u * TAU);
        add(
            f.origin,
            add(scale(f.x, radii[0] * c), scale(f.y, radii[1] * s)),
        )
    };
    // The second derivative in radians is bounded by the largest semi-axis.
    // Limit the turn independently so a coarse ring still preserves its chart.
    let step = (libm::sqrt(8.0 * tolerance / radii[0].max(radii[1])) / TAU).min(0.125);
    let n = libm::ceil(sweep / step) as usize;
    if n == 0 || n >= max {
        return None;
    }
    let first = ends[0].unwrap_or_else(|| at(start));
    let mut points = Vec::with_capacity(n + 1);
    points.push(first);
    for i in 1..n {
        points.push(at(start + sweep * i as f64 / n as f64));
    }
    points.push(ends[1].unwrap_or(first));
    Some(points)
}

pub(super) fn sample_span(
    at: &mut impl FnMut(f64) -> Xyz,
    interval: [f64; 2],
    ends: [Xyz; 2],
    tolerance: f64,
    out: &mut impl FnMut(f64, Xyz) -> Option<()>,
) -> Option<()> {
    fn visit(
        at: &mut impl FnMut(f64) -> Xyz,
        [a, b]: [f64; 2],
        [p, q]: [Xyz; 2],
        tolerance: f64,
        out: &mut impl FnMut(f64, Xyz) -> Option<()>,
        depth: u32,
    ) -> Option<()> {
        let mid = a.midpoint(b);
        let m = at(mid);
        let error = |t, x| norm(sub(x, add(p, scale(sub(q, p), t))));
        if error(0.5, m) <= tolerance
            && error(0.25, at(a.midpoint(mid))) <= tolerance
            && error(0.75, at(mid.midpoint(b))) <= tolerance
        {
            return out(b, q);
        }
        if depth >= 32 || mid <= a || mid >= b {
            return None;
        }
        visit(at, [a, mid], [p, m], tolerance, out, depth + 1)?;
        visit(at, [mid, b], [m, q], tolerance, out, depth + 1)
    }
    visit(at, interval, ends, tolerance, out, 0)
}

struct Intersection<'a> {
    supports: [Projector<'a>; 2],
    seeds: [Option<[f64; 2]>; 2],
}

impl<'a> Intersection<'a> {
    fn new(a: Projector<'a>, b: Projector<'a>) -> Self {
        Self {
            supports: [a, b],
            seeds: [None, None],
        }
    }

    // Newton correction in the plane of the normals; it does not move along
    // the local tangent. Neither this correction nor the trace claims a bound.
    fn correct(&mut self, mut p: Xyz, accuracy: f64) -> Option<(Xyz, Xyz)> {
        let mut seeds = self.seeds;
        for _ in 0..20 {
            let mut normals = [[0.0; 3]; 2];
            let mut residuals = [0.0; 2];
            for (i, s) in self.supports.iter_mut().enumerate() {
                let uv = s.project(p, seeds[i]);
                seeds[i] = Some(uv);
                let uv = s.canonical(uv);
                let (q, [du, dv]) = s.jet(uv);
                normals[i] = unit(cross(du, dv))?;
                residuals[i] = dot(normals[i], sub(q, p));
            }
            let [a, b] = normals;
            if residuals.iter().all(|v| v.abs() <= accuracy) {
                self.seeds = seeds;
                return Some((p, unit(cross(a, b))?));
            }
            let [x, y] = solve_gram(a, b, residuals)?;
            p = add(p, add(scale(a, x), scale(b, y)));
        }
        None
    }

    fn trace(
        mut self,
        ends: [Option<Xyz>; 2],
        witness: Xyz,
        sense: f64,
        tolerance: f64,
        max: usize,
    ) -> Option<Vec<Xyz>> {
        let start = ends[0].unwrap_or(witness);
        let end = ends[1].unwrap_or(start);
        let accuracy = tolerance * 0.001;
        let (mut p, mut tangent) = self.correct(start, accuracy)?;
        tangent = scale(tangent, sense);
        let mut result = vec![start];
        let mut step = (norm(sub(end, start)) * 0.125).max(tolerance);
        let mut departed = start != end;
        for _ in 0..max {
            let target = sub(end, p);
            // Join the retained endpoint only on the forward side of the
            // branch, and only when the remaining chord meets the visual goal.
            if departed && norm(target) <= step * 1.5 && dot(target, tangent) >= -accuracy {
                let middle = add(p, scale(target, 0.5));
                if let Some((q, t)) = self.correct(middle, accuracy)
                    && norm(sub(q, middle)) <= tolerance * 0.5
                    && dot(scale(t, sense), tangent) > 0.8
                {
                    result.push(end);
                    return (result.len() <= max).then_some(result);
                }
            }
            let mut next = None;
            for _ in 0..32 {
                let guess = add(p, scale(tangent, step));
                if let Some((q, t)) = self.correct(guess, accuracy) {
                    let t = scale(t, sense);
                    if dot(t, tangent) >= 0.98 && norm(sub(q, guess)) <= tolerance * 0.5 {
                        next = Some((q, t));
                        break;
                    }
                }
                step *= 0.5;
                if step <= accuracy {
                    return None;
                }
            }
            (p, tangent) = next?;
            if departed && norm(sub(p, end)) <= accuracy {
                result.push(end);
                return (result.len() <= max).then_some(result);
            }
            if result.len() == max {
                return None;
            }
            result.push(p);
            departed |= norm(sub(p, start)) > tolerance * 2.0;
            step *= 1.5;
        }
        None
    }
}
