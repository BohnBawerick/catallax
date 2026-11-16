//! Numerical single-edge fillets between retained spline faces.
//!
//! The source and every continuation remain numerical bodies. The fitted
//! carrier, support pcurves and cap curves carry sampled checks, not exact
//! incidence or a whole-body approximation bound.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{BSplineCurve, Curve, NumericalFilletSpine, NumericalSplineFillet, Point3, RationalAtlas, SplineSupport, Surface};
use catallax_topology::{Sense, numerical::{NumericalBody, NumericalDefinition, NumericalEdge, NumericalFace, NumericalNativeTrim, NumericalTrim, NumericalVertex}};

fn bad(detail: impl Into<String>) -> Refusal { Refusal::new("fillet.numerical-body", detail) }
fn q(n: i64) -> Rational { Rational::from_i64(n) }
fn flip(s: Sense) -> Sense { if s == Sense::Same { Sense::Reversed } else { Sense::Same } }
fn scalar(value: &Real, dec: &Decider) -> Result<Rational, Refusal> { Ok(dec.enclose_at(value, 128)?.mid()) }
fn rational(value: f64) -> Result<Rational, Refusal> { Rational::from_f64_exact(value).ok_or_else(|| bad("non-finite numerical fillet parameter")) }
fn point_at(curve: &Curve, t: &Rational, dec: &Decider) -> Result<[Rational; 3], Refusal> {
    match curve {
        Curve::BSpline(curve) => Ok(curve.evaluate_box([t, t], 128, dec)?.map(|v| v.mid())),
        Curve::Line(line) => {
            let p = line.origin.plus(&line.dir.scaled(&Real::from_rational(t.clone())));
            Ok([scalar(&p.x, dec)?, scalar(&p.y, dec)?, scalar(&p.z, dec)?])
        }
        _ => Err(bad("a numerical spline fillet needs native line or spline border curves")),
    }
}
fn uv_at(curve: &BSplineCurve, t: &Rational, dec: &Decider) -> Result<[Rational; 2], Refusal> {
    let p = curve.evaluate_box([t, t], 128, dec)?;
    Ok([p[0].mid(), p[1].mid()])
}
fn plane_uv(surface: &Surface, p: &[Rational; 3], dec: &Decider) -> Result<[Rational; 2], Refusal> {
    let p = p.each_ref().map(|v| Ball::from_rational(v, 128));
    Ok(surface.parameters_of(&p, dec)?.ok_or_else(|| bad("end-plane parameters are unresolved"))?.map(|v| v.mid()))
}
fn atlas(face: &NumericalFace, dec: &mut Decider) -> Result<RationalAtlas, Refusal> {
    match &face.surface {
        Surface::BSpline(SplineSupport::RationalAtlas(atlas)) => Ok(atlas.clone()),
        Surface::BSpline(SplineSupport::RationalPatch(patch)) => {
            let patch = patch.extended(patch.domain().clone(), dec)?;
            let [atlas] = <[_; 1]>::try_from(RationalAtlas::partition(patch.carrier().clone(), [false, false], dec)?)
                .map_err(|_| bad("a spline fillet support must be one connected logical face"))?;
            Ok(atlas)
        }
        _ => Err(bad("the selected edge must join two retained rational spline faces")),
    }
}
fn plane(face: &NumericalFace, dec: &Decider) -> Result<([Rational; 3], Rational), Refusal> {
    let Surface::Plane(plane) = &face.surface else { return Err(bad("this spline fillet requires planar end faces")); };
    let normal = [scalar(&plane.normal.x, dec)?, scalar(&plane.normal.y, dec)?, scalar(&plane.normal.z, dec)?];
    let origin = [scalar(&plane.origin.x, dec)?, scalar(&plane.origin.y, dec)?, scalar(&plane.origin.z, dec)?];
    let offset = (0..3).fold(q(0), |sum, i| sum + &normal[i] * &origin[i]);
    Ok((normal, offset))
}
fn walked(edge: &NumericalEdge, sense: Sense) -> [usize; 2] {
    if sense == Sense::Same { edge.vertices } else { [edge.vertices[1], edge.vertices[0]] }
}
fn coedge_uv(trim: &NumericalTrim) -> [[Rational; 2]; 2] {
    if trim.sense == Sense::Same { trim.uv.clone() } else { [trim.uv[1].clone(), trim.uv[0].clone()] }
}
fn along(uv: &[[Rational; 2]; 2]) -> Result<usize, Refusal> {
    let changes = [(&uv[1][0] - &uv[0][0]).abs(), (&uv[1][1] - &uv[0][1]).abs()];
    if changes.iter().all(|v| v.sign() == Sign::Zero) {
        return Err(bad("selected spline edge has no nonempty parameter graph"));
    }
    Ok(usize::from(changes[1] >= changes[0]))
}
fn directions(
    atlas: &RationalAtlas, trim: &NumericalTrim, sense: Sense,
    tolerance: &Rational, dec: &mut Decider,
) -> Result<[[f64; 3]; 2], Refusal> {
    let at = [
        (&trim.uv[0][0] + &trim.uv[1][0]).mul_pow2(-1),
        (&trim.uv[0][1] + &trim.uv[1][1]).mul_pow2(-1),
    ];
    // The full tensor source can include a sharp neighboring logical face.
    // Its derivative hull is not either face's one-sided normal.
    let chart = atlas.charts().iter().find(|chart| {
        (0..2).all(|axis| at[axis] >= chart.patch().domain()[axis][0]
            && at[axis] <= chart.patch().domain()[axis][1])
    }).ok_or_else(|| bad("selected edge has no retained support chart"))?;
    let at_real = at.each_ref().map(|value| Real::from_rational(value.clone()));
    let partials = chart.patch().derivatives_at(&at_real, dec)?
        .map(|value| [value.x.to_f64(), value.y.to_f64(), value.z.to_f64()]);
    let [du, dv] = partials;
    let point = chart.patch().point_at(&at_real, dec)?;
    let free = 1 - along(&coedge_uv(trim))?;
    let bounds = domain(atlas);
    let low = (&at[free] - &bounds[free][0]).abs() <= (&at[free] - &bounds[free][1]).abs();
    let mut border = at.clone();
    border[free] = bounds[free][usize::from(!low)].clone();
    let border_real = border.each_ref().map(|value| Real::from_rational(value.clone()));
    let border_point = chart.patch().point_at(&border_real, dec)?;
    let differences = [
        &point.x - &border_point.x,
        &point.y - &border_point.y,
        &point.z - &border_point.z,
    ];
    let mut miss = q(0);
    for difference in &differences {
        let gap = dec.enclose_at(difference, 128)?;
        miss = miss + gap.lo().abs().max(gap.hi().abs());
    }
    if miss > *tolerance {
        return Err(bad("selected edge is not on a support border within the source tolerance"));
    }
    let inward_parameter = if low { 1.0 } else { -1.0 };
    let sign = if sense == Sense::Same { 1.0 } else { -1.0 };
    let into = partials[free].map(|value| inward_parameter * value);
    let outward = [
        sign * (du[1] * dv[2] - du[2] * dv[1]),
        sign * (du[2] * dv[0] - du[0] * dv[2]),
        sign * (du[0] * dv[1] - du[1] * dv[0]),
    ];
    dec.charge_work(64)?;
    Ok([into, outward])
}

fn end_face(d: &NumericalDefinition, vertex: usize, sides: [usize; 2]) -> Result<usize, Refusal> {
    let others: Vec<_> = d.faces.iter().enumerate().filter(|(i, face)| {
        !sides.contains(i) && face.loops.iter().flatten().any(|trim| d.edges[trim.edge].vertices.contains(&vertex))
    }).map(|(i, _)| i).collect();
    let [face] = <[_; 1]>::try_from(others).map_err(|_| bad("a spline fillet end must meet exactly one other face"))?;
    Ok(face)
}
fn neighbour(d: &NumericalDefinition, face: usize, edge: usize, vertex: usize) -> Result<usize, Refusal> {
    let found: Vec<_> = d.faces[face].loops.iter().flatten().filter(|trim| trim.edge != edge && d.edges[trim.edge].vertices.contains(&vertex)).map(|trim| trim.edge).collect();
    let [edge] = <[_; 1]>::try_from(found).map_err(|_| bad("a spline fillet end needs one adjacent border edge per side"))?;
    Ok(edge)
}
fn project(edge: &NumericalEdge, point: &[Rational; 3], tolerance: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    let mut bounds = edge.parameters.clone();
    bounds.sort();
    let target = point.each_ref().map(Rational::to_f64_nearest);
    let (mut lo, mut hi) = (bounds[0].to_f64_nearest(), bounds[1].to_f64_nearest());
    let mut t = 0.5 * (lo + hi);
    for _ in 0..48 {
        dec.check_deadline()?;
        let at = rational(t)?;
        let (value, tangent) = match &edge.curve {
            Curve::BSpline(curve) => {
                let jet = curve.jet_box([&at, &at], 64, dec)?;
                (jet.value.map(|v| v.mid().to_f64_nearest()), jet.first.map(|v| v.mid().to_f64_nearest()))
            }
            Curve::Line(line) => (point_at(&edge.curve, &at, dec)?.map(|v| v.to_f64_nearest()), [line.dir.x.to_f64(), line.dir.y.to_f64(), line.dir.z.to_f64()]),
            _ => return Err(bad("a fillet border is not a native line or spline")),
        };
        dec.charge_work(32)?;
        let squared = tangent.iter().map(|v| v * v).sum::<f64>();
        if !squared.is_finite() || squared <= f64::MIN_POSITIVE { return Err(bad("a fillet border has a singular tangent")); }
        let gradient = (0..3).map(|i| (value[i] - target[i]) * tangent[i]).sum::<f64>();
        if gradient < 0.0 { lo = t; } else { hi = t; }
        let next = (t - gradient / squared).clamp(lo, hi);
        if (next - t).abs() <= 1e-14 * t.abs().max(1.0) { t = next; break; }
        t = next;
    }
    let parameter = rational(t)?;
    let actual = point_at(&edge.curve, &parameter, dec)?;
    let error = (0..3).fold(q(0), |sum, i| sum + (&actual[i] - &point[i]).abs());
    if error > *tolerance { return Err(bad(format!("fillet rail misses its end border by {error} mm against {tolerance} mm"))); }
    if parameter <= bounds[0] || parameter >= bounds[1] { return Err(bad("fillet rail does not end strictly inside its adjacent border")); }
    Ok(parameter)
}
fn domain(atlas: &RationalAtlas) -> [[Rational; 2]; 2] {
    let mut result = atlas.charts()[0].patch().domain().clone();
    for chart in &atlas.charts()[1..] {
        for (axis, side) in result.iter_mut().enumerate() {
            side[0] = side[0].clone().min(chart.patch().domain()[axis][0].clone());
            side[1] = side[1].clone().max(chart.patch().domain()[axis][1].clone());
        }
    }
    result
}

/// Numerical fillet result, retaining the original source for radius edits.
#[derive(Clone, Debug)]
pub struct NumericalFilletBuilt {
    source: Arc<NumericalBody>,
    body: NumericalBody,
    edge: String,
    key: String,
    radius: Rational,
    tolerance_mm: Rational,
    fillet: NumericalSplineFillet,
    constrained_parameter_controls: usize,
}
impl NumericalFilletBuilt {
    /// Complete admitted numerical body, suitable for the numerical STEP writer.
    pub fn body(&self) -> &NumericalBody { &self.body }
    /// Numerical carrier and its separate sampled rolling-ball checks.
    pub fn fillet(&self) -> &NumericalSplineFillet { &self.fillet }
    /// Radius used for this build, in millimetres.
    pub fn radius(&self) -> &Rational { &self.radius }
    /// Pcurve controls constrained to the physical support box or end borders.
    /// These constraints do not move the three-dimensional rail curves.
    pub fn constrained_parameter_controls(&self) -> usize { self.constrained_parameter_controls }
    /// Rebuilds from the same original source and preserves operation-derived keys.
    pub fn edit_radius(&self, radius: Rational, dec: &mut Decider) -> Result<Self, Refusal> {
        guarded(|| build_impl(self.source.clone(), &self.edge, &self.key, radius, self.tolerance_mm.clone(), Some(&self.body), dec))
    }
}

/// Rounds a selected native edge of a numerical body. Source conversion and
/// approximate STEP remain explicit caller operations. No exact model is returned.
pub fn build(
    source: NumericalBody,
    edge: &str,
    key: &str,
    radius: Rational,
    tolerance_mm: Rational,
    dec: &mut Decider,
) -> Result<NumericalFilletBuilt, Refusal> {
    guarded(|| build_impl(Arc::new(source), edge, key, radius, tolerance_mm, None, dec))
}

fn guarded<T>(run: impl FnOnce() -> Result<T, Refusal>) -> Result<T, Refusal> {
    if catallax_numbers::operation_budget_active() { return run(); }
    let (result, work) = catallax_numbers::with_operation_budget(
        catallax_topology::numerical::DEFAULT_WORK, run,
    );
    if work.exhausted {
        Err(Refusal::new("numbers.operation-budget", "numerical fillet exhausted its operation allowance"))
    } else {
        result
    }
}

type Bindings = BTreeMap<(String, String), BSplineCurve>;
type Parameters = BTreeMap<(usize, usize), [Rational; 2]>;

fn constrained_rail(
    curve: &BSplineCurve,
    interval: &[Rational; 2],
    endpoints: &[[Rational; 2]; 2],
    bounds: &[[Rational; 2]; 2],
    dec: &mut Decider,
) -> Result<(BSplineCurve, usize), Refusal> {
    let forward = interval[0] < interval[1];
    let (range, ends) = if forward { (interval.clone(), endpoints.clone()) } else { ([interval[1].clone(), interval[0].clone()], [endpoints[1].clone(), endpoints[0].clone()]) };
    let degree = curve.degree();
    let mut controls = Vec::new();
    let mut weights = Vec::new();
    let mut knots = vec![range[0].clone(); degree + 1];
    let mut changed = 0;
    for (span, old) in curve.spans() {
        let lo = old[0].clone().max(range[0].clone());
        let hi = old[1].clone().min(range[1].clone());
        if lo >= hi { continue; }
        let data = curve.span_bezier_controls_on(span, [&lo, &hi], dec)?;
        let first = controls.is_empty();
        for (index, h) in data.into_iter().enumerate() {
            if !first && index == 0 { continue; }
            let weight = scalar(&h[3], dec)?;
            let mut uv = [scalar(&(&h[0] / &h[3]), dec)?, scalar(&(&h[1] / &h[3]), dec)?];
            for axis in 0..2 {
                let next = uv[axis].clone().max(bounds[axis][0].clone()).min(bounds[axis][1].clone());
                changed += usize::from(next != uv[axis]);
                uv[axis] = next;
            }
            controls.push(Point3::new(uv[0].clone().into(), uv[1].clone().into(), Real::zero()));
            weights.push(Real::from_rational(weight));
        }
        if hi < range[1] { knots.extend(std::iter::repeat_n(hi, degree)); }
    }
    if controls.len() < degree + 1 { return Err(bad("a finite support pcurve has no clipped span")); }
    let last = controls.len() - 1;
    for (index, uv) in [(0, &ends[0]), (last, &ends[1])] {
        let next = Point3::new(uv[0].clone().into(), uv[1].clone().into(), Real::zero());
        changed += usize::from(controls[index].x.as_rational() != Some(uv[0].clone()) || controls[index].y.as_rational() != Some(uv[1].clone()));
        controls[index] = next;
    }
    knots.extend(std::iter::repeat_n(range[1].clone(), degree + 1));
    Ok((BSplineCurve::new(degree, knots, controls, weights, dec)?, changed))
}
fn parameter_line(range: [&Rational; 2], side: usize, dec: &mut Decider) -> Result<BSplineCurve, Refusal> {
    BSplineCurve::new(1, vec![range[0].clone(), range[0].clone(), range[1].clone(), range[1].clone()],
        range.map(|tau| Point3::new(tau.clone().into(), q(side as i64).into(), Real::zero())).to_vec(), vec![Real::one(); 2], dec)
}
fn new_trim(edge: usize, sense: Sense, uv: [[Rational; 2]; 2]) -> NumericalTrim { NumericalTrim { edge, sense, uv } }
fn compact(d: &mut NumericalDefinition) {
    let used: BTreeSet<_> = d.faces.iter().flat_map(|face| face.loops.iter().flatten().map(|trim| trim.edge)).collect();
    let mut edge_map = vec![usize::MAX; d.edges.len()];
    let old = std::mem::take(&mut d.edges);
    for (index, edge) in old.into_iter().enumerate() {
        if used.contains(&index) { edge_map[index] = d.edges.len(); d.edges.push(edge); }
    }
    for trim in d.faces.iter_mut().flat_map(|face| face.loops.iter_mut().flatten()) { trim.edge = edge_map[trim.edge]; }
    let used: BTreeSet<_> = d.edges.iter().flat_map(|edge| edge.vertices).collect();
    let mut vertex_map = vec![usize::MAX; d.vertices.len()];
    let old = std::mem::take(&mut d.vertices);
    for (index, vertex) in old.into_iter().enumerate() {
        if used.contains(&index) { vertex_map[index] = d.vertices.len(); d.vertices.push(vertex); }
    }
    for edge in &mut d.edges { edge.vertices = edge.vertices.map(|v| vertex_map[v]); }
}

fn at_stage(stage: &str, mut refusal: Refusal) -> Refusal {
    refusal.detail = format!("numerical spline fillet {stage}: {}", refusal.detail);
    refusal
}

#[allow(clippy::too_many_arguments)]
fn build_impl(
    source: Arc<NumericalBody>, edge_key: &str, key: &str, radius: Rational,
    tolerance_mm: Rational, previous: Option<&NumericalBody>, dec: &mut Decider,
) -> Result<NumericalFilletBuilt, Refusal> {
    if key.is_empty() || radius.sign() != Sign::Positive || tolerance_mm.sign() != Sign::Positive || tolerance_mm > Rational::parse("0.001").unwrap() {
        return Err(bad("a numerical fillet needs an operation key, positive radius and stated tolerance at most 1 micrometre"));
    }
    let original = source.definition();
    let selected = original.edges.iter().position(|edge| edge.key == edge_key).ok_or_else(|| bad("selected numerical edge does not exist"))?;
    let input = &original.edges[selected];
    let uses: Vec<_> = original.faces.iter().enumerate().flat_map(|(face, f)| f.loops.iter().enumerate().flat_map(move |(boundary, row)| row.iter().enumerate().filter_map(move |(index, trim)| (trim.edge == selected).then_some((face, boundary, index))))).collect();
    let [(first, row_a, at_a), (second, row_b, at_b)] = <[_; 2]>::try_from(uses).map_err(|_| bad("selected numerical edge does not have exactly two uses"))?;
    if row_a != 0 || row_b != 0 { return Err(bad("this numerical spline fillet requires an edge on both outer boundaries")); }
    let sides = [first, second];
    let trims = [original.faces[first].loops[0][at_a].clone(), original.faces[second].loops[0][at_b].clone()];
    if trims[0].sense == trims[1].sense { return Err(bad("selected edge uses do not oppose each other")); }
    let uv = trims.each_ref().map(coedge_uv);
    let axes = [along(&uv[0])?, along(&uv[1])?];
    let physical = [atlas(&original.faces[first], dec)?, atlas(&original.faces[second], dec)?];
    let turns = [
        directions(&physical[0], &trims[0], original.faces[first].sense, &source.options().tolerance_mm, dec)?,
        directions(&physical[1], &trims[1], original.faces[second].sense, &source.options().tolerance_mm, dec)?,
    ];
    let convex = |a: usize, b: usize| {
        (0..3).map(|i| turns[a][0][i] * turns[b][1][i]).sum::<f64>() < 0.0
    };
    if !convex(0, 1) || !convex(1, 0) {
        return Err(bad("this numerical spline fillet requires a convex selected edge"));
    }
    let extended = [
        NumericalFilletSpine::extended_support(&physical[0], axes[0], &Rational::from_ratio(1, 4).unwrap(), dec)
            .map_err(|error| at_stage("first support extension", error))?,
        NumericalFilletSpine::extended_support(&physical[1], axes[1], &Rational::from_ratio(1, 4).unwrap(), dec)
            .map_err(|error| at_stage("second support extension", error))?,
    ];
    let mut range = [uv[0][0][axes[0]].clone(), uv[0][1][axes[0]].clone()];
    range.sort();
    let margin = (&range[1] - &range[0]).mul_pow2(-3);
    range = [&range[0] - &margin, &range[1] + &margin];
    let share = (&range[0] - &uv[0][0][axes[0]]) / (&uv[0][1][axes[0]] - &uv[0][0][axes[0]]);
    let second_seed: [Rational; 2] = std::array::from_fn(|i| &uv[1][0][i] + &(&(&uv[1][1][i] - &uv[1][0][i]) * &share));
    let first_seed = &uv[0][0][1 - axes[0]]
        + &(&(&uv[0][1][1 - axes[0]] - &uv[0][0][1 - axes[0]]) * &share);
    let seed = [first_seed, second_seed[0].clone(), second_seed[1].clone()];
    let distances = sides.map(|face| if original.faces[face].sense == Sense::Same { -radius.clone() } else { radius.clone() });
    let spine = NumericalFilletSpine::trace(&extended, distances, axes[0], range, seed, tolerance_mm.clone(), dec)
        .map_err(|error| at_stage("centre reconstruction", error))?;
    let fillet = NumericalSplineFillet::fit(spine, dec)
        .map_err(|error| at_stage("carrier reconstruction", error))?;
    let caps = [end_face(original, input.vertices[0], sides)?, end_face(original, input.vertices[1], sides)?];
    if caps[0] == caps[1] { return Err(bad("a numerical spline fillet needs two distinct end faces")); }
    let planes = [plane(&original.faces[caps[0]], dec)?, plane(&original.faces[caps[1]], dec)?];
    let cuts = [
        fillet.plane_cut(planes[0].0.clone(), planes[0].1.clone(), dec)
            .map_err(|error| at_stage("first end-plane cut", error))?,
        fillet.plane_cut(planes[1].0.clone(), planes[1].1.clone(), dec)
            .map_err(|error| at_stage("second end-plane cut", error))?,
    ];
    let mut d = original.clone();
    d.name = format!("{} / {key}", original.name);
    let mut bindings: Bindings = source.native_trims().iter().map(|binding| ((binding.face.clone(), binding.edge.clone()), binding.curve.clone())).collect();
    let mut parameters: Parameters = BTreeMap::new();
    for (face, f) in original.faces.iter().enumerate() {
        for trim in f.loops.iter().flatten() {
            for (vertex, uv) in walked(&original.edges[trim.edge], trim.sense).into_iter().zip(&trim.uv) { parameters.insert((face, vertex), uv.clone()); }
        }
    }
    let mut rail_vertices = [[0usize; 2]; 2];
    for side in 0..2 {
        for end in 0..2 {
            let tau = &cuts[end].tau_ends()[side];
            let curve = fillet.rail(side).ok_or_else(|| bad("fitted numerical rail is missing"))?;
            let point = curve.evaluate_box([tau, tau], 128, dec)?.map(|v| v.mid());
            rail_vertices[side][end] = d.vertices.len();
            d.vertices.push(NumericalVertex { key: format!("{key}/rail/{}/end/{}", original.faces[sides[side]].key, original.faces[caps[end]].key), position: point });
        }
    }
    let mut rail_uv: [[[Rational; 2]; 2]; 2] = std::array::from_fn(|_| std::array::from_fn(|_| std::array::from_fn(|_| q(0))));
    for side in 0..2 {
        for end in 0..2 {
            let border = neighbour(original, sides[side], selected, input.vertices[end])?;
            let vertex = rail_vertices[side][end];
            let t = project(&original.edges[border], &d.vertices[vertex].position, &source.options().tolerance_mm, dec)?;
            let slot = usize::from(original.edges[border].vertices[1] == input.vertices[end]);
            d.edges[border].vertices[slot] = vertex;
            d.edges[border].parameters[slot] = t.clone();
            let binding = bindings.get(&(original.faces[sides[side]].key.clone(), original.edges[border].key.clone())).ok_or_else(|| bad("a spline side border needs its native parameter curve"))?;
            let uv = uv_at(binding, &t, dec)?;
            rail_uv[side][end] = uv.clone();
            parameters.insert((sides[side], vertex), uv);
            parameters.insert((caps[end], vertex), plane_uv(&original.faces[caps[end]].surface, &d.vertices[vertex].position, dec)?);
        }
    }
    let mut rails = [0; 2];
    let mut constrained_parameter_controls = 0;
    for side in 0..2 {
        let edge_key = format!("{key}/rail/{}", original.faces[sides[side]].key);
        let interval = [cuts[0].tau_ends()[side].clone(), cuts[1].tau_ends()[side].clone()];
        let (curve, changed) = constrained_rail(fillet.parameter_rail(side).ok_or_else(|| bad("support pcurve is missing"))?, &interval, &rail_uv[side], &domain(&physical[side]), dec)?;
        constrained_parameter_controls += changed;
        bindings.insert((original.faces[sides[side]].key.clone(), edge_key.clone()), curve);
        rails[side] = d.edges.len();
        d.edges.push(NumericalEdge { key: edge_key, curve: Curve::BSpline(fillet.rail(side).unwrap().clone()), vertices: rail_vertices[side], parameters: interval });
    }
    d.faces[first].loops[0][at_a].edge = rails[0];
    d.faces[second].loops[0][at_b].edge = rails[1];
    let mut cut_edges = [0; 2];
    for end in 0..2 {
        cut_edges[end] = d.edges.len();
        d.edges.push(NumericalEdge { key: format!("{key}/cap/{}", original.faces[caps[end]].key), curve: Curve::BSpline(cuts[end].curve().clone()), vertices: [rail_vertices[0][end], rail_vertices[1][end]], parameters: [q(0), q(1)] });
        let face = caps[end];
        let loops = std::mem::take(&mut d.faces[face].loops);
        for boundary in loops {
            let mut replacement = Vec::with_capacity(boundary.len() + 1);
            for (i, trim) in boundary.iter().enumerate() {
                replacement.push(trim.clone());
                let a = walked(&d.edges[trim.edge], trim.sense)[1];
                let next = &boundary[(i + 1) % boundary.len()];
                let b = walked(&d.edges[next.edge], next.sense)[0];
                if a != b {
                    let endpoints = d.edges[cut_edges[end]].vertices;
                    let sense = if [a, b] == endpoints { Sense::Same } else if [b, a] == endpoints { Sense::Reversed } else { return Err(bad("end-plane loop has a gap unrelated to the selected fillet")); };
                    replacement.push(new_trim(cut_edges[end], sense, [parameters[&(face, a)].clone(), parameters[&(face, b)].clone()]));
                }
            }
            d.faces[face].loops.push(replacement);
        }
    }
    let face_index = d.faces.len();
    let face_key = format!("{key}/blend/{edge_key}");
    let mut sampled_carriers = source.sampled_carriers().clone();
    sampled_carriers.insert(face_key.clone());
    for side in 0..2 {
        let curve = fillet.rail(side).unwrap();
        bindings.insert((face_key.clone(), d.edges[rails[side]].key.clone()), parameter_line(curve.domain(), side, dec)?);
        for end in 0..2 { parameters.insert((face_index, rail_vertices[side][end]), [cuts[end].tau_ends()[side].clone(), q(side as i64)]); }
    }
    for end in 0..2 { bindings.insert((face_key.clone(), d.edges[cut_edges[end]].key.clone()), cuts[end].parameters().clone()); }
    let first_use = (rails[0], flip(trims[0].sense));
    let reaches = usize::from(first_use.1 == Sense::Same);
    let path = [first_use, (cut_edges[reaches], Sense::Same), (rails[1], flip(trims[1].sense)), (cut_edges[1 - reaches], Sense::Reversed)];
    let boundary = path.into_iter().map(|(edge, sense)| {
        let ends = walked(&d.edges[edge], sense);
        new_trim(edge, sense, ends.map(|vertex| parameters[&(face_index, vertex)].clone()))
    }).collect();
    let carrier = fillet.surface().source();
    let native = carrier.domain();
    let tau = (native[0][0] + native[0][1]).mul_pow2(-1);
    let half = Rational::from_ratio(1, 2).unwrap();
    let jet = carrier.jet_box([&tau, &tau], [&half, &half], 64, dec)?;
    let du = jet.first[0].each_ref().map(|v| v.mid().to_f64_nearest());
    let dv = jet.first[1].each_ref().map(|v| v.mid().to_f64_nearest());
    let normal = [du[1] * dv[2] - du[2] * dv[1], du[2] * dv[0] - du[0] * dv[2], du[0] * dv[1] - du[1] * dv[0]];
    let centre = fillet.spine().centre_curve().evaluate_box([&tau, &tau], 128, dec)?;
    let direction = (0..3).map(|i| normal[i] * (jet.value[i].mid().to_f64_nearest() - centre[i].mid().to_f64_nearest())).sum::<f64>();
    if !direction.is_finite() || direction == 0.0 { return Err(bad("numerical fillet material orientation is unresolved")); }
    d.faces.push(NumericalFace { key: face_key, surface: Surface::BSpline(SplineSupport::RationalAtlas(fillet.surface().clone())), sense: if direction > 0.0 { Sense::Same } else { Sense::Reversed }, loops: vec![boundary] });
    let mut native_trims = Vec::new();
    for (face, f) in d.faces.iter_mut().enumerate() {
        for (boundary, row) in f.loops.iter_mut().enumerate() {
            for (use_index, trim) in row.iter_mut().enumerate() {
                let edge = &d.edges[trim.edge];
                trim.uv = walked(edge, trim.sense).map(|vertex| parameters[&(face, vertex)].clone());
                if let Some(curve) = bindings.get(&(f.key.clone(), edge.key.clone())) {
                    native_trims.push(NumericalNativeTrim { face: f.key.clone(), boundary, use_index, edge: edge.key.clone(), curve: curve.clone() });
                }
            }
        }
    }
    compact(&mut d);
    let body = previous.unwrap_or(&source).rebuild_with_sampled_carriers(d, native_trims, sampled_carriers, source.options().clone(), dec)
        .map_err(|error| at_stage("body native admission", error))?;
    Ok(NumericalFilletBuilt { source, body, edge: edge_key.into(), key: key.into(), radius, tolerance_mm, fillet, constrained_parameter_controls })
}
