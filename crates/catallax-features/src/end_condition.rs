//! Definition-native first-hit arrangements for bounded planar limits.
//! The arrangement is shared data, not a B-rep evaluation: field compilation
//! consumes its halfspaces directly, while the B-rep constructs its boundary.

use crate::{
    Profile,
    profile::{Arc, sign},
};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{
    BooleanKind, FaceRef, Operation, Origin, Recipe, Role, StepId, ThroughSide, Turn,
};
use catallax_shapes::{Plane, Point3, RPoint3, RVector3, Surface, Vector3};
use catallax_topology::{BodyId, Model};
use std::collections::{BTreeMap, BTreeSet};

type P = [Real; 2];
fn r(n: i64) -> Real {
    Real::from_i64(n)
}
fn bad(code: &'static str, detail: &str) -> Refusal {
    Refusal::new(code, detail)
}
fn sqrt(value: &Real) -> Real {
    value
        .rational_ref()
        .and_then(crate::extrude::sqrt)
        .map(Real::from)
        .unwrap_or_else(|| value.sqrt())
}
fn cross(a: &P, b: &P) -> Real {
    &a[0] * &b[1] - &a[1] * &b[0]
}
fn sub(a: &P, b: &P) -> P {
    [&a[0] - &b[0], &a[1] - &b[1]]
}
fn dot(a: &P, b: &P) -> Real {
    &a[0] * &b[0] + &a[1] * &b[1]
}
fn midpoint(a: &P, b: &P) -> P {
    std::array::from_fn(|i| &(&a[i] + &b[i]) / &r(2))
}
fn same(a: &P, b: &P, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(&(&a[0] - &b[0]), dec)? == Sign::Zero && sign(&(&a[1] - &b[1]), dec)? == Sign::Zero)
}
fn centroid(poly: &[P]) -> P {
    std::array::from_fn(|i| {
        poly.iter().fold(Real::zero(), |s, p| s + p[i].clone()) / Real::from_i64(poly.len() as i64)
    })
}
fn qpoint(p: &[Rational; 2]) -> P {
    p.clone().map(Real::from)
}
fn zero_vec() -> RVector3 {
    RVector3::zero()
}
fn vector(v: &catallax_recipe::Vector3) -> RVector3 {
    let [x, y, z] = v.0.each_ref().map(|n| n.value().clone());
    RVector3::new(x, y, z)
}
fn point(p: &catallax_recipe::Point3) -> Point3 {
    let [x, y, z] = p.0.each_ref().map(|n| Real::from(n.value().clone()));
    Point3::new(x, y, z)
}

#[derive(Clone, Debug)]
struct Line2([Real; 3]);
impl Line2 {
    fn at(&self, p: &P) -> Real {
        &self.0[0] * &p[0] + &self.0[1] * &p[1] + self.0[2].clone()
    }
    fn through(a: &P, b: &P) -> Self {
        let d = sub(b, a);
        Self([-&d[1], d[0].clone(), &d[1] * &a[0] - &d[0] * &a[1]])
    }
    fn normalised(mut self, dec: &mut Decider) -> Result<Option<Self>, Refusal> {
        let mut pivot = None;
        for c in &self.0[..2] {
            if sign(c, dec)? != Sign::Zero {
                pivot = Some(c.clone());
                break;
            }
        }
        let Some(pivot) = pivot else { return Ok(None) };
        for c in &mut self.0 {
            *c = &*c / &pivot;
        }
        Ok(Some(self))
    }
}
fn push_line(lines: &mut Vec<Line2>, line: Line2, dec: &mut Decider) -> Result<(), Refusal> {
    let Some(line) = line.normalised(dec)? else {
        return Ok(());
    };
    for old in lines.iter() {
        let mut equal = true;
        for i in 0..3 {
            if sign(&(&old.0[i] - &line.0[i]), dec)? != Sign::Zero {
                equal = false;
                break;
            }
        }
        if equal {
            return Ok(());
        }
    }
    lines.push(line);
    Ok(())
}
fn clip2(poly: &[P], line: &Line2, negative: bool, dec: &mut Decider) -> Result<Vec<P>, Refusal> {
    let mut out = Vec::with_capacity(poly.len() + 1);
    for i in 0..poly.len() {
        let a = &poly[i];
        let b = &poly[(i + 1) % poly.len()];
        let va = line.at(a);
        let vb = line.at(b);
        let sa = sign(&va, dec)?;
        let sb = sign(&vb, dec)?;
        let keep = |s| {
            if negative {
                s != Sign::Positive
            } else {
                s != Sign::Negative
            }
        };
        if keep(sa) {
            out.push(a.clone());
        }
        if (sa == Sign::Negative && sb == Sign::Positive)
            || (sa == Sign::Positive && sb == Sign::Negative)
        {
            let t = &va / &(&va - &vb);
            out.push(std::array::from_fn(|j| &a[j] + &(&t * &(&b[j] - &a[j]))));
        }
    }
    if out.len() < 3 {
        return Ok(vec![]);
    }
    let area = out.iter().enumerate().fold(Real::zero(), |s, (i, a)| {
        s + cross(a, &out[(i + 1) % out.len()])
    });
    if sign(&area, dec)? == Sign::Zero {
        out.clear();
    }
    Ok(out)
}
fn inside_polygon(poly: &[P], p: &P, strict: bool, dec: &mut Decider) -> Result<bool, Refusal> {
    for i in 0..poly.len() {
        let s = sign(
            &cross(
                &sub(&poly[(i + 1) % poly.len()], &poly[i]),
                &sub(p, &poly[i]),
            ),
            dec,
        )?;
        if s == Sign::Negative || (strict && s == Sign::Zero) {
            return Ok(false);
        }
    }
    Ok(true)
}
fn ordered(values: &mut Vec<Real>, dec: &mut Decider) -> Result<(), Refusal> {
    for i in 1..values.len() {
        let mut j = i;
        while j > 0 && sign(&(&values[j] - &values[j - 1]), dec)? == Sign::Negative {
            values.swap(j, j - 1);
            j -= 1;
        }
    }
    let mut unique = Vec::with_capacity(values.len());
    for v in values.drain(..) {
        if unique.last().is_none() || sign(&(&v - unique.last().unwrap()), dec)? != Sign::Zero {
            unique.push(v);
        }
    }
    *values = unique;
    Ok(())
}
fn line_circle(
    a: &P,
    d: &P,
    c: &P,
    radius_sq: &Real,
    dec: &mut Decider,
) -> Result<Vec<P>, Refusal> {
    let v = sub(a, c);
    let aa = dot(d, d);
    let b = dot(&v, d);
    let disc = &b.square() - &(&aa * &(&dot(&v, &v) - radius_sq));
    let s = sign(&disc, dec)?;
    if s == Sign::Negative {
        return Ok(vec![]);
    }
    let root = if s == Sign::Zero {
        Real::zero()
    } else {
        sqrt(&disc)
    };
    let mut out = Vec::with_capacity(2);
    for k in if s == Sign::Zero {
        &[1][..]
    } else {
        &[-1, 1][..]
    } {
        let t = &(-&b + &r(*k) * &root) / &aa;
        out.push(std::array::from_fn(|i| &a[i] + &(&d[i] * &t)));
    }
    Ok(out)
}
#[derive(Clone, Copy, PartialEq, Eq)]
enum RegionLocation {
    Inside,
    Outside,
    Boundary,
}
fn profile_location(
    profile: &Profile,
    p: &P,
    dec: &mut Decider,
) -> Result<RegionLocation, Refusal> {
    for e in profile.loops.iter().flat_map(|b| &b.segments) {
        let on = if let Some((c, rs, _)) = &e.circle {
            sign(
                &(&dot(&sub(p, &qpoint(c)), &sub(p, &qpoint(c))) - &Real::from(rs.clone())),
                dec,
            )? == Sign::Zero
                && e.contains(p, dec)?
        } else {
            e.contains(p, dec)?
        };
        if on {
            return Ok(RegionLocation::Boundary);
        }
    }
    let mut slope = 0i64;
    let d = loop {
        let d = [r(1), r(slope)];
        let mut good = true;
        for e in profile.loops.iter().flat_map(|b| &b.segments) {
            if !e.full {
                for q in [&e.start, &e.end] {
                    if sign(&cross(&sub(&qpoint(q), p), &d), dec)? == Sign::Zero {
                        good = false;
                        break;
                    }
                }
            }
            if let Some((c, rs, _)) = &e.circle {
                let det = cross(&sub(&qpoint(c), p), &d);
                if sign(
                    &(&det.square() - &(&Real::from(rs.clone()) * &dot(&d, &d))),
                    dec,
                )? == Sign::Zero
                {
                    good = false;
                }
            }
            if !good {
                break;
            }
        }
        if good {
            break d;
        }
        slope = slope.checked_add(1).ok_or_else(|| {
            bad(
                "features.precision",
                "profile ray direction index exhausted",
            )
        })?;
    };
    let mut crossings = 0usize;
    for e in profile.loops.iter().flat_map(|b| &b.segments) {
        if let Some((c, rs, _)) = &e.circle {
            for hit in line_circle(p, &d, &qpoint(c), &Real::from(rs.clone()), dec)? {
                if sign(&dot(&sub(&hit, p), &d), dec)? == Sign::Positive && e.contains(&hit, dec)? {
                    crossings += 1;
                }
            }
        } else {
            let a = qpoint(&e.start);
            let edge = sub(&qpoint(&e.end), &a);
            let det = cross(&d, &edge);
            if sign(&det, dec)? == Sign::Zero {
                continue;
            }
            let v = sub(&a, p);
            let t = &cross(&v, &edge) / &det;
            let u = &cross(&v, &d) / &det;
            if sign(&t, dec)? == Sign::Positive
                && sign(&u, dec)? == Sign::Positive
                && sign(&(&u - &r(1)), dec)? == Sign::Negative
            {
                crossings += 1;
            }
        }
    }
    Ok(if crossings % 2 == 1 {
        RegionLocation::Inside
    } else {
        RegionLocation::Outside
    })
}
fn arc_middle(e: &Arc, a: &P, b: &P, full: bool, dec: &mut Decider) -> Result<P, Refusal> {
    let Some((c, rs, turn)) = &e.circle else {
        return Ok(midpoint(a, b));
    };
    let c = qpoint(c);
    let va = sub(a, &c);
    let vb = sub(b, &c);
    let mut v = [&va[0] + &vb[0], &va[1] + &vb[1]];
    if full {
        v = [-&va[0], -&va[1]];
    } else if sign(&dot(&v, &v), dec)? == Sign::Zero {
        v = if *turn == Turn::Ccw {
            [-&va[1], va[0].clone()]
        } else {
            [va[1].clone(), -&va[0]]
        };
    } else {
        let order = sign(&cross(&va, &vb), dec)?;
        if (order == Sign::Negative && *turn == Turn::Ccw)
            || (order == Sign::Positive && *turn == Turn::Cw)
        {
            v = v.map(|x| -x);
        }
    }
    let factor = sqrt(&(&Real::from(rs.clone()) / &dot(&v, &v)));
    Ok(std::array::from_fn(|i| &c[i] + &(&v[i] * &factor)))
}
fn boundary_in_cell(profile: &Profile, poly: &[P], dec: &mut Decider) -> Result<bool, Refusal> {
    for e in profile.loops.iter().flat_map(|b| &b.segments) {
        let a = qpoint(&e.start);
        let b = qpoint(&e.end);
        let mut cuts = vec![a.clone()];
        if !e.full {
            cuts.push(b.clone());
        }
        for i in 0..poly.len() {
            let l = Line2::through(&poly[i], &poly[(i + 1) % poly.len()]);
            let d = sub(&poly[(i + 1) % poly.len()], &poly[i]);
            if let Some((c, rs, _)) = &e.circle {
                for p in line_circle(&poly[i], &d, &qpoint(c), &Real::from(rs.clone()), dec)? {
                    if e.contains(&p, dec)? {
                        cuts.push(p);
                    }
                }
            } else {
                let va = l.at(&a);
                let vb = l.at(&b);
                if sign(&(&va - &vb), dec)? != Sign::Zero {
                    let t = &va / &(&va - &vb);
                    if sign(&t, dec)? != Sign::Negative
                        && sign(&(&t - &r(1)), dec)? != Sign::Positive
                    {
                        cuts.push(std::array::from_fn(|j| &a[j] + &(&t * &(&b[j] - &a[j]))));
                    }
                }
            }
        }
        if let Some((c, _, turn)) = &e.circle {
            let c = qpoint(c);
            let start = sub(&a, &c);
            let flip = if *turn == Turn::Ccw { 1 } else { -1 };
            let half = |v: &P, dec: &mut Decider| -> Result<bool, Refusal> {
                let s = sign(&(&cross(&start, v) * &r(flip)), dec)?;
                Ok(s == Sign::Negative
                    || (s == Sign::Zero && sign(&dot(&start, v), dec)? == Sign::Negative))
            };
            for i in 1..cuts.len() {
                let mut j = i;
                while j > 0 {
                    let v = sub(&cuts[j], &c);
                    let w = sub(&cuts[j - 1], &c);
                    let hv = half(&v, dec)?;
                    let hw = half(&w, dec)?;
                    let before = if hv != hw {
                        !hv
                    } else {
                        sign(&(&cross(&v, &w) * &r(flip)), dec)? == Sign::Positive
                    };
                    if !before {
                        break;
                    }
                    cuts.swap(j, j - 1);
                    j -= 1;
                }
            }
        } else {
            let d = sub(&b, &a);
            for i in 1..cuts.len() {
                let mut j = i;
                while j > 0 && sign(&dot(&sub(&cuts[j], &cuts[j - 1]), &d), dec)? == Sign::Negative
                {
                    cuts.swap(j, j - 1);
                    j -= 1;
                }
            }
        }
        let mut unique = Vec::new();
        for p in cuts {
            if unique.last().is_none() || !same(&p, unique.last().unwrap(), dec)? {
                unique.push(p);
            }
        }
        let count = if e.full {
            unique.len()
        } else {
            unique.len().saturating_sub(1)
        };
        for i in 0..count {
            let p = arc_middle(
                e,
                &unique[i],
                &unique[(i + 1) % unique.len()],
                e.full && unique.len() == 1,
                dec,
            )?;
            if inside_polygon(poly, &p, true, dec)? {
                return Ok(true);
            }
        }
    }
    Ok(false)
}
fn profile_intersects(profile: &Profile, poly: &[P], dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(
        profile_location(profile, &centroid(poly), dec)? == RegionLocation::Inside
            || boundary_in_cell(profile, poly, dec)?,
    )
}
fn material_on_segment(
    profile: &Profile,
    a: &P,
    b: &P,
    dec: &mut Decider,
) -> Result<Option<P>, Refusal> {
    let d = sub(b, a);
    let norm = dot(&d, &d);
    let mut cuts = vec![r(0), r(1)];
    for e in profile.loops.iter().flat_map(|b| &b.segments) {
        if let Some((c, rs, _)) = &e.circle {
            for p in line_circle(a, &d, &qpoint(c), &Real::from(rs.clone()), dec)? {
                if e.contains(&p, dec)? {
                    cuts.push(&dot(&sub(&p, a), &d) / &norm);
                }
            }
        } else {
            let p = qpoint(&e.start);
            let v = sub(&qpoint(&e.end), &p);
            let det = cross(&d, &v);
            if sign(&det, dec)? != Sign::Zero {
                let u = &cross(&sub(&p, a), &d) / &det;
                if sign(&u, dec)? != Sign::Negative && sign(&(&u - &r(1)), dec)? != Sign::Positive {
                    cuts.push(&cross(&sub(&p, a), &v) / &det);
                }
            }
        }
    }
    ordered(&mut cuts, dec)?;
    for pair in cuts.windows(2) {
        let t = &(&pair[0] + &pair[1]) / &r(2);
        if sign(&t, dec)? != Sign::Positive || sign(&(&t - &r(1)), dec)? != Sign::Negative {
            continue;
        }
        let p = std::array::from_fn(|i| &a[i] + &(&d[i] * &t));
        if profile_location(profile, &p, dec)? == RegionLocation::Inside {
            return Ok(Some(p));
        }
    }
    Ok(None)
}

mod linear;

/// One cell of the exact first-hit height arrangement.
#[derive(Clone, Debug)]
pub struct CapCell {
    /// Counterclockwise polygon in the authored sketch coordinates.
    pub footprint: Vec<[Real; 2]>,
    /// The affine limiting support, oriented outward from the extrusion.
    pub limit: Plane,
    /// Which vertical walls can contribute interior discontinuity curtains.
    pub curtains: Vec<bool>,
    /// False only for construction fill proved disjoint from the profile interior.
    pub active: bool,
}

/// Definition-native reduction of a first-hit extrusion to planar height cells.
#[derive(Clone, Debug)]
pub struct Capped {
    pub(crate) profile: Profile,
    /// Raw authored sweep direction; line parameters do not normalise it.
    pub direction: RVector3,
    /// Height cells, including exterior construction fill for boundary assembly.
    pub cells: Vec<CapCell>,
    /// Certified result sphere about world zero.
    pub radius: Rational,
    /// Certified sphere for each active column before profile trimming.
    pub column_radius: Rational,
    extent: Rational,
}
#[derive(Clone)]
struct Hit {
    plane: usize,
}

impl Capped {
    /// Resolves trimmed limits and proves positive, uniformly bounded first hits.
    pub fn compile(
        recipe: &Recipe,
        target: &StepId,
        profile: &StepId,
        direction: Option<&catallax_recipe::Vector3>,
        limit: Option<&FaceRef>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        recipe
            .check()
            .map_err(|e| Refusal::new(e.code(), e.to_string()))?;
        let profile = linear::read_profile(recipe, profile, dec)?;
        let normal = crate::profile::rational_vector(
            &profile.pose,
            &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
        );
        let direction = direction.map(vector).unwrap_or_else(|| normal.clone());
        if direction.is_zero() {
            return Err(bad(
                "extrude.zero-direction",
                "extrusion direction must not be zero",
            ));
        }
        let projection = direction.dot(&normal);
        if projection.is_zero() {
            return Err(bad(
                "extrude.zero-thickness",
                "extrusion direction lies in the sketch plane",
            ));
        }
        let bound = linear::recipe_bound(recipe, target, dec)?;
        let profile_radius = crate::evaluate::profile_radius(&profile);
        let radius = bound.clone().max(profile_radius);
        let origin = crate::profile::rational_point(&profile.pose, &RPoint3::origin());
        let norm_origin = origin.x.abs() + origin.y.abs() + origin.z.abs();
        let norm_normal = normal.x.abs() + normal.y.abs() + normal.z.abs();
        let extent = ((&bound + &norm_origin) * norm_normal) / projection.abs() + Rational::one();
        let mut result = Self {
            profile,
            direction,
            cells: vec![],
            radius,
            column_radius: bound.clone(),
            extent,
        };
        if result.profile.is_empty() {
            return Ok(result);
        }
        let mut scene = linear::Scene::compile(recipe, target, dec)?;
        let boundary = scene.boundary(dec)?;
        let selected = if let Some(reference) = limit {
            boundary.resolve(&scene, reference, dec)?
        } else {
            (0..boundary.groups.len()).collect()
        };
        let mut bounds: Option<[[Rational; 2]; 2]> = None;
        for edge in result.profile.loops.iter().flat_map(|b| &b.segments) {
            for p in [&edge.start, &edge.end] {
                let b = bounds.get_or_insert_with(|| p.each_ref().map(|q| [q.clone(), q.clone()]));
                for i in 0..2 {
                    b[i][0] = b[i][0].clone().min(p[i].clone());
                    b[i][1] = b[i][1].clone().max(p[i].clone());
                }
            }
            if let Some((c, rs, _)) = &edge.circle {
                let radius = Ball::from_rational(rs, 128)
                    .sqrt()
                    .ok_or_else(|| bad("sketch.arc-radius", "profile radius is not positive"))?
                    .hi();
                let b = bounds.as_mut().unwrap();
                for i in 0..2 {
                    b[i][0] = b[i][0].clone().min(&c[i] - &radius);
                    b[i][1] = b[i][1].clone().max(&c[i] + &radius);
                }
            }
        }
        let bounds = bounds
            .unwrap()
            .map(|[lo, hi]| [lo - Rational::one(), hi + Rational::one()]);
        let start = vec![
            [
                Real::from(bounds[0][0].clone()),
                Real::from(bounds[1][0].clone()),
            ],
            [
                Real::from(bounds[0][1].clone()),
                Real::from(bounds[1][0].clone()),
            ],
            [
                Real::from(bounds[0][1].clone()),
                Real::from(bounds[1][1].clone()),
            ],
            [
                Real::from(bounds[0][0].clone()),
                Real::from(bounds[1][1].clone()),
            ],
        ];
        let mut column_radius = bound;
        for q in &start {
            let p = result.world(q);
            let mut upper = Rational::zero();
            for x in [&p.x, &p.y, &p.z] {
                let b = if let Some(q) = x.rational_ref() {
                    Ball::from_rational(q, 128)
                } else {
                    dec.enclose_budgeted_at(x, 128)?
                };
                upper = upper + b.lo().abs().max(b.hi().abs());
            }
            column_radius = column_radius.max(upper);
        }
        result.column_radius = column_radius;
        let mut lines = Vec::new();
        let mut heights = Vec::new();
        let mut physical = BTreeSet::new();
        for patch in &boundary.patches {
            physical.insert(patch.plane);
        }
        for plane in physical {
            let p = &scene.planes[plane];
            let divisor = p.normal.dot(&result.direction.to_real());
            if sign(&divisor, dec)? == Sign::Zero {
                continue;
            }
            let at = |q: &P| result.world(q).minus(&p.origin).dot(&p.normal);
            let c = -at(&[r(0), r(0)]) / divisor.clone();
            let a = &(-at(&[r(1), r(0)]) / divisor.clone()) - &c;
            let b = &(-at(&[r(0), r(1)]) / divisor) - &c;
            let height = Line2([a, b, c]);
            push_line(&mut lines, height.clone(), dec)?;
            heights.push(height);
        }
        for i in 0..heights.len() {
            for other in &heights[..i] {
                push_line(
                    &mut lines,
                    Line2(std::array::from_fn(|j| &heights[i].0[j] - &other.0[j])),
                    dec,
                )?;
            }
        }
        for group in &boundary.groups {
            let mut edges: BTreeMap<(usize, usize), usize> = BTreeMap::new();
            for i in group {
                let patch = &boundary.patches[*i];
                for j in 0..patch.vertices.len() {
                    let a = patch.vertices[j];
                    let b = patch.vertices[(j + 1) % patch.vertices.len()];
                    *edges.entry((a.min(b), a.max(b))).or_default() += 1;
                }
            }
            for ((a, b), count) in edges {
                if count == 1 {
                    push_line(
                        &mut lines,
                        Line2::through(
                            &result.project(&boundary.points[a]),
                            &result.project(&boundary.points[b]),
                        ),
                        dec,
                    )?;
                }
            }
        }
        let mut polygons = vec![start];
        for line in &lines {
            let mut split = Vec::new();
            for poly in polygons {
                let mut negative = false;
                let mut positive = false;
                for p in &poly {
                    match sign(&line.at(p), dec)? {
                        Sign::Negative => negative = true,
                        Sign::Positive => positive = true,
                        Sign::Zero => (),
                    }
                }
                if negative && positive {
                    for side in [true, false] {
                        let p = clip2(&poly, line, side, dec)?;
                        if !p.is_empty() {
                            split.push(p);
                        }
                    }
                } else {
                    split.push(poly);
                }
            }
            polygons = split;
        }
        let forward = projection.sign() == Sign::Positive;
        let fill = result.start_plane().origin.plus(
            &result
                .direction
                .to_real()
                .scaled(&Real::from(result.extent.clone())),
        );
        let fill_plane = Plane {
            origin: fill,
            normal: normal
                .scaled(&Rational::from_i64(if forward { 1 } else { -1 }))
                .into(),
            ref_dir: result.profile.pose.apply_vector(
                &RVector3::new(Rational::one(), Rational::zero(), Rational::zero()).to_real(),
            ),
        };
        for poly in polygons {
            let active = profile_intersects(&result.profile, &poly, dec)?;
            let limit = if active {
                let hit = result
                    .hit(
                        &scene,
                        &boundary,
                        &selected,
                        limit.is_some(),
                        &centroid(&poly),
                        dec,
                    )?
                    .ok_or_else(|| {
                        bad(
                            "extrude.no-finite-cap",
                            "a profile ray has no positive finite limiting boundary",
                        )
                    })?;
                let mut plane = scene.planes[hit.plane].clone();
                if sign(&plane.normal.dot(&result.direction.to_real()), dec)? == Sign::Negative {
                    plane.normal = plane.normal.scaled(&r(-1));
                }
                plane
            } else {
                fill_plane.clone()
            };
            result.cells.push(CapCell {
                curtains: vec![false; poly.len()],
                footprint: poly,
                limit,
                active,
            });
        }
        let mut vertices: Vec<P> = Vec::new();
        let mut edges: BTreeMap<(usize, usize), Vec<(usize, usize)>> = BTreeMap::new();
        for (i, cell) in result.cells.iter().enumerate() {
            let mut ids = Vec::with_capacity(cell.footprint.len());
            for p in &cell.footprint {
                ids.push(point_index(&mut vertices, p, dec)?);
            }
            for j in 0..ids.len() {
                let a = ids[j];
                let b = ids[(j + 1) % ids.len()];
                edges.entry((a.min(b), a.max(b))).or_default().push((i, j));
            }
        }
        for q in &vertices {
            if profile_location(&result.profile, q, dec)? == RegionLocation::Inside
                && result
                    .hit(&scene, &boundary, &selected, limit.is_some(), q, dec)?
                    .is_none()
            {
                return Err(bad(
                    "extrude.no-finite-cap",
                    "an interior profile ray on a first-hit branch boundary has no positive cap",
                ));
            }
        }
        for ((a, b), uses) in edges {
            if let Some(q) = material_on_segment(&result.profile, &vertices[a], &vertices[b], dec)?
            {
                if result
                    .hit(&scene, &boundary, &selected, limit.is_some(), &q, dec)?
                    .is_none()
                {
                    return Err(bad(
                        "extrude.no-finite-cap",
                        "an interior profile ray along a first-hit branch boundary has no positive cap",
                    ));
                }
                for (i, j) in uses {
                    result.cells[i].curtains[j] = true;
                }
            }
        }
        Ok(result)
    }
    /// An exact world-space profile-plane point.
    pub fn world(&self, q: &[Real; 2]) -> Point3 {
        self.profile
            .pose
            .apply_point(&Point3::new(q[0].clone(), q[1].clone(), r(0)))
    }
    fn project(&self, p: &Point3) -> P {
        let n = self.profile.pose.apply_vector(
            &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()).to_real(),
        );
        let origin = self.profile.pose.apply_point(&Point3::origin());
        let t = &p.minus(&origin).dot(&n) / &self.direction.to_real().dot(&n);
        let p = p.plus(&self.direction.to_real().scaled(&(-t)));
        let local = self.profile.pose.inverse().apply_point(&p);
        [local.x, local.y]
    }
    fn start_plane(&self) -> Plane {
        let mut normal = self.profile.pose.apply_vector(
            &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()).to_real(),
        );
        let raw = self.direction.to_real();
        if normal
            .dot(&raw)
            .rational_ref()
            .is_some_and(|q| q.sign() == Sign::Positive)
        {
            normal = normal.scaled(&r(-1));
        }
        Plane {
            origin: self.profile.pose.apply_point(&Point3::origin()),
            normal,
            ref_dir: self.profile.pose.apply_vector(
                &RVector3::new(Rational::one(), Rational::zero(), Rational::zero()).to_real(),
            ),
        }
    }
    fn height(&self, cell: &CapCell, q: &P) -> Real {
        -self
            .world(q)
            .minus(&cell.limit.origin)
            .dot(&cell.limit.normal)
            / self.direction.to_real().dot(&cell.limit.normal)
    }
    /// Oriented halfspaces and real names of one active convex height cell.
    pub fn cell_planes(
        &self,
        cell: &CapCell,
        dec: &mut Decider,
    ) -> Result<Vec<(Plane, Option<Role>)>, Refusal> {
        let mut planes = vec![
            (self.start_plane(), Some(Role::Start)),
            (cell.limit.clone(), Some(Role::End)),
        ];
        let normal = self.profile.pose.apply_vector(
            &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()).to_real(),
        );
        let forward = sign(&normal.dot(&self.direction.to_real()), dec)? == Sign::Positive;
        for i in 0..cell.footprint.len() {
            let a = self.world(&cell.footprint[i]);
            let b = self.world(&cell.footprint[(i + 1) % cell.footprint.len()]);
            let mut normal = b.minus(&a).cross(&self.direction.to_real());
            if !forward {
                normal = normal.scaled(&r(-1));
            }
            planes.push((
                Plane {
                    origin: a,
                    normal,
                    ref_dir: self.direction.to_real(),
                },
                cell.curtains[i].then_some(Role::End),
            ));
        }
        Ok(planes)
    }
    fn hit(
        &self,
        scene: &linear::Scene,
        boundary: &linear::Boundary,
        selected: &BTreeSet<usize>,
        to: bool,
        q: &P,
        dec: &mut Decider,
    ) -> Result<Option<Hit>, Refusal> {
        let origin = self.world(q);
        let direction = self.direction.to_real();
        let mut roots = vec![r(0)];
        let mut physical = BTreeSet::new();
        for patch in &boundary.patches {
            physical.insert(patch.plane);
        }
        for plane in &physical {
            let p = &scene.planes[*plane];
            let divisor = p.normal.dot(&direction);
            if sign(&divisor, dec)? != Sign::Zero {
                let t = -origin.minus(&p.origin).dot(&p.normal) / divisor;
                if sign(&t, dec)? == Sign::Positive {
                    roots.push(t);
                }
            }
        }
        ordered(&mut roots, dec)?;
        if to {
            let t = if roots.len() > 1 {
                &roots[1] / &r(2)
            } else {
                r(1)
            };
            if boundary.on_selected(scene, selected, &origin.plus(&direction.scaled(&t)), dec)? {
                return Err(bad(
                    "extrude.no-finite-cap",
                    "positive hits accumulate at zero without a least positive cap",
                ));
            }
        }
        for i in 1..roots.len() {
            let p = origin.plus(&direction.scaled(&roots[i]));
            let hit = if to {
                boundary.on_selected(scene, selected, &p, dec)?
            } else {
                let before =
                    origin.plus(&direction.scaled(&(&(&roots[i - 1] + &roots[i]) / &r(2))));
                let after = origin.plus(&direction.scaled(&if i + 1 < roots.len() {
                    &(&roots[i] + &roots[i + 1]) / &r(2)
                } else {
                    &roots[i] + &r(1)
                }));
                let a = boundary.location(scene, &before, dec)?;
                let b = boundary.location(scene, &after, dec)?;
                matches!(
                    (a, b),
                    (RegionLocation::Inside, RegionLocation::Outside)
                        | (RegionLocation::Outside, RegionLocation::Inside)
                )
            };
            if hit {
                for plane in &physical {
                    let s = &scene.planes[*plane];
                    if sign(&s.normal.dot(&direction), dec)? != Sign::Zero
                        && sign(&p.minus(&s.origin).dot(&s.normal), dec)? == Sign::Zero
                    {
                        return Ok(Some(Hit { plane: *plane }));
                    }
                }
            }
        }
        Ok(None)
    }
}
fn point_index(points: &mut Vec<P>, p: &P, dec: &mut Decider) -> Result<usize, Refusal> {
    for (i, old) in points.iter().enumerate() {
        if same(old, p, dec)? {
            return Ok(i);
        }
    }
    let i = points.len();
    points.push(p.clone());
    Ok(i)
}

fn level_node(
    levels: &[Vec<Real>],
    offsets: &[usize],
    vertex: usize,
    height: &Real,
    dec: &mut Decider,
) -> Result<usize, Refusal> {
    for (i, value) in levels[vertex].iter().enumerate() {
        if sign(&(value - height), dec)? == Sign::Zero {
            return Ok(offsets[vertex] + i);
        }
    }
    Err(bad(
        "features.limit-arrangement",
        "a curtain endpoint has no canonical height vertex",
    ))
}
fn boundary_face(
    model: &mut Model,
    mut ids: Vec<usize>,
    plane: Plane,
    origins: Vec<Origin>,
    points: &[Point3],
    vertices: &[catallax_topology::VertexId],
    edges: &mut BTreeMap<(usize, usize), catallax_topology::EdgeId>,
) -> Option<catallax_topology::FaceId> {
    ids.dedup();
    if ids.len() > 1 && ids.first() == ids.last() {
        ids.pop();
    }
    if ids.len() < 3 {
        return None;
    }
    let mut path = Vec::with_capacity(ids.len());
    for i in 0..ids.len() {
        let a = ids[i];
        let b = ids[(i + 1) % ids.len()];
        if a == b {
            continue;
        }
        let key = (a.min(b), a.max(b));
        let edge = *edges.entry(key).or_insert_with(|| {
            let curve =
                model
                    .geom_mut()
                    .add_curve(catallax_shapes::Curve::Line(catallax_shapes::Line {
                        origin: points[key.0].clone(),
                        dir: points[key.1].minus(&points[key.0]),
                    }));
            model.add_edge(catallax_topology::Edge {
                curve,
                sense: catallax_topology::Sense::Same,
                start: Some(vertices[key.0]),
                end: Some(vertices[key.1]),
            })
        });
        path.push((
            edge,
            if a < b {
                catallax_topology::Sense::Same
            } else {
                catallax_topology::Sense::Reversed
            },
        ));
    }
    Some(crate::extrude::face(
        model,
        Surface::Plane(plane),
        catallax_topology::Sense::Same,
        vec![path],
        None,
        origins,
    ))
}

impl Capped {
    /// Constructs the complete height-envelope boundary, then trims it by the
    /// authored profile. Construction fill and caps cannot survive that trim.
    pub fn build(
        &self,
        model: &mut Model,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<BodyId, Refusal> {
        use catallax_topology::{Body, BodyKind, Lump, Shell, Vertex};
        if self.profile.is_empty() {
            return Ok(catallax_topology::empty(model));
        }
        let mut footprint = Vec::new();
        let mut cells = Vec::with_capacity(self.cells.len());
        let mut edge_uses: BTreeMap<(usize, usize), Vec<(usize, usize)>> = BTreeMap::new();
        for (i, cell) in self.cells.iter().enumerate() {
            let mut ids = Vec::with_capacity(cell.footprint.len());
            for p in &cell.footprint {
                ids.push(point_index(&mut footprint, p, dec)?);
            }
            for j in 0..ids.len() {
                let a = ids[j];
                let b = ids[(j + 1) % ids.len()];
                edge_uses
                    .entry((a.min(b), a.max(b)))
                    .or_default()
                    .push((i, j));
            }
            cells.push(ids);
        }
        let mut levels = vec![vec![r(0)]; footprint.len()];
        for (cell, ids) in self.cells.iter().zip(&cells) {
            for vertex in ids {
                let height = self.height(cell, &footprint[*vertex]);
                if sign(&height, dec)? == Sign::Negative {
                    return Err(bad(
                        "features.limit-arrangement",
                        "a positive first-hit branch became negative on its closure",
                    ));
                }
                levels[*vertex].push(height);
            }
        }
        let mut offsets = Vec::with_capacity(footprint.len());
        let mut points = Vec::new();
        let mut vertices = Vec::new();
        for (i, values) in levels.iter_mut().enumerate() {
            ordered(values, dec)?;
            offsets.push(points.len());
            let base = self.world(&footprint[i]);
            for height in values {
                let p = base.plus(&self.direction.to_real().scaled(height));
                let point = model
                    .geom_mut()
                    .add_point(catallax_shapes::PointDef::Explicit(p.clone()));
                vertices.push(model.add_vertex(Vertex { point }));
                points.push(p);
            }
        }
        let start = self.start_plane();
        let normal = self.profile.pose.apply_vector(
            &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()).to_real(),
        );
        let forward = sign(&normal.dot(&self.direction.to_real()), dec)? == Sign::Positive;
        let names = |role: Option<Role>| {
            role.into_iter()
                .map(|role| Origin {
                    step: step.clone(),
                    role,
                })
                .collect()
        };
        let mut edges = BTreeMap::new();
        let mut faces = Vec::new();
        for (cell, ids) in self.cells.iter().zip(&cells) {
            let mut bottom = Vec::with_capacity(ids.len());
            let mut top = Vec::with_capacity(ids.len());
            for id in ids {
                bottom.push(level_node(&levels, &offsets, *id, &r(0), dec)?);
                top.push(level_node(
                    &levels,
                    &offsets,
                    *id,
                    &self.height(cell, &footprint[*id]),
                    dec,
                )?);
            }
            if forward {
                bottom.reverse();
            } else {
                top.reverse();
            }
            if let Some(face) = boundary_face(
                model,
                bottom,
                start.clone(),
                names(Some(Role::Start)),
                &points,
                &vertices,
                &mut edges,
            ) {
                faces.push(face);
            }
            if let Some(face) = boundary_face(
                model,
                top,
                cell.limit.clone(),
                names(cell.active.then_some(Role::End)),
                &points,
                &vertices,
                &mut edges,
            ) {
                faces.push(face);
            }
        }
        for (_, uses) in edge_uses {
            if uses.len() > 2 {
                return Err(bad(
                    "features.limit-arrangement",
                    "a footprint edge has more than two incident cells",
                ));
            }
            let (i, j) = uses[0];
            let a = cells[i][j];
            let b = cells[i][(j + 1) % cells[i].len()];
            let first = [
                self.height(&self.cells[i], &footprint[a]),
                self.height(&self.cells[i], &footprint[b]),
            ];
            let second = if uses.len() == 2 {
                let other = uses[1].0;
                [
                    self.height(&self.cells[other], &footprint[a]),
                    self.height(&self.cells[other], &footprint[b]),
                ]
            } else {
                [r(0), r(0)]
            };
            let order = sign(
                &(&(&first[0] + &first[1]) - &(&second[0] + &second[1])),
                dec,
            )?;
            if order == Sign::Zero {
                continue;
            }
            let (owner, edge, a, b, lo, hi) = if order == Sign::Positive {
                (i, j, a, b, second, first)
            } else {
                let (other, k) = uses[1];
                (
                    other,
                    k,
                    b,
                    a,
                    [first[1].clone(), first[0].clone()],
                    [second[1].clone(), second[0].clone()],
                )
            };
            let lo_a = level_node(&levels, &offsets, a, &lo[0], dec)?;
            let lo_b = level_node(&levels, &offsets, b, &lo[1], dec)?;
            let hi_a = level_node(&levels, &offsets, a, &hi[0], dec)?;
            let hi_b = level_node(&levels, &offsets, b, &hi[1], dec)?;
            let mut polygon = vec![lo_a, lo_b];
            for id in lo_b + 1..=hi_b {
                polygon.push(id);
            }
            polygon.push(hi_a);
            for id in (lo_a..hi_a).rev() {
                polygon.push(id);
            }
            if !forward {
                polygon.reverse();
            }
            let pa = self.world(&footprint[a]);
            let pb = self.world(&footprint[b]);
            let mut normal = pb.minus(&pa).cross(&self.direction.to_real());
            if !forward {
                normal = normal.scaled(&r(-1));
            }
            let plane = Plane {
                origin: pa,
                normal,
                ref_dir: self.direction.to_real(),
            };
            if let Some(face) = boundary_face(
                model,
                polygon,
                plane,
                names(self.cells[owner].curtains[edge].then_some(Role::End)),
                &points,
                &vertices,
                &mut edges,
            ) {
                faces.push(face);
            }
        }
        let shell = model.add_shell(Shell { faces });
        let lump = model.add_lump(Lump {
            shells: vec![shell],
        });
        let envelope = model.add_body(Body {
            kind: BodyKind::Solid,
            lumps: vec![lump],
            marks: Default::default(),
        });
        let offset = self.direction.scaled(&self.extent);
        let mut profile = self.profile.clone();
        profile.pose.translation = profile.pose.translation.minus(&offset.to_real());
        let tool = crate::extrude::extrude_vector(
            model,
            &profile,
            step,
            &offset.scaled(&Rational::from_i64(3)),
            [None, None],
        )?;
        Ok(catallax_boolean::boolean(
            model,
            envelope,
            tool,
            catallax_boolean::BoolOp::Intersect,
            dec,
        )
        .map_err(|e| e.refusal)?
        .body)
    }
}
