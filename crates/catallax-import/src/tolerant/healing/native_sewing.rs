//! Native sewing uses measured finite samples, never an exact-image claim.
use std::sync::Arc;

use crate::geometry_source::{AngleUnit, ParsedStep, Topology};
use super::super::native::{NativeCurve, NativeSample};
use super::{sewing, *};

const STATIONS: usize = 16;
type Samples = [NativeSample; STATIONS + 1];
type Curves = BTreeMap<u64, Arc<NativeCurve>>;
fn norm(point: &[f64; 3]) -> f64 { point.iter().fold(0.0, |result, value| result.hypot(*value)) }
fn distance(a: &[f64; 3], b: &[f64; 3]) -> f64 { (0..3).fold(0.0, |result, axis| result.hypot(a[axis] - b[axis])) }
fn point(value: &Vec3) -> Result<[f64; 3], Error> {
    let point = value.each_ref().map(Rational::to_f64_nearest);
    if point.iter().all(|coordinate| coordinate.is_finite()) { Ok(point) } else { Err(failure("import.native-sampling-range", "source coordinate exceeds binary64 sampling range")) }
}

fn sample_edge(parsed: &ParsedStep, scale: &Rational, angles: &AngleUnit, curves: &mut Curves, body: &TolerantBody, edge_index: usize, dec: &mut Decider) -> Result<Samples, Error> {
    let edge = &body.edges[edge_index];
    let endpoints = [point(&body.vertices[edge.vertices[0]].source.point)?, point(&body.vertices[edge.vertices[1]].source.point)?];
    let source_curve = match edge.name.strip_prefix('#').and_then(|name| name.parse().ok()) {
        Some(entity) => match parsed.topology(entity, dec)? { Topology::Edge { curve, .. } => curve, _ => None }, None => None,
    };
    if let Some(entity) = source_curve {
        let curve = match curves.entry(entity) {
            std::collections::btree_map::Entry::Occupied(entry) => entry.into_mut(),
            std::collections::btree_map::Entry::Vacant(entry) => entry.insert(Arc::new(NativeCurve::from_source(parsed, entity, scale, angles.clone(), dec)?)),
        };
        let a = curve.project(endpoints[0], None, dec)?;
        let b = curve.project(endpoints[1], None, dec)?;
        let range = curve.finite_range([a.location.parameters[0], b.location.parameters[0]], edge.curve_same_sense, edge.vertices[0] == edge.vertices[1])?;
        let evaluate = |station: usize, dec: &mut Decider| {
            let parameter = range[0] + (range[1] - range[0]) * station as f64 / STATIONS as f64;
            curve.sample(curve.location(parameter)?, dec)
        };
        let first = evaluate(0, dec)?;
        let mut samples = [first; STATIONS + 1];
        for (station, value) in samples.iter_mut().enumerate().skip(1) { *value = evaluate(station, dec)?; }
        return Ok(samples);
    }
    let Curve::Line(origin, direction) = &edge.curve else {
        return Err(failure("import.tolerant-broken-reference", "free boundary has no source curve"));
    };
    let origin = point(origin)?; let direction = point(direction)?;
    let length: f64 = direction.iter().map(|value| value * value).sum();
    if length <= 0.0 || !length.is_finite() { return Err(failure("import.native-sampling-range", "source line has no finite sampled direction")); }
    let parameters = endpoints.map(|point| (0..3).map(|axis| (point[axis] - origin[axis]) * direction[axis]).sum::<f64>() / length);
    let sample = |station: usize| {
        let parameter = parameters[0] + station as f64 / STATIONS as f64 * (parameters[1] - parameters[0]);
        let position = std::array::from_fn(|axis| origin[axis] + parameter * direction[axis]);
        NativeSample { point_mm: position, uncertainty_mm: 16.0 * f64::EPSILON * norm(&position) }
    };
    let first = sample(0);
    let mut samples = [first; STATIONS + 1];
    for (station, value) in samples.iter_mut().enumerate().skip(1) { dec.charge_work(1)?; *value = sample(station); }
    Ok(samples)
}

pub(super) fn sew(body_index: usize, body: &mut TolerantBody, parsed: &ParsedStep, scale: &Rational, angles: &AngleUnit, curves: &mut Curves, plan: &mut Plan, dec: &mut Decider) -> Result<(), Error> {
    let mut samples = BTreeMap::<usize, Samples>::new();
    let join_shells = !body.declared_solid;
    sewing::sew_with(body_index, body, plan, dec, true, join_shells, |body, a, b, dec| {
        let cap_squared = arithmetic(&body.options.maximum_mm, &body.options.maximum_mm, '*', dec)?;
        let left = &body.edges[a]; let right = &body.edges[b];
        let mut result = None;
        for reversed in [false, true] {
            let mut gaps = [Rational::zero(), Rational::zero()];
            let mut endpoints_fit = true;
            for endpoint in 0..2 {
                let other = if reversed { 1 - endpoint } else { endpoint };
                gaps[endpoint] = distance_squared(&body.vertices[left.vertices[endpoint]].source.point, &body.vertices[right.vertices[other]].source.point, dec)?;
                if compare(&gaps[endpoint], &cap_squared, dec)?.is_gt() { endpoints_fit = false; break; }
            }
            if !endpoints_fit { continue; }
            for edge in [a, b] {
                if let std::collections::btree_map::Entry::Vacant(entry) = samples.entry(edge) {
                    entry.insert(sample_edge(parsed, scale, angles, curves, body, edge, dec)?);
                }
            }
            let first = &samples[&a]; let second = &samples[&b];
            let mut maximum: f64 = 0.0;
            for station in 0..=STATIONS {
                dec.charge_work(1)?;
                let other = if reversed { STATIONS - station } else { station };
                let residual = distance(&first[station].point_mm, &second[other].point_mm);
                maximum = maximum.max(residual + first[station].uncertainty_mm + second[other].uncertainty_mm + 16.0 * f64::EPSILON * residual.max(1.0));
            }
            let maximum = Rational::from_f64_exact(maximum).ok_or_else(|| failure("import.native-sampling-range", "nonfinite native sew residual"))?;
            if compare(&maximum, &body.options.maximum_mm, dec)?.is_gt() { continue; }
            if result.is_some() { return Ok(None); }
            result = Some(sewing::Match { other: b, reversed, gaps, check: "17-finite-native-curve-sample-correspondence", maximum_sample_mm: Some(maximum) });
        }
        Ok(result)
    })?;
    if join_shells { merge_shells(body, plan, dec)?; }
    Ok(())
}

fn merge_shells(body: &mut TolerantBody, plan: &mut Plan, dec: &mut Decider) -> Result<(), Error> {
    fn root(parents: &[usize], mut value: usize) -> usize { while parents[value] != value { value = parents[value]; } value }
    let uses = usage(body, dec)?;
    let mut parents: Vec<_> = (0..body.shells.len()).collect();
    for edge in &uses {
        if let [a, b] = edge.as_slice() {
            let a = root(&parents, a.shell); let b = root(&parents, b.shell);
            if a != b { parents[a.max(b)] = a.min(b); }
        }
    }
    if parents.iter().enumerate().all(|(index, &parent)| index == parent) { return Ok(()); }
    let mut combined = BTreeMap::<usize, (super::super::TolerantShell, Vec<String>)>::new();
    for (index, shell) in std::mem::take(&mut body.shells).into_iter().enumerate() {
        dec.charge_work(shell.faces.len() as u64 + 1)?;
        let parent = root(&parents, index);
        match combined.entry(parent) {
            std::collections::btree_map::Entry::Vacant(entry) => { let name = shell.name.clone(); entry.insert((shell, vec![name])); }
            std::collections::btree_map::Entry::Occupied(mut entry) => { entry.get_mut().1.push(shell.name); entry.get_mut().0.faces.extend(shell.faces); }
        }
    }
    for (_, (shell, names)) in combined {
        if names.len() > 1 { plan.changes.push(HealingChange::ShellsSewn { body: body.name.clone(), source_shells: names, retained_shell: shell.name.clone() }); }
        body.shells.push(shell);
    }
    Ok(())
}
