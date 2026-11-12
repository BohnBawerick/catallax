//! Region growing: local normals, plane-driven growth into small regions,
//! merging of neighbouring regions that one plane, cylinder or cone explains,
//! and a final assignment of every sample to the surface that explains it
//! best. Samples no surface explains stay unassigned; the caller reports them.

use catallax_numbers::{Decider, Refusal};

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use super::fit::{self, Surface, SurfaceKind};
use super::linalg::{V3, dot, eigen, moments, moments_checked, unit};
use super::par_map;

/// Neighbours used for each sample's normal, itself included.
pub(super) const K_NORMAL: usize = 16;
/// Neighbours that connect a sample to the region graph.
const K_GRAPH: usize = 10;
/// Growth admits a sample within this many noise deviations of the surface.
pub(super) const GROWTH_SIGMAS: f64 = 4.0;
/// A merge keeps each side's root-mean-square residual within this many
/// noise deviations.
const MERGE_RMS_SIGMAS: f64 = 1.5;
/// Largest angle between a sample's normal and a growing plane's normal.
const GROWTH_ANGLE_DEGREES: f64 = 30.0;
/// Fewest samples a fitted surface must explain to be kept.
pub(super) const MIN_SURFACE: usize = 30;
/// Regions within this size ratio may merge into a surface of a new kind.
const SIZE_RATIO: usize = 20;

/// Samples with their neighbourhood structure.
pub(super) struct Cloud<'a> {
    pub points: &'a [V3],
    pub normals: Vec<V3>,
    curvature: Vec<f64>,
    centres: Vec<V3>,
    /// Symmetric neighbour lists in ascending order.
    pub graph: Vec<Vec<u32>>,
    /// Estimated noise deviation normal to the surface, in millimetres.
    pub sigma: f64,
}

impl<'a> Cloud<'a> {
    /// Normals, curvature and noise from the `K_NORMAL` nearest samples.
    /// `noise` overrides the estimate when given.
    pub(super) fn new(
        points: &'a [V3],
        noise: Option<f64>,
        workers: usize,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let tree = super::knn::KdTree::new(points, dec)?;
        dec.charge_work(points.len() as u64)?;
        let mut local = Vec::with_capacity(points.len());
        for start in (0..points.len()).step_by(256) {
            let end = (start + 256).min(points.len());
            let near: Vec<_> = (start..end)
                .map(|i| tree.nearest(i, K_NORMAL, dec))
                .collect::<Result<_, _>>()?;
            dec.charge_work(((end - start) * K_NORMAL * 8) as u64)?;
            local.extend(par_map(workers, near.len(), |i| {
                let neighbours = &near[i];
                let (mean, cov, n) = moments(neighbours.iter().map(|&j| points[j as usize]));
                let (values, vectors) = eigen(cov);
                let total = values[0] + values[1] + values[2];
                let normal = unit(vectors[0]).unwrap_or([0.0, 0.0, 1.0]);
                let curvature = if total > 0.0 { values[0] / total } else { 0.0 };
                let variance = if n > 3 {
                    values[0].max(0.0) * n as f64 / (n - 3) as f64
                } else {
                    0.0
                };
                (normal, curvature, mean, variance, neighbours.clone())
            }));
            dec.check_deadline()?;
        }
        dec.charge_work(points.len() as u64)?;
        let mut graph = vec![Vec::new(); points.len()];
        for (i, entry) in local.iter().enumerate() {
            for &j in entry.4.iter().skip(1).take(K_GRAPH) {
                dec.charge_work(1)?;
                graph[i].push(j);
                graph[j as usize].push(i as u32);
            }
        }
        for list in &mut graph {
            super::charge_sort(list.len(), dec)?;
            list.sort_unstable();
            list.dedup();
            dec.check_deadline()?;
        }
        dec.charge_work(local.len() as u64)?;
        let mut variances: Vec<f64> = local.iter().map(|e| e.3).collect();
        super::charge_sort(variances.len(), dec)?;
        variances.sort_by(f64::total_cmp);
        dec.check_deadline()?;
        let estimate = variances.get(variances.len() / 2).map_or(0.0, |v| v.sqrt());
        Ok(Self {
            points,
            normals: local.iter().map(|e| e.0).collect(),
            curvature: local.iter().map(|e| e.1).collect(),
            centres: local.iter().map(|e| e.2).collect(),
            graph,
            sigma: noise.unwrap_or(estimate).max(f64::MIN_POSITIVE),
        })
    }

    fn eps(&self) -> f64 {
        GROWTH_SIGMAS * self.sigma
    }

    /// Whether `surface` explains these samples: root-mean-square residual
    /// within the merge bound and at most one in a hundred beyond growth.
    fn explains(
        &self,
        surface: &Surface,
        members: &[u32],
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let eps = self.eps();
        let mut sum = 0.0;
        let mut far = 0usize;
        for &i in members {
            dec.charge_work(1)?;
            let r = surface.residual(self.points[i as usize]);
            sum += r * r;
            far += usize::from(r.abs() > eps);
        }
        let rms = (sum / members.len().max(1) as f64).sqrt();
        Ok(rms <= MERGE_RMS_SIGMAS * self.sigma && far * 100 <= members.len())
    }

    /// The samples `surface` explains within the growth bound, when they are
    /// at least nine in ten of `members` and their root-mean-square residual
    /// is within the merge bound.
    fn core(
        &self,
        surface: &Surface,
        members: &[u32],
        dec: &mut Decider,
    ) -> Result<Option<Vec<u32>>, Refusal> {
        let eps = self.eps();
        let mut inliers = Vec::new();
        for &i in members {
            dec.charge_work(1)?;
            if surface.residual(self.points[i as usize]).abs() <= eps {
                inliers.push(i);
            }
        }
        Ok(
            (inliers.len() * 10 >= members.len() * 9 && self.explains(surface, &inliers, dec)?)
                .then_some(inliers),
        )
    }

    /// Whether a cone's radius changes visibly along its samples: by more
    /// than the growth bound, or a cylinder would explain them as well.
    fn tapers_enough(
        &self,
        surface: &Surface,
        members: &[u32],
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let Surface::Cone {
            point, axis, slope, ..
        } = *surface
        else {
            return Ok(true);
        };
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for &i in members {
            dec.charge_work(1)?;
            let t = dot(super::linalg::sub(self.points[i as usize], point), axis);
            lo = lo.min(t);
            hi = hi.max(t);
        }
        Ok(slope.abs() * (hi - lo) > self.eps())
    }

    /// Whether a cylinder or cone bends visibly over its samples: its
    /// sagitta across their angular span exceeds the growth bound, so a plane
    /// could not explain them as well. Planes pass.
    fn curved_enough(
        &self,
        surface: &Surface,
        members: &[u32],
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let Some((point, axis)) = surface.axis() else {
            return Ok(true);
        };
        let (u, v) = super::linalg::frame(axis);
        let mut mean = [0.0f64; 2];
        dec.charge_work(members.len() as u64)?;
        let mut directions = Vec::with_capacity(members.len());
        for &i in members {
            dec.charge_work(1)?;
            let q = super::linalg::sub(self.points[i as usize], point);
            let d = [dot(q, u), dot(q, v)];
            let n = d[0].hypot(d[1]).max(f64::MIN_POSITIVE);
            directions.push([d[0] / n, d[1] / n]);
        }
        for d in &directions {
            dec.charge_work(1)?;
            mean[0] += d[0];
            mean[1] += d[1];
        }
        let length = mean[0].hypot(mean[1]);
        if length <= f64::EPSILON * members.len() as f64 {
            return Ok(true);
        }
        let mean = [mean[0] / length, mean[1] / length];
        let mut half_span_cos = 1.0_f64;
        for d in &directions {
            dec.charge_work(1)?;
            half_span_cos = half_span_cos.min(d[0] * mean[0] + d[1] * mean[1]);
        }
        let half_span_cos = half_span_cos.max(-1.0);
        let mut radius = 0.0_f64;
        for &i in members {
            dec.charge_work(1)?;
            let q = super::linalg::sub(self.points[i as usize], point);
            radius = radius.max(surface.radius_at(dot(q, axis)).unwrap_or(0.0).abs());
        }
        Ok(radius * (1.0 - half_span_cos) > self.eps())
    }

    /// The simplest surface that explains each part separately.
    fn fit_union(&self, parts: [&[u32]; 2], dec: &mut Decider) -> Result<Option<Surface>, Refusal> {
        dec.charge_work((parts[0].len() + parts[1].len()) as u64)?;
        let union: Vec<u32> = parts[0].iter().chain(parts[1]).copied().collect();
        for kind in [SurfaceKind::Plane, SurfaceKind::Cylinder, SurfaceKind::Cone] {
            let surface = refit(self, kind, &union, dec)?;
            if let Some(s) = surface
                && self.curved_enough(&s, &union, dec)?
                && self.tapers_enough(&s, &union, dec)?
                && self.explains(&s, parts[0], dec)?
                && self.explains(&s, parts[1], dec)?
            {
                return Ok(Some(s));
            }
        }
        Ok(None)
    }
}

/// Refits a surface of the same kind on its samples.
pub(super) fn refit(
    cloud: &Cloud<'_>,
    kind: SurfaceKind,
    members: &[u32],
    dec: &mut Decider,
) -> Result<Option<Surface>, Refusal> {
    match kind {
        SurfaceKind::Plane => fit::plane(cloud.points, members, dec),
        SurfaceKind::Cylinder => fit::cylinder(cloud.points, &cloud.normals, members, dec),
        SurfaceKind::Cone => fit::cone(cloud.points, &cloud.normals, members, dec),
    }
}

/// Plane-driven growth from the flattest unassigned sample, repeatedly. Every
/// sample ends in exactly one small or large region.
fn grow_planes(cloud: &Cloud<'_>, dec: &mut Decider) -> Result<Vec<Vec<u32>>, Refusal> {
    const NONE: u32 = u32::MAX;
    let n = cloud.points.len();
    dec.charge_work(n as u64)?;
    super::charge_sort(n, dec)?;
    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_by(|&a, &b| {
        cloud.curvature[a as usize]
            .total_cmp(&cloud.curvature[b as usize])
            .then(a.cmp(&b))
    });
    dec.check_deadline()?;
    let cos = GROWTH_ANGLE_DEGREES.to_radians().cos();
    let eps = cloud.eps();
    let mut label = vec![NONE; n];
    let mut regions = Vec::new();
    for &seed in &order {
        dec.charge_work(1)?;
        if label[seed as usize] != NONE {
            continue;
        }
        let id = regions.len() as u32;
        let seed_normal = cloud.normals[seed as usize];
        let mut normal = seed_normal;
        let mut offset = dot(normal, cloud.centres[seed as usize]);
        let mut members = vec![seed];
        label[seed as usize] = id;
        let mut queue = VecDeque::from([seed]);
        let mut next_refit = 16;
        while let Some(i) = queue.pop_front() {
            for &j in &cloud.graph[i as usize] {
                dec.charge_work(1)?;
                let ju = j as usize;
                if label[ju] == NONE
                    && (dot(normal, cloud.points[ju]) - offset).abs() <= eps
                    && dot(normal, cloud.normals[ju]).abs() >= cos
                {
                    label[ju] = id;
                    members.push(j);
                    queue.push_back(j);
                }
            }
            if members.len() >= next_refit {
                next_refit *= 2;
                let (mean, cov, _) =
                    moments_checked(members.iter().map(|&k| cloud.points[k as usize]), dec)?;
                let (values, vectors) = eigen(cov);
                // Refit only a patch, not a line of samples, and never turn
                // away from the seed's own normal.
                if let Some(m) = unit(vectors[0])
                    && values[1] >= 4.0 * values[0].max(0.0)
                    && dot(m, seed_normal).abs() >= cos
                {
                    normal = m;
                    offset = dot(normal, mean);
                }
            }
        }
        regions.push(members);
    }
    Ok(regions)
}

struct Region {
    members: Vec<u32>,
    surface: Option<Surface>,
    /// Bumped when the surface changes kind or the region doubles, so that
    /// failed merges are retried only after a real change.
    version: u32,
    size_at_version: usize,
    alive: bool,
}

/// Regions after merging, and the surfaces kept.
pub(super) struct Segmentation {
    /// Per sample: index into `surfaces`, or `None` when nothing explains it.
    pub assignment: Vec<Option<u32>>,
    /// Kept surfaces with their samples, in a deterministic order.
    pub surfaces: Vec<(Surface, Vec<u32>)>,
}

/// Region growing, merging and final assignment.
pub(super) fn segment(cloud: &Cloud<'_>, dec: &mut Decider) -> Result<Segmentation, Refusal> {
    let n = cloud.points.len();
    let mut regions = Vec::new();
    for members in grow_planes(cloud, dec)? {
        let surface = match fit::plane(cloud.points, &members, dec)? {
            Some(surface) if cloud.explains(&surface, &members, dec)? => Some(surface),
            _ => None,
        };
        regions.push(Region {
            surface,
            size_at_version: members.len(),
            members,
            version: 0,
            alive: true,
        });
    }
    dec.charge_work(n as u64)?;
    let mut owner = vec![0; n];
    for (r, region) in regions.iter().enumerate() {
        for &i in &region.members {
            dec.charge_work(1)?;
            owner[i as usize] = r as u32;
        }
    }
    let mut failed = BTreeSet::new();
    loop {
        let mut changed = false;
        dec.charge_work(regions.len() as u64)?;
        let mut order: Vec<u32> = (0..regions.len() as u32)
            .filter(|&r| regions[r as usize].alive)
            .collect();
        super::charge_sort(order.len(), dec)?;
        order.sort_by(|&a, &b| {
            regions[b as usize]
                .members
                .len()
                .cmp(&regions[a as usize].members.len())
                .then(a.cmp(&b))
        });
        dec.check_deadline()?;
        for &a in &order {
            dec.charge_work(1)?;
            if !regions[a as usize].alive {
                continue;
            }
            loop {
                let mut shared: BTreeMap<u32, usize> = BTreeMap::new();
                for &i in &regions[a as usize].members {
                    for &j in &cloud.graph[i as usize] {
                        dec.charge_work(1)?;
                        let b = owner[j as usize];
                        if b != a {
                            *shared.entry(b).or_default() += 1;
                        }
                    }
                }
                let mut neighbours: Vec<(usize, u32)> =
                    shared.into_iter().map(|(b, c)| (c, b)).collect();
                super::charge_sort(neighbours.len(), dec)?;
                neighbours.sort_by(|x, y| y.0.cmp(&x.0).then(x.1.cmp(&y.1)));
                dec.check_deadline()?;
                let mut merged = false;
                for (_, b) in neighbours {
                    dec.charge_work(1)?;
                    let key = (
                        a,
                        regions[a as usize].version,
                        b,
                        regions[b as usize].version,
                    );
                    if failed.contains(&key) {
                        continue;
                    }
                    let parts = [
                        regions[a as usize].members.as_slice(),
                        regions[b as usize].members.as_slice(),
                    ];
                    let (big, small) = if parts[0].len() >= parts[1].len() {
                        (0, 1)
                    } else {
                        (1, 0)
                    };
                    let comparable = parts[small].len() >= MIN_SURFACE
                        && parts[small].len() * SIZE_RATIO >= parts[big].len();
                    let attempt = if comparable || parts[big].len() < MIN_SURFACE {
                        cloud.fit_union(parts, dec)?
                    } else {
                        let r = if big == 0 { a } else { b };
                        match regions[r as usize].surface.as_ref() {
                            Some(surface) if cloud.explains(surface, parts[small], dec)? => {
                                dec.charge_work((parts[0].len() + parts[1].len()) as u64)?;
                                let union =
                                    parts[0].iter().chain(parts[1]).copied().collect::<Vec<_>>();
                                Some(
                                    refit(cloud, surface.kind(), &union, dec)?
                                        .unwrap_or_else(|| surface.clone()),
                                )
                            }
                            _ => None,
                        }
                    };
                    match attempt {
                        Some(surface) => {
                            let moved = std::mem::take(&mut regions[b as usize].members);
                            for &i in &moved {
                                dec.charge_work(1)?;
                                owner[i as usize] = a;
                            }
                            regions[b as usize].alive = false;
                            let region = &mut regions[a as usize];
                            region.members.extend(moved);
                            let kind_changed =
                                region.surface.as_ref().map(Surface::kind) != Some(surface.kind());
                            region.surface = Some(surface);
                            if kind_changed || region.members.len() >= 2 * region.size_at_version {
                                region.version += 1;
                                region.size_at_version = region.members.len();
                            }
                            merged = true;
                            changed = true;
                            break;
                        }
                        None => {
                            failed.insert(key);
                        }
                    }
                }
                if !merged {
                    break;
                }
            }
        }
        if !changed {
            break;
        }
    }
    let mut kept: Vec<_> = regions
        .into_iter()
        .filter(|r| r.alive && r.members.len() >= MIN_SURFACE)
        .filter_map(|r| r.surface.map(|s| (s, r.members)))
        .collect();
    let mut found = assign(cloud, &mut kept, dec)?;
    while merge_kept(cloud, &mut found.surfaces, dec)? {
        found = assign(cloud, &mut found.surfaces, dec)?;
    }
    let mut recovered = false;
    for group in unexplained(cloud, &found.assignment, dec)? {
        if group.len() >= MIN_SURFACE
            && let Some(surface) = cloud.fit_union([&group, &[]], dec)?
        {
            found.surfaces.push((surface, group));
            recovered = true;
        }
    }
    if recovered {
        found = assign(cloud, &mut found.surfaces, dec)?;
        while merge_kept(cloud, &mut found.surfaces, dec)? {
            found = assign(cloud, &mut found.surfaces, dec)?;
        }
    }
    Ok(found)
}

/// Joins kept surfaces that one surface explains, connected or not: a
/// smaller one onto a larger one's surface of any kind, or two of one kind
/// onto a common fit. Returns whether anything joined.
fn merge_kept(
    cloud: &Cloud<'_>,
    kept: &mut Vec<(Surface, Vec<u32>)>,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    super::charge_sort(kept.len(), dec)?;
    let mut order: Vec<usize> = (0..kept.len()).collect();
    order.sort_by(|&a, &b| kept[b].1.len().cmp(&kept[a].1.len()).then(a.cmp(&b)));
    dec.check_deadline()?;
    for (x, &a) in order.iter().enumerate() {
        for &b in &order[x + 1..] {
            dec.charge_work(1)?;
            let (big, small) = (&kept[a], &kept[b]);
            let joined = if cloud.explains(&big.0, &small.1, dec)? {
                dec.charge_work((big.1.len() + small.1.len()) as u64)?;
                let union: Vec<u32> = big.1.iter().chain(&small.1).copied().collect();
                match refit(cloud, big.0.kind(), &union, dec)? {
                    Some(surface)
                        if cloud.explains(&surface, &big.1, dec)?
                            && cloud.explains(&surface, &small.1, dec)? =>
                    {
                        Some(surface)
                    }
                    _ => None,
                }
            } else if small.0.kind() == big.0.kind() && small.1.len() * SIZE_RATIO >= big.1.len() {
                cloud.fit_union([&big.1, &small.1], dec)?
            } else {
                None
            };
            let (joined, moved) = match joined {
                Some(s) => (Some(s), std::mem::take(&mut kept[b].1)),
                None => match cloud.core(&kept[a].0, &kept[b].1, dec)? {
                    Some(inliers) => {
                        dec.charge_work((kept[a].1.len() + inliers.len()) as u64)?;
                        let union = kept[a]
                            .1
                            .iter()
                            .chain(&inliers)
                            .copied()
                            .collect::<Vec<_>>();
                        let surface = refit(cloud, kept[a].0.kind(), &union, dec)?
                            .unwrap_or_else(|| kept[a].0.clone());
                        (Some(surface), inliers)
                    }
                    None => (None, Vec::new()),
                },
            };
            if let Some(surface) = joined {
                kept[a].1.extend(moved);
                kept[a].0 = surface;
                kept.remove(b);
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// Grows kept surfaces over unassigned neighbours within the growth bound,
/// moves boundary samples to the surface with the smallest residual, refits,
/// and drops surfaces left with too few samples. Every round updates all
/// samples from the previous round's state, so the result does not depend on
/// the worker count.
pub(super) fn assign(
    cloud: &Cloud<'_>,
    kept: &mut [(Surface, Vec<u32>)],
    dec: &mut Decider,
) -> Result<Segmentation, Refusal> {
    let n = cloud.points.len();
    let eps = cloud.eps();
    dec.charge_work(n as u64)?;
    let mut assignment = vec![None; n];
    for (s, (_, members)) in kept.iter().enumerate() {
        for &i in members {
            dec.charge_work(1)?;
            assignment[i as usize] = Some(s as u32);
        }
    }
    for _ in 0..3 {
        for round in 0..64 {
            let before = &assignment;
            dec.charge_work(n as u64)?;
            let mut next = Vec::with_capacity(n);
            for i in 0..n {
                dec.charge_work(1)?;
                if round > 0 && before[i].is_some() {
                    next.push(before[i]);
                    continue;
                }
                let mut best =
                    before[i].map(|s| (kept[s as usize].0.residual(cloud.points[i]).abs(), s));
                for &j in &cloud.graph[i] {
                    dec.charge_work(1)?;
                    if let Some(s) = before[j as usize] {
                        let r = kept[s as usize].0.residual(cloud.points[i]).abs();
                        if r <= eps
                            && best.is_none_or(|(b, t)| r.total_cmp(&b).then(s.cmp(&t)).is_lt())
                        {
                            best = Some((r, s));
                        }
                    }
                }
                next.push(best.map(|(_, s)| s));
            }
            dec.charge_work(n as u64)?;
            let changed = next != assignment;
            assignment = next;
            if !changed {
                break;
            }
        }
        dec.charge_work(kept.len() as u64)?;
        let mut members = vec![Vec::new(); kept.len()];
        for (i, s) in assignment.iter().enumerate() {
            dec.charge_work(1)?;
            if let Some(s) = s {
                members[*s as usize].push(i as u32);
            }
        }
        for (k, list) in kept.iter_mut().zip(members) {
            let refitted = match refit(cloud, k.0.kind(), &list, dec)? {
                Some(s) if cloud.tapers_enough(&s, &list, dec)? => Some(s),
                Some(_) => refit(cloud, SurfaceKind::Cylinder, &list, dec)?,
                None => None,
            };
            if let Some(s) = refitted {
                k.0 = s;
            }
            k.1 = list;
        }
        for (i, s) in assignment.iter_mut().enumerate() {
            dec.charge_work(1)?;
            if let Some(t) = *s
                && kept[t as usize].0.residual(cloud.points[i]).abs() > eps
            {
                *s = None;
            }
        }
    }
    dec.charge_work(kept.len() as u64)?;
    let mut renumber = vec![None; kept.len()];
    let mut surfaces = Vec::new();
    for (s, (surface, _)) in kept.iter().enumerate() {
        let mut list = Vec::new();
        for i in 0..n as u32 {
            dec.charge_work(1)?;
            if assignment[i as usize] == Some(s as u32) {
                list.push(i);
            }
        }
        if list.len() >= MIN_SURFACE {
            renumber[s] = Some(surfaces.len() as u32);
            surfaces.push((surface.clone(), list));
        }
    }
    dec.charge_work(n as u64)?;
    let assignment = assignment
        .into_iter()
        .map(|s| s.and_then(|s| renumber[s as usize]))
        .collect();
    Ok(Segmentation {
        assignment,
        surfaces,
    })
}

/// Connected groups of unassigned samples through the neighbour graph, each
/// in ascending sample order, the groups ordered by their first sample.
pub(super) fn unexplained(
    cloud: &Cloud<'_>,
    assignment: &[Option<u32>],
    dec: &mut Decider,
) -> Result<Vec<Vec<u32>>, Refusal> {
    dec.charge_work(assignment.len() as u64)?;
    let mut seen = vec![false; assignment.len()];
    let mut groups = Vec::new();
    for start in 0..assignment.len() {
        dec.charge_work(1)?;
        if seen[start] || assignment[start].is_some() {
            continue;
        }
        seen[start] = true;
        let mut group = vec![start as u32];
        let mut queue = vec![start as u32];
        while let Some(i) = queue.pop() {
            for &j in &cloud.graph[i as usize] {
                dec.charge_work(1)?;
                if !seen[j as usize] && assignment[j as usize].is_none() {
                    seen[j as usize] = true;
                    group.push(j);
                    queue.push(j);
                }
            }
        }
        super::charge_sort(group.len(), dec)?;
        group.sort_unstable();
        dec.check_deadline()?;
        groups.push(group);
    }
    Ok(groups)
}
