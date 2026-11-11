//! Region growing: local normals, plane-driven growth into small regions,
//! merging of neighbouring regions that one plane, cylinder or cone explains,
//! and a final assignment of every sample to the surface that explains it
//! best. Samples no surface explains stay unassigned; the caller reports them.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use super::fit::{self, Surface, SurfaceKind};
use super::linalg::{V3, dot, eigen, moments, unit};
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
    pub(super) fn new(points: &'a [V3], noise: Option<f64>, workers: usize) -> Self {
        let tree = super::knn::KdTree::new(points);
        let local = par_map(workers, points.len(), |i| {
            let near = tree.nearest(i, K_NORMAL);
            let (mean, cov, n) = moments(near.iter().map(|&j| points[j as usize]));
            let (values, vectors) = eigen(cov);
            let total = values[0] + values[1] + values[2];
            let normal = unit(vectors[0]).unwrap_or([0.0, 0.0, 1.0]);
            let curvature = if total > 0.0 { values[0] / total } else { 0.0 };
            let variance = if n > 3 {
                values[0].max(0.0) * n as f64 / (n - 3) as f64
            } else {
                0.0
            };
            (normal, curvature, mean, variance, near)
        });
        let mut graph = vec![Vec::new(); points.len()];
        for (i, entry) in local.iter().enumerate() {
            for &j in entry.4.iter().skip(1).take(K_GRAPH) {
                graph[i].push(j);
                graph[j as usize].push(i as u32);
            }
        }
        for list in &mut graph {
            list.sort_unstable();
            list.dedup();
        }
        let mut variances: Vec<f64> = local.iter().map(|e| e.3).collect();
        variances.sort_by(f64::total_cmp);
        let estimate = variances.get(variances.len() / 2).map_or(0.0, |v| v.sqrt());
        Self {
            points,
            normals: local.iter().map(|e| e.0).collect(),
            curvature: local.iter().map(|e| e.1).collect(),
            centres: local.iter().map(|e| e.2).collect(),
            graph,
            sigma: noise.unwrap_or(estimate).max(f64::MIN_POSITIVE),
        }
    }

    fn eps(&self) -> f64 {
        GROWTH_SIGMAS * self.sigma
    }

    /// Whether `surface` explains these samples: root-mean-square residual
    /// within the merge bound and at most one in a hundred beyond growth.
    fn explains(&self, surface: &Surface, members: &[u32]) -> bool {
        let eps = self.eps();
        let mut sum = 0.0;
        let mut far = 0usize;
        for &i in members {
            let r = surface.residual(self.points[i as usize]);
            sum += r * r;
            far += usize::from(r.abs() > eps);
        }
        let rms = (sum / members.len().max(1) as f64).sqrt();
        rms <= MERGE_RMS_SIGMAS * self.sigma && far * 100 <= members.len()
    }

    /// The samples `surface` explains within the growth bound, when they are
    /// at least nine in ten of `members` and their root-mean-square residual
    /// is within the merge bound.
    fn core(&self, surface: &Surface, members: &[u32]) -> Option<Vec<u32>> {
        let eps = self.eps();
        let inliers: Vec<u32> = members
            .iter()
            .copied()
            .filter(|&i| surface.residual(self.points[i as usize]).abs() <= eps)
            .collect();
        (inliers.len() * 10 >= members.len() * 9 && self.explains(surface, &inliers))
            .then_some(inliers)
    }

    /// Whether a cone's radius changes visibly along its samples: by more
    /// than the growth bound, or a cylinder would explain them as well.
    fn tapers_enough(&self, surface: &Surface, members: &[u32]) -> bool {
        let Surface::Cone {
            point, axis, slope, ..
        } = *surface
        else {
            return true;
        };
        let (lo, hi) = members
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &i| {
                let t = dot(super::linalg::sub(self.points[i as usize], point), axis);
                (lo.min(t), hi.max(t))
            });
        slope.abs() * (hi - lo) > self.eps()
    }

    /// Whether a cylinder or cone bends visibly over its samples: its
    /// sagitta across their angular span exceeds the growth bound, so a plane
    /// could not explain them as well. Planes pass.
    fn curved_enough(&self, surface: &Surface, members: &[u32]) -> bool {
        let Some((point, axis)) = surface.axis() else {
            return true;
        };
        let (u, v) = super::linalg::frame(axis);
        let mut mean = [0.0f64; 2];
        let directions: Vec<[f64; 2]> = members
            .iter()
            .map(|&i| {
                let q = super::linalg::sub(self.points[i as usize], point);
                let d = [dot(q, u), dot(q, v)];
                let n = d[0].hypot(d[1]).max(f64::MIN_POSITIVE);
                [d[0] / n, d[1] / n]
            })
            .collect();
        for d in &directions {
            mean[0] += d[0];
            mean[1] += d[1];
        }
        let length = mean[0].hypot(mean[1]);
        if length <= f64::EPSILON * members.len() as f64 {
            return true;
        }
        let mean = [mean[0] / length, mean[1] / length];
        let half_span_cos = directions
            .iter()
            .map(|d| d[0] * mean[0] + d[1] * mean[1])
            .fold(1.0f64, f64::min)
            .max(-1.0);
        let radius = members
            .iter()
            .map(|&i| {
                let q = super::linalg::sub(self.points[i as usize], point);
                surface.radius_at(dot(q, axis)).unwrap_or(0.0).abs()
            })
            .fold(0.0f64, f64::max);
        radius * (1.0 - half_span_cos) > self.eps()
    }

    /// The simplest surface that explains each part separately.
    fn fit_union(&self, parts: [&[u32]; 2]) -> Option<Surface> {
        let union: Vec<u32> = parts[0].iter().chain(parts[1]).copied().collect();
        for kind in [SurfaceKind::Plane, SurfaceKind::Cylinder, SurfaceKind::Cone] {
            let surface = match kind {
                SurfaceKind::Plane => fit::plane(self.points, &union),
                SurfaceKind::Cylinder => fit::cylinder(self.points, &self.normals, &union),
                SurfaceKind::Cone => fit::cone(self.points, &self.normals, &union),
            };
            if let Some(s) = surface
                && self.curved_enough(&s, &union)
                && self.tapers_enough(&s, &union)
                && self.explains(&s, parts[0])
                && self.explains(&s, parts[1])
            {
                return Some(s);
            }
        }
        None
    }
}

/// Refits a surface of the same kind on its samples.
pub(super) fn refit(cloud: &Cloud<'_>, kind: SurfaceKind, members: &[u32]) -> Option<Surface> {
    match kind {
        SurfaceKind::Plane => fit::plane(cloud.points, members),
        SurfaceKind::Cylinder => fit::cylinder(cloud.points, &cloud.normals, members),
        SurfaceKind::Cone => fit::cone(cloud.points, &cloud.normals, members),
    }
}

/// Plane-driven growth from the flattest unassigned sample, repeatedly. Every
/// sample ends in exactly one small or large region.
fn grow_planes(cloud: &Cloud<'_>, work: &mut u64) -> Vec<Vec<u32>> {
    const NONE: u32 = u32::MAX;
    let n = cloud.points.len();
    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_by(|&a, &b| {
        cloud.curvature[a as usize]
            .total_cmp(&cloud.curvature[b as usize])
            .then(a.cmp(&b))
    });
    let cos = GROWTH_ANGLE_DEGREES.to_radians().cos();
    let eps = cloud.eps();
    let mut label = vec![NONE; n];
    let mut regions = Vec::new();
    for &seed in &order {
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
                *work += 1;
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
                let (mean, cov, _) = moments(members.iter().map(|&k| cloud.points[k as usize]));
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
    regions
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
    /// Deterministic work units spent.
    pub work: u64,
}

/// Region growing, merging and final assignment.
pub(super) fn segment(cloud: &Cloud<'_>, workers: usize) -> Segmentation {
    let mut work = 0u64;
    let n = cloud.points.len();
    let mut regions: Vec<Region> = grow_planes(cloud, &mut work)
        .into_iter()
        .map(|members| Region {
            surface: fit::plane(cloud.points, &members).filter(|s| cloud.explains(s, &members)),
            size_at_version: members.len(),
            members,
            version: 0,
            alive: true,
        })
        .collect();
    let mut owner: Vec<u32> = vec![0; n];
    for (r, region) in regions.iter().enumerate() {
        for &i in &region.members {
            owner[i as usize] = r as u32;
        }
    }
    // Greedy merging, largest regions first, each with its neighbours in
    // order of shared graph edges. Failed pairs are not retried until one
    // side changes.
    let mut failed: BTreeSet<(u32, u32, u32, u32)> = BTreeSet::new();
    loop {
        let mut changed = false;
        let mut order: Vec<u32> = (0..regions.len() as u32)
            .filter(|&r| regions[r as usize].alive)
            .collect();
        order.sort_by(|&a, &b| {
            regions[b as usize]
                .members
                .len()
                .cmp(&regions[a as usize].members.len())
                .then(a.cmp(&b))
        });
        for &a in &order {
            if !regions[a as usize].alive {
                continue;
            }
            loop {
                let mut shared: BTreeMap<u32, usize> = BTreeMap::new();
                for &i in &regions[a as usize].members {
                    for &j in &cloud.graph[i as usize] {
                        let b = owner[j as usize];
                        if b != a {
                            *shared.entry(b).or_default() += 1;
                        }
                    }
                }
                let mut neighbours: Vec<(usize, u32)> =
                    shared.into_iter().map(|(b, c)| (c, b)).collect();
                neighbours.sort_by(|x, y| y.0.cmp(&x.0).then(x.1.cmp(&y.1)));
                let mut merged = false;
                for (_, b) in neighbours {
                    let key = (
                        a,
                        regions[a as usize].version,
                        b,
                        regions[b as usize].version,
                    );
                    if failed.contains(&key) {
                        continue;
                    }
                    work += (regions[a as usize].members.len() + regions[b as usize].members.len())
                        as u64;
                    let parts = [
                        regions[a as usize].members.as_slice(),
                        regions[b as usize].members.as_slice(),
                    ];
                    // A region much smaller than its neighbour joins it only on
                    // the neighbour's own surface. Regions of comparable size,
                    // or two small ones, may agree on any kind.
                    let (big, small) = if parts[0].len() >= parts[1].len() {
                        (0, 1)
                    } else {
                        (1, 0)
                    };
                    let comparable = parts[small].len() >= MIN_SURFACE
                        && parts[small].len() * SIZE_RATIO >= parts[big].len();
                    let attempt = if comparable || parts[big].len() < MIN_SURFACE {
                        cloud.fit_union(parts)
                    } else {
                        let owner_region = if big == 0 { a } else { b };
                        regions[owner_region as usize]
                            .surface
                            .as_ref()
                            .filter(|s| cloud.explains(s, parts[small]))
                            .map(|s| {
                                let union: Vec<u32> =
                                    parts[0].iter().chain(parts[1]).copied().collect();
                                refit(cloud, s.kind(), &union).unwrap_or_else(|| s.clone())
                            })
                    };
                    match attempt {
                        Some(surface) => {
                            let moved = std::mem::take(&mut regions[b as usize].members);
                            for &i in &moved {
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
    let mut kept: Vec<(Surface, Vec<u32>)> = regions
        .into_iter()
        .filter(|r| r.alive && r.members.len() >= MIN_SURFACE)
        .filter_map(|r| r.surface.map(|s| (s, r.members)))
        .collect();
    let mut found = assign(cloud, &mut kept, workers, &mut work);
    while merge_kept(cloud, &mut found.surfaces, &mut work) {
        found = assign(cloud, &mut found.surfaces, workers, &mut work);
    }
    // One recovery pass: a connected group of unexplained samples that one
    // surface explains as a whole, missed because its growth strips were too
    // narrow to merge, becomes a surface.
    let mut recovered = false;
    for group in unexplained(cloud, &found.assignment) {
        if group.len() >= MIN_SURFACE
            && let Some(surface) = cloud.fit_union([&group, &[]])
        {
            work += group.len() as u64;
            found.surfaces.push((surface, group));
            recovered = true;
        }
    }
    if recovered {
        found = assign(cloud, &mut found.surfaces, workers, &mut work);
        while merge_kept(cloud, &mut found.surfaces, &mut work) {
            found = assign(cloud, &mut found.surfaces, workers, &mut work);
        }
    }
    found.work = work;
    found
}

/// Joins kept surfaces that one surface explains, connected or not: a
/// smaller one onto a larger one's surface of any kind, or two of one kind
/// onto a common fit. Returns whether anything joined.
fn merge_kept(cloud: &Cloud<'_>, kept: &mut Vec<(Surface, Vec<u32>)>, work: &mut u64) -> bool {
    let mut order: Vec<usize> = (0..kept.len()).collect();
    order.sort_by(|&a, &b| kept[b].1.len().cmp(&kept[a].1.len()).then(a.cmp(&b)));
    for (x, &a) in order.iter().enumerate() {
        for &b in &order[x + 1..] {
            *work += (kept[a].1.len() + kept[b].1.len()) as u64;
            let (big, small) = (&kept[a], &kept[b]);
            let joined = if cloud.explains(&big.0, &small.1) {
                let union: Vec<u32> = big.1.iter().chain(&small.1).copied().collect();
                refit(cloud, big.0.kind(), &union)
                    .filter(|s| cloud.explains(s, &big.1) && cloud.explains(s, &small.1))
            } else if small.0.kind() == big.0.kind() && small.1.len() * SIZE_RATIO >= big.1.len() {
                cloud.fit_union([&big.1, &small.1])
            } else {
                None
            };
            // A smaller surface whose samples the larger one mostly explains
            // gives those samples up; the rest are assigned again.
            let (joined, moved) = match joined {
                Some(s) => (Some(s), std::mem::take(&mut kept[b].1)),
                None => match cloud.core(&kept[a].0, &kept[b].1) {
                    Some(inliers) => {
                        let union: Vec<u32> = kept[a].1.iter().chain(&inliers).copied().collect();
                        let surface = refit(cloud, kept[a].0.kind(), &union)
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
                return true;
            }
        }
    }
    false
}

/// Grows kept surfaces over unassigned neighbours within the growth bound,
/// moves boundary samples to the surface with the smallest residual, refits,
/// and drops surfaces left with too few samples. Every round updates all
/// samples from the previous round's state, so the result does not depend on
/// the worker count.
pub(super) fn assign(
    cloud: &Cloud<'_>,
    kept: &mut [(Surface, Vec<u32>)],
    workers: usize,
    work: &mut u64,
) -> Segmentation {
    let n = cloud.points.len();
    let eps = cloud.eps();
    let mut assignment: Vec<Option<u32>> = vec![None; n];
    for (s, (_, members)) in kept.iter().enumerate() {
        for &i in members {
            assignment[i as usize] = Some(s as u32);
        }
    }
    for _ in 0..3 {
        // Expansion over unassigned samples, then one boundary exchange.
        for round in 0..64 {
            let surfaces: Vec<&Surface> = kept.iter().map(|k| &k.0).collect();
            let before = &assignment;
            let next = par_map(workers, n, |i| {
                if round > 0 && before[i].is_some() {
                    return before[i];
                }
                let mut best: Option<(f64, u32)> =
                    before[i].map(|s| (surfaces[s as usize].residual(cloud.points[i]).abs(), s));
                for &j in &cloud.graph[i] {
                    if let Some(s) = before[j as usize] {
                        let r = surfaces[s as usize].residual(cloud.points[i]).abs();
                        if r <= eps
                            && best.is_none_or(|(b, t)| r.total_cmp(&b).then(s.cmp(&t)).is_lt())
                        {
                            best = Some((r, s));
                        }
                    }
                }
                best.map(|(_, s)| s)
            });
            *work += (n * K_GRAPH) as u64;
            let changed = next != assignment;
            assignment = next;
            if !changed {
                break;
            }
        }
        let mut members: Vec<Vec<u32>> = vec![Vec::new(); kept.len()];
        for (i, s) in assignment.iter().enumerate() {
            if let Some(s) = s {
                members[*s as usize].push(i as u32);
            }
        }
        for (k, list) in kept.iter_mut().zip(members) {
            *work += list.len() as u64;
            // A cone that no longer tapers visibly is a cylinder.
            let refitted = refit(cloud, k.0.kind(), &list).and_then(|s| {
                if cloud.tapers_enough(&s, &list) {
                    Some(s)
                } else {
                    refit(cloud, SurfaceKind::Cylinder, &list)
                }
            });
            if let Some(s) = refitted {
                k.0 = s;
            }
            k.1 = list;
        }
        // A sample beyond the growth bound of its refitted surface is released.
        for (i, s) in assignment.iter_mut().enumerate() {
            if let Some(t) = *s
                && kept[t as usize].0.residual(cloud.points[i]).abs() > eps
            {
                *s = None;
            }
        }
    }
    // Drop surfaces left too small, renumbering the rest in order.
    let mut renumber = vec![None; kept.len()];
    let mut surfaces = Vec::new();
    for (s, (surface, _)) in kept.iter().enumerate() {
        let list: Vec<u32> = (0..n as u32)
            .filter(|&i| assignment[i as usize] == Some(s as u32))
            .collect();
        if list.len() >= MIN_SURFACE {
            renumber[s] = Some(surfaces.len() as u32);
            surfaces.push((surface.clone(), list));
        }
    }
    let assignment = assignment
        .into_iter()
        .map(|s| s.and_then(|s| renumber[s as usize]))
        .collect();
    Segmentation {
        assignment,
        surfaces,
        work: *work,
    }
}

/// Connected groups of unassigned samples through the neighbour graph, each
/// in ascending sample order, the groups ordered by their first sample.
pub(super) fn unexplained(cloud: &Cloud<'_>, assignment: &[Option<u32>]) -> Vec<Vec<u32>> {
    let mut seen = vec![false; assignment.len()];
    let mut groups = Vec::new();
    for start in 0..assignment.len() {
        if seen[start] || assignment[start].is_some() {
            continue;
        }
        seen[start] = true;
        let mut group = vec![start as u32];
        let mut queue = vec![start as u32];
        while let Some(i) = queue.pop() {
            for &j in &cloud.graph[i as usize] {
                if !seen[j as usize] && assignment[j as usize].is_none() {
                    seen[j as usize] = true;
                    group.push(j);
                    queue.push(j);
                }
            }
        }
        group.sort_unstable();
        groups.push(group);
    }
    groups
}
