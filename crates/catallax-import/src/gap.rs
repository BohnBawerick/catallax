//! Gap measurement on recognised data, with the repair policy table.
//!
//! A gap is a set of straight boundary edges, used once and never matched by an
//! opposite use, that pair up with each other within the search radius. Distances
//! are exact squares of rational coordinates; the decimal bounds shown are a
//! certified enclosure of their square roots.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Integer, Rational};

use crate::exact::{cross, dot, isqrt, sub};
use crate::{
    Curve, Disposition, Face, GapKind, ImportRefusal, OrientedEdge, Proposal, ProposalSet, Solid,
    StepFile, Surface, Vec3, sha256_hex,
};

/// Largest gap repaired without asking, in millimetres.
pub(crate) fn auto_cap() -> Rational {
    Rational::parse(crate::AUTO_CAP_MM).expect("decimal cap")
}

/// Options of [`propose`].
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct ImportOptions {
    /// Farthest pair of boundary edges that is still reported as a gap, in millimetres.
    /// Farther unmatched boundary edges are an open boundary, which is refused.
    pub search_radius_mm: Rational,
}

impl Default for ImportOptions {
    fn default() -> Self {
        Self {
            search_radius_mm: Rational::one(),
        }
    }
}

fn not_solid(detail: impl Into<String>) -> ImportRefusal {
    ImportRefusal::NotSolid {
        detail: detail.into(),
    }
}

/// One use of a straight edge in a bound.
#[derive(Clone, Debug)]
pub(crate) struct Use {
    pub face: usize,
    pub edge: String,
    pub from: Vec3,
    pub to: Vec3,
}

/// Two opposite free uses that would be one seam.
#[derive(Clone, Debug)]
pub(crate) struct Pair {
    pub a: Use,
    pub b: Use,
    pub d2: Rational,
}

/// A measured gap, before the policy table.
#[derive(Clone, Debug)]
pub(crate) struct Gap {
    pub solid: usize,
    pub pairs: Vec<Pair>,
    pub kind: GapKind,
    pub d2: Rational,
    /// Squared length of the shortest straight edge of the faces involved.
    pub feature2: Option<Rational>,
    pub entities: Vec<String>,
}

fn is_seam(edges: &[OrientedEdge], i: usize) -> bool {
    edges
        .iter()
        .enumerate()
        .any(|(j, o)| j != i && o.edge == edges[i].edge && o.forward != edges[i].forward)
}

fn dist2(a: &Vec3, b: &Vec3) -> Rational {
    let d = sub(a, b);
    dot(&d, &d)
}

fn normal(face: &Face) -> Option<Vec3> {
    match &face.surface {
        Surface::Plane(p) => Some(if face.same_sense {
            p.axis.clone()
        } else {
            p.axis.clone().map(|c| -c)
        }),
        _ => None,
    }
}

/// Measures the gaps of every solid, in file order.
pub(crate) fn measure(file: &StepFile, radius: &Rational) -> Result<Vec<Gap>, ImportRefusal> {
    let mut gaps = Vec::new();
    for (si, solid) in file.solids.iter().enumerate() {
        gaps.extend(measure_solid(si, solid, radius)?);
    }
    Ok(gaps)
}

fn measure_solid(si: usize, solid: &Solid, radius: &Rational) -> Result<Vec<Gap>, ImportRefusal> {
    // Directed uses by unordered endpoint pair, for straight edges.
    let mut segments: BTreeMap<(Vec3, Vec3), Vec<Use>> = BTreeMap::new();
    let mut rings: BTreeMap<String, (i32, i32)> = BTreeMap::new();
    let mut ring_faces: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for (fi, face) in solid.faces.iter().enumerate() {
        for bound in &face.loops {
            for (i, e) in bound.edges.iter().enumerate() {
                if is_seam(&bound.edges, i) {
                    continue;
                }
                let travel = e.forward == bound.orientation;
                let (from, to) = if travel {
                    (e.start.point.clone(), e.end.point.clone())
                } else {
                    (e.end.point.clone(), e.start.point.clone())
                };
                if matches!(e.curve, Curve::Line(..)) {
                    let key = (from.clone().min(to.clone()), from.clone().max(to.clone()));
                    segments.entry(key).or_default().push(Use {
                        face: fi,
                        edge: e.edge.clone(),
                        from,
                        to,
                    });
                } else {
                    ring_faces
                        .entry(e.edge.clone())
                        .or_default()
                        .insert(face.name.clone());
                    let counts = rings.entry(e.edge.clone()).or_default();
                    if travel {
                        counts.0 += 1;
                    } else {
                        counts.1 += 1;
                    }
                }
            }
        }
    }
    if let Some((name, _)) = rings.iter().find(|(_, c)| **c != (1, 1)) {
        return Err(not_solid(format!(
            "gap at curved edge {name} (faces {}): curved gaps are outside the first path, only straight seams between planar faces are repaired",
            ring_faces[name]
                .iter()
                .cloned()
                .collect::<Vec<_>>()
                .join(", ")
        )));
    }
    let mut free: Vec<Use> = Vec::new();
    for uses in segments.into_values() {
        let (mut forward, mut backward): (Vec<_>, Vec<_>) =
            uses.into_iter().partition(|u| u.from < u.to);
        while !forward.is_empty() && !backward.is_empty() {
            forward.pop();
            backward.pop();
        }
        free.extend(forward);
        free.extend(backward);
    }
    if free.is_empty() {
        return Ok(Vec::new());
    }
    let limit2 = radius * radius;
    let mut candidates: Vec<(Rational, usize, usize)> = Vec::new();
    for i in 0..free.len() {
        for j in i + 1..free.len() {
            if free[i].face == free[j].face {
                continue;
            }
            let (a, b) = (
                dist2(&free[i].from, &free[j].to),
                dist2(&free[i].to, &free[j].from),
            );
            let d2 = a.clone().max(b);
            if d2 <= limit2 {
                candidates.push((d2, i, j));
            }
        }
    }
    candidates.sort();
    let mut taken = vec![false; free.len()];
    let mut pairs: Vec<Pair> = Vec::new();
    for (d2, i, j) in candidates {
        if !taken[i] && !taken[j] {
            taken[i] = true;
            taken[j] = true;
            pairs.push(Pair {
                a: free[i].clone(),
                b: free[j].clone(),
                d2,
            });
        }
    }
    // Join free uses that share a point, or are paired, into one gap.
    let mut group: Vec<usize> = (0..free.len()).collect();
    fn root(g: &mut [usize], mut x: usize) -> usize {
        while g[x] != x {
            g[x] = g[g[x]];
            x = g[x];
        }
        x
    }
    let union = |group: &mut Vec<usize>, a: usize, b: usize| {
        let (x, y) = (root(group, a), root(group, b));
        group[x] = y;
    };
    let mut by_point: BTreeMap<Vec3, usize> = BTreeMap::new();
    for (k, u) in free.iter().enumerate() {
        for point in [&u.from, &u.to] {
            match by_point.get(point) {
                Some(&other) => union(&mut group, k, other),
                None => {
                    by_point.insert(point.clone(), k);
                }
            }
        }
    }
    let index_of = |u: &Use, free: &[Use]| {
        free.iter()
            .position(|f| f.face == u.face && f.edge == u.edge && f.from == u.from)
            .expect("pair member is a free use")
    };
    for p in &pairs {
        let (a, b) = (index_of(&p.a, &free), index_of(&p.b, &free));
        union(&mut group, a, b);
    }
    let mut grouped: BTreeMap<usize, (Vec<Pair>, Vec<Use>)> = BTreeMap::new();
    for (k, u) in free.iter().enumerate() {
        let r = root(&mut group, k);
        let entry = grouped.entry(r).or_default();
        if !taken[k] {
            entry.1.push(u.clone());
        }
    }
    for p in pairs {
        let r = root(&mut group, index_of(&p.a, &free));
        grouped.entry(r).or_default().0.push(p);
    }
    let mut gaps = Vec::new();
    for (pairs, unmatched) in grouped.into_values() {
        let mut kind = if unmatched.is_empty() {
            GapKind::SeamClose
        } else {
            GapKind::Opening
        };
        let mut faces: BTreeSet<usize> = unmatched.iter().map(|u| u.face).collect();
        for p in &pairs {
            faces.insert(p.a.face);
            faces.insert(p.b.face);
            let (Some(na), Some(nb)) = (
                normal(&solid.faces[p.a.face]),
                normal(&solid.faces[p.b.face]),
            ) else {
                return Err(not_solid(format!(
                    "gap between faces {} and {} (edges {}, {}) touches a curved face: curved gaps are outside the first path, only planar seams are repaired",
                    solid.faces[p.a.face].name, solid.faces[p.b.face].name, p.a.edge, p.b.edge
                )));
            };
            // Faces that would fold onto each other close an opening, not a seam.
            if crate::exact::is_zero(&cross(&na, &nb))
                && dot(&na, &nb).sign() == catallax_numbers::Sign::Negative
            {
                kind = GapKind::Opening;
            }
        }
        let d2 = pairs
            .iter()
            .map(|p| p.d2.clone())
            .max()
            .unwrap_or_else(Rational::zero);
        let mut feature2: Option<Rational> = None;
        let mut entities = BTreeSet::new();
        for &f in &faces {
            let face = &solid.faces[f];
            entities.insert(face.name.clone());
            for bound in &face.loops {
                for e in &bound.edges {
                    if matches!(e.curve, Curve::Line(..)) {
                        let l2 = dist2(&e.start.point, &e.end.point);
                        feature2 = Some(match feature2 {
                            Some(m) if m <= l2 => m,
                            _ => l2,
                        });
                    }
                }
            }
        }
        for p in &pairs {
            entities.insert(p.a.edge.clone());
            entities.insert(p.b.edge.clone());
        }
        for u in &unmatched {
            entities.insert(u.edge.clone());
        }
        let mut entities: Vec<String> = entities.into_iter().collect();
        entities.sort_by_key(|n| n.trim_start_matches('#').parse::<u64>().unwrap_or(u64::MAX));
        gaps.push(Gap {
            solid: si,
            pairs,
            kind,
            d2,
            feature2,
            entities,
        });
    }
    gaps.sort_by(|a, b| {
        let key = |g: &Gap| {
            g.entities
                .first()
                .and_then(|n| n.trim_start_matches('#').parse::<u64>().ok())
        };
        (a.solid, key(a)).cmp(&(b.solid, key(b)))
    });
    Ok(gaps)
}

const SCALE: i64 = 1_000_000_000_000;

/// Decimal enclosure of `sqrt(q)` at 1e-12 mm: exact when `q` is a perfect square there.
pub(crate) fn sqrt_enclosure(q: &Rational) -> (Rational, Rational) {
    let s = Rational::from_i64(SCALE);
    let scaled = &(q * &s) * &s;
    let k = isqrt(&scaled.floor());
    let lo = Rational::new(k.clone(), Integer::from_i64(SCALE)).expect("nonzero");
    let exact = Rational::from_integer(&k * &k) == scaled;
    let hi = if exact {
        lo.clone()
    } else {
        Rational::new(&k + &Integer::one(), Integer::from_i64(SCALE)).expect("nonzero")
    };
    (lo, hi)
}

pub(crate) fn decimal(q: &Rational) -> String {
    q.to_decimal_exact().expect("terminating decimal")
}

/// The policy table for one gap: automatic or ask, and why.
pub(crate) fn disposition(gap: &Gap) -> (Disposition, String) {
    let cap = auto_cap();
    match gap.kind {
        GapKind::Opening => (
            Disposition::Ask,
            "the faces fold onto each other, so closing would close an opening".into(),
        ),
        GapKind::BodyGap => (
            Disposition::Ask,
            "bodies are never merged automatically".into(),
        ),
        GapKind::SeamClose => {
            if gap.d2 > &cap * &cap {
                return (
                    Disposition::Ask,
                    "gap is above the 0.01 mm automatic cap".into(),
                );
            }
            match &gap.feature2 {
                Some(l2) if &Rational::from_i64(100) * &gap.d2 <= *l2 => (
                    Disposition::Auto,
                    "boundary gap within the automatic cap, under a tenth of the nearest feature"
                        .into(),
                ),
                Some(_) => (
                    Disposition::Ask,
                    "the move would exceed a tenth of the nearest feature size".into(),
                ),
                None => (
                    Disposition::Ask,
                    "the nearest feature size could not be established".into(),
                ),
            }
        }
    }
}

/// Question id for a gap of the source with hash `sha256`.
pub(crate) fn question_id(sha256: &str, id: &str) -> String {
    format!("repair:{}:{id}", &sha256[..12])
}

pub(crate) fn proposals(file: &StepFile, source_sha: &str, gaps: &[Gap]) -> Vec<Proposal> {
    let _ = file;
    gaps.iter()
        .enumerate()
        .map(|(i, gap)| {
            let id = format!("g{i}");
            let (disposition, reason) = disposition(gap);
            let (lo, hi) = sqrt_enclosure(&gap.d2);
            let tolerance = if disposition == Disposition::Auto {
                crate::AUTO_CAP_MM.to_string()
            } else {
                decimal(&hi)
            };
            Proposal {
                question: question_id(source_sha, &id),
                id,
                entities: gap.entities.clone(),
                kind: gap.kind,
                distance_min_mm: decimal(&lo),
                distance_max_mm: decimal(&hi),
                tolerance_mm: tolerance,
                disposition,
                reason,
            }
        })
        .collect()
}

/// Reads the source, measures every gap and applies the policy table. Nothing moves.
pub fn propose(source: &[u8], options: &ImportOptions) -> Result<ProposalSet, ImportRefusal> {
    let file = crate::read_step(source)?;
    let sha = sha256_hex(source);
    let gaps = measure(&file, &options.search_radius_mm)?;
    Ok(ProposalSet {
        format: crate::PROPOSALS_FORMAT.to_string(),
        proposals: proposals(&file, &sha, &gaps),
        source_sha256: sha,
    })
}
