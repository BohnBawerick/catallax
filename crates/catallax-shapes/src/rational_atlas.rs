//! Exact joined rational charts. Numerical chart cuts are not model edges.

use std::{collections::VecDeque, sync::Arc};

use catallax_numbers::{Decider, Rational, Refusal};

mod query;
mod section;
use crate::rational_patch::MAX_PROOF_CELLS;
use crate::{BSplineSurface, RationalPatch, SplineBoundary};
pub use section::AtlasSection;

type Domain = [[Rational; 2]; 2];

/// One certified chart and its exact G1 identifications within a logical face.
#[derive(Clone, Debug)]
pub struct AtlasChart {
    patch: RationalPatch,
    neighbors: [Option<usize>; 4],
}

impl AtlasChart {
    /// The unchanged exact finite chart, including its embedding certificate.
    pub fn patch(&self) -> &RationalPatch {
        &self.patch
    }

    /// A joined chart on this side, or a physical face boundary.
    pub fn neighbor(&self, side: SplineBoundary) -> Option<usize> {
        self.neighbors[side_index(side)]
    }
}

#[derive(Debug)]
struct Data {
    source: BSplineSurface,
    charts: Vec<AtlasChart>,
    periodic: [bool; 2],
}

/// One connected G1 logical face, covered by regular exact rational charts.
/// This proves chart regularity and gluing, not global body admissibility.
#[derive(Clone, Debug)]
pub struct RationalAtlas {
    data: Arc<Data>,
}

#[derive(Clone)]
struct Band {
    source_span: usize,
    bounds: [Rational; 2],
}

type CertifiedCover = (Vec<Band>, Vec<Band>, Vec<RationalPatch>);

fn side_index(side: SplineBoundary) -> usize {
    match side {
        SplineBoundary::UMin => 0,
        SplineBoundary::UMax => 1,
        SplineBoundary::VMin => 2,
        SplineBoundary::VMax => 3,
    }
}

fn unproven(detail: &str) -> Refusal {
    Refusal::new("shapes.rational-atlas-unproven", detail)
}

fn bands(source: &BSplineSurface, axis: usize) -> Vec<Band> {
    let degree = source.degrees()[axis];
    let knots = source.knots()[axis];
    (degree..source.control_count()[axis])
        .filter(|&i| knots[i] < knots[i + 1])
        .map(|i| Band {
            source_span: i,
            bounds: [knots[i].clone(), knots[i + 1].clone()],
        })
        .collect()
}

fn domain(u: &Band, v: &Band) -> Domain {
    [u.bounds.clone(), v.bounds.clone()]
}

fn split_band(bands: &mut Vec<Band>, index: usize) {
    let old = &mut bands[index];
    let middle = (&old.bounds[0] + &old.bounds[1]).mul_pow2(-1);
    let high = Band {
        source_span: old.source_span,
        bounds: [middle.clone(), old.bounds[1].clone()],
    };
    old.bounds[1] = middle;
    bands.insert(index + 1, high);
}

fn cover(source: &BSplineSurface, dec: &mut Decider) -> Result<CertifiedCover, Refusal> {
    let mut u = bands(source, 0);
    let mut v = bands(source, 1);
    let mut cells: Vec<Option<RationalPatch>> = vec![None; u.len() * v.len()];
    let widths = source.domain().map(|[lo, hi]| hi - lo);
    loop {
        dec.check_deadline()?;
        let mut failed = None;
        'rows: for (i, ub) in u.iter().enumerate() {
            for (j, vb) in v.iter().enumerate() {
                let slot = &mut cells[i * v.len() + j];
                if slot.is_some() {
                    continue;
                }
                let bounds = domain(ub, vb);
                let carrier = source.bezier_patch(
                    [ub.source_span, vb.source_span],
                    bounds.each_ref().map(|x| x.each_ref()),
                    dec,
                )?;
                match RationalPatch::new_cell(carrier, dec) {
                    Ok(patch) => *slot = Some(patch),
                    Err(error) if error.code == "shapes.rational-patch-unproven" => {
                        failed = Some((i, j));
                        break 'rows;
                    }
                    Err(error) => return Err(error),
                }
            }
        }
        let Some((i, j)) = failed else {
            return Ok((
                u,
                v,
                cells
                    .into_iter()
                    .map(|x| x.expect("complete cover"))
                    .collect(),
            ));
        };
        let axis = usize::from(
            &(&v[j].bounds[1] - &v[j].bounds[0]) * &widths[0]
                > &(&u[i].bounds[1] - &u[i].bounds[0]) * &widths[1],
        );
        let count = if axis == 0 {
            (u.len() + 1) * v.len()
        } else {
            u.len() * (v.len() + 1)
        };
        if count > MAX_PROOF_CELLS || dec.policy().max_precision_bits == 0 {
            return Err(unproven(
                "the complete regular chart cover could not be proved",
            ));
        }
        let old_u = u.len();
        let old_v = v.len();
        if axis == 0 {
            split_band(&mut u, i);
        } else {
            split_band(&mut v, j);
        }
        let mut next = vec![None; u.len() * v.len()];
        for a in 0..u.len() {
            for b in 0..v.len() {
                let old_a = if axis == 0 && a > i { a - 1 } else { a };
                let old_b = if axis == 1 && b > j { b - 1 } else { b };
                debug_assert!(old_a < old_u && old_b < old_v);
                if let Some(patch) = &cells[old_a * old_v + old_b] {
                    let bounds = domain(&u[a], &v[b]);
                    next[a * v.len() + b] = Some(if patch.domain() == &bounds {
                        patch.clone()
                    } else {
                        patch.restricted(bounds, dec)?
                    });
                }
            }
        }
        cells = next;
    }
}

impl RationalAtlas {
    /// Covers the complete tensor domain and partitions only at genuine non-G1
    /// joins. `periodic` identifies the exact end rows of each parameter; it never
    /// closes a geometric gap. Every accepted identification is proved C0/G1.
    pub fn partition(
        source: BSplineSurface,
        periodic: [bool; 2],
        dec: &mut Decider,
    ) -> Result<Vec<Self>, Refusal> {
        let (u, v, patches) = cover(&source, dec)?;
        let mut links = vec![[None; 4]; patches.len()];
        for (i, ub) in u.iter().enumerate() {
            for (j, vb) in v.iter().enumerate() {
                let here = i * v.len() + j;
                for axis in 0..2 {
                    let (next, wrap, same_span, knot) = if axis == 0 {
                        if i + 1 == u.len() && !periodic[0] {
                            continue;
                        }
                        let ni = (i + 1) % u.len();
                        (
                            ni * v.len() + j,
                            ni == 0,
                            ub.source_span == u[ni].source_span,
                            &ub.bounds[1],
                        )
                    } else {
                        if j + 1 == v.len() && !periodic[1] {
                            continue;
                        }
                        let nj = (j + 1) % v.len();
                        (
                            i * v.len() + nj,
                            nj == 0,
                            vb.source_span == v[nj].source_span,
                            &vb.bounds[1],
                        )
                    };
                    let multiplicity = source.knots()[axis].iter().filter(|k| *k == knot).count();
                    let smooth = if !wrap && (same_span || multiplicity < source.degrees()[axis]) {
                        true
                    } else {
                        patches[here].seam(
                            &patches[next],
                            axis,
                            wrap || multiplicity > source.degrees()[axis],
                            dec,
                        )?
                    };
                    if smooth {
                        links[here][2 * axis + 1] = Some(next);
                        links[next][2 * axis] = Some(here);
                    }
                }
            }
        }
        let mut seen = vec![false; patches.len()];
        let mut result = Vec::new();
        for start in 0..patches.len() {
            if seen[start] {
                continue;
            }
            let mut queue = VecDeque::from([start]);
            let mut members = Vec::new();
            seen[start] = true;
            while let Some(index) = queue.pop_front() {
                members.push(index);
                for next in links[index].iter().flatten() {
                    if !seen[*next] {
                        seen[*next] = true;
                        queue.push_back(*next);
                    }
                }
            }
            members.sort_unstable();
            let mut local = vec![usize::MAX; patches.len()];
            for (index, global) in members.iter().enumerate() {
                local[*global] = index;
            }
            let charts = members
                .into_iter()
                .map(|global| AtlasChart {
                    patch: patches[global].clone(),
                    neighbors: links[global].map(|next| next.map(|index| local[index])),
                })
                .collect();
            result.push(Self {
                data: Arc::new(Data {
                    source: source.clone(),
                    charts,
                    periodic,
                }),
            });
        }
        Ok(result)
    }

    /// The retained full tensor definition. A component can cover only a subset.
    pub fn source(&self) -> &BSplineSurface {
        &self.data.source
    }
    /// Complete regular chart cover of this logical face.
    pub fn charts(&self) -> &[AtlasChart] {
        &self.data.charts
    }
    /// Whether the source's u end rows are identified by an exact position proof.
    pub fn periodic_u(&self) -> bool {
        self.data.periodic[0]
    }

    /// Whether the source's v end rows are identified by an exact position proof.
    pub fn periodic_v(&self) -> bool {
        self.data.periodic[1]
    }

    pub(crate) fn has_periodic_join(&self, axis: usize) -> bool {
        let first = self.source().domain()[axis][0];
        let side = [SplineBoundary::UMin, SplineBoundary::VMin][axis];
        self.data.periodic[axis]
            && self.charts().iter().any(|chart| {
                &chart.patch().domain()[axis][0] == first && chart.neighbor(side).is_some()
            })
    }

    /// Exact analytic reduction when all Euclidean poles are rational and
    /// coplanar. The clamped corner tangents fix the natural chart orientation;
    /// the positive control hull proves the complete face lies on this plane.
    pub fn rational_plane(&self) -> Option<crate::Plane> {
        let carrier = self.charts().first()?.patch().carrier();
        let nv = carrier.control_count()[1];
        let origin = carrier.controls()[0].as_rational()?;
        let u = carrier.controls()[nv].as_rational()?.minus(&origin);
        let v = carrier.controls()[1].as_rational()?.minus(&origin);
        let normal = u.cross(&v);
        if normal.is_zero() {
            return None;
        }
        for chart in self.charts() {
            for point in chart.patch().carrier().controls() {
                if !normal.dot(&point.as_rational()?.minus(&origin)).is_zero() {
                    return None;
                }
            }
        }
        Some(crate::Plane {
            origin: origin.to_real(),
            normal: normal.to_real(),
            ref_dir: u.to_real(),
        })
    }
}
