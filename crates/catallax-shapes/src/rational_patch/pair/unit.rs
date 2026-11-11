//! One first-patch cell with its cover of second-patch cells: the cover's
//! exclusion, each unit's certificate and border points, or a split.

use super::*;
use face::{Endpoint, Hit, classify, contact, excluded, face_excluded, overlap, solve, strict};
use net::{interval_cross, interval_dot, less};
use std::collections::{BTreeMap, VecDeque};

/// Most splits of one first-patch cell before it refuses.
const MAX_DEPTH: u32 = 48;
/// Most splits of one second-patch cover cell.
const MAX_COVER_DEPTH: u32 = 64;

/// A first-patch cell and the second-patch cells that may still meet it.
/// Cells of the second patch outside the cover were excluded from an
/// ancestor of this cell, so they hold no meeting point with it.
#[derive(Clone, Debug)]
pub(super) struct Node {
    pub(super) cell: Arc<Cell>,
    pub(super) cover: Arc<[Arc<Cell>]>,
}

impl Node {
    pub(super) fn describe(&self, cells: &[Arc<Cell>]) -> String {
        let interval = |r: &[Rational; 2]| format!("[{}, {}]", r[0], r[1]);
        let rect = &self.cell.rect;
        let mut text = format!(
            "first span {} cell u {} v {} after {} splits",
            self.cell.span,
            interval(&rect[0]),
            interval(&rect[1]),
            self.cell.depth
        );
        if let Some(first) = cells.first() {
            let hull = hull(cells);
            text.push_str(&format!(
                ", against {} cells of second span {} in s {} t {}",
                cells.len(),
                first.span,
                interval(&hull[0]),
                interval(&hull[1])
            ));
        }
        text
    }
}

/// One accepted unit: a first-patch cell against one knot span of the
/// second patch.
#[derive(Debug)]
pub(super) struct Unit {
    pub(super) second_span: usize,
    pub(super) hull: Domain,
    pub(super) certificate: PairCertificate,
    pub(super) ends: Option<[Endpoint; 2]>,
}

/// What one cell proved.
#[derive(Debug)]
pub(super) enum Step {
    /// Every cover cell was excluded: no meeting point.
    Excluded,
    /// Every unit settled, each with no arc or one arc.
    Settled(Vec<Unit>),
    /// These cells cover it.
    Split(Vec<Node>),
    /// Still unsettled at the split limit, or holding an exact tangent contact.
    Unresolved(Refusal),
}

/// A cell's result, the cover cells it kept and its face systems.
#[derive(Debug)]
pub(super) struct Done {
    pub(super) step: Step,
    pub(super) cover: usize,
    pub(super) systems: usize,
}

/// The cover cells that may meet `alpha`, each split until it is no larger
/// than `size`. Order is deterministic: breadth first in cover order.
pub(super) fn refine(
    alpha: &Cell,
    cover: &[Arc<Cell>],
    size: &Rational,
    dec: &mut Decider,
) -> Result<Vec<Arc<Cell>>, Refusal> {
    let mut queue: VecDeque<Arc<Cell>> = cover.iter().cloned().collect();
    let mut kept = Vec::new();
    while let Some(beta) = queue.pop_front() {
        if excluded(alpha, &beta, dec)? {
            continue;
        }
        if beta.depth < MAX_COVER_DEPTH && less(size, &beta.extent, dec)? {
            // The second patch splits at two and four fifths; the first in
            // thirds. The rules never share a split value and never split at
            // a half, so cell faces of the two patches meet only on knot
            // lines, and a symmetric centre is never on a face.
            for child in beta.second_children(dec)? {
                queue.push_back(child.clone());
            }
            continue;
        }
        kept.push(beta);
    }
    Ok(kept)
}

/// The parameter box round some cells of one knot span.
fn hull(cells: &[Arc<Cell>]) -> Domain {
    let mut out = cells[0].rect.clone();
    for cell in &cells[1..] {
        for (range, [lo, hi]) in out.iter_mut().zip(&cell.rect) {
            if *lo < range[0] {
                range[0] = lo.clone();
            }
            if *hi > range[1] {
                range[1] = hi.clone();
            }
        }
    }
    out
}

/// A component of a field that keeps one strict sign on every cell.
fn common_axis(fields: &[[[Rational; 2]; 2]]) -> Option<[i8; 2]> {
    for (component, direction) in [(0, [1, 0]), (1, [0, 1])] {
        for sign in [Sign::Positive, Sign::Negative] {
            let holds = fields.iter().all(|field| {
                let range = &field[component];
                if sign == Sign::Positive {
                    range[0].sign() == Sign::Positive
                } else {
                    range[1].sign() == Sign::Negative
                }
            });
            if holds {
                let d: [i8; 2] = if sign == Sign::Positive {
                    direction
                } else {
                    direction.map(|x| -x)
                };
                return Some(d);
            }
        }
    }
    None
}

/// The first certificate that holds over `alpha` against every cell.
fn certificate(
    alpha: &Cell,
    cells: &[Arc<Cell>],
    dec: &mut Decider,
) -> Result<Option<PairCertificate>, Refusal> {
    let mut fields = Vec::with_capacity(cells.len());
    for beta in cells {
        fields.push([
            interval_dot(&alpha.derived().tu_range, &beta.derived().normal_range, dec)?,
            interval_dot(&alpha.derived().tv_range, &beta.derived().normal_range, dec)?,
        ]);
    }
    if let Some(direction) = common_axis(&fields) {
        return Ok(Some(PairCertificate::FieldFirst(direction)));
    }
    fields.clear();
    for beta in cells {
        fields.push([
            interval_dot(&beta.derived().tu_range, &alpha.derived().normal_range, dec)?,
            interval_dot(&beta.derived().tv_range, &alpha.derived().normal_range, dec)?,
        ]);
    }
    if let Some(direction) = common_axis(&fields) {
        return Ok(Some(PairCertificate::FieldSecond(direction)));
    }
    'component: for component in 0..3 {
        let mut sign = None;
        for beta in cells {
            let range = interval_cross(
                &alpha.derived().normal_range,
                &beta.derived().normal_range,
                component,
                dec,
            )?;
            if !strict(&range) {
                continue 'component;
            }
            let here = range[0].sign() == Sign::Positive;
            if sign.is_some_and(|s| s != here) {
                continue 'component;
            }
            sign = Some(here);
        }
        return Ok(Some(PairCertificate::NormalCone(component as u8)));
    }
    Ok(None)
}

/// Proves one first-patch cell against its cover, or splits it.
pub(super) fn process(ctx: &Context, node: &Node, dec: &mut Decider) -> Result<Done, Refusal> {
    dec.check_deadline()?;
    let alpha = &node.cell;
    // The writer refined the cover against this cell before the level ran.
    let cover: Vec<Arc<Cell>> = node.cover.to_vec();
    let mut systems = 0;
    if cover.is_empty() {
        return Ok(Done {
            step: Step::Excluded,
            cover: 0,
            systems,
        });
    }
    let mut groups: BTreeMap<usize, Vec<Arc<Cell>>> = BTreeMap::new();
    for beta in &cover {
        groups.entry(beta.span).or_default().push(beta.clone());
    }
    let mut units = Vec::with_capacity(groups.len());
    let mut failing = None;
    for (span, cells) in &groups {
        let Some(certificate) = certificate(alpha, cells, dec)? else {
            failing = Some(cells.as_slice());
            break;
        };
        let ends = match border(ctx, alpha, *span, cells, &mut systems, dec)? {
            Some(ends) => ends,
            None => {
                failing = Some(cells.as_slice());
                break;
            }
        };
        let ends = match ends.len() {
            0 => None,
            2 if ends[0].start != ends[1].start => {
                let mut ends = ends;
                if !ends[0].start {
                    ends.swap(0, 1);
                }
                Some(ends.try_into().expect("two ends"))
            }
            _ => {
                failing = Some(cells.as_slice());
                break;
            }
        };
        units.push(Unit {
            second_span: *span,
            hull: hull(cells),
            certificate,
            ends,
        });
    }
    let kept = cover.len();
    let Some(failing) = failing else {
        return Ok(Done {
            step: Step::Settled(units),
            cover: kept,
            systems,
        });
    };
    let step = split(ctx, node, cover, failing, dec)?;
    Ok(Done {
        step,
        cover: kept,
        systems,
    })
}

/// Splits the cell into thirds along its longer side; the children inherit
/// the refined cover. A cell holding an exact tangent contact, or at the
/// split limit, refuses and names itself.
fn split(
    ctx: &Context,
    node: &Node,
    cover: Vec<Arc<Cell>>,
    failing: &[Arc<Cell>],
    dec: &mut Decider,
) -> Result<Step, Refusal> {
    let alpha = &node.cell;
    let hull = hull(failing);
    let rect = [
        alpha.rect[0].clone(),
        alpha.rect[1].clone(),
        hull[0].clone(),
        hull[1].clone(),
    ];
    let spans = [alpha.span, failing[0].span];
    if alpha.depth >= 1
        && let Some(refusal) = contact(ctx, spans, &rect, &node.describe(failing), dec)?
    {
        return Ok(Step::Unresolved(refusal));
    }
    if alpha.depth >= MAX_DEPTH {
        return Ok(Step::Unresolved(Refusal::new(
            "intersect.spline-pair-unresolved",
            format!(
                "no exclusion, field or normal-cone certificate with two transversal border points settled the cell: {}",
                node.describe(failing)
            ),
        )));
    }
    let axis = alpha.split_axis(dec)?;
    let thirds = [
        Rational::from_ratio(1, 3).expect("nonzero"),
        Rational::from_ratio(2, 3).expect("nonzero"),
    ];
    let children = alpha.split(axis, &thirds, dec)?;
    let cover: Arc<[Arc<Cell>]> = cover.into();
    Ok(Step::Split(
        children
            .into_iter()
            .map(|cell| Node {
                cell: Arc::new(cell),
                cover: cover.clone(),
            })
            .collect(),
    ))
}

/// Border points of one unit, `alpha` against the cells of one second-patch
/// knot span, where branches cross; `None` when the border does not settle
/// at this size. Faces are the four sides of `alpha` and the sides of the
/// knot span; splits of the second patch inside a span are not faces.
fn border(
    ctx: &Context,
    alpha: &Cell,
    span: usize,
    cells: &[Arc<Cell>],
    systems: &mut usize,
    dec: &mut Decider,
) -> Result<Option<Vec<Endpoint>>, Refusal> {
    let domain = &ctx.sides[1].spans[span].root.rect;
    let unit = [
        alpha.rect[0].clone(),
        alpha.rect[1].clone(),
        domain[0].clone(),
        domain[1].clone(),
    ];
    let spans = [alpha.span, span];
    let mut hits: Vec<Hit> = Vec::new();
    // Sides of the first-patch cell, against the cover cells they may meet.
    for edge in 0..4 {
        let mut near = Vec::new();
        for beta in cells {
            if !face_excluded(alpha, edge, beta, dec)? {
                near.push(beta.clone());
            }
        }
        if near.is_empty() {
            continue;
        }
        let around = hull(&near);
        let search = [
            alpha.rect[0].clone(),
            alpha.rect[1].clone(),
            around[0].clone(),
            around[1].clone(),
        ];
        let (axis, side) = (edge / 2, edge % 2);
        let value = alpha.rect[axis][side].clone();
        *systems += 1;
        let Some(found) = solve(ctx, spans, &search, axis, &value, dec)? else {
            return Ok(None);
        };
        if let Some(hit) = found
            && !add_hit(&mut hits, hit, dec)?
        {
            return Ok(None);
        }
    }
    // Sides of the second patch's knot span, where the cover reaches them.
    for edge in 0..4 {
        let (b_axis, side) = (edge / 2, edge % 2);
        let value = domain[b_axis][side].clone();
        let mut near = Vec::new();
        for beta in cells {
            if beta.rect[b_axis][side] == value && !face_excluded(beta, edge, alpha, dec)? {
                near.push(beta.clone());
            }
        }
        if near.is_empty() {
            continue;
        }
        let around = hull(&near);
        let search = [
            alpha.rect[0].clone(),
            alpha.rect[1].clone(),
            around[0].clone(),
            around[1].clone(),
        ];
        *systems += 1;
        let Some(found) = solve(ctx, spans, &search, 2 + b_axis, &value, dec)? else {
            return Ok(None);
        };
        if let Some(hit) = found
            && !add_hit(&mut hits, hit, dec)?
        {
            return Ok(None);
        }
    }
    let mut ends = Vec::with_capacity(hits.len());
    for hit in hits {
        match classify(ctx, &unit, hit, dec)? {
            Some(Some(end)) => ends.push(end),
            Some(None) => (),
            None => return Ok(None),
        }
    }
    Ok(Some(ends))
}

/// Adds a border point unless the same corner point was already found on
/// another face. `false` when two points cannot be told apart.
fn add_hit(hits: &mut Vec<Hit>, hit: Hit, dec: &mut Decider) -> Result<bool, Refusal> {
    for old in hits.iter() {
        if old.fixed == hit.fixed && old.bounds == hit.bounds {
            return Ok(true);
        }
        if overlap(&old.bounds, &hit.bounds, dec)? {
            return Ok(false);
        }
    }
    hits.push(hit);
    Ok(true)
}
