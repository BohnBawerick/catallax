//! The level-by-level search. Each level's cells are independent proofs:
//! helper threads prove them ahead on every core, and the request's own
//! thread adopts their results and charges in cell order, so the answer, a
//! refusal and the work count are those of one thread.

use super::*;
use catallax_numbers::scoped_try_collect;
use join::Resolved;
use unit::{Done, Node, Step};

/// The complete search: every accepted arc in a fixed order, and the counts.
pub(super) fn search(
    ctx: &Context,
    reuse_proofs: bool,
    dec: &mut Decider,
) -> Result<(Vec<Resolved>, PatchPairSearch), Refusal> {
    let mut record = PatchPairSearch::default();
    let mut cache = reuse_proofs.then(reuse::Cache::default);
    let cover: Arc<[Arc<Cell>]> = ctx.sides[1]
        .spans
        .iter()
        .map(|span| span.root.clone())
        .collect();
    let mut level: Vec<Node> = ctx.sides[0]
        .spans
        .iter()
        .map(|span| Node {
            cell: span.root.clone(),
            cover: cover.clone(),
        })
        .collect();
    let mut arcs = Vec::new();
    while !level.is_empty() {
        record.levels += 1;
        record.cells += level.len();
        if record.cells > MAX_CELLS {
            return Err(Refusal::new(
                "intersect.spline-pair-unresolved",
                format!(
                    "the cell search reached its unchanged proof allowance of {MAX_CELLS} cells at {}",
                    level[0].describe(&[])
                ),
            ));
        }
        // The writer refines every cell's cover first, in order. Children
        // of second-patch cells are shared and built once, so their charge
        // falls on the first cell in this order that needs them.
        for node in &mut level {
            let cover = unit::refine(ctx, &node.cell, &node.cover, &node.cell.extent, dec)
                .map_err(|e| named(e, node))?;
            // Tangent and normal nets, and the span pair's equations, only
            // for cells that are proved against something, so excluded
            // cells never pay for them.
            if !cover.is_empty() {
                node.cell.ensure_derived(dec).map_err(|e| named(e, node))?;
                for beta in &cover {
                    beta.ensure_derived(dec).map_err(|e| named(e, node))?;
                    ctx.ensure_equations([node.cell.span, beta.span], dec)
                        .map_err(|e| named(e, node))?;
                }
            }
            node.cover = cover.into();
        }
        let snapshot = cache.as_ref().map(reuse::Cache::snapshot);
        let done = run_level(ctx, &level, snapshot.as_ref(), dec)?;
        drop(snapshot);
        let mut next = Vec::new();
        for (
            node,
            Processed {
                done:
                    Done {
                        step,
                        cover,
                        systems,
                    },
                requests,
            },
        ) in level.into_iter().zip(done)
        {
            if let Some(cache) = &mut cache {
                cache.commit(requests);
            }
            record.depth = record.depth.max(node.cell.depth);
            record.face_systems += systems;
            record.cover += cover;
            match step {
                Step::Excluded => record.excluded += 1,
                Step::Settled(units) => {
                    for unit in units {
                        match unit.certificate {
                            PairCertificate::FieldFirst(_) => record.field_first += 1,
                            PairCertificate::FieldSecond(_) => record.field_second += 1,
                            PairCertificate::NormalCone(_) => record.normal_cone += 1,
                        }
                        if let Some(ends) = unit.ends {
                            record.arcs += 1;
                            arcs.push(Resolved {
                                first: node.cell.clone(),
                                second_span: unit.second_span,
                                hull: unit.hull,
                                certificate: unit.certificate,
                                ends,
                            });
                        }
                    }
                }
                Step::Split(children) => {
                    record.splits += 1;
                    next.extend(children);
                }
                Step::Unresolved(refusal) => return Err(refusal),
            }
        }
        level = next;
    }
    Ok((arcs, record))
}

/// An allowance or deadline refusal keeps its code and names the cell pair
/// it stopped in.
fn named(mut refusal: Refusal, node: &Node) -> Refusal {
    if refusal.code == "numbers.operation-budget" || refusal.code == "numbers.deadline" {
        refusal.detail = format!("{}; in the cell {}", refusal.detail, node.describe(&[]));
    }
    refusal
}

struct Processed {
    done: Done,
    requests: Vec<Arc<reuse::Certificate>>,
}

fn process(
    ctx: &Context,
    node: &Node,
    snapshot: Option<&reuse::Snapshot>,
    dec: &mut Decider,
) -> Result<Processed, Refusal> {
    let (done, requests) = reuse::scope(snapshot, || unit::process(ctx, node, dec));
    done.map(|done| Processed { done, requests })
}

/// Proves one level's pairs, in order.
fn run_level(
    ctx: &Context,
    pairs: &[Node],
    snapshot: Option<&reuse::Snapshot>,
    dec: &mut Decider,
) -> Result<Vec<Processed>, Refusal> {
    scoped_try_collect(
        pairs.len(),
        |index, dec| process(ctx, &pairs[index], snapshot, dec).map_err(|error| (index, error)),
        dec,
    )
    .map_err(|(index, error)| named(error, &pairs[index]))
}
