use super::*;
use std::collections::{HashMap, hash_map::Entry};

impl SpanPolynomials {
    pub(super) fn new(
        homogeneous: &[MPoly; 4],
        rotation: &[[Rational; 3]; 3],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let homogeneous = rotated(homogeneous, rotation, dec)?;
        let mut tangents = Vec::with_capacity(2);
        for axis in 0..2 {
            let dw = homogeneous[3].partial(axis);
            let mut parts = Vec::with_capacity(3);
            for coordinate in &homogeneous[..3] {
                let left = product(&coordinate.partial(axis), &homogeneous[3], dec)?;
                let right = product(coordinate, &dw, dec)?;
                parts.push(difference(&left, &right, dec)?);
            }
            tangents.push(parts.try_into().expect("three components"));
        }
        let tangents: [[MPoly; 3]; 2] = tangents.try_into().expect("two tangents");
        let [tu, tv] = &tangents;
        let mut normal = Vec::with_capacity(3);
        for i in 0..3 {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            let left = product(&tu[j], &tv[k], dec)?;
            let right = product(&tu[k], &tv[j], dec)?;
            normal.push(difference(&left, &right, dec)?);
        }
        let normal: [MPoly; 3] = normal.try_into().expect("three components");
        Ok(Self {
            homogeneous,
            tangents,
            normal,
        })
    }
}

struct Group<'a> {
    homogeneous: &'a [MPoly; 4],
    domain: Domain,
    count: usize,
}

// Singleton polynomials move into their cell without another allocation.
#[allow(clippy::large_enum_variant)]
enum Prepared {
    Single(Option<SpanPolynomials>),
    Shared {
        polynomials: Arc<SpanPolynomials>,
        nets: net::SpanNets,
    },
}

impl Side {
    fn cached(
        source: Source<'_>,
        homogeneous: Vec<[MPoly; 4]>,
        rotation: &[[Rational; 3]; 3],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        // Borrow the retained input arrays as keys: HashMap verifies full exact
        // polynomial equality, without cloned keys or fingerprint authority.
        let mut keys: HashMap<&[MPoly; 4], usize> = HashMap::new();
        let mut groups: Vec<Group<'_>> = Vec::new();
        let mut entries = Vec::with_capacity(homogeneous.len());
        for ((chart, local_span, span), h) in source.spans().zip(&homogeneous) {
            dec.check_deadline()?;
            let group = match keys.entry(h) {
                Entry::Vacant(entry) => {
                    let index = groups.len();
                    groups.push(Group {
                        homogeneous: h,
                        domain: span.domain.clone(),
                        count: 1,
                    });
                    entry.insert(index);
                    index
                }
                Entry::Occupied(entry) => {
                    let index = *entry.get();
                    let group = &mut groups[index];
                    group.count += 1;
                    for axis in 0..2 {
                        if net::less(&span.domain[axis][0], &group.domain[axis][0], dec)? {
                            group.domain[axis][0] = span.domain[axis][0].clone();
                        }
                        if net::less(&group.domain[axis][1], &span.domain[axis][1], dec)? {
                            group.domain[axis][1] = span.domain[axis][1].clone();
                        }
                    }
                    index
                }
            };
            entries.push((chart, local_span, span, group));
        }
        let mut prepared = Vec::with_capacity(groups.len());
        for group in groups {
            let polynomials = SpanPolynomials::new(group.homogeneous, rotation, dec)?;
            prepared.push(if group.count == 1 {
                Prepared::Single(Some(polynomials))
            } else {
                let nets = net::SpanNets::new(&group.domain, &polynomials.homogeneous, dec)?;
                Prepared::Shared {
                    polynomials: Arc::new(polynomials),
                    nets,
                }
            });
        }
        let mut spans = Vec::with_capacity(entries.len());
        for (index, (chart, local_span, span, group)) in entries.into_iter().enumerate() {
            dec.check_deadline()?;
            let (polynomials, root) = match &mut prepared[group] {
                Prepared::Single(value) => {
                    let polynomials = value.take().expect("a singleton group has one chart");
                    let root =
                        Cell::span(index, &span.domain, &polynomials.homogeneous, None, dec)?;
                    (PolynomialStorage::Owned(polynomials), root)
                }
                Prepared::Shared { polynomials, nets } => (
                    PolynomialStorage::Shared(polynomials.clone()),
                    nets.restricted_cell(index, &span.domain, dec)?,
                ),
            };
            spans.push(SideSpan {
                chart,
                local_span,
                physical: source.physical(chart, span),
                polynomials,
                offset: None,
                root: Arc::new(root),
            });
        }
        Ok(Self {
            periodic: source.periodic(),
            distance: None,
            spans,
        })
    }
}

pub(super) fn meeting_cached(
    patches: [Source<'_>; 2],
    dec: &mut Decider,
) -> Result<PatchPairIntersection, Refusal> {
    let [this, other] = patches;
    let budget = catallax_numbers::operation_budget_active();
    let start = if budget {
        catallax_numbers::current_operation_work().used
    } else {
        dec.stats().work
    };
    let spans = [rational_spans(this)?, rational_spans(other)?];
    let rotation = upright(this, &spans[0], dec)?;
    let [first, second] = spans;
    let sides = [
        Side::cached(this, first, &rotation, dec)?,
        Side::cached(other, second, &rotation, dec)?,
    ];
    finish(patches, [None, None], sides, rotation, (budget, start), dec)
}
