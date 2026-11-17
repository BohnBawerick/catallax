//! Contact setup for OCCT-style Cliff intersections. The actual fitted rails
//! remain authoritative; a failed finite-domain solve never extrapolates them.

use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_topology::{Model, VertexId};

use super::super::super::rolling::RolledChain;
use super::{BITS, StripEnd, bad, gap, strip};

pub(crate) fn common_contact(
    model: &Model,
    chains: &[RolledChain],
    vertex: VertexId,
    ends: [StripEnd; 2],
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    if tolerance.sign() != Sign::Positive {
        return Err(bad("Cliff contact tolerance must be positive"));
    }
    let sources = [
        strip(model, chains, vertex, ends[0])?,
        strip(model, chains, vertex, ends[1])?,
    ];
    let common = sources[0]
        .source
        .faces
        .iter()
        .enumerate()
        .find_map(|(i, face)| {
            sources[1]
                .source
                .faces
                .iter()
                .position(|other| other == face)
                .map(|j| [i, j])
        })
        .ok_or_else(|| bad("Cliff strips have no common source support"))?;
    let rails = [
        sources[0]
            .fillet
            .rail(common[0])
            .ok_or_else(|| bad("first Cliff contact rail is missing"))?,
        sources[1]
            .fillet
            .rail(common[1])
            .ok_or_else(|| bad("second Cliff contact rail is missing"))?,
    ];
    let domains = rails.map(|rail| rail.domain().map(Rational::to_f64_nearest));
    if domains.iter().flatten().any(|value| !value.is_finite()) {
        return Err(bad("Cliff rail domain is outside binary64 proposal range"));
    }
    let mut parameters = ends.map(|end| end.source_end as f64);
    let limit = tolerance.mul_pow2(-3);
    for _ in 0..32 {
        dec.charge_work(1)?;
        let at = |side: usize| -> Result<Rational, Refusal> {
            let q = Rational::from_f64_exact(parameters[side])
                .ok_or_else(|| bad("non-finite Cliff contact proposal"))?;
            let domain = rails[side].domain();
            Ok(q.max(domain[0].clone()).min(domain[1].clone()))
        };
        let tau = [at(0)?, at(1)?];
        let jets = [
            rails[0].jet_box([&tau[0], &tau[0]], BITS, dec)?,
            rails[1].jet_box([&tau[1], &tau[1]], BITS, dec)?,
        ];
        let error = gap(&jets[0].value, &jets[1].value)?;
        if error <= limit {
            return Ok(tau);
        }
        let residual: [f64; 3] = std::array::from_fn(|i| {
            jets[0].value[i].mid().to_f64_nearest() - jets[1].value[i].mid().to_f64_nearest()
        });
        let a = jets[0]
            .first
            .each_ref()
            .map(|value| value.mid().to_f64_nearest());
        let b = jets[1]
            .first
            .each_ref()
            .map(|value| value.mid().to_f64_nearest());
        let dot = |a: &[f64; 3], b: &[f64; 3]| (0..3).map(|i| a[i] * b[i]).sum::<f64>();
        let aa = dot(&a, &a);
        let bb = dot(&b, &b);
        let ab = -dot(&a, &b);
        let determinant = aa * bb - ab * ab;
        if !determinant.is_finite() || determinant <= 0.0 {
            return Err(bad(
                "Cliff contact rails have an unresolved intersection direction",
            ));
        }
        let rhs = [-dot(&a, &residual), dot(&b, &residual)];
        let step = [
            (bb * rhs[0] - ab * rhs[1]) / determinant,
            (aa * rhs[1] - ab * rhs[0]) / determinant,
        ];
        if step.iter().any(|value| !value.is_finite()) {
            return Err(bad("Cliff contact solve produced a non-finite step"));
        }
        let next: [f64; 2] =
            std::array::from_fn(|i| (parameters[i] + step[i]).clamp(domains[i][0], domains[i][1]));
        if (0..2).all(|i| next[i].to_bits() == parameters[i].to_bits()) {
            return Err(bad(format!(
                "common rail contact stalled in finite strip domains: tau={parameters:?}, domains={domains:?}, world separation={} mm, points={:?}",
                error.to_f64_nearest(),
                jets.each_ref().map(|jet| jet
                    .value
                    .each_ref()
                    .map(|value| value.mid().to_f64_nearest())),
            )));
        }
        parameters = next;
    }
    Err(bad(format!(
        "common rail contact did not converge in finite strip domains: tau={parameters:?}, domains={domains:?}"
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{BudgetPolicy, with_operation_budget};
    use catallax_recipe::{Recipe, Role};

    #[test]
    fn owner_vertex_cylindrical_cliff_common_contact_is_retained() {
        let ((), work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::new(BudgetPolicy {
                max_work: 100_000_000,
                ..Default::default()
            });
            let recipe =
                Recipe::parse(include_str!("../tests/fixtures/owner-curved-pocket.json")).unwrap();
            let built = crate::evaluate(&recipe, &mut dec).unwrap();
            let model = &built.model;
            let selected: Vec<_> = model.edges_of_body(built.result).into_iter().filter(|edge| {
                model.edge_name(built.result, *edge).is_some_and(|name| {
                    name.between().iter().any(|origin| origin.step.as_str() == "blank"
                        && matches!(&origin.role, Role::Side { segment, .. } if segment.as_str() == "right-round"))
                        && name.between().iter().any(|origin| origin.step.as_str() == "pocket-tool"
                            && matches!(&origin.role, Role::Side { segment, .. } if ["far", "right"].contains(&segment.as_str())))
                })
            }).collect();
            assert_eq!(selected.len(), 2);
            let tolerance = Rational::parse("0.0000001").unwrap();
            let chains = super::super::super::super::rolling::construct(
                model,
                built.result,
                &selected,
                &Rational::one(),
                &tolerance,
                &mut dec,
            )
            .unwrap();
            assert_eq!(chains.len(), 2);
            let a = model.edge(chains[0].strips[0].source.edge);
            let b = model.edge(chains[1].strips[0].source.edge);
            let vertex = [a.start, a.end]
                .into_iter()
                .flatten()
                .find(|vertex| b.start == Some(*vertex) || b.end == Some(*vertex))
                .unwrap();
            let ends = [a, b].map(|edge| usize::from(edge.start != Some(vertex)));
            let ends = std::array::from_fn(|i| StripEnd {
                chain: i,
                strip: 0,
                source_end: ends[i],
            });
            let tau = common_contact(model, &chains, vertex, ends, &tolerance, &mut dec).unwrap();
            for side in 0..2 {
                let rail = chains[side].strips[0].fillet.rail(0).unwrap();
                let point = rail
                    .evaluate_box([&tau[side], &tau[side]], BITS, &dec)
                    .unwrap();
                for (axis, expected) in [(0, "33.75"), (1, "-25.665")] {
                    let expected = Rational::parse(expected).unwrap();
                    assert!(
                        (&point[axis].lo() - &expected)
                            .abs()
                            .max((&point[axis].hi() - &expected).abs())
                            <= tolerance
                    );
                }
            }
        });
        assert!(!work.exhausted);
        println!(
            "owner cylindrical Cliff common-contact work={}; actual retained rails, no extrapolation",
            work.used
        );
    }
}
