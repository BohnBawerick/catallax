#![allow(missing_docs)]
mod common;
use catallax_numbers::*;
use common::*;

fn fact(kind: FactKind, a: u32, b: u32, at: Option<u32>) -> Fact<u32> {
    Fact {
        kind,
        a,
        b,
        at,
        origin: FactOrigin::Construction("test fixture"),
    }
}

blind_test!(fact_symmetry_transitivity_and_direction, {
    for kind in [
        FactKind::Coincident,
        FactKind::Coaxial,
        FactKind::Concentric,
        FactKind::Parallel,
    ] {
        let mut book = FactBook::new();
        assert!(book.holds(kind, 99, 99));
        assert!(book.record(fact(kind, 1, 2, None)));
        assert!(book.holds(kind, 2, 1));
        assert!(book.record(fact(kind, 2, 3, None)));
        assert!(book.holds(kind, 1, 3));
        assert!(!book.holds(kind, 1, 4));
    }
    let mut book = FactBook::new();
    book.record(fact(FactKind::Tangent, 1, 2, None));
    book.record(fact(FactKind::Tangent, 2, 3, None));
    assert!(book.holds(FactKind::Tangent, 2, 1));
    assert!(!book.holds(FactKind::Tangent, 1, 3));
    assert!(!book.holds(FactKind::Tangent, 1, 1));
    book.record(fact(FactKind::Incident, 10, 11, None));
    book.record(fact(FactKind::Incident, 11, 12, None));
    book.record(fact(FactKind::Coincident, 12, 13, None));
    book.record(fact(FactKind::Coincident, 9, 10, None));
    assert!(book.holds(FactKind::Incident, 10, 12));
    assert!(book.holds(FactKind::Incident, 9, 13));
    assert!(!book.holds(FactKind::Incident, 12, 10));
    assert!(!book.holds(FactKind::Incident, 99, 99));
    book.record(fact(FactKind::Incident, 12, 10, None));
    assert!(book.holds(FactKind::Incident, 10, 10));
    for kind in [
        FactKind::Coaxial,
        FactKind::Concentric,
        FactKind::Parallel,
        FactKind::Tangent,
    ] {
        let mut b = FactBook::new();
        b.record(fact(FactKind::Coincident, 1, 2, None));
        b.record(fact(kind, 2, 3, None));
        assert!(!b.holds(kind, 1, 3));
    }
});

blind_test!(
    facts_deduplicate_ignore_origin_and_preserve_recording_order,
    {
        let mut book = FactBook::new();
        assert!(book.is_empty());
        let first = fact(FactKind::Tangent, 1, 2, Some(8));
        assert!(book.record(first));
        assert!(!book.record(Fact {
            a: 2,
            b: 1,
            origin: FactOrigin::Proof(ZeroProof::Algebraic),
            ..first
        }));
        let second = fact(FactKind::Tangent, 1, 2, Some(9));
        assert!(book.record(second));
        let third = fact(FactKind::Incident, 2, 1, None);
        assert!(book.record(third));
        let fourth = fact(FactKind::Incident, 1, 2, None);
        assert!(book.record(fourth));
        assert!(!book.record(fact(FactKind::Coaxial, 1, 1, None)));
        assert_eq!(book.len(), 4);
        assert_eq!(
            book.iter().copied().collect::<Vec<_>>(),
            vec![first, second, third, fourth]
        );
        assert_eq!(book.facts_about(8), vec![&first]);
        assert_eq!(book.facts_about(2), vec![&first, &second, &third, &fourth]);
        assert!(book.facts_about(99).is_empty());
    }
);

blind_test!(hundred_thousand_facts_and_thousand_link_chains, {
    let mut book = FactBook::new();
    for i in 0..100_000 {
        assert!(book.record(fact(FactKind::Coaxial, i, i + 1, None)));
    }
    assert_eq!(book.len(), 100_000);
    for start in (0..=99_000).step_by(1000) {
        assert!(book.holds(FactKind::Coaxial, start, start + 1000));
        assert!(book.holds(FactKind::Coaxial, start + 1000, start));
    }
    assert!(!book.holds(FactKind::Coaxial, 0, 100_001));
});

blind_test!(
    forked_fact_closures_keep_independent_unions_and_incidence,
    {
        let mut original = FactBook::new();
        for i in 0..300 {
            original.record(fact(FactKind::Incident, i, i + 1, None));
        }
        let recorded = original.iter().copied().collect::<Vec<_>>();
        let mut branch = original.clone();
        branch.record(fact(FactKind::Coincident, 0, 300, None));
        branch.record(fact(FactKind::Tangent, 4, 9, Some(20)));
        assert!(branch.holds(FactKind::Incident, 300, 1));
        assert!(!original.holds(FactKind::Incident, 300, 1));
        assert!(!original.holds(FactKind::Tangent, 4, 9));
        original.record(fact(FactKind::Parallel, 5, 7, None));
        assert!(!branch.holds(FactKind::Parallel, 5, 7));
        assert!(branch.holds(FactKind::Tangent, 4, 9));
        assert_eq!(
            original.iter().take(300).copied().collect::<Vec<_>>(),
            recorded
        );
        assert_eq!(
            branch.iter().take(300).copied().collect::<Vec<_>>(),
            recorded
        );
    }
);

blind_test!(generated_fact_equivalence_matches_graph_reachability, {
    let mut rng = Rng(0xfac7);
    let mut book = FactBook::new();
    let mut edges = [[false; 32]; 32];
    for (i, row) in edges.iter_mut().enumerate() {
        row[i] = true;
    }
    for _ in 0..128 {
        let a = (rng.next() % 32) as usize;
        let b = (rng.next() % 32) as usize;
        book.record(fact(FactKind::Parallel, a as u32, b as u32, None));
        edges[a][b] = true;
        edges[b][a] = true;
    }
    // Independent transitive closure, no FactBook answers in the oracle.
    for k in 0..32 {
        for i in 0..32 {
            for j in 0..32 {
                edges[i][j] |= edges[i][k] && edges[k][j];
            }
        }
    }
    for (i, row) in edges.iter().enumerate() {
        for (j, &expected) in row.iter().enumerate() {
            assert_eq!(book.holds(FactKind::Parallel, i as u32, j as u32), expected);
        }
    }
});

blind_test!(
    incidence_quotient_matches_mixed_directed_graph_after_every_union,
    {
        use std::collections::{BTreeSet, VecDeque};
        let mut rng = Rng(0xc01ac1de);
        let mut book = FactBook::new();
        let mut recorded = Vec::new();
        for turn in 0..96 {
            let a = (rng.next() % 16) as u32;
            let b = (rng.next() % 16) as u32;
            let kind = if turn % 3 == 0 {
                FactKind::Coincident
            } else {
                FactKind::Incident
            };
            let next = fact(kind, a, b, (turn % 4 == 0).then_some(100));
            if book.record(next) {
                recorded.push(next);
            }
            for a in 0..17 {
                for b in 0..17 {
                    let mut seen = BTreeSet::from([(a, false)]);
                    let mut queue = VecDeque::from([(a, false)]);
                    let mut expected = false;
                    while let Some((key, used)) = queue.pop_front() {
                        for fact in &recorded {
                            let next = if fact.a == key {
                                Some(fact.b)
                            } else if fact.kind == FactKind::Coincident && fact.b == key {
                                Some(fact.a)
                            } else {
                                None
                            };
                            if let Some(next) = next {
                                let used = used || fact.kind == FactKind::Incident;
                                expected |= next == b && used;
                                if seen.insert((next, used)) {
                                    queue.push_back((next, used));
                                }
                            }
                        }
                    }
                    assert_eq!(
                        book.holds(FactKind::Incident, a, b),
                        expected,
                        "insertion {turn}, {a} -> {b}"
                    );
                }
            }
        }
    }
);

blind_test!(marks_worst_by_kind_order_and_answer_mapping, {
    let mut first = MarkLog::new();
    assert!(first.is_empty());
    assert!(first.worst(BoundKind::Distance).is_none());
    let ne = NotExact {
        bound: r("1/10"),
        code: "numbers.precision-limit",
        detail: "oracle case".into(),
    };
    let value = Mark::from_not_exact(&ne, BoundKind::Value, ne.bound.clone(), "unconverted");
    let distance = Mark::from_not_exact(&ne, BoundKind::Distance, r("1/5"), "converted by caller");
    let angle = Mark::new(
        "numbers.budget-exhausted",
        BoundKind::Angle,
        r("100"),
        "angle",
    );
    assert_eq!(distance.code, ne.code);
    assert_eq!(distance.bound, r("1/5"));
    first.push(value.clone());
    first.push(distance.clone());
    let mut second = MarkLog::new();
    second.push(angle.clone());
    second.push(Mark::new(
        "numbers.precision-limit",
        BoundKind::Distance,
        r("1/3"),
        "larger distance",
    ));
    first.extend(&second);
    assert_eq!(first.len(), 4);
    assert_eq!(first.worst(BoundKind::Distance), Some(&r("1/3")));
    assert_eq!(first.worst(BoundKind::Angle), Some(&r("100")));
    assert_eq!(first.worst(BoundKind::Value), Some(&r("1/10")));
    assert_eq!(
        first.iter().take(3).cloned().collect::<Vec<_>>(),
        vec![value, distance, angle]
    );
    assert_eq!(second.len(), 2);
    let exact = Answer::exact(7);
    assert!(exact.is_exact());
    assert_eq!(exact.map(|x| x * 2).value, 14);
    let marked = Answer::with_marks(7, first.clone()).map(|x| x.to_string());
    assert!(!marked.is_exact());
    assert_eq!(marked.value, "7");
    assert_eq!(marked.marks, first);
});

blind_test!(outcome_mapping_preserves_uncertainty_and_refusals, {
    let variants: Vec<Outcome<i32>> = vec![
        Outcome::Proven(3),
        Outcome::ProvenCoincident(Coincidence {
            proof: ZeroProof::Algebraic,
        }),
        Outcome::NotExact(NotExact {
            bound: r("1/7"),
            code: "numbers.precision-limit",
            detail: "bounded".into(),
        }),
        Outcome::Refused(Refusal::new("numbers.deadline", "expired")),
    ];
    for outcome in variants {
        let mut called = false;
        let mapped = outcome.clone().map(|v| {
            called = true;
            v + 1
        });
        assert_eq!(called, outcome.is_proven());
        match outcome {
            Outcome::Proven(v) => assert_eq!(mapped.proven(), Some(&(v + 1))),
            Outcome::ProvenCoincident(c) => {
                assert_eq!(mapped, Outcome::ProvenCoincident(c));
                assert!(mapped.is_coincident());
            }
            Outcome::NotExact(ne) => {
                assert_eq!(mapped, Outcome::NotExact(ne));
                assert!(mapped.is_not_exact());
            }
            Outcome::Refused(e) => {
                assert_eq!(mapped, Outcome::Refused(e));
                assert!(mapped.is_refused());
            }
        }
    }
});
