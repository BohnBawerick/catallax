//! Certified trace continuity, orientation, Taylor identities and arc lengths.
mod common;
use common::curves::*;
use common::*;

fn ball(q: &Rational, bits: u32) -> Ball {
    Ball::from_rational(q, bits)
}
fn overlap(a: &Ball, b: &Ball) {
    assert!(
        a.lo() <= b.hi() && b.lo() <= a.hi(),
        "disjoint enclosures {a:?}, {b:?}"
    );
}
fn point_overlap(actual: &[Ball; 3], expected: &Point3) {
    for (a, e) in actual.iter().zip(xyz(expected)) {
        if let Some(q) = e.as_rational() {
            assert!(a.contains(&q));
        } else {
            let bound = Decider::default().enclose(&e, &r("1e-40")).unwrap();
            overlap(a, &bound);
        }
    }
}
fn dot(a: &[Ball; 3], b: &[Ball; 3]) -> Ball {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}
fn tangent_bound(t: &Vector3, bits: u32) -> [Ball; 3] {
    [&t.x, &t.y, &t.z].map(|x| {
        let b = Decider::default()
            .enclose(x, &n(2).pow(-(bits as i32)).unwrap())
            .unwrap();
        b.with_precision(bits)
    })
}
fn polynomial_jet_identity(equation: &MPoly, jets: &[[Ball; 3]], bits: u32) {
    let size = jets.len();
    let zero = ball(&n(0), bits);
    let one = ball(&n(1), bits);
    let mut result = vec![zero.clone(); size];
    for (monomial, coefficient) in equation.terms() {
        let mut product = vec![zero.clone(); size];
        product[0] = one.clone();
        for (axis, exponent) in monomial.exponents().iter().enumerate() {
            for _ in 0..*exponent {
                let mut next = vec![zero.clone(); size];
                for k in 0..size {
                    for j in 0..=k {
                        next[k] = &next[k] + &(&product[j] * &jets[k - j][axis]);
                    }
                }
                product = next;
            }
        }
        let coefficient = ball(&coefficient, bits);
        for k in 0..size {
            result[k] = &result[k] + &(&coefficient * &product[k]);
        }
    }
    // f(C(t))=0 identically, so every normalized derivative must enclose zero.
    // This catches higher-order coefficients incorrectly returned as zero.
    for coefficient in result {
        assert!(coefficient.contains_zero());
    }
}
fn check_trace(
    g: &GeomStore,
    c: CurveId,
    from: Option<PointId>,
    to: Option<PointId>,
    start: &Point3,
    end: &Point3,
    target: Ball,
) -> String {
    let before = fingerprint(g);
    let mut d = Decider::default();
    let trace = exact(trace(g, c, from, to, &mut d));
    assert_eq!(fingerprint(g), before);
    assert!(!trace.pieces().is_empty());
    let pieces = trace.pieces();
    for bits in [64, 1000] {
        let first = &pieces[0];
        let last = pieces.last().unwrap();
        point_overlap(
            &first
                .evaluate(&ball(&first.range().0, bits), &Decider::default())
                .unwrap(),
            start,
        );
        point_overlap(
            &last
                .evaluate(&ball(&last.range().1, bits), &Decider::default())
                .unwrap(),
            end,
        );
        if from.is_some() && from == to {
            let witness = g.curve(c).parameter_origin(g);
            let witness =
                xyz(&witness).map(|x| Decider::default().enclose(&x, &r("1e-40")).unwrap());
            assert!(
                pieces.iter().any(|piece| {
                    let (a, b) = piece.range();
                    let image = piece
                        .evaluate(&ball(&a, bits).hull(&ball(&b, bits)), &Decider::default())
                        .unwrap();
                    image
                        .iter()
                        .zip(&witness)
                        .all(|(a, b)| a.lo() <= b.hi() && b.lo() <= a.hi())
                }),
                "whole loop must pass through its witness"
            );
        }
        for (i, piece) in pieces.iter().enumerate() {
            let (a, b) = piece.range();
            assert!(a < b);
            if i > 0 {
                let previous = &pieces[i - 1];
                let lhs = previous
                    .evaluate(&ball(&previous.range().1, bits), &Decider::default())
                    .unwrap();
                let rhs = piece
                    .evaluate(&ball(&a, bits), &Decider::default())
                    .unwrap();
                for (l, r) in lhs.iter().zip(rhs) {
                    overlap(l, &r);
                }
            }
            // Rational sub-balls, including both ends and interior intervals.
            let mut rng = Rng(0x044043 + i as u64);
            for sample in 0..10 {
                let (lo, hi) = match sample {
                    0 => (a.clone(), a.clone()),
                    1 => (b.clone(), b.clone()),
                    _ => {
                        let j = (rng.next() % 64) as i64;
                        let width = (&b - &a) / n(64);
                        (&a + &(&width * &n(j)), &a + &(width * n(j + 1)))
                    }
                };
                let t = ball(&lo, bits).hull(&ball(&hi, bits));
                let p = piece.evaluate(&t, &Decider::default()).unwrap();
                let v = piece.derivative(&t, &Decider::default()).unwrap();
                for x in p.iter().chain(v.iter()) {
                    assert_eq!(x.precision(), bits);
                }
                for equation in g.curve(c).equations(g) {
                    assert!(equation.eval_ball(&p).contains_zero());
                }
                for order in [4, 16] {
                    let jets = piece.taylor(&t, order, &Decider::default()).unwrap();
                    assert_eq!(jets.len(), order + 1);
                    for j in &jets {
                        for x in j {
                            assert_eq!(x.precision(), bits);
                        }
                    }
                    for axis in 0..3 {
                        overlap(&jets[0][axis], &p[axis]);
                        overlap(&jets[1][axis], &v[axis]);
                    }
                    for equation in g.curve(c).equations(g) {
                        polynomial_jet_identity(&equation, &jets, bits);
                    }
                    // Taylor's theorem provides another point enclosure from midpoint
                    // jets and an interval remainder, independent of evaluate(t).
                    let mid = ball(&((&lo + &hi) / n(2)), bits);
                    let centre = piece.taylor(&mid, order, &Decider::default()).unwrap();
                    let delta = &t - &mid;
                    for axis in 0..3 {
                        let mut enclosure = jets[order][axis].clone();
                        for k in (0..order).rev() {
                            enclosure = &(&enclosure * &delta) + &centre[k][axis];
                        }
                        overlap(&enclosure, &p[axis]);
                    }
                }
            }
        }
        // At externally specified arc endpoints we have exact PointIds. Internal
        // joins expose only balls, so their continuity is checked by enclosure.
        for (id, piece, t, which) in [
            (from, first, first.range().0, 1usize),
            (to, last, last.range().1, 0usize),
        ] {
            if let Some(id) = id {
                let tangent = tangent_at(g, c, id, &mut Decider::default()).unwrap();
                let velocity = piece
                    .derivative(&ball(&t, bits), &Decider::default())
                    .unwrap();
                assert!(
                    dot(&velocity, &tangent_bound(&tangent[which], bits)).lo() > n(0),
                    "wrong endpoint direction"
                );
            }
        }
    }
    // Validated interval quadrature: taylor(t,1)[1] encloses C' throughout
    // each subinterval. Its norm times the interval width bounds that length.
    // No parametrisation or subdivision count is prescribed by the contract.
    let bits = 128;
    let mut length = ball(&n(0), bits);
    for piece in pieces {
        let (a, b) = piece.range();
        let step = (&b - &a) / n(128);
        for j in 0..128 {
            let lo = &a + &(&step * &n(j));
            let hi = &a + &(&step * &n(j + 1));
            let t = ball(&lo, bits).hull(&ball(&hi, bits));
            let jets = piece.taylor(&t, 1, &Decider::default()).unwrap();
            let v = &jets[1];
            let speed = (v[0].square() + v[1].square() + v[2].square())
                .sqrt()
                .unwrap();
            length = length + speed * ball(&step, bits);
        }
    }
    overlap(&length, &target);
    eprintln!("trace length enclosure [{}, {}]", length.lo(), length.hi());
    format!("{trace:?}; {:?}", d.stats())
}
blind_test!(cases_40_41_circle_and_line_traces, {
    let mut g = GeomStore::new();
    let k = g.add_curve(k5());
    let l = g.add_curve(line28());
    let p = explicit(&mut g, pi([5, 0, 0]));
    let q = explicit(&mut g, pi([0, 5, 0]));
    let a = explicit(&mut g, pi([-5, 0, 0]));
    for (c, from, to, start, end, length) in [
        (
            k,
            None,
            None,
            pi([5, 0, 0]).to_real(),
            pi([5, 0, 0]).to_real(),
            Ball::pi(128) * ball(&n(10), 128),
        ),
        (
            k,
            Some(p),
            Some(q),
            g.position(p),
            g.position(q),
            Ball::pi(128) * ball(&r("5/2"), 128),
        ),
        (
            l,
            Some(a),
            Some(p),
            g.position(a),
            g.position(p),
            ball(&n(10), 128),
        ),
    ] {
        let first = check_trace(&g, c, from, to, &start, &end, length.clone());
        let second = check_trace(&g, c, from, to, &start, &end, length);
        assert_eq!(first, second);
    }
});
blind_test!(cases_42_43_quartic_loops_and_singular_endpoints, {
    let mut g = GeomStore::new();
    let c25 = loop25().add(&mut g);
    let c24 = viviani(true).add(&mut g);
    let node = explicit(&mut g, pi([-1, 0, 0]));
    // Independent rational-parameter integrands and numerical cross-checks are
    // in verify_math.py. Twelve decimal places suffice for these enclosure checks.
    let l24 = ball(&r("7.640395578055"), 128).hull(&ball(&r("7.640395578056"), 128));
    let l25 = ball(&r("9.538833617686"), 128).hull(&ball(&r("9.538833617687"), 128));
    let w25 = Point3::new(rr(1), rr(0), rr(3).sqrt());
    let w24 = pi([1, 0, 2]).to_real();
    check_trace(&g, c25, None, None, &w25, &w25, l25);
    check_trace(&g, c24, None, None, &w24, &w24, l24.clone());
    check_trace(
        &g,
        c24,
        Some(node),
        Some(node),
        &g.position(node),
        &g.position(node),
        l24,
    );
});
blind_test!(trace_budget_failure_has_no_marks_or_mutation, {
    let mut g = GeomStore::new();
    let c = g.add_curve(k5());
    let p = explicit(&mut g, pi([5, 0, 0]));
    let q = explicit(&mut g, pi([0, 5, 0]));
    for mut d in [tiny(), starved(), expired()] {
        let before = fingerprint(&g);
        match trace(&g, c, Some(p), Some(q), &mut d) {
            Ok(a) => {
                assert!(a.marks.is_empty());
                let first = a.value.pieces().first().unwrap();
                let last = a.value.pieces().last().unwrap();
                point_overlap(
                    &first
                        .evaluate(&ball(&first.range().0, 128), &Decider::default())
                        .unwrap(),
                    &g.position(p),
                );
                point_overlap(
                    &last
                        .evaluate(&ball(&last.range().1, 128), &Decider::default())
                        .unwrap(),
                    &g.position(q),
                );
            }
            Err(e) => assert!(
                [
                    "numbers.deadline",
                    "numbers.budget-exhausted",
                    "numbers.precision-limit",
                    "numbers.divisor-not-separated",
                    "numbers.root-not-certified"
                ]
                .contains(&e.code)
            ),
        }
        assert_eq!(fingerprint(&g), before);
    }
});
