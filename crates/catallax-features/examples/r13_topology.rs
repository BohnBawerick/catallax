//! R13 exact, bounded applicability controls; not either published SSI algorithm.
use catallax_numbers::Rational as Q;

type Poly = Vec<Q>;
fn q(n: i64) -> Q {
    Q::from_i64(n)
}
fn r(s: &str) -> Q {
    Q::parse(s).unwrap()
}
fn add(a: &Poly, b: &Poly) -> Poly {
    (0..a.len().max(b.len()))
        .map(|i| a.get(i).cloned().unwrap_or_default() + b.get(i).cloned().unwrap_or_default())
        .collect()
}
fn scale(a: &Poly, c: &Q) -> Poly {
    a.iter().map(|v| v * c).collect()
}
fn mul(a: &Poly, b: &Poly) -> Poly {
    let mut out = vec![q(0); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] = &out[i + j] + &(x * y);
        }
    }
    out
}
fn eval(p: &Poly, x: &Q) -> Q {
    p.iter().rev().fold(q(0), |v, c| &(&v * x) + c)
}
fn derivative(p: &Poly) -> Poly {
    p.iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| c * &q(i as i64))
        .collect()
}
fn choose(n: usize, k: usize) -> i64 {
    (0..k).fold(1, |a, i| a * (n - i) as i64 / (i + 1) as i64)
}
// Cox-de Boor on one nonempty knot span. Exact polynomial identity, not fitting.
fn basis(i: usize, degree: usize, span: usize, knots: &[i64]) -> Poly {
    if degree == 0 {
        return vec![q(i64::from(i == span))];
    }
    let mut result = vec![q(0)];
    if knots[i + degree] != knots[i] {
        let a = scale(
            &vec![q(-knots[i]), q(1)],
            &(&q(1) / &q(knots[i + degree] - knots[i])),
        );
        result = add(&result, &mul(&a, &basis(i, degree - 1, span, knots)));
    }
    if knots[i + degree + 1] != knots[i + 1] {
        let a = scale(
            &vec![q(knots[i + degree + 1]), q(-1)],
            &(&q(1) / &q(knots[i + degree + 1] - knots[i + 1])),
        );
        result = add(&result, &mul(&a, &basis(i + 1, degree - 1, span, knots)));
    }
    result
}
// Power polynomial on [a,b] -> Bernstein coefficients. Their convex hull
// encloses the entire polynomial, including between sample locations.
fn bernstein(p: &Poly, a: &Q, b: &Q) -> Poly {
    let mut shifted = vec![q(0)];
    let mut power = vec![q(1)];
    for c in p {
        shifted = add(&shifted, &scale(&power, c));
        power = mul(&power, &vec![a.clone(), b - a]);
    }
    let n = p.len() - 1;
    (0..=n)
        .map(|k| {
            (0..=k).fold(q(0), |v, j| {
                v + &shifted[j] * &(&q(choose(k, j)) / &q(choose(n, j)))
            })
        })
        .collect()
}
fn positive_inside(p: &Poly, a: &Q, b: &Q) -> bool {
    let c = bernstein(p, a, b);
    c.iter().all(|v| *v >= q(0)) && c.iter().any(|v| *v > q(0))
}
#[derive(Clone)]
struct Hill {
    height: Vec<Poly>,
    x: Vec<Poly>,
}
impl Hill {
    fn new() -> Self {
        // Exact public controls: patch_fixtures::height_patch(true), unit weights,
        // degree (3,3), knots [0,0,0,0,1,2,3,3,3,3], z_ij=1 for i,j in {2,3}.
        let knots = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3];
        let height = (3..6)
            .map(|s| add(&basis(2, 3, s, &knots), &basis(3, 3, s, &knots)))
            .collect();
        let x = (3..6)
            .map(|s| {
                (0..6).fold(vec![q(0)], |p, i| {
                    add(&p, &scale(&basis(i, 3, s, &knots), &q(i as i64)))
                })
            })
            .collect();
        Self { height, x }
    }
    fn h(&self, t: &Q) -> Q {
        let span = if *t < q(1) {
            0
        } else if *t < q(2) {
            1
        } else {
            2
        };
        eval(&self.height[span], t)
    }
    fn certify_shape(&self) {
        assert_eq!(self.h(&q(0)), q(0));
        assert_eq!(self.h(&q(3)), q(0));
        assert_eq!(self.h(&r("3/2")), r("15/16"));
        assert_eq!(eval(&self.x[1], &r("3/2")), r("5/2"));
        for join in [1, 2] {
            for pieces in [&self.height, &self.x] {
                assert_eq!(
                    eval(&pieces[join - 1], &q(join as i64)),
                    eval(&pieces[join], &q(join as i64))
                );
            }
        }
        for (s, p) in self.height.iter().enumerate() {
            assert!(positive_inside(p, &q(s as i64), &q(s as i64 + 1)));
            assert!(positive_inside(
                &derivative(&self.x[s]),
                &q(s as i64),
                &q(s as i64 + 1)
            ));
        }
        for (s, a, b, sign) in [
            (0, q(0), q(1), 1),
            (1, q(1), r("3/2"), 1),
            (1, r("3/2"), q(2), -1),
            (2, q(2), q(3), -1),
        ] {
            assert!(positive_inside(
                &scale(&derivative(&self.height[s]), &q(sign)),
                &a,
                &b
            ));
        }
        // Exact signs on the whole central square's four edges establish an
        // inward gradient homotopy to the radial field, hence index +1.
        assert!(eval(&derivative(&self.height[1]), &q(1)) > q(0));
        assert!(eval(&derivative(&self.height[1]), &q(2)) < q(0));
        assert_eq!(eval(&derivative(&self.height[1]), &r("3/2")), q(0));
    }
}
#[derive(Debug, PartialEq, Eq)]
enum Topology {
    OneInteriorLoop,
    IsolatedContact,
    Disjoint,
    Unresolved,
}
fn level(h: &Hill, c: &Q) -> Topology {
    h.certify_shape();
    let max = &h.h(&r("3/2")) * &h.h(&r("3/2"));
    if *c > max || *c < q(0) {
        Topology::Disjoint
    } else if *c == max {
        Topology::IsolatedContact
    } else if *c > q(0) {
        Topology::OneInteriorLoop
    } else {
        Topology::Unresolved
    } // boundary level, not covered by the interior argument
}
fn experiment() -> Vec<String> {
    let h = Hill::new();
    h.certify_shape();
    let mut rows = Vec::new();
    for (name, c, want) in [
        ("public-hill", r("3/10"), Topology::OneInteriorLoop),
        (
            "public-tiny",
            r("225/256") - r("1/1000000"),
            Topology::OneInteriorLoop,
        ),
        ("public-tangent", r("225/256"), Topology::IsolatedContact),
        (
            "public-above",
            r("225/256") + r("1/1000000"),
            Topology::Disjoint,
        ),
        ("boundary-control", q(0), Topology::Unresolved),
    ] {
        let actual = level(&h, &c);
        assert_eq!(actual, want);
        rows.push(format!(
            "{name}: normal-rule=unresolved gradient-index=+1 topology={actual:?}"
        ));
    }
    let mut sampled_max = q(0);
    for i in 0..=17 {
        for j in 0..=17 {
            sampled_max = sampled_max.max(h.h(&(&q(3 * i) / &q(17))) * h.h(&(&q(3 * j) / &q(17))));
        }
    }
    assert!(sampled_max < r("225/256") - r("1/1000000"));
    rows.push(format!("sampling-counterexample: 324 exact grid samples below tiny level; sampled-max={sampled_max}; loop still exists"));
    // Original control F=x, graph normal (-1,0,1) vs plane (0,0,1):
    // constant, nonparallel normals; a boundary-to-boundary intersection x=0.
    // No closed loop is not a disjointness verdict.
    rows.push(
        "linear-control: normal-rule=no-closed-loop topology=open-line gradient-index=0".into(),
    );
    // F=(x^2-1)^2+y^2. grad=(4x(x^2-1),2y), roots (-1,0),(0,0),(1,0).
    // Hessian diag(12x^2-4,2). Rectangle [-1/2,3/2] x [-1,1]
    // encloses saddle (0,0), min (1,0): degree -1+1=0, despite two zeros.
    let critical = [q(0), q(1)];
    let mut index = 0;
    for x in critical {
        let x2 = &x * &x;
        assert_eq!(&(&q(4) * &x) * &(&x2 - &q(1)), q(0));
        let det = &(&q(12) * &x2 - q(4)) * &q(2);
        index += if det > q(0) { 1 } else { -1 };
    }
    assert_eq!(index, 0);
    rows.push("index-cancellation-control: index=0 critical-points=2 absence=unresolved".into());
    // F=x^2 has an entire critical line; isolated-zero index machinery alone
    // is inapplicable. This is not an implementation of the 2026 subdivision.
    rows.push("tangent-line-control: F=x^2 critical-set=line isolated-index=inapplicable".into());
    rows
}
fn main() {
    let expected = experiment();
    for workers in [1, 2, 4] {
        std::thread::scope(|s| {
            let handles: Vec<_> = (0..workers).map(|_| s.spawn(experiment)).collect();
            for h in handles {
                assert_eq!(h.join().unwrap(), expected);
            }
        });
    }
    for row in expected {
        println!("{row}");
    }
    println!(
        "workers=1/2/4 identical; fixed controls; no full paper implementation or production proof"
    );
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bernstein_detects_interior_sign_change() {
        // Positive endpoints but negative centre: endpoint-only tests are unsound.
        assert!(!positive_inside(&vec![r("3/4"), q(-4), q(4)], &q(0), &q(1)));
        assert!(positive_inside(&vec![q(0), q(1), q(-1)], &q(0), &q(1)));
    }
    #[test]
    fn all_exact_controls() {
        assert_eq!(experiment().len(), 9);
    }
}
