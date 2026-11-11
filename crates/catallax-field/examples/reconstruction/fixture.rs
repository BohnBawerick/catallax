//! Original MIT-licensed procedural controls. Never called by the fitting code.
use catallax_numbers::Rational;
use catallax_recipe::Number;
fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
/// A rectangular plate with a faceted through-bore. Exact terminating-decimal
/// circle points follow repeated rotations by (cos,sin)=(0.96,0.28).
pub fn fixture(top_radius: i64) -> Vec<u8> {
    let mut inner = Vec::new();
    let mut outer = Vec::new();
    for quadrant in 0..4 {
        let mut u = [q(1), q(0)];
        for k in 0..6 {
            let rotate = |mut p: [Rational; 2]| {
                for _ in 0..quadrant {
                    p = [-&p[1], p[0].clone()];
                }
                p
            };
            inner.push(rotate(u.clone()));
            let p = if k <= 3 {
                [q(5), q([0, 2, 4, 5][k as usize])]
            } else {
                [q(if k == 4 { 3 } else { 1 }), q(5)]
            };
            outer.push(rotate(p));
            u = [
                (&u[0] * &q(24) - &u[1] * &q(7))
                    .checked_div(&q(25))
                    .unwrap(),
                (&u[0] * &q(7) + &u[1] * &q(24))
                    .checked_div(&q(25))
                    .unwrap(),
            ];
        }
    }
    let n = inner.len();
    let mut points = Vec::new();
    for z in [0, 6] {
        for (ring, radius) in [(&outer, 1), (&inner, if z == 0 { 2 } else { top_radius })] {
            for p in ring {
                points.push([&p[0] * &q(radius), &p[1] * &q(radius), q(z)]);
            }
        }
    }
    let mut triangles = Vec::new();
    let mut quad = |a, b, c, d| {
        triangles.push([a, b, c]);
        triangles.push([a, c, d]);
    };
    for i in 0..n {
        let j = (i + 1) % n;
        quad(i, j, 2 * n + j, 2 * n + i); // outer wall, outward
        quad(n + j, n + i, 3 * n + i, 3 * n + j); // bore wall, inward
        quad(i, n + i, n + j, j); // lower cap, downward
        quad(2 * n + j, 3 * n + j, 3 * n + i, 2 * n + i); // upper cap, upward
    }
    let mut text = String::from("# Original Catallax M23 control, MIT; units mm\n");
    for p in points {
        text.push_str(&format!(
            "v {} {} {}\n",
            Number::new(p[0].clone()),
            Number::new(p[1].clone()),
            Number::new(p[2].clone())
        ));
    }
    for t in triangles {
        text.push_str(&format!("f {} {} {}\n", t[0] + 1, t[1] + 1, t[2] + 1));
    }
    text.into_bytes()
}
