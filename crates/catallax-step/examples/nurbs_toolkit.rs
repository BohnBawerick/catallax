//! Builds each phase 2.2 surface, checks interpolation and exports real STEP files.
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_shapes::spline_tools::{
    RationalPoint, RationalSplineCurve as Curve, RationalSplineSurface as Surface,
};
use catallax_step::write_spline_surfaces_step;
use serde_json::json;
use std::{fs, path::Path};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn p(x: i64, y: i64, z: i64) -> RationalPoint {
    [x, y, z].map(Rational::from_i64)
}
fn line(a: RationalPoint, b: RationalPoint, dec: &Decider) -> Curve {
    Curve::new(
        1,
        vec![q("0"), q("0"), q("1"), q("1")],
        vec![a, b],
        vec![q("1"); 2],
        dec,
    )
    .unwrap()
}
fn law(u: &Rational, v: &Rational) -> RationalPoint {
    [
        u * &q("10"),
        v * &q("8"),
        u * &(&q("1") - u) * q("4") + v * &(&q("1") - v) * q("3") + u * v * q("2"),
    ]
}
fn section(t: &Rational, u_curve: bool, dec: &Decider) -> Curve {
    let eval = |s: &Rational| if u_curve { law(s, t) } else { law(t, s) };
    let a = eval(&q("0"));
    let m = eval(&q("1/2"));
    let b = eval(&q("1"));
    let control = std::array::from_fn(|i| &m[i] * &q("2") - (&a[i] + &b[i]) * q("1/2"));
    Curve::new(
        2,
        vec![q("0"), q("0"), q("0"), q("1"), q("1"), q("1")],
        vec![a, control, b],
        vec![q("1"); 3],
        dec,
    )
    .unwrap()
}
fn export(out: &Path, name: &str, surfaces: &[Surface], dec: &mut Decider) -> serde_json::Value {
    let mut bytes = Vec::new();
    let report =
        write_spline_surfaces_step(surfaces, &q("1/1000000"), 17, dec, &mut bytes).unwrap();
    fs::write(out.join(format!("{name}.step")), &bytes).unwrap();
    json!({"name":name,"surfaces":report.surfaces,"patches":report.patches,"bytes":bytes.len(),
        "exact":report.exact,"tolerance_mm":report.tolerance_mm.to_string(),"whole_surface_rounding_bound_mm":report.max_deviation_mm.to_string()})
}
fn run(out: &Path) {
    fs::create_dir_all(out).unwrap();
    let mut dec = Decider::default();
    let nodes = [q("0"), q("1/3"), q("1")];
    let points: Vec<_> = nodes
        .iter()
        .flat_map(|u| nodes.iter().map(move |v| law(u, v)))
        .collect();
    let grid = Surface::through_grid([&nodes, &nodes], [2, 2], &points, &mut dec).unwrap();
    for u in &nodes {
        for v in &nodes {
            assert_eq!(grid.point_at(u, v, &dec).unwrap(), law(u, v));
        }
    }
    let profiles: Vec<_> = nodes.iter().map(|v| section(v, true, &dec)).collect();
    let guides: Vec<_> = nodes.iter().map(|u| section(u, false, &dec)).collect();
    let gordon = Surface::gordon(&profiles, &nodes, &guides, &nodes, &mut dec).unwrap();
    let skin = Surface::guided_skin(
        &profiles,
        &nodes,
        &[guides[1].clone()],
        &[nodes[1].clone()],
        &mut dec,
    )
    .unwrap();
    for u in [q("0"), q("1/5"), q("1/3"), q("4/5"), q("1")] {
        for v in [q("0"), q("2/7"), q("1/3"), q("1")] {
            assert_eq!(gordon.point_at(&u, &v, &dec).unwrap(), law(&u, &v));
            // The single guide preserves every complete profile, not only its guide intersection.
            for (i, node) in nodes.iter().enumerate() {
                assert_eq!(
                    skin.point_at(&u, node, &dec).unwrap(),
                    profiles[i].point_at(&u, &dec).unwrap()
                );
            }
        }
    }
    let bottom = Curve::new(
        2,
        vec![q("0"), q("0"), q("0"), q("1"), q("1"), q("1")],
        vec![p(0, 0, 0), p(4, 0, 4), p(8, 0, 0)],
        vec![q("1"), q("2"), q("1")],
        &dec,
    )
    .unwrap()
    .insert_knot(&q("1/3"), 1, &dec)
    .unwrap();
    let top = Curve::new(
        2,
        vec![q("0"), q("0"), q("0"), q("1"), q("1"), q("1")],
        vec![p(0, 6, 0), p(4, 6, 4), p(8, 6, 0)],
        vec![q("1"), q("2"), q("1")],
        &dec,
    )
    .unwrap()
    .raise_degree(3, &dec)
    .unwrap();
    let coons = Surface::coons(
        &[bottom.clone(), top],
        &[
            line(p(0, 0, 0), p(0, 6, 0), &dec),
            line(p(8, 0, 0), p(8, 6, 0), &dec),
        ],
        &mut dec,
    )
    .unwrap();
    for u in [q("0"), q("1/7"), q("1/3"), q("3/4"), q("1")] {
        for v in [q("0"), q("2/5"), q("1")] {
            let mut expected = bottom.point_at(&u, &dec).unwrap();
            expected[1] = &v * &q("6");
            assert_eq!(coons.point_at(&u, &v, &dec).unwrap(), expected);
        }
    }
    let vertices = [p(0, 0, 0), p(6, 0, 1), p(9, 4, 2), p(4, 8, 1), p(-2, 4, 0)];
    let boundaries: Vec<_> = (0..vertices.len())
        .map(|i| {
            line(
                vertices[i].clone(),
                vertices[(i + 1) % vertices.len()].clone(),
                &dec,
            )
        })
        .collect();
    let fill = Surface::n_sided_fill(&boundaries, &p(3, 3, 2), &mut dec).unwrap();
    for i in 0..fill.len() {
        let previous = (i + fill.len() - 1) % fill.len();
        for t in [q("0"), q("1/3"), q("1")] {
            assert_eq!(
                fill[i].point_at(&t, &q("0"), &dec).unwrap(),
                boundaries[i].point_at(&(&t * &q("1/2")), &dec).unwrap()
            );
            assert_eq!(
                fill[i].point_at(&t, &q("1"), &dec).unwrap(),
                fill[previous].point_at(&q("1"), &t, &dec).unwrap()
            );
        }
    }
    // The exact editing API must retain a weighted curve at knots and inside spans.
    let raised = bottom.raise_degree(4, &dec).unwrap();
    let split = bottom.split(&q("2/5"), &dec).unwrap();
    let joined =
        Curve::joined(&[(split[0].clone(), true), (split[1].clone(), true)], &dec).unwrap();
    for t in [q("0"), q("1/7"), q("1/3"), q("2/5"), q("3/4"), q("1")] {
        assert_eq!(
            raised.point_at(&t, &dec).unwrap(),
            bottom.point_at(&t, &dec).unwrap()
        );
        let mapped = if t <= q("2/5") {
            &t / &q("2/5")
        } else {
            q("1") + (&t - &q("2/5")) / q("3/5")
        };
        assert_eq!(
            joined.point_at(&mapped, &dec).unwrap(),
            bottom.point_at(&t, &dec).unwrap()
        );
    }
    let edited = coons
        .insert_knot(1, &q("2/5"), 1, &dec)
        .unwrap()
        .raise_degree(1, 3, &dec)
        .unwrap();
    let pieces = edited.split(0, &q("1/2"), &dec).unwrap();
    assert_eq!(
        pieces[0].point_at(&q("1/2"), &q("2/5"), &dec).unwrap(),
        coons.point_at(&q("1/2"), &q("2/5"), &dec).unwrap()
    );
    let mut reports = vec![
        export(out, "grid", &[grid], &mut dec),
        export(out, "coons", &[coons], &mut dec),
        export(out, "gordon", &[gordon], &mut dec),
        export(out, "guided-skin", &[skin], &mut dec),
        export(out, "n-sided-fill", &fill, &mut dec),
    ];
    // Decimal rounding is checked by an independent rational control bound.
    reports.push(export(out, "edited-coons", &pieces, &mut dec));
    let report = json!({"format":"catallax-nurbs-toolkit/1","checks":"exact interpolation, complete rational boundaries, guide profiles, fill seams, curve edits, surface edits","exports":reports});
    fs::write(
        out.join("report.json"),
        serde_json::to_vec_pretty(&report).unwrap(),
    )
    .unwrap();
    println!("{report}");
}
fn main() {
    let out = std::env::args()
        .nth(1)
        .expect("usage: nurbs_toolkit OUTPUT_DIRECTORY");
    let (_, work) = with_operation_budget(100_000_000, || run(Path::new(&out)));
    assert!(!work.exhausted);
    println!("work={}", work.used);
}
