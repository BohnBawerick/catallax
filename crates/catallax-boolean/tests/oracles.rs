//! Implementation-independent checks of the exact test oracles and fixtures.
mod common;
use common::*;

#[test]
fn exact_primitive_oracle_resolves_tiny_separations() {
    let ball = ball_at(5, [0; 3]);
    assert_eq!(ball.locate(&p("3", "4", "0")), Location::On);
    assert_eq!(ball.locate(&p("3", "4", "1e-30")), Location::Outside);
    assert_eq!(
        ball.locate(&RPoint3::new(r(3), r(4) - q("1e-30"), r(0))),
        Location::Inside
    );
    let cylinder = cyl_at(2, 5, [0; 3], false);
    assert_eq!(cylinder.locate(&p("6/5", "8/5", "5")), Location::On);
    assert_eq!(cylinder.locate(&p("2", "2", "1")), Location::Outside);
    let c = case("09-overlap");
    let point = RPoint3::new(r(1) - q("1e-30") / r(2), q("1/2"), q("1/2"));
    assert_eq!(
        set_rule(BoolOp::Intersect, c.a.locate(&point), c.b.locate(&point)),
        Some(Location::Inside)
    );
    let c = case("13-bite");
    let point = RPoint3::new(r(5), r(5), r(10) - q("1e-31"));
    assert_eq!(
        set_rule(BoolOp::Intersect, c.a.locate(&point), c.b.locate(&point)),
        Some(Location::Inside)
    );
}

#[test]
fn certified_normal_step_stays_inside_the_thin_overlap_cell() {
    let c = case("09-overlap");
    // On the right face of the thin intersection, the next plane is only
    // 10^-30 away. The oracle must derive a step smaller than that gap.
    let point = p("1", "1/2", "1/2");
    let normal = v(1, 0, 0);
    let step = safe_step(&c.a, &c.b, &point, &normal);
    assert!(step > r(0) && step < q("1e-30"));
    let inward = point.plus(&normal.scaled(&(-&step)));
    let outward = point.plus(&normal.scaled(&step));
    assert_eq!(
        set_rule(BoolOp::Intersect, c.a.locate(&inward), c.b.locate(&inward)),
        Some(Location::Inside)
    );
    assert_eq!(
        set_rule(
            BoolOp::Intersect,
            c.a.locate(&outward),
            c.b.locate(&outward)
        ),
        Some(Location::Outside)
    );
}

#[test]
fn architect_contact_fallbacks_stay_on_their_primitive_faces() {
    for (id, point) in [
        ("12", p("5", "5", "10")),
        ("19", p("5", "5", "10")),
        ("20", p("5", "5", "10")),
        ("21", p("5", "5", "10")),
        ("26", p("5", "0", "5")),
        ("29", p("0", "2", "0")),
        ("34", p("1", "0", "0")),
        ("37", p("5", "0", "0")),
        ("38", p("5", "0", "0")),
    ] {
        let c = case(id);
        let mut m = Model::new();
        let a = c.a.build(&mut m);
        let b = c.b.build(&mut m);
        for (body, source, other) in [(a, &c.a, &c.b), (b, &c.b, &c.a)] {
            assert_eq!(source.locate(&point), Location::On);
            assert_eq!(other.locate(&point), Location::On);
            let faces: Vec<_> = m
                .faces_of_body(body)
                .into_iter()
                .filter(|&f| {
                    implicit_normal(m.geom().surface(m.face(f).surface), &point)
                        .0
                        .is_zero()
                })
                .collect();
            assert_eq!(faces.len(), 1);
            let f = faces[0];
            let q = contact_witness(&m, f, source, other, &point);
            assert_eq!(source.locate(&q), Location::On);
            assert_ne!(other.locate(&q), Location::On);
            assert_eq!(
                source
                    .constraints(&q)
                    .iter()
                    .filter(|x| x.is_zero())
                    .count(),
                1
            );
            assert_eq!(
                implicit_normal(m.geom().surface(m.face(f).surface), &q).0,
                r(0)
            );
        }
    }
}

#[test]
fn all_numbered_operands_pass_the_finished_topology_checker() {
    // This runs before pieces 4/5 exist. It validates the actual operands fed
    // to ignored tests, including their quaternion poses and 10^-30 offsets.
    for &id in IDS {
        let c = case(id);
        let mut m = Model::new();
        let a = c.a.build(&mut m);
        let b = c.b.build(&mut m);
        assert!(m.check(a).is_empty(), "{id} A");
        assert!(m.check(b).is_empty(), "{id} B");
    }
}
