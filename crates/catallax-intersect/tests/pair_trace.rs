use catallax_intersect::TracePiece;
use catallax_numbers::{
    Ball, Decider, Outcome, Rational, Real, Refusal, Sign, with_operation_budget,
};
use catallax_shapes::{BSplineSurface, Point3, RationalPatch};

fn q(n: i64, d: i64) -> Rational {
    Rational::parse(&format!("{n}/{d}")).unwrap()
}

// Both inputs curve. Their meeting is 13 x² + 7 y² = 10,
// z = x² + y², over the positive square [1/5, 4/5]².
fn patches(dec: &mut Decider) -> [RationalPatch; 2] {
    let linear = [q(1, 5), q(1, 2), q(4, 5)];
    let square = [q(1, 25), q(4, 25), q(16, 25)];
    std::array::from_fn(|side| {
        let mut points = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                let z = if side == 0 {
                    &square[i] + &square[j]
                } else {
                    &Rational::one() + &(&q(3, 10) * &(&square[j] - &square[i]))
                };
                points.push(Point3::new(
                    linear[i].clone().into(),
                    linear[j].clone().into(),
                    z.into(),
                ));
            }
        }
        let knots = [0, 0, 0, 1, 1, 1].map(Rational::from_i64).to_vec();
        let surface = BSplineSurface::new(
            [2, 2],
            [knots.clone(), knots],
            [3, 3],
            points,
            vec![Real::one(); 9],
            dec,
        )
        .unwrap();
        RationalPatch::new(surface, dec).unwrap()
    })
}

fn zero(value: &Real, dec: &mut Decider) -> Result<(), Refusal> {
    match dec.sign(value) {
        Outcome::Proven(Sign::Zero) => Ok(()),
        Outcome::Refused(refusal) => Err(refusal),
        Outcome::NotExact(answer) => Err(Refusal::from(answer)),
        other => panic!("the retained graph point must satisfy its defining equation: {other:?}"),
    }
}

#[test]
fn curved_pair_exact_incidence_and_whole_panel_second_jets() {
    let mut dec = Decider::default();
    let [a, b] = patches(&mut dec);
    let (result, work) = with_operation_budget(20_000_000, || -> Result<(), Refusal> {
        let mut dec = Decider::default();
        let meeting = a.patch_meeting(&b, &mut dec)?;
        assert_eq!(meeting.branches().len(), 1);
        assert!(!meeting.branches()[0].is_closed());
        let pieces = meeting.trace_branch(0, &mut dec)?;
        for (index, piece) in pieces.into_iter().enumerate() {
            let piece = TracePiece::from_patch_pair(piece);
            let parameter =
                Ball::from_rational(&q(49, 100), 128).hull(&Ball::from_rational(&q(51, 100), 128));
            let jets = piece.taylor(&parameter, 2, &dec)?;
            let [[x, y, z], [dx, dy, dz], [xx, yy, zz]]: [[Ball; 3]; 3] = jets.try_into().unwrap();
            let scale = |n: i64, value: &Ball| &Ball::from_rational(&q(n, 1), 128) * value;
            let ellipse = &(&scale(13, &(&x * &x)) + &scale(7, &(&y * &y)))
                - &Ball::from_rational(&q(10, 1), 128);
            let velocity = &scale(26, &(&x * &dx)) + &scale(14, &(&y * &dy));
            let curvature = &scale(13, &(&scale(2, &(&x * &xx)) + &(&dx * &dx)))
                + &scale(7, &(&scale(2, &(&y * &yy)) + &(&dy * &dy)));
            let height = &z - &(&(&x * &x) + &(&y * &y));
            let height_velocity = &dz - &scale(2, &(&(&x * &dx) + &(&y * &dy)));
            let height_curvature = &zz
                - &(&(&scale(2, &(&x * &xx)) + &(&dx * &dx))
                    + &(&scale(2, &(&y * &yy)) + &(&dy * &dy)));
            for residual in [
                ellipse,
                velocity,
                curvature,
                height,
                height_velocity,
                height_curvature,
            ] {
                assert!(
                    residual.contains_zero(),
                    "independent graph identity on panel {index}: {residual:?}"
                );
            }
            assert!(
                &x.hi() - &x.lo() < q(1, 10),
                "a useful panel enclosure, not the entire chart"
            );
            assert!(
                &y.hi() - &y.lo() < q(1, 10),
                "a useful panel enclosure, not the entire chart"
            );
            println!("piece {index}: Taylor work {}", dec.stats().work);
            let point = piece.evaluate_exact(&q(1, 3), &mut dec)?;
            println!("piece {index}: point work {}", dec.stats().work);
            let ellipse = &(&Real::from_i64(13) * &point.x.square())
                + &(&Real::from_i64(7) * &point.y.square());
            zero(&(&ellipse - &Real::from_i64(10)), &mut dec)?;
            zero(
                &(&point.z - &(&point.x.square() + &point.y.square())),
                &mut dec,
            )?;
        }
        Ok(())
    });
    println!("four-parameter graph: {result:?}; {} work", work.used);
    assert!(
        !work.exhausted,
        "the graph proof stays inside the unchanged allowance"
    );
    result.unwrap();
}
