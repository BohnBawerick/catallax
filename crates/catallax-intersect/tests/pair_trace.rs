//! Exact incidence and whole-panel jets of a retained four-parameter graph.

use catallax_intersect::TracePiece;
use catallax_numbers::{Ball, Decider, Outcome, Rational, Real, Refusal, with_operation_budget};
use catallax_shapes::{BSplineSurface, Point3, RationalPatch};

fn q(n: i64, d: i64) -> Rational {
    Rational::parse(&format!("{n}/{d}")).unwrap()
}

// Both inputs curve. Their meeting is 13 x² + 7 y² = 10,
// z = x² + y², over the positive square [1/5, 4/5]².
fn patches(dec: &mut Decider, swapped: bool) -> [RationalPatch; 2] {
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
                let [x, y] = if side == 0 && swapped { [j, i] } else { [i, j] };
                points.push(Point3::new(
                    linear[x].clone().into(),
                    linear[y].clone().into(),
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
    let separated = dec.stats().coincident_separation;
    match dec.sign(value) {
        Outcome::ProvenCoincident(_) => {
            assert_eq!(
                separated,
                dec.stats().coincident_separation,
                "exact defining-equation authority, not coincident separation"
            );
            Ok(())
        }
        Outcome::Refused(refusal) => Err(refusal),
        Outcome::NotExact(answer) => Err(Refusal::from(answer)),
        other => panic!("the retained graph point must satisfy its defining equation: {other:?}"),
    }
}

fn check_panel(piece: &TracePiece, index: usize, dec: &Decider) -> Result<Vec<String>, Refusal> {
    let parameter =
        Ball::from_rational(&q(49, 100), 128).hull(&Ball::from_rational(&q(51, 100), 128));
    let jets = piece.taylor(&parameter, 2, dec)?;
    let [[x, y, z], [dx, dy, dz], [xx, yy, zz]]: [[Ball; 3]; 3] = jets.try_into().unwrap();
    let scale = |n: i64, value: &Ball| &Ball::from_rational(&q(n, 1), 128) * value;
    let ellipse =
        &(&scale(13, &(&x * &x)) + &scale(7, &(&y * &y))) - &Ball::from_rational(&q(10, 1), 128);
    let velocity = &scale(26, &(&x * &dx)) + &scale(14, &(&y * &dy));
    let curvature = &scale(13, &(&scale(2, &(&x * &xx)) + &(&dx * &dx)))
        + &scale(7, &(&scale(2, &(&y * &yy)) + &(&dy * &dy)));
    let height = &z - &(&(&x * &x) + &(&y * &y));
    let height_velocity = &dz - &scale(2, &(&(&x * &dx) + &(&y * &dy)));
    let height_curvature = &zz
        - &(&(&scale(2, &(&x * &xx)) + &(&dx * &dx)) + &(&scale(2, &(&y * &yy)) + &(&dy * &dy)));
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
    Ok([[x, y, z], [dx, dy, dz], [xx, yy, zz]]
        .into_iter()
        .flat_map(|xyz| xyz.map(|value| format!("{}..{}", value.lo(), value.hi())))
        .collect())
}

#[test]
fn row_pivot_graph_jets_and_work_are_thread_independent() {
    catallax_numbers::set_run_ahead_wait(true);
    for swapped in [false, true] {
        let mut seen = None;
        for threads in [Some(1), Some(2), Some(4), None] {
            catallax_numbers::set_threads(threads.and_then(std::num::NonZeroUsize::new));
            let mut dec = Decider::default();
            let [a, b] = patches(&mut dec, swapped);
            let (result, work) =
                with_operation_budget(20_000_000, || -> Result<Vec<String>, Refusal> {
                    let mut dec = Decider::default();
                    let meeting = a.patch_meeting(&b, &mut dec)?;
                    let [branch] = meeting.branches() else {
                        panic!("one positive-quadrant arc");
                    };
                    assert!(!branch.is_closed());
                    let mut pieces = meeting.trace_branch(0, &mut dec)?.into_iter();
                    let first = pieces.next().expect("the proven arc has a graph");
                    let mut records = Vec::new();
                    for (index, piece) in std::iter::once(first).chain(pieces).enumerate() {
                        records.extend(check_panel(
                            &TracePiece::from_patch_pair(piece),
                            index,
                            &dec,
                        )?);
                    }
                    Ok(records)
                });
            assert!(
                !work.exhausted,
                "{swapped:?}/{threads:?}: {result:?}, {} work",
                work.used
            );
            let now = (result.unwrap(), work.used);
            println!(
                "graph jets swapped={swapped}/{threads:?}: {} work",
                work.used
            );
            if let Some(first) = &seen {
                assert_eq!(first, &now);
            } else {
                seen = Some(now);
            }
        }
    }
    catallax_numbers::set_threads(None);
    catallax_numbers::set_run_ahead_wait(false);
}

#[test]
fn curved_pair_exact_incidence_and_whole_panel_second_jets() {
    catallax_numbers::set_run_ahead_wait(true);
    let mut seen = None;
    for threads in [Some(1), Some(2), Some(4), None] {
        catallax_numbers::set_threads(threads.and_then(std::num::NonZeroUsize::new));
        let mut dec = Decider::default();
        let [a, b] = patches(&mut dec, false);
        let (result, work) =
            with_operation_budget(20_000_000, || -> Result<Vec<[u64; 6]>, Refusal> {
                let mut dec = Decider::default();
                let meeting = a.patch_meeting(&b, &mut dec)?;
                assert_eq!(meeting.branches().len(), 1);
                assert!(!meeting.branches()[0].is_closed());
                let pieces = meeting.trace_branch(0, &mut dec)?;
                let mut point_bounds = Vec::with_capacity(pieces.len());
                for (index, piece) in pieces.into_iter().enumerate() {
                    let piece = TracePiece::from_patch_pair(piece);
                    let _ = check_panel(&piece, index, &dec)?;
                    let point = piece.evaluate_exact(&q(1, 3), &mut dec)?;
                    let ellipse = &(&Real::from_i64(13) * &point.x.square())
                        + &(&Real::from_i64(7) * &point.y.square());
                    zero(&(&ellipse - &Real::from_i64(10)), &mut dec)?;
                    zero(
                        &(&point.z - &(&point.x.square() + &point.y.square())),
                        &mut dec,
                    )?;
                    point_bounds.push([
                        point.x.approx().lo().to_bits(),
                        point.x.approx().hi().to_bits(),
                        point.y.approx().lo().to_bits(),
                        point.y.approx().hi().to_bits(),
                        point.z.approx().lo().to_bits(),
                        point.z.approx().hi().to_bits(),
                    ]);
                }
                Ok(point_bounds)
            });
        println!("exact graph at {threads:?}: {result:?}; {} work", work.used);
        assert!(
            !work.exhausted,
            "the graph proof stays inside the unchanged allowance"
        );
        let now = (result.unwrap(), work.used);
        if let Some(first) = &seen {
            assert_eq!(
                first, &now,
                "exact graph point bounds and work at {threads:?}"
            );
        } else {
            seen = Some(now);
        }
    }
    catallax_numbers::set_threads(None);
    catallax_numbers::set_run_ahead_wait(false);
}
