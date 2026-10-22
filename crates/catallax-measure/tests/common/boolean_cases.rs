use super::*;
use catallax_boolean::{BoolOp, boolean};
use catallax_topology::{BodyId, Model, block, cylinder, sphere};

#[derive(Clone, Copy, Debug)]
pub enum BooleanCase {
    BlockUnion,
    BlockIntersection,
    BlockDifference,
    ThroughHole,
    ThroughIntersection,
    FourLumps,
    Steinmetz,
    CylinderBall,
    CrossingHoles,
}
pub const BOOLEAN_CASES: [BooleanCase; 9] = [
    BooleanCase::BlockUnion,
    BooleanCase::BlockIntersection,
    BooleanCase::BlockDifference,
    BooleanCase::ThroughHole,
    BooleanCase::ThroughIntersection,
    BooleanCase::FourLumps,
    BooleanCase::Steinmetz,
    BooleanCase::CylinderBall,
    BooleanCase::CrossingHoles,
];
pub struct BooleanFixture {
    pub model: Model,
    pub body: BodyId,
    pub volume: Expr,
    // Case 41 specifies volume only. No guessed area or centroid assertion.
    pub mass: Option<Expected>,
}
fn translation(x: i64, y: i64, z: i64) -> Pose {
    Pose::new(
        Rotation::identity(),
        (v(x, y, z)).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap()
}
fn xpose(x: i64, y: i64, z: i64) -> Pose {
    Pose::new(
        Rotation::from_quaternion(r(1), r(0), r(1), r(0)).unwrap(),
        (v(x, y, z)).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap()
}
pub fn build_boolean(case: BooleanCase) -> BooleanFixture {
    use BooleanCase::*;
    let mut m = Model::new();
    let (a, b, op, volume, area, centroid) = match case {
        BlockUnion | BlockIntersection | BlockDifference => {
            let a = block(&mut m, [r(2), r(2), r(2)], &Pose::identity()).unwrap();
            let b = block(&mut m, [r(2), r(2), r(2)], &translation(1, 1, 1)).unwrap();
            let (op, vol, ar, c) = match case {
                BlockUnion => (BoolOp::Unite, rat(15, 1), rat(42, 1), rat(3, 2)),
                BlockIntersection => (BoolOp::Intersect, rat(1, 1), rat(6, 1), rat(3, 2)),
                _ => (BoolOp::Subtract, rat(7, 1), rat(24, 1), rat(13, 14)),
            };
            (a, b, op, vol, Some(ar), [c.clone(), c.clone(), c])
        }
        ThroughHole | ThroughIntersection | FourLumps | CrossingHoles => {
            let a = block(&mut m, [r(10), r(10), r(10)], &Pose::identity()).unwrap();
            let radius = match case {
                FourLumps => 5,
                CrossingHoles => 2,
                _ => 3,
            };
            let b = cylinder(&mut m, r(radius), r(12), &translation(5, 5, -1)).unwrap();
            let (op, vol, ar) = match case {
                ThroughHole => (
                    BoolOp::Subtract,
                    rat(1000, 1) - rat(90, 1) * pi(),
                    Some(rat(600, 1) + rat(42, 1) * pi()),
                ),
                ThroughIntersection => (
                    BoolOp::Intersect,
                    rat(90, 1) * pi(),
                    Some(rat(78, 1) * pi()),
                ),
                FourLumps => (
                    BoolOp::Subtract,
                    rat(1000, 1) - rat(250, 1) * pi(),
                    Some(rat(600, 1) + rat(50, 1) * pi()),
                ),
                _ => (
                    BoolOp::Subtract,
                    rat(1000, 1) - rat(80, 1) * pi() + rat(128, 3),
                    None,
                ),
            };
            (a, b, op, vol, ar, [rat(5, 1), rat(5, 1), rat(5, 1)])
        }
        Steinmetz => {
            let a = cylinder(&mut m, r(1), r(4), &translation(0, 0, -2)).unwrap();
            let b = cylinder(&mut m, r(1), r(4), &xpose(-2, 0, 0)).unwrap();
            (
                a,
                b,
                BoolOp::Intersect,
                rat(16, 3),
                Some(rat(16, 1)),
                [rat(0, 1), rat(0, 1), rat(0, 1)],
            )
        }
        CylinderBall => {
            let a = cylinder(&mut m, r(3), r(20), &translation(0, 0, -10)).unwrap();
            let b = sphere(&mut m, r(5), &Pose::identity()).unwrap();
            // Central cylinder: 72*pi volume, 48*pi side area. Two height-1
            // sphere caps add 28*pi/3 volume and 20*pi surface area.
            (
                a,
                b,
                BoolOp::Intersect,
                rat(244, 3) * pi(),
                Some(rat(68, 1) * pi()),
                [rat(0, 1), rat(0, 1), rat(0, 1)],
            )
        }
    };
    let result = boolean(&mut m, a, b, op, &mut Decider::default()).unwrap();
    assert!(
        result.marks.is_empty(),
        "exact cross-check construction unexpectedly marked"
    );
    let mut body = result.body;
    if matches!(case, CrossingHoles) {
        let c = cylinder(&mut m, r(2), r(12), &xpose(-1, 5, 5)).unwrap();
        let result = boolean(&mut m, body, c, BoolOp::Subtract, &mut Decider::default()).unwrap();
        assert!(result.marks.is_empty());
        body = result.body;
    }
    assert!(m.check(body).is_empty());
    // Meeting vertices do not impose rational exactness, even on all-plane results.
    let mass = area.map(|area| Expected {
        volume: volume.clone(),
        area,
        centroid: Some(centroid),
        exact: false,
    });
    BooleanFixture {
        model: m,
        body,
        volume,
        mass,
    }
}
