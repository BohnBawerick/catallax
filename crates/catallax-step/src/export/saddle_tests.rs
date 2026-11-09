use super::*;
use catallax_numbers::{
    prefer_reduced_angle_series, reuse_operation_enclosures, with_operation_budget,
};
use catallax_shapes::{BallSide, Cylinder, RollingSaddleDefinition, Vector3};

fn point(x: i64, y: i64, z: i64) -> Point3 {
    Point3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
}
fn vector(x: i64, y: i64, z: i64) -> Vector3 {
    Vector3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
}

#[test]
fn crossing_table_keys_and_reflections_match_independent_cold_tables() {
    let tolerance = Rational::from_ratio(1, TABLE_RADIANS).unwrap();
    let mut canonical: Option<(TableKey, ParameterTable)> = None;
    for (y, z) in [(4, 4), (-4, 4), (4, -4), (-4, -4)] {
        for row in 0..2 {
            let ((key, table), work) = with_operation_budget(1_000_000, || {
                reuse_operation_enclosures();
                prefer_reduced_angle_series();
                let mut dec = Decider::default();
                let wide = Cylinder {
                    origin: Point3::origin(),
                    axis: vector(0, 0, 1),
                    radius_sq: Real::from_i64(25),
                    ref_dir: vector(1, 0, 0),
                };
                let narrow = Cylinder {
                    axis: vector(0, 1, 0),
                    ..wide.clone()
                };
                let definition = RollingSaddleDefinition::crossing(
                    &wide,
                    &narrow,
                    &Real::from(Rational::from_ratio(1, 2).unwrap()),
                    [BallSide::Outside; 2],
                    &point(3, y, z),
                    &mut dec,
                )
                .unwrap();
                let band = RingBand::Saddle(RollingSaddle::fit(definition, &mut dec).unwrap());
                let (key, negative) = band.crossing_key(row, &tolerance).unwrap();
                let cold = band.table(row, &tolerance, &dec).unwrap();
                let normalized = reflect_table(&cold, negative, &dec).unwrap();
                let restored = reflect_table(&normalized, negative, &dec).unwrap();
                assert_eq!(cold.knots, restored.knots);
                assert_eq!(cold.angle, restored.angle);
                assert_eq!(cold.height, restored.height);
                assert_eq!(cold.bound, restored.bound);
                assert_eq!(cold.height_bound, restored.height_bound);
                assert!(band.crossing_key(row, &tolerance.mul_pow2(-1)).unwrap().0 != key);
                (key, normalized)
            });
            assert!(!work.exhausted);
            if let Some((expected_key, expected)) = &canonical {
                assert_eq!(&key, expected_key);
                assert_eq!(table.knots, expected.knots);
                assert_eq!(table.angle, expected.angle);
                assert_eq!(table.height, expected.height);
                assert_eq!(table.bound, expected.bound);
                assert_eq!(table.height_bound, expected.height_bound);
            } else {
                canonical = Some((key, table));
            }
        }
    }
    let table = canonical.unwrap().1;
    let mut expired = Decider::default();
    expired.set_deadline(Some(
        std::time::Instant::now() - std::time::Duration::from_secs(1),
    ));
    assert_eq!(
        reflect_table(&table, true, &expired).unwrap_err().code,
        "numbers.deadline"
    );
    let (result, work) =
        with_operation_budget(0, || reflect_table(&table, true, &Decider::default()));
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    assert!(work.exhausted);
}
