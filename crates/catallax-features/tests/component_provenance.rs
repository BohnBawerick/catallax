//! Geometry-only shell subdivision does not relax named extraction contracts.

use std::collections::BTreeSet;

use catallax_features::draft::{
    DraftDefinition,
    arrangement::{Patch, Skin},
};
use catallax_numbers::{Decider, Rational, Real, with_operation_budget};
use catallax_shapes::{Plane, RPoint3, Vector3};

fn cube() -> Skin {
    let vertices = [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
    ]
    .map(|xyz| {
        let [x, y, z] = xyz.map(Rational::from_i64);
        RPoint3::new(x, y, z).to_real()
    })
    .to_vec();
    let patches = [
        ([-1, 0, 0], [0, 4, 7, 3]),
        ([1, 0, 0], [1, 2, 6, 5]),
        ([0, -1, 0], [0, 1, 5, 4]),
        ([0, 1, 0], [3, 7, 6, 2]),
        ([0, 0, -1], [0, 3, 2, 1]),
        ([0, 0, 1], [4, 5, 6, 7]),
    ]
    .into_iter()
    .enumerate()
    .map(|(i, (normal, cycle))| {
        let [x, y, z] = normal.map(Real::from_i64);
        Patch {
            plane: Plane {
                origin: vertices[cycle[0]].clone(),
                normal: Vector3::new(x, y, z),
                ref_dir: if normal[2] == 0 {
                    Vector3::new(Real::zero(), Real::zero(), Real::one())
                } else {
                    Vector3::new(Real::one(), Real::zero(), Real::zero())
                },
            },
            loops: vec![cycle.to_vec()],
            origins: vec![
                serde_json::from_value(serde_json::json!({
                    "step":"stock","role":{"side":{"segment":format!("face-{i}")}}
                }))
                .unwrap(),
            ],
        }
    })
    .collect();
    Skin { vertices, patches }
}

#[test]
fn named_draft_and_partition_reject_unnamed_boundaries() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut source = cube();
        source.patches[0].origins.clear();
        let mut dec = Decider::default();
        assert_eq!(
            source.components(&mut dec).unwrap_err().code,
            "draft.nonclosing-boundary"
        );
        let definition = DraftDefinition::new(
            Vector3::new(Real::zero(), Real::zero(), Real::one()),
            source.patches[4].plane.clone(),
            &Rational::zero(),
            &mut dec,
        )
        .unwrap();
        assert_eq!(
            source
                .draft(&BTreeSet::new(), &definition, &mut dec)
                .unwrap_err()
                .code,
            "draft.nonclosing-boundary"
        );
        let partition = source.partition(
            |point, dec| {
                for coordinate in [&point.x, &point.y, &point.z] {
                    let value = dec
                        .rational_value(coordinate)?
                        .expect("rational cube cell witness");
                    if value <= Rational::zero() || value >= Rational::one() {
                        return Ok(false);
                    }
                }
                Ok(true)
            },
            |_, _, _| Ok(Vec::new()),
            &mut dec,
        );
        assert_eq!(partition.unwrap_err().code, "draft.nonclosing-boundary");
    });
    assert!(!work.exhausted, "{work:?}");
}
