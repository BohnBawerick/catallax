#![allow(missing_docs)]
mod common;
use catallax_shapes::*;
use common::*;
use std::hash::{Hash, Hasher};

fn hash(x: &impl Hash) -> u64 {
    let mut h = std::collections::hash_map::DefaultHasher::new();
    x.hash(&mut h);
    h.finish()
}

blind_test!(cases_01_to_09_quaternions, {
    for w in ["1", "2", "-3", "1/7"] {
        assert_eq!(rotation([w, "0", "0", "0"]), Rotation::identity());
    }
    let cases = [
        (
            ["1", "0", "0", "1"],
            [["0", "-1", "0"], ["1", "0", "0"], ["0", "0", "1"]],
        ),
        (
            ["1", "1", "1", "1"],
            [["0", "0", "1"], ["1", "0", "0"], ["0", "1", "0"]],
        ),
        (
            ["2", "1", "0", "0"],
            [["1", "0", "0"], ["0", "3/5", "-4/5"], ["0", "4/5", "3/5"]],
        ),
        (
            ["1", "2", "3", "4"],
            [
                ["-2/3", "2/15", "11/15"],
                ["2/3", "-1/3", "2/3"],
                ["1/3", "14/15", "2/15"],
            ],
        ),
        (
            ["0", "0", "0", "1"],
            [["-1", "0", "0"], ["0", "-1", "0"], ["0", "0", "1"]],
        ),
    ];
    for (q, m) in cases {
        let rot = rotation(q);
        assert_eq!(rot.matrix(), &m.map(|row| row.map(r)));
    }
    let a = rotation(["1", "2", "3", "4"]);
    for c in ["-1", "5", "1/7"] {
        let [w, x, y, z] = ["1", "2", "3", "4"].map(|x| r(x) * r(c));
        let b = Rotation::from_quaternion(w, x, y, z).unwrap();
        assert_eq!(a, b);
        assert_eq!(hash(&a), hash(&b));
    }
    assert_eq!(
        a.then(&rotation(["2", "-1", "0", "5"])),
        rotation(["-16", "-12", "20", "10"])
    );
    rejected(
        Rotation::from_quaternion(r("0"), r("0"), r("0"), r("0")),
        "shapes.zero-quaternion",
    );
    assert_eq!(
        rotation(["1", "0", "0", "1"]).apply(&v("1", "0", "0")),
        v("0", "1", "0")
    );
});

blind_test!(cases_07_08_generated_exact_rotation_laws, {
    let mut rng = Rng(0x20708);
    for _ in 0..256 {
        // Full 64-bit numerators and denominators, not f64-generated inputs.
        let a: [Q; 4] = std::array::from_fn(|_| Q::from((rng.next() as i64, rng.next() | 1)));
        let b: [Q; 4] = std::array::from_fn(|_| Q::from((rng.next() as i64, rng.next() | 1)));
        let build = |q: &[Q; 4]| {
            let [w, x, y, z] = q.clone().map(|x| eq(&x));
            Rotation::from_quaternion(w, x, y, z).unwrap()
        };
        let rot = build(&a);
        let matrix = quaternion_matrix(&a);
        assert_eq!(rot.matrix(), &matrix.clone().map(|row| row.map(|x| eq(&x))));
        for i in 0..3 {
            for j in 0..3 {
                let product: Q = (0..3)
                    .map(|k| oq(&rot.matrix()[i][k]) * oq(&rot.matrix()[j][k]))
                    .sum();
                assert_eq!(product, q(i64::from(i == j), 1));
                assert_eq!(rot.inverse().matrix()[i][j], rot.matrix()[j][i]);
            }
        }
        let m = rot.matrix().clone().map(|row| row.map(|x| oq(&x)));
        assert_eq!(dot(&m[0], &cross(&m[1], &m[2])), q(1, 1));
        assert_eq!(
            rot.inverse(),
            build(&[a[0].clone(), -a[1].clone(), -a[2].clone(), -a[3].clone()])
        );
        assert_eq!(rot.then(&rot.inverse()), Rotation::identity());
        assert_eq!(rot.then(&build(&b)), build(&hamilton(&b, &a)));
    }
});

blind_test!(case_10_pose_order_and_vectors, {
    let pose = Pose::new(rotation(["1", "0", "0", "1"]), v("1", "2", "3"));
    assert_eq!(pose.apply_point(&p("1", "0", "0")), p("1", "3", "3"));
    assert_eq!(
        pose.inverse().apply_point(&p("1", "3", "3")),
        p("1", "0", "0")
    );
    point_eq(
        &pose.apply_real_point(&p("1", "0", "0").to_real()),
        &p("1", "3", "3").to_real(),
    );
    let mut rng = Rng(0x210);
    for _ in 0..200 {
        let a = rng.pose();
        let b = rng.pose();
        let point = rp(std::array::from_fn(|_| rng.rational()));
        assert_eq!(
            a.then(&b).apply_point(&point),
            b.apply_point(&a.apply_point(&point))
        );
        assert_eq!(a.then(&a.inverse()), Pose::identity());
        assert_eq!(a.inverse().apply_point(&a.apply_point(&point)), point);
        let vector = RVector3::new(point.x.clone(), point.y.clone(), point.z.clone());
        assert_eq!(a.apply_vector(&vector), a.rotation.apply(&vector));
        vector_eq(
            &a.rotation.apply_real(&vector.to_real()),
            &a.apply_vector(&vector),
        );
    }
});

blind_test!(rational_and_real_vector_identities, {
    let a = v("1/3", "-2", "5");
    let b = v("3", "4", "-1");
    assert_eq!(a.dot(&b), r("-12"));
    assert_eq!(a.cross(&b), v("-18", "46/3", "22/3"));
    assert_eq!(a.norm_squared(), r("262/9"));
    assert_eq!(a.plus(&b), v("10/3", "2", "4"));
    assert_eq!(a.minus(&b), v("-8/3", "-6", "6"));
    assert_eq!(a.scaled(&r("-3")), v("-1", "6", "-15"));
    assert!(RVector3::zero().is_zero());
    assert!(!a.is_zero());
    let origin = RPoint3::origin();
    assert_eq!(origin, p("0", "0", "0"));
    let point = p("1/3", "-2", "5");
    assert_eq!(point.minus(&origin), a);
    assert_eq!(origin.plus(&a), point);
    let ar = a.to_real();
    let br: Vector3 = b.clone().into();
    vector_eq(&ar.cross(&br), &v("-18", "46/3", "22/3"));
    vector_eq(&ar.plus(&br), &v("10/3", "2", "4"));
    vector_eq(&ar.minus(&br), &v("-8/3", "-6", "6"));
    vector_eq(&ar.scaled(&real("-3")), &v("-1", "6", "-15"));
    real_eq(&ar.dot(&br), &real("-12"));
    real_eq(&ar.norm_squared(), &real("262/9"));
    let pr: Point3 = point.into();
    vector_eq(&pr.minus(&origin.to_real()), &a);
    point_eq(&origin.to_real().plus(&ar), &pr);
});
