#![allow(missing_docs)]
mod common;
use catallax_numbers::{Decider, Real};
use catallax_shapes::*;
use common::*;

blind_test!(cases_01_to_09_quaternions, {
    for w in ["1", "2", "-3", "1/7"] {
        assert!(
            rotation([w, "0", "0", "0"])
                .same_rotation(&Rotation::identity(), &mut Decider::default())
                .unwrap()
        );
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
        let expected = m.map(|row| row.map(r));
        assert_eq!(
            rot.rational_matrix().unwrap(),
            expected.each_ref().map(|row| row.each_ref())
        );
    }
    let a = rotation(["1", "2", "3", "4"]);
    for c in ["-1", "5", "1/7"] {
        let [w, x, y, z] = ["1", "2", "3", "4"].map(|x| r(x) * r(c));
        let b = Rotation::from_quaternion(w, x, y, z).unwrap();
        assert!(a.same_rotation(&b, &mut Decider::default()).unwrap());
    }
    assert!(
        a.then(&rotation(["2", "-1", "0", "5"]))
            .same_rotation(
                &rotation(["-16", "-12", "20", "10"]),
                &mut Decider::default()
            )
            .unwrap()
    );
    rejected(
        Rotation::from_quaternion(r("0"), r("0"), r("0"), r("0")),
        "shapes.zero-quaternion",
    );
    assert_eq!(
        rotation(["1", "0", "0", "1"])
            .apply(&v("1", "0", "0").to_real())
            .as_rational()
            .unwrap(),
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
        let expected = matrix.each_ref().map(|row| row.each_ref().map(eq));
        assert_eq!(
            rot.rational_matrix().unwrap(),
            expected.each_ref().map(|row| row.each_ref())
        );
        for i in 0..3 {
            for j in 0..3 {
                let product: Q = (0..3)
                    .map(|k| {
                        oq(rot.rational_matrix().unwrap()[i][k])
                            * oq(rot.rational_matrix().unwrap()[j][k])
                    })
                    .sum();
                assert_eq!(product, q(i64::from(i == j), 1));
                assert_eq!(
                    rot.inverse().rational_matrix().unwrap()[i][j],
                    rot.rational_matrix().unwrap()[j][i]
                );
            }
        }
        let m = rot.rational_matrix().unwrap().map(|row| row.map(oq));
        assert_eq!(dot(&m[0], &cross(&m[1], &m[2])), q(1, 1));
        let mut dec = Decider::default();
        assert!(
            rot.inverse()
                .same_rotation(
                    &build(&[a[0].clone(), -a[1].clone(), -a[2].clone(), -a[3].clone()]),
                    &mut dec
                )
                .unwrap()
        );
        assert!(
            rot.then(&rot.inverse())
                .same_rotation(&Rotation::identity(), &mut dec)
                .unwrap()
        );
        assert!(
            rot.then(&build(&b))
                .same_rotation(&build(&hamilton(&b, &a)), &mut dec)
                .unwrap()
        );
    }
});

blind_test!(case_10_pose_order_and_vectors, {
    let pose = Pose::new(
        rotation(["1", "0", "0", "1"]),
        v("1", "2", "3").to_real(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(
        pose.apply_point(&p("1", "0", "0").to_real())
            .as_rational()
            .unwrap(),
        p("1", "3", "3")
    );
    assert_eq!(
        pose.inverse()
            .apply_point(&p("1", "3", "3").to_real())
            .as_rational()
            .unwrap(),
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
        point_eq(
            &a.then(&b).apply_point(&point.to_real()),
            &b.apply_point(&a.apply_point(&point.to_real())),
        );
        assert!(
            a.then(&a.inverse())
                .same_motion(&Pose::identity(), &mut Decider::default())
                .unwrap()
        );
        point_eq(
            &a.inverse().apply_point(&a.apply_point(&point.to_real())),
            &point.to_real(),
        );
        let vector = RVector3::new(point.x.clone(), point.y.clone(), point.z.clone());
        real_eq(
            &a.apply_vector(&vector.to_real()).norm_squared(),
            &Real::from(vector.norm_squared()),
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
