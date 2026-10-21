use super::*;

#[derive(Clone, Debug)]
pub struct PairCase {
    pub name: String,
    pub a: Surface,
    pub b: Surface,
}
pub fn pairs() -> Vec<PairCase> {
    let mut out = vec![];
    let mut add = |name: &str, a: Surface, b: Surface| {
        out.push(PairCase {
            name: name.into(),
            a,
            b,
        })
    };
    let o = || pi([0, 0, 0]);
    let z5 = || zc(5);
    let s5 = || sphere(o(), n(5));
    let eps = r("1e-30");
    add(
        "01 PP1",
        plane(pi([0, 0, 1]), vi([0, 0, 2])),
        axis_plane(0, n(3)),
    );
    add(
        "02 PP2",
        axis_plane(2, n(0)),
        plane(p(["0", "0", "1e-30"]), vi([0, 0, -3])),
    );
    add(
        "02 PP3 opposite",
        axis_plane(2, n(0)),
        plane(pi([5, 7, 0]), vi([0, 0, -3])),
    );
    add(
        "02 PP3 same",
        axis_plane(2, n(0)),
        plane(pi([5, 7, 0]), vi([0, 0, 4])),
    );
    add("03 PC1", axis_plane(2, n(2)), z5());
    let tilted = Surface::Plane(Plane::new(o(), vi([0, 3, 4]), vi([1, 0, 0])).unwrap());
    add("04 PC2", tilted, z5());
    add("05 PC3 rational", axis_plane(0, n(3)), z5());
    add("06 PC3 irrational", axis_plane(0, n(1)), z5());
    add("07 PC4", axis_plane(0, n(5)), z5());
    add("07 PC5 near", axis_plane(0, &n(5) + &eps), z5());
    add("07 PC3 near", axis_plane(0, &n(5) - &eps), z5());
    add("08 PC3 seam", axis_plane(1, n(0)), z5());
    add("09 PS1", axis_plane(2, n(3)), s5());
    add("10 PS2 pole", axis_plane(2, n(5)), s5());
    add("10 PS3 near", axis_plane(2, &n(5) + &eps), s5());
    add("10 PS1 near", axis_plane(2, &n(5) - &eps), s5());
    add("11 PS1 seam", axis_plane(1, n(0)), s5());
    add("12 SS3", s5(), sphere(pi([6, 0, 0]), n(5)));
    add("13 SS4", s5(), sphere(pi([10, 0, 0]), n(5)));
    add("13 SS5", s5(), sphere(pi([2, 0, 0]), n(3)));
    add(
        "14 SS6 outside",
        s5(),
        sphere(RPoint3::new(&n(10) + &eps, n(0), n(0)), n(5)),
    );
    add("14 SS6 nested", s5(), sphere(pi([1, 0, 0]), n(3)));
    add("14 SS2", s5(), sphere(o(), n(4)));
    add("14 SS1", s5(), s5());
    add(
        "15 CP3 opposite axes",
        z5(),
        cyl(pi([8, 0, 0]), vi([0, 0, -2]), n(5)),
    );
    add("16 CP4", z5(), cyl(pi([10, 0, 0]), vi([0, 0, 1]), n(5)));
    add("16 CP5", z5(), cyl(pi([2, 0, 0]), vi([0, 0, 1]), n(3)));
    add("16 CP2", z5(), cyl(o(), vi([0, 0, 3]), n(2)));
    add("16 CP1", z5(), cyl(o(), vi([0, 0, 3]), n(5)));
    for (name, y) in [
        ("17 CC5", n(0)),
        ("18 CC3", n(2)),
        ("19 CC4", n(1)),
        ("20 CC2", n(3)),
        ("20 CC1 near", &n(3) + &eps),
        ("20 CC3 near", &n(3) - &eps),
    ] {
        add(
            name,
            zc(2),
            cyl(RPoint3::new(n(0), y, n(0)), vi([1, 0, 0]), n(1)),
        );
    }
    add("21 CC6 right angle", zc(1), cyl(o(), vi([1, 0, 0]), n(1)));
    add(
        "22 CC6 rational ratio",
        zc(1),
        cyl(o(), vi([3, 0, 4]), n(1)),
    );
    add(
        "22 CC6 irrational ratio",
        zc(1),
        cyl(o(), vi([1, 1, 0]), n(1)),
    );
    for (name, radius) in [
        ("23 CS1 rational", 13),
        ("23 CS1 irrational", 6),
        ("23 CS2", 5),
        ("23 CS3", 4),
    ] {
        add(name, z5(), sphere(pi([0, 0, 1]), n(radius)));
    }
    add("24 CS7 Viviani", zc(1), sphere(pi([1, 0, 0]), n(2)));
    add("25 CS6", zc(1), sphere(pi([2, 0, 0]), n(2)));
    add("26 CS8", zc(1), sphere(pi([1, 0, 0]), n(3)));
    add("27 CS5 outside", zc(1), sphere(pi([3, 0, 0]), n(2)));
    add("27 CS5 inside", zc(3), sphere(pi([1, 0, 0]), n(2)));
    for d in [eps.clone(), -&eps] {
        add(
            &format!("27 outside offset {d}"),
            zc(1),
            sphere(RPoint3::new(&n(3) + &d, n(0), n(0)), n(2)),
        );
        // Case 27's original prose reversed this pair. The technical review resolved it:
        // x²+y²=9 and (x-d)²+y²+z²=4 imply z²=-5-d²+2dx.
        // Its maximum for x in [-3,3] is 4-(3-d)². At d=1+eps
        // it is 4eps-eps²>0 (CS6); at 1-eps it is -4eps-eps²<0 (CS4).
        add(
            &format!("27 inside offset {d}"),
            zc(3),
            sphere(RPoint3::new(&n(1) + &d, n(0), n(0)), n(2)),
        );
    }
    // Configurations not assigned their own numbered example.
    add(
        "CP3 irrational same axes",
        z5(),
        cyl(pi([1, 0, 7]), vi([0, 0, 2]), n(5)),
    );
    add("CP6 nested", z5(), cyl(pi([1, 0, 2]), vi([0, 0, -3]), n(2)));
    add(
        "CP6 outside",
        z5(),
        cyl(pi([11, 0, 0]), vi([0, 0, 1]), n(5)),
    );
    add(
        "CC3 equal radii skew",
        zc(1),
        cyl(pi([0, 1, 0]), vi([3, 0, 4]), n(1)),
    );
    add(
        "CC5 oblique nonunit",
        zc(3),
        cyl(pi([0, 1, 7]), vi([3, 0, 4]), n(1)),
    );
    add(
        "CC4 oblique nonunit",
        zc(3),
        cyl(pi([0, 2, 7]), vi([3, 0, 4]), n(1)),
    );
    add("CS6 irrational offset", zc(2), sphere(pi([1, 1, 3]), n(2)));
    add(
        "CC6 both axis lengths irrational ratio rational",
        cyl(o(), vi([1, 1, 0]), n(1)),
        cyl(o(), vi([1, 0, 1]), n(1)),
    );
    for delta in [eps.clone(), -&eps] {
        let shifted = |x: i64| RPoint3::new(&n(x) + &delta, n(0), n(0));
        let label = |base: &str| format!("{base} near {delta}");
        add(&label("SS4"), s5(), sphere(shifted(10), n(5)));
        add(&label("SS5"), s5(), sphere(shifted(2), n(3)));
        add(&label("SS1 radius"), s5(), sphere(o(), &n(5) + &delta));
        add(&label("SS1 centre"), s5(), sphere(shifted(0), n(5)));
        add(&label("CP4"), z5(), cyl(shifted(10), vi([0, 0, 1]), n(5)));
        add(&label("CP5"), z5(), cyl(shifted(2), vi([0, 0, 1]), n(3)));
        add(
            &label("CP1 radius"),
            z5(),
            cyl(o(), vi([0, 0, 1]), &n(5) + &delta),
        );
        add(
            &label("CP1 axis"),
            z5(),
            cyl(shifted(0), vi([0, 0, 1]), n(5)),
        );
        add(
            &label("CC4"),
            zc(2),
            cyl(
                RPoint3::new(n(0), &n(1) + &delta, n(0)),
                vi([1, 0, 0]),
                n(1),
            ),
        );
        add(
            &label("CC6 radius"),
            zc(1),
            cyl(o(), vi([1, 0, 0]), &n(1) + &delta),
        );
        add(
            &label("CC6 skew"),
            zc(1),
            cyl(RPoint3::new(n(0), delta.clone(), n(0)), vi([1, 0, 0]), n(1)),
        );
        add(&label("CS2"), z5(), sphere(o(), &n(5) + &delta));
        add(&label("CS7"), zc(1), sphere(pi([1, 0, 0]), &n(2) + &delta));
        add(
            &label("PP3 angle"),
            axis_plane(2, n(0)),
            plane(o(), RVector3::new(delta.clone(), n(0), n(1))),
        );
        add(
            &label("PC parallel"),
            plane(o(), RVector3::new(n(1), n(0), delta.clone())),
            z5(),
        );
        add(
            &label("PC perpendicular"),
            plane(o(), RVector3::new(n(0), delta, n(1))),
            z5(),
        );
    }
    out
}
