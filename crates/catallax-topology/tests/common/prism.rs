use super::*;

pub fn prism(n: usize) -> Expected {
    assert!(n >= 3);
    // Increasing samples of y=x^2 followed by the closing chord are convex and
    // counterclockwise. Each turn along the parabola has positive determinant.
    let mut d = Expected {
        points: Vec::with_capacity(2 * n),
        edges: Vec::with_capacity(3 * n),
        faces: Vec::with_capacity(n + 2),
        interior: p(0, 0, 0),
        chi: 2,
    };
    for z in [0, 1] {
        for i in 0..n {
            let i = i as i64;
            d.points.push(p(i, i * i, z));
        }
    }
    for layer in 0..2 {
        for i in 0..n {
            let a = layer * n + i;
            let b = layer * n + (i + 1) % n;
            d.edges.push((
                Curve::Line(Line {
                    origin: d.points[a].clone(),
                    dir: d.points[b].minus(&d.points[a]),
                }),
                Some(a),
                Some(b),
            ));
        }
    }
    for i in 0..n {
        d.edges.push((
            Curve::Line(Line {
                origin: d.points[i].clone(),
                dir: v(0, 0, 1),
            }),
            Some(i),
            Some(n + i),
        ));
    }
    d.faces.push((
        plane(p(0, 0, 0), v(0, 0, -1), v(1, 0, 0)),
        Sense::Same,
        vec![Boundary::Cycle(
            (0..n).rev().map(|i| (i, Sense::Reversed)).collect(),
        )],
    ));
    d.faces.push((
        plane(p(0, 0, 1), v(0, 0, 1), v(1, 0, 0)),
        Sense::Same,
        vec![Boundary::Cycle(
            (0..n).map(|i| (n + i, Sense::Same)).collect(),
        )],
    ));
    for i in 0..n {
        let j = (i + 1) % n;
        let direction = d.points[j].minus(&d.points[i]);
        d.faces.push((
            plane(
                d.points[i].clone(),
                direction.cross(&v(0, 0, 1)),
                v(0, 0, 1),
            ),
            Sense::Same,
            vec![Boundary::Cycle(vec![
                (i, Sense::Same),
                (2 * n + j, Sense::Same),
                (n + i, Sense::Reversed),
                (2 * n + i, Sense::Reversed),
            ])],
        ));
    }
    // Mean of the polygon vertices lies strictly inside this convex polygon.
    let sum = d.points[..n].iter().fold(RVector3::zero(), |s, p| {
        s.plus(&p.minus(&RPoint3::origin()))
    });
    d.interior = RPoint3::origin()
        .plus(&sum.scaled(&q(1, n as i64)))
        .plus(&RVector3::new(r(0), r(0), q(1, 2)));
    d
}
