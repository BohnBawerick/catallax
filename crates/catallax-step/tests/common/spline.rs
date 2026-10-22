use super::*;
#[derive(Debug)]
pub struct Spline {
    pub degree: usize,
    pub points: Vec<[Rational; 3]>,
    pub knots: Vec<Rational>,
}
impl Spline {
    pub fn read(s: &Step, id: usize) -> Self {
        let e = s.e(id, "B_SPLINE_CURVE_WITH_KNOTS");
        assert_eq!(e.args.len(), 9);
        assert_eq!(e.args[0], "''");
        assert_eq!(e.args[3], ".UNSPECIFIED.");
        assert_eq!(e.args[5], ".F.");
        assert_eq!(e.args[8], ".UNSPECIFIED.");
        let degree: usize = e.args[1].parse().unwrap();
        assert!(degree > 0);
        let cps = list(&e.args[2]);
        let points = cps.iter().map(|p| s.xyz(reference(p))).collect::<Vec<_>>();
        assert!(points.len() > degree);
        let mult = list(&e.args[6])
            .iter()
            .map(|m| m.parse::<usize>().unwrap())
            .collect::<Vec<_>>();
        let knots = list(&e.args[7])
            .iter()
            .map(|k| number(k))
            .collect::<Vec<_>>();
        assert_eq!(knots.len(), mult.len());
        assert!(knots.len() >= 2);
        assert!(knots.windows(2).all(|k| k[0] < k[1]));
        assert!(mult.iter().all(|&m| m > 0));
        assert_eq!(mult[0], degree + 1);
        assert_eq!(*mult.last().unwrap(), degree + 1);
        let knots = knots
            .into_iter()
            .zip(mult)
            .flat_map(|(k, m)| std::iter::repeat_n(k, m))
            .collect::<Vec<_>>();
        assert_eq!(knots.len(), points.len() + degree + 1);
        Self {
            degree,
            points,
            knots,
        }
    }
    pub fn at_span(&self, span: usize, t: &Rational) -> [Rational; 3] {
        // de Boor recurrence in Q; there is no floating-point evaluation.
        let p = self.degree;
        let mut d = self.points[span - p..=span].to_vec();
        for level in 1..=p {
            for j in (level..=p).rev() {
                let i = span - p + j;
                let den = &self.knots[i + p + 1 - level] - &self.knots[i];
                assert!(den > n(0));
                let alpha = (t - &self.knots[i]) / den;
                d[j] = std::array::from_fn(|axis| {
                    &(&n(1) - &alpha) * &d[j - 1][axis] + &alpha * &d[j][axis]
                });
            }
        }
        d[p].clone()
    }
    pub fn samples(&self, mut f: impl FnMut([Rational; 3])) {
        for span in self.degree..self.points.len() {
            let lo = &self.knots[span];
            let hi = &self.knots[span + 1];
            if lo == hi {
                continue;
            }
            // Exactly 1,000 parameters per nonzero knot span, both ends included.
            for j in 0..1000 {
                let t = lo + &((hi - lo) * Rational::from_ratio(j, 999).unwrap());
                f(self.at_span(span, &t));
            }
        }
    }
}
pub fn check_splines(s: &Step, d: &Rational, supports: impl Fn(&[Rational; 3], &Rational)) {
    check_splines_by_curve(s, d, |_, p, d| supports(p, d));
}
pub fn check_splines_by_curve(
    s: &Step,
    d: &Rational,
    supports: impl Fn(usize, &[Rational; 3], &Rational),
) {
    let mut owners = BTreeSet::new();
    let mut points = BTreeSet::new();
    for edge in s.ids("EDGE_CURVE") {
        let e = &s.entities[&edge];
        let id = reference(&e.args[3]);
        if s.entities[&id].kind != "B_SPLINE_CURVE_WITH_KNOTS" {
            continue;
        }
        assert!(owners.insert(id), "one spline per edge, never per CurveId");
        assert_eq!(e.args[4], ".T.");
        let curve = &s.entities[&id];
        let cp = list(&curve.args[2]);
        let start = s.vertex(reference(&e.args[1]));
        let end = s.vertex(reference(&e.args[2]));
        assert_eq!(s.xyz_text(reference(cp[0])), s.xyz_text(start));
        assert_eq!(s.xyz_text(reference(cp[cp.len() - 1])), s.xyz_text(end));
        assert_eq!(
            curve.args[4],
            if e.args[1] == e.args[2] { ".T." } else { ".F." }
        );
        let local_points = cp.into_iter().map(reference).collect::<BTreeSet<_>>();
        for p in local_points {
            assert!(points.insert(p), "control point shared by splines");
            assert!(
                !s.ids("VERTEX_POINT").iter().any(|&v| s.vertex(v) == p),
                "spline owns its points"
            );
        }
        let spline = Spline::read(s, id);
        spline.samples(|p| supports(id, &p, d));
    }
    assert_eq!(owners.len(), s.count("B_SPLINE_CURVE_WITH_KNOTS"));
}
pub fn crossing_cylinders(p: &[Rational; 3], d: &Rational) {
    radial_distance(sq(&p[0]) + sq(&p[1]), &n(2), d);
    radial_distance(sq(&p[1]) + sq(&p[2]), &n(1), d);
}
pub fn offset_sphere(p: &[Rational; 3], d: &Rational) {
    radial_distance(sq(&p[0]) + sq(&p[1]), &n(1), d);
    radial_distance(sq(&(&p[0] - &n(2))) + sq(&p[1]) + sq(&p[2]), &n(2), d);
}

// Exact support distances, independent of the fitted curve and of piece 2's
// implicit-value-to-distance estimates.
pub fn on_exact_support(surface: &Surface, p: &[Rational; 3], d: &Rational) {
    match surface {
        Surface::Plane(s) => {
            let v = sub(p, &real_xyz(&s.origin));
            let normal = real_vec_xyz(&s.normal);
            assert!(sq(&dot(&v, &normal)) <= sq(d) * dot(&normal, &normal));
        }
        Surface::Cylinder(s) => {
            let v = sub(p, &real_xyz(&s.origin));
            let axis = real_vec_xyz(&s.axis);
            let radial2 = dot(&v, &v) - sq(&dot(&v, &axis)) / dot(&axis, &axis);
            radial_distance(radial2, s.radius.rational_ref().unwrap(), d);
        }
        Surface::Sphere(s) => {
            radial_distance(
                distance_sq(p, &real_xyz(&s.center)),
                s.radius.rational_ref().unwrap(),
                d,
            );
        }
        _ => panic!("piece 4 does not define branches on this support"),
    }
}
