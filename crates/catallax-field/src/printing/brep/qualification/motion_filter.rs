//! A filter for the existing swept separating-axis test. Every direction is a
//! proposal. Outward intervals for relative positions certify the whole move;
//! an uncertain pair goes through the Bernstein test.

use catallax_numbers::Interval;

type IV = [Interval; 3];
type V = [f64; 3];

pub(super) struct Filter {
    start: Vec<IV>,
    end: Vec<IV>,
}

fn point(x: f64) -> Interval {
    Interval::point(x).expect("finite filter coordinate")
}
fn sub(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] - b[i])
}
fn dot(a: V, b: V) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn cross(a: V, b: V) -> V {
    std::array::from_fn(|i| {
        let [j, k] = [(i + 1) % 3, (i + 2) % 3];
        a[j] * b[k] - a[k] * b[j]
    })
}
fn cross_i(a: IV, b: IV) -> IV {
    std::array::from_fn(|i| {
        let [j, k] = [(i + 1) % 3, (i + 2) % 3];
        a[j] * b[k] - a[k] * b[j]
    })
}
fn dot_i(a: IV, b: IV) -> Interval {
    (0..3).fold(point(0.0), |v, i| v + a[i] * b[i])
}
fn unit(v: V) -> V {
    let length = crate::metric::software_sqrt(dot(v, v));
    if length > 0.0 {
        v.map(|x| x / length)
    } else {
        v
    }
}

impl Filter {
    pub(super) fn new(start: Vec<IV>, end: &[[f32; 3]]) -> Self {
        let end = end.iter().map(|p| p.map(|x| point(f64::from(x)))).collect();
        Self { start, end }
    }

    fn difference(&self, a: usize, b: usize) -> IV {
        std::array::from_fn(|i| {
            (self.start[a][i] - self.start[b][i]).hull(&(self.end[a][i] - self.end[b][i]))
        })
    }

    pub(super) fn orientation(&self, [a, b, c]: [usize; 3]) -> bool {
        let initial = |v: usize| std::array::from_fn(|i| self.start[v][i] - self.start[a][i]);
        let start = cross_i(initial(b), initial(c));
        let swept = cross_i(self.difference(b, a), self.difference(c, a));
        dot_i(start, swept).lo() > 0.0
    }

    pub(super) fn separated(&self, a: [usize; 3], b: [usize; 3]) -> bool {
        let mut shared = [0; 3];
        let mut count = 0;
        for i in a {
            if b.contains(&i) {
                shared[count] = i;
                count += 1;
            }
        }
        match count {
            2 => {
                let [p, q] = [shared[0], shared[1]];
                let edge = self.difference(q, p);
                let u = self.difference(a.into_iter().find(|&i| i != p && i != q).unwrap(), p);
                let v = self.difference(b.into_iter().find(|&i| i != p && i != q).unwrap(), p);
                let normal = cross_i(edge, u);
                // Non-coplanar planes meet only on the common line. For
                // coplanar faces, opposite projections across it suffice.
                !dot_i(normal, v).contains_zero() || dot_i(normal, cross_i(edge, v)).hi() < 0.0
            }
            1 => self.separate_wedges(a, b, shared[0]),
            0 => self.separate_disjoint(a, b),
            _ => false,
        }
    }

    fn separate_wedges(&self, a: [usize; 3], b: [usize; 3], anchor: usize) -> bool {
        // Both wedges are based at the same moving vertex. Negate the second
        // wedge's rays to look for one direction strictly positive on all four.
        let mut rays = a
            .into_iter()
            .filter(|i| *i != anchor)
            .map(|i| self.difference(i, anchor))
            .chain(
                b.into_iter()
                    .filter(|i| *i != anchor)
                    .map(|i| self.difference(anchor, i)),
            );
        let rays: [IV; 4] = std::array::from_fn(|_| rays.next().unwrap());
        let guesses = rays.map(|v| unit(v.map(|v| v.mid())));
        let test = |dir: V| {
            dir.iter().all(|v| v.is_finite())
                && rays
                    .iter()
                    .all(|ray| dot_i(*ray, dir.map(point)).lo() > 0.0)
        };
        let centre = std::array::from_fn(|i| guesses.iter().map(|v| v[i]).sum());
        if test(centre) || guesses.iter().any(|v| test(*v)) {
            return true;
        }

        // If the origin is outside the convex hull of the unit rays, the
        // nearest point is on a vertex, edge or triangle of that hull. Test
        // the projections on those spans without trusting their rounding.
        for i in 0..4 {
            let u = guesses[i];
            for j in i + 1..4 {
                let e = sub(guesses[j], u);
                let length = dot(e, e);
                if length > 0.0 {
                    let t = (-dot(u, e) / length).clamp(0.0, 1.0);
                    if test(std::array::from_fn(|k| u[k] + t * e[k])) {
                        return true;
                    }
                }
                for k in j + 1..4 {
                    let n = cross(e, sub(guesses[k], u));
                    let sense = dot(n, u);
                    if test(n.map(|v| v * sense)) {
                        return true;
                    }
                }
            }
        }
        false
    }

    fn separate_disjoint(&self, a: [usize; 3], b: [usize; 3]) -> bool {
        let av = a.map(|i| self.start[i].map(|v| v.mid()));
        let bv = b.map(|i| self.start[i].map(|v| v.mid()));
        let edges = |v: [V; 3]| std::array::from_fn::<_, 3, _>(|i| sub(v[(i + 1) % 3], v[i]));
        let ae = edges(av);
        let be = edges(bv);
        let an = cross(ae[0], ae[1]);
        let bn = cross(be[0], be[1]);
        let relative: [[IV; 3]; 3] = a.map(|i| b.map(|j| self.difference(j, i)));
        let test = |u: V| {
            if !u.iter().all(|v| v.is_finite()) {
                return false;
            }
            let axis = u.map(point);
            let dots = relative.map(|row| row.map(|d| dot_i(d, axis)));
            dots.iter().flatten().all(|v| v.lo() > 0.0)
                || dots.iter().flatten().all(|v| v.hi() < 0.0)
        };
        let delta = std::array::from_fn(|i| (0..3).map(|j| bv[j][i] - av[j][i]).sum());
        test(delta)
            || test(an)
            || test(bn)
            || ae.iter().any(|a| be.iter().any(|b| test(cross(*a, *b))))
            || ae.iter().any(|e| test(cross(an, *e)))
            || be.iter().any(|e| test(cross(bn, *e)))
    }
}
