//! Full-text 2026 core tests. A nearest grid node is only a candidate closest point.
//! Yang and Jia, doi:10.1145/3807948, Eqs. 3-5 and Definition 3.
use super::{Meter, P, Pair, Q, Result, bases, certify_hill, patches, q, r, refuse, shared_chart};
use catallax_numbers::Refusal;

#[derive(Clone)]
struct Jet {
    point: [Q; 3],
    du: [Q; 3],
    dv: [Q; 3],
    normal: [Q; 3],
}
fn dot(a: &[Q; 3], b: &[Q; 3], m: &mut Meter) -> Result<Q> {
    let mut value = q(0);
    for (x, y) in a.iter().zip(b) {
        let xy = m.mul(x, y)?;
        value = m.add(&value, &xy)?;
    }
    Ok(value)
}
fn cross(a: &[Q; 3], b: &[Q; 3], m: &mut Meter) -> Result<[Q; 3]> {
    let mut out = [q(0), q(0), q(0)];
    for (i, x) in out.iter_mut().enumerate() {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        let left = m.mul(&a[j], &b[k])?;
        let right = m.mul(&a[k], &b[j])?;
        *x = m.sub(&left, &right)?;
    }
    Ok(out)
}
fn difference(a: &[Q; 3], b: &[Q; 3], m: &mut Meter) -> Result<[Q; 3]> {
    Ok([
        m.sub(&a[0], &b[0])?,
        m.sub(&a[1], &b[1])?,
        m.sub(&a[2], &b[2])?,
    ])
}
fn eq4(a: &Jet, b: &Jet, m: &mut Meter) -> Result<bool> {
    let delta = difference(&a.point, &b.point, m)?;
    let residuals = [
        dot(&delta, &a.du, m)?,
        dot(&delta, &a.dv, m)?,
        dot(&a.normal, &b.du, m)?,
        dot(&a.normal, &b.dv, m)?,
    ];
    let mut zero = true;
    for residual in residuals {
        zero &= m.cmp(&residual, &q(0))?.is_eq();
    }
    Ok(zero)
}
struct Surface {
    xyz: Vec<[P; 3]>,
    du: Vec<[P; 3]>,
    dv: Vec<[P; 3]>,
}
impl Surface {
    fn new(net: &super::Net, b: &[Vec<P>], m: &mut Meter) -> Result<Self> {
        let xyz = patches(net, b, m)?;
        let mut du = Vec::new();
        let mut dv = Vec::new();
        for p in &xyz {
            du.push([
                p[0].derivative(0, m)?,
                p[1].derivative(0, m)?,
                p[2].derivative(0, m)?,
            ]);
            dv.push([
                p[0].derivative(1, m)?,
                p[1].derivative(1, m)?,
                p[2].derivative(1, m)?,
            ]);
        }
        Ok(Self { xyz, du, dv })
    }
    fn jet(&self, u: &Q, v: &Q, m: &mut Meter) -> Result<Jet> {
        let span = |t: &Q, m: &mut Meter| -> Result<usize> {
            if m.cmp(t, &q(1))?.is_lt() {
                Ok(0)
            } else if m.cmp(t, &q(2))?.is_lt() {
                Ok(1)
            } else {
                Ok(2)
            }
        };
        let (i, j) = (span(u, m)?, span(v, m)?);
        let local_u = m.sub(u, &q(i as i64))?;
        let local_v = m.sub(v, &q(j as i64))?;
        let evaluate = |p: &[P; 3], m: &mut Meter| -> Result<[Q; 3]> {
            Ok([
                p[0].eval(&local_u, &local_v, m)?,
                p[1].eval(&local_u, &local_v, m)?,
                p[2].eval(&local_u, &local_v, m)?,
            ])
        };
        let point = evaluate(&self.xyz[3 * i + j], m)?;
        let du = evaluate(&self.du[3 * i + j], m)?;
        let dv = evaluate(&self.dv[3 * i + j], m)?;
        let normal = cross(&du, &dv, m)?;
        let square = dot(&normal, &normal, m)?;
        if !m.cmp(&square, &q(0))?.is_gt() {
            return Err(refuse("r13.degenerate-distance-normal"));
        }
        Ok(Jet {
            point,
            du,
            dv,
            normal,
        })
    }
}

// Equation 5, including its asymmetric treatment of coordinate axes.
fn quadrant(x: &Q, y: &Q, m: &mut Meter) -> Result<u8> {
    let (x, y) = (m.cmp(x, &q(0))?, m.cmp(y, &q(0))?);
    Ok(if x.is_eq() && y.is_eq() {
        0
    } else if !x.is_lt() && y.is_gt() {
        1
    } else if x.is_lt() && !y.is_lt() {
        2
    } else if !x.is_gt() && y.is_lt() {
        3
    } else {
        4
    })
}
#[derive(Clone, Copy)]
enum EdgeRule {
    Printed,
    Cyclic,
}
fn valid_edge(a: u8, b: u8, rule: EdgeRule) -> bool {
    if a == 0 || b == 0 {
        return true;
    }
    let d = a.abs_diff(b);
    match rule {
        EdgeRule::Printed => d % 4 <= 1,
        EdgeRule::Cyclic => d.min(4 - d) <= 1,
    }
}
fn cell_index(corners: [u8; 4]) -> i32 {
    for offset in 0..4 {
        if (0..4).all(|i| corners[i] == [1, 2, 3, 4][(i + offset) % 4]) {
            return 1;
        }
        if (0..4).all(|i| corners[i] == [1, 4, 3, 2][(i + offset) % 4]) {
            return -1;
        }
    }
    0
}
#[derive(Clone)]
struct Sample {
    quadrant: u8,
    tied_targets: usize,
    boundary_target: bool,
    eq4_zero: bool,
}
struct Target {
    u: Q,
    v: Q,
    jet: Jet,
}
struct Grid<'a> {
    source: &'a Surface,
    targets: Vec<Target>,
    cache: Vec<(Q, Q, Sample)>,
}
impl<'a> Grid<'a> {
    fn new(source: &'a Surface, target: &Surface, m: &mut Meter) -> Result<Self> {
        let mut targets = Vec::new();
        for i in 0..=4 {
            for j in 0..=4 {
                let (u, v) = (m.div(&q(3 * i), &q(4))?, m.div(&q(3 * j), &q(4))?);
                let jet = target.jet(&u, &v, m)?;
                targets.push(Target { u, v, jet });
            }
        }
        Ok(Self {
            source,
            targets,
            cache: Vec::new(),
        })
    }
    fn sample(&mut self, u: &Q, v: &Q, m: &mut Meter) -> Result<Sample> {
        // A small linear cache makes every rational key comparison chargeable.
        for (cu, cv, sample) in &self.cache {
            if m.cmp(cu, u)?.is_eq() && m.cmp(cv, v)?.is_eq() {
                return Ok(sample.clone());
            }
        }
        let a = self.source.jet(u, v, m)?;
        let mut best: Option<(usize, Q)> = None;
        let mut tied_targets = 0;
        for (i, target) in self.targets.iter().enumerate() {
            let delta = difference(&a.point, &target.jet.point, m)?;
            let square = dot(&delta, &delta, m)?;
            let order = match &best {
                Some((_, old)) => m.cmp(&square, old)?,
                None => std::cmp::Ordering::Less,
            };
            if order.is_lt() {
                best = Some((i, square));
                tied_targets = 1;
            } else if order.is_eq() {
                tied_targets += 1;
            }
        }
        let target = &self.targets[best.unwrap().0];
        // Dropping the positive normal length preserves quadrants and zeros.
        let x = dot(&target.jet.normal, &a.du, m)?;
        let y = dot(&target.jet.normal, &a.dv, m)?;
        let quadrant = quadrant(&x, &y, m)?;
        let mut boundary_target = false;
        for t in [&target.u, &target.v] {
            boundary_target |= m.cmp(t, &q(0))?.is_eq() || m.cmp(t, &q(3))?.is_eq();
        }
        let eq4_zero = quadrant == 0 && eq4(&a, &target.jet, m)?;
        let sample = Sample {
            quadrant,
            tied_targets,
            boundary_target,
            eq4_zero,
        };
        self.cache.push((u.clone(), v.clone(), sample.clone()));
        Ok(sample)
    }
    fn cell(
        &mut self,
        u: [&Q; 2],
        v: [&Q; 2],
        depth: usize,
        rule: EdgeRule,
        stats: &mut Stats,
        m: &mut Meter,
    ) -> Result<()> {
        let corners = [
            self.sample(u[0], v[0], m)?.quadrant,
            self.sample(u[1], v[0], m)?.quadrant,
            self.sample(u[1], v[1], m)?.quadrant,
            self.sample(u[0], v[1], m)?.quadrant,
        ];
        let invalid: [bool; 4] =
            std::array::from_fn(|i| !valid_edge(corners[i], corners[(i + 1) % 4], rule));
        if invalid.iter().all(|x| !x) || depth == 2 {
            stats.leaves += 1;
            stats.invalid_edges += invalid.iter().filter(|x| **x).count();
            if invalid.iter().all(|x| !x) {
                match cell_index(corners) {
                    1 => stats.positive += 1,
                    -1 => stats.negative += 1,
                    _ => (),
                }
            }
            return Ok(());
        }
        stats.splits += 1;
        let um = m.add(u[0], u[1])?;
        let um = m.div(&um, &q(2))?;
        let vm = m.add(v[0], v[1])?;
        let vm = m.div(&vm, &q(2))?;
        let adjacent = (0..4).any(|i| invalid[i] && invalid[(i + 1) % 4]);
        let us = if adjacent || invalid[0] || invalid[2] {
            vec![[u[0], &um], [&um, u[1]]]
        } else {
            vec![u]
        };
        let vs = if adjacent || invalid[1] || invalid[3] {
            vec![[v[0], &vm], [&vm, v[1]]]
        } else {
            vec![v]
        };
        for uu in us {
            for vv in &vs {
                self.cell(uu, *vv, depth + 1, rule, stats, m)?;
            }
        }
        Ok(())
    }
    fn run(&mut self, rule: EdgeRule, m: &mut Meter) -> Result<Stats> {
        let mut stats = Stats::default();
        for i in 0..4 {
            for j in 0..4 {
                let u = [m.div(&q(3 * i), &q(4))?, m.div(&q(3 * (i + 1)), &q(4))?];
                let v = [m.div(&q(3 * j), &q(4))?, m.div(&q(3 * (j + 1)), &q(4))?];
                self.cell([&u[0], &u[1]], [&v[0], &v[1]], 0, rule, &mut stats, m)?;
            }
        }
        Ok(stats)
    }
}
#[derive(Default)]
struct Stats {
    leaves: usize,
    invalid_edges: usize,
    positive: usize,
    negative: usize,
    splits: usize,
}

fn level(h: &[P], t: &Q, c: &Q, m: &mut Meter) -> Result<Q> {
    let span = if m.cmp(t, &q(1))?.is_lt() {
        0
    } else if m.cmp(t, &q(2))?.is_lt() {
        1
    } else {
        2
    };
    let local = m.sub(t, &q(span as i64))?;
    let height = h[span].eval(&local, &q(0), m)?;
    let product = m.mul(&r("15/16"), &height)?;
    m.sub(&product, c)
}
fn root_box(h: &[P], mut lo: Q, mut hi: Q, c: &Q, m: &mut Meter) -> Result<String> {
    let left = level(h, &lo, c, m)?;
    let left_sign = m.cmp(&left, &q(0))?;
    let right = level(h, &hi, c, m)?;
    let right_sign = m.cmp(&right, &q(0))?;
    if left_sign.is_eq() || right_sign.is_eq() || left_sign == right_sign {
        return Err(refuse("r13.start-not-bracketed"));
    }
    // Fixed isolation rounds, not a geometric tolerance or an acceptance threshold.
    for _ in 0..12 {
        let sum = m.add(&lo, &hi)?;
        let mid = m.div(&sum, &q(2))?;
        let value = level(h, &mid, c, m)?;
        let sign = m.cmp(&value, &q(0))?;
        if sign.is_eq() {
            return Ok(format!("[{mid},{mid}]"));
        }
        if sign == left_sign {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(format!("[{lo},{hi}]"))
}
fn peak_witness(
    pair: &Pair,
    b: &[Vec<P>],
    a: &Surface,
    target: &Surface,
    m: &mut Meter,
) -> Result<String> {
    certify_hill(b, m)?;
    let c = m.mul(&q(2), &pair[1][0][0][2])?;
    for (i, (a_row, b_row)) in pair[0].iter().zip(&pair[1]).enumerate() {
        for (j, (a_point, b_point)) in a_row.iter().zip(b_row).enumerate() {
            let height = q(i64::from((2..=3).contains(&i) && (2..=3).contains(&j)));
            m.equal(&a_point[2], &height)?;
            let z = m.add(&height, &c)?;
            let z = m.div(&z, &q(2))?;
            m.equal(&b_point[2], &z)?;
        }
    }
    let mut h = Vec::new();
    let mut curvature = q(0);
    let mut min_speed = q(100);
    for span in b {
        let height = span[2].add(&span[3], m)?;
        let range = height.derivative(0, m)?.derivative(0, m)?.bounds(m)?;
        for x in range {
            if m.cmp(&x.abs(), &curvature)?.is_gt() {
                curvature = x.abs();
            }
        }
        let mut x = P::zero(1, 1);
        for (i, p) in span.iter().enumerate() {
            x = x.add(&p.scale(&q(i as i64), m)?, m)?;
        }
        let speed = x.derivative(0, m)?.bounds(m)?;
        if m.cmp(&speed[0], &min_speed)?.is_lt() {
            min_speed = speed[0].clone();
        }
        h.push(height);
    }
    if !m.cmp(&min_speed, &q(0))?.is_gt() {
        return Err(refuse("r13.graph-speed-not-positive"));
    }
    let speed_square = m.mul(&min_speed, &min_speed)?;
    let denominator = m.mul(&q(4), &speed_square)?;
    let numerator = m.mul(&r("15/16"), &curvature)?;
    let k = m.div(&numerator, &denominator)?;
    let p = a.jet(&r("3/2"), &r("3/2"), m)?;
    let closest = target.jet(&r("3/2"), &r("3/2"), m)?;
    if !eq4(&p, &closest, m)? {
        return Err(refuse("r13.critical-residual-not-zero"));
    }
    for jet in [&p, &closest] {
        m.equal(&jet.du[2], &q(0))?;
        m.equal(&jet.dv[2], &q(0))?;
        if !m.cmp(&jet.normal[2], &q(0))?.is_gt() {
            return Err(refuse("r13.graph-normal-not-positive"));
        }
    }
    let gap = m.sub(&p.point[2], &closest.point[2])?;
    let e = m.sub(&q(0), &gap)?;
    let margin = if m.cmp(&e, &q(0))?.is_gt() {
        let ek = m.mul(&e, &k)?;
        let twice = m.mul(&q(2), &ek)?;
        m.sub(&q(1), &twice)?
    } else {
        q(1)
    };
    if !m.cmp(&margin, &q(0))?.is_gt() {
        return Err(refuse("r13.closest-peak-not-certified"));
    }
    // Height drop <= K * horizontal distance squared gives
    // distance^2 - peak_distance^2 >= margin * horizontal distance squared.
    let starts = match (m.cmp(&c, &q(0))?, m.cmp(&c, &r("225/256"))?) {
        (_, std::cmp::Ordering::Equal) => "starts=1;contact=exact".into(),
        (std::cmp::Ordering::Greater, std::cmp::Ordering::Less) => {
            let left = root_box(&h, q(0), r("3/2"), &c, m)?;
            let right = root_box(&h, r("3/2"), q(3), &c, m)?;
            format!("starts=4;root-boxes={left},{right}")
        }
        (_, std::cmp::Ordering::Greater) => "starts=0;iso-level=above-maximum".into(),
        _ => return Err(refuse("r13.start-family-not-covered")),
    };
    Ok(format!(
        "peak=exact;closest=unique;gap={gap};nearest-margin={margin};{starts}"
    ))
}

pub(super) fn experiment(pair: &Pair, m: &mut Meter) -> Result<String> {
    shared_chart(pair, m)?;
    let b = bases(m)?;
    let source = Surface::new(&pair[0], &b, m)?;
    let target = Surface::new(&pair[1], &b, m)?;
    let mut grid = Grid::new(&source, &target, m)?;
    let printed = grid.run(EdgeRule::Printed, m)?;
    let cyclic = grid.run(EdgeRule::Cyclic, m)?;
    let ties = grid
        .cache
        .iter()
        .filter(|(_, _, x)| x.tied_targets > 1)
        .count();
    let boundary = grid
        .cache
        .iter()
        .filter(|(_, _, x)| x.boundary_target)
        .count();
    let exact_grid = grid.cache.iter().filter(|(_, _, x)| x.eq4_zero).count();
    let peak = match peak_witness(pair, &b, &source, &target, m) {
        Ok(proof) => proof,
        Err(refusal) if refusal.code.starts_with("numbers.") => return Err(refusal),
        Err(refusal) => format!("peak=Refused({})", refusal.code),
    };
    Err(Refusal::new(
        "r13.distance-field-coverage-unproved",
        format!(
            "samples={};tied-grid-nearest={ties};boundary-grid-nearest={boundary};grid-eq4-zeros={exact_grid};printed-leaves={};printed-invalid={};cyclic-leaves={};cyclic-invalid={};cyclic-index-plus={};cyclic-index-minus={};cyclic-splits={};{peak}",
            grid.cache.len(),
            printed.leaves,
            printed.invalid_edges,
            cyclic.leaves,
            cyclic.invalid_edges,
            cyclic.positive,
            cyclic.negative,
            cyclic.splits
        ),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn quadrant_axes_and_published_edge_discrepancy() {
        let mut m = Meter::default();
        for (x, y, expected) in [(0, 0, 0), (0, 1, 1), (-1, 0, 2), (0, -1, 3), (1, 0, 4)] {
            assert_eq!(quadrant(&q(x), &q(y), &mut m).unwrap(), expected);
        }
        assert!(!valid_edge(4, 1, EdgeRule::Printed));
        assert!(valid_edge(4, 1, EdgeRule::Cyclic));
        assert!(!valid_edge(1, 3, EdgeRule::Cyclic));
        assert_eq!(cell_index([1, 2, 3, 4]), 1);
        assert_eq!(cell_index([1, 4, 3, 2]), -1);
        assert_eq!(cell_index([1, 2, 1, 2]), 0);
    }
    #[test]
    fn critical_equations_reject_nonperpendicular_distance() {
        let mut m = Meter::default();
        let a = Jet {
            point: [q(0), q(0), q(1)],
            du: [q(1), q(0), q(0)],
            dv: [q(0), q(1), q(0)],
            normal: [q(0), q(0), q(1)],
        };
        let mut b = a.clone();
        b.point[2] = q(0);
        assert!(eq4(&a, &b, &mut m).unwrap());
        b.point[0] = q(1);
        assert!(!eq4(&a, &b, &mut m).unwrap());
    }
    #[test]
    fn stationary_pair_can_fail_the_nearest_point_test() {
        let mut pair = super::super::fixtures()
            .into_iter()
            .find(|c| c.name == "loop/base")
            .unwrap()
            .pair;
        for (a_row, b_row) in pair[0].clone().iter().zip(&mut pair[1]) {
            for (a, b) in a_row.iter().zip(b_row) {
                b[2] = (&a[2] + &q(100)) / q(2);
            }
        }
        let mut m = Meter::default();
        let b = bases(&mut m).unwrap();
        let a = Surface::new(&pair[0], &b, &mut m).unwrap();
        let target = Surface::new(&pair[1], &b, &mut m).unwrap();
        let source_peak = a.jet(&r("3/2"), &r("3/2"), &mut m).unwrap();
        let target_peak = target.jet(&r("3/2"), &r("3/2"), &mut m).unwrap();
        let target_boundary = target.jet(&q(0), &r("3/2"), &mut m).unwrap();
        assert!(eq4(&source_peak, &target_peak, &mut m).unwrap());
        let peak_delta = difference(&source_peak.point, &target_peak.point, &mut m).unwrap();
        let boundary_delta =
            difference(&source_peak.point, &target_boundary.point, &mut m).unwrap();
        let peak_square = dot(&peak_delta, &peak_delta, &mut m).unwrap();
        let boundary_square = dot(&boundary_delta, &boundary_delta, &mut m).unwrap();
        assert!(boundary_square < peak_square);
        assert_eq!(
            peak_witness(&pair, &b, &a, &target, &mut m)
                .unwrap_err()
                .code,
            "r13.closest-peak-not-certified"
        );
    }
}
