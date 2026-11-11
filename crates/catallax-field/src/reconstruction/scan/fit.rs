//! Least-squares plane, cylinder and cone candidates in binary64.
//!
//! Every fit is a candidate: geometric least squares on the given samples,
//! with no claim about samples not given or about the scanned object.

use super::linalg::{
    V3, add, canonical, dot, eigen, frame, moments, norm, scale, solve, sub, unit,
};

/// Kinds of fitted surface, simplest first.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SurfaceKind {
    /// A plane.
    Plane,
    /// A circular cylinder.
    Cylinder,
    /// A circular cone.
    Cone,
}

/// A fitted analytic surface, in millimetres.
#[derive(Clone, Debug, PartialEq)]
pub enum Surface {
    /// Points `x` with `normal . x = offset`; `normal` is unit length.
    Plane {
        /// Unit normal, in its canonical sign.
        normal: V3,
        /// Signed distance of the plane from the origin along `normal`.
        offset: f64,
    },
    /// Points at distance `radius` from the line through `point` along `axis`.
    Cylinder {
        /// A point on the axis, the foot of the fitted samples' centroid.
        point: V3,
        /// Unit axis direction, in its canonical sign.
        axis: V3,
        /// Radius.
        radius: f64,
    },
    /// Points whose distance from the axis is `radius + slope * t`, where `t`
    /// is the signed axial distance from `point`.
    Cone {
        /// A point on the axis, the foot of the fitted samples' centroid.
        point: V3,
        /// Unit axis direction, in its canonical sign.
        axis: V3,
        /// Radius at `point`.
        radius: f64,
        /// Radius change per unit length along `axis`: the tangent of the half-angle.
        slope: f64,
    },
}

impl Surface {
    /// The kind of surface.
    pub fn kind(&self) -> SurfaceKind {
        match self {
            Surface::Plane { .. } => SurfaceKind::Plane,
            Surface::Cylinder { .. } => SurfaceKind::Cylinder,
            Surface::Cone { .. } => SurfaceKind::Cone,
        }
    }

    /// Signed distance from `p`, positive away from a plane along its normal
    /// and away from the axis for a cylinder or cone. Exact for a plane and a
    /// cylinder; for a cone it is the distance to the generating line, which is
    /// the surface distance away from the apex.
    pub fn residual(&self, p: V3) -> f64 {
        match *self {
            Surface::Plane { normal, offset } => dot(normal, p) - offset,
            Surface::Cylinder {
                point,
                axis,
                radius,
            } => {
                let q = sub(p, point);
                norm(sub(q, scale(axis, dot(q, axis)))) - radius
            }
            Surface::Cone {
                point,
                axis,
                radius,
                slope,
            } => {
                let q = sub(p, point);
                let t = dot(q, axis);
                let rho = norm(sub(q, scale(axis, t)));
                (rho - radius - slope * t) / (1.0 + slope * slope).sqrt()
            }
        }
    }

    /// The axis of a cylinder or cone.
    pub fn axis(&self) -> Option<(V3, V3)> {
        match *self {
            Surface::Plane { .. } => None,
            Surface::Cylinder { point, axis, .. } | Surface::Cone { point, axis, .. } => {
                Some((point, axis))
            }
        }
    }

    /// The radius at axial distance `t` from the axis point.
    pub fn radius_at(&self, t: f64) -> Option<f64> {
        match *self {
            Surface::Plane { .. } => None,
            Surface::Cylinder { radius, .. } => Some(radius),
            Surface::Cone { radius, slope, .. } => Some(radius + slope * t),
        }
    }
}

/// Root-mean-square residual of `surface` over the samples.
pub(super) fn rms(surface: &Surface, points: &[V3], members: &[u32]) -> f64 {
    let sum: f64 = members
        .iter()
        .map(|&i| surface.residual(points[i as usize]).powi(2))
        .sum();
    (sum / members.len().max(1) as f64).sqrt()
}

/// Least-squares plane through the samples, by their covariance.
pub(super) fn plane(points: &[V3], members: &[u32]) -> Option<Surface> {
    if members.len() < 3 {
        return None;
    }
    let (mean, cov, _) = moments(members.iter().map(|&i| points[i as usize]));
    let (values, vectors) = eigen(cov);
    if values[1].is_nan() || values[1] <= 0.0 {
        return None;
    }
    let normal = canonical(unit(vectors[0])?);
    Some(Surface::Plane {
        normal,
        offset: dot(normal, mean),
    })
}

fn centroid(points: &[V3], members: &[u32]) -> V3 {
    moments(members.iter().map(|&i| points[i as usize])).0
}

/// Moves the axis point to the foot of the samples' centroid.
fn recentre(c: V3, a: V3, mean: V3) -> (V3, f64) {
    let shift = dot(sub(mean, c), a);
    (add(c, scale(a, shift)), shift)
}

/// Levenberg-Marquardt on `N` parameters. `eval` returns the residuals and
/// Jacobian rows at the current state; `step` applies an increment.
fn minimise<S: Clone, const N: usize>(
    mut state: S,
    eval: impl Fn(&S, &mut Vec<(f64, [f64; N])>),
    step: impl Fn(&S, &[f64; N]) -> Option<S>,
) -> Option<S> {
    let mut rows = Vec::new();
    let cost = |s: &S, rows: &mut Vec<(f64, [f64; N])>| {
        rows.clear();
        eval(s, rows);
        rows.iter().map(|r| r.0 * r.0).sum::<f64>()
    };
    let mut current = cost(&state, &mut rows);
    if !current.is_finite() {
        return None;
    }
    let mut lambda = 1e-3;
    for _ in 0..200 {
        let mut jtj = [[0.0; N]; N];
        let mut jtf = [0.0; N];
        for (f, j) in &rows {
            for a in 0..N {
                jtf[a] -= j[a] * f;
                for b in 0..N {
                    jtj[a][b] += j[a] * j[b];
                }
            }
        }
        let mut improved = false;
        while lambda < 1e12 {
            let mut damped = jtj;
            for a in 0..N {
                damped[a][a] += lambda * jtj[a][a].max(1e-12);
            }
            let Some(delta) = solve(damped, jtf) else {
                lambda *= 4.0;
                continue;
            };
            let Some(next) = step(&state, &delta) else {
                lambda *= 4.0;
                continue;
            };
            let mut next_rows = Vec::new();
            let next_cost = cost(&next, &mut next_rows);
            if next_cost.is_finite() && next_cost < current {
                let gain = current - next_cost;
                state = next;
                rows = next_rows;
                lambda = (lambda / 3.0).max(1e-12);
                improved = gain > current * 1e-14;
                current = next_cost;
                break;
            }
            lambda *= 4.0;
        }
        if !improved {
            break;
        }
    }
    Some(state)
}

/// Axis estimate from sample normals: the direction they are most nearly
/// perpendicular to.
fn normal_axis(normals: &[V3], members: &[u32]) -> Option<V3> {
    let mut m = [[0.0; 3]; 3];
    for &i in members {
        let n = normals[i as usize];
        for a in 0..3 {
            for b in 0..3 {
                m[a][b] += n[a] * n[b];
            }
        }
    }
    let (_, vectors) = eigen(m);
    unit(vectors[0])
}

/// Algebraic circle through the samples projected along `axis`.
fn circle(points: &[V3], members: &[u32], mean: V3, axis: V3) -> Option<(V3, f64)> {
    let (u, v) = frame(axis);
    let mut m = [[0.0; 3]; 3];
    let mut b = [0.0; 3];
    for &i in members {
        let q = sub(points[i as usize], mean);
        let (x, y) = (dot(q, u), dot(q, v));
        let row = [x, y, 1.0];
        let rhs = -(x * x + y * y);
        for a in 0..3 {
            b[a] += row[a] * rhs;
            for c in 0..3 {
                m[a][c] += row[a] * row[c];
            }
        }
    }
    let [d, e, f] = solve(m, b)?;
    let (cx, cy) = (-d / 2.0, -e / 2.0);
    let r2 = cx * cx + cy * cy - f;
    (r2 > 0.0).then(|| (add(mean, add(scale(u, cx), scale(v, cy))), r2.sqrt()))
}

#[derive(Clone)]
struct Axial {
    c: V3,
    a: V3,
    r: f64,
    k: f64,
}

fn axial_rows<const N: usize>(
    s: &Axial,
    points: &[V3],
    members: &[u32],
    rows: &mut Vec<(f64, [f64; N])>,
) {
    let (u, v) = frame(s.a);
    let w = (1.0 + s.k * s.k).sqrt();
    for &i in members {
        let q = sub(points[i as usize], s.c);
        let t = dot(q, s.a);
        let radial = sub(q, scale(s.a, t));
        let rho = norm(radial).max(1e-300);
        let f = (rho - s.r - s.k * t) / w;
        let (pu, pv) = (dot(radial, u) / rho, dot(radial, v) / rho);
        // Tilting the axis toward u or v, moving its point along u or v,
        // changing the radius, and (for a cone) the slope.
        let mut j = [0.0; N];
        j[0] = (-t * pu - s.k * dot(q, u)) / w;
        j[1] = (-t * pv - s.k * dot(q, v)) / w;
        j[2] = -pu / w;
        j[3] = -pv / w;
        j[4] = -1.0 / w;
        if N > 5 {
            j[5] = -t / w - f * s.k / (w * w);
        }
        rows.push((f, j));
    }
}

fn axial_step<const N: usize>(s: &Axial, d: &[f64; N], mean: V3) -> Option<Axial> {
    let (u, v) = frame(s.a);
    let a = unit(add(s.a, add(scale(u, d[0]), scale(v, d[1]))))?;
    let c = add(s.c, add(scale(u, d[2]), scale(v, d[3])));
    let k = if N > 5 { s.k + d[5] } else { s.k };
    let (c, shift) = recentre(c, a, mean);
    let r = s.r + d[4] + k * shift;
    (r > 0.0 && r.is_finite() && k.is_finite()).then_some(Axial { c, a, r, k })
}

/// Least-squares cylinder, started from the normals' axis and an algebraic circle.
pub(super) fn cylinder(points: &[V3], normals: &[V3], members: &[u32]) -> Option<Surface> {
    if members.len() < 8 {
        return None;
    }
    let mean = centroid(points, members);
    let axis = normal_axis(normals, members)?;
    let (centre, radius) = circle(points, members, mean, axis)?;
    let (c, _) = recentre(centre, axis, mean);
    let start = Axial {
        c,
        a: axis,
        r: radius,
        k: 0.0,
    };
    let fit = minimise::<Axial, 5>(
        start,
        |s, rows| axial_rows(s, points, members, rows),
        |s, d| axial_step(s, d, mean),
    )?;
    let axis = canonical(fit.a);
    Some(Surface::Cylinder {
        point: fit.c,
        axis,
        radius: fit.r,
    })
}

/// Least-squares cone, started from a fitted cylinder and the slope of the
/// samples' radial distance along its axis.
pub(super) fn cone(points: &[V3], normals: &[V3], members: &[u32]) -> Option<Surface> {
    if members.len() < 12 {
        return None;
    }
    let mean = centroid(points, members);
    let Surface::Cylinder {
        point,
        axis,
        radius,
    } = cylinder(points, normals, members)?
    else {
        unreachable!("cylinder fit returns a cylinder");
    };
    let (mut st, mut stt, mut sr, mut str) = (0.0, 0.0, 0.0, 0.0);
    for &i in members {
        let q = sub(points[i as usize], point);
        let t = dot(q, axis);
        let rho = norm(sub(q, scale(axis, t)));
        st += t;
        stt += t * t;
        sr += rho;
        str += t * rho;
    }
    let n = members.len() as f64;
    let denominator = n * stt - st * st;
    let slope = if denominator > 0.0 {
        (n * str - st * sr) / denominator
    } else {
        0.0
    };
    let start = Axial {
        c: point,
        a: axis,
        r: radius,
        k: slope,
    };
    let fit = minimise::<Axial, 6>(
        start,
        |s, rows| axial_rows(s, points, members, rows),
        |s, d| axial_step(s, d, mean),
    )?;
    let (axis, slope) = if dot(canonical(fit.a), fit.a) > 0.0 {
        (fit.a, fit.k)
    } else {
        (scale(fit.a, -1.0), -fit.k)
    };
    Some(Surface::Cone {
        point: fit.c,
        axis,
        radius: fit.r,
        slope,
    })
}
