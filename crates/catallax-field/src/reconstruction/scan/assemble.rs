//! From fitted surfaces to an ordinary recipe: a block stock from three
//! perpendicular pairs of planes, cylinders and cones parallel to a stock axis
//! as bosses (with their cap planes) or bores, and cylinders on stock edges as
//! convex fillets. Any fitted surface this family cannot place is returned as
//! unplaced, with its reason, never dropped.

use catallax_numbers::{Decider, Rational, Refusal};
use catallax_recipe::{
    BooleanKind, EdgeRef, EdgeSelection, ElementId, Number, Operation, OriginPattern, Pose, Recipe,
    Role, Step, StepId, Vector3,
};

use super::fit::Surface;
use super::linalg::{V3, cross, dot, norm, scale, sub, unit};
use super::{SurfaceUse, Unplaced};

/// Planes or axes closer than this to parallel or perpendicular are taken as such.
const ANGLE_DEGREES: f64 = 1.0;
/// Lengths in the recipe are rounded to this many millimetres.
const LENGTH_GRID: i64 = 1000;
/// Through cuts run this far past the faces they open, so that no cutting
/// face coincides with a stock face. Clamped to half the stock's extent.
const OVERLAP_MM: f64 = 1.0;
/// Quaternion components in the recipe are rounded to this grid.
const ROTATION_GRID: i64 = 1_000_000;

fn id(s: &str) -> StepId {
    StepId::new(s).expect("static step id")
}

fn grid(x: f64) -> Rational {
    Rational::from_ratio((x * LENGTH_GRID as f64).round() as i64, LENGTH_GRID)
        .expect("nonzero grid")
}

fn length(x: f64) -> Number {
    Number::new(grid(x))
}

fn quaternion_number(x: f64) -> Number {
    Number::new(
        Rational::from_ratio((x * ROTATION_GRID as f64).round() as i64, ROTATION_GRID)
            .expect("nonzero grid"),
    )
}

/// The rotation's quaternion `[w, x, y, z]` from its matrix columns.
fn quaternion(cols: [V3; 3]) -> [f64; 4] {
    let m = |r: usize, c: usize| cols[c][r];
    let trace = m(0, 0) + m(1, 1) + m(2, 2);
    let q = if trace > 0.0 {
        let s = (trace + 1.0).sqrt() * 2.0;
        [
            0.25 * s,
            (m(2, 1) - m(1, 2)) / s,
            (m(0, 2) - m(2, 0)) / s,
            (m(1, 0) - m(0, 1)) / s,
        ]
    } else if m(0, 0) > m(1, 1) && m(0, 0) > m(2, 2) {
        let s = (1.0 + m(0, 0) - m(1, 1) - m(2, 2)).sqrt() * 2.0;
        [
            (m(2, 1) - m(1, 2)) / s,
            0.25 * s,
            (m(0, 1) + m(1, 0)) / s,
            (m(0, 2) + m(2, 0)) / s,
        ]
    } else if m(1, 1) > m(2, 2) {
        let s = (1.0 + m(1, 1) - m(0, 0) - m(2, 2)).sqrt() * 2.0;
        [
            (m(0, 2) - m(2, 0)) / s,
            (m(0, 1) + m(1, 0)) / s,
            0.25 * s,
            (m(1, 2) + m(2, 1)) / s,
        ]
    } else {
        let s = (1.0 + m(2, 2) - m(0, 0) - m(1, 1)).sqrt() * 2.0;
        [
            (m(1, 0) - m(0, 1)) / s,
            (m(0, 2) + m(2, 0)) / s,
            (m(1, 2) + m(2, 1)) / s,
            0.25 * s,
        ]
    };
    if q[0] < 0.0 { q.map(|x| -x) } else { q }
}

/// Columns of the rotation of a quaternion of any nonzero length.
fn columns(q: [f64; 4]) -> [V3; 3] {
    let [w, x, y, z] = q;
    let n = w * w + x * x + y * y + z * z;
    [
        [
            (w * w + x * x - y * y - z * z) / n,
            2.0 * (x * y + w * z) / n,
            2.0 * (x * z - w * y) / n,
        ],
        [
            2.0 * (x * y - w * z) / n,
            (w * w - x * x + y * y - z * z) / n,
            2.0 * (y * z + w * x) / n,
        ],
        [
            2.0 * (x * z + w * y) / n,
            2.0 * (y * z - w * x) / n,
            (w * w - x * x - y * y + z * z) / n,
        ],
    ]
}

/// Exact quaternion product `a b`.
fn compose(a: &[Number; 4], b: [i64; 4]) -> [Number; 4] {
    let a = a.clone().map(|n| n.value().clone());
    let b = b.map(Rational::from);
    let m = |i: usize, j: usize| &a[i] * &b[j];
    [
        &(&m(0, 0) - &m(1, 1)) - &(&m(2, 2) + &m(3, 3)),
        &(&m(0, 1) + &m(1, 0)) + &(&m(2, 3) - &m(3, 2)),
        &(&m(0, 2) - &m(1, 3)) + &(&m(2, 0) + &m(3, 1)),
        &(&m(0, 3) + &m(1, 2)) + &(&m(3, 0) - &m(2, 1)),
    ]
    .map(Number::new)
}

/// The quaternion that turns local `+z` onto stock axis `k`.
const TURN_TO: [[i64; 4]; 3] = [[1, 0, 1, 0], [1, -1, 0, 0], [1, 0, 0, 0]];

fn side_role(axis: usize, high: bool) -> Role {
    match (axis, high) {
        (2, false) => Role::Start,
        (2, true) => Role::End,
        _ => Role::Side {
            segment: ElementId::new(match (axis, high) {
                (0, false) => "x_min",
                (0, true) => "x_max",
                (1, false) => "y_min",
                _ => "y_max",
            })
            .expect("static segment id"),
            path: None,
        },
    }
}

/// A stock edge, as two (axis, high side) faces.
type Edge = [(usize, bool); 2];
/// A round's radius, its sample weight or surface, and its stock edge.
type RoundedEdge = (f64, usize, Edge);

/// Whether the samples go more than half way round the axis, as a bore's
/// or a boss's do and a round's never does.
fn wraps(
    surface: &Surface,
    members: &[u32],
    points: &[V3],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let Some((point, axis)) = surface.axis() else {
        return Ok(false);
    };
    let (u, v) = super::linalg::frame(axis);
    let mut angles: Vec<f64> = members
        .iter()
        .map(|&i| {
            dec.charge_work(1)?;
            let q = sub(points[i as usize], point);
            Ok(dot(q, v).atan2(dot(q, u)))
        })
        .collect::<Result<_, Refusal>>()?;
    if angles.len() < 3 || angles.iter().any(|a| !a.is_finite()) {
        return Ok(false);
    }
    super::charge_sort(angles.len(), dec)?;
    angles.sort_by(f64::total_cmp);
    dec.check_deadline()?;
    dec.charge_work(angles.len() as u64)?;
    let seam_gap = angles[0] + std::f64::consts::TAU - angles[angles.len() - 1];
    let largest_gap = angles
        .windows(2)
        .map(|a| a[1] - a[0])
        .fold(seam_gap, f64::max);
    dec.check_deadline()?;
    Ok(largest_gap < std::f64::consts::PI)
}

/// What the assembly made.
pub(super) struct Assembly {
    pub recipe: Recipe,
    pub uses: Vec<SurfaceUse>,
    pub unplaced: Vec<Unplaced>,
    /// Fillet cylinders refitted with their axis held one radius inside both
    /// stock faces, as the recipe's fillet defines them.
    pub refined: Vec<(usize, Surface)>,
}

/// The radius of a round tangent to two stock faces that best fits the
/// samples: the axis sits `r` inside both faces, so only `r` is free.
fn tangent_radius(
    frame: &Frame,
    points: &[V3],
    members: &[u32],
    edge: Edge,
    start: f64,
    dec: &mut Decider,
) -> Result<f64, Refusal> {
    let corner = edge.map(|(o, high)| if high { frame.hi[o] } else { frame.lo[o] });
    let inward = edge.map(|(_, high)| if high { -1.0 } else { 1.0 });
    let local: Vec<[f64; 2]> = members
        .iter()
        .map(|&i| {
            dec.charge_work(1)?;
            let l = frame.local(points[i as usize]);
            Ok([l[edge[0].0], l[edge[1].0]])
        })
        .collect::<Result<_, Refusal>>()?;
    let mut r = start;
    for _ in 0..50 {
        dec.charge_work(1)?;
        let (mut jtj, mut jtf) = (0.0, 0.0);
        for q in &local {
            dec.charge_work(1)?;
            let d = [
                q[0] - (corner[0] + inward[0] * r),
                q[1] - (corner[1] + inward[1] * r),
            ];
            let rho = d[0].hypot(d[1]).max(f64::MIN_POSITIVE);
            let f = rho - r;
            let j = -(d[0] * inward[0] + d[1] * inward[1]) / rho - 1.0;
            jtj += j * j;
            jtf += j * f;
        }
        if jtj.is_nan() || jtj <= 0.0 {
            break;
        }
        let step = jtf / jtj;
        let next = r - step;
        if next.is_nan() || next <= 0.0 || !next.is_finite() {
            break;
        }
        r = next;
        if step.abs() <= 1e-12 * r {
            break;
        }
    }
    Ok(r)
}

struct Frame {
    axes: [V3; 3],
    rotation: [Number; 4],
    /// The exact rotation matrix of `rotation`, by rows.
    matrix: [[Rational; 3]; 3],
    /// The stock's exact placed corner and its exact sizes.
    origin: [Rational; 3],
    size: [Rational; 3],
    lo: [f64; 3],
    hi: [f64; 3],
}

/// The exact rotation matrix, by rows, of a rational quaternion of any
/// nonzero length.
fn exact_matrix(q: &[Number; 4]) -> [[Rational; 3]; 3] {
    let [w, x, y, z] = q.clone().map(|n| n.value().clone());
    let m = |a: &Rational, b: &Rational| a * b;
    let n = &(&m(&w, &w) + &m(&x, &x)) + &(&m(&y, &y) + &m(&z, &z));
    let two = Rational::from(2i64);
    let rows = [
        [
            &(&m(&w, &w) + &m(&x, &x)) - &(&m(&y, &y) + &m(&z, &z)),
            &two * &(&m(&x, &y) - &m(&w, &z)),
            &two * &(&m(&x, &z) + &m(&w, &y)),
        ],
        [
            &two * &(&m(&x, &y) + &m(&w, &z)),
            &(&m(&w, &w) - &m(&x, &x)) + &(&m(&y, &y) - &m(&z, &z)),
            &two * &(&m(&y, &z) - &m(&w, &x)),
        ],
        [
            &two * &(&m(&x, &z) - &m(&w, &y)),
            &two * &(&m(&y, &z) + &m(&w, &x)),
            &(&m(&w, &w) - &m(&x, &x)) - &(&m(&y, &y) - &m(&z, &z)),
        ],
    ];
    rows.map(|row| row.map(|v| v.checked_div(&n).expect("nonzero quaternion")))
}

impl Frame {
    fn local(&self, p: V3) -> V3 {
        self.axes.map(|a| dot(a, p))
    }

    fn world(&self, l: V3) -> V3 {
        (0..3).fold([0.0; 3], |s, k| {
            [
                s[0] + self.axes[k][0] * l[k],
                s[1] + self.axes[k][1] * l[k],
                s[2] + self.axes[k][2] * l[k],
            ]
        })
    }

    /// The pose at an exact offset from the stock's corner in its own axes,
    /// so tools placed on a stock face lie exactly in that face.
    fn pose(&self, turn: [i64; 4], offset: &[Rational; 3]) -> Pose {
        let translation = std::array::from_fn(|r| {
            (0..3).fold(self.origin[r].clone(), |s, c| {
                &s + &(&self.matrix[r][c] * &offset[c])
            })
        });
        Pose {
            rotation: compose(&self.rotation, turn),
            translation: Vector3(translation.map(Number::new)),
        }
    }

    /// The exact offset of a level along axis `k` from the stock's corner.
    fn offset(&self, k: usize, level: f64) -> Rational {
        grid(level - self.lo[k])
    }
}

fn refuse(code: &'static str, detail: impl Into<String>) -> Refusal {
    Refusal::new(code, detail)
}

/// Groups plane normals into parallel classes, largest support first.
fn classes(
    surfaces: &[(Surface, Vec<u32>)],
    dec: &mut Decider,
) -> Result<Vec<(V3, usize, Vec<usize>)>, Refusal> {
    dec.charge_work(surfaces.len() as u64)?;
    let cos = ANGLE_DEGREES.to_radians().cos();
    let mut planes: Vec<usize> = (0..surfaces.len())
        .filter(|&s| matches!(surfaces[s].0, Surface::Plane { .. }))
        .collect();
    super::charge_sort(planes.len(), dec)?;
    planes.sort_by(|&a, &b| {
        surfaces[b]
            .1
            .len()
            .cmp(&surfaces[a].1.len())
            .then(a.cmp(&b))
    });
    dec.check_deadline()?;
    let mut out: Vec<(V3, usize, Vec<usize>)> = Vec::new();
    for s in planes {
        let Surface::Plane { normal, .. } = surfaces[s].0 else {
            unreachable!()
        };
        let mut matched = None;
        for (i, class) in out.iter().enumerate() {
            dec.charge_work(1)?;
            if dot(class.0, normal).abs() >= cos {
                matched = Some(i);
                break;
            }
        }
        match matched {
            Some(i) => {
                out[i].1 += surfaces[s].1.len();
                out[i].2.push(s);
            }
            None => out.push((normal, surfaces[s].1.len(), vec![s])),
        }
    }
    Ok(out)
}

/// Builds the recipe for `scan`'s fitted surfaces, or refuses when no block
/// stock can be found.
pub(super) fn assemble(
    scan: &Operation,
    surfaces: &[(Surface, Vec<u32>)],
    points: &[V3],
    sigma: f64,
    dec: &mut Decider,
) -> Result<Assembly, Refusal> {
    dec.charge_work(surfaces.len() as u64)?;
    let eps = super::segment::GROWTH_SIGMAS * sigma;
    let sin = ANGLE_DEGREES.to_radians().sin();
    let cos = ANGLE_DEGREES.to_radians().cos();
    let mut uses = vec![SurfaceUse::Unplaced; surfaces.len()];
    let mut unplaced = Vec::new();

    // Three mutually perpendicular plane classes with the most support.
    let classes = classes(surfaces, dec)?;
    let mut best: Option<([usize; 3], usize)> = None;
    for a in 0..classes.len() {
        dec.charge_work(1)?;
        for b in a + 1..classes.len() {
            dec.charge_work(1)?;
            for c in b + 1..classes.len() {
                dec.charge_work(1)?;
                let [x, y, z] = [&classes[a], &classes[b], &classes[c]];
                let perpendicular = dot(x.0, y.0).abs() <= sin
                    && dot(x.0, z.0).abs() <= sin
                    && dot(y.0, z.0).abs() <= sin;
                let support = x.1 + y.1 + z.1;
                if perpendicular && best.is_none_or(|(_, s)| support > s) {
                    best = Some(([a, b, c], support));
                }
            }
        }
    }
    let Some((chosen, _)) = best else {
        return Err(refuse(
            "reconstruction.no-stock",
            "no three mutually perpendicular plane directions were fitted",
        ));
    };
    // Local axes follow the world axes they are nearest, as a right-handed frame.
    let mut directions: Vec<V3> = chosen.iter().map(|&c| classes[c].0).collect();
    let mut axes = [[0.0; 3]; 3];
    for k in 0..2 {
        dec.charge_work(1)?;
        let (i, sign) = directions
            .iter()
            .enumerate()
            .map(|(i, d)| (i, d[k]))
            .max_by(|a, b| a.1.abs().total_cmp(&b.1.abs()).then(b.0.cmp(&a.0)))
            .expect("three directions");
        axes[k] = scale(directions.remove(i), sign.signum());
    }
    axes[1] = unit(sub(axes[1], scale(axes[0], dot(axes[0], axes[1])))).expect("perpendicular");
    axes[2] = cross(axes[0], axes[1]);
    let q = quaternion(axes).map(|x| (x * ROTATION_GRID as f64).round() / ROTATION_GRID as f64);
    let axes = columns(q);
    let rotation = q.map(quaternion_number);
    let class_axis = |c: usize| (0..3).find(|&k| dot(axes[k], classes[c].0).abs() >= cos);

    // Caps: planes across a cylinder or cone axis whose samples lie within
    // its radius, as a boss top or a bore bottom.
    let mut cap_of: Vec<Option<usize>> = vec![None; surfaces.len()];
    for (s, (surface, members)) in surfaces.iter().enumerate() {
        dec.charge_work(1)?;
        let Surface::Plane { normal, .. } = surface else {
            continue;
        };
        for (t, (other, _)) in surfaces.iter().enumerate() {
            dec.charge_work(1)?;
            let Some((point, axis)) = other.axis() else {
                continue;
            };
            if dot(*normal, axis).abs() < cos {
                continue;
            }
            let mut inside = true;
            for &i in members {
                dec.charge_work(1)?;
                let q = sub(points[i as usize], point);
                let along = dot(q, axis);
                let radial = norm(sub(q, scale(axis, along)));
                if !other.radius_at(along).is_some_and(|r| radial <= r + eps) {
                    inside = false;
                    break;
                }
            }
            if inside {
                cap_of[s] = Some(t);
                break;
            }
        }
    }

    // Stock faces: the outermost non-cap planes of each chosen class.
    let mut lo = [0.0; 3];
    let mut hi = [0.0; 3];
    for &c in &chosen {
        dec.charge_work(1)?;
        let k = class_axis(c).ok_or_else(|| {
            refuse(
                "reconstruction.no-stock",
                "fitted plane classes do not fit the orthogonalized stock frame",
            )
        })?;
        let mut offsets: Vec<(f64, usize)> = classes[c]
            .2
            .iter()
            .filter(|&&s| cap_of[s].is_none())
            .map(|&s| {
                let Surface::Plane { normal, offset } = surfaces[s].0 else {
                    unreachable!()
                };
                (offset * dot(normal, axes[k]).signum(), s)
            })
            .collect();
        super::charge_sort(offsets.len(), dec)?;
        offsets.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        dec.check_deadline()?;
        let (Some(&(low, a)), Some(&(high, b))) = (offsets.first(), offsets.last()) else {
            return Err(refuse(
                "reconstruction.no-stock",
                format!("stock axis {k} has no face plane"),
            ));
        };
        if high - low <= 2.0 * eps {
            return Err(refuse(
                "reconstruction.no-stock",
                format!("stock axis {k} has one face plane, not two"),
            ));
        }
        lo[k] = low;
        hi[k] = high;
        uses[a] = SurfaceUse::StockFace;
        uses[b] = SurfaceUse::StockFace;
    }
    let matrix = exact_matrix(&rotation);
    let corner = (0..3).fold([0.0; 3], |s, k| {
        [
            s[0] + axes[k][0] * lo[k],
            s[1] + axes[k][1] * lo[k],
            s[2] + axes[k][2] * lo[k],
        ]
    });
    let frame = Frame {
        axes,
        rotation,
        matrix,
        origin: corner.map(grid),
        size: std::array::from_fn(|k| grid(hi[k] - lo[k])),
        lo,
        hi,
    };
    let block = Operation::Block {
        size: frame.size.clone().map(Number::new),
        pose: frame.pose([1, 0, 0, 0], &[0, 0, 0].map(Rational::from)),
    };

    let mut steps = vec![
        Step {
            id: id("scan"),
            op: scan.clone(),
        },
        Step {
            id: id("stock"),
            op: block,
        },
    ];
    let mut current = id("stock");
    let mut bosses = Vec::new();
    let mut bores = Vec::new();
    let mut fillets: Vec<RoundedEdge> = Vec::new();
    let mut refined = Vec::new();
    for (s, (surface, members)) in surfaces.iter().enumerate() {
        dec.charge_work(1)?;
        let Some((point, axis)) = surface.axis() else {
            continue;
        };
        let Some(k) = (0..3).find(|&k| dot(frame.axes[k], axis).abs() >= cos) else {
            unplaced.push(Unplaced {
                surface: s,
                samples: members.len(),
                reason: "axis is not parallel to a stock axis",
            });
            continue;
        };
        let foot = frame.local(point);
        let along: Vec<f64> = members
            .iter()
            .map(|&i| {
                dec.charge_work(1)?;
                Ok(frame.local(points[i as usize])[k])
            })
            .collect::<Result<_, Refusal>>()?;
        dec.charge_work((along.len() * 2) as u64)?;
        let t_min = along.iter().copied().fold(f64::INFINITY, f64::min);
        let t_max = along.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let radius_at = |level: f64| {
            let t = (level - foot[k]) * dot(frame.axes[k], axis);
            surface.radius_at(t).expect("axial surface")
        };
        let others = [(k + 1) % 3, (k + 2) % 3];
        let r_mid = radius_at((t_min + t_max) / 2.0);
        let inside_face = others
            .iter()
            .all(|&o| foot[o] - r_mid > frame.lo[o] && foot[o] + r_mid < frame.hi[o]);
        // A convex fillet: the axis lies one radius inside two stock faces
        // and the samples face the edge between them.
        if let Surface::Cylinder { radius, .. } = *surface {
            let tol = (eps).max(0.1 * radius);
            let sides: Vec<Option<bool>> = others
                .iter()
                .map(|&o| {
                    if (foot[o] - (frame.lo[o] + radius)).abs() <= tol {
                        Some(false)
                    } else if (foot[o] - (frame.hi[o] - radius)).abs() <= tol {
                        Some(true)
                    } else {
                        None
                    }
                })
                .collect();
            if let [Some(a), Some(b)] = sides[..] {
                let outward =
                    members
                        .iter()
                        .try_fold(0usize, |count, &i| -> Result<usize, Refusal> {
                            dec.charge_work(1)?;
                            let l = frame.local(points[i as usize]);
                            let outside = others.iter().zip([a, b]).all(|(&o, high)| {
                                let d = l[o] - foot[o];
                                if high { d >= -eps } else { d <= eps }
                            });
                            Ok(count + usize::from(outside))
                        })?;
                if outward * 3 >= members.len() * 2 {
                    let edge = [(others[0], a), (others[1], b)];
                    let r = tangent_radius(&frame, points, members, edge, radius, dec)?;
                    let mut foot = [0.0; 3];
                    for (o, high) in edge {
                        dec.charge_work(1)?;
                        foot[o] = if high {
                            frame.hi[o] - r
                        } else {
                            frame.lo[o] + r
                        };
                    }
                    refined.push((
                        s,
                        Surface::Cylinder {
                            point: frame.world(foot),
                            axis: super::linalg::canonical(frame.axes[k]),
                            radius: r,
                        },
                    ));
                    uses[s] = SurfaceUse::Fillet;
                    fillets.push((r, s, edge));
                    continue;
                }
            }
        }
        let mut cap = None;
        for (c, candidate) in cap_of.iter().enumerate() {
            dec.charge_work(1)?;
            if *candidate == Some(s) {
                cap = Some(c);
                break;
            }
        }
        let cap_level = cap.map(|c| {
            let Surface::Plane { normal, offset } = surfaces[c].0 else {
                unreachable!()
            };
            offset * dot(normal, frame.axes[k]).signum()
        });
        let above = t_min >= frame.hi[k] - eps;
        let below = t_max <= frame.lo[k] + eps;
        let within = t_min >= frame.lo[k] - eps && t_max <= frame.hi[k] + eps;
        if !wraps(surface, members, points, dec)? {
            unplaced.push(Unplaced {
                surface: s,
                samples: members.len(),
                reason: "a boss or bore needs samples more than half way round its axis",
            });
            continue;
        }
        if (above || below) && inside_face {
            let Some(level) = cap_level else {
                unplaced.push(Unplaced {
                    surface: s,
                    samples: members.len(),
                    reason: "boss without a fitted cap plane",
                });
                continue;
            };
            // A boss stands exactly on its stock face.
            let height = grid(if above {
                level - frame.hi[k]
            } else {
                frame.lo[k] - level
            });
            let (base, top) = if above {
                (frame.size[k].clone(), &frame.size[k] + &height)
            } else {
                (-&height, Rational::zero())
            };
            if height.to_f64_nearest() <= 2.0 * eps {
                unplaced.push(Unplaced {
                    surface: s,
                    samples: members.len(),
                    reason: "boss cap does not lie beyond the stock face",
                });
                continue;
            }
            uses[s] = SurfaceUse::Boss;
            uses[cap.expect("cap level")] = SurfaceUse::Cap;
            bosses.push((s, k, foot, base, top));
        } else if within && inside_face {
            let overlap = grid(OVERLAP_MM.min((frame.hi[k] - frame.lo[k]) / 2.0));
            let (base, top) = match cap_level {
                // A blind bore opens on the face its samples reach.
                Some(level) if (t_max - frame.hi[k]).abs() < (t_min - frame.lo[k]).abs() => {
                    (frame.offset(k, level), &frame.size[k] + &overlap)
                }
                Some(level) => (-&overlap, frame.offset(k, level)),
                None => (-&overlap, &frame.size[k] + &overlap),
            };
            uses[s] = SurfaceUse::Bore;
            if let Some(c) = cap {
                uses[c] = SurfaceUse::Cap;
            }
            bores.push((s, k, foot, base, top));
        } else {
            unplaced.push(Unplaced {
                surface: s,
                samples: members.len(),
                reason: "neither a boss, a bore nor a stock-edge fillet",
            });
        }
    }
    let axial = |s: usize, k: usize, foot: V3, base: Rational, top: Rational| {
        let surface = &surfaces[s].0;
        let axis = surface.axis().expect("axial").1;
        let r = |offset: &Rational| {
            let level = frame.lo[k] + offset.to_f64_nearest();
            let t = (level - foot[k]) * dot(frame.axes[k], axis);
            surface.radius_at(t).expect("axial")
        };
        let mut at: [Rational; 3] = std::array::from_fn(|o| frame.offset(o, foot[o]));
        at[k] = base.clone();
        let pose = frame.pose(TURN_TO[k], &at);
        let height = Number::new(&top - &base);
        match surface {
            Surface::Cone { .. } => Operation::Cone {
                base_radius: length(r(&base)),
                top_radius: length(r(&top)),
                height,
                pose,
            },
            _ => Operation::Cylinder {
                radius: length(r(&base)),
                height,
                pose,
            },
        }
    };
    for (n, (s, k, foot, base, top)) in bosses.into_iter().enumerate() {
        dec.charge_work(1)?;
        let tool = id(&format!("boss-{}", n + 1));
        let joined = id(&format!("with-boss-{}", n + 1));
        steps.push(Step {
            id: tool.clone(),
            op: axial(s, k, foot, base, top),
        });
        steps.push(Step {
            id: joined.clone(),
            op: Operation::Boolean {
                kind: BooleanKind::Union,
                target: current,
                tool,
            },
        });
        current = joined;
    }
    for (n, (s, k, foot, base, top)) in bores.into_iter().enumerate() {
        dec.charge_work(1)?;
        let tool = id(&format!("bore-{}", n + 1));
        let cut = id(&format!("with-bore-{}", n + 1));
        steps.push(Step {
            id: tool.clone(),
            op: axial(s, k, foot, base, top),
        });
        steps.push(Step {
            id: cut.clone(),
            op: Operation::Boolean {
                kind: BooleanKind::Subtract,
                target: current,
                tool,
            },
        });
        current = cut;
    }
    // One radius per stock edge, weighted by samples; edges whose radii
    // agree within the noise share one step, so they round together.
    let mut edges: std::collections::BTreeMap<Edge, (f64, usize)> =
        std::collections::BTreeMap::new();
    for (radius, s, edge) in fillets {
        dec.charge_work(1)?;
        let mut key = edge;
        key.sort();
        let weight = surfaces[s].1.len();
        let entry = edges.entry(key).or_default();
        entry.0 += radius * weight as f64;
        entry.1 += weight;
    }
    let mut edges: Vec<RoundedEdge> = edges
        .into_iter()
        .map(|(key, (sum, weight))| (sum / weight as f64, weight, key))
        .collect();
    super::charge_sort(edges.len(), dec)?;
    edges.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.2.cmp(&b.2)));
    let mut groups: Vec<Vec<RoundedEdge>> = Vec::new();
    for f in edges {
        dec.charge_work(1)?;
        match groups.last_mut() {
            Some(g) if (f.0 - g[0].0).abs() <= eps.max(0.02 * g[0].0) => g.push(f),
            _ => groups.push(vec![f]),
        }
    }
    for (n, group) in groups.iter().enumerate() {
        dec.charge_work(1)?;
        let weight: usize = group.iter().map(|f| f.1).sum();
        let radius = group.iter().map(|f| f.0 * f.1 as f64).sum::<f64>() / weight as f64;
        let edges = group
            .iter()
            .map(|f| EdgeRef {
                between: f.2.map(|(axis, high)| OriginPattern {
                    step: id("stock"),
                    role: Some(side_role(axis, high)),
                }),
                touching: Vec::new(),
                near: None,
            })
            .collect();
        let step = id(&format!("fillet-{}", n + 1));
        steps.push(Step {
            id: step.clone(),
            op: Operation::Fillet {
                body: current,
                edges: EdgeSelection::Refs(edges),
                radius: length(radius),
                overflow: None,
                corner: None,
            },
        });
        current = step;
    }
    for (s, cap) in cap_of.iter().enumerate() {
        dec.charge_work(1)?;
        if cap.is_some() && uses[s] == SurfaceUse::Unplaced {
            unplaced.push(Unplaced {
                surface: s,
                samples: surfaces[s].1.len(),
                reason: "cap of a surface that was not placed",
            });
        }
    }
    for (s, (surface, members)) in surfaces.iter().enumerate() {
        dec.charge_work(1)?;
        if matches!(surface, Surface::Plane { .. })
            && uses[s] == SurfaceUse::Unplaced
            && cap_of[s].is_none()
        {
            unplaced.push(Unplaced {
                surface: s,
                samples: members.len(),
                reason: "plane is neither an outer stock face nor a cap",
            });
        }
    }
    super::charge_sort(unplaced.len(), dec)?;
    unplaced.sort_by_key(|u| u.surface);
    dec.check_deadline()?;
    let mut recipe = Recipe::new("scan reconstruction", steps, current);
    recipe.format = catallax_recipe::FORMAT_V4.into();
    Ok(Assembly {
        recipe,
        uses,
        unplaced,
        refined,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn antipodal_samples_cover_exactly_half_a_circle() {
        let surface = Surface::Cylinder {
            point: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
            radius: 1.0,
        };
        let points: Vec<_> = (0..60)
            .map(|i| [if i < 30 { 1.0 } else { -1.0 }, 0.0, i as f64])
            .collect();
        assert!(
            !wraps(
                &surface,
                &(0..60).collect::<Vec<_>>(),
                &points,
                &mut Decider::default()
            )
            .unwrap()
        );
    }

    #[test]
    fn a_density_biased_120_degree_cylinder_is_unplaced() {
        let size = [30.0, 20.0, 8.0];
        let mut points = Vec::new();
        let mut surfaces = Vec::new();
        for axis in 0..3 {
            for high in [false, true] {
                let mut members = Vec::new();
                for u in [0.0, 1.0] {
                    for v in [0.0, 1.0] {
                        let mut p = [0.0; 3];
                        p[axis] = if high { size[axis] } else { 0.0 };
                        p[(axis + 1) % 3] = u * size[(axis + 1) % 3];
                        p[(axis + 2) % 3] = v * size[(axis + 2) % 3];
                        members.push(points.len() as u32);
                        points.push(p);
                    }
                }
                let mut normal = [0.0; 3];
                normal[axis] = 1.0;
                surfaces.push((
                    Surface::Plane {
                        normal,
                        offset: if high { size[axis] } else { 0.0 },
                    },
                    members,
                ));
            }
        }
        let mut members = Vec::new();
        for i in 0..120 {
            let angle = if i < 90 { 0.0 } else { 120.0_f64.to_radians() };
            let angle = angle + (i % 30) as f64 * 1e-6;
            members.push(points.len() as u32);
            points.push([
                15.0 + 4.0 * angle.cos(),
                10.0 + 4.0 * angle.sin(),
                8.0 * i as f64 / 120.0,
            ]);
        }
        surfaces.push((
            Surface::Cylinder {
                point: [15.0, 10.0, 0.0],
                axis: [0.0, 0.0, 1.0],
                radius: 4.0,
            },
            members,
        ));
        let scan = Operation::ScanImport {
            file: catallax_recipe::FileRef {
                name: "partial-cylinder.obj".into(),
                sha256: "0".repeat(64),
            },
            pose: Pose::identity(),
            scale: length(1.0),
        };
        let assembly = assemble(&scan, &surfaces, &points, 0.01, &mut Decider::default())
            .expect("stock planes fit");
        assert_eq!(assembly.uses[6], SurfaceUse::Unplaced);
        assert!(
            assembly
                .unplaced
                .iter()
                .any(|u| { u.surface == 6 && u.reason.contains("more than half way round") })
        );
    }
}
