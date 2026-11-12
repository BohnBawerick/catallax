//! Seeded noisy point clouds sampled from analytic part surfaces, for scan
//! fitting controls. Independent of the fitter: each part is written out
//! patch by patch from its analytic definition, and every sample keeps the
//! name of the patch it came from as ground truth for scoring only.
//!
//! The generator uses only addition, multiplication, division, square roots
//! and the standard sine and cosine, from a SplitMix64 stream. Noise is the
//! Irwin-Hall sum of twelve uniforms, so it needs no logarithm.
#![allow(dead_code)]

use std::f64::consts::{FRAC_PI_2, TAU};

type V3 = [f64; 3];
type Sample = Box<dyn Fn(f64, f64) -> Option<(V3, V3)>>;

/// One analytic patch: `sample` maps `(u, v)` uniform on the unit square to a
/// point and its outward unit normal, or `None` outside the patch, so that
/// accepted points are uniform by area over the patch.
pub struct Patch {
    pub name: &'static str,
    /// Area of the parameter domain before rejection, in square millimetres.
    pub area: f64,
    pub sample: Sample,
}

pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform on `[0, 1)` with 53 random bits.
    pub fn uniform(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Approximately standard normal: the sum of twelve uniforms, less six.
    pub fn normal(&mut self) -> f64 {
        (0..12).map(|_| self.uniform()).sum::<f64>() - 6.0
    }
}

/// `count` samples, each moved along its normal by noise of deviation
/// `sigma`, with the index of its patch.
pub fn sample(patches: &[Patch], count: usize, sigma: f64, seed: u64) -> Vec<(V3, usize)> {
    let total: f64 = patches.iter().map(|p| p.area).sum();
    let mut rng = Rng::new(seed);
    let mut out = Vec::with_capacity(count);
    while out.len() < count {
        let mut pick = rng.uniform() * total;
        let mut k = 0;
        while k + 1 < patches.len() && pick >= patches[k].area {
            pick -= patches[k].area;
            k += 1;
        }
        let (u, v) = (rng.uniform(), rng.uniform());
        if let Some((p, n)) = (patches[k].sample)(u, v) {
            let e = sigma * rng.normal();
            out.push((std::array::from_fn(|i| p[i] + e * n[i]), k));
        }
    }
    out
}

/// Wavefront OBJ of the samples as isolated vertices: a point cloud.
pub fn obj(samples: &[(V3, usize)]) -> String {
    let mut s = String::from("# seeded noisy point cloud; millimetres\n");
    for (p, _) in samples {
        s.push_str(&format!("v {:.6} {:.6} {:.6}\n", p[0], p[1], p[2]));
    }
    s
}

/// A rectangle on `z = level` (or `x`, `y`) with an optional circular cut.
fn plane(
    name: &'static str,
    axis: usize,
    level: f64,
    outward: f64,
    lo: [f64; 2],
    hi: [f64; 2],
    cut: Option<([f64; 2], f64)>,
) -> Patch {
    let (a, b) = ((axis + 1) % 3, (axis + 2) % 3);
    Patch {
        name,
        area: (hi[0] - lo[0]) * (hi[1] - lo[1]),
        sample: Box::new(move |u, v| {
            let s = lo[0] + u * (hi[0] - lo[0]);
            let t = lo[1] + v * (hi[1] - lo[1]);
            if let Some((c, r)) = cut
                && (s - c[0]).hypot(t - c[1]) < r
            {
                return None;
            }
            let mut p = [0.0; 3];
            p[axis] = level;
            p[a] = s;
            p[b] = t;
            let mut n = [0.0; 3];
            n[axis] = outward;
            Some((p, n))
        }),
    }
}

/// A disk on `z = level` facing `outward` along z.
fn disk(name: &'static str, centre: [f64; 2], radius: f64, level: f64, outward: f64) -> Patch {
    Patch {
        name,
        area: 4.0 * radius * radius,
        sample: Box::new(move |u, v| {
            let (x, y) = ((2.0 * u - 1.0) * radius, (2.0 * v - 1.0) * radius);
            (x.hypot(y) <= radius)
                .then(|| ([centre[0] + x, centre[1] + y, level], [0.0, 0.0, outward]))
        }),
    }
}

/// A z-axis cylinder or cone wall between two levels; `outward` is +1 for a
/// boss (material inside) and -1 for a bore (material outside).
fn wall(
    name: &'static str,
    centre: [f64; 2],
    radii: [f64; 2],
    levels: [f64; 2],
    outward: f64,
) -> Patch {
    let height = levels[1] - levels[0];
    let rmax = radii[0].max(radii[1]);
    let slant = (height * height + (radii[1] - radii[0]).powi(2)).sqrt();
    Patch {
        name,
        area: TAU * rmax * slant,
        sample: Box::new(move |u, v| {
            let r = radii[0] + v * (radii[1] - radii[0]);
            // Accept in proportion to the local circumference.
            if u * rmax > r {
                return None;
            }
            let theta = TAU * (u * rmax / r);
            let (c, s) = (theta.cos(), theta.sin());
            let z = levels[0] + v * height;
            let dr = radii[1] - radii[0];
            let nz = -dr / slant;
            let nr = height / slant;
            Some((
                [centre[0] + r * c, centre[1] + r * s, z],
                [outward * nr * c, outward * nr * s, outward * nz],
            ))
        }),
    }
}

/// The quarter-cylinder round of radius `r` on a block's top edge, mitred
/// against the rounds of the two neighbouring top edges. `axis` 0 runs along
/// x on the side `y = side`, axis 1 along y on the side `x = side`.
fn top_round(
    name: &'static str,
    axis: usize,
    side: f64,
    outward: f64,
    span: [f64; 2],
    top: f64,
    r: f64,
) -> Patch {
    Patch {
        name,
        area: (span[1] - span[0]) * FRAC_PI_2 * r,
        sample: Box::new(move |u, v| {
            let phi = v * FRAC_PI_2;
            let inset = r - r * phi.sin();
            let s = span[0] + u * (span[1] - span[0]);
            if s < span[0] + inset || s > span[1] - inset {
                return None;
            }
            let across = side - outward * inset;
            let z = top - r + r * phi.cos();
            let (sin, cos) = (phi.sin(), phi.cos());
            Some(if axis == 0 {
                ([s, across, z], [0.0, outward * sin, cos])
            } else {
                ([across, s, z], [outward * sin, 0.0, cos])
            })
        }),
    }
}

/// The concave round of radius `r` where a z-axis boss of radius `boss` meets
/// a plane `z = level`: a quarter of a torus.
fn junction(name: &'static str, centre: [f64; 2], boss: f64, level: f64, r: f64) -> Patch {
    let major = boss + r;
    Patch {
        name,
        area: TAU * major * FRAC_PI_2 * r,
        sample: Box::new(move |u, v| {
            let alpha = v * FRAC_PI_2;
            let rho = major - r * alpha.cos();
            if u * major > rho {
                return None;
            }
            let theta = TAU * (u * major / rho);
            let (c, s) = (theta.cos(), theta.sin());
            let z = level + r - r * alpha.sin();
            Some((
                [centre[0] + rho * c, centre[1] + rho * s, z],
                [alpha.cos() * c, alpha.cos() * s, alpha.sin()],
            ))
        }),
    }
}

/// A spherical dome of base radius `a` and height `a / 2` standing on the
/// plane `z = level`: a surface no plane, cylinder or cone explains.
pub fn dome(name: &'static str, centre: [f64; 2], a: f64, level: f64) -> Patch {
    let h = a / 2.0;
    let big = (a * a + h * h) / (2.0 * h);
    let floor = big / (big - h);
    Patch {
        name,
        area: 4.0 * a * a * floor,
        sample: Box::new(move |u, v| {
            let (x, y) = ((2.0 * u - 1.0) * a, (2.0 * v - 1.0) * a);
            let rho2 = x * x + y * y;
            if rho2 > a * a {
                return None;
            }
            let root = (big * big - rho2).sqrt();
            // Accept in proportion to the area element over the disk.
            let weight = (big / root) / floor;
            let w = (u * 7919.0 + v * 104_729.0).fract();
            if w > weight {
                return None;
            }
            let z = root - (big - h);
            Some((
                [centre[0] + x, centre[1] + y, level + z],
                [x / big, y / big, root / big],
            ))
        }),
    }
}

/// A 30 x 20 x 8 mm plate with a through bore at (15, 10): radius 4, or a
/// cone from radius 3 at the bottom to 5 at the top.
pub fn plate_with_bore(cone: bool) -> Vec<Patch> {
    let c = [15.0, 10.0];
    let radii = if cone { [3.0, 5.0] } else { [4.0, 4.0] };
    vec![
        plane(
            "bottom",
            2,
            0.0,
            -1.0,
            [0.0, 0.0],
            [30.0, 20.0],
            Some((c, radii[0])),
        ),
        plane(
            "top",
            2,
            8.0,
            1.0,
            [0.0, 0.0],
            [30.0, 20.0],
            Some((c, radii[1])),
        ),
        plane("x_min", 0, 0.0, -1.0, [0.0, 0.0], [20.0, 8.0], None),
        plane("x_max", 0, 30.0, 1.0, [0.0, 0.0], [20.0, 8.0], None),
        plane("y_min", 1, 0.0, -1.0, [0.0, 0.0], [8.0, 30.0], None),
        plane("y_max", 1, 20.0, 1.0, [0.0, 0.0], [8.0, 30.0], None),
        wall("bore", c, radii, [0.0, 8.0], -1.0),
    ]
}

/// P3's block, boss and top-edge rounds: block 40 x 30 x 10, boss radius 8
/// and height 12 at (20, 15), the four top edges rounded at radius 1. With
/// `junction_round` the boss meets the block through P3's radius-2 concave
/// round; without it the junction is sharp.
pub fn boss_block(junction_round: bool) -> Vec<Patch> {
    let c = [20.0, 15.0];
    let (cut, wall_from) = if junction_round {
        (10.0, 12.0)
    } else {
        (8.0, 10.0)
    };
    let mut patches = vec![
        plane("bottom", 2, 0.0, -1.0, [0.0, 0.0], [40.0, 30.0], None),
        plane(
            "top",
            2,
            10.0,
            1.0,
            [1.0, 1.0],
            [39.0, 29.0],
            Some((c, cut)),
        ),
        plane("x_min", 0, 0.0, -1.0, [0.0, 0.0], [30.0, 9.0], None),
        plane("x_max", 0, 40.0, 1.0, [0.0, 0.0], [30.0, 9.0], None),
        plane("y_min", 1, 0.0, -1.0, [0.0, 0.0], [9.0, 40.0], None),
        plane("y_max", 1, 30.0, 1.0, [0.0, 0.0], [9.0, 40.0], None),
        top_round("round_y_min", 0, 0.0, -1.0, [0.0, 40.0], 10.0, 1.0),
        top_round("round_y_max", 0, 30.0, 1.0, [0.0, 40.0], 10.0, 1.0),
        top_round("round_x_min", 1, 0.0, -1.0, [0.0, 30.0], 10.0, 1.0),
        top_round("round_x_max", 1, 40.0, 1.0, [0.0, 30.0], 10.0, 1.0),
        wall("boss", c, [8.0, 8.0], [wall_from, 22.0], 1.0),
        disk("boss_top", c, 8.0, 22.0, 1.0),
    ];
    if junction_round {
        patches.push(junction("junction", c, 8.0, 10.0, 2.0));
    }
    patches
}

/// Total area of the patches after rejection, by deterministic quadrature on
/// an `n` by `n` grid; for reports only.
pub fn area(patches: &[Patch], n: usize) -> Vec<f64> {
    patches
        .iter()
        .map(|p| {
            let mut hits = 0usize;
            for i in 0..n {
                for j in 0..n {
                    let (u, v) = ((i as f64 + 0.5) / n as f64, (j as f64 + 0.5) / n as f64);
                    hits += usize::from((p.sample)(u, v).is_some());
                }
            }
            p.area * hits as f64 / (n * n) as f64
        })
        .collect()
}
