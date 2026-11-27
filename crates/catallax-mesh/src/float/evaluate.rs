use std::f64::consts::TAU;

use catallax_shapes::{ParamKind, Point3, Surface, Vector3};

use super::{DisplayFrame, DisplaySurface, Uv, Xyz, spline};

pub(super) fn finite(x: f64) -> Option<f64> {
    x.is_finite().then_some(if x == 0.0 { 0.0 } else { x })
}

pub(super) fn positive(x: f64) -> Option<f64> {
    (x > 0.0 && x.is_finite()).then_some(x)
}

pub(super) fn scalar(x: &catallax_numbers::Real) -> f64 {
    x.display_value().unwrap_or(f64::NAN)
}

pub(super) fn point(p: &Point3) -> Option<Xyz> {
    let [x, y, z] = [&p.x, &p.y, &p.z].map(|v| finite(scalar(v)));
    Some([x?, y?, z?])
}

fn vector(p: &Vector3) -> Xyz {
    [&p.x, &p.y, &p.z].map(scalar)
}

pub(super) fn sub(a: Xyz, b: Xyz) -> Xyz {
    std::array::from_fn(|i| a[i] - b[i])
}

pub(super) fn add(a: Xyz, b: Xyz) -> Xyz {
    std::array::from_fn(|i| a[i] + b[i])
}

pub(super) fn scale(a: Xyz, s: f64) -> Xyz {
    a.map(|x| x * s)
}

pub(super) fn dot(a: Xyz, b: Xyz) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(super) fn cross(a: Xyz, b: Xyz) -> Xyz {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}

pub(super) fn norm(v: Xyz) -> f64 {
    libm::sqrt(dot(v, v))
}

pub(super) fn unit(v: Xyz) -> Option<Xyz> {
    let largest = positive(v[0].abs().max(v[1].abs()).max(v[2].abs()))?;
    let v = v.map(|x| x / largest);
    Some(scale(v, 1.0 / positive(norm(v))?))
}

pub(super) fn angle(y: f64, x: f64) -> f64 {
    libm::atan2(y, x) / TAU
}

pub(super) fn align(x: f64, near: f64, period: f64) -> f64 {
    x + libm::round((near - x) / period) * period
}

pub(super) fn solve_gram(a: Xyz, b: Xyz, rhs: Uv) -> Option<Uv> {
    let aa = dot(a, a);
    let bb = dot(b, b);
    let ab = dot(a, b);
    let det = aa * bb - ab * ab;
    (det > aa * bb * 1e-14 && det.is_finite()).then(|| {
        [
            (bb * rhs[0] - ab * rhs[1]) / det,
            (aa * rhs[1] - ab * rhs[0]) / det,
        ]
    })
}

pub(super) fn frame(
    p: &Point3,
    axis: &Vector3,
    reference: &Vector3,
    origin: Xyz,
) -> Option<DisplayFrame> {
    let z = unit(vector(axis))?;
    let y = unit(cross(z, unit(vector(reference))?))?;
    Some(DisplayFrame {
        origin: sub(point(p)?, origin),
        x: cross(y, z),
        y,
        z,
    })
}

pub(super) fn anchor(s: &Surface) -> Option<Xyz> {
    point(match s {
        Surface::Plane(s) => &s.origin,
        Surface::Cylinder(s) => &s.origin,
        Surface::Cone(s) => &s.apex,
        Surface::Sphere(s) => &s.center,
        Surface::Torus(s) => &s.center,
        Surface::TorusSheet(s) => &s.center,
        Surface::EllipticCylinder(s) => &s.origin,
        Surface::BSpline(s) => s.carrier().controls().first()?,
        _ => return None,
    })
}

pub(super) fn describe(s: &Surface, origin: Xyz) -> Option<DisplaySurface> {
    use DisplaySurface as D;
    let f = |p, a, r| frame(p, a, r, origin);
    let sqrt = |x| positive(libm::sqrt(x));
    Some(match s {
        Surface::Plane(s) => D::Plane(f(&s.origin, &s.normal, &s.ref_dir)?),
        Surface::Cylinder(s) => D::Cylinder {
            frame: f(&s.origin, &s.axis, &s.ref_dir)?,
            radius: sqrt(scalar(&s.radius_sq))?,
        },
        Surface::Cone(s) => D::Cone {
            frame: f(&s.apex, &s.axis, &s.ref_dir)?,
            slope: positive(scalar(&s.slope))?,
        },
        Surface::Sphere(s) => D::Sphere {
            frame: f(&s.center, &s.axis, &s.ref_dir)?,
            radius: positive(scalar(&s.radius))?,
        },
        Surface::Torus(s) => D::Torus {
            frame: f(&s.center, &s.axis, &s.ref_dir)?,
            major: positive(scalar(&s.major))?,
            minor: positive(scalar(&s.minor))?,
        },
        Surface::TorusSheet(s) => D::Torus {
            frame: f(&s.center, &s.axis, &s.ref_dir)?,
            major: positive(scalar(&s.major))?,
            minor: positive(scalar(&s.minor))?,
        },
        Surface::EllipticCylinder(s) => D::EllipticCylinder {
            frame: f(&s.origin, &s.axis, &s.major_dir)?,
            radii: [sqrt(scalar(&s.major_sq))?, sqrt(scalar(&s.minor_sq))?],
        },
        Surface::BSpline(s) => D::BSpline(spline::net(s.carrier(), origin)?),
        _ => return None,
    })
}

pub(super) struct Projector<'a> {
    surface: &'a DisplaySurface,
    spline: Option<spline::Patch<'a>>,
    pub(super) domain: [Option<[f64; 2]>; 2],
    pub(super) periods: [Option<f64>; 2],
}

impl<'a> Projector<'a> {
    pub(super) fn new(surface: &'a DisplaySurface, source: Option<&Surface>) -> Option<Self> {
        use DisplaySurface as D;
        let mut spline = None;
        let (domain, periods) = match surface {
            D::Plane(_) => ([None, None], [None, None]),
            D::Cylinder { .. } | D::Cone { .. } | D::EllipticCylinder { .. } => {
                ([Some([0.0, 1.0]), None], [Some(1.0), None])
            }
            D::Sphere { .. } => ([Some([0.0, 1.0]), Some([-0.25, 0.25])], [Some(1.0), None]),
            D::Torus { major, minor, .. } if matches!(source, Some(Surface::TorusSheet(_))) => {
                let end = libm::acos(-major / minor) / TAU;
                ([Some([0.0, 1.0]), Some([-end, end])], [Some(1.0), None])
            }
            D::Torus { .. } => ([Some([0.0, 1.0]); 2], [Some(1.0); 2]),
            D::BSpline(net) => {
                let patch = spline::Patch::new(net)?;
                let domain = patch.bounds;
                let kinds = match source {
                    Some(Surface::BSpline(s)) => s.param_kinds(),
                    _ => [ParamKind::Finite; 2],
                };
                let periods = std::array::from_fn(|i| {
                    (kinds[i] == ParamKind::PeriodicFinite || patch.closed(i))
                        .then(|| domain[i][1] - domain[i][0])
                });
                spline = Some(patch);
                (domain.map(Some), periods)
            }
            D::Grid(_) => return None,
        };
        Some(Self {
            surface,
            spline,
            domain,
            periods,
        })
    }

    pub(super) fn canonical(&self, mut uv: Uv) -> Uv {
        for (i, coord) in uv.iter_mut().enumerate() {
            if let Some(period) = self.periods[i] {
                let base = self.domain[i].unwrap()[0];
                // Keep the principal upper endpoint to evaluate a clamped spline.
                if *coord < base || *coord > base + period {
                    *coord -= libm::floor((*coord - base) / period) * period;
                }
            }
        }
        uv
    }

    pub(super) fn project(&mut self, p: Xyz, previous: Option<Uv>) -> Uv {
        use DisplaySurface as D;
        let previous = previous.map(|uv| self.canonical(uv));
        if let Some(patch) = &mut self.spline {
            return patch.project(p, previous, self.periods);
        }
        let frame = match self.surface {
            D::Plane(frame)
            | D::Cylinder { frame, .. }
            | D::Cone { frame, .. }
            | D::Sphere { frame, .. }
            | D::Torus { frame, .. }
            | D::EllipticCylinder { frame, .. } => frame,
            _ => unreachable!("only direct carriers have a projector"),
        };
        let delta = sub(p, frame.origin);
        let [x, y, z] = [frame.x, frame.y, frame.z].map(|axis| dot(delta, axis));
        let u = if x == 0.0 && y == 0.0 {
            previous.map_or(0.0, |p| p[0])
        } else {
            angle(y, x)
        };
        match self.surface {
            D::Plane(_) => [x, y],
            D::Cylinder { .. } | D::Cone { .. } => [u, z],
            D::Sphere { .. } => [u, angle(z, libm::sqrt(x * x + y * y))],
            D::Torus { major, .. } => [u, angle(z, libm::sqrt(x * x + y * y) - major)],
            D::EllipticCylinder { radii, .. } => [angle(y / radii[1], x / radii[0]), z],
            _ => unreachable!(),
        }
    }

    pub(super) fn position(&mut self, uv: Uv) -> Xyz {
        if let Some(patch) = &mut self.spline {
            patch.at(uv)
        } else {
            self.jet(uv).0
        }
    }

    // Closest-parameter surface point and its two parametric tangents. The
    // analytic projection is a chart inverse on the carrier; off it, only the
    // tangents and normal component of the displacement are used.
    pub(super) fn jet(&mut self, uv: Uv) -> (Xyz, [Xyz; 2]) {
        use DisplaySurface as D;
        if let Some(patch) = &mut self.spline {
            return patch.jet(uv);
        }
        let [u, v] = uv;
        if let D::Plane(f) = self.surface {
            return (add(f.origin, add(scale(f.x, u), scale(f.y, v))), [f.x, f.y]);
        }
        let (s, c) = libm::sincos(u * TAU);
        let circle = |f: &DisplayFrame, rx: f64, ry: f64| {
            (
                add(scale(f.x, rx * c), scale(f.y, ry * s)),
                add(scale(f.x, -rx * s), scale(f.y, ry * c)),
            )
        };
        match self.surface {
            D::Cylinder { frame: f, radius } => {
                let (r, t) = circle(f, *radius, *radius);
                (add(f.origin, add(r, scale(f.z, v))), [t, f.z])
            }
            D::EllipticCylinder { frame: f, radii } => {
                let (r, t) = circle(f, radii[0], radii[1]);
                (add(f.origin, add(r, scale(f.z, v))), [t, f.z])
            }
            D::Cone { frame: f, slope } => {
                let (r, t) = circle(f, *slope, *slope);
                let dv = add(f.z, r);
                (add(f.origin, scale(dv, v)), [t, dv])
            }
            D::Sphere { frame: f, radius } => {
                let (r, t) = circle(f, *radius, *radius);
                let (sv, cv) = libm::sincos(v * TAU);
                (
                    add(f.origin, add(scale(r, cv), scale(f.z, radius * sv))),
                    [t, add(scale(r, -sv), scale(f.z, radius * cv))],
                )
            }
            D::Torus {
                frame: f,
                major,
                minor,
            } => {
                let (r, t) = circle(f, 1.0, 1.0);
                let (sv, cv) = libm::sincos(v * TAU);
                (
                    add(
                        f.origin,
                        add(scale(r, major + minor * cv), scale(f.z, minor * sv)),
                    ),
                    [t, add(scale(r, -sv), scale(f.z, cv))],
                )
            }
            _ => unreachable!(),
        }
    }
}
