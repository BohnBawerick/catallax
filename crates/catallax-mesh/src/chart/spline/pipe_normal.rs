use catallax_numbers::{Decider, Interval, Outcome, Real, Refusal, StrictSign};
use catallax_shapes::{RollingPipe, Vector3};

use super::{failure, point};

/// A positive multiple of the exact radial normal of a regular rolling pipe.
#[derive(Clone, Copy)]
pub(super) struct PipeNormal {
    plane: [Interval; 3],
    wall: [[Interval; 3]; 3],
    contact_dot: [Interval; 2],
}

impl PipeNormal {
    pub(super) fn new(pipe: &RollingPipe, dec: &mut Decider) -> Result<Self, Refusal> {
        let law = pipe.law();
        let wall_scale = &law.wall / &law.scale;
        let first = law.normal.scaled(&law.reach);
        let across = law.across.scaled(&wall_scale);
        let side = law.side.scaled(&wall_scale);
        let axial = law.axis.scaled(&-(&law.slope * &wall_scale));
        // At u = v = 0 the centre moves along side, and the arc moves
        // toward last. The triple product gives the chart's orientation,
        // including after a reflection. Regularity keeps that sign constant.
        let last = across.plus(&axial);
        let orientation = match dec.sign(&law.side.cross(&last).dot(&first)) {
            Outcome::Proven(StrictSign::Positive) => Real::one(),
            Outcome::Proven(StrictSign::Negative) => -&Real::one(),
            Outcome::ProvenCoincident(_) => {
                return Err(failure("rolling pipe has no chart orientation"));
            }
            Outcome::NotExact(reason) => return Err(reason.into()),
            Outcome::Refused(reason) => return Err(reason),
        };
        let enclose = |value: &Real| {
            dec.enclose_at(value, crate::chart::BITS)
                .map(|ball| ball.to_interval())
        };
        let vector = |v: Vector3| {
            let [x, y, z] = [&v.x, &v.y, &v.z].map(enclose);
            Ok::<_, Refusal>([x?, y?, z?])
        };
        let signed_wall = &wall_scale * &law.reach;
        Ok(Self {
            plane: vector(first.scaled(&orientation))?,
            wall: [
                vector(across.scaled(&orientation))?,
                vector(side.scaled(&orientation))?,
                vector(axial.scaled(&orientation))?,
            ],
            contact_dot: [
                enclose(&(&signed_wall * &law.tilt))?,
                enclose(&-(&(&signed_wall * &law.slope) * &law.along))?,
            ],
        })
    }

    pub(super) fn bounds(
        &self,
        quarter: usize,
        [u, v]: [Interval; 2],
        dec: &mut Decider,
    ) -> Result<[Interval; 3], Refusal> {
        dec.charge_work(42)?;
        // On a quarter turn, cos falls and sin rises as the half-angle
        // tangent runs from zero to one.
        let turn = |t: f64| {
            let t = point(t);
            let square = t.square();
            let denominator = point(1.0) + square;
            [
                (point(1.0) - square) / denominator,
                point(2.0) * t / denominator,
            ]
        };
        let [a, b] = [turn(u.lo()), turn(u.hi())];
        let c = a[0].hull(&b[0]);
        let s = a[1].hull(&b[1]);
        let [c, s] = match quarter {
            0 => [c, s],
            1 => [-s, c],
            2 => [-c, -s],
            _ => [s, -c],
        };
        let dot = self.contact_dot[0] * c + self.contact_dot[1];
        // If D = (2-v) F + v L, the contact ray is
        // 2 (F·D) D / |D|² - F. Expand its numerator before bounding it.
        let plane = point(4.0) * (point(1.0) - v);
        let wall = point(2.0) * v * (point(2.0) - v * (point(1.0) - dot));
        Ok(std::array::from_fn(|i| {
            plane * self.plane[i]
                + wall * (c * self.wall[0][i] + s * self.wall[1][i] + self.wall[2][i])
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{Ball, Rational};
    use catallax_shapes::{Cone, Cylinder, Plane, Point3};

    #[test]
    fn radial_normal_agrees_with_the_exact_carrier_on_both_pipe_laws() {
        let vector = |xyz: [i64; 3]| {
            let [x, y, z] = xyz.map(Real::from_i64);
            Vector3::new(x, y, z)
        };
        let plane = Plane {
            origin: Point3::new(Real::zero(), Real::zero(), Real::from_i64(3)),
            normal: vector([3, 0, 4]),
            ref_dir: vector([0, 1, 0]),
        };
        let cylinder = Cylinder {
            origin: Point3::origin(),
            axis: vector([0, 0, 1]),
            radius_sq: Real::from_i64(36),
            ref_dir: vector([1, 0, 0]),
        };
        let cone = Cone {
            apex: Point3::new(Real::zero(), Real::zero(), Real::from_i64(8)),
            axis: vector([0, 0, -1]),
            slope: Real::from(Rational::from_ratio(3, 4).unwrap()),
            ref_dir: vector([1, 0, 0]),
        };
        let mut dec = Decider::default();
        let radius = Real::from(Rational::from_ratio(1, 2).unwrap());
        let pipes = [
            RollingPipe::new(&plane, &cylinder, [true; 2], radius.clone(), true, &mut dec).unwrap(),
            RollingPipe::conical(&plane, &cone, [true; 2], radius, true, &mut dec).unwrap(),
        ];
        for pipe in pipes {
            let field = PipeNormal::new(&pipe, &mut dec).unwrap();
            for quarter in 0..4 {
                for v in [0.0, 0.25, 0.75, 1.0] {
                    let t = 0.375;
                    let radial = field
                        .bounds(quarter, [point(t), point(v)], &mut dec)
                        .unwrap();
                    let [u, v] = [quarter as f64 + t, v]
                        .map(|x| Ball::from_rational(&crate::chart::q(x), 128));
                    let [du, dv] = pipe
                        .carrier()
                        .derivatives(&u, &v, &dec)
                        .unwrap()
                        .map(|v| v.map(|x| x.to_interval()));
                    let normal: [_; 3] = std::array::from_fn(|i| {
                        let a = (i + 1) % 3;
                        let b = (i + 2) % 3;
                        du[a] * dv[b] - du[b] * dv[a]
                    });
                    let dot = (0..3).fold(point(0.0), |sum, i| sum + radial[i] * normal[i]);
                    assert!(dot.lo() > 0.0);
                    for i in 0..3 {
                        let a = (i + 1) % 3;
                        let b = (i + 2) % 3;
                        assert!((radial[a] * normal[b] - radial[b] * normal[a]).contains(0.0));
                    }
                }
            }
        }
    }
}
