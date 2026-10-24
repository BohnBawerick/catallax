//! Exact meridian/parallel trim rectangles, including seam-crossing long arcs.

use catallax_intersect::CurveForm;
use catallax_numbers::{Ball, Decider, Outcome, Rational, Real, Refusal, StrictSign};
use catallax_shapes::{Point3, Torus, Vector3};
use catallax_topology::{FaceId, Loop, Model, Sense};

use crate::integrate::{Stop, Work};

pub(crate) struct Domain {
    coordinates: [Real; 8],
    u: Sweep,
    v: Sweep,
    pub(crate) orientation: i64,
}

struct Sweep {
    cosine: Real,
    sine: Real,
    direction: i64,
    wrap: bool,
}

fn sign(value: &Real, dec: &mut Decider) -> Result<Option<i64>, Refusal> {
    Ok(match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Some(1),
        Outcome::Proven(StrictSign::Negative) => Some(-1),
        Outcome::ProvenCoincident(_) => Some(0),
        Outcome::NotExact(_) => None,
        Outcome::Refused(error) => return Err(error),
    })
}

fn zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(value, dec)? == Some(0))
}

impl Sweep {
    fn new(
        from: &[Real],
        to: &[Real],
        direction: i64,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        let cosine = &from[0] * &to[0] + &from[1] * &to[1];
        let sine = (&from[0] * &to[1] - &from[1] * &to[0]) * Real::from_i64(direction);
        let Some(c) = sign(&cosine, dec)? else {
            return Ok(None);
        };
        let Some(s) = sign(&sine, dec)? else {
            return Ok(None);
        };
        if c > 0 && s == 0 {
            return Ok(None);
        }
        Ok(Some(Self {
            cosine,
            sine,
            direction,
            wrap: c >= 0 && s < 0,
        }))
    }

    fn evaluate(&self, work: &mut Work) -> Result<Ball, Stop> {
        let cosine = work.real(&self.cosine)?;
        let sine = work.real(&self.sine)?;
        let mut angle = crate::integrate::angle(&cosine, &sine, work)?;
        if self.wrap {
            let pi = work.pi()?;
            let two = work.rational(&Rational::from_i64(2))?;
            let period = work.mul(&pi, &two)?;
            angle = work.add(&angle, &period)?;
        }
        let direction = work.rational(&Rational::from_i64(self.direction))?;
        work.mul(&angle, &direction)
    }
}

impl Domain {
    pub(crate) fn extras(&self, work: &mut Work) -> Result<[Ball; 10], Stop> {
        Ok([
            work.real(&self.coordinates[0])?,
            work.real(&self.coordinates[1])?,
            work.real(&self.coordinates[2])?,
            work.real(&self.coordinates[3])?,
            work.real(&self.coordinates[4])?,
            work.real(&self.coordinates[5])?,
            work.real(&self.coordinates[6])?,
            work.real(&self.coordinates[7])?,
            self.u.evaluate(work)?,
            self.v.evaluate(work)?,
        ])
    }
}

struct Frame<'a> {
    torus: &'a Torus,
    axis: Vector3,
    first: Vector3,
    second: Vector3,
}

impl<'a> Frame<'a> {
    fn new(torus: &'a Torus) -> Self {
        let unit = |v: &Vector3| v.scaled(&(&Real::one() / &v.norm_squared().sqrt()));
        let axis = unit(&torus.axis);
        let first = unit(&torus.ref_dir);
        let second = axis.cross(&first);
        Self {
            torus,
            axis,
            first,
            second,
        }
    }

    fn coordinates(&self, point: &Point3, dec: &mut Decider) -> Result<Option<[Real; 4]>, Refusal> {
        let q = point.minus(&self.torus.center);
        let x = q.dot(&self.first);
        let y = q.dot(&self.second);
        let rho = (x.square() + y.square()).sqrt();
        if sign(&rho, dec)? != Some(1) {
            return Ok(None);
        }
        let cv = &(&rho - &self.torus.major) / &self.torus.minor;
        let sv = &q.dot(&self.axis) / &self.torus.minor;
        if !zero(&(cv.square() + sv.square() - Real::one()), dec)? {
            return Ok(None);
        }
        Ok(Some([&x / &rho, &y / &rho, cv, sv]))
    }

    fn direction(
        &self,
        circle: &CurveForm,
        at: &[Real; 4],
        dec: &mut Decider,
    ) -> Result<Option<(usize, i64)>, Refusal> {
        let CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ..
        } = circle
        else {
            return Ok(None);
        };
        let offset = center.minus(&self.torus.center);
        if zero(&normal.cross(&self.axis).norm_squared(), dec)? {
            let expected = self.axis.scaled(&(&self.torus.minor * &at[3]));
            let radial = &self.torus.major + &(&self.torus.minor * &at[2]);
            if !zero(&offset.minus(&expected).norm_squared(), dec)?
                || !zero(&(radius_sq - &radial.square()), dec)?
            {
                return Ok(None);
            }
            let Some(direction) = sign(&normal.dot(&self.axis), dec)? else {
                return Ok(None);
            };
            return Ok((direction != 0).then_some((0, direction)));
        }
        let expected = self
            .first
            .scaled(&at[0])
            .plus(&self.second.scaled(&at[1]))
            .scaled(&self.torus.major);
        let meridian_normal = offset.cross(&self.axis);
        if !zero(&offset.minus(&expected).norm_squared(), dec)?
            || !zero(&(radius_sq - &self.torus.minor.square()), dec)?
            || !zero(&normal.cross(&meridian_normal).norm_squared(), dec)?
        {
            return Ok(None);
        }
        let Some(direction) = sign(&normal.dot(&meridian_normal), dec)? else {
            return Ok(None);
        };
        Ok((direction != 0).then_some((1, direction)))
    }
}

pub(crate) fn prepare(
    model: &Model,
    face: FaceId,
    torus: &Torus,
    dec: &mut Decider,
) -> Result<Option<Domain>, Refusal> {
    let [lp] = model.face(face).loops.as_slice() else {
        return Ok(None);
    };
    let Loop::Coedges(coedges) = model.loop_(*lp) else {
        return Ok(None);
    };
    let [a, b, c, d] = coedges.as_slice() else {
        return Ok(None);
    };
    let coedges = [a, b, c, d].map(|id| model.coedge(*id));
    let edges = coedges.map(|c| model.edge(c.edge));
    let starts = std::array::from_fn::<_, 4, _>(|i| {
        if coedges[i].sense == Sense::Same {
            edges[i].start
        } else {
            edges[i].end
        }
    });
    let [Some(a), Some(b), Some(c), Some(d)] = starts else {
        return Ok(None);
    };
    for i in 0..4 {
        let end = if coedges[i].sense == Sense::Same {
            edges[i].end
        } else {
            edges[i].start
        };
        if end != starts[(i + 1) % 4] {
            return Ok(None);
        }
    }
    let frame = Frame::new(torus);
    let Some(a) = frame.coordinates(&model.geom().position(model.vertex(a).point), dec)? else {
        return Ok(None);
    };
    let Some(b) = frame.coordinates(&model.geom().position(model.vertex(b).point), dec)? else {
        return Ok(None);
    };
    let Some(c) = frame.coordinates(&model.geom().position(model.vertex(c).point), dec)? else {
        return Ok(None);
    };
    let Some(d) = frame.coordinates(&model.geom().position(model.vertex(d).point), dec)? else {
        return Ok(None);
    };
    let coordinates = [a, b, c, d];
    let mut directions = [(0, 0); 4];
    for i in 0..4 {
        let circle = catallax_intersect::curve_form(model.geom(), edges[i].curve, dec)?;
        let Some((axis, direction)) = frame.direction(&circle, &coordinates[i], dec)? else {
            return Ok(None);
        };
        let sense = if (edges[i].sense == Sense::Reversed) != (coedges[i].sense == Sense::Reversed)
        {
            -1
        } else {
            1
        };
        directions[i] = (axis, direction * sense);
        let fixed = if axis == 0 { 2 } else { 0 };
        for (from, to) in coordinates[i][fixed..fixed + 2]
            .iter()
            .zip(&coordinates[(i + 1) % 4][fixed..fixed + 2])
        {
            if !zero(&(from - to), dec)? {
                return Ok(None);
            }
        }
    }
    if directions[0].0 == directions[1].0
        || directions[1].0 == directions[2].0
        || directions[2].0 == directions[3].0
        || directions[3].0 == directions[0].0
        || directions[0].1 != -directions[2].1
        || directions[1].1 != -directions[3].1
    {
        return Ok(None);
    }
    let u_edge = if directions[0].0 == 0 { 0 } else { 1 };
    let v_edge = 1 - u_edge;
    let Some(u) = Sweep::new(
        &coordinates[u_edge][..2],
        &coordinates[u_edge + 1][..2],
        directions[u_edge].1,
        dec,
    )?
    else {
        return Ok(None);
    };
    let Some(v) = Sweep::new(
        &coordinates[v_edge][2..],
        &coordinates[v_edge + 1][2..],
        directions[v_edge].1,
        dec,
    )?
    else {
        return Ok(None);
    };
    let [[a, b, c, d], [e, f, g, h], [i, j, k, l], _] = coordinates;
    let coordinates = if u_edge == 0 {
        [a, b, e, f, g, h, k, l]
    } else {
        [e, f, i, j, c, d, g, h]
    };
    Ok(Some(Domain {
        coordinates,
        u,
        v,
        orientation: if u_edge == 0 { 1 } else { -1 },
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Accuracy, Enclosure, mass_properties};
    use catallax_shapes::{Pose, RVector3, Rotation};

    fn b(n: i64, d: i64) -> Ball {
        Ball::from_rational(&Rational::from_ratio(n, d).unwrap(), 512)
    }

    fn contains(actual: &Enclosure, expected: &Ball) {
        assert!(
            actual.lo <= expected.lo() && actual.hi >= expected.hi(),
            "{actual:?} excludes {expected:?}"
        );
        assert!(Accuracy::default().reached(actual));
    }

    #[test]
    fn rectangular_torus_sectors_preserve_long_arcs_and_posed_inertia() {
        for long in [false, true] {
            let pose = if long {
                Pose::from_rational(
                    Rotation::from_quaternion(1.into(), 2.into(), 3.into(), 4.into()).unwrap(),
                    RVector3::new(7.into(), (-5).into(), 11.into()),
                )
            } else {
                Pose::identity()
            };
            let (model, body) = crate::fixtures::torus_rectangle(long, &pose);
            let result =
                mass_properties(&model, body, &Accuracy::default(), &mut Decider::default())
                    .unwrap();
            assert!(result.marks.is_empty(), "{:?}", result.marks);
            let mass = result.value;
            let pi = Ball::pi(512);
            let length = &pi * &b(if long { 3 } else { 1 }, 2);
            let sign = b(if long { -1 } else { 1 }, 1);
            // Integrate cylindrical horizontal sections of the solid, rather
            // than a torus flux or a boundary-form antiderivative.
            let volume = &length * &(b(29, 6) + &pi * &b(3, 4));
            let area = &length * &(&pi * &b(3, 2) + b(27, 2)) + b(6, 1) + &pi * &b(1, 2);
            let first_xy = b(11, 1) + &pi * &b(37, 16);
            let first = [&first_xy * &sign, first_xy, &length * &b(27, 8)];
            let center = first.each_ref().map(|v| v.checked_div(&volume).unwrap());
            let radial = b(1763, 60) + &pi * &b(117, 16);
            let xx = &(&length * &radial) * &b(1, 2);
            let xy = &radial * &b(1, 2);
            let xz = &b(499, 60) * &sign;
            let yz = b(499, 60);
            let zz = &length * &(b(47, 30) + &pi * &b(3, 16));
            let raw = [
                [xx.clone(), xy.clone(), xz.clone()],
                [xy, xx, yz.clone()],
                [xz, yz, zz],
            ];
            let centered = std::array::from_fn::<_, 3, _>(|i| {
                std::array::from_fn::<_, 3, _>(|j| {
                    &raw[i][j] - &(&first[i] * &first[j]).checked_div(&volume).unwrap()
                })
            });
            let tensor = std::array::from_fn::<_, 3, _>(|i| {
                std::array::from_fn::<_, 3, _>(|j| {
                    if i == j {
                        &centered[(i + 1) % 3][(i + 1) % 3] + &centered[(i + 2) % 3][(i + 2) % 3]
                    } else {
                        -&centered[i][j]
                    }
                })
            });
            contains(&mass.volume, &volume);
            contains(&mass.area, &area);
            let matrix = pose.rotation.rational_matrix().unwrap();
            let translation = [
                &pose.translation.x,
                &pose.translation.y,
                &pose.translation.z,
            ];
            let actual_center = mass.centroid.unwrap();
            let actual_tensor = mass.inertia.unwrap();
            for i in 0..3 {
                let expected = (0..3).fold(
                    Ball::from_rational(translation[i].rational_ref().unwrap(), 512),
                    |sum, k| sum + &Ball::from_rational(matrix[i][k], 512) * &center[k],
                );
                contains(&actual_center[i], &expected);
                for j in 0..3 {
                    let expected = (0..3).fold(b(0, 1), |sum, k| {
                        sum + (0..3).fold(b(0, 1), |sum, l| {
                            sum + &tensor[k][l]
                                * &Ball::from_rational(&(matrix[i][k] * matrix[j][l]), 512)
                        })
                    });
                    contains(&actual_tensor[i][j], &expected);
                }
            }
        }
    }

    #[test]
    fn torus_rectangle_reversal_keeps_the_complement_and_material_sense() {
        use catallax_topology::{Coedge, Face};
        let (mut model, body) = crate::fixtures::torus_rectangle(false, &Pose::identity());
        let face = model
            .faces_of_body(body)
            .into_iter()
            .find(|face| {
                matches!(
                    model.geom().surface(model.face(*face).surface),
                    catallax_shapes::Surface::Torus(_)
                )
            })
            .unwrap();
        let original = model.face(face).clone();
        let origins = model.face_origins(face).to_vec();
        let Loop::Coedges(coedges) = model.loop_(original.loops[0]) else {
            unreachable!()
        };
        let reversed = coedges
            .iter()
            .rev()
            .map(|id| {
                let c = model.coedge(*id);
                Coedge {
                    edge: c.edge,
                    sense: if c.sense == Sense::Same {
                        Sense::Reversed
                    } else {
                        Sense::Same
                    },
                }
            })
            .collect::<Vec<_>>();
        let reversed = reversed.into_iter().map(|c| model.add_coedge(c)).collect();
        let lp = model.add_loop(Loop::Coedges(reversed));
        let pi = Ball::pi(512);
        let patch = &(&pi * &pi) * &b(3, 4) + &pi * &b(1, 2);
        let full = &(&pi * &pi) * &b(12, 1);
        for sense in [Sense::Same, Sense::Reversed] {
            let face = model.add_named_face(
                Face {
                    surface: original.surface,
                    sense,
                    loops: vec![lp],
                },
                origins.clone(),
            );
            let answer =
                crate::face_area(&model, face, &Accuracy::default(), &mut Decider::default())
                    .unwrap();
            assert!(answer.marks.is_empty(), "{:?}", answer.marks);
            let expected = if sense == Sense::Same {
                &full - &patch
            } else {
                patch.clone()
            };
            contains(&answer.value, &expected);
        }
    }
}
