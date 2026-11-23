//! Boundary exclusion for imported rectangular patch shells.

use super::*;
use crate::import_spline::bilinear;
use catallax_shapes::{SplineSupport, Surface};

type Affine<T> = [T; 4];

fn evaluate(form: &Affine<Interval>, at: &[Interval; 3]) -> Interval {
    form[0] * at[0] + form[1] * at[1] + form[2] * at[2] + form[3]
}

fn affine(axis: V, origin: &V) -> Affine<Rational> {
    let offset = -dot(&axis, origin);
    let [x, y, z] = axis;
    [x, y, z, offset]
}

fn strict(range: Interval) -> bool {
    range.lo() > 0.0 || range.hi() < 0.0
}

struct Coordinates {
    exact: [Affine<Rational>; 3],
    intervals: [Affine<Interval>; 3],
}

impl Coordinates {
    fn new(patch: &PatchFace) -> Option<Self> {
        let surface = Surface::BSpline(SplineSupport::RationalPatch(patch.support.clone()));
        let (a, b, c, d) = bilinear(&surface).ok()??.basis();
        let a = point(&a);
        let scale = b.cross(&c).dot(&d);
        if scale.is_zero() {
            return None;
        }
        let inverse = Rational::one() / scale;
        // In the basis (b, c, d), the patch is (u, v, u*v) on the unit square.
        let exact = [c.cross(&d), d.cross(&b), b.cross(&c)].map(|v| {
            let v = v.scaled(&inverse);
            affine([v.x, v.y, v.z], &a)
        });
        let intervals = exact
            .each_ref()
            .map(|row| row.each_ref().map(Interval::from_rational));
        Some(Self { exact, intervals })
    }

    fn misses(&self, at: &[Interval; 3]) -> bool {
        let [u, v, uv] = self.intervals.each_ref().map(|row| evaluate(row, at));
        [u, v].iter().any(|s| s.lo() > 1.0 || s.hi() < 0.0) || strict(uv - u * v)
    }

    fn contains(&self, at: &V) -> bool {
        let [u, v, uv] = self
            .exact
            .each_ref()
            .map(|[x, y, z, offset]| &(&(x * &at[0]) + &(y * &at[1])) + &(&(z * &at[2]) + offset));
        [&u, &v]
            .iter()
            .all(|&s| s.sign() != Sign::Negative && *s <= Rational::one())
            && uv == u * v
    }
}

enum FaceBound {
    Polygon {
        plane: Affine<Interval>,
        bounds: [Interval; 3],
    },
    Bilinear(Box<Coordinates>),
    Patch(Vec<[Interval; 3]>),
}

pub(crate) struct Boundary {
    inverse: [[Interval; 3]; 3],
    translation: [Interval; 3],
    faces: Vec<FaceBound>,
}

impl Boundary {
    pub(crate) fn new(shell: &Shell, dec: &mut Decider) -> Result<Option<Self>, Refusal> {
        if !shell.rectangular_import {
            return Ok(None);
        }
        let mut faces = Vec::with_capacity(shell.faces.len());
        for face in &shell.faces {
            dec.charge_work(1)?;
            faces.push(match face {
                Face::Polygon(p) => FaceBound::Polygon {
                    plane: affine(p.normal.clone(), &p.vertices[0])
                        .each_ref()
                        .map(Interval::from_rational),
                    bounds: p.rings().flatten().fold(
                        p.vertices[0].each_ref().map(Interval::from_rational),
                        |bounds, point| {
                            std::array::from_fn(|i| {
                                bounds[i].hull(&Interval::from_rational(&point[i]))
                            })
                        },
                    ),
                },
                Face::Patch(p) => match Coordinates::new(p) {
                    Some(uv) => FaceBound::Bilinear(Box::new(uv)),
                    None => FaceBound::Patch(p.distance.boxes()),
                },
            });
        }
        Ok(Some(Self {
            inverse: shell
                .inverse_rotation
                .each_ref()
                .map(|row| row.each_ref().map(Interval::from_rational)),
            translation: shell.translation.each_ref().map(Interval::from_rational),
            faces,
        }))
    }

    pub(crate) fn disjoint(
        &self,
        world: &[Interval; 3],
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let delta: [Interval; 3] = std::array::from_fn(|i| world[i] - self.translation[i]);
        let local = self
            .inverse
            .map(|[a, b, c]| a * delta[0] + b * delta[1] + c * delta[2]);
        for face in &self.faces {
            dec.charge_work(1)?;
            let misses = match face {
                FaceBound::Polygon { plane, bounds } => {
                    disjoint(bounds, &local) || strict(evaluate(plane, &local))
                }
                FaceBound::Bilinear(uv) => uv.misses(&local),
                FaceBound::Patch(boxes) => boxes.iter().all(|b| disjoint(b, &local)),
            };
            if !misses {
                return Ok(false);
            }
        }
        Ok(true)
    }

    pub(crate) fn location(
        &self,
        shell: &Shell,
        point: &V,
        boundary_excluded: bool,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        let local = shell.local(point);
        if (0..3)
            .any(|i| local[i] < shell.local_bounds[i][0] || local[i] > shell.local_bounds[i][1])
        {
            return Ok(Location::Outside);
        }
        if !boundary_excluded {
            for (face, bound) in shell.faces.iter().zip(&self.faces) {
                dec.charge_work(1)?;
                match (face, bound) {
                    (Face::Polygon(p), FaceBound::Polygon { .. })
                        if dot(&p.normal, &sub(&local, &p.vertices[0])).is_zero()
                            && p.planar(&local) != Planar::Outside =>
                    {
                        return Ok(Location::Boundary);
                    }
                    (_, FaceBound::Bilinear(uv)) if uv.contains(&local) => {
                        return Ok(Location::Boundary);
                    }
                    (_, FaceBound::Patch(_)) => return shell.location(point, dec),
                    _ => (),
                }
            }
        }
        for direction in DIRECTIONS {
            let direction = direction.map(|(n, m)| Rational::from_ratio(n, m).expect("nonzero"));
            if let Some(winding) = shell.winding(&local, &direction, dec)? {
                return Ok(if winding == 0 {
                    Location::Outside
                } else {
                    Location::Inside
                });
            }
        }
        Ok(Location::Indeterminate)
    }
}

fn disjoint(a: &[Interval; 3], b: &[Interval; 3]) -> bool {
    (0..3).any(|i| a[i].hi() < b[i].lo() || b[i].hi() < a[i].lo())
}
