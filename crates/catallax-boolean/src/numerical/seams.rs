//! Ordered generator roots own cone seam meetings across cylinder quarters.

use super::*;
use catallax_numbers::RealPoly;
use std::collections::BTreeMap;

#[derive(Clone)]
pub(super) struct Owner {
    pub uv: [Rational; 2],
    pub xyz: [Rational; 3],
}

pub(super) struct Roots<'a> {
    source: &'a BSplineSurface,
    cylinder: &'a Cylinder,
    cache: BTreeMap<Rational, Vec<Owner>>,
}

impl<'a> Roots<'a> {
    pub fn new(intersection: &'a NumericalNativeIntersection) -> Result<Self, Refusal> {
        let Surface::Cylinder(cylinder) = &intersection.windows().originals()[1] else {
            return Err(bad("finite seam ownership lost its cylinder"));
        };
        Ok(Self {
            source: &intersection.windows().sources()[0],
            cylinder,
            cache: BTreeMap::new(),
        })
    }

    pub fn owner(
        &mut self,
        u: Rational,
        xyz: &[Rational; 3],
        tolerance: &Rational,
        dec: &mut Decider,
    ) -> Result<(usize, Owner), Refusal> {
        if !self.cache.contains_key(&u) {
            let roots = self.build(&u, dec)?;
            self.cache.insert(u.clone(), roots);
        }
        let roots = &self.cache[&u];
        let mut distances = roots
            .iter()
            .enumerate()
            .map(|(i, root)| {
                dec.charge_work(1)?;
                Ok((
                    distance(&events::box_point(xyz), &events::box_point(&root.xyz))?,
                    i,
                ))
            })
            .collect::<Result<Vec<_>, Refusal>>()?;
        distances.sort();
        let Some((gap, index)) = distances.first() else {
            return Err(bad("finite cone seam has no isolated cylinder root"));
        };
        if gap > tolerance
            || distances
                .get(1)
                .is_some_and(|(next, _)| gap * &q(2) >= *next)
        {
            return Err(bad(
                "finite cone seam needs separated root ownership refinement",
            ));
        }
        let owner = roots[*index].clone();
        if owner.uv[1] < q(0) || owner.uv[1] > q(1) {
            return Err(bad("owned cone seam root is outside the finite band"));
        }
        Ok((*index, owner))
    }

    fn build(&self, u: &Rational, dec: &mut Decider) -> Result<Vec<Owner>, Refusal> {
        let endpoints = [q(0), q(1)].map(|v| {
            self.source
                .evaluate_box([u, u], [&v, &v], BITS, dec)
                .map(|point| point.map(|coordinate| Real::from(coordinate.mid())))
        });
        let [a, b] = endpoints;
        let [a, b] = [a?, b?].map(|p| Point3::new(p[0].clone(), p[1].clone(), p[2].clone()));
        let direction = b.minus(&a).cross(&self.cylinder.axis);
        let offset = a.minus(&self.cylinder.origin).cross(&self.cylinder.axis);
        let polynomial = RealPoly::from_terms(
            1,
            vec![
                (vec![2], direction.norm_squared()),
                (vec![1], Real::from(2) * offset.dot(&direction)),
                (
                    vec![0],
                    offset.norm_squared()
                        - &self.cylinder.radius_sq * &self.cylinder.axis.norm_squared(),
                ),
            ],
        );
        polynomial
            .real_roots(dec)?
            .into_iter()
            .map(|root| {
                let v = interval(&root, dec)?.mid();
                let point = a.plus(&b.minus(&a).scaled(&Real::from(v.clone())));
                Ok(Owner {
                    uv: [u.clone(), v],
                    xyz: [
                        interval(&point.x, dec)?.mid(),
                        interval(&point.y, dec)?.mid(),
                        interval(&point.z, dec)?.mid(),
                    ],
                })
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn close_generator_roots_keep_distinct_owners_inside_the_physical_tolerance() {
        let mut dec = Decider::default();
        let carrier = Cylinder {
            origin: Point3::new(0.into(), 11.into(), 0.into()),
            axis: Vector3::new(0.into(), 1.into(), 0.into()),
            radius_sq: 1.into(),
            ref_dir: Vector3::new(1.into(), 0.into(), 0.into()),
        };
        let tolerance = Rational::parse("1e-7").unwrap();
        let source = NumericalCylinderSource::new(
            &carrier,
            [q(0), q(1)],
            NumericalOptions {
                tolerance_mm: tolerance.clone(),
            },
            &mut dec,
        )
        .unwrap();
        let tool = Cylinder {
            origin: Point3::new(
                0.into(),
                Rational::from_ratio(23, 2).unwrap().into(),
                0.into(),
            ),
            axis: Vector3::new(1.into(), 0.into(), 0.into()),
            radius_sq: Rational::parse("1e-18").unwrap().into(),
            ref_dir: Vector3::new(0.into(), 1.into(), 0.into()),
        };
        let mut roots = Roots {
            source: source.source(),
            cylinder: &tool,
            cache: BTreeMap::new(),
        };
        let delta = Rational::parse("1e-9").unwrap();
        let middle = Rational::from_ratio(23, 2).unwrap();
        let lower = roots
            .owner(q(0), &[q(1), &middle - &delta, q(0)], &tolerance, &mut dec)
            .unwrap();
        let upper = roots
            .owner(q(0), &[q(1), &middle + &delta, q(0)], &tolerance, &mut dec)
            .unwrap();
        assert_eq!((lower.0, upper.0), (0, 1));
        assert!(lower.1.uv[1] < upper.1.uv[1]);
        assert!(
            roots
                .owner(q(0), &[q(1), middle, q(0)], &tolerance, &mut dec)
                .is_err()
        );
    }
}
