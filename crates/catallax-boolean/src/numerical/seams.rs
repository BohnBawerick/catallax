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
                    xyz: [&point.x, &point.y, &point.z]
                        .map(|coordinate| interval(coordinate, dec).map(|b| b.mid()))
                        .into_iter()
                        .collect::<Result<Vec<_>, Refusal>>()?
                        .try_into()
                        .map_err(|_| unreachable!())?,
                })
            })
            .collect()
    }
}
