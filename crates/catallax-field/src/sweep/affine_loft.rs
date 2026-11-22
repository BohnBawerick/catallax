//! Two sections related by an affine map in space.

use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_recipe::{LoftSection, StepId};
use catallax_shapes::{BSplineSurface, RPoint3, RVector3, RationalPatch};
use std::collections::BTreeMap;

use crate::Plan;
use crate::geometry::Leaf;
use crate::patch_shell::{Face, PatchFace, Polygon, Shell};
use crate::profile::{P, Profile, cross, rational_point, sub};

use super::stations::knots;

fn unit(axis: usize) -> RVector3 {
    let [x, y, z] = std::array::from_fn(|i| {
        if i == axis {
            Rational::one()
        } else {
            Rational::zero()
        }
    });
    RVector3::new(x, y, z)
}

fn coordinates(v: RVector3) -> [Rational; 3] {
    [v.x, v.y, v.z]
}

/// Each corner's symmetric Jacobian is positive definite by Gershgorin.
/// The Jacobian is affine in (x,y,t), so convexity of that matrix cone
/// proves positive definiteness over the box. Integrating along any
/// segment in the box gives strict monotonicity and hence injectivity.
fn embedded(
    source: &[P],
    a: &RVector3,
    b: &RVector3,
    c: &RVector3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    if c.z.is_zero() {
        return Ok(false);
    }
    let limits = |axis: usize| {
        [
            source.iter().map(|q| &q[axis]).min().expect("polygon"),
            source.iter().map(|q| &q[axis]).max().expect("polygon"),
        ]
    };
    let da = a.minus(&unit(0));
    let db = b.minus(&unit(1));
    for t in [Rational::zero(), Rational::one()] {
        let u = coordinates(unit(0).plus(&da.scaled(&t)));
        let v = coordinates(unit(1).plus(&db.scaled(&t)));
        for x in limits(0) {
            for y in limits(1) {
                dec.charge_work(1)?;
                // The third source coordinate is c.z * t.
                let w = coordinates(
                    c.plus(&da.scaled(x))
                        .plus(&db.scaled(y))
                        .scaled(&(&Rational::one() / &c.z)),
                );
                let columns = [&u, &v, &w];
                for (i, diagonal) in columns.iter().enumerate() {
                    let radius = (0..3)
                        .filter(|&j| j != i)
                        .map(|j| (&columns[j][i] + &columns[i][j]).abs().mul_pow2(-1))
                        .fold(Rational::zero(), |sum, bound| sum + bound);
                    if diagonal[i] <= radius {
                        return Ok(false);
                    }
                }
            }
        }
    }
    Ok(true)
}

impl Plan {
    pub(super) fn affine_loft(
        &mut self,
        sections: &[LoftSection],
        profiles: &BTreeMap<StepId, Profile>,
        dec: &mut Decider,
    ) -> Result<Option<usize>, Refusal> {
        let [first, last] = sections else {
            return Ok(None);
        };
        let ordered = |section: &LoftSection| -> Option<Vec<P>> {
            let [edges] = profiles[&section.profile].loops.as_slice() else {
                return None;
            };
            if edges.len() < 3 || edges.iter().any(|e| e.circle.is_some() || e.ids.len() != 1) {
                return None;
            }
            let start = edges.iter().position(|e| e.ids[0] == section.start[0])?;
            Some(
                (0..edges.len())
                    .map(|i| edges[(start + i) % edges.len()].start.clone())
                    .collect(),
            )
        };
        let (Some(bottom), Some(top)) = (ordered(first), ordered(last)) else {
            return Ok(None);
        };
        if top.len() != bottom.len() {
            return Ok(None);
        }
        let lower = &profiles[&first.profile];
        let upper = &profiles[&last.profile];
        let inverse = lower.pose.inverse();
        let relative = upper.pose.then(&inverse);
        if crate::profile::rational_vector(&relative, &unit(2))
            .cross(&unit(2))
            .is_zero()
        {
            return Ok(None);
        }
        let origin = RPoint3::new(bottom[0][0].clone(), bottom[0][1].clone(), Rational::zero());
        let source: Vec<_> = bottom.iter().map(|p| sub(p, &bottom[0])).collect();
        let target: Vec<_> = top
            .iter()
            .map(|p| {
                rational_point(
                    &relative,
                    &RPoint3::new(p[0].clone(), p[1].clone(), Rational::zero()),
                )
                .minus(&origin)
            })
            .collect();
        let Some((j, determinant)) = (2..source.len()).find_map(|j| {
            let area = cross(&source[1], &source[j]);
            (!area.is_zero()).then_some((j, area))
        }) else {
            return Ok(None);
        };
        let c = &target[0];
        let d = target[1].minus(c);
        let e = target[j].minus(c);
        let inverse_det = &Rational::one() / &determinant;
        let a = d
            .scaled(&source[j][1])
            .minus(&e.scaled(&source[1][1]))
            .scaled(&inverse_det);
        let b = e
            .scaled(&source[1][0])
            .minus(&d.scaled(&source[j][0]))
            .scaled(&inverse_det);
        for (q, p) in source.iter().zip(&target) {
            dec.charge_work(1)?;
            if c.plus(&a.scaled(&q[0])).plus(&b.scaled(&q[1])) != *p {
                return Ok(None);
            }
        }
        if !embedded(&source, &a, &b, c, dec)? {
            return Ok(None);
        }
        let area = bottom
            .iter()
            .zip(bottom.iter().cycle().skip(1))
            .fold(Rational::zero(), |sum, (a, b)| sum + cross(a, b));
        let outward = area.sign() == c.z.sign();
        let mut bottom: Vec<_> = bottom.iter().map(|p| lower.world(p)).collect();
        let mut top: Vec<_> = top.iter().map(|p| upper.world(p)).collect();
        let mut faces = Vec::new();
        for i in 0..bottom.len() {
            let next = (i + 1) % bottom.len();
            let controls = [&bottom[i], &top[i], &bottom[next], &top[next]].map(RPoint3::to_real);
            let carrier = BSplineSurface::new(
                [1, 1],
                [knots(2), knots(2)],
                [2, 2],
                controls.into(),
                vec![Real::one(); 4],
                dec,
            )?;
            let patch = RationalPatch::new(carrier, dec)?;
            faces.push(Face::Patch(Box::new(PatchFace::new(&patch, outward, dec)?)));
        }
        if outward {
            bottom.reverse();
        } else {
            top.reverse();
        }
        for cap in [bottom, top] {
            faces.push(Face::Polygon(Box::new(Polygon::new(cap)?)));
        }
        self.leaf(Leaf::Shell(Box::new(Shell::new(faces)?)), vec![], dec)
            .map(Some)
    }
}
