//! Straight normal sweeps of one planar cap of a polygonal extrusion, independent of B-rep topology.

use super::{FastProfile, Feature, Sweep, bad, face};
use crate::{
    Face, Location, geometry,
    profile::{self, Edge},
    support::Support,
};
use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_recipe::{EdgeName, ElementId, Origin, Role, StepId};
use catallax_shapes::{Isometry, Plane, Point3, RVector3, Surface, Vector3};

/// A planar cap swept along its own normal, with the roles that evaluate its faces.
pub(crate) struct CapSweep {
    pub feature: Feature,
    pub faces: Vec<Face>,
    pub roles: Vec<Role>,
}

fn unit(sign: i64) -> Rational {
    Rational::from_i64(sign)
}

impl Feature {
    /// Whether the profile is one polygon with optional holes.
    fn simple_polygon(&self, dec: &mut Decider) -> Result<bool, Refusal> {
        if self
            .loops
            .iter()
            .flatten()
            .any(|e| e.circle.is_some() || e.full)
        {
            return Ok(false);
        }
        let mut outer = 0;
        for (i, edges) in self.loops.iter().enumerate() {
            let start = profile::real(&edges[0].start);
            let mut depth = 0;
            for (j, other) in self.loops.iter().enumerate() {
                if i != j
                    && profile::location(std::slice::from_ref(other), &start, dec)?
                        == Location::Inside
                {
                    depth += 1;
                }
            }
            if depth % 2 == 0 {
                outer += 1;
            }
        }
        Ok(outer == 1)
    }

    /// Sweeps the start or end cap of this extrusion by `thickness` along that cap's
    /// outward normal (`front`) or against it. `neighbour` names the source face across a
    /// profile edge. A shape outside the plain polygonal extrusion reports none.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn cap_sweep(
        &self,
        end: bool,
        thickness: &Rational,
        front: bool,
        cap_origins: &[Origin],
        neighbour: &dyn Fn(&[ElementId]) -> Option<Origin>,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<CapSweep>, Refusal> {
        let Sweep::Extrude {
            vector,
            through: false,
        } = &self.sweep
        else {
            return Ok(None);
        };
        if self.fillet.is_some()
            || self.chamfer.is_some()
            || self.loops.is_empty()
            || !self.simple_polygon(dec)?
        {
            return Ok(None);
        }
        let up = if vector.z.sign() == Sign::Positive {
            1
        } else {
            -1
        };
        // Local +z points away from the cap when the sweep runs along it.
        let outward = if end { up } else { -up };
        let along = if front { outward } else { -outward };
        let height = thickness * &unit(along);
        let (inverse, base_z) = if end {
            (
                self.inverse
                    .then(&Isometry::translation(vector.scaled(&unit(-1)))),
                Rational::zero(),
            )
        } else {
            (self.inverse.clone(), Rational::zero())
        };
        let placement = inverse.inverse();
        let z = |v: i64| RVector3::new(Rational::zero(), Rational::zero(), unit(v)).to_real();
        let at = |x: &Rational, y: &Rational, h: &Rational| -> Point3 {
            placement.apply_point(
                &catallax_shapes::RPoint3::new(x.clone(), y.clone(), h.clone()).to_real(),
            )
        };
        let x_axis = placement.apply_vector(
            &RVector3::new(Rational::one(), Rational::zero(), Rational::zero()).to_real(),
        );
        let origin = at(&Rational::zero(), &Rational::zero(), &base_z);
        let far = at(&Rational::zero(), &Rational::zero(), &height);
        let mut faces = Vec::new();
        let mut roles = Vec::new();
        let cap = |faces: &mut Vec<Face>,
                   roles: &mut Vec<Role>,
                   origin: Point3,
                   normal: Vector3,
                   names: Vec<Origin>,
                   role: Role,
                   dec: &mut Decider|
         -> Result<(), Refusal> {
            for name in names {
                faces.push(geometry::face(
                    Surface::Plane(Plane {
                        origin: origin.clone(),
                        normal: normal.clone(),
                        ref_dir: x_axis.clone(),
                    }),
                    name,
                    dec,
                )?);
                roles.push(role.clone());
            }
            Ok(())
        };
        cap(
            &mut faces,
            &mut roles,
            origin,
            placement.apply_vector(&z(-along)),
            cap_origins.to_vec(),
            Role::Start,
            dec,
        )?;
        cap(
            &mut faces,
            &mut roles,
            far,
            placement.apply_vector(&z(along)),
            cap_origins
                .iter()
                .map(|of| Origin {
                    step: step.clone(),
                    role: Role::Offset {
                        of: Box::new(of.clone()),
                    },
                })
                .collect(),
            Role::End,
            dec,
        )?;
        let first = cap_origins
            .first()
            .ok_or_else(|| bad("distance.planar-domain"))?;
        let rise = placement.apply_vector(&z(1));
        for Edge {
            ids, start, end, ..
        } in self.loops.iter().flatten()
        {
            let other = neighbour(ids).ok_or_else(|| bad("distance.planar-domain"))?;
            let name =
                EdgeName::new(first.clone(), other).map_err(|_| bad("distance.planar-domain"))?;
            // Outer loops run counter-clockwise and holes clockwise, so material lies left.
            let right = RVector3::new(&end[1] - &start[1], &start[0] - &end[0], Rational::zero());
            faces.push(face(
                Surface::Plane(Plane {
                    origin: at(&start[0], &start[1], &base_z),
                    normal: placement.apply_vector(&right.to_real()),
                    ref_dir: rise.clone(),
                }),
                step,
                Role::Rim { edge: name },
                dec,
            )?);
            roles.push(Role::Side {
                segment: ids[0].clone(),
                path: None,
            });
        }
        let mut constraints = Vec::new();
        let mut supports = Vec::new();
        for face in &faces {
            if !constraints
                .iter()
                .any(|c: &geometry::Constraint| c.polynomial.same_definition(&face.polynomial))
            {
                constraints.push(geometry::constraint(face.surface.clone())?);
                supports.push(Support::new(&face.surface));
            }
        }
        let mut feature = self.clone();
        feature.sweep = Sweep::Extrude {
            vector: RVector3::new(Rational::zero(), Rational::zero(), height.clone()),
            through: false,
        };
        feature.fast = FastProfile {
            matrix: geometry::fast_matrix(&inverse),
            translation: geometry::fast_translation(&inverse),
            sweep: [0.0, 0.0, height.to_f64_nearest(), 0.0, 0.0],
            valid: false,
            ..self.fast.clone()
        };
        feature.fast.valid = feature.fast.finite();
        feature.inverse = inverse;
        feature.constraints = constraints;
        feature.supports = supports;
        Ok(Some(CapSweep {
            feature,
            faces,
            roles,
        }))
    }
}
