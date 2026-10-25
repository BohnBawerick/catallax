//! Exact finite cap parameters over the independent rational extrusion profile.
//!
//! The rational sweep direction remains a projective carrier. Only its finite
//! parameter interval needs a radical; profile projection cancels that scale.
use super::*;

impl Feature {
    pub(super) fn algebraic_extrude(
        profile: &Profile,
        step: &StepId,
        distance: &Rational,
        direction: &RVector3,
        dec: &mut Decider,
    ) -> Result<Option<(Self, Vec<Face>)>, Refusal> {
        let end = &Real::one() / &profile::positive_sqrt(&direction.norm_squared(), dec)?;
        let sweep = direction.scaled(distance);
        let input = catallax_recipe::Vector3([
            sweep.x.clone().into(),
            sweep.y.clone().into(),
            sweep.z.clone().into(),
        ]);
        // Reuse only the independent through-profile definition and its rational
        // analytic reductions. No B-rep body or topology is built here.
        let Some((tool, mut faces)) = Self::extrude(profile, step, None, Some(&input), false, dec)?
        else {
            return Ok(None);
        };
        let normal = profile::rational_vector(&profile.pose, &z()).to_real();
        let reference = profile::rational_vector(
            &profile.pose,
            &RVector3::new(Rational::one(), Rational::zero(), Rational::zero()),
        )
        .to_real();
        let origin = profile::rational_point(&profile.pose, &RPoint3::origin()).to_real();
        let offset = sweep.to_real().scaled(&end);
        for (origin, role) in [
            (origin.clone(), Role::Start),
            (origin.plus(&offset), Role::End),
        ] {
            faces.push(face(
                Surface::Plane(Plane {
                    origin,
                    normal: normal.clone(),
                    ref_dir: reference.clone(),
                }),
                step,
                role,
                dec,
            )?);
        }
        let Sweep::Extrude { vector, .. } = tool.sweep else {
            unreachable!()
        };
        let mut result = Self::new(
            profile,
            Sweep::Extrude {
                vector,
                through: false,
            },
            &faces,
            Rational::zero(),
            Vec::new(),
        )?;
        result.algebraic_end = Some(end.clone());
        let enclosed = geometry::enclose(&end, 192, dec)?;
        let approximate = enclosed.mid().to_f64_nearest();
        for value in &mut result.fast.sweep[..3] {
            *value *= approximate;
        }
        result.fast.valid = result.fast.finite();
        Ok(Some((result, faces)))
    }

    pub(super) fn extrusion_end(&self) -> Real {
        self.algebraic_end.clone().unwrap_or_else(Real::one)
    }

    pub(super) fn extrusion_vector(&self, vector: &RVector3) -> Vector3 {
        if let Some(end) = &self.algebraic_end {
            vector.to_real().scaled(end)
        } else {
            vector.to_real()
        }
    }

    pub(super) fn algebraic_display_box(
        &self,
        dec: &mut Decider,
    ) -> Result<Option<[Ball; 3]>, Refusal> {
        let Sweep::Extrude {
            vector,
            through: false,
        } = &self.sweep
        else {
            return Ok(None);
        };
        let profile = self.sample_profile_bounds();
        let mut bounds = [
            profile[0].clone(),
            profile[1].clone(),
            [Rational::zero(), Rational::zero()],
        ];
        let v = self.extrusion_vector(vector);
        for (i, v) in [&v.x, &v.y, &v.z].into_iter().enumerate() {
            let b = geometry::enclose(v, 192, dec)?;
            bounds[i][0] = &bounds[i][0] + &b.lo().min(Rational::zero());
            bounds[i][1] = &bounds[i][1] + &b.hi().max(Rational::zero());
        }
        Ok(Some(geometry::display_box(
            &self.inverse.inverse(),
            &bounds,
            dec,
        )?))
    }
}

impl Feature {
    /// An upright reference prism with the same independently validated profile.
    /// The finite sweep's invertible affine image is exactly this feature.
    pub(crate) fn algebraic_display_source(
        &self,
        step: &StepId,
        height: &Rational,
        dec: &mut Decider,
    ) -> Result<Option<(Self, Vec<Face>)>, Refusal> {
        if self.algebraic_end.is_none() || self.fillet.is_some() || self.chamfer.is_some() {
            return Ok(None);
        }
        let profile = Profile {
            pose: catallax_shapes::Pose::identity(),
            loops: self.loops.clone(),
        };
        Self::extrude(&profile, step, Some(height), None, false, dec)
    }

    pub(crate) fn algebraic_display_affine(&self, height: &Rational) -> (Point3, [Vector3; 3]) {
        let Sweep::Extrude { vector, .. } = &self.sweep else {
            unreachable!()
        };
        let placement = self.inverse.inverse();
        let sweep = self
            .extrusion_vector(vector)
            .scaled(&Real::from(height.recip().expect("positive source height")));
        let x = RVector3::new(Rational::one(), Rational::zero(), Rational::zero()).to_real();
        let y = RVector3::new(Rational::zero(), Rational::one(), Rational::zero()).to_real();
        (
            placement.apply_point(&Point3::origin()),
            [x, y, sweep].map(|v| placement.apply_vector(&v)),
        )
    }

    /// Frobenius norm bounds the operator norm of [e_x,e_y,sweep/height]. The final
    /// placement is an isometry, including reflection, and does not enlarge it.
    pub(crate) fn algebraic_display_lipschitz(
        &self,
        height: &Rational,
        dec: &mut Decider,
    ) -> Result<Rational, Refusal> {
        let Sweep::Extrude { vector, .. } = &self.sweep else {
            unreachable!()
        };
        let physical = self.extrusion_vector(vector);
        let squared =
            &Real::from_i64(2) + &(&physical.norm_squared() / &Real::from(height * height));
        Ok(geometry::enclose(&squared.sqrt(), 192, dec)?.hi())
    }

    /// A rational reference height near the physical sweep length keeps affine
    /// error amplification bounded for both long and short extrusions.
    pub(crate) fn algebraic_display_height(
        &self,
        dec: &mut Decider,
    ) -> Result<Option<Rational>, Refusal> {
        if self.algebraic_end.is_none() || self.fillet.is_some() || self.chamfer.is_some() {
            return Ok(None);
        }
        let Sweep::Extrude { vector, .. } = &self.sweep else {
            unreachable!()
        };
        let norm = self.extrusion_vector(vector).norm_squared().sqrt();
        let upper = geometry::enclose(&norm, 192, dec)?.hi();
        Ok(Some(
            Rational::from_integer(upper.floor()).max(Rational::one()),
        ))
    }

    pub(crate) fn algebraic_display_reversing(&self) -> bool {
        let Sweep::Extrude { vector, .. } = &self.sweep else {
            unreachable!()
        };
        (vector.z.sign() == Sign::Negative) ^ self.inverse.reversing()
    }
}
