//! Exact eliminants enumerate finite preimages; original equations keep branches.

use super::*;

mod affine;
mod plane;

fn univariate(poly: &RealPoly, axis: usize, dec: &mut Decider) -> Result<RealPoly, Refusal> {
    let mut terms = Vec::new();
    for (monomial, coefficient) in poly.terms() {
        let exponents = monomial.exponents();
        if exponents
            .iter()
            .enumerate()
            .any(|(i, &power)| i != axis && power != 0)
        {
            if exact_sign(&coefficient, dec)? == Sign::Zero {
                continue;
            }
            return Err(unproven("an eliminant retained another parameter"));
        }
        terms.push((vec![exponents[axis]], coefficient));
    }
    Ok(RealPoly::from_terms(1, terms))
}

fn within(value: &Real, bounds: &[Rational; 2], dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(
        exact_sign(&(value - &Real::from_rational(bounds[0].clone())), dec)? != Sign::Negative
            && exact_sign(&(&Real::from_rational(bounds[1].clone()) - value), dec)?
                != Sign::Negative,
    )
}

fn restrict(poly: &RealPoly, value: &Real, dec: &mut Decider) -> Result<RealPoly, Refusal> {
    univariate(
        &poly.substitute(0, &RealPoly::constant(2, value.clone())),
        1,
        dec,
    )
}

pub(super) fn roots(
    equations: &[RealPoly; 2],
    domain: &Domain,
    dec: &mut Decider,
) -> Result<Vec<[Real; 2]>, Refusal> {
    dec.check_deadline()?;
    let bits = dec.policy().max_precision_bits.min(128);
    if bits == 0 {
        return Err(unproven(
            "no precision is available for a finite root query",
        ));
    }
    let enclosure = box_balls(domain, bits);
    for equation in equations {
        if !equation.eval_ball(&enclosure, dec)?.contains_zero() {
            return Ok(Vec::new());
        }
        if equation.total_degree(dec)? == Some(0) {
            return Ok(Vec::new());
        }
    }
    let [f, g] = equations;
    let degree_v = |poly: &RealPoly, dec: &mut Decider| -> Result<u32, Refusal> {
        let mut degree = 0;
        for (monomial, coefficient) in poly.terms() {
            if monomial.exponents()[1] > degree && exact_sign(&coefficient, dec)? != Sign::Zero {
                degree = monomial.exponents()[1];
            }
        }
        Ok(degree)
    };
    if degree_v(f, dec)? == 0 && degree_v(g, dec)? == 0 {
        let common = univariate(f, 0, dec)?.gcd(&univariate(g, 0, dec)?, dec)?;
        if common.is_zero(dec)? {
            return Err(unproven(
                "the finite root query has a nonisolated parameter family",
            ));
        }
        for value in common.real_roots(dec)? {
            if within(&value, &domain[0], dec)? {
                return Err(unproven(
                    "the finite root query contains a whole parameter row",
                ));
            }
        }
        return Ok(Vec::new());
    }
    let resultant = f.resultant(g, 1, dec)?;
    if resultant.is_zero(dec)? {
        return Err(unproven(
            "the finite root query has no nonzero completeness eliminant",
        ));
    }
    let eliminant = univariate(&resultant, 0, dec)?;
    let mut found = Vec::new();
    for u in eliminant.real_roots(dec)? {
        if !within(&u, &domain[0], dec)? {
            continue;
        }
        let a = restrict(f, &u, dec)?;
        let b = restrict(g, &u, dec)?;
        let common = a.gcd(&b, dec)?;
        if common.is_zero(dec)? {
            return Err(unproven(
                "an eliminant root contains a nonisolated parameter row",
            ));
        }
        for v in common.real_roots(dec)? {
            if !within(&v, &domain[1], dec)? {
                continue;
            }
            let at = [u.clone(), v];
            if equations.iter().try_fold(true, |equal, equation| {
                Ok::<_, Refusal>(equal && exact_sign(&equation.eval_real(&at), dec)? == Sign::Zero)
            })? {
                found.push(at);
            }
        }
    }
    Ok(found)
}

fn projected(span: &Span, point: &Point3, projection: &[Vector3; 2]) -> [RealPoly; 2] {
    let values = [&point.x, &point.y, &point.z];
    let residuals: [RealPoly; 3] =
        std::array::from_fn(|i| &span.coordinates[i] - &span.coordinates[3].scale(values[i]));
    projection.each_ref().map(|direction| {
        let coefficients = [&direction.x, &direction.y, &direction.z];
        (0..3).fold(RealPoly::zero(2), |sum, i| {
            sum + residuals[i].scale(coefficients[i])
        })
    })
}

impl RationalPatch {
    fn projected_parameter(
        &self,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Option<[Real; 2]>, Refusal> {
        math::validate_point(point, dec)?;
        let mut found: Option<[Real; 2]> = None;
        for span in &self.data.spans {
            let equations = projected(span, point, &self.data.projection);
            for at in roots(&equations, &span.domain, dec)? {
                if let Some(previous) = &found {
                    for axis in 0..2 {
                        if exact_sign(&(&at[axis] - &previous[axis]), dec)? != Sign::Zero {
                            return Err(unproven(
                                "a checked embedding acquired distinct projected preimages",
                            ));
                        }
                    }
                } else {
                    found = Some(at);
                }
            }
        }
        Ok(found)
    }

    /// The unique exact finite preimage, or proven absence. Projection roots are
    /// checked against all three world coordinates, never accepted by distance.
    pub fn parameters_at(
        &self,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Option<[Real; 2]>, Refusal> {
        match affine::parameters(self, point, dec)? {
            affine::Preimage::Point(at) => return Ok(Some(at)),
            affine::Preimage::Absent => return Ok(None),
            affine::Preimage::Unsupported => {}
        }
        let Some(at) = self.projected_parameter(point, dec)? else {
            return Ok(None);
        };
        let actual = self.point_at(&at, dec)?;
        Ok((exact_sign(&actual.minus(point).norm_squared(), dec)? == Sign::Zero).then_some(at))
    }

    /// Natural chart normal at a genuine finite point. At a crease, incompatible
    /// one-sided normals refuse rather than choosing an arbitrary knot branch.
    pub fn normal_at(&self, point: &Point3, dec: &mut Decider) -> Result<Vector3, Refusal> {
        let at = self
            .parameters_at(point, dec)?
            .ok_or_else(|| domain_error("normal point is not on the finite patch"))?;
        let mut normal: Option<Vector3> = None;
        for span in &self.data.spans {
            if !contains(&span.domain, &at, dec)? {
                continue;
            }
            let [u, v] = span.derivatives(&at);
            let next = u.cross(&v);
            if let Some(previous) = &normal {
                if exact_sign(&previous.cross(&next).norm_squared(), dec)? != Sign::Zero
                    || exact_sign(&previous.dot(&next), dec)? != Sign::Positive
                {
                    return Err(domain_error(
                        "the finite knot join has incompatible one-sided normals",
                    ));
                }
            } else {
                normal = Some(next);
            }
        }
        normal.ok_or_else(|| domain_error("normal has no finite span"))
    }

    /// Oriented side only over the checked finite projection image. An outside
    /// projection has no side and returns a domain refusal, not a global sign.
    pub fn side_at(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        let at = self.projected_parameter(point, dec)?.ok_or_else(|| {
            domain_error("the side point lies outside the finite projection image")
        })?;
        let actual = self.point_at(&at, dec)?;
        let normal = self.data.projection[0].cross(&self.data.projection[1]);
        exact_sign(&point.minus(&actual).dot(&normal), dec)
    }

    /// Complete transverse intersections of a whole line with every finite span.
    /// Eliminant candidates are checked in the original equations. Tangencies
    /// and contained curve families retain a refusal, never a false miss.
    pub fn line_hits(
        &self,
        origin: &Point3,
        direction: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vec<Point3>, Refusal> {
        math::validate_point(origin, dec)?;
        math::validate_vector(direction, dec)?;
        let coordinates = [&direction.x, &direction.y, &direction.z];
        let mut pivot = None;
        for (axis, value) in coordinates.iter().enumerate() {
            if exact_sign(value, dec)? != Sign::Zero {
                pivot = Some(axis);
                break;
            }
        }
        let pivot = pivot.ok_or_else(|| domain_error("a line direction must be nonzero"))?;
        let others = [[1, 2], [0, 2], [0, 1]][pivot];
        let values = [&origin.x, &origin.y, &origin.z];
        let mut hits: Vec<(Real, Point3)> = Vec::new();
        for span in &self.data.spans {
            let residuals: [RealPoly; 3] = std::array::from_fn(|i| {
                &span.coordinates[i] - &span.coordinates[3].scale(values[i])
            });
            let equations = std::array::from_fn(|i| {
                &residuals[others[i]].scale(coordinates[pivot])
                    - &residuals[pivot].scale(coordinates[others[i]])
            });
            for at in roots(&equations, &span.domain, dec)? {
                let point = span.point(&at);
                let [u, v] = span.derivatives(&at);
                if exact_sign(&u.cross(&v).dot(direction), dec)? == Sign::Zero {
                    return Err(unproven("the line is tangent to the finite rational patch"));
                }
                let xyz = [&point.x, &point.y, &point.z];
                let station = &(xyz[pivot] - values[pivot]) / coordinates[pivot];
                let expected = origin.plus(&direction.scaled(&station));
                if exact_sign(&expected.minus(&point).norm_squared(), dec)? != Sign::Zero {
                    return Err(unproven(
                        "a projected line root failed its original three-dimensional equations",
                    ));
                }
                let mut index = hits.len();
                let mut duplicate = false;
                for (i, (previous, _)) in hits.iter().enumerate() {
                    match exact_sign(&(&station - previous), dec)? {
                        Sign::Negative => {
                            index = i;
                            break;
                        }
                        Sign::Zero => {
                            duplicate = true;
                            break;
                        }
                        Sign::Positive => (),
                    }
                }
                if !duplicate {
                    hits.insert(index, (station, point));
                }
            }
        }
        Ok(hits.into_iter().map(|(_, point)| point).collect())
    }

    /// Conservative parameter cover of every preimage of a spatial box. A cover
    /// is not incidence. `None` proves the whole finite patch misses the box.
    pub fn parameters_of_box(
        &self,
        point: &[Ball; 3],
        dec: &Decider,
    ) -> Result<Option<[Ball; 2]>, Refusal> {
        let bits = point.iter().map(Ball::precision).max().unwrap_or(1);
        if bits > dec.policy().max_precision_bits {
            return Err(unproven("inverse box precision exceeds caller policy"));
        }
        let mut pending: Vec<_> = self
            .data
            .spans
            .iter()
            .map(|span| (span.domain.clone(), 0u32))
            .collect();
        let mut cover: Option<Domain> = None;
        let mut visited = 0usize;
        while let Some((cell, depth)) = pending.pop() {
            dec.check_deadline()?;
            visited += 1;
            let image =
                self.carrier()
                    .evaluate_box(cell[0].each_ref(), cell[1].each_ref(), bits, dec)?;
            if (0..3).any(|i| image[i].hi() < point[i].lo() || image[i].lo() > point[i].hi()) {
                continue;
            }
            if depth >= 12 || visited >= MAX_PROOF_CELLS {
                match &mut cover {
                    Some(cover) => {
                        for axis in 0..2 {
                            cover[axis][0] = cover[axis][0].clone().min(cell[axis][0].clone());
                            cover[axis][1] = cover[axis][1].clone().max(cell[axis][1].clone());
                        }
                    }
                    None => cover = Some(cell),
                }
            } else {
                let [low, high] = split(&cell, &self.data.widths);
                pending.push((high, depth + 1));
                pending.push((low, depth + 1));
            }
        }
        Ok(cover.map(|cover| box_balls(&cover, bits)))
    }
}
