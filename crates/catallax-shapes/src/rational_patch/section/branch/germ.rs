//! Exact outgoing acceleration from the retained patch and cutting plane.

use super::*;

fn missing(detail: &str) -> Refusal {
    Refusal::new("shapes.patch-section-germ", detail)
}

fn add(a: &Real, b: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
    reserve(a, b, false, dec)?;
    Ok(a + b)
}
fn sub(a: &Real, b: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
    add(a, &(-b), dec)
}
fn mul(a: &Real, b: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
    reserve(a, b, true, dec)?;
    Ok(a * b)
}
fn div(a: &Real, b: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
    reserve(a, b, true, dec)?;
    Ok(a / b)
}
fn reserve(a: &Real, b: &Real, product: bool, dec: &mut Decider) -> Result<(), Refusal> {
    if let (Some(a), Some(b)) = (a.rational_ref(), b.rational_ref()) {
        dec.reserve_exact_arithmetic(a, b, product)?;
    } else {
        dec.charge_work(1)?;
    }
    Ok(())
}
fn xyz(v: &Vector3) -> [&Real; 3] {
    [&v.x, &v.y, &v.z]
}
fn vector(values: Vec<Real>) -> Vector3 {
    let [x, y, z] = values.try_into().expect("three Cartesian components");
    Vector3::new(x, y, z)
}
fn dot(a: &Vector3, b: &Vector3, dec: &mut Decider) -> Result<Real, Refusal> {
    let mut result = Real::zero();
    for (a, b) in xyz(a).into_iter().zip(xyz(b)) {
        result = add(&result, &mul(a, b, dec)?, dec)?;
    }
    Ok(result)
}
fn scale(v: &Vector3, s: &Real, dec: &mut Decider) -> Result<Vector3, Refusal> {
    xyz(v)
        .into_iter()
        .map(|v| mul(v, s, dec))
        .collect::<Result<Vec<_>, _>>()
        .map(vector)
}
fn plus(a: &Vector3, b: &Vector3, dec: &mut Decider) -> Result<Vector3, Refusal> {
    xyz(a)
        .into_iter()
        .zip(xyz(b))
        .map(|(a, b)| add(a, b, dec))
        .collect::<Result<Vec<_>, _>>()
        .map(vector)
}
fn equal(a: &Vector3, b: &Vector3, dec: &mut Decider) -> Result<bool, Refusal> {
    for (a, b) in xyz(a).into_iter().zip(xyz(b)) {
        let difference = sub(a, b, dec)?;
        let sign = if difference.as_rational().is_none()
            && let Some(value) = dec.rational_identity(&difference)?
        {
            value.sign()
        } else {
            exact_sign(&difference, dec)?
        };
        if sign != Sign::Zero {
            return Ok(false);
        }
    }
    Ok(true)
}

fn value(
    poly: &RealPoly,
    at: &[Real; 2],
    order: [u32; 2],
    dec: &mut Decider,
) -> Result<Real, Refusal> {
    let mut sum = Real::zero();
    for (monomial, coefficient) in poly.terms() {
        if monomial
            .exponents()
            .iter()
            .zip(order)
            .any(|(&power, derivative)| power < derivative)
        {
            continue;
        }
        let mut term = coefficient;
        for axis in 0..2 {
            let mut power = monomial.exponents()[axis];
            for _ in 0..order[axis] {
                term = mul(&term, &Real::from_i64(i64::from(power)), dec)?;
                power -= 1;
            }
            for _ in 0..power {
                term = mul(&term, &at[axis], dec)?;
            }
        }
        sum = add(&sum, &term, dec)?;
    }
    Ok(sum)
}

struct Jet {
    point: Vector3,
    first: [Vector3; 2],
    second: [Vector3; 3],
}

fn position(span: &Span, at: &[Real; 2], dec: &mut Decider) -> Result<Vector3, Refusal> {
    let weight = value(&span.coordinates[3], at, [0, 0], dec)?;
    if exact_sign(&weight, dec)? != Sign::Positive {
        return Err(missing(
            "the retained patch needs a positive homogeneous weight",
        ));
    }
    span.coordinates[..3]
        .iter()
        .map(|p| div(&value(p, at, [0, 0], dec)?, &weight, dec))
        .collect::<Result<Vec<_>, _>>()
        .map(vector)
}

fn in_cell(domain: &Domain, at: &[Real; 2], dec: &mut Decider) -> Result<bool, Refusal> {
    for axis in 0..2 {
        let lower = sub(
            &at[axis],
            &Real::from_rational(domain[axis][0].clone()),
            dec,
        )?;
        let upper = sub(
            &Real::from_rational(domain[axis][1].clone()),
            &at[axis],
            dec,
        )?;
        if exact_sign(&lower, dec)? == Sign::Negative || exact_sign(&upper, dec)? == Sign::Negative
        {
            return Ok(false);
        }
    }
    Ok(true)
}
impl Jet {
    fn new(span: &Span, at: &[Real; 2], dec: &mut Decider) -> Result<Self, Refusal> {
        let mut h = Vec::new();
        for order in [[0, 0], [1, 0], [0, 1], [2, 0], [1, 1], [0, 2]] {
            h.push(
                span.coordinates
                    .iter()
                    .map(|p| value(p, at, order, dec))
                    .collect::<Result<Vec<_>, _>>()?,
            );
        }
        let w = &h[0][3];
        if exact_sign(w, dec)? != Sign::Positive {
            return Err(missing(
                "the retained patch needs a positive homogeneous weight",
            ));
        }
        let mut point = Vec::new();
        let mut first = [Vec::new(), Vec::new()];
        let mut second = [Vec::new(), Vec::new(), Vec::new()];
        for axis in 0..3 {
            let p = div(&h[0][axis], w, dec)?;
            let mut d = Vec::new();
            for (parameter, values) in first.iter_mut().enumerate() {
                let correction = mul(&p, &h[1 + parameter][3], dec)?;
                let derivative = div(&sub(&h[1 + parameter][axis], &correction, dec)?, w, dec)?;
                values.push(derivative.clone());
                d.push(derivative);
            }
            for (index, (a, b)) in [(0, 0), (0, 1), (1, 1)].into_iter().enumerate() {
                // H_ab = W_ab P + W_a P_b + W_b P_a + W P_ab.
                let correction = mul(&p, &h[3 + index][3], dec)?;
                let numerator = sub(&h[3 + index][axis], &correction, dec)?;
                let correction = mul(&d[b], &h[1 + a][3], dec)?;
                let numerator = sub(&numerator, &correction, dec)?;
                let correction = mul(&d[a], &h[1 + b][3], dec)?;
                second[index].push(div(&sub(&numerator, &correction, dec)?, w, dec)?);
            }
            point.push(p);
        }
        Ok(Self {
            point: vector(point),
            first: first.map(vector),
            second: second.map(vector),
        })
    }

    fn velocity(&self, tangent: &Vector3, dec: &mut Decider) -> Result<Option<[Real; 2]>, Refusal> {
        let [u, v] = &self.first;
        let uu = dot(u, u, dec)?;
        let uv = dot(u, v, dec)?;
        let vv = dot(v, v, dec)?;
        let determinant = sub(&mul(&uu, &vv, dec)?, &mul(&uv, &uv, dec)?, dec)?;
        if exact_sign(&determinant, dec)? != Sign::Positive {
            return Err(missing(
                "section acceleration needs a regular patch tangent frame",
            ));
        }
        let tu = dot(tangent, u, dec)?;
        let tv = dot(tangent, v, dec)?;
        let a = div(
            &sub(&mul(&tu, &vv, dec)?, &mul(&tv, &uv, dec)?, dec)?,
            &determinant,
            dec,
        )?;
        let b = div(
            &sub(&mul(&tv, &uu, dec)?, &mul(&tu, &uv, dec)?, dec)?,
            &determinant,
            dec,
        )?;
        let reconstructed = plus(&scale(u, &a, dec)?, &scale(v, &b, dec)?, dec)?;
        Ok(equal(&reconstructed, tangent, dec)?.then_some([a, b]))
    }

    fn acceleration(
        &self,
        velocity: &[Real; 2],
        normal: &Vector3,
        tangent: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vector3, Refusal> {
        let [a, b] = velocity;
        let aa = mul(a, a, dec)?;
        let ab = mul(&Real::from_i64(2), &mul(a, b, dec)?, dec)?;
        let bb = mul(b, b, dec)?;
        let remainder = plus(
            &plus(
                &scale(&self.second[0], &aa, dec)?,
                &scale(&self.second[1], &ab, dec)?,
                dec,
            )?,
            &scale(&self.second[2], &bb, dec)?,
            dec,
        )?;
        let [u, v] = &self.first;
        let nu = dot(normal, u, dec)?;
        let nv = dot(normal, v, dec)?;
        let tu = dot(tangent, u, dec)?;
        let tv = dot(tangent, v, dec)?;
        let determinant = sub(&mul(&nu, &tv, dec)?, &mul(&nv, &tu, dec)?, dec)?;
        if exact_sign(&determinant, dec)? == Sign::Zero {
            return Err(missing("the cutting plane needs a regular section tangent"));
        }
        // The plane and constant-speed conditions give n.A = t.A = 0.
        // A = R + S_u u'' + S_v v'' determines both parameter accelerations.
        let nr = -dot(normal, &remainder, dec)?;
        let tr = -dot(tangent, &remainder, dec)?;
        let a = div(
            &sub(&mul(&nr, &tv, dec)?, &mul(&nv, &tr, dec)?, dec)?,
            &determinant,
            dec,
        )?;
        let b = div(
            &sub(&mul(&nu, &tr, dec)?, &mul(&nr, &tu, dec)?, dec)?,
            &determinant,
            dec,
        )?;
        let acceleration = plus(
            &remainder,
            &plus(&scale(u, &a, dec)?, &scale(v, &b, dec)?, dec)?,
            dec,
        )?;
        if exact_sign(&dot(normal, &acceleration, dec)?, dec)? != Sign::Zero
            || exact_sign(&dot(tangent, &acceleration, dec)?, dec)? != Sign::Zero
        {
            return Err(missing(
                "the section acceleration equations have not been proved",
            ));
        }
        Ok(acceleration)
    }
}

impl PatchSectionBranch {
    /// Exact normal acceleration for the supplied nonzero outgoing velocity.
    /// Plane sections use the retained homogeneous patch derivatives.
    /// Every retained chart preimage is checked before the outgoing cell sides
    /// select the result at knots, seams and cell ends.
    /// Other cutting surfaces and unproved selectors refuse by name.
    pub fn normal_acceleration(
        &self,
        point: &Point3,
        tangent: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vector3, Refusal> {
        let SectionSurface::Plane(plane) = &self.surface else {
            return Err(missing(
                "this cutting surface needs an exact section-germ Hessian query",
            ));
        };
        if exact_sign(&dot(tangent, tangent, dec)?, dec)? != Sign::Positive {
            return Err(missing(
                "section acceleration needs a nonzero outgoing tangent",
            ));
        }
        if exact_sign(&dot(&plane.normal, tangent, dec)?, dec)? != Sign::Zero {
            return Err(missing("the outgoing tangent is not in the section plane"));
        }
        let position = Vector3::new(point.x.clone(), point.y.clone(), point.z.clone());
        let origin = Vector3::new(
            plane.origin.x.clone(),
            plane.origin.y.clone(),
            plane.origin.z.clone(),
        );
        if exact_sign(
            &sub(
                &dot(&plane.normal, &position, dec)?,
                &dot(&plane.normal, &origin, dec)?,
                dec,
            )?,
            dec,
        )? != Sign::Zero
        {
            return Err(missing(
                "the acceleration point is outside the section plane",
            ));
        }
        let mut candidates = Vec::new();
        for (chart, patch) in self.charts.iter().enumerate() {
            if self.cells.iter().all(|cell| cell.chart != chart) {
                continue;
            }
            // A certified patch is injective. An exact cell end therefore
            // supplies its unique chart preimage without a second inversion.
            let mut known = None;
            'ends: for cell in self.cells.iter().filter(|cell| cell.chart == chart) {
                for at in &cell.ends {
                    if equal(&self::position(self.span(cell), at, dec)?, &position, dec)? {
                        known = Some(at.clone());
                        break 'ends;
                    }
                }
            }
            let at = match known {
                Some(at) => Some(at),
                None => patch.parameters_at(point, dec)?,
            };
            let Some(at) = at else { continue };
            for (index, cell) in self
                .cells
                .iter()
                .enumerate()
                .filter(|(_, cell)| cell.chart == chart)
            {
                if !in_cell(&cell.domain, &at, dec)? {
                    continue;
                }
                let free = 1 - cell.solved;
                let end = if exact_sign(&sub(&at[free], &cell.ends[0][free], dec)?, dec)?
                    == Sign::Zero
                {
                    Some(0)
                } else if exact_sign(&sub(&at[free], &cell.ends[1][free], dec)?, dec)? == Sign::Zero
                {
                    Some(1)
                } else {
                    None
                };
                candidates.push((index, at.clone(), end));
            }
        }
        let mut result = None;
        for (index, at, end) in candidates {
            let cell = &self.cells[index];
            let jet = Jet::new(self.span(cell), &at, dec)?;
            if !equal(&jet.point, &position, dec)? {
                return Err(missing(
                    "the retained cell join needs exact point correspondence",
                ));
            }
            let Some(velocity) = jet.velocity(tangent, dec)? else {
                continue;
            };
            let free = 1 - cell.solved;
            let traversal = sub(&cell.ends[1][free], &cell.ends[0][free], dec)?;
            let along = exact_sign(&mul(&velocity[free], &traversal, dec)?, dec)?;
            if along == Sign::Zero {
                return Err(missing(
                    "the outgoing tangent needs a regular free-parameter direction",
                ));
            }
            if matches!(
                (end, along),
                (Some(0), Sign::Negative) | (Some(1), Sign::Positive)
            ) {
                continue;
            }
            let acceleration = jet.acceleration(&velocity, &plane.normal, tangent, dec)?;
            if let Some(previous) = &result {
                if !equal(previous, &acceleration, dec)? {
                    return Err(missing(
                        "the tangent selects unequal cell-side accelerations; a branch selector is missing",
                    ));
                }
            } else {
                result = Some(acceleration);
            }
        }
        result.ok_or_else(|| missing("the tangent does not select a regular outgoing cell side"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BSplineSurface, Sphere};
    use catallax_numbers::with_operation_budget;

    fn q(n: i64, d: i64) -> Real {
        Rational::from_ratio(n, d).unwrap().into()
    }
    fn point(x: i64, y: i64, z: i64) -> Point3 {
        Point3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
    }
    fn plane(axis: usize, coordinate: Real) -> SectionSurface {
        let mut origin = [Real::zero(), Real::zero(), Real::zero()];
        origin[axis] = coordinate;
        let mut normal = [Real::zero(), Real::zero(), Real::zero()];
        normal[axis] = Real::one();
        let [x, y, z] = origin;
        let [a, b, c] = normal;
        SectionSurface::Plane(Plane {
            origin: Point3::new(x, y, z),
            normal: Vector3::new(a, b, c),
            ref_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
        })
    }
    fn knots(values: &[i64]) -> Vec<Rational> {
        values.iter().copied().map(Rational::from_i64).collect()
    }
    fn section(
        patch: &RationalPatch,
        cut: &SectionSurface,
        dec: &mut Decider,
    ) -> PatchSectionBranch {
        let found = super::super::search(patch, cut, dec).unwrap();
        assert_eq!(found.branches().len(), 1);
        found.branches()[0].clone()
    }
    fn assert_equal(a: &Vector3, b: &Vector3, dec: &mut Decider) {
        assert!(equal(a, b, dec).unwrap(), "{a:?} != {b:?}");
    }
    fn circle_patch(dec: &mut Decider) -> RationalPatch {
        let controls = [(q(1, 1), q(0, 1)), (q(1, 1), q(1, 2)), (q(3, 5), q(4, 5))]
            .into_iter()
            .flat_map(|(x, y)| {
                [
                    Point3::new(x.clone(), y.clone(), q(0, 1)),
                    Point3::new(x, y, q(1, 1)),
                ]
            })
            .collect();
        RationalPatch::new(
            BSplineSurface::new(
                [2, 1],
                [knots(&[0, 0, 0, 1, 1, 1]), knots(&[0, 0, 1, 1])],
                [3, 2],
                controls,
                vec![q(1, 1), q(1, 1), q(1, 1), q(1, 1), q(5, 4), q(5, 4)],
                dec,
            )
            .unwrap(),
            dec,
        )
        .unwrap()
    }

    #[test]
    fn rational_circle_section_uses_weight_derivatives_and_outgoing_end_sides() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let patch = circle_patch(&mut dec);
            let branch = section(&patch, &plane(2, q(1, 2)), &mut dec);
            for (u, speed) in [(q(0, 1), 2), (q(1, 2), 2), (q(1, 2), -2), (q(1, 1), -2)] {
                let p = patch.point_at(&[u, q(1, 2)], &mut dec).unwrap();
                let tangent = Vector3::new(
                    -&p.y * Real::from_i64(speed),
                    &p.x * &Real::from_i64(speed),
                    Real::zero(),
                );
                let acceleration = branch.normal_acceleration(&p, &tangent, &mut dec).unwrap();
                assert_equal(
                    &acceleration,
                    &Vector3::new(
                        &p.x * &Real::from_i64(-4),
                        &p.y * &Real::from_i64(-4),
                        Real::zero(),
                    ),
                    &mut dec,
                );
            }
        });
        assert!(!work.exhausted, "{} units", work.used);
    }

    #[test]
    fn closed_circle_cell_joins_and_reversed_cover_have_the_same_acceleration() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let mut controls = Vec::new();
            for i in 0..3 {
                for j in 0..3 {
                    let square = [1, -1, 1];
                    controls.push(point(i - 1, j - 1, square[i as usize] + square[j as usize]));
                }
            }
            let patch = RationalPatch::new(
                BSplineSurface::new(
                    [2, 2],
                    [knots(&[-1, -1, -1, 1, 1, 1]), knots(&[-1, -1, -1, 1, 1, 1])],
                    [3, 3],
                    controls,
                    vec![Real::one(); 9],
                    &mut dec,
                )
                .unwrap(),
                &mut dec,
            )
            .unwrap();
            let branch = section(&patch, &plane(2, q(1, 3)), &mut dec);
            assert!(branch.is_closed());
            let mut points = vec![Point3::new(
                q(1, 3),
                &Real::from_i64(2).sqrt() * &q(1, 3),
                q(1, 3),
            )];
            points.push(branch.endpoints()[0].clone());
            points.extend(branch.witnesses().take(3));
            for branch in [branch.clone(), branch.reversed()] {
                for p in &points {
                    for speed in [1, -2] {
                        let tangent = Vector3::new(
                            -&p.y * Real::from_i64(speed),
                            &p.x * &Real::from_i64(speed),
                            Real::zero(),
                        );
                        let acceleration =
                            branch.normal_acceleration(p, &tangent, &mut dec).unwrap();
                        let scale = Real::from_i64(-speed * speed);
                        assert_equal(
                            &acceleration,
                            &Vector3::new(&p.x * &scale, &p.y * &scale, Real::zero()),
                            &mut dec,
                        );
                    }
                }
            }
        });
        assert!(!work.exhausted, "{} units", work.used);
    }

    #[test]
    fn joined_cells_select_their_different_second_derivatives() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let controls = [q(0, 1), q(1, 2), q(1, 1), q(3, 2), q(2, 1)]
                .into_iter()
                .enumerate()
                .flat_map(|(i, x)| {
                    let z = Real::from_i64(i64::from(i == 4));
                    [
                        Point3::new(x.clone(), Real::zero(), z.clone()),
                        Point3::new(x, Real::one(), z),
                    ]
                })
                .collect();
            let patch = RationalPatch::new(
                BSplineSurface::new(
                    [2, 1],
                    [knots(&[0, 0, 0, 1, 1, 2, 2, 2]), knots(&[0, 0, 1, 1])],
                    [5, 2],
                    controls,
                    vec![Real::one(); 10],
                    &mut dec,
                )
                .unwrap(),
                &mut dec,
            )
            .unwrap();
            let branch = section(&patch, &plane(1, q(1, 2)), &mut dec);
            let p = Point3::new(Real::one(), q(1, 2), Real::zero());
            for branch in [branch.clone(), branch.reversed()] {
                for (direction, expected) in [(1, 2), (-1, 0)] {
                    let tangent =
                        Vector3::new(Real::from_i64(direction), Real::zero(), Real::zero());
                    let acceleration = branch.normal_acceleration(&p, &tangent, &mut dec).unwrap();
                    assert_equal(
                        &acceleration,
                        &Vector3::new(Real::zero(), Real::zero(), Real::from_i64(expected)),
                        &mut dec,
                    );
                }
            }
        });
        assert!(!work.exhausted, "{} units", work.used);
    }

    #[test]
    fn invalid_directions_and_zero_budget_do_not_publish_an_acceleration() {
        let ((branch, p), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let patch = circle_patch(&mut dec);
            let branch = section(&patch, &plane(2, q(1, 2)), &mut dec);
            let p = patch.point_at(&[q(1, 2), q(1, 2)], &mut dec).unwrap();
            for tangent in [
                Vector3::new(Real::zero(), Real::zero(), Real::one()),
                Vector3::new(Real::zero(), Real::zero(), Real::zero()),
                Vector3::new(Real::one(), Real::zero(), Real::zero()),
            ] {
                assert_eq!(
                    branch
                        .normal_acceleration(&p, &tangent, &mut dec)
                        .unwrap_err()
                        .code,
                    "shapes.patch-section-germ"
                );
            }
            (branch, p)
        });
        assert!(!work.exhausted, "{} units", work.used);
        let (answer, stopped) = with_operation_budget(0, || {
            branch.normal_acceleration(
                &p,
                &Vector3::new(-&p.y, p.x.clone(), Real::zero()),
                &mut Decider::default(),
            )
        });
        assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
        assert!(stopped.exhausted);
    }

    #[test]
    fn a_nonplane_section_names_the_missing_hessian_query() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let patch = RationalPatch::new(
                BSplineSurface::new(
                    [1, 1],
                    [knots(&[0, 0, 1, 1]), knots(&[0, 0, 1, 1])],
                    [2, 2],
                    vec![
                        point(0, 0, 0),
                        point(0, 1, 0),
                        point(1, 0, 0),
                        point(1, 1, 0),
                    ],
                    vec![Real::one(); 4],
                    &mut dec,
                )
                .unwrap(),
                &mut dec,
            )
            .unwrap();
            let branch = section(
                &patch,
                &SectionSurface::Sphere(Sphere {
                    center: Point3::new(q(1, 2), q(1, 2), Real::zero()),
                    radius: q(1, 3),
                    axis: Vector3::new(Real::zero(), Real::zero(), Real::one()),
                    ref_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
                }),
                &mut dec,
            );
            assert_eq!(
                branch
                    .normal_acceleration(
                        &branch.endpoints()[0],
                        &Vector3::new(Real::one(), Real::zero(), Real::zero()),
                        &mut dec
                    )
                    .unwrap_err()
                    .detail,
                "this cutting surface needs an exact section-germ Hessian query"
            );
        });
        assert!(!work.exhausted, "{} units", work.used);
    }
}
