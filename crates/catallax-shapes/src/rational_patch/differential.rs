//! Charged homogeneous derivatives shared by retained intersection germs.

use super::*;

fn missing(detail: &str) -> Refusal {
    Refusal::new("shapes.patch-section-germ", detail)
}
pub(super) fn add(a: &Real, b: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
    reserve(a, b, false, dec)?;
    Ok(a + b)
}
pub(super) fn sub(a: &Real, b: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
    add(a, &(-b), dec)
}
pub(super) fn mul(a: &Real, b: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
    reserve(a, b, true, dec)?;
    Ok(a * b)
}
pub(super) fn div(a: &Real, b: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
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
pub(super) fn vector(values: Vec<Real>) -> Vector3 {
    let [x, y, z] = values.try_into().expect("three Cartesian components");
    Vector3::new(x, y, z)
}
pub(super) fn dot(a: &Vector3, b: &Vector3, dec: &mut Decider) -> Result<Real, Refusal> {
    let mut result = Real::zero();
    for (a, b) in xyz(a).into_iter().zip(xyz(b)) {
        result = add(&result, &mul(a, b, dec)?, dec)?;
    }
    Ok(result)
}
pub(super) fn scale(v: &Vector3, s: &Real, dec: &mut Decider) -> Result<Vector3, Refusal> {
    xyz(v)
        .into_iter()
        .map(|v| mul(v, s, dec))
        .collect::<Result<Vec<_>, _>>()
        .map(vector)
}
pub(super) fn plus(a: &Vector3, b: &Vector3, dec: &mut Decider) -> Result<Vector3, Refusal> {
    xyz(a)
        .into_iter()
        .zip(xyz(b))
        .map(|(a, b)| add(a, b, dec))
        .collect::<Result<Vec<_>, _>>()
        .map(vector)
}
pub(super) fn equal(a: &Vector3, b: &Vector3, dec: &mut Decider) -> Result<bool, Refusal> {
    for (a, b) in xyz(a).into_iter().zip(xyz(b)) {
        if exact_sign(&sub(a, b, dec)?, dec)? != Sign::Zero {
            return Ok(false);
        }
    }
    Ok(true)
}

pub(super) fn value(
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

pub(super) struct Jet {
    pub(super) point: Vector3,
    first: [Vector3; 2],
    second: [Vector3; 3],
}
impl Jet {
    pub(super) fn new(span: &Span, at: &[Real; 2], dec: &mut Decider) -> Result<Self, Refusal> {
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

    pub(super) fn velocity(
        &self,
        tangent: &Vector3,
        dec: &mut Decider,
    ) -> Result<Option<[Real; 2]>, Refusal> {
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

    pub(super) fn acceleration(
        &self,
        velocity: &[Real; 2],
        normal: &Vector3,
        tangent: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vector3, Refusal> {
        let remainder = self.second_along(velocity, dec)?;
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

    pub(super) fn normal(&self, dec: &mut Decider) -> Result<Vector3, Refusal> {
        let [a, b] = &self.first;
        Ok(Vector3::new(
            sub(&mul(&a.y, &b.z, dec)?, &mul(&a.z, &b.y, dec)?, dec)?,
            sub(&mul(&a.z, &b.x, dec)?, &mul(&a.x, &b.z, dec)?, dec)?,
            sub(&mul(&a.x, &b.y, dec)?, &mul(&a.y, &b.x, dec)?, dec)?,
        ))
    }

    pub(super) fn second_along(
        &self,
        velocity: &[Real; 2],
        dec: &mut Decider,
    ) -> Result<Vector3, Refusal> {
        let [a, b] = velocity;
        let aa = mul(a, a, dec)?;
        let ab = mul(&Real::from_i64(2), &mul(a, b, dec)?, dec)?;
        let bb = mul(b, b, dec)?;
        plus(
            &plus(
                &scale(&self.second[0], &aa, dec)?,
                &scale(&self.second[1], &ab, dec)?,
                dec,
            )?,
            &scale(&self.second[2], &bb, dec)?,
            dec,
        )
    }
}
