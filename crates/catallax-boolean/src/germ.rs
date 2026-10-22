//! Infinitesimal radial order from exact implicit derivatives. Equal tangents
//! and curvatures are separated by higher coefficients, never nearby samples.

use crate::geometry::{oriented, product, refusal, sign};
use catallax_intersect::tangent_at;
use catallax_numbers::{Decider, Real, RealPoly, Refusal, Sign};
use catallax_shapes::{Curve, CurveId, Point3, PointId, Surface, SurfaceId, Vector3};
use catallax_topology::{EdgeId, Model, Sense};
use std::cmp::Ordering;

#[derive(Clone)]
pub(crate) struct Germ {
    pub(crate) key: usize,
    pub(crate) tangent: Vector3,
    pub(crate) curvature: Vector3,
    source: Source,
    point: PointId,
}

#[derive(Clone, Copy)]
enum Source {
    Curve(CurveId),
    Surface(SurfaceId),
}

fn gradient(p: &RealPoly, point: &Point3) -> Vector3 {
    let x = [point.x.clone(), point.y.clone(), point.z.clone()];
    Vector3::new(
        p.partial(0).eval_real(&x),
        p.partial(1).eval_real(&x),
        p.partial(2).eval_real(&x),
    )
}
fn hessian(p: &RealPoly, point: &Point3, v: &Vector3) -> Vector3 {
    let x = [point.x.clone(), point.y.clone(), point.z.clone()];
    let t = [&v.x, &v.y, &v.z];
    let components: [Real; 3] = std::array::from_fn(|i| {
        (0..3).fold(Real::zero(), |sum, j| {
            sum + p.partial(i).partial(j).eval_real(&x) * t[j].clone()
        })
    });
    let [x, y, z] = components;
    Vector3::new(x, y, z)
}
fn solve(
    rows: [&Vector3; 3],
    rhs: [Real; 3],
    dec: &mut Decider,
    code: &'static str,
) -> Result<Vector3, Refusal> {
    let [a, b, c] = rows;
    let determinant = a.dot(&b.cross(c));
    if sign(&determinant, dec, code)? == Sign::Zero {
        return Err(refusal(
            code,
            "second-order equations do not determine an outgoing germ",
        ));
    }
    Ok(b.cross(c)
        .scaled(&rhs[0])
        .plus(&c.cross(a).scaled(&rhs[1]))
        .plus(&a.cross(b).scaled(&rhs[2]))
        .scaled(&(Real::one() / determinant)))
}

pub(crate) fn curvature(
    model: &Model,
    curve: CurveId,
    point: PointId,
    tangent: &Vector3,
    dec: &mut Decider,
    code: &'static str,
) -> Result<Vector3, Refusal> {
    let line = match model.geom().curve(curve) {
        Curve::Line(_) => true,
        Curve::Intersection(i) => {
            i.surfaces
                .iter()
                .all(|s| matches!(model.geom().surface(*s), Surface::Plane(_)))
                || matches!(
                    catallax_intersect::curve_form(model.geom(), curve, dec)?,
                    catallax_intersect::CurveForm::Line { .. }
                )
        }
        _ => false,
    };
    if line {
        // tangent_at has certified the direction; affine carriers have zero
        // second derivatives, so their normal acceleration is identically zero.
        return Ok(Vector3::new(Real::zero(), Real::zero(), Real::zero()));
    }
    let position = model.geom().position(point);
    let equations = model.geom().curve(curve).equations(model.geom());
    let mut first = None;
    for p in &equations {
        let g = gradient(p, &position);
        if sign(&g.norm_squared(), dec, code)? != Sign::Zero {
            first = Some((p, g));
            break;
        }
    }
    let Some((f, a)) = first else {
        return Err(refusal(code, "curve gradients vanish at the vertex"));
    };
    for g in &equations {
        let b = gradient(g, &position);
        if sign(&a.cross(&b).norm_squared(), dec, code)? != Sign::Zero {
            return solve(
                [&a, &b, tangent],
                [
                    -tangent.dot(&hessian(f, &position, tangent)),
                    -tangent.dot(&hessian(g, &position, tangent)),
                    Real::zero(),
                ],
                dec,
                code,
            );
        }
    }
    // At a node of two quadrics, the common normal fixes normal curvature.
    // The cubic term of g-lambda*f fixes its remaining transverse component.
    let aa = a.norm_squared();
    for g in &equations {
        let b = gradient(g, &position);
        let lambda = b.dot(&a) / aa.clone();
        let difference =
            hessian(g, &position, tangent).minus(&hessian(f, &position, tangent).scaled(&lambda));
        if sign(&a.cross(&difference).dot(tangent), dec, code)? != Sign::Zero {
            return solve(
                [&a, &difference, tangent],
                [
                    -tangent.dot(&hessian(f, &position, tangent)),
                    Real::zero(),
                    Real::zero(),
                ],
                dec,
                code,
            );
        }
    }
    Err(refusal(
        code,
        "second-order branch information does not separate tangent germs",
    ))
}

struct Jet {
    hessians: [[Vector3; 3]; 2],
    dual: [Vector3; 2],
    coefficients: Vec<Vector3>,
}
impl Jet {
    fn new(
        model: &Model,
        germ: &Germ,
        reference: &Vector3,
        normal: &Vector3,
        dec: &mut Decider,
        code: &'static str,
    ) -> Result<Self, Refusal> {
        let point = model.geom().position(germ.point);
        let (equations, gradients) = match germ.source {
            Source::Curve(curve) => {
                let equations = model.geom().curve(curve).equations(model.geom());
                let mut independent = None;
                'pair: for (i, a) in equations.iter().enumerate() {
                    let a = gradient(a, &point);
                    for b in &equations[i + 1..] {
                        let b_gradient = gradient(b, &point);
                        if sign(&a.cross(&b_gradient).norm_squared(), dec, code)? != Sign::Zero {
                            independent = Some((
                                [Some(equations[i].clone()), Some(b.clone())],
                                [a, b_gradient],
                            ));
                            break 'pair;
                        }
                    }
                }
                independent
                    .ok_or_else(|| refusal(code, "higher-order germs require a regular curve"))?
            }
            Source::Surface(surface) => {
                let equation = model.geom().surface(surface).implicit();
                let gradient = gradient(&equation, &point);
                ([Some(equation), None], [gradient, normal.clone()])
            }
        };
        let determinant = gradients[0].dot(&gradients[1].cross(reference));
        if sign(&determinant, dec, code)? == Sign::Zero {
            return Err(refusal(code, "germ parameter is not transverse"));
        }
        let inverse = Real::one() / determinant;
        let dual = [
            gradients[1].cross(reference).scaled(&inverse),
            reference.cross(&gradients[0]).scaled(&inverse),
        ];
        let position = [point.x, point.y, point.z];
        let hessians = equations.map(|equation| {
            std::array::from_fn(|i| {
                let row: [Real; 3] = std::array::from_fn(|j| {
                    equation.as_ref().map_or_else(Real::zero, |equation| {
                        equation.partial(i).partial(j).eval_real(&position)
                    })
                });
                let [x, y, z] = row;
                Vector3::new(x, y, z)
            })
        });
        let parameter = reference.dot(&germ.tangent);
        let mut coefficients = Vec::with_capacity(9);
        coefficients.push(Vector3::new(Real::zero(), Real::zero(), Real::zero()));
        coefficients.push(germ.tangent.scaled(&(&Real::one() / &parameter)));
        coefficients.push(
            germ.curvature
                .scaled(&(Real::one() / (Real::from_i64(2) * parameter.square()))),
        );
        Ok(Self {
            hessians,
            dual,
            coefficients,
        })
    }

    fn next(&mut self) -> Vector3 {
        let n = self.coefficients.len();
        let rhs: [Real; 2] = self.hessians.each_ref().map(|hessian| {
            let mut sum = Real::zero();
            for i in 1..n {
                let a = &self.coefficients[i];
                let b = &self.coefficients[n - i];
                let product = Vector3::new(hessian[0].dot(b), hessian[1].dot(b), hessian[2].dot(b));
                sum = sum + a.dot(&product);
            }
            -sum / Real::from_i64(2)
        });
        let next = self.dual[0]
            .scaled(&rhs[0])
            .plus(&self.dual[1].scaled(&rhs[1]));
        self.coefficients.push(next.clone());
        next
    }
}

fn higher_order(
    model: &Model,
    a: &Germ,
    b: &Germ,
    normal: &Vector3,
    dec: &mut Decider,
    code: &'static str,
) -> Result<Ordering, Refusal> {
    let mut left = Jet::new(model, a, &a.tangent, normal, dec, code)?;
    let mut right = Jet::new(model, b, &a.tangent, normal, dec, code)?;
    // Quadric sections of one quadric have contact order at most eight
    // unless they share a component. Both use t = tangent . (p-p0).
    for _ in 3..=8 {
        let difference = right.next().minus(&left.next());
        match sign(&normal.dot(&a.tangent.cross(&difference)), dec, code)? {
            Sign::Positive => return Ok(Ordering::Less),
            Sign::Negative => return Ok(Ordering::Greater),
            Sign::Zero => {}
        }
    }
    Err(refusal(code, "tangent germs agree through eighth order"))
}

pub(crate) fn outgoing(
    model: &Model,
    edge: EdgeId,
    sense: Sense,
    point: PointId,
    key: usize,
    dec: &mut Decider,
    code: &'static str,
) -> Result<Germ, Refusal> {
    let e = model.edge(edge);
    let travel = product(e.sense, sense);
    let ts = tangent_at(model.geom(), e.curve, point, dec)?;
    let tangent = if travel == Sense::Same {
        ts[1].clone()
    } else {
        oriented(&ts[0], Sense::Reversed)
    };
    let curvature = curvature(model, e.curve, point, &tangent, dec, code)?;
    Ok(Germ {
        key,
        tangent,
        curvature,
        source: Source::Curve(e.curve),
        point,
    })
}

fn half(
    reference: &Vector3,
    v: &Vector3,
    normal: &Vector3,
    dec: &mut Decider,
    code: &'static str,
) -> Result<bool, Refusal> {
    match sign(&normal.dot(&reference.cross(v)), dec, code)? {
        Sign::Positive => Ok(false),
        Sign::Negative => Ok(true),
        Sign::Zero => Ok(sign(&reference.dot(v), dec, code)? == Sign::Negative),
    }
}

pub(crate) fn compare(
    model: &Model,
    a: &Germ,
    b: &Germ,
    reference: &Vector3,
    normal: &Vector3,
    dec: &mut Decider,
    code: &'static str,
) -> Result<Ordering, Refusal> {
    if a.key == b.key {
        return Ok(Ordering::Equal);
    }
    let ah = half(reference, &a.tangent, normal, dec, code)?;
    let bh = half(reference, &b.tangent, normal, dec, code)?;
    if ah != bh {
        return Ok(ah.cmp(&bh));
    }
    match sign(&normal.dot(&a.tangent.cross(&b.tangent)), dec, code)? {
        Sign::Positive => Ok(Ordering::Less),
        Sign::Negative => Ok(Ordering::Greater),
        Sign::Zero => {
            if sign(&a.tangent.dot(&b.tangent), dec, code)? != Sign::Positive {
                return Err(refusal(code, "radial germ has no strict tangent direction"));
            }
            let difference = b
                .curvature
                .scaled(&a.tangent.norm_squared())
                .minus(&a.curvature.scaled(&b.tangent.norm_squared()));
            match sign(&normal.dot(&a.tangent.cross(&difference)), dec, code)? {
                Sign::Positive => Ok(Ordering::Less),
                Sign::Negative => Ok(Ordering::Greater),
                Sign::Zero => higher_order(model, a, b, normal, dec, code),
            }
        }
    }
}

pub(crate) fn sort(
    model: &Model,
    germs: &mut [Germ],
    normal: &Vector3,
    dec: &mut Decider,
    code: &'static str,
) -> Result<(), Refusal> {
    if germs.len() < 2 {
        return Ok(());
    }
    let reference = germs[0].tangent.clone();
    for i in 1..germs.len() {
        let mut j = i;
        while j > 0
            && compare(
                model,
                &germs[j],
                &germs[j - 1],
                &reference,
                normal,
                dec,
                code,
            )? == Ordering::Less
        {
            germs.swap(j, j - 1);
            j -= 1;
        }
    }
    Ok(())
}

pub(crate) fn surface_germ(
    model: &Model,
    surface: catallax_shapes::SurfaceId,
    point: PointId,
    tangent: Vector3,
    key: usize,
    dec: &mut Decider,
) -> Result<Germ, Refusal> {
    let s = model.geom().surface(surface);
    let p = model.geom().position(point);
    let normal = s.gradient(&p);
    let n2 = normal.norm_squared();
    if sign(&n2, dec, "boolean.merge-ambiguous")? == Sign::Zero {
        return Err(refusal(
            "boolean.merge-ambiguous",
            "surface normal vanished on a sewn edge",
        ));
    }
    let acceleration = -tangent.dot(&hessian(&s.implicit(), &p, &tangent)) / n2;
    let curvature = normal.scaled(&acceleration);
    Ok(Germ {
        key,
        tangent,
        curvature,
        source: Source::Surface(surface),
        point,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::Rational;
    use catallax_shapes::{Circle, Ellipse, PointDef, RPoint3, RVector3};

    fn q(value: &str) -> Rational {
        Rational::parse(value).unwrap()
    }

    fn germ(model: &Model, curve: CurveId, point: PointId, key: usize, speed: i64) -> Germ {
        let tangent = Vector3::new(Real::from_i64(speed), Real::zero(), Real::zero());
        let curvature = curvature(
            model,
            curve,
            point,
            &tangent,
            &mut Decider::default(),
            "boolean.imprint-ambiguous",
        )
        .unwrap();
        Germ {
            key,
            tangent,
            curvature,
            source: Source::Curve(curve),
            point,
        }
    }

    #[test]
    fn osculating_conics_are_ordered_by_third_and_fourth_coefficients() {
        let mut model = Model::new();
        let normal = RVector3::new(q("0"), q("0"), q("1"));
        let circle = model.geom_mut().add_curve(Curve::Circle(Circle {
            center: RPoint3::new(q("0"), q("1"), q("0")).into(),
            normal: normal.clone().into(),
            radius_sq: q("1").into(),
            ref_dir: RVector3::new(q("1"), q("0"), q("0")).into(),
        }));
        // x² + 3xy/2 + y² - 2y = 0 gives y = x²/2 + 3x³/8 + O(x⁴).
        let tilted = model.geom_mut().add_curve(Curve::Ellipse(Ellipse {
            center: RPoint3::new(q("-12/7"), q("16/7"), q("0")).into(),
            normal: normal.clone().into(),
            major_dir: RVector3::new(q("1"), q("-1"), q("0")).into(),
            major_sq: q("64/7").into(),
            minor_sq: q("64/49").into(),
        }));
        // The upright ellipse gives y = x²/2 + x⁴/16 + O(x⁶);
        // the circle gives y = x²/2 + x⁴/8 + O(x⁶).
        let upright = model.geom_mut().add_curve(Curve::Ellipse(Ellipse {
            center: RPoint3::new(q("0"), q("2"), q("0")).into(),
            normal: normal.clone().into(),
            major_dir: RVector3::new(q("0"), q("1"), q("0")).into(),
            major_sq: q("4").into(),
            minor_sq: q("2").into(),
        }));
        let point = model.geom_mut().add_point(PointDef::Explicit(
            (RPoint3::new(q("0"), q("0"), q("0"))).into(),
        ));
        let normal = normal.to_real();
        let a = germ(&model, circle, point, 0, 2);
        for (curve, want) in [(tilted, Ordering::Less), (upright, Ordering::Greater)] {
            let b = germ(&model, curve, point, 1, 3);
            let mut dec = Decider::default();
            assert_eq!(
                compare(
                    &model,
                    &a,
                    &b,
                    &a.tangent,
                    &normal,
                    &mut dec,
                    "boolean.imprint-ambiguous"
                )
                .unwrap(),
                want,
            );
            assert_eq!(
                compare(
                    &model,
                    &b,
                    &a,
                    &a.tangent,
                    &normal,
                    &mut dec,
                    "boolean.imprint-ambiguous"
                )
                .unwrap(),
                want.reverse(),
            );
        }
        let duplicate = germ(&model, circle, point, 1, 3);
        assert_eq!(
            compare(
                &model,
                &a,
                &duplicate,
                &a.tangent,
                &normal,
                &mut Decider::default(),
                "boolean.imprint-ambiguous",
            )
            .unwrap_err()
            .code,
            "boolean.imprint-ambiguous",
        );
    }
}
