//! Rational polyhedral integration, without Ball work or Decider calls.

use catallax_numbers::{Integer, Rational};
use catallax_shapes::{Curve, PointDef, RPoint3, RVector3, Surface};
use catallax_topology::{FaceId, Loop, Model, Sense};

pub(crate) struct Planar {
    pub(crate) volume: Rational,
    pub(crate) moments: [Rational; 3],
    pub(crate) area_numerator: Rational,
    pub(crate) normal_squared: Rational,
}

pub(crate) fn components(p: &RPoint3) -> [&Rational; 3] {
    [&p.x, &p.y, &p.z]
}

pub(crate) fn direction(v: &RVector3) -> [&Rational; 3] {
    [&v.x, &v.y, &v.z]
}

pub(crate) fn real_components(p: &catallax_shapes::Point3) -> [&catallax_numbers::Real; 3] {
    [&p.x, &p.y, &p.z]
}
pub(crate) fn real_direction(v: &catallax_shapes::Vector3) -> [&catallax_numbers::Real; 3] {
    [&v.x, &v.y, &v.z]
}

pub(crate) fn sign(sense: Sense) -> Rational {
    Rational::from_i64(if sense == Sense::Same { 1 } else { -1 })
}

pub(crate) fn face(
    model: &Model,
    id: FaceId,
    anchor: &RPoint3,
    property: crate::surface::Property,
) -> Option<Planar> {
    let face = model.face(id);
    let Surface::Plane(plane) = model.geom().surface(face.surface) else {
        return None;
    };
    let normal = plane.normal.as_rational()?;
    let origin = plane.origin.as_rational()?.minus(anchor);
    let mut result = Planar::new(if property == crate::surface::Property::Area {
        normal.norm_squared()
    } else {
        Rational::one()
    });
    for boundary in &face.loops {
        let Loop::Coedges(coedges) = model.loop_(*boundary) else {
            return None;
        };
        for id in coedges {
            let coedge = model.coedge(*id);
            let edge = model.edge(coedge.edge);
            if !matches!(model.geom().curve(edge.curve), Curve::Line(_)) {
                return None;
            }
            let (a, b) = (edge.start?, edge.end?);
            let PointDef::Explicit(a) = model.geom().point(model.vertex(a).point) else {
                return None;
            };
            let PointDef::Explicit(b) = model.geom().point(model.vertex(b).point) else {
                return None;
            };
            let (a, b) = if coedge.sense == Sense::Same {
                (a, b)
            } else {
                (b, a)
            };
            let a = a.as_rational()?.minus(anchor);
            let b = b.as_rational()?.minus(anchor);
            result.segment(&origin, &a, &b, &normal, property);
        }
    }
    result.area_numerator = result.area_numerator * sign(face.sense) / Rational::from_i64(2);
    Some(result)
}

pub(crate) fn sqrt(r: &Rational) -> Option<Rational> {
    fn integer_sqrt(n: &Integer) -> Option<Integer> {
        if n.is_zero() {
            return Some(Integer::zero());
        }
        let mut x = Integer::one().mul_pow2(n.bit_length().div_ceil(2));
        loop {
            let quotient = n.div_rem_floor(&x)?.0;
            let next = (&x + &quotient).div_rem_floor(&Integer::from_i64(2))?.0;
            if next >= x {
                return (&x * &x == *n).then_some(x);
            }
            x = next;
        }
    }
    Rational::new(integer_sqrt(r.numerator())?, integer_sqrt(r.denominator())?)
}

pub(crate) fn line(
    surface: &Surface,
    geom: &catallax_shapes::GeomStore,
    from: Option<catallax_shapes::PointId>,
    to: Option<catallax_shapes::PointId>,
    anchor: &RPoint3,
    property: crate::surface::Property,
) -> Option<[Rational; 5]> {
    let PointDef::Explicit(from) = geom.point(from?) else {
        return None;
    };
    let PointDef::Explicit(to) = geom.point(to?) else {
        return None;
    };
    if let Surface::Cone(cone) = surface
        && [from, to].into_iter().any(|p| {
            [&p.x, &p.y, &p.z]
                .into_iter()
                .zip([&cone.apex.x, &cone.apex.y, &cone.apex.z])
                .all(|(a, b)| a.same_definition(b))
        })
    {
        // Along a straight generator, axis . (q x dq) is identically zero.
        // Every cone boundary form therefore has zero integral. Evaluating its
        // radial-unit-vector factors at the apex would introduce a removable 0/0.
        return Some(std::array::from_fn(|_| Rational::zero()));
    }
    if matches!(
        property,
        crate::surface::Property::SecondMoment(_, _) | crate::surface::Property::FaceCentroid
    ) {
        return None;
    }
    let Surface::Plane(plane) = surface else {
        return None;
    };
    let origin = plane.origin.as_rational()?;
    let normal = plane.normal.as_rational()?;
    let reference = plane.ref_dir.as_rational()?;
    let from = from.as_rational()?;
    let to = to.as_rational()?;
    let reference_squared = reference.norm_squared();
    let second = normal.cross(&reference);
    let a = from.minus(&origin);
    let d = to.minus(&from);
    let u0 = a.dot(&reference) / reference_squared.clone();
    let du = d.dot(&reference) / reference_squared.clone();
    let v0 = a.dot(&second) / second.norm_squared();
    let dv = d.dot(&second) / second.norm_squared();
    let vv = [v0, dv, Rational::zero(), Rational::zero()];
    let origin = origin.minus(anchor);
    let average_v = &vv[0] + &(&vv[1] / &Rational::from_i64(2));
    let area = if property == crate::surface::Property::Area {
        -sqrt(&normal.norm_squared())? * reference_squared.clone() * du.clone() * average_v.clone()
    } else {
        Rational::zero()
    };
    let volume = if property != crate::surface::Property::Area {
        -reference_squared.clone() * normal.dot(&origin) * du.clone() * average_v
            / Rational::from_i64(3)
    } else {
        Rational::zero()
    };
    let moments = if property == crate::surface::Property::Centroid {
        let first = direction(&reference);
        let second = direction(&second);
        let normal = direction(&normal);
        let origin = direction(&origin);
        let moments: [Rational; 3] = std::array::from_fn(|i| {
            let base = [
                origin[i] + &(&u0 * first[i]),
                &du * first[i],
                Rational::zero(),
                Rational::zero(),
            ];
            let term1 = multiply(&multiply(&base, &base), &vv);
            let term2 = multiply(&base, &multiply(&vv, &vv));
            let term3 = multiply(&multiply(&vv, &vv), &vv);
            let integral = (0..4).fold(Rational::zero(), |sum, k| {
                let coefficient = &term1[k]
                    + &(&term2[k] * second[i])
                    + &term3[k] * &(second[i] * second[i] / Rational::from_i64(3));
                sum + coefficient / Rational::from_i64((k + 1) as i64)
            });
            -reference_squared.clone() * normal[i].clone() * du.clone() * integral
                / Rational::from_i64(2)
        });
        moments
    } else {
        std::array::from_fn(|_| Rational::zero())
    };
    Some([
        area,
        volume,
        moments[0].clone(),
        moments[1].clone(),
        moments[2].clone(),
    ])
}

fn multiply(a: &[Rational; 4], b: &[Rational; 4]) -> [Rational; 4] {
    std::array::from_fn(|k| (0..=k).fold(Rational::zero(), |sum, j| sum + &a[j] * &b[k - j]))
}

impl Planar {
    fn new(normal_squared: Rational) -> Self {
        Self {
            volume: Rational::zero(),
            moments: std::array::from_fn(|_| Rational::zero()),
            area_numerator: Rational::zero(),
            normal_squared,
        }
    }

    fn segment(
        &mut self,
        origin: &RVector3,
        a: &RVector3,
        b: &RVector3,
        normal: &RVector3,
        property: crate::surface::Property,
    ) {
        let cross = a.cross(b);
        if property == crate::surface::Property::FaceCentroid {
            let area = normal.dot(&a.minus(origin).cross(&b.minus(origin)));
            let coordinates = [direction(origin), direction(a), direction(b)];
            for (i, moment) in self.moments.iter_mut().enumerate() {
                *moment = &*moment
                    + &(&area * &(&(coordinates[0][i] + coordinates[1][i]) + coordinates[2][i])
                        / Rational::from_i64(6));
            }
            self.area_numerator = &self.area_numerator + &normal.dot(&cross);
            return;
        }
        if property != crate::surface::Property::Area {
            self.volume = &self.volume + &(origin.dot(&cross) / Rational::from_i64(6));
        }
        if property == crate::surface::Property::Centroid {
            // Use the same field as curved faces, F_i = p_i^2 e_i / 2.
            let area_vector = a.minus(origin).cross(&b.minus(origin));
            let o = direction(origin);
            let aa = direction(a);
            let bb = direction(b);
            let area = direction(&area_vector);
            for i in 0..3 {
                let square_integral = o[i] * o[i]
                    + aa[i] * aa[i]
                    + bb[i] * bb[i]
                    + o[i] * aa[i]
                    + o[i] * bb[i]
                    + aa[i] * bb[i];
                self.moments[i] =
                    &self.moments[i] + &(area[i] * &square_integral / Rational::from_i64(24));
            }
        }
        if let crate::surface::Property::SecondMoment(i, j) = property {
            let area_vector = a.minus(origin).cross(&b.minus(origin));
            let o = direction(origin);
            let aa = direction(a);
            let bb = direction(b);
            let integral = (&(o[i] + aa[i]) + bb[i]) * (&(o[j] + aa[j]) + bb[j])
                + o[i] * o[j]
                + aa[i] * aa[j]
                + bb[i] * bb[j];
            self.moments[0] =
                &self.moments[0] + &(origin.dot(&area_vector) * integral / Rational::from_i64(120));
        }
        if property == crate::surface::Property::Area {
            self.area_numerator = &self.area_numerator + &normal.dot(&cross);
        }
    }
}

struct PlaneGroup<'a> {
    plane: &'a catallax_shapes::Plane,
    sense: Sense,
    edges: std::collections::BTreeMap<catallax_topology::EdgeId, i64>,
}

struct EdgeUse<'a> {
    first: &'a catallax_shapes::Plane,
    second: Option<&'a catallax_shapes::Plane>,
    count: u8,
}

pub(crate) fn body(
    model: &Model,
    faces: &[FaceId],
    anchor: &RPoint3,
    property: crate::surface::Property,
) -> Option<[Rational; 5]> {
    if property == crate::surface::Property::FaceCentroid {
        return None;
    }
    use std::collections::BTreeMap;
    let mut groups = BTreeMap::new();
    let mut uses = BTreeMap::new();
    for id in faces {
        let face = model.face(*id);
        let Surface::Plane(plane) = model.geom().surface(face.surface) else {
            return None;
        };
        let normal = plane.normal.as_rational()?;
        let origin = plane.origin.as_rational()?;
        let coefficient = direction(&normal)
            .into_iter()
            .find(|r| !r.is_zero())
            .unwrap();
        let canonical = normal.scaled(&coefficient.recip()?);
        let offset = &canonical.x * &origin.x + &canonical.y * &origin.y + &canonical.z * &origin.z;
        let positive =
            (coefficient.sign() == catallax_numbers::Sign::Positive) == (face.sense == Sense::Same);
        let key = ([canonical.x, canonical.y, canonical.z, offset], positive);
        let group = groups.entry(key).or_insert_with(|| PlaneGroup {
            plane,
            sense: face.sense,
            edges: BTreeMap::new(),
        });
        for boundary in &face.loops {
            match model.loop_(*boundary) {
                Loop::Vertex(vertex) => {
                    if !matches!(
                        model.geom().point(model.vertex(*vertex).point),
                        PointDef::Explicit(_)
                    ) {
                        return None;
                    }
                }
                Loop::Coedges(coedges) => {
                    for id in coedges {
                        let coedge = model.coedge(*id);
                        let edge = model.edge(coedge.edge);
                        for vertex in [edge.start, edge.end].into_iter().flatten() {
                            if !matches!(
                                model.geom().point(model.vertex(vertex).point),
                                PointDef::Explicit(_)
                            ) {
                                return None;
                            }
                        }
                        *group.edges.entry(coedge.edge).or_insert(0) +=
                            if coedge.sense == Sense::Same { 1 } else { -1 };
                        let usage = uses.entry(coedge.edge).or_insert(EdgeUse {
                            first: plane,
                            second: None,
                            count: 0,
                        });
                        usage.count += 1;
                        if usage.count == 2 {
                            usage.second = Some(plane);
                        }
                        if usage.count > 2 {
                            return None;
                        }
                    }
                }
            }
        }
    }
    // This reduction applies to a whole closed boundary, not a face query.
    if uses.values().any(|usage| usage.count != 2) {
        return None;
    }
    let mut sums: [Rational; 5] = std::array::from_fn(|_| Rational::zero());
    for group in groups.into_values() {
        let normal = group.plane.normal.as_rational()?;
        let mut integral = Planar::new(if property == crate::surface::Property::Area {
            normal.norm_squared()
        } else {
            Rational::one()
        });
        let origin = group.plane.origin.as_rational()?.minus(anchor);
        for (id, count) in group.edges {
            if count == 0 {
                continue;
            }
            let edge = model.edge(id);
            if !matches!(model.geom().curve(edge.curve), Curve::Line(_)) {
                let usage = &uses[&id];
                // An uncancelled edge on two distinct planes is a straight
                // segment, even if its stored carrier is an intersection.
                if usage
                    .first
                    .normal
                    .cross(&usage.second?.normal)
                    .norm_squared()
                    .as_rational()?
                    .is_zero()
                {
                    return None;
                }
            }
            let PointDef::Explicit(a) = model.geom().point(model.vertex(edge.start?).point) else {
                return None;
            };
            let PointDef::Explicit(b) = model.geom().point(model.vertex(edge.end?).point) else {
                return None;
            };
            let (a, b) = if count > 0 { (a, b) } else { (b, a) };
            integral.segment(
                &origin,
                &a.as_rational()?.minus(anchor),
                &b.as_rational()?.minus(anchor),
                &normal,
                property,
            );
        }
        if property == crate::surface::Property::Area {
            sums[0] = &sums[0]
                + &(integral.area_numerator * sign(group.sense)
                    / (Rational::from_i64(2) * sqrt(&integral.normal_squared)?));
        } else {
            sums[1] = &sums[1] + &integral.volume;
            for (sum, moment) in sums[2..].iter_mut().zip(integral.moments) {
                *sum = &*sum + &moment;
            }
        }
    }
    Some(sums)
}
