//! Exact setback contacts and ruled carriers for elementary bevels.

use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{
    Circle, Cone, Curve, Cylinder, Isometry, Line, Plane, Point3, Surface, Vector3,
};

use crate::profile::sign;

/// Exact finite circular domains used by the independent chamfer builders.
#[path = "chamfer_arc.rs"]
pub mod arc;

/// The immutable generating definition of a plane or conical chamfer strip.
#[derive(Clone, Debug)]
pub struct ChamferDefinition {
    /// Persistent bevel origin, independent of the setbacks.
    pub origin: Origin,
    /// Input supports, in the same order as the setbacks and contacts.
    pub supports: [Surface; 2],
    /// Input material orientations relative to the support normals.
    pub same_sense: [bool; 2],
    /// Original sharp edge carrier.
    pub edge: Curve,
    /// Positive Euclidean distances from the edge on the corresponding supports.
    pub distances: [Real; 2],
    /// Convex strips remove material; concave strips add it.
    pub convex: bool,
    /// Exact ruled carrier joining corresponding setback points.
    pub surface: Surface,
    /// Material orientation relative to the ruled carrier's normal.
    pub surface_same_sense: bool,
    /// Contact curves on the corresponding supports.
    pub contacts: [Curve; 2],
}

pub(crate) fn bad(detail: &str) -> Refusal {
    Refusal::new("chamfer.unsupported-contact", detail)
}
pub(crate) fn root(x: &Real) -> Real {
    x.as_rational()
        .and_then(|r| crate::extrude::sqrt(&r))
        .map(Real::from)
        .unwrap_or_else(|| x.sqrt())
}
pub(crate) fn unit(v: &Vector3) -> Vector3 {
    v.scaled(&(&Real::one() / &root(&v.norm_squared())))
}
fn normal(v: &Vector3, same: bool) -> Vector3 {
    unit(v).scaled(&Real::from_i64(if same { 1 } else { -1 }))
}

impl ChamferDefinition {
    /// Constructs exact contacts and their analytic reduction. Trimming and
    /// whole-body validity must be established separately before publication.
    pub fn new(
        origin: Origin,
        supports: [Surface; 2],
        same_sense: [bool; 2],
        edge: Curve,
        distances: [Real; 2],
        convex: bool,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        for distance in &distances {
            if sign(distance, dec)? != Sign::Positive {
                return Err(Refusal::new(
                    "chamfer.non-positive-distance",
                    "chamfer distances must be positive",
                ));
            }
        }
        let k = Real::from_i64(if convex { 1 } else { -1 });
        let (surface, surface_same_sense, contacts) = match (&supports[0], &supports[1], &edge) {
            (Surface::Plane(a), Surface::Plane(b), Curve::Line(e)) => {
                if sign(&e.dir.norm_squared(), dec)? == Sign::Zero {
                    return Err(bad("zero edge direction"));
                }
                for plane in [a, b] {
                    if sign(&plane.normal.norm_squared(), dec)? == Sign::Zero
                        || sign(&plane.normal.dot(&e.dir), dec)? != Sign::Zero
                        || sign(&e.origin.minus(&plane.origin).dot(&plane.normal), dec)?
                            != Sign::Zero
                    {
                        return Err(bad("line must lie on both regular input planes"));
                    }
                }
                let na = normal(&a.normal, same_sense[0]);
                let nb = normal(&b.normal, same_sense[1]);
                let dot = na.dot(&nb);
                let sine_sq = Real::one() - dot.square();
                if sign(&sine_sq, dec)? != Sign::Positive {
                    return Err(bad("parallel planes have no sharp edge"));
                }
                let factor = -k / root(&sine_sq);
                let a = e.origin.plus(
                    &nb.minus(&na.scaled(&dot))
                        .scaled(&(&factor * &distances[0])),
                );
                let b = e.origin.plus(
                    &na.minus(&nb.scaled(&dot))
                        .scaled(&(&factor * &distances[1])),
                );
                let surface = Surface::Plane(Plane {
                    origin: a.clone(),
                    normal: na.scaled(&distances[0]).plus(&nb.scaled(&distances[1])),
                    ref_dir: e.dir.clone(),
                });
                (
                    surface,
                    true,
                    [a, b].map(|origin| {
                        Curve::Line(Line {
                            origin,
                            dir: e.dir.clone(),
                        })
                    }),
                )
            }
            (Surface::Plane(p), Surface::Cylinder(c), Curve::Circle(e)) => {
                circular(p, c, e, same_sense, &distances, &k, dec)?
            }
            (Surface::Cylinder(c), Surface::Plane(p), Curve::Circle(e)) => {
                let (surface, sense, [a, b]) = circular(
                    p,
                    c,
                    e,
                    [same_sense[1], same_sense[0]],
                    &[distances[1].clone(), distances[0].clone()],
                    &k,
                    dec,
                )?;
                (surface, sense, [b, a])
            }
            _ => {
                return Err(bad(
                    "elementary bevels require a plane-plane line or perpendicular plane-cylinder circle",
                ));
            }
        };
        Ok(Self {
            origin,
            supports,
            same_sense,
            edge,
            distances,
            convex,
            surface,
            surface_same_sense,
            contacts,
        })
    }

    /// Applies an exact isometry while retaining support order and material sense.
    pub fn transform(&mut self, map: &Isometry) -> Result<(), Refusal> {
        let curve = |c: &Curve| {
            map.curve(c)
                .ok_or_else(|| bad("contact has no exact isometry"))
        };
        self.edge = curve(&self.edge)?;
        self.contacts = [curve(&self.contacts[0])?, curve(&self.contacts[1])?];
        self.supports = self.supports.each_ref().map(|s| map.surface(s));
        self.surface = map.surface(&self.surface);
        Ok(())
    }

    /// Certifies that the complete affected region misses another support.
    /// For a finite line run, `vertices` must include both original endpoints
    /// and all four end contacts. Circle runs use their full annular section.
    /// Failure to separate is conservative and requires further trimming work.
    pub fn avoids_support(
        &self,
        surface: &Surface,
        vertices: &[Point3],
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        if let Curve::Line(_) = &self.edge {
            if vertices.len() < 6 {
                return Ok(false);
            }
            if let Surface::Plane(plane) = surface {
                let mut side = None;
                for p in vertices {
                    let s = sign(&p.minus(&plane.origin).dot(&plane.normal), dec)?;
                    if s == Sign::Zero || side.is_some_and(|old| old != s) {
                        return Ok(false);
                    }
                    side = Some(s);
                }
                return Ok(true);
            }
        }
        if let (Curve::Circle(edge), Surface::Cylinder(c)) = (&self.edge, surface)
            && sign(&edge.normal.cross(&c.axis).norm_squared(), dec)? == Sign::Zero
            && sign(
                &edge.center.minus(&c.origin).cross(&c.axis).norm_squared(),
                dec,
            )? == Sign::Zero
        {
            let [Curve::Circle(a), Curve::Circle(b)] = &self.contacts else {
                return Ok(false);
            };
            let radii = [&edge.radius_sq, &a.radius_sq, &b.radius_sq];
            let [a, b, c] = radii.map(|r| sign(&(r - &c.radius_sq), dec));
            let [a, b, c] = [a?, b?, c?];
            return Ok(a != Sign::Zero && a == b && b == c);
        }
        let Some(bounds) = self.region_box(vertices, dec)? else {
            return Ok(false);
        };
        let Some(polynomial) = surface.implicit() else {
            return Ok(false);
        };
        let mut value = Ball::from_rational(&Rational::zero(), 128);
        for (monomial, coefficient) in polynomial.terms() {
            let mut term = dec.enclose_budgeted_at(&coefficient, 128)?;
            for (exponent, coordinate) in monomial.exponents().iter().zip(&bounds) {
                if *exponent != 0 {
                    term = &term * &coordinate.pow_u32(*exponent);
                }
            }
            value = &value + &term;
        }
        Ok(value.lo().sign() == Sign::Positive || value.hi().sign() == Sign::Negative)
    }

    /// Conservative box of the entire swept affected region, including its
    /// original edge. Used only to prove separation, never to accept overlap.
    pub fn region_box(
        &self,
        vertices: &[Point3],
        dec: &mut Decider,
    ) -> Result<Option<[Ball; 3]>, Refusal> {
        let coordinates = |p: &Point3| [p.x.clone(), p.y.clone(), p.z.clone()];
        let mut points = vertices.to_vec();
        let mut extent = [Real::zero(), Real::zero(), Real::zero()];
        if let Curve::Circle(edge) = &self.edge {
            let [Curve::Circle(a), Curve::Circle(b)] = &self.contacts else {
                return Ok(None);
            };
            points.extend([edge.center.clone(), a.center.clone(), b.center.clone()]);
            let mut radius_sq = edge.radius_sq.clone();
            for r in [&a.radius_sq, &b.radius_sq] {
                if sign(&(r - &radius_sq), dec)? == Sign::Positive {
                    radius_sq = r.clone();
                }
            }
            let n = [
                edge.normal.x.clone(),
                edge.normal.y.clone(),
                edge.normal.z.clone(),
            ];
            extent = n.map(|v| {
                root(&(&radius_sq * &(Real::one() - v.square() / edge.normal.norm_squared())))
            });
        }
        let Some(first) = points.first() else {
            return Ok(None);
        };
        let mut bounds = coordinates(first)
            .map(|v| dec.enclose_budgeted_at(&v, 128))
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?;
        for point in &points[1..] {
            for (i, value) in coordinates(point).iter().enumerate() {
                bounds[i] = bounds[i].hull(&dec.enclose_budgeted_at(value, 128)?);
            }
        }
        for i in 0..3 {
            let e = dec.enclose_budgeted_at(&extent[i], 128)?;
            bounds[i] = (&bounds[i] - &e).hull(&(&bounds[i] + &e));
        }
        Ok(Some(bounds.try_into().expect("three coordinates")))
    }

    /// Proves disjoint affected regions from finite hulls or coaxial annuli.
    /// Finite line vertices have the same coverage requirement as `region_box`.
    pub fn avoids_region(
        &self,
        other: &Self,
        vertices: &[Point3],
        other_vertices: &[Point3],
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        if let (Some(a), Some(b)) = (
            self.region_box(vertices, dec)?,
            other.region_box(other_vertices, dec)?,
        ) && (0..3).any(|i| a[i].hi() < b[i].lo() || b[i].hi() < a[i].lo())
        {
            return Ok(true);
        }
        let (Curve::Circle(a), Curve::Circle(b)) = (&self.edge, &other.edge) else {
            return Ok(false);
        };
        if sign(&a.normal.cross(&b.normal).norm_squared(), dec)? != Sign::Zero
            || sign(
                &a.center.minus(&b.center).cross(&a.normal).norm_squared(),
                dec,
            )? != Sign::Zero
        {
            return Ok(false);
        }
        let radial_range =
            |definition: &Self, dec: &mut Decider| -> Result<(Real, Real), Refusal> {
                let Curve::Circle(edge) = &definition.edge else {
                    return Err(bad("expected circular region"));
                };
                let (mut lo, mut hi) = (edge.radius_sq.clone(), edge.radius_sq.clone());
                for contact in &definition.contacts {
                    let Curve::Circle(c) = contact else {
                        return Err(bad("expected circular contact"));
                    };
                    if sign(&(&c.radius_sq - &lo), dec)? == Sign::Negative {
                        lo = c.radius_sq.clone();
                    }
                    if sign(&(&c.radius_sq - &hi), dec)? == Sign::Positive {
                        hi = c.radius_sq.clone();
                    }
                }
                Ok((lo, hi))
            };
        let (al, ah) = radial_range(self, dec)?;
        let (bl, bh) = radial_range(other, dec)?;
        Ok(
            sign(&(&bl - &ah), dec)? == Sign::Positive
                || sign(&(&al - &bh), dec)? == Sign::Positive,
        )
    }
}

fn circular(
    plane: &Plane,
    cylinder: &Cylinder,
    edge: &Circle,
    senses: [bool; 2],
    distances: &[Real; 2],
    k: &Real,
    dec: &mut Decider,
) -> Result<(Surface, bool, [Curve; 2]), Refusal> {
    if sign(&plane.normal.norm_squared(), dec)? == Sign::Zero
        || sign(&cylinder.axis.norm_squared(), dec)? == Sign::Zero
        || sign(&cylinder.radius_sq, dec)? != Sign::Positive
        || sign(&edge.normal.norm_squared(), dec)? == Sign::Zero
        || sign(&edge.ref_dir.norm_squared(), dec)? == Sign::Zero
        || sign(&plane.normal.cross(&cylinder.axis).norm_squared(), dec)? != Sign::Zero
        || sign(&edge.normal.cross(&cylinder.axis).norm_squared(), dec)? != Sign::Zero
        || sign(&edge.ref_dir.dot(&edge.normal), dec)? != Sign::Zero
        || sign(
            &edge
                .center
                .minus(&cylinder.origin)
                .cross(&cylinder.axis)
                .norm_squared(),
            dec,
        )? != Sign::Zero
        || sign(&edge.center.minus(&plane.origin).dot(&plane.normal), dec)? != Sign::Zero
        || sign(&(&edge.radius_sq - &cylinder.radius_sq), dec)? != Sign::Zero
    {
        return Err(bad(
            "circle must lie on a regular perpendicular plane and cylinder",
        ));
    }
    let n = normal(&plane.normal, senses[0]);
    let s = Real::from_i64(if senses[1] { 1 } else { -1 });
    let radius = root(&cylinder.radius_sq);
    let contact_radius = &radius - &(&(k * &s) * &distances[0]);
    if sign(&contact_radius, dec)? != Sign::Positive {
        return Err(Refusal::new(
            "chamfer.unsupported-overflow",
            "a radial contact reaches the axis and needs continuation",
        ));
    }
    let apex_height = &(&(&s * &distances[1]) * &contact_radius) / &distances[0];
    let surface = Surface::Cone(Cone {
        apex: edge.center.plus(&n.scaled(&apex_height)),
        axis: n.scaled(&(-s)),
        slope: &distances[0] / &distances[1],
        ref_dir: edge.ref_dir.clone(),
    });
    let circle = |center, radius_sq| {
        Curve::Circle(Circle {
            center,
            radius_sq,
            normal: edge.normal.clone(),
            ref_dir: edge.ref_dir.clone(),
        })
    };
    Ok((
        surface,
        senses[1],
        [
            circle(edge.center.clone(), contact_radius.square()),
            circle(
                edge.center.plus(&n.scaled(&(-(k * &distances[1])))),
                cylinder.radius_sq.clone(),
            ),
        ],
    ))
}
