//! Exact source coordinates and native carriers, without topology admission.
use super::*;
use catallax_numbers::{Outcome, Rational, Real, Refusal, StrictSign, ZeroProof};
use catallax_shapes::{self as shapes, Point3, Vector3};
use shapes::spline_tools::{RationalSplineCurve, RationalSplineSurface};

/// A source declaration error or a native geometry construction refusal.
#[derive(Clone, Debug)]
pub enum MappingError {
    /// Invalid or missing STEP declaration.
    Source(ImportRefusal),
    /// Invalid geometry or a carrier requiring finite numerical mapping.
    Geometry(Refusal),
}
impl From<ImportRefusal> for MappingError {
    fn from(value: ImportRefusal) -> Self {
        Self::Source(value)
    }
}
impl From<Refusal> for MappingError {
    fn from(value: Refusal) -> Self {
        Self::Geometry(value)
    }
}
fn invalid(detail: &str) -> MappingError {
    Refusal::new("import.geometry-mapping", detail).into()
}

/// Normalized right-handed source placement, in the caller's Cartesian units.
#[derive(Clone, Debug)]
pub struct Frame {
    /// Placement origin.
    pub origin: Point3,
    /// Orthonormal x, y and z directions.
    pub axes: [Vector3; 3],
}

fn unit(v: Vector3, dec: &mut Decider) -> Result<Vector3, MappingError> {
    let positive = match dec.sign(&v.norm_squared()) {
        Outcome::Proven(StrictSign::Positive) => true,
        Outcome::Proven(_) | Outcome::ProvenCoincident(_) => false,
        Outcome::NotExact(n) => return Err(Refusal::from(n).into()),
        Outcome::Refused(r) => return Err(r.into()),
    };
    if !positive {
        return Err(invalid("zero placement direction"));
    }
    Ok(v.scaled(&(&Real::one() / &v.norm_squared().sqrt())))
}
fn vector(v: [Rational; 3]) -> Vector3 {
    let [x, y, z] = v.map(Real::from_rational);
    Vector3::new(x, y, z)
}

impl ParsedStep {
    /// Cartesian coordinates with explicit per-coordinate scales. For a pcurve,
    /// these may mix angle and length units; no global length factor is assumed.
    pub fn cartesian_point(
        &self,
        id: EntityId,
        scales: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<[Rational; 3], MappingError> {
        let Geometry::Point(values) = self.geometry(id, dec)? else {
            return Err(invalid("expected CARTESIAN_POINT"));
        };
        Ok(std::array::from_fn(|i| {
            values
                .get(i)
                .map_or_else(Rational::zero, |v| v * &scales[i])
        }))
    }

    /// Resolve and normalize a two- or three-dimensional DIRECTION.
    pub fn unit_direction(&self, id: EntityId, dec: &mut Decider) -> Result<Vector3, MappingError> {
        let Geometry::Direction(values) = self.geometry(id, dec)? else {
            return Err(invalid("expected DIRECTION"));
        };
        unit(
            vector(std::array::from_fn(|i| {
                values.get(i).cloned().unwrap_or_else(Rational::zero)
            })),
            dec,
        )
    }

    /// Resolve an Axis2 placement using OCCT's projected reference direction.
    /// A missing reference uses its Cartesian x direction unless parallel to z.
    pub fn axis_placement(
        &self,
        id: EntityId,
        length_mm: &Rational,
        dec: &mut Decider,
    ) -> Result<Frame, MappingError> {
        let (origin, z, x) = match self.geometry(id, dec)? {
            Geometry::Axis2Placement3d(p, z, x) => (p, z, x),
            Geometry::Axis2Placement2d(p, x) => (p, None, x),
            _ => return Err(invalid("expected AXIS2_PLACEMENT")),
        };
        let p = self.cartesian_point(origin, &std::array::from_fn(|_| length_mm.clone()), dec)?;
        let [px, py, pz] = p.map(Real::from_rational);
        let z = match z {
            Some(id) => self.unit_direction(id, dec)?,
            None => vector([Rational::zero(), Rational::zero(), Rational::one()]),
        };
        let missing_reference = x.is_none();
        let mut x = match x {
            Some(id) => self.unit_direction(id, dec)?,
            None => vector([Rational::one(), Rational::zero(), Rational::zero()]),
        };
        if missing_reference
            && x.cross(&z)
                .norm_squared()
                .rational_ref()
                .is_some_and(Rational::is_zero)
            && x.dot(&z).rational_ref().is_some()
        {
            x = vector([Rational::zero(), Rational::one(), Rational::zero()]);
        }
        let x = unit(x.minus(&z.scaled(&x.dot(&z))), dec)?;
        let y = z.cross(&x);
        Ok(Frame {
            origin: Point3::new(px, py, pz),
            axes: [x, y, z],
        })
    }

    /// Exact spline control data, including nonclamped active domains.
    /// Coordinate scales are explicit so the same reader supports 2D pcurves.
    pub fn spline_curve(
        &self,
        id: EntityId,
        scales: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<RationalSplineCurve, MappingError> {
        let Geometry::Curve(curve) = self.geometry(id, dec)? else {
            return Err(invalid("expected a spline curve"));
        };
        let (degree, knots, ids, weights) = match curve {
            Curve::BSpline(s) => (
                s.degree,
                expanded_knots(&s.knots, s.degree, s.controls.len(), dec)?,
                s.controls,
                s.weights,
            ),
            Curve::Polyline(ids) => {
                let n = ids.len();
                if n < 2 {
                    return Err(invalid("polyline requires at least two points"));
                }
                let knots = (0..n + 2)
                    .map(|i| Rational::from_integer((i.saturating_sub(1).min(n - 1) as u64).into()))
                    .collect();
                (1, knots, ids, None)
            }
            _ => return Err(invalid("expected a spline curve or polyline")),
        };
        let controls = ids
            .into_iter()
            .map(|id| self.cartesian_point(id, scales, dec))
            .collect::<Result<Vec<_>, _>>()?;
        let weights = weights.unwrap_or_else(|| vec![Rational::one(); controls.len()]);
        Ok(RationalSplineCurve::new(
            degree, knots, controls, weights, dec,
        )?)
    }

    /// Exact surface control data. Regularity and finite trim admission are separate.
    pub fn spline_surface(
        &self,
        id: EntityId,
        length_mm: &Rational,
        dec: &mut Decider,
    ) -> Result<RationalSplineSurface, MappingError> {
        let Geometry::Surface(Surface::BSpline(s)) = self.geometry(id, dec)? else {
            return Err(invalid("expected a spline surface"));
        };
        let count = [s.controls.len(), s.controls.first().map_or(0, Vec::len)];
        if s.controls.iter().any(|row| row.len() != count[1]) {
            return Err(invalid("ragged spline control grid"));
        }
        let knots = [
            expanded_knots(&s.knots[0], s.degree[0], count[0], dec)?,
            expanded_knots(&s.knots[1], s.degree[1], count[1], dec)?,
        ];
        let scales = std::array::from_fn(|_| length_mm.clone());
        let controls = s
            .controls
            .into_iter()
            .flatten()
            .map(|id| self.cartesian_point(id, &scales, dec))
            .collect::<Result<Vec<_>, _>>()?;
        let weights = match s.weights {
            None => vec![Rational::one(); controls.len()],
            Some(rows) => {
                if rows.len() != count[0] || rows.iter().any(|r| r.len() != count[1]) {
                    return Err(invalid("spline weight grid differs from controls"));
                }
                rows.into_iter().flatten().collect()
            }
        };
        Ok(RationalSplineSurface::new(
            s.degree, knots, count, controls, weights, dec,
        )?)
    }

    /// Exact finite quadratic carrier for a parabola in its original STEP
    /// parameter. Algebraic placement coordinates remain exact Real values.
    pub fn parabola_segment(
        &self,
        id: EntityId,
        interval: [Rational; 2],
        length_mm: &Rational,
        dec: &mut Decider,
    ) -> Result<shapes::BSplineCurve, MappingError> {
        let Geometry::Curve(Curve::Parabola(placement, focal)) = self.geometry(id, dec)? else {
            return Err(invalid("expected PARABOLA"));
        };
        if focal <= Rational::zero() || interval[0] >= interval[1] {
            return Err(invalid(
                "parabola requires positive focal distance and increasing parameters",
            ));
        }
        let frame = self.axis_placement(placement, length_mm, dec)?;
        let [a, b] = &interval;
        let divisor = &Rational::from_i64(4) * &focal;
        let xs = [
            (a * a) / divisor.clone(),
            (a * b) / divisor.clone(),
            (b * b) / divisor,
        ];
        let ys = [a.clone(), (a + b) / Rational::from_i64(2), b.clone()];
        let controls = xs
            .into_iter()
            .zip(ys)
            .map(|(x, y)| {
                frame
                    .origin
                    .plus(&frame.axes[0].scaled(&Real::from_rational(&x * length_mm)))
                    .plus(&frame.axes[1].scaled(&Real::from_rational(&y * length_mm)))
            })
            .collect();
        Ok(shapes::BSplineCurve::new(
            2,
            vec![
                a.clone(),
                a.clone(),
                a.clone(),
                b.clone(),
                b.clone(),
                b.clone(),
            ],
            controls,
            vec![Real::one(); 3],
            dec,
        )?)
    }

    /// Map an unwrapped curve to an exact native carrier. Trims and replicas
    /// require their own parameter/placement handling and are never discarded.
    pub fn native_curve(
        &self,
        id: EntityId,
        length_mm: &Rational,
        dec: &mut Decider,
    ) -> Result<shapes::Curve, MappingError> {
        let Geometry::Curve(curve) = self.geometry(id, dec)? else {
            return Err(invalid("expected CURVE"));
        };
        Ok(match curve {
            Curve::Line(p, v) => {
                let Geometry::Vector(direction, magnitude) = self.geometry(v, dec)? else {
                    return Err(invalid("expected VECTOR"));
                };
                if magnitude <= Rational::zero() {
                    return Err(invalid("line VECTOR magnitude must be positive"));
                }
                let p =
                    self.cartesian_point(p, &std::array::from_fn(|_| length_mm.clone()), dec)?;
                let [x, y, z] = p.map(Real::from_rational);
                let dir = self
                    .unit_direction(direction, dec)?
                    .scaled(&Real::from_rational(&magnitude * length_mm));
                shapes::Curve::Line(shapes::Line::new(Point3::new(x, y, z), dir, dec)?)
            }
            Curve::Circle(placement, r) => {
                if r <= Rational::zero() {
                    return Err(invalid("circle radius must be positive"));
                }
                let f = self.axis_placement(placement, length_mm, dec)?;
                shapes::Curve::Circle(shapes::Circle::new(
                    f.origin,
                    f.axes[2].clone(),
                    Real::from_rational(&r * length_mm).square(),
                    f.axes[0].clone(),
                    dec,
                )?)
            }
            Curve::Ellipse(placement, a, b) => {
                if a <= Rational::zero() || b <= Rational::zero() {
                    return Err(invalid("ellipse axes must be positive"));
                }
                let f = self.axis_placement(placement, length_mm, dec)?;
                if a == b {
                    shapes::Curve::Circle(shapes::Circle::new(
                        f.origin,
                        f.axes[2].clone(),
                        Real::from_rational(&a * length_mm).square(),
                        f.axes[0].clone(),
                        dec,
                    )?)
                } else {
                    let (major, minor, direction) = if a > b {
                        (a, b, f.axes[0].clone())
                    } else {
                        (b, a, f.axes[1].scaled(&Real::from_i64(-1)))
                    };
                    shapes::Curve::Ellipse(shapes::Ellipse::new(
                        f.origin,
                        f.axes[2].clone(),
                        direction,
                        Real::from_rational(&major * length_mm).square(),
                        Real::from_rational(&minor * length_mm).square(),
                        dec,
                    )?)
                }
            }
            Curve::Offset3d(parent, distance, _, reference) => {
                let raw = self.geometry(parent, dec)?;
                if !matches!(
                    raw,
                    Geometry::Curve(Curve::Line(_, _) | Curve::Circle(_, _))
                ) {
                    return Err(invalid("offset curve requires finite numerical mapping"));
                }
                let reference = self.unit_direction(reference, dec)?;
                let distance = Real::from_rational(&distance * length_mm);
                match self.native_curve(parent, length_mm, dec)? {
                    shapes::Curve::Line(line) => {
                        let normal = unit(line.dir.cross(&reference), dec)?;
                        shapes::Curve::Line(shapes::Line::new(
                            line.origin.plus(&normal.scaled(&(-distance))),
                            line.dir,
                            dec,
                        )?)
                    }
                    shapes::Curve::Circle(circle) => {
                        match dec.sign(&circle.normal.cross(&reference).norm_squared()) {
                            Outcome::ProvenCoincident(p)
                                if p.proof != ZeroProof::SeparationBound => {}
                            _ => {
                                return Err(invalid(
                                    "circle offset reference is not proved parallel to its normal",
                                ));
                            }
                        }
                        let direction = match dec.sign(&circle.normal.dot(&reference)) {
                            Outcome::Proven(StrictSign::Positive) => Real::one(),
                            Outcome::Proven(StrictSign::Negative) => Real::from_i64(-1),
                            _ => {
                                return Err(invalid(
                                    "circle offset reference orientation is undecided",
                                ));
                            }
                        };
                        let radius = &circle.radius_sq.sqrt() - &(&distance * &direction);
                        if !matches!(dec.sign(&radius), Outcome::Proven(StrictSign::Positive)) {
                            return Err(invalid(
                                "circle offset radius is not proved positive; signed phase requires finite mapping",
                            ));
                        }
                        shapes::Curve::Circle(shapes::Circle::new(
                            circle.center,
                            circle.normal,
                            radius.square(),
                            circle.ref_dir,
                            dec,
                        )?)
                    }
                    _ => unreachable!("selected line or circle"),
                }
            }
            Curve::BSpline(_) | Curve::Polyline(_) => shapes::Curve::BSpline(
                self.spline_curve(id, &std::array::from_fn(|_| length_mm.clone()), dec)?
                    .to_curve(dec)?,
            ),
            _ => {
                return Err(invalid(
                    "curve requires explicit wrapper or finite numerical mapping",
                ));
            }
        })
    }

    /// Map an unwrapped elementary surface to its complete exact native carrier.
    /// A spindle torus is not reduced to its outer sheet. Spline control data
    /// remain available through `spline_surface` without a regularity certificate.
    pub fn native_surface(
        &self,
        id: EntityId,
        length_mm: &Rational,
        dec: &mut Decider,
    ) -> Result<shapes::Surface, MappingError> {
        let Geometry::Surface(surface) = self.geometry(id, dec)? else {
            return Err(invalid("expected SURFACE"));
        };
        if let Surface::Offset(parent, distance, _) = &surface {
            let Geometry::Surface(base) = self.geometry(*parent, dec)? else {
                return Err(invalid("expected offset basis surface"));
            };
            let placement = match &base {
                Surface::Plane(p) | Surface::Cylinder(p, _) | Surface::Sphere(p, _) => *p,
                _ => return Err(invalid("offset surface requires finite numerical mapping")),
            };
            let f = self.axis_placement(placement, length_mm, dec)?;
            let distance = Real::from_rational(distance * length_mm);
            return Ok(match base {
                Surface::Plane(_) => shapes::Surface::Plane(shapes::Plane::new(
                    f.origin.plus(&f.axes[2].scaled(&distance)),
                    f.axes[2].clone(),
                    f.axes[0].clone(),
                    dec,
                )?),
                Surface::Cylinder(_, r) => shapes::Surface::Cylinder(shapes::Cylinder::new(
                    f.origin,
                    f.axes[2].clone(),
                    Real::from_rational(&r * length_mm) + distance,
                    f.axes[0].clone(),
                    dec,
                )?),
                Surface::Sphere(_, r) => shapes::Surface::Sphere(shapes::Sphere::new(
                    f.origin,
                    Real::from_rational(&r * length_mm) + distance,
                    f.axes[2].clone(),
                    f.axes[0].clone(),
                    dec,
                )?),
                _ => unreachable!("selected analytic offset"),
            });
        }
        let placement = match &surface {
            Surface::Plane(p)
            | Surface::Cylinder(p, _)
            | Surface::Sphere(p, _)
            | Surface::Torus(p, _, _) => *p,
            _ => return Err(invalid("surface requires finite or numerical mapping")),
        };
        let f = self.axis_placement(placement, length_mm, dec)?;
        Ok(match surface {
            Surface::Plane(_) => shapes::Surface::Plane(shapes::Plane::new(
                f.origin,
                f.axes[2].clone(),
                f.axes[0].clone(),
                dec,
            )?),
            Surface::Cylinder(_, r) => shapes::Surface::Cylinder(shapes::Cylinder::new(
                f.origin,
                f.axes[2].clone(),
                Real::from_rational(&r * length_mm),
                f.axes[0].clone(),
                dec,
            )?),
            Surface::Sphere(_, r) => shapes::Surface::Sphere(shapes::Sphere::new(
                f.origin,
                Real::from_rational(&r * length_mm),
                f.axes[2].clone(),
                f.axes[0].clone(),
                dec,
            )?),
            Surface::Torus(_, a, b) => shapes::Surface::Torus(shapes::Torus::new(
                f.origin,
                f.axes[2].clone(),
                Real::from_rational(&a.abs() * length_mm),
                Real::from_rational(&b.abs() * length_mm),
                f.axes[0].clone(),
                dec,
            )?),
            _ => unreachable!("selected elementary carrier"),
        })
    }
}

fn expanded_knots(
    data: &KnotData,
    degree: usize,
    count: usize,
    dec: &mut Decider,
) -> Result<Vec<Rational>, MappingError> {
    let len = count
        .checked_add(degree)
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| invalid("knot count overflow"))?;
    dec.charge_work(len as u64)?;
    if degree == 0 || degree >= count {
        return Err(invalid("invalid spline degree"));
    }
    match data {
        KnotData::Explicit(mult, values, _) => {
            if mult.len() != values.len() {
                return Err(invalid("knot and multiplicity counts differ"));
            }
            let mut out = Vec::with_capacity(len);
            for (m, v) in mult.iter().zip(values) {
                if *m == 0 || out.len().checked_add(*m).is_none_or(|n| n > len) {
                    return Err(invalid("invalid expanded knot count"));
                }
                out.extend(std::iter::repeat_n(v.clone(), *m));
            }
            if out.len() != len {
                return Err(invalid(
                    "expanded knot count differs from degree and controls",
                ));
            }
            Ok(out)
        }
        _ => {
            if matches!(data, KnotData::Bezier) && count != degree + 1 {
                return Err(invalid("Bezier control count differs from degree"));
            }
            Ok((0..len)
                .map(|i| {
                    Rational::from_integer(
                        (if matches!(data, KnotData::Uniform) {
                            i
                        } else {
                            i.saturating_sub(degree).min(count - degree)
                        } as u64)
                            .into(),
                    )
                })
                .collect())
        }
    }
}
