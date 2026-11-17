//! Source parameter charts for numerical fitting, with interval enclosures.
use super::expression::{self as ex, Expr, Triple};
use super::*;
use catallax_numbers::{Ball, Rational, Refusal};
use catallax_shapes::{Point3, Vector3};
use std::collections::BTreeSet;

/// Plane-angle unit of the containing representation.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum AngleUnit {
    /// Written angular parameters are radians.
    Radians,
    /// Written angular parameters are degrees.
    Degrees,
    /// Exact radians per written angular unit from a representation context.
    ScaledRadians(Rational),
}
impl AngleUnit {
    fn factor(&self) -> Expr {
        match self {
            Self::Radians => Expr::integer(1),
            Self::Degrees => Expr::pi().div(&Expr::integer(180)),
            Self::ScaledRadians(scale) => Expr::rational(scale.clone()),
        }
    }
}

/// Units of a chart's evaluated coordinates.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CoordinateSpace {
    /// Physical Cartesian coordinates in millimetres.
    Millimetres,
    /// Raw surface parameters, which may mix angular and length units.
    SurfaceParameters,
}

/// One source carrier chart. It is not a finite face, edge or validity proof.
/// Wrapper instances remain listed so topology admission must apply their trims.
#[derive(Clone, Debug)]
pub struct ParameterChart {
    source: EntityId,
    dimension: usize,
    domain: [[Option<Rational>; 2]; 2],
    wrappers: Vec<EntityId>,
    coordinates: Triple,
    coordinate_space: CoordinateSpace,
}
impl ParameterChart {
    /// Original geometry instance requested by the caller.
    pub fn source(&self) -> EntityId {
        self.source
    }
    /// One for a curve or two for a surface.
    pub fn dimension(&self) -> usize {
        self.dimension
    }
    /// Evaluated coordinate units, independent of curve parameter units.
    pub fn coordinate_space(&self) -> CoordinateSpace {
        self.coordinate_space
    }
    /// Basis parameter limits. `None` means unbounded, not an inferred trim.
    pub fn domain(&self) -> &[[Option<Rational>; 2]; 2] {
        &self.domain
    }
    /// Retained trim, composite-segment or replica instances, outermost last.
    pub fn wrappers(&self) -> &[EntityId] {
        &self.wrappers
    }
    /// Enclose the carrier at source parameters. This is a whole-interval
    /// enclosure, not a sampled error estimate or a source-to-body certificate.
    pub fn evaluate(
        &self,
        parameters: &[Ball],
        dec: &mut Decider,
    ) -> Result<[Ball; 3], MappingError> {
        self.evaluate_derivative(parameters, &[], dec)
    }
    /// Enclose a mixed parameter derivative. Repeated axes request higher order.
    pub fn evaluate_derivative(
        &self,
        parameters: &[Ball],
        axes: &[usize],
        dec: &mut Decider,
    ) -> Result<[Ball; 3], MappingError> {
        if parameters.len() != self.dimension || axes.iter().any(|i| *i >= self.dimension) {
            return Err(error("wrong parameter dimension"));
        }
        for (parameter, limits) in parameters.iter().zip(&self.domain) {
            if limits[0].as_ref().is_some_and(|lo| parameter.lo() < *lo)
                || limits[1].as_ref().is_some_and(|hi| parameter.hi() > *hi)
            {
                return Err(error("parameter interval is outside its source knot span"));
            }
        }
        let bits = parameters.iter().map(Ball::precision).min().unwrap_or(0);
        self.evaluate_expressions(parameters, axes, bits, dec)
    }

    /// Enclose an exact rational parameter box, checking its source-domain
    /// limits before outward arithmetic rounding of the parameters.
    pub fn evaluate_box(
        &self,
        bounds: &[[Rational; 2]],
        precision: u32,
        dec: &mut Decider,
    ) -> Result<[Ball; 3], MappingError> {
        self.evaluate_box_derivative(bounds, &[], precision, dec)
    }

    /// Enclose a mixed derivative over an exact rational parameter box.
    /// Polynomial extensions only bound rounding around the validated box;
    /// requests outside the source knot span still refuse.
    pub fn evaluate_box_derivative(
        &self,
        bounds: &[[Rational; 2]],
        axes: &[usize],
        precision: u32,
        dec: &mut Decider,
    ) -> Result<[Ball; 3], MappingError> {
        if bounds.len() != self.dimension || axes.iter().any(|i| *i >= self.dimension) {
            return Err(error("wrong parameter dimension"));
        }
        if precision == 0 || precision > dec.policy().max_precision_bits {
            return Err(Refusal::new(
                "numbers.precision-limit",
                "parameter evaluation exceeds the caller's precision allowance",
            )
            .into());
        }
        for (range, limits) in bounds.iter().zip(&self.domain) {
            if range[0] > range[1]
                || limits[0].as_ref().is_some_and(|lo| range[0] < *lo)
                || limits[1].as_ref().is_some_and(|hi| range[1] > *hi)
            {
                return Err(error("parameter interval is outside its source knot span"));
            }
        }
        let parameters = bounds
            .iter()
            .map(|[a, b]| {
                Ball::from_rational(a, precision).hull(&Ball::from_rational(b, precision))
            })
            .collect::<Vec<_>>();
        self.evaluate_expressions(&parameters, axes, precision, dec)
    }

    fn evaluate_expressions(
        &self,
        parameters: &[Ball],
        axes: &[usize],
        bits: u32,
        dec: &mut Decider,
    ) -> Result<[Ball; 3], MappingError> {
        let mut coordinates = self.coordinates.clone();
        for axis in axes {
            coordinates = ex::derivative(&coordinates, *axis, dec)?;
        }
        Ok([
            coordinates[0].evaluate(parameters, bits, dec)?,
            coordinates[1].evaluate(parameters, bits, dec)?,
            coordinates[2].evaluate(parameters, bits, dec)?,
        ])
    }
}
fn error(detail: &str) -> MappingError {
    Refusal::new("import.geometry-parameter-mapping", detail).into()
}
fn real_vector(v: &Vector3) -> Triple {
    [
        Expr::constant(v.x.clone()),
        Expr::constant(v.y.clone()),
        Expr::constant(v.z.clone()),
    ]
}
fn real_point(p: &Point3) -> Triple {
    [
        Expr::constant(p.x.clone()),
        Expr::constant(p.y.clone()),
        Expr::constant(p.z.clone()),
    ]
}
fn domain() -> [[Option<Rational>; 2]; 2] {
    [[None, None], [None, None]]
}
fn chart(id: EntityId, dimension: usize, coordinates: Triple) -> ParameterChart {
    ParameterChart {
        source: id,
        dimension,
        domain: domain(),
        wrappers: Vec::new(),
        coordinates,
        coordinate_space: CoordinateSpace::Millimetres,
    }
}
fn frame_point(f: &Frame, local: &Triple) -> Triple {
    let mut out = real_point(&f.origin);
    for (axis, value) in f.axes.iter().zip(local) {
        out = ex::add(&out, &ex::scale(&real_vector(axis), value));
    }
    out
}

impl ParsedStep {
    /// Read differentiable source curve carriers in their written parameter
    /// units. Spline spans and composite segments remain separate charts.
    /// Cartesian-only trims remain declarations for the topology mapper.
    pub fn curve_charts(
        &self,
        id: EntityId,
        length_mm: &Rational,
        angles: AngleUnit,
        dec: &mut Decider,
    ) -> Result<Vec<ParameterChart>, MappingError> {
        Compiler {
            source: self,
            length_mm: length_mm.clone(),
            angles,
            active: BTreeSet::new(),
        }
        .curve(id, 0, dec)
    }
    /// Read complete source surface carriers, including offsets, sweeps and
    /// spindle tori. Rectangular and curve bounds remain explicit wrappers.
    pub fn surface_charts(
        &self,
        id: EntityId,
        length_mm: &Rational,
        angles: AngleUnit,
        dec: &mut Decider,
    ) -> Result<Vec<ParameterChart>, MappingError> {
        Compiler {
            source: self,
            length_mm: length_mm.clone(),
            angles,
            active: BTreeSet::new(),
        }
        .surface(id, dec)
    }
}
struct Compiler<'a> {
    source: &'a ParsedStep,
    length_mm: Rational,
    angles: AngleUnit,
    active: BTreeSet<EntityId>,
}
impl Compiler<'_> {
    fn enter(&mut self, id: EntityId, dec: &mut Decider) -> Result<(), MappingError> {
        dec.charge_work(1)?;
        if self.active.len() >= 128 || !self.active.insert(id) {
            return Err(error(
                "cyclic or excessively deep intrinsic geometry references",
            ));
        }
        Ok(())
    }
    fn frame(&self, id: EntityId, dec: &mut Decider) -> Result<Frame, MappingError> {
        self.source.axis_placement(id, &self.length_mm, dec)
    }
    fn point(&self, id: EntityId, dec: &mut Decider) -> Result<Triple, MappingError> {
        Ok(self
            .source
            .cartesian_point(id, &std::array::from_fn(|_| self.length_mm.clone()), dec)?
            .map(Expr::rational))
    }
    fn direction(&self, id: EntityId, dec: &mut Decider) -> Result<Triple, MappingError> {
        Ok(real_vector(&self.source.unit_direction(id, dec)?))
    }
    fn length(&self, v: &Rational) -> Expr {
        Expr::rational(v * &self.length_mm)
    }
    fn torus_coordinates(
        &self,
        placement: EntityId,
        major: &Expr,
        minor: &Expr,
        angles: [&Expr; 2],
        dec: &mut Decider,
    ) -> Result<Triple, MappingError> {
        let [a, b] = angles;
        let radial = major.add(&minor.mul(&b.cos()));
        Ok(frame_point(
            &self.frame(placement, dec)?,
            &[
                radial.mul(&a.cos()),
                radial.mul(&a.sin()),
                minor.mul(&b.sin()),
            ],
        ))
    }
    fn curve(
        &mut self,
        id: EntityId,
        axis: usize,
        dec: &mut Decider,
    ) -> Result<Vec<ParameterChart>, MappingError> {
        self.enter(id, dec)?;
        let Geometry::Curve(curve) = self.source.geometry(id, dec)? else {
            return Err(error("expected CURVE"));
        };
        let t = Expr::parameter(axis);
        let zero = Expr::integer(0);
        let mut out = match curve {
            Curve::Line(p, v) => {
                let Geometry::Vector(direction, magnitude) = self.source.geometry(v, dec)? else {
                    return Err(error("expected VECTOR"));
                };
                vec![chart(
                    id,
                    1,
                    ex::add(
                        &self.point(p, dec)?,
                        &ex::scale(
                            &self.direction(direction, dec)?,
                            &t.mul(&self.length(&magnitude)),
                        ),
                    ),
                )]
            }
            Curve::Circle(p, r) => {
                let a = t.mul(&self.angles.factor());
                vec![chart(
                    id,
                    1,
                    frame_point(
                        &self.frame(p, dec)?,
                        &[
                            self.length(&r).mul(&a.cos()),
                            self.length(&r).mul(&a.sin()),
                            zero,
                        ],
                    ),
                )]
            }
            Curve::Ellipse(p, a, b) => {
                let t = t.mul(&self.angles.factor());
                vec![chart(
                    id,
                    1,
                    frame_point(
                        &self.frame(p, dec)?,
                        &[
                            self.length(&a).mul(&t.cos()),
                            self.length(&b).mul(&t.sin()),
                            zero,
                        ],
                    ),
                )]
            }
            Curve::Hyperbola(p, a, b) => {
                let positive = t.exp();
                let negative = t.neg().exp();
                let cosh = positive.add(&negative).div(&Expr::integer(2));
                let sinh = positive.sub(&negative).div(&Expr::integer(2));
                vec![chart(
                    id,
                    1,
                    frame_point(
                        &self.frame(p, dec)?,
                        &[self.length(&a).mul(&cosh), self.length(&b).mul(&sinh), zero],
                    ),
                )]
            }
            Curve::Parabola(p, f) => {
                let t = t.mul(&Expr::rational(self.length_mm.clone()));
                vec![chart(
                    id,
                    1,
                    frame_point(
                        &self.frame(p, dec)?,
                        &[
                            t.pow(2).div(&Expr::integer(4).mul(&self.length(&f))),
                            t,
                            zero,
                        ],
                    ),
                )]
            }
            Curve::BSpline(_) | Curve::Polyline(_) => {
                let native = self
                    .source
                    .spline_curve(id, &std::array::from_fn(|_| self.length_mm.clone()), dec)?
                    .to_curve(dec)?;
                let mut charts = Vec::new();
                for (span, interval) in native.spans() {
                    let polys = native.span_polynomials(span, dec)?;
                    let p = polynomials(&polys, std::slice::from_ref(&t), dec)?;
                    let mut c = chart(id, 1, p);
                    c.domain[axis] = [Some(interval[0].clone()), Some(interval[1].clone())];
                    charts.push(c);
                }
                charts
            }
            Curve::Trimmed(parent, _, _, _) => self.curve(parent, axis, dec)?,
            Curve::SurfaceCurve(parent, _, _, _) => self.curve(parent, axis, dec)?,
            Curve::Replica(parent, transform) => {
                let mut charts = self.curve(parent, axis, dec)?;
                for c in &mut charts {
                    c.coordinates = self.transform(&c.coordinates, transform, dec)?;
                }
                charts
            }
            Curve::Offset3d(parent, distance, _, direction) => {
                let mut charts = self.curve(parent, axis, dec)?;
                let normal = self.direction(direction, dec)?;
                for c in &mut charts {
                    let tangent = ex::derivative(&c.coordinates, axis, dec)?;
                    // OCCT negates STEP distance before its tangent x reference convention.
                    let offset = ex::scale(
                        &ex::unit(&ex::cross(&tangent, &normal)),
                        &self.length(&distance).neg(),
                    );
                    c.coordinates = ex::add(&c.coordinates, &offset);
                }
                charts
            }
            Curve::Composite(segments, _) => {
                let mut charts = Vec::new();
                for segment in segments {
                    let Geometry::CompositeSegment(_, _, parent) =
                        self.source.geometry(segment, dec)?
                    else {
                        return Err(error("expected COMPOSITE_CURVE_SEGMENT"));
                    };
                    let mut part = self.curve(parent, axis, dec)?;
                    for c in &mut part {
                        c.wrappers.push(segment);
                    }
                    charts.extend(part);
                }
                charts
            }
            Curve::Pcurve(_, representation) => {
                let rep = self.source.representation(representation)?;
                let mut charts = Vec::new();
                let length = std::mem::replace(&mut self.length_mm, Rational::one());
                let angle = std::mem::replace(&mut self.angles, AngleUnit::Radians);
                for item in rep.items {
                    charts.extend(self.curve(item, axis, dec)?);
                }
                self.length_mm = length;
                self.angles = angle;
                for c in &mut charts {
                    c.coordinate_space = CoordinateSpace::SurfaceParameters;
                }
                charts
            }
        };
        for c in &mut out {
            if c.source != id {
                c.wrappers.push(id);
            }
            c.source = id;
        }
        self.active.remove(&id);
        Ok(out)
    }
    fn surface(
        &mut self,
        id: EntityId,
        dec: &mut Decider,
    ) -> Result<Vec<ParameterChart>, MappingError> {
        self.enter(id, dec)?;
        let Geometry::Surface(surface) = self.source.geometry(id, dec)? else {
            return Err(error("expected SURFACE"));
        };
        let u = Expr::parameter(0);
        let v = Expr::parameter(1);
        let a = u.mul(&self.angles.factor());
        let b = v.mul(&self.angles.factor());
        let zero = Expr::integer(0);
        let mut out = match surface {
            Surface::Plane(p) => vec![chart(
                id,
                2,
                frame_point(
                    &self.frame(p, dec)?,
                    &[
                        u.mul(&Expr::rational(self.length_mm.clone())),
                        v.mul(&Expr::rational(self.length_mm.clone())),
                        zero,
                    ],
                ),
            )],
            Surface::Cylinder(p, r) => vec![chart(
                id,
                2,
                frame_point(
                    &self.frame(p, dec)?,
                    &[
                        self.length(&r).mul(&a.cos()),
                        self.length(&r).mul(&a.sin()),
                        v.mul(&Expr::rational(self.length_mm.clone())),
                    ],
                ),
            )],
            Surface::Cone(p, r, angle) => {
                let angle = Expr::rational(angle).mul(&self.angles.factor());
                let height = v.mul(&Expr::rational(self.length_mm.clone()));
                let radius = self
                    .length(&r)
                    .add(&height.mul(&angle.sin().div(&angle.cos())));
                vec![chart(
                    id,
                    2,
                    frame_point(
                        &self.frame(p, dec)?,
                        &[radius.mul(&a.cos()), radius.mul(&a.sin()), height],
                    ),
                )]
            }
            Surface::Sphere(p, r) => {
                let radial = self.length(&r).mul(&b.cos());
                vec![chart(
                    id,
                    2,
                    frame_point(
                        &self.frame(p, dec)?,
                        &[
                            radial.mul(&a.cos()),
                            radial.mul(&a.sin()),
                            self.length(&r).mul(&b.sin()),
                        ],
                    ),
                )]
            }
            Surface::Torus(p, major, minor) => vec![chart(
                id,
                2,
                self.torus_coordinates(
                    p,
                    &self.length(&major.abs()),
                    &self.length(&minor.abs()),
                    [&a, &b],
                    dec,
                )?,
            )],
            Surface::BSpline(_) => {
                let native = self
                    .source
                    .spline_surface(id, &self.length_mm, dec)?
                    .to_surface(dec)?;
                let mut charts = Vec::new();
                for span in native.span_indices() {
                    let p = polynomials(
                        &native.span_polynomials(span, dec)?,
                        &[u.clone(), v.clone()],
                        dec,
                    )?;
                    let mut c = chart(id, 2, p);
                    for axis in 0..2 {
                        let knots = native.knots()[axis];
                        c.domain[axis] = [
                            Some(knots[span[axis]].clone()),
                            Some(knots[span[axis] + 1].clone()),
                        ];
                    }
                    charts.push(c);
                }
                charts
            }
            Surface::Extrusion(profile, vector) => {
                let Geometry::Vector(direction, magnitude) = self.source.geometry(vector, dec)?
                else {
                    return Err(error("expected extrusion VECTOR"));
                };
                if magnitude == Rational::zero() {
                    return Err(error("extrusion VECTOR is zero"));
                }
                let direction = self.direction(direction, dec)?;
                let distance = v
                    .mul(&Expr::rational(self.length_mm.clone()))
                    .mul(&Expr::integer(if magnitude < Rational::zero() {
                        -1
                    } else {
                        1
                    }));
                let mut charts = self.curve(profile, 0, dec)?;
                for c in &mut charts {
                    if c.coordinate_space != CoordinateSpace::Millimetres {
                        return Err(error("extrusion profile is a parameter-space curve"));
                    }
                    c.dimension = 2;
                    c.coordinates = ex::add(&c.coordinates, &ex::scale(&direction, &distance));
                }
                charts
            }
            Surface::Revolution(profile, placement) => {
                let Geometry::Axis1Placement(origin, axis) =
                    self.source.geometry(placement, dec)?
                else {
                    return Err(error("expected AXIS1_PLACEMENT"));
                };
                let origin = self.point(origin, dec)?;
                let axis = match axis {
                    Some(id) => self.direction(id, dec)?,
                    None => [zero.clone(), zero.clone(), Expr::integer(1)],
                };
                let mut charts = self.curve(profile, 1, dec)?;
                for c in &mut charts {
                    if c.coordinate_space != CoordinateSpace::Millimetres {
                        return Err(error("revolution profile is a parameter-space curve"));
                    }
                    let relative = ex::sub(&c.coordinates, &origin);
                    let axial = ex::scale(&axis, &ex::dot(&axis, &relative));
                    let radial = ex::sub(&relative, &axial);
                    c.coordinates = ex::add(
                        &origin,
                        &ex::add(
                            &axial,
                            &ex::add(
                                &ex::scale(&radial, &a.cos()),
                                &ex::scale(&ex::cross(&axis, &radial), &a.sin()),
                            ),
                        ),
                    );
                    c.dimension = 2;
                }
                charts
            }
            Surface::RectangularTrim(parent, _, _) | Surface::CurveBounded(parent, _, _) => {
                self.surface(parent, dec)?
            }
            Surface::Replica(parent, transform) => {
                let mut charts = self.surface(parent, dec)?;
                for c in &mut charts {
                    c.coordinates = self.transform(&c.coordinates, transform, dec)?;
                }
                charts
            }
            Surface::Offset(parent, distance, _) => {
                let ring = match self.source.geometry(parent, dec)? {
                    Geometry::Surface(Surface::Torus(p, major, minor)) => {
                        let major = major.abs();
                        let minor = minor.abs();
                        (major > minor && !minor.is_zero()).then_some((p, major, minor))
                    }
                    _ => None,
                };
                if let Some((p, major, minor)) = ring {
                    // The regular basis normal is cos(v)*radial + sin(v)*axis.
                    // Keep signed r+d and the original UV even across a collapse;
                    // no equivalent checked torus or embedding is asserted.
                    vec![chart(
                        parent,
                        2,
                        self.torus_coordinates(
                            p,
                            &self.length(&major),
                            &self.length(&(&minor + &distance)),
                            [&a, &b],
                            dec,
                        )?,
                    )]
                } else {
                    let mut charts = self.surface(parent, dec)?;
                    for c in &mut charts {
                        let du = ex::derivative(&c.coordinates, 0, dec)?;
                        let dv = ex::derivative(&c.coordinates, 1, dec)?;
                        c.coordinates = ex::add(
                            &c.coordinates,
                            &ex::scale(&ex::unit(&ex::cross(&du, &dv)), &self.length(&distance)),
                        );
                    }
                    charts
                }
            }
            Surface::RectangularComposite(rows) => {
                let mut charts = Vec::new();
                for patch in rows.into_iter().flatten() {
                    let Geometry::SurfacePatch(parent, _, _) = self.source.geometry(patch, dec)?
                    else {
                        return Err(error("expected SURFACE_PATCH"));
                    };
                    let mut part = self.surface(parent, dec)?;
                    for c in &mut part {
                        c.wrappers.push(patch);
                    }
                    charts.extend(part);
                }
                charts
            }
        };
        for c in &mut out {
            if c.source != id {
                c.wrappers.push(id);
            }
            c.source = id;
        }
        self.active.remove(&id);
        Ok(out)
    }
    fn transform(
        &self,
        p: &Triple,
        id: EntityId,
        dec: &mut Decider,
    ) -> Result<Triple, MappingError> {
        let Geometry::Transformation(t) = self.source.geometry(id, dec)? else {
            return Err(error("expected Cartesian transformation"));
        };
        let x = match t.axes[0] {
            Some(id) => self.direction(id, dec)?,
            None => [Expr::integer(1), Expr::integer(0), Expr::integer(0)],
        };
        let origin = self.point(t.origin, dec)?;
        if t.dimension == Some(2) {
            if t.axes[1].is_some() || t.scales[0].as_ref().is_some_and(|s| *s != Rational::one()) {
                return Err(error(
                    "2D replica scale and second-axis interpretation require explicit mapping",
                ));
            }
            let y = [x[1].neg(), x[0].clone(), Expr::integer(0)];
            return Ok(ex::add(
                &origin,
                &ex::add(&ex::scale(&x, &p[0]), &ex::scale(&y, &p[1])),
            ));
        }
        let y = match t.axes[1] {
            Some(id) => self.direction(id, dec)?,
            None => [Expr::integer(0), Expr::integer(1), Expr::integer(0)],
        };
        let z = match t.axes[2] {
            Some(id) => self.direction(id, dec)?,
            None => ex::unit(&ex::cross(&x, &y)),
        };
        let x = ex::unit(&ex::sub(&x, &ex::scale(&z, &ex::dot(&x, &z))));
        let y = ex::cross(&z, &x);
        let scale = t.scales[0].clone().unwrap_or_else(Rational::one);
        if scale == Rational::zero() {
            return Err(error("Cartesian transformation scale is zero"));
        }
        let rotated = ex::add(
            &ex::add(&ex::scale(&x, &p[0]), &ex::scale(&y, &p[1])),
            &ex::scale(&z, &p[2]),
        );
        Ok(ex::scale(
            &ex::add(&origin, &rotated),
            &Expr::integer(1).div(&Expr::rational(scale)),
        ))
    }
}
fn polynomials(
    p: &[catallax_numbers::RealPoly; 4],
    vars: &[Expr],
    dec: &mut Decider,
) -> Result<Triple, MappingError> {
    let denominator = Expr::polynomial(&p[3], vars, dec)?;
    Ok([
        Expr::polynomial(&p[0], vars, dec)?.div(&denominator),
        Expr::polynomial(&p[1], vars, dec)?.div(&denominator),
        Expr::polynomial(&p[2], vars, dec)?.div(&denominator),
    ])
}
