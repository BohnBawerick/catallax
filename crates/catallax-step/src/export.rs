use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::rc::Rc;

use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Answer, Answered, BoundKind, Decider, Mark, Rational, Real, Refusal};
use catallax_shapes::{CurveId, Point3, Surface, SurfaceId, Vector3};
use catallax_topology::{BodyId, BodyKind, EdgeId, FaceId, Loop, Model, Sense, ShellId, VertexId};

use crate::graph::{Graph, Id, boolean, refs, string, text};
use crate::number::{self, Number, Point};
use crate::{StepOptions, StepReport};

type Bounds = [[Rational; 2]; 3];

struct Vertex {
    id: Id,
    point: Point,
}

struct Edge {
    id: Id,
    bounds: Bounds,
}

#[derive(Default)]
struct Topology {
    vertices: BTreeMap<VertexId, Rc<Vertex>>,
    edges: BTreeMap<EdgeId, Rc<Edge>>,
}

struct Curve {
    id: Id,
    error: Rational,
    line_origin: Option<[Rational; 3]>,
    bounds: Option<Bounds>,
    // Semi-axis ratio bounds the inverse parameter map of an ellipse.
    aspect: Rational,
    diameter: Rational,
}

struct Placement {
    id: Id,
    origin: Point,
    error: Rational,
}

struct SurfaceGeometry {
    id: Id,
    origin: Bounds,
    base_error: Rational,
    extent_error: Rational,
}

pub(crate) struct Exporter<'a> {
    model: &'a Model,
    options: &'a StepOptions,
    dec: &'a mut Decider,
    graph: Graph,
    curves: BTreeMap<CurveId, Option<Rc<Curve>>>,
    surfaces: BTreeMap<SurfaceId, SurfaceGeometry>,
    report: StepReport,
}

impl<'a> Exporter<'a> {
    pub(crate) fn new(model: &'a Model, options: &'a StepOptions, dec: &'a mut Decider) -> Self {
        Self {
            model,
            options,
            dec,
            graph: Graph::default(),
            curves: BTreeMap::new(),
            surfaces: BTreeMap::new(),
            report: StepReport::default(),
        }
    }

    fn digits(&self) -> u32 {
        self.options.significant_digits.max(1)
    }

    fn deviation(&mut self, error: Rational) {
        self.report.max_deviation = self.report.max_deviation.clone().max(error);
    }

    fn point_entity(&mut self, point: &Point) -> Id {
        self.graph.add(
            "CARTESIAN_POINT",
            vec![
                text("''"),
                text(format!(
                    "({},{},{})",
                    point.xyz[0].text, point.xyz[1].text, point.xyz[2].text
                )),
            ],
        )
    }

    fn direction_entity(&mut self, point: &Point) -> Id {
        self.graph.add(
            "DIRECTION",
            vec![
                text("''"),
                text(format!(
                    "({},{},{})",
                    point.xyz[0].text, point.xyz[1].text, point.xyz[2].text
                )),
            ],
        )
    }

    fn placement(
        &mut self,
        origin: &Point3,
        axis: &Vector3,
        reference: &Vector3,
    ) -> Result<Placement, Refusal> {
        let p = number::point(origin, self.digits(), self.dec)?;
        let z = number::direction(axis, self.digits(), self.dec)?;
        let x = number::direction(reference, self.digits(), self.dec)?;
        // STEP normalizes the axis, then projects and normalizes the reference.
        // For a unit u, |normalize(w)-u| <= 2|w-u|. Gram-Schmidt gives
        // |dx|+|dy|+|dz| <= 8|reference-e1|+12|axis-e3|.
        let error = x.error() * Rational::from_i64(8) + z.error() * Rational::from_i64(12);
        let location = self.point_entity(&p);
        let axis = self.direction_entity(&z);
        let reference = self.direction_entity(&x);
        let id = self.graph.add(
            "AXIS2_PLACEMENT_3D",
            vec![text("''"), location.into(), axis.into(), reference.into()],
        );
        Ok(Placement {
            id,
            origin: p,
            error,
        })
    }

    fn vertex_point(&mut self, point: Point) -> Rc<Vertex> {
        self.deviation(point.error());
        let p = self.point_entity(&point);
        let id = self.graph.add("VERTEX_POINT", vec![text("''"), p.into()]);
        self.report.vertices += 1;
        Rc::new(Vertex { id, point })
    }

    fn added_vertex(&mut self, point: &Point3) -> Result<Rc<Vertex>, Refusal> {
        let point = number::point(point, self.digits(), self.dec)?;
        self.report.added_vertices += 1;
        Ok(self.vertex_point(point))
    }

    fn vertex(&mut self, id: VertexId, topology: &mut Topology) -> Result<Rc<Vertex>, Refusal> {
        if let Some(vertex) = topology.vertices.get(&id) {
            return Ok(vertex.clone());
        }
        let point_id = self.model.vertex(id).point;
        // The public enclosure call is required even for explicit rational points.
        let [x, y, z] = self.model.geom().enclose_point(
            point_id,
            &self
                .options
                .curve_tolerance
                .abs()
                .max(Rational::from_ratio(1, 1_000_000_000).unwrap()),
            self.dec,
        )?;
        let position = self.model.geom().position(point_id);
        let point = Point {
            xyz: [
                number::enclosed(&position.x, x, self.digits(), self.dec)?,
                number::enclosed(&position.y, y, self.digits(), self.dec)?,
                number::enclosed(&position.z, z, self.digits(), self.dec)?,
            ],
        };
        let vertex = self.vertex_point(point);
        topology.vertices.insert(id, vertex.clone());
        Ok(vertex)
    }

    fn analytic_curve(&mut self, id: CurveId) -> Result<Option<Rc<Curve>>, Refusal> {
        if let Some(curve) = self.curves.get(&id) {
            return Ok(curve.clone());
        }
        let form = curve_form(self.model.geom(), id, self.dec)?;
        let curve = self.curve_geometry(&form)?.map(Rc::new);
        self.curves.insert(id, curve.clone());
        Ok(curve)
    }

    fn curve_geometry(&mut self, form: &CurveForm) -> Result<Option<Curve>, Refusal> {
        let curve = match form {
            CurveForm::Line { origin, dir } => {
                let p = number::point(origin, self.digits(), self.dec)?;
                let d = number::direction(dir, self.digits(), self.dec)?;
                let location = self.point_entity(&p);
                let direction = self.direction_entity(&d);
                let vector = self
                    .graph
                    .add("VECTOR", vec![text("''"), direction.into(), text("1.")]);
                let id = self
                    .graph
                    .add("LINE", vec![text("''"), location.into(), vector.into()]);
                Curve {
                    id,
                    error: p.error(),
                    line_origin: Some(p.xyz.map(|n| n.value)),
                    bounds: None,
                    aspect: d.error() * Rational::from_i64(2),
                    diameter: Rational::zero(),
                }
            }
            CurveForm::Circle {
                center,
                normal,
                radius_sq,
                ref_dir,
            } => {
                let radius = number::real(&number::sqrt_real(radius_sq), self.digits(), self.dec)?;
                self.ellipse(center, normal, ref_dir, radius, None)?
            }
            CurveForm::Ellipse {
                center,
                normal,
                major_dir,
                major_sq,
                minor_sq,
            } => {
                let major = number::real(&number::sqrt_real(major_sq), self.digits(), self.dec)?;
                let minor = number::real(&number::sqrt_real(minor_sq), self.digits(), self.dec)?;
                self.ellipse(center, normal, major_dir, major, Some(minor))?
            }
            CurveForm::Parabola {
                vertex,
                normal,
                axis_dir,
                focal,
            } => {
                let focal = number::real(focal, self.digits(), self.dec)?;
                self.open_conic("PARABOLA", vertex, normal, axis_dir, vec![focal])?
            }
            CurveForm::Hyperbola {
                center,
                normal,
                axis_dir,
                transverse_sq,
                conjugate_sq,
                branch,
            } => {
                let major =
                    number::real(&number::sqrt_real(transverse_sq), self.digits(), self.dec)?;
                let minor =
                    number::real(&number::sqrt_real(conjugate_sq), self.digits(), self.dec)?;
                if *branch == catallax_numbers::StrictSign::Negative {
                    self.open_conic(
                        "HYPERBOLA",
                        center,
                        &normal.scaled(&Real::from_i64(-1)),
                        &axis_dir.scaled(&Real::from_i64(-1)),
                        vec![major, minor],
                    )?
                } else {
                    self.open_conic("HYPERBOLA", center, normal, axis_dir, vec![major, minor])?
                }
            }
            CurveForm::Branch => return Ok(None),
            _ => return Err(unsupported()),
        };
        Ok(Some(curve))
    }

    fn open_conic(
        &mut self,
        name: &'static str,
        center: &Point3,
        normal: &Vector3,
        reference: &Vector3,
        dimensions: Vec<Number>,
    ) -> Result<Curve, Refusal> {
        let placement = self.placement(center, normal, reference)?;
        let mut aspect = placement.error;
        let mut fields = vec![text("''"), placement.id.into()];
        for dimension in dimensions {
            let lower = &dimension.value - &dimension.error;
            if lower <= Rational::zero() {
                return Err(number::budget(
                    "rounded conic dimension has no positive lower bound",
                ));
            }
            aspect = aspect + &dimension.error / &lower;
            fields.push(text(dimension.text));
        }
        let error = placement.origin.error() * (Rational::one() + aspect.clone());
        Ok(Curve {
            id: self.graph.add(name, fields),
            error,
            line_origin: Some(placement.origin.xyz.map(|n| n.value)),
            bounds: None,
            aspect,
            diameter: Rational::zero(),
        })
    }

    fn ellipse(
        &mut self,
        center: &Point3,
        axis: &Vector3,
        reference: &Vector3,
        major: Number,
        minor: Option<Number>,
    ) -> Result<Curve, Refusal> {
        if major.value <= Rational::zero() {
            return Err(number::budget("rounded radius is not positive"));
        }
        let placement = self.placement(center, axis, reference)?;
        let radius_upper = &major.value.abs() + &major.error;
        let mut error =
            placement.origin.error() + &radius_upper * &placement.error + major.error.clone();
        let mut fields = vec![text("''"), placement.id.into(), text(major.text.clone())];
        let (name, aspect) = if let Some(minor) = minor {
            if minor.value <= Rational::zero() || major.value < minor.value {
                return Err(number::budget(
                    "rounded ellipse axes are not ordered and positive",
                ));
            }
            let aspect = &radius_upper / &minor.value;
            error = &error + &minor.error;
            fields.push(text(minor.text));
            ("ELLIPSE", aspect)
        } else {
            ("CIRCLE", &radius_upper / &major.value)
        };
        let id = self.graph.add(name, fields);
        let bounds = std::array::from_fn(|i| {
            let n = &placement.origin.xyz[i];
            [
                &(&n.value - &n.error) - &radius_upper,
                &(&n.value + &n.error) + &radius_upper,
            ]
        });
        Ok(Curve {
            id,
            error,
            line_origin: None,
            bounds: Some(bounds),
            aspect,
            diameter: radius_upper * Rational::from_i64(2),
        })
    }

    fn edge_entity(&mut self, start: Id, end: Id, curve: Id, same: bool) -> Id {
        self.report.edges += 1;
        self.graph.add(
            "EDGE_CURVE",
            vec![
                text("''"),
                start.into(),
                end.into(),
                curve.into(),
                boolean(same),
            ],
        )
    }

    fn edge(&mut self, id: EdgeId, topology: &mut Topology) -> Result<Rc<Edge>, Refusal> {
        if let Some(edge) = topology.edges.get(&id) {
            return Ok(edge.clone());
        }
        let source = self.model.edge(id).clone();
        let (start, end) = match (source.start, source.end) {
            (Some(a), Some(b)) => (self.vertex(a, topology)?, self.vertex(b, topology)?),
            (None, None) => {
                if !catallax_intersect::curve_is_closed(self.model.geom(), source.curve, self.dec)?
                {
                    return Err(unsupported());
                }
                let origin = self
                    .model
                    .geom()
                    .curve(source.curve)
                    .parameter_origin(self.model.geom());
                let v = self.added_vertex(&origin)?;
                (v.clone(), v)
            }
            _ => panic!("edge has only one endpoint"),
        };
        let (curve, bounds, same) = if let Some(curve) = self.analytic_curve(source.curve)? {
            let bounds = curve
                .bounds
                .clone()
                .unwrap_or_else(|| endpoint_bounds(&start.point, &end.point));
            let error = if curve.bounds.is_none() {
                let origin = curve.line_origin.as_ref().expect("line origin");
                let extent = start.point.extent(origin).max(end.point.extent(origin));
                (&curve.error + &(&extent * &curve.aspect)) * Rational::from_i64(2)
                    + start.point.error().max(end.point.error())
            } else {
                arc_error(&curve, &start, &end)
            };
            self.deviation(error);
            (curve.id, bounds, source.sense == Sense::Same)
        } else {
            let from = source.start.map(|v| self.model.vertex(v).point);
            let to = source.end.map(|v| self.model.vertex(v).point);
            let (from, to) = if source.sense == Sense::Same {
                (from, to)
            } else {
                (to, from)
            };
            let trace =
                catallax_intersect::trace(self.model.geom(), source.curve, from, to, self.dec)?
                    .value;
            let fit = crate::fit::fit(
                &trace,
                &start.point,
                &end.point,
                source.sense == Sense::Reversed,
                self.options,
                self.dec,
            )?;
            self.deviation(fit.error);
            let points = fit
                .points
                .iter()
                .map(|p| self.point_entity(p))
                .collect::<Vec<_>>();
            let knots = fit
                .knots
                .iter()
                .map(|k| number::rational(k, self.digits()).text)
                .collect::<Vec<_>>()
                .join(",");
            let multiplicities = fit
                .multiplicities
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",");
            let curve = self.graph.add(
                "B_SPLINE_CURVE_WITH_KNOTS",
                vec![
                    text("''"),
                    text(fit.degree.to_string()),
                    refs(points),
                    text(".UNSPECIFIED."),
                    boolean(start.id == end.id),
                    boolean(false),
                    text(format!("({multiplicities})")),
                    text(format!("({knots})")),
                    text(".UNSPECIFIED."),
                ],
            );
            self.report.approximated_curves += 1;
            (curve, fit.bounds, true)
        };
        let edge = Rc::new(Edge {
            id: self.edge_entity(start.id, end.id, curve, same),
            bounds,
        });
        topology.edges.insert(id, edge.clone());
        Ok(edge)
    }

    fn oriented(&mut self, edge: Id, same: bool) -> Id {
        self.graph.add(
            "ORIENTED_EDGE",
            vec![text("''"), text("*"), text("*"), edge.into(), boolean(same)],
        )
    }

    fn added_circle(
        &mut self,
        center: &Point3,
        axis: &Vector3,
        reference: &Vector3,
        radius: &Real,
        start: &Vertex,
        end: &Vertex,
    ) -> Result<Id, Refusal> {
        let curve = self
            .curve_geometry(&CurveForm::Circle {
                center: center.clone(),
                normal: axis.clone(),
                radius_sq: radius * radius,
                ref_dir: reference.clone(),
            })?
            .unwrap();
        self.deviation(arc_error(&curve, start, end));
        self.report.added_edges += 1;
        Ok(self.edge_entity(start.id, end.id, curve.id, true))
    }

    fn loopless(&mut self, surface: &Surface) -> Result<Id, Refusal> {
        let (edges, orientations) = match surface {
            Surface::Sphere(s) => {
                let axis = number::unit(&s.axis);
                let low = number::point_plus(&s.center, &axis, &(-&s.radius));
                let high = number::point_plus(&s.center, &axis, &s.radius);
                let a = self.added_vertex(&low)?;
                let b = self.added_vertex(&high)?;
                let normal = s.axis.cross(&s.ref_dir).scaled(&Real::from_i64(-1));
                let edge = self.added_circle(&s.center, &normal, &s.ref_dir, &s.radius, &a, &b)?;
                (vec![edge, edge], vec![true, false])
            }
            Surface::Torus(s) => {
                let radial = number::unit(&s.ref_dir);
                let center = number::point_plus(&s.center, &radial, &s.major);
                let radius = &s.major + &s.minor;
                let origin = number::point_plus(&s.center, &radial, &radius);
                let vertex = self.added_vertex(&origin)?;
                let normal = s.axis.cross(&s.ref_dir).scaled(&Real::from_i64(-1));
                let tube =
                    self.added_circle(&center, &normal, &s.ref_dir, &s.minor, &vertex, &vertex)?;
                let equator =
                    self.added_circle(&s.center, &s.axis, &s.ref_dir, &radius, &vertex, &vertex)?;
                (
                    vec![tube, equator, tube, equator],
                    vec![true, true, false, false],
                )
            }
            _ => return Err(unsupported()),
        };
        let oriented = edges
            .into_iter()
            .zip(orientations)
            .map(|(edge, same)| self.oriented(edge, same))
            .collect::<Vec<_>>();
        Ok(self
            .graph
            .add("EDGE_LOOP", vec![text("''"), refs(oriented)]))
    }

    fn surface(&mut self, id: SurfaceId) -> Result<(), Refusal> {
        if self.surfaces.contains_key(&id) {
            return Ok(());
        }
        let surface = self.model.geom().surface(id);
        if let Surface::EllipticCylinder(s) = surface {
            let major = number::real(&number::sqrt_real(&s.major_sq), self.digits(), self.dec)?;
            let minor = number::real(&number::sqrt_real(&s.minor_sq), self.digits(), self.dec)?;
            let ellipse = self.ellipse(&s.origin, &s.axis, &s.major_dir, major, Some(minor))?;
            let axis = number::direction(&s.axis, self.digits(), self.dec)?;
            let direction = self.direction_entity(&axis);
            let vector = self
                .graph
                .add("VECTOR", vec![text("''"), direction.into(), text("1.")]);
            let entity = self.graph.add(
                "SURFACE_OF_LINEAR_EXTRUSION",
                vec![text("''"), ellipse.id.into(), vector.into()],
            );
            let origin = self.origin_bounds(&s.origin)?;
            self.surfaces.insert(
                id,
                SurfaceGeometry {
                    id: entity,
                    origin,
                    base_error: ellipse.error,
                    extent_error: axis.error() * Rational::from_i64(2),
                },
            );
            return Ok(());
        }
        let (origin, axis, reference, name, sizes) = match surface {
            Surface::Plane(s) => (&s.origin, &s.normal, &s.ref_dir, "PLANE", vec![]),
            Surface::Cylinder(s) => (
                &s.origin,
                &s.axis,
                &s.ref_dir,
                "CYLINDRICAL_SURFACE",
                vec![s.radius()],
            ),
            Surface::Cone(s) => (&s.apex, &s.axis, &s.ref_dir, "CONICAL_SURFACE", vec![]),
            Surface::Sphere(s) => (
                &s.center,
                &s.axis,
                &s.ref_dir,
                "SPHERICAL_SURFACE",
                vec![s.radius.clone()],
            ),
            Surface::Torus(s) => (
                &s.center,
                &s.axis,
                &s.ref_dir,
                "TOROIDAL_SURFACE",
                vec![s.major.clone(), s.minor.clone()],
            ),
            _ => return Err(unsupported()),
        };
        let placement = self.placement(origin, axis, reference)?;
        let mut base_error = placement.origin.error();
        let mut extent_error = placement.error.clone();
        let mut fields = vec![text("''"), placement.id.into()];
        for size in &sizes {
            let n = number::real(size, self.digits(), self.dec)?;
            base_error = &base_error + &n.error;
            fields.push(text(n.text));
        }
        match surface {
            Surface::Sphere(s) => {
                base_error = base_error + &number::upper(&s.radius, self.dec)? * &extent_error;
                extent_error = Rational::zero();
            }
            Surface::Torus(s) => {
                base_error =
                    base_error + &number::upper(&(&s.major + &s.minor), self.dec)? * &extent_error;
                extent_error = Rational::zero();
            }
            Surface::Cylinder(s) => {
                base_error = base_error + &number::upper(&s.radius(), self.dec)? * &extent_error
            }
            Surface::Cone(s) => {
                let p = number::work_precision(
                    self.dec.policy(),
                    self.digits().saturating_mul(4).saturating_add(32),
                    2,
                )
                .ok_or_else(|| number::budget("no work available for the cone's arctangent"))?;
                let angle = self.dec.enclose_at(&s.slope, p)?.atan_checked(self.dec)?;
                // Rounding down keeps the written half-angle strictly acute,
                // including at one digit and for nearly flat cones.
                let n =
                    number::certified(number::floor_ball(&angle, self.digits()), self.digits())?;
                extent_error = (&extent_error + &n.error)
                    * (Rational::one() + number::upper(&s.slope, self.dec)?);
                fields.push(text("0."));
                fields.push(text(n.text));
            }
            _ => {}
        }
        let geometry = SurfaceGeometry {
            id: self.graph.add(name, fields),
            origin: self.origin_bounds(origin)?,
            base_error,
            extent_error,
        };
        self.surfaces.insert(id, geometry);
        Ok(())
    }

    fn origin_bounds(&mut self, point: &Point3) -> Result<Bounds, Refusal> {
        let mut result = std::array::from_fn(|_| [Rational::zero(), Rational::zero()]);
        for (output, value) in result.iter_mut().zip([&point.x, &point.y, &point.z]) {
            if let Some(value) = value.rational_ref() {
                *output = [value.clone(), value.clone()];
            } else {
                let value = number::real(value, self.digits(), self.dec)?;
                *output = [&value.value - &value.error, &value.value + &value.error];
            }
        }
        Ok(result)
    }

    fn face(&mut self, id: FaceId, reverse: bool, topology: &mut Topology) -> Result<Id, Refusal> {
        let face = self.model.face(id).clone();
        let mut bounds = Vec::new();
        let mut boundary_extent: Option<Bounds> = None;
        for loop_id in &face.loops {
            let loop_entity = match self.model.loop_(*loop_id).clone() {
                Loop::Vertex(v) => {
                    let vertex = self.vertex(v, topology)?;
                    include_bounds(&mut boundary_extent, &point_bounds(&vertex.point));
                    self.graph
                        .add("VERTEX_LOOP", vec![text("''"), vertex.id.into()])
                }
                Loop::Coedges(coedges) => {
                    let mut oriented = Vec::new();
                    for coedge in coedges {
                        let coedge = self.model.coedge(coedge);
                        let edge = self.edge(coedge.edge, topology)?;
                        include_bounds(&mut boundary_extent, &edge.bounds);
                        oriented.push(self.oriented(edge.id, coedge.sense == Sense::Same));
                    }
                    self.graph
                        .add("EDGE_LOOP", vec![text("''"), refs(oriented)])
                }
            };
            bounds.push(self.graph.add(
                "FACE_BOUND",
                vec![text("''"), loop_entity.into(), boolean(!reverse)],
            ));
        }
        if face.loops.is_empty() {
            let loop_entity = self.loopless(self.model.geom().surface(face.surface))?;
            bounds.push(self.graph.add(
                "FACE_BOUND",
                vec![text("''"), loop_entity.into(), boolean(!reverse)],
            ));
        }
        self.surface(face.surface)?;
        let surface = &self.surfaces[&face.surface];
        let extent = boundary_extent.map_or_else(Rational::zero, |bounds| {
            (0..3).fold(Rational::zero(), |s, i| {
                s + (&bounds[i][0] - &surface.origin[i][1])
                    .abs()
                    .max((&bounds[i][1] - &surface.origin[i][0]).abs())
            })
        });
        let error = &surface.base_error + &(&surface.extent_error * &extent);
        let surface_id = surface.id;
        self.deviation(error);
        self.report.faces += 1;
        Ok(self.graph.add(
            "ADVANCED_FACE",
            vec![
                text("''"),
                refs(bounds),
                surface_id.into(),
                boolean((face.sense == Sense::Same) != reverse),
            ],
        ))
    }

    fn shell(
        &mut self,
        id: ShellId,
        reverse: bool,
        topology: &mut Topology,
    ) -> Result<Id, Refusal> {
        let mut faces = Vec::new();
        for face in &self.model.shell(id).faces {
            faces.push(self.face(*face, reverse, topology)?);
        }
        Ok(self
            .graph
            .add("CLOSED_SHELL", vec![text("''"), refs(faces)]))
    }

    pub(crate) fn build(mut self, bodies: &[BodyId]) -> Answered<(String, StepReport)> {
        let mut solids = Vec::new();
        for id in bodies {
            let body = self.model.body(*id);
            match body.kind {
                BodyKind::Empty => continue,
                BodyKind::Solid => {}
                _ => return Err(Refusal::new("step.not-solid", "body is not a solid")),
            }
            for lump in &body.lumps {
                let lump = self.model.lump(*lump);
                let mut topology = Topology::default();
                let outer = self.shell(lump.shells[0], false, &mut topology)?;
                let solid = if lump.shells.len() == 1 {
                    self.graph
                        .add("MANIFOLD_SOLID_BREP", vec![text("''"), outer.into()])
                } else {
                    let mut voids = Vec::new();
                    for shell in &lump.shells[1..] {
                        let shell = self.shell(*shell, true, &mut topology)?;
                        voids.push(self.graph.add(
                            "ORIENTED_CLOSED_SHELL",
                            vec![text("''"), text("*"), shell.into(), boolean(false)],
                        ));
                    }
                    self.graph.add(
                        "BREP_WITH_VOIDS",
                        vec![text("''"), outer.into(), refs(voids)],
                    )
                };
                solids.push(solid);
                self.report.solids += 1;
            }
        }
        let uncertainty = self
            .options
            .curve_tolerance
            .clone()
            .max(self.report.max_deviation.clone());
        let uncertainty = number::decimal(&uncertainty.to_decimal_bounds(self.digits()).1);
        let mut output = self.graph.write(&solids, |numbers| {
            header(&self.options.product_name, numbers, &uncertainty)
        });
        output.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
        let mut answer = Answer::exact((output, self.report));
        if answer.value.1.max_deviation > self.options.curve_tolerance {
            answer.marks.push(Mark::new(
                "step.tolerance-not-reached",
                BoundKind::Distance,
                answer.value.1.max_deviation.clone(),
                "decimal rounding or spline fitting exceeded the requested tolerance",
            ));
        }
        Ok(answer)
    }
}

fn unsupported() -> Refusal {
    Refusal::new(
        "step.unsupported-geometry",
        "geometry has no specified STEP form",
    )
}

fn point_bounds(point: &Point) -> Bounds {
    point
        .xyz
        .each_ref()
        .map(|n| [&n.value - &n.error, &n.value + &n.error])
}

fn endpoint_bounds(start: &Point, end: &Point) -> Bounds {
    let mut bounds = Some(point_bounds(start));
    include_bounds(&mut bounds, &point_bounds(end));
    bounds.unwrap()
}

fn include_bounds(bounds: &mut Option<Bounds>, other: &Bounds) {
    if let Some(bounds) = bounds {
        for i in 0..3 {
            bounds[i][0] = bounds[i][0].clone().min(other[i][0].clone());
            bounds[i][1] = bounds[i][1].clone().max(other[i][1].clone());
        }
    } else {
        *bounds = Some(other.clone());
    }
}

pub(crate) fn header(product_name: &str, solids: &[usize], uncertainty: &str) -> String {
    let name = string(product_name);
    let filename = string(&format!("{product_name}.step"));
    let mut items = String::from("#17");
    for id in solids {
        write!(items, ",#{id}").unwrap();
    }
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('catallax export'),'2;1');\nFILE_NAME({filename},'1970-01-01T00:00:00',(''),(''),'catallax','catallax','');\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF {{ 1 0 10303 442 1 1 4 }}'));\nENDSEC;\nDATA;\n#1=APPLICATION_CONTEXT('managed model based 3d engineering');\n#2=APPLICATION_PROTOCOL_DEFINITION('international standard','ap242_managed_model_based_3d_engineering',2014,#1);\n#3=PRODUCT_CONTEXT('',#1,'mechanical');\n#4=PRODUCT({name},{name},'',(#3));\n#5=PRODUCT_DEFINITION_FORMATION('','',#4);\n#6=PRODUCT_DEFINITION_CONTEXT('part definition',#1,'design');\n#7=PRODUCT_DEFINITION('design','',#5,#6);\n#8=PRODUCT_DEFINITION_SHAPE('','',#7);\n#9=SHAPE_DEFINITION_REPRESENTATION(#8,#16);\n#10=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\n#11=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));\n#12=(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT());\n#13=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({uncertainty}),#10,'distance_accuracy_value','maximum deviation from the exact model');\n#14=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#13))GLOBAL_UNIT_ASSIGNED_CONTEXT((#10,#11,#12))REPRESENTATION_CONTEXT('',''));\n#15=PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#4));\n#16=ADVANCED_BREP_SHAPE_REPRESENTATION({name},({items}),#14);\n#17=AXIS2_PLACEMENT_3D('',#18,#19,#20);\n#18=CARTESIAN_POINT('',(0.,0.,0.));\n#19=DIRECTION('',(0.,0.,1.));\n#20=DIRECTION('',(1.,0.,0.));\n"
    )
}

fn arc_error(curve: &Curve, start: &Vertex, end: &Vertex) -> Rational {
    if start.id == end.id {
        return curve.error.clone();
    }
    let movement = start.point.error().max(end.point.error()) + curve.error.clone();
    if movement.is_zero() {
        return Rational::zero();
    }
    // The inverse written ellipse scales by at most 1 / written_minor.
    // `aspect` bounds exact_major / written_minor and written_major /
    // written_minor, so both endpoint order and endpoint displacement use
    // the same certificate. Without separation, cover the whole carrier.
    let separated = (0..3).any(|i| {
        (&start.point.xyz[i].value - &end.point.xyz[i].value).abs()
            > &start.point.xyz[i].error
                + &end.point.xyz[i].error
                + &movement * &curve.aspect * Rational::from_i64(4)
    });
    if separated {
        &curve.error + &(&movement * &curve.aspect * Rational::from_i64(4))
    } else {
        &curve.error + &curve.diameter
    }
}
