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

#[path = "export/charts.rs"]
mod charts;

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
    ring_vertices: BTreeMap<EdgeId, Rc<Vertex>>,
    chart: Option<charts::ChartTopology>,
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

struct CurveProduct {
    curve: Curve,
    written: Option<crate::pcurve::Curve>,
}

struct Placement {
    id: Id,
    origin: Point,
    error: Rational,
    written: Option<Rc<crate::pcurve::Placement>>,
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
    face_errors: BTreeMap<FaceId, Rational>,
    edge_errors: BTreeMap<EdgeId, Rational>,
    bands: crate::band::Plan,
    band_period: Option<Rc<Number>>,
    chart_state: Option<charts::State>,
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
            face_errors: BTreeMap::new(),
            edge_errors: BTreeMap::new(),
            bands: crate::band::Plan::default(),
            band_period: None,
            chart_state: None,
        }
    }

    pub(crate) fn regularizations(
        mut self,
        records: Vec<crate::StepRegularization>,
        face_errors: BTreeMap<FaceId, Rational>,
        edge_errors: BTreeMap<EdgeId, Rational>,
    ) -> Self {
        self.report.regularizations = records;
        self.face_errors = face_errors;
        self.edge_errors = edge_errors;
        self
    }

    pub(crate) fn charts(mut self, plan: crate::chart::Plan) -> Self {
        self.report.added_edges += plan.added_edges;
        self.report.added_vertices += plan.added_vertices;
        self.chart_state = Some(charts::State::new(plan));
        self
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
        capture: bool,
    ) -> Result<Placement, Refusal> {
        let p = number::point(origin, self.digits(), self.dec)?;
        let z = number::direction(axis, self.digits(), self.dec)?;
        let x = number::direction(reference, self.digits(), self.dec)?;
        // STEP normalizes the axis, then projects and normalizes the reference.
        // For a unit u, |normalize(w)-u| <= 2|w-u|. Gram-Schmidt gives
        // |dx|+|dy|+|dz| <= 8|reference-e1|+12|axis-e3|.
        let error = x.error() * Rational::from_i64(8) + z.error() * Rational::from_i64(12);
        let written = if capture {
            Some(Rc::new(crate::pcurve::Placement::new(
                p.xyz.each_ref().map(|n| n.value.clone()),
                z.xyz.each_ref().map(|n| n.value.clone()),
                x.xyz.each_ref().map(|n| n.value.clone()),
            )?))
        } else {
            None
        };
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
            written,
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
        let point = match self
            .chart_state
            .as_mut()
            .and_then(|state| state.plan.written_vertices.remove(&id))
        {
            Some(point) => point,
            None => number::vertex(self.model, id, self.options, self.dec)?,
        };
        let vertex = self.vertex_point(point);
        topology.vertices.insert(id, vertex.clone());
        Ok(vertex)
    }

    fn analytic_curve(&mut self, id: CurveId) -> Result<Option<Rc<Curve>>, Refusal> {
        if let Some(curve) = self.curves.get(&id) {
            return Ok(curve.clone());
        }
        let form = if let Some(form) = self.bands.curves.get(&id) {
            form.clone()
        } else {
            curve_form(self.model.geom(), id, self.dec)?
        };
        let capture = self
            .chart_state
            .as_ref()
            .is_some_and(|state| state.plan.curves.contains(&id));
        let curve = match self.curve_geometry(&form, capture)? {
            Some(product) => {
                if let Some(written) = product.written {
                    self.chart_state
                        .as_mut()
                        .expect("captured chart geometry")
                        .analytic
                        .insert(
                            id,
                            Rc::new(charts::Analytic {
                                written,
                                exact: form,
                            }),
                        );
                }
                Some(Rc::new(product.curve))
            }
            None => None,
        };
        self.curves.insert(id, curve.clone());
        Ok(curve)
    }

    fn curve_geometry(
        &mut self,
        form: &CurveForm,
        capture: bool,
    ) -> Result<Option<CurveProduct>, Refusal> {
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
                let written = capture.then(|| crate::pcurve::Curve::Line {
                    origin: p.xyz.each_ref().map(|n| n.value.clone()),
                    direction: d.xyz.each_ref().map(|n| n.value.clone()),
                });
                CurveProduct {
                    curve: Curve {
                        id,
                        error: p.error(),
                        line_origin: Some(p.xyz.map(|n| n.value)),
                        bounds: None,
                        aspect: d.error() * Rational::from_i64(2),
                        diameter: Rational::zero(),
                    },
                    written,
                }
            }
            CurveForm::Circle {
                center,
                normal,
                radius_sq,
                ref_dir,
            } => {
                let radius = number::real(&number::sqrt_real(radius_sq), self.digits(), self.dec)?;
                self.ellipse(center, normal, ref_dir, radius, None, capture)?
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
                self.ellipse(center, normal, major_dir, major, Some(minor), capture)?
            }
            CurveForm::Parabola {
                vertex,
                normal,
                axis_dir,
                focal,
            } => {
                let focal = number::real(focal, self.digits(), self.dec)?;
                CurveProduct {
                    curve: self.open_conic("PARABOLA", vertex, normal, axis_dir, vec![focal])?,
                    written: None,
                }
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
                    CurveProduct {
                        curve: self.open_conic(
                            "HYPERBOLA",
                            center,
                            &normal.scaled(&Real::from_i64(-1)),
                            &axis_dir.scaled(&Real::from_i64(-1)),
                            vec![major, minor],
                        )?,
                        written: None,
                    }
                } else {
                    CurveProduct {
                        curve: self.open_conic(
                            "HYPERBOLA",
                            center,
                            normal,
                            axis_dir,
                            vec![major, minor],
                        )?,
                        written: None,
                    }
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
        let placement = self.placement(center, normal, reference, false)?;
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
        capture: bool,
    ) -> Result<CurveProduct, Refusal> {
        if major.value <= Rational::zero() {
            return Err(number::budget("rounded radius is not positive"));
        }
        let placement = self.placement(center, axis, reference, capture)?;
        let written = placement
            .written
            .map(|placement| crate::pcurve::Curve::Ellipse {
                placement,
                major: major.value.clone(),
                minor: minor
                    .as_ref()
                    .map_or_else(|| major.value.clone(), |n| n.value.clone()),
            });
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
        Ok(CurveProduct {
            curve: Curve {
                id,
                error,
                line_origin: None,
                bounds: Some(bounds),
                aspect,
                diameter: radius_upper * Rational::from_i64(2),
            },
            written,
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
                let origin = self.bands.rings.get(&id).map_or_else(
                    || {
                        self.model
                            .geom()
                            .curve(source.curve)
                            .parameter_origin(self.model.geom())
                    },
                    |ring| ring.origin.clone(),
                );
                let v = self.added_vertex(&origin)?;
                if self.bands.rings.contains_key(&id) {
                    topology.ring_vertices.insert(id, v.clone());
                }
                (v.clone(), v)
            }
            _ => panic!("edge has only one endpoint"),
        };
        let capture = self.has_chart_edge(id, topology);
        let mut geometry = None;
        let (curve, bounds, same, error) = if let Some(curve) = self.analytic_curve(source.curve)? {
            if capture {
                let captured = self
                    .chart_state
                    .as_ref()
                    .expect("chart edge")
                    .analytic
                    .get(&source.curve)
                    .ok_or_else(unsupported)?
                    .clone();
                geometry = Some(charts::Geometry::Analytic(captured));
            }
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
            (
                curve.id,
                bounds,
                source.sense == Sense::Same,
                error + self.edge_errors.get(&id).cloned().unwrap_or_default(),
            )
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
            let (fit, spans) = if capture {
                let (fit, spans) = crate::fit::fit_with_ranges(
                    &trace,
                    &start.point,
                    &end.point,
                    source.sense == Sense::Reversed,
                    self.options,
                    self.dec,
                )?;
                (fit, Some(spans))
            } else {
                (
                    crate::fit::fit(
                        &trace,
                        &start.point,
                        &end.point,
                        source.sense == Sense::Reversed,
                        self.options,
                        self.dec,
                    )?,
                    None,
                )
            };
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
            if let Some(spans) = spans {
                geometry = Some(charts::Geometry::spline(&fit, trace, spans));
            }
            self.report.approximated_curves += 1;
            (curve, fit.bounds, true, fit.error)
        };
        let curve = if let Some(ring) = self.bands.rings.get(&id).cloned() {
            self.ring_curve(curve, &bounds, &ring, error)?
        } else if let Some(geometry) = geometry.as_ref() {
            self.chart_edge(
                id,
                charts::EdgeInput {
                    curve,
                    geometry,
                    bounds: &bounds,
                    error: &error,
                    same,
                },
                topology,
            )?
        } else {
            self.deviation(error);
            curve
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

    fn pcurve_line(
        &mut self,
        surface: Id,
        context: Id,
        origin: [&str; 2],
        direction: &'static str,
    ) -> Id {
        let point = self.graph.add(
            "CARTESIAN_POINT",
            vec![text("''"), text(format!("({},{})", origin[0], origin[1]))],
        );
        let direction = self
            .graph
            .add("DIRECTION", vec![text("''"), text(direction)]);
        let vector = self
            .graph
            .add("VECTOR", vec![text("''"), direction.into(), text("1.")]);
        let line = self
            .graph
            .add("LINE", vec![text("''"), point.into(), vector.into()]);
        let representation = self.graph.add(
            "DEFINITIONAL_REPRESENTATION",
            vec![text("''"), refs([line]), context.into()],
        );
        self.graph.add(
            "PCURVE",
            vec![text("''"), surface.into(), representation.into()],
        )
    }

    fn parameter_context(&mut self) -> Id {
        self.graph.add(
            "GEOMETRIC_REPRESENTATION_CONTEXT",
            vec![text("''"), text("''"), text("2")],
        )
    }

    fn surface_error(&self, id: SurfaceId, bounds: &Bounds) -> Rational {
        let surface = &self.surfaces[&id];
        let extent = (0..3).fold(Rational::zero(), |sum, i| {
            sum + (&bounds[i][0] - &surface.origin[i][1])
                .abs()
                .max((&bounds[i][1] - &surface.origin[i][0]).abs())
        });
        &surface.base_error + &(&surface.extent_error * &extent)
    }

    fn band_period(&mut self) -> Result<Rc<Number>, Refusal> {
        if let Some(period) = &self.band_period {
            return Ok(period.clone());
        }
        let precision = number::work_precision(
            self.dec.policy(),
            self.digits().saturating_mul(4).saturating_add(32),
            2,
        )
        .ok_or_else(|| number::budget("no work available for the band's period"))?;
        let pi = catallax_numbers::Ball::pi_checked(precision, self.dec)?;
        let period = Rc::new(number::certified(
            number::rounded(&pi.lo().mul_pow2(1), &pi.hi().mul_pow2(1), self.digits()),
            self.digits(),
        )?);
        self.band_period = Some(period.clone());
        Ok(period)
    }

    fn ring_curve(
        &mut self,
        curve: Id,
        bounds: &Bounds,
        ring: &crate::band::Ring,
        error: Rational,
    ) -> Result<Id, Refusal> {
        let context = self.parameter_context();
        let period = self.band_period()?;
        let mut pcurves = Vec::with_capacity(ring.associations.len());
        for association in &ring.associations {
            self.surface(association.surface)?;
            let v = number::real(&association.v, self.digits(), self.dec)?;
            let surface = self.surfaces[&association.surface].id;
            pcurves.push(self.pcurve_line(
                surface,
                context,
                [
                    if association.increasing {
                        "0."
                    } else {
                        &period.text
                    },
                    &v.text,
                ],
                if association.increasing {
                    "(1.,0.)"
                } else {
                    "(-1.,0.)"
                },
            ));
            let axial_error = &v.error * &number::upper(&association.v_scale, self.dec)?;
            let angular_error = if association.increasing {
                Rational::zero()
            } else {
                &period.error * &number::upper(&ring.radius, self.dec)?
            };
            self.deviation(
                &error
                    + &(&axial_error
                        + &angular_error
                        + self.surface_error(association.surface, bounds)),
            );
        }
        Ok(self.graph.add(
            "SURFACE_CURVE",
            vec![text("''"), curve.into(), refs(pcurves), text(".CURVE_3D.")],
        ))
    }

    fn band_loop(
        &mut self,
        band: &crate::band::Band,
        topology: &mut Topology,
    ) -> Result<(Id, Bounds), Refusal> {
        let low = self.edge(band.edges[0], topology)?;
        let high = self.edge(band.edges[1], topology)?;
        let mut bounds = Some(low.bounds.clone());
        include_bounds(&mut bounds, &high.bounds);
        let bounds = bounds.unwrap();
        let a = topology.ring_vertices[&band.edges[0]].clone();
        let b = topology.ring_vertices[&band.edges[1]].clone();
        let origin = &self.bands.rings[&band.edges[0]].origin;
        let form = CurveForm::Line {
            origin: origin.clone(),
            dir: band.seam_direction.clone(),
        };
        let curve = self.curve_geometry(&form, false)?.unwrap().curve;
        let period = self.band_period()?;
        let v = number::real(&band.low_v, self.digits(), self.dec)?;
        self.surface(band.surface)?;
        let surface = self.surfaces[&band.surface].id;
        let context = self.parameter_context();
        let pcurves = [&period.text[..], "0."]
            .map(|u| self.pcurve_line(surface, context, [u, &v.text], "(0.,1.)"));
        let seam_curve = self.graph.add(
            "SEAM_CURVE",
            vec![
                text("''"),
                curve.id.into(),
                refs(pcurves),
                text(".CURVE_3D."),
            ],
        );
        let seam = self.edge_entity(a.id, b.id, seam_curve, true);
        self.report.added_edges += 1;
        let origin = curve.line_origin.as_ref().unwrap();
        let extent = a.point.extent(origin).max(b.point.extent(origin));
        let line_error = (&curve.error + &(&extent * &curve.aspect)) * Rational::from_i64(2)
            + a.point.error().max(b.point.error());
        let parameter_error = &number::upper(&band.radius_bound, self.dec)? * &period.error
            + &v.error * &number::upper(&band.v_scale, self.dec)?;
        self.deviation(line_error + parameter_error + self.surface_error(band.surface, &bounds));
        let uses = [
            (high.id, !band.forward[1]),
            (seam, false),
            (low.id, band.forward[0]),
            (seam, true),
        ];
        let uses = uses.map(|(edge, same)| self.oriented(edge, same));
        Ok((
            self.graph.add("EDGE_LOOP", vec![text("''"), refs(uses)]),
            bounds,
        ))
    }

    fn added_seam(
        &mut self,
        form: &CurveForm,
        vertex: &Vertex,
        surface: Id,
        period: &str,
        along_u: bool,
        parameter_error: &Rational,
    ) -> Result<Id, Refusal> {
        let curve = self.curve_geometry(form, false)?.unwrap().curve;
        self.deviation(&arc_error(&curve, vertex, vertex) + parameter_error);
        let context = self.graph.add(
            "GEOMETRIC_REPRESENTATION_CONTEXT",
            vec![text("''"), text("''"), text("2")],
        );
        let pcurves = [period, "0."].map(|offset| {
            let point = self.graph.add(
                "CARTESIAN_POINT",
                vec![
                    text("''"),
                    text(if along_u {
                        format!("(0.,{offset})")
                    } else {
                        format!("({offset},0.)")
                    }),
                ],
            );
            let direction = self.graph.add(
                "DIRECTION",
                vec![
                    text("''"),
                    text(if along_u { "(1.,0.)" } else { "(0.,1.)" }),
                ],
            );
            let vector = self
                .graph
                .add("VECTOR", vec![text("''"), direction.into(), text("1.")]);
            let line = self
                .graph
                .add("LINE", vec![text("''"), point.into(), vector.into()]);
            let representation = self.graph.add(
                "DEFINITIONAL_REPRESENTATION",
                vec![text("''"), refs([line]), context.into()],
            );
            self.graph.add(
                "PCURVE",
                vec![text("''"), surface.into(), representation.into()],
            )
        });
        let seam = self.graph.add(
            "SEAM_CURVE",
            vec![
                text("''"),
                curve.id.into(),
                refs(pcurves),
                text(".CURVE_3D."),
            ],
        );
        self.report.added_edges += 1;
        Ok(self.edge_entity(vertex.id, vertex.id, seam, true))
    }

    fn loopless(&mut self, id: SurfaceId) -> Result<Id, Refusal> {
        self.surface(id)?;
        let surface = self.surfaces[&id].id;
        match self.model.geom().surface(id) {
            Surface::Sphere(s) => {
                // A sole vertex_loop denotes the surface's complete natural
                // domain, including its seam and both singular pole edges.
                let axis = number::unit(&s.axis);
                let low = number::point_plus(&s.center, &axis, &(-&s.radius));
                let vertex = self.added_vertex(&low)?;
                Ok(self
                    .graph
                    .add("VERTEX_LOOP", vec![text("''"), vertex.id.into()]))
            }
            Surface::Torus(s) => {
                let radial = number::unit(&s.ref_dir);
                let center = number::point_plus(&s.center, &radial, &s.major);
                let radius = &s.major + &s.minor;
                let origin = number::point_plus(&s.center, &radial, &radius);
                let vertex = self.added_vertex(&origin)?;
                let normal = s.axis.cross(&s.ref_dir).scaled(&Real::from_i64(-1));
                let precision = number::work_precision(
                    self.dec.policy(),
                    self.digits().saturating_mul(4).saturating_add(32),
                    2,
                )
                .ok_or_else(|| number::budget("no work available for the torus's period"))?;
                let pi = catallax_numbers::Ball::pi_checked(precision, self.dec)?;
                let period = number::certified(
                    number::rounded(&pi.lo().mul_pow2(1), &pi.hi().mul_pow2(1), self.digits()),
                    self.digits(),
                )?;
                // The surface map is Lipschitz in either angle with constant
                // R+r. Include decimal period rounding in each seam's bound.
                let parameter_error = &number::upper(&radius, self.dec)? * &period.error;
                let tube = self.added_seam(
                    &CurveForm::Circle {
                        center,
                        normal,
                        radius_sq: &s.minor * &s.minor,
                        ref_dir: s.ref_dir.clone(),
                    },
                    &vertex,
                    surface,
                    &period.text,
                    false,
                    &parameter_error,
                )?;
                let equator = self.added_seam(
                    &CurveForm::Circle {
                        center: s.center.clone(),
                        normal: s.axis.clone(),
                        radius_sq: &radius * &radius,
                        ref_dir: s.ref_dir.clone(),
                    },
                    &vertex,
                    surface,
                    &period.text,
                    true,
                    &parameter_error,
                )?;
                let oriented = [
                    (equator, false),
                    (tube, false),
                    (equator, true),
                    (tube, true),
                ]
                .map(|(edge, same)| self.oriented(edge, same));
                Ok(self
                    .graph
                    .add("EDGE_LOOP", vec![text("''"), refs(oriented)]))
            }
            _ => Err(unsupported()),
        }
    }

    fn surface(&mut self, id: SurfaceId) -> Result<(), Refusal> {
        if self.surfaces.contains_key(&id) {
            return Ok(());
        }
        let surface = self.model.geom().surface(id);
        if let Surface::EllipticCylinder(s) = surface {
            let major = number::real(&number::sqrt_real(&s.major_sq), self.digits(), self.dec)?;
            let minor = number::real(&number::sqrt_real(&s.minor_sq), self.digits(), self.dec)?;
            let ellipse = self
                .ellipse(&s.origin, &s.axis, &s.major_dir, major, Some(minor), false)?
                .curve;
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
        let capture = self
            .chart_state
            .as_ref()
            .is_some_and(|state| state.plan.surfaces.contains(&id));
        let placement = if let Some(reference) = self.bands.references.get(&id).cloned() {
            self.placement(origin, axis, &reference, capture)?
        } else {
            self.placement(origin, axis, reference, capture)?
        };
        let mut base_error = placement.origin.error();
        let mut extent_error = placement.error.clone();
        let mut fields = vec![text("''"), placement.id.into()];
        let mut parameters = None;
        let mut written = Vec::with_capacity(sizes.len());
        for size in &sizes {
            let n = number::real(size, self.digits(), self.dec)?;
            base_error = &base_error + &n.error;
            written.push(n.value);
            fields.push(text(n.text));
        }
        if capture && !written.is_empty() {
            parameters = Some(match (surface, written.as_slice()) {
                (Surface::Cylinder(_), [radius]) => crate::pcurve::SurfaceParameters::Cylinder {
                    radius: radius.clone(),
                },
                (Surface::Sphere(_), [radius]) => crate::pcurve::SurfaceParameters::Sphere {
                    radius: radius.clone(),
                },
                (Surface::Torus(_), [major, minor]) => crate::pcurve::SurfaceParameters::Torus {
                    major: major.clone(),
                    minor: minor.clone(),
                },
                _ => unreachable!(
                    "captured surface sizes are not cylindrical, spherical or toroidal"
                ),
            });
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
                if capture {
                    parameters = Some(crate::pcurve::SurfaceParameters::Cone {
                        angle: n.value.clone(),
                    });
                }
                fields.push(text("0."));
                fields.push(text(n.text));
            }
            _ => {}
        }
        if let Some(parameters) = parameters {
            self.chart_state
                .as_mut()
                .expect("captured chart surface")
                .surfaces
                .insert(
                    id,
                    Rc::new(crate::pcurve::Surface {
                        placement: placement.written.expect("captured chart placement"),
                        parameters,
                    }),
                );
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
        if let Some(band) = self.bands.faces.get(&id).cloned() {
            let (loop_entity, extent) = self.band_loop(&band, topology)?;
            boundary_extent = Some(extent);
            bounds.push(self.graph.add(
                "FACE_BOUND",
                // Keep the chart walk CCW and reverse the bound together with
                // the face, as the native writer does for cylindrical holes.
                vec![
                    text("''"),
                    loop_entity.into(),
                    boolean((face.sense == Sense::Same) != reverse),
                ],
            ));
        } else {
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
        }
        if face.loops.is_empty() {
            let loop_entity = self.loopless(face.surface)?;
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
        let error = &surface.base_error
            + &(&surface.extent_error * &extent)
            + self.face_errors.get(&id).cloned().unwrap_or_default();
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
            if self
                .chart_state
                .as_ref()
                .is_some_and(|state| state.plan.faces.contains_key(face))
            {
                self.chart_faces(*face, reverse, topology, &mut faces)?;
            } else {
                faces.push(self.face(*face, reverse, topology)?);
            }
        }
        Ok(self
            .graph
            .add("CLOSED_SHELL", vec![text("''"), refs(faces)]))
    }

    pub(crate) fn build(self, bodies: &[BodyId]) -> Answered<(String, StepReport)> {
        self.build_parts(bodies, None)
    }

    pub(crate) fn build_stand_in(
        self,
        envelope: BodyId,
        bodies: &[BodyId],
    ) -> Answered<(String, StepReport)> {
        self.build_parts(bodies, Some(envelope))
    }

    fn build_parts(
        mut self,
        bodies: &[BodyId],
        envelope: Option<BodyId>,
    ) -> Answered<(String, StepReport)> {
        if self.chart_state.is_none() {
            self.bands = crate::band::Plan::prepare(
                self.model,
                bodies.iter().copied().chain(envelope),
                self.dec,
            )?;
            // Charts are an optional encoding. Only an exhausted request
            // allowance or deadline may turn their failure into a refusal;
            // otherwise the scene keeps its legacy encoding. Charts that miss
            // the requested tolerance would turn an exact answer inexact, so
            // they are not used either. A trim passing close to a parameter
            // pole can need that many pcurve spans, so such a scene is tried
            // once more with every sphere on its perpendicular frame.
            if self.bands.faces.is_empty() {
                for reframe in [false, true] {
                    let charted = crate::partial::prepare(
                        self.model,
                        bodies.iter().copied().chain(envelope),
                        self.options,
                        reframe,
                        self.dec,
                    )
                    .and_then(|prepared| {
                        let Some(prepared) = prepared else {
                            return Ok(None);
                        };
                        let (face_errors, edge_errors) =
                            prepared.errors(self.face_errors.clone(), self.edge_errors.clone());
                        let envelope = envelope.map(|_| prepared.bodies[bodies.len()]);
                        Exporter::new(&prepared.model, self.options, self.dec)
                            .regularizations(
                                self.report.regularizations.clone(),
                                face_errors,
                                edge_errors,
                            )
                            .charts(prepared.plan)
                            .build_parts(&prepared.bodies[..bodies.len()], envelope)
                            .map(Some)
                    });
                    match charted {
                        Ok(Some(answer))
                            if answer.value.1.max_deviation <= self.options.curve_tolerance =>
                        {
                            return Ok(answer);
                        }
                        Ok(Some(_)) => {}
                        Ok(None) => break,
                        Err(refusal)
                            if !matches!(
                                refusal.code,
                                "numbers.operation-budget" | "numbers.deadline"
                            ) =>
                        {
                            break;
                        }
                        Err(refusal) => return Err(refusal),
                    }
                }
            }
        }
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
                let mut topology = self.topology(&lump.shells);
                let outer = self.shell(lump.shells[0], false, &mut topology)?;
                let solid = if lump.shells.len() == 1 {
                    self.graph.add(
                        "MANIFOLD_SOLID_BREP",
                        vec![
                            text(if envelope.is_some() {
                                "'unit-density mass surrogate'"
                            } else {
                                "''"
                            }),
                            outer.into(),
                        ],
                    )
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
                self.finish_chart_edges(&topology)?;
                solids.push(solid);
                self.report.solids += 1;
            }
        }
        if let Some(envelope) = envelope {
            let body = self.model.body(envelope);
            if body.kind != BodyKind::Solid {
                return Err(Refusal::new(
                    "step.not-solid",
                    "the source envelope must be a solid before surface-only export",
                ));
            }
            let mut shells = Vec::new();
            for lump in &body.lumps {
                let mut topology = self.topology(&self.model.lump(*lump).shells);
                for (index, shell) in self.model.lump(*lump).shells.iter().enumerate() {
                    shells.push(self.shell(*shell, index != 0, &mut topology)?);
                }
                self.finish_chart_edges(&topology)?;
            }
            solids.push(self.graph.add(
                "SHELL_BASED_SURFACE_MODEL",
                vec![text("'massless field envelope'"), refs(shells)],
            ));
        }
        let uncertainty = self
            .options
            .curve_tolerance
            .clone()
            .max(self.report.max_deviation.clone());
        let uncertainty = number::decimal(&uncertainty.to_decimal_bounds(self.digits()).1);
        let mut output = self.graph.write(&solids, |numbers| {
            if envelope.is_some() {
                stand_in_header(&self.options.product_name, numbers, &uncertainty)
            } else {
                header(&self.options.product_name, numbers, &uncertainty)
            }
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
    representation_header(
        product_name,
        solids,
        uncertainty,
        "ADVANCED_BREP_SHAPE_REPRESENTATION",
    )
}

fn stand_in_header(product_name: &str, items: &[usize], uncertainty: &str) -> String {
    representation_header(product_name, items, uncertainty, "SHAPE_REPRESENTATION")
}

fn representation_header(
    product_name: &str,
    solids: &[usize],
    uncertainty: &str,
    representation: &str,
) -> String {
    let name = string(product_name);
    let filename = string(&format!("{product_name}.step"));
    let mut items = String::from("#17");
    for id in solids {
        write!(items, ",#{id}").unwrap();
    }
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('catallax export'),'2;1');\nFILE_NAME({filename},'1970-01-01T00:00:00',(''),(''),'catallax','catallax','');\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF {{ 1 0 10303 442 1 1 4 }}'));\nENDSEC;\nDATA;\n#1=APPLICATION_CONTEXT('managed model based 3d engineering');\n#2=APPLICATION_PROTOCOL_DEFINITION('international standard','ap242_managed_model_based_3d_engineering',2014,#1);\n#3=PRODUCT_CONTEXT('',#1,'mechanical');\n#4=PRODUCT({name},{name},'',(#3));\n#5=PRODUCT_DEFINITION_FORMATION('','',#4);\n#6=PRODUCT_DEFINITION_CONTEXT('part definition',#1,'design');\n#7=PRODUCT_DEFINITION('design','',#5,#6);\n#8=PRODUCT_DEFINITION_SHAPE('','',#7);\n#9=SHAPE_DEFINITION_REPRESENTATION(#8,#16);\n#10=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\n#11=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));\n#12=(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT());\n#13=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({uncertainty}),#10,'distance_accuracy_value','maximum deviation from the exact model');\n#14=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#13))GLOBAL_UNIT_ASSIGNED_CONTEXT((#10,#11,#12))REPRESENTATION_CONTEXT('',''));\n#15=PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#4));\n#16={representation}({name},({items}),#14);\n#17=AXIS2_PLACEMENT_3D('',#18,#19,#20);\n#18=CARTESIAN_POINT('',(0.,0.,0.));\n#19=DIRECTION('',(0.,0.,1.));\n#20=DIRECTION('',(1.,0.,0.));\n"
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
