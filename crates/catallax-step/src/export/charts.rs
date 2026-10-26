use super::*;
use crate::{fit, pcurve, uv, wire::Side};
use catallax_intersect::CurveTrace;
use catallax_topology::CoedgeId;
use std::collections::BTreeSet;

pub(super) struct Analytic {
    pub written: pcurve::Curve,
    pub exact: CurveForm,
}

pub(super) enum Geometry {
    Analytic(Rc<Analytic>),
    Spline {
        curves: Vec<pcurve::Curve>,
        trace: CurveTrace,
        spans: Vec<fit::TraceSpan>,
    },
}

impl Geometry {
    pub(super) fn spline(fit: &fit::Fit, trace: CurveTrace, spans: Vec<fit::TraceSpan>) -> Self {
        let curves = (0..fit.knots.len() - 1)
            .map(|i| {
                let points = std::array::from_fn(|j| {
                    fit.points[i * fit::DEGREE + j]
                        .xyz
                        .each_ref()
                        .map(|n| n.value.clone())
                });
                pcurve::Curve::Polynomial {
                    powers: fit::bezier_power(points),
                    range: [fit.knots[i].clone(), fit.knots[i + 1].clone()],
                }
            })
            .collect();
        Self::Spline {
            curves,
            trace,
            spans,
        }
    }
}
pub(super) struct EdgeInput<'a> {
    pub curve: Id,
    pub geometry: &'a Geometry,
    pub bounds: &'a Bounds,
    pub error: &'a Rational,
    pub same: bool,
}

pub(super) struct State {
    pub plan: crate::chart::Plan,
    pub analytic: BTreeMap<CurveId, Rc<Analytic>>,
    projections: BTreeMap<(SurfaceId, u32), Rc<pcurve::Projection>>,
    pub surfaces: BTreeMap<SurfaceId, Rc<pcurve::Surface>>,
}

impl State {
    pub(super) fn new(plan: crate::chart::Plan) -> Self {
        Self {
            plan,
            analytic: BTreeMap::new(),
            surfaces: BTreeMap::new(),
            projections: BTreeMap::new(),
        }
    }
}

struct Association {
    fit: Option<uv::Fit>,
    id: Id,
    representation: Id,
    surface: Id,
    same: bool,
    seam_side: Option<Side>,
}

struct Wrapper {
    id: Id,
    curve: Id,
    coedges: [Option<CoedgeId>; 2],
}
struct UvVertex {
    point: [Rational; 2],
    id: Id,
}

impl UvVertex {
    fn clone(&self) -> Self {
        Self {
            point: self.point.clone(),
            id: self.id,
        }
    }
}

pub(super) struct ChartTopology {
    faces: BTreeSet<FaceId>,
    associations: BTreeMap<CoedgeId, Association>,
    wrappers: Vec<Wrapper>,
    vertices: BTreeMap<(FaceId, VertexId, Side), UvVertex>,
    context: Option<Id>,
}

impl ChartTopology {
    pub(super) fn new(faces: impl Iterator<Item = FaceId>) -> Self {
        Self {
            faces: faces.collect(),
            associations: BTreeMap::new(),
            wrappers: Vec::new(),
            vertices: BTreeMap::new(),
            context: None,
        }
    }
}

impl Exporter<'_> {
    pub(super) fn topology(&self, shells: &[ShellId]) -> super::Topology {
        let chart = self.chart_state.as_ref().and_then(|state| {
            let topology = ChartTopology::new(
                shells
                    .iter()
                    .flat_map(|shell| self.model.shell(*shell).faces.iter().copied())
                    .filter(|face| state.plan.faces.contains_key(face)),
            );
            (!topology.faces.is_empty()).then_some(topology)
        });
        super::Topology {
            chart,
            ..super::Topology::default()
        }
    }

    pub(super) fn has_chart_edge(&self, edge: EdgeId, topology: &super::Topology) -> bool {
        let (Some(state), Some(topology)) = (&self.chart_state, &topology.chart) else {
            return false;
        };
        state
            .plan
            .uses
            .get(&edge)
            .is_some_and(|uses| uses.iter().any(|(face, _)| topology.faces.contains(face)))
    }

    fn projection(
        &mut self,
        surface: SurfaceId,
        precision: u32,
    ) -> Result<Rc<pcurve::Projection>, Refusal> {
        let state = self.chart_state.as_ref().expect("chart projection");
        if let Some(projection) = state.projections.get(&(surface, precision)) {
            return Ok(projection.clone());
        }
        self.surface(surface)?;
        let written = self
            .chart_state
            .as_ref()
            .expect("chart projection")
            .surfaces[&surface]
            .clone();
        let projection = Rc::new(pcurve::Projection::new(
            written,
            self.model.geom().surface(surface),
            precision,
            self.digits(),
            self.dec,
        )?);
        self.chart_state
            .as_mut()
            .unwrap()
            .projections
            .insert((surface, precision), projection.clone());
        Ok(projection)
    }

    fn uv_vertex(
        &mut self,
        face: FaceId,
        vertex: VertexId,
        side: Side,
        projection: &pcurve::Projection,
        topology: &mut super::Topology,
    ) -> Result<UvVertex, Refusal> {
        let key = (face, vertex, side);
        if let Some(vertex) = topology.chart.as_ref().unwrap().vertices.get(&key) {
            return Ok(vertex.clone());
        }
        let written = self.vertex(vertex, topology)?;
        let point = projection.vertex(
            &self.model.geom().position(self.model.vertex(vertex).point),
            &written.point.xyz.each_ref().map(|n| n.value.clone()),
            side,
            self.dec,
        )?;
        let coordinates = point
            .each_ref()
            .map(|n| number::rational(n, self.digits()).text);
        let id = self.graph.add(
            "CARTESIAN_POINT",
            vec![
                text("''"),
                text(format!("({},{})", coordinates[0], coordinates[1])),
            ],
        );
        let vertex = UvVertex { point, id };
        topology
            .chart
            .as_mut()
            .unwrap()
            .vertices
            .insert(key, vertex.clone());
        Ok(vertex)
    }

    fn uv_representation(&mut self, fit: &uv::Fit, ends: [Id; 2], context: Id) -> Id {
        let mut points = Vec::with_capacity(fit.points.len());
        for (i, point) in fit.points.iter().enumerate() {
            let id = if i == 0 {
                ends[0]
            } else if i + 1 == fit.points.len() {
                ends[1]
            } else {
                self.graph.add(
                    "CARTESIAN_POINT",
                    vec![
                        text("''"),
                        text(format!("({},{})", point[0].text, point[1].text)),
                    ],
                )
            };
            points.push(id);
        }
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
                text(fit::DEGREE.to_string()),
                refs(points),
                text(".UNSPECIFIED."),
                boolean(ends[0] == ends[1]),
                boolean(false),
                text(format!("({multiplicities})")),
                text(format!("({knots})")),
                text(".UNSPECIFIED."),
            ],
        );
        self.graph.add(
            "DEFINITIONAL_REPRESENTATION",
            vec![text("''"), refs([curve]), context.into()],
        )
    }

    pub(super) fn chart_edge(
        &mut self,
        id: EdgeId,
        input: EdgeInput<'_>,
        topology: &mut super::Topology,
    ) -> Result<Id, Refusal> {
        let state = self.chart_state.as_ref().expect("chart edge state");
        let active = &topology.chart.as_ref().expect("chart edge topology").faces;
        let mut uses = [None; 2];
        let mut count = 0;
        for &(face, coedge) in &state.plan.uses[&id] {
            if !active.contains(&face) {
                continue;
            }
            if count == 2 {
                return Err(Refusal::new(
                    "step.not-solid",
                    "chart edge has more than two face uses",
                ));
            }
            uses[count] = Some((face, coedge));
            count += 1;
        }
        assert!(count != 0, "chart edge has no active face");
        let source = self.model.edge(id);
        let (Some(start), Some(end)) = (source.start, source.end) else {
            return Err(unsupported());
        };
        let vertices = if input.same {
            [start, end]
        } else {
            [end, start]
        };
        let source_natural = source.sense == Sense::Same;
        let exact = vertices.map(|v| self.model.geom().position(self.model.vertex(v).point));
        let written = [
            self.vertex(vertices[0], topology)?
                .point
                .xyz
                .each_ref()
                .map(|n| n.value.clone()),
            self.vertex(vertices[1], topology)?
                .point
                .xyz
                .each_ref()
                .map(|n| n.value.clone()),
        ];
        let curves = match input.geometry {
            Geometry::Analytic(curve) => std::slice::from_ref(&curve.written),
            Geometry::Spline { curves, .. } => curves.as_slice(),
        };
        let mut reservation = Some(uv::Budget::new(curves.len(), self.digits(), self.dec)?);
        let precision = reservation.as_ref().unwrap().precision;
        let first = uses[0].unwrap().0;
        let first_projection = self.projection(self.model.face(first).surface, precision)?;
        let analytic_range = match input.geometry {
            Geometry::Analytic(curve) => Some(curve.written.analytic_range(
                [&written[0], &written[1]],
                vertices[0] == vertices[1],
                self.digits(),
                &first_projection.pi,
                self.dec,
            )?),
            Geometry::Spline { .. } => None,
        };
        let bounds: Bounds = std::array::from_fn(|i| {
            [
                &input.bounds[i][0] - input.error,
                &input.bounds[i][1] + input.error,
            ]
        });
        let context = if let Some(context) = topology.chart.as_ref().unwrap().context {
            context
        } else {
            let context = self.parameter_context();
            topology.chart.as_mut().unwrap().context = Some(context);
            context
        };
        let mut coedges = [None; 2];
        let mut pcurves = [None; 2];
        for (index, use_) in uses.into_iter().take(count).enumerate() {
            let (face, coedge) = use_.unwrap();
            let definition = &self.chart_state.as_ref().unwrap().plan.faces[&face];
            let mut sides = definition.sides[&coedge];
            let seam = definition.seams.contains(&id);
            if input.same != source_natural {
                sides.reverse();
            }
            let surface_id = self.model.face(face).surface;
            let projection = if index == 0 {
                first_projection.clone()
            } else {
                self.projection(surface_id, precision)?
            };
            let start = self.uv_vertex(face, vertices[0], sides[0], &projection, topology)?;
            let end = self.uv_vertex(face, vertices[1], sides[1], &projection, topology)?;
            let mut pieces = Vec::with_capacity(curves.len());
            for (piece, curve) in curves.iter().enumerate() {
                let range = match curve {
                    pcurve::Curve::Polynomial { range, .. } => range.clone(),
                    _ => analytic_range.as_ref().unwrap().clone(),
                };
                let guide = if seam {
                    if sides[0] != sides[1] || sides[0] == Side::Interior {
                        return Err(unsupported());
                    }
                    pcurve::Guide::Fixed(if sides[0] == Side::Low {
                        catallax_numbers::Ball::from_rational(&Rational::zero(), precision)
                    } else {
                        &projection.pi
                            * &catallax_numbers::Ball::from_rational(
                                &Rational::from_i64(2),
                                precision,
                            )
                    })
                } else {
                    match input.geometry {
                        Geometry::Analytic(analytic) => projection.analytic_guide(
                            &analytic.exact,
                            [&exact[0], &exact[1]],
                            &range,
                            self.dec,
                        )?,
                        Geometry::Spline { trace, spans, .. } => projection
                            .trace_guide(&trace.pieces()[spans[piece].piece], &spans[piece]),
                    }
                };
                pieces.push(pcurve::Piece {
                    curve: curve.local(&projection.written.placement, precision),
                    guide,
                    projection: &projection,
                    range,
                    endpoints: [
                        (piece == 0).then(|| start.point.clone()),
                        (piece + 1 == curves.len()).then(|| end.point.clone()),
                    ],
                });
            }
            let tolerance = &self.options.curve_tolerance
                / &(&projection.metric(&bounds, None) * &Rational::from_i64(16));
            let budget = if let Some(budget) = reservation.take() {
                budget
            } else {
                uv::Budget::new(curves.len(), self.digits(), self.dec)?
            };
            assert_eq!(
                budget.precision, precision,
                "changed chart working precision"
            );
            let fit = uv::fit(&pieces, self.digits(), &tolerance, budget, self.dec)?;
            let geometry = &self.surfaces[&surface_id];
            let extent = (0..3).fold(Rational::zero(), |sum, i| {
                sum + (&bounds[i][0] - &geometry.origin[i][1])
                    .abs()
                    .max((&bounds[i][1] - &geometry.origin[i][0]).abs())
            });
            let error = &(input.error * &Rational::from_i64(2))
                + &geometry.base_error
                + &geometry.extent_error * &extent
                + projection.metric(&bounds, Some(&fit)) * fit.error.clone()
                + self.face_errors.get(&face).cloned().unwrap_or_default();
            self.deviation(error);
            let surface = self.surfaces[&surface_id].id;
            let representation = self.uv_representation(&fit, [start.id, end.id], context);
            let pcurve = self.graph.add(
                "PCURVE",
                vec![text("''"), surface.into(), representation.into()],
            );
            topology.chart.as_mut().unwrap().associations.insert(
                coedge,
                Association {
                    fit: Some(fit),
                    id: pcurve,
                    representation,
                    surface,
                    same: input.same,
                    seam_side: seam.then_some(sides[0]),
                },
            );
            coedges[index] = Some(coedge);
            pcurves[index] = Some(pcurve);
        }
        let wrapper = self.graph.add(
            "SURFACE_CURVE",
            vec![
                text("''"),
                input.curve.into(),
                refs(pcurves.into_iter().flatten()),
                text(".CURVE_3D."),
            ],
        );
        topology.chart.as_mut().unwrap().wrappers.push(Wrapper {
            id: wrapper,
            curve: input.curve,
            coedges,
        });
        Ok(wrapper)
    }

    pub(super) fn finish_chart_edges(&mut self, topology: &super::Topology) -> Result<(), Refusal> {
        let Some(topology) = &topology.chart else {
            return Ok(());
        };
        for wrapper in &topology.wrappers {
            let first = &topology.associations[&wrapper.coedges[0].unwrap()];
            let second = wrapper.coedges[1].map(|id| &topology.associations[&id]);
            let seam = second.is_some_and(|second| first.surface == second.surface);
            let mut ids = [Some(first.id), second.map(|s| s.id)];
            if seam {
                let second = second.unwrap();
                if first.seam_side == Some(Side::Low) && second.seam_side == Some(Side::High) {
                    ids.swap(0, 1);
                } else if first.seam_side != Some(Side::High) || second.seam_side != Some(Side::Low)
                {
                    return Err(unsupported());
                }
            }
            self.graph.replace(
                wrapper.id,
                if seam { "SEAM_CURVE" } else { "SURFACE_CURVE" },
                vec![
                    text("''"),
                    wrapper.curve.into(),
                    refs(ids.into_iter().flatten()),
                    text(".CURVE_3D."),
                ],
            );
        }
        Ok(())
    }

    pub(super) fn chart_faces(
        &mut self,
        id: FaceId,
        reverse: bool,
        topology: &mut super::Topology,
        faces: &mut Vec<Id>,
    ) -> Result<(), Refusal> {
        let definition = self.model.face(id).clone();
        let mut wires = Vec::with_capacity(definition.loops.len());
        let mut loop_entities = Vec::with_capacity(definition.loops.len());
        let mut boundary_extent = None;
        for loop_id in &definition.loops {
            let Loop::Coedges(source) = self.model.loop_(*loop_id) else {
                return Err(unsupported());
            };
            let mut wire: Vec<_> = source
                .iter()
                .map(|coedge| (*coedge, self.model.coedge(*coedge).sense == Sense::Same))
                .collect();
            if definition.sense == Sense::Reversed {
                wire.reverse();
                for (_, sense) in &mut wire {
                    *sense = !*sense;
                }
            }
            let mut oriented = Vec::with_capacity(wire.len());
            for &(coedge, sense) in &wire {
                let edge = self.edge(self.model.coedge(coedge).edge, topology)?;
                include_bounds(&mut boundary_extent, &edge.bounds);
                oriented.push(self.oriented(edge.id, sense));
            }
            loop_entities.push(
                self.graph
                    .add("EDGE_LOOP", vec![text("''"), refs(oriented)]),
            );
            wires.push(wire);
        }
        let groups = {
            let topology = topology.chart.as_ref().unwrap();
            let curves: Vec<Vec<_>> = wires
                .iter()
                .map(|wire| {
                    wire.iter()
                        .map(|(coedge, sense)| {
                            let association = &topology.associations[coedge];
                            crate::uv_loop::Curve {
                                fit: association.fit.as_ref().expect("pending chart face"),
                                forward: *sense == association.same,
                            }
                        })
                        .collect()
                })
                .collect();
            crate::uv_loop::groups(&curves, self.dec)?
        };
        self.surface(definition.surface)?;
        let geometry = &self.surfaces[&definition.surface];
        let extent = boundary_extent.map_or_else(Rational::zero, |bounds| {
            (0..3).fold(Rational::zero(), |sum, i| {
                sum + (&bounds[i][0] - &geometry.origin[i][1])
                    .abs()
                    .max((&bounds[i][1] - &geometry.origin[i][0]).abs())
            })
        });
        let error = &geometry.base_error
            + &(&geometry.extent_error * &extent)
            + self.face_errors.get(&id).cloned().unwrap_or_default();
        let original = geometry.id;
        self.deviation(error);
        let same = (definition.sense == Sense::Same) != reverse;
        for (index, group) in groups.into_iter().enumerate() {
            let surface = if index == 0 {
                original
            } else {
                self.graph.duplicate(original)
            };
            let mut bounds = Vec::with_capacity(group.holes.len() + 1);
            for loop_index in std::iter::once(group.outer).chain(group.holes) {
                bounds.push(self.graph.add(
                    "FACE_BOUND",
                    vec![text("''"), loop_entities[loop_index].into(), boolean(same)],
                ));
                for &(coedge, _) in &wires[loop_index] {
                    let association = topology
                        .chart
                        .as_mut()
                        .unwrap()
                        .associations
                        .get_mut(&coedge)
                        .unwrap();
                    association.surface = surface;
                    self.graph.replace(
                        association.id,
                        "PCURVE",
                        vec![
                            text("''"),
                            surface.into(),
                            association.representation.into(),
                        ],
                    );
                    association.fit = None;
                }
            }
            faces.push(self.graph.add(
                "ADVANCED_FACE",
                vec![text("''"), refs(bounds), surface.into(), boolean(same)],
            ));
            self.report.faces += 1;
        }
        Ok(())
    }
}
