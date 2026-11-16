//! Numerical through-hole assembly on a collar with two outer bevels.

use std::collections::{BTreeMap, BTreeSet};

use catallax_boolean::numerical::{
    NumericalCylinderSource, NumericalFacePartition, clip_cone_cylinder_face,
};
use catallax_intersect::numerical_native_cylinder_intersection;
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{BSplineCurve, BSplineSurface, Circle, Curve, Cylinder, Point3, Surface};
use catallax_topology::{
    BodyId, EdgeId, FaceId, Loop, Model, Sense,
    numerical::{
        self, NumericalBody, NumericalDefinition, NumericalEdge, NumericalNativeTrim,
        NumericalOptions, NumericalVertex,
    },
};

use super::{
    CylinderSpan, CylinderWallBoundary, append_cylinder_wall, append_planar_cap, cap_curve,
    cylinder_ring_curve,
};

/// A numerical Hole body with the original face identities behind its names.
#[derive(Clone, Debug)]
pub struct NumericalHole {
    body: NumericalBody,
    source_faces: BTreeMap<String, FaceId>,
    tool_faces: BTreeSet<String>,
}
impl NumericalHole {
    /// Checked body, retaining its preset numerical tolerance and native curves.
    pub fn body(&self) -> &NumericalBody {
        &self.body
    }
    /// Output faces descended from each original stock face.
    pub fn source_faces(&self) -> &BTreeMap<String, FaceId> {
        &self.source_faces
    }
    /// New faces supplied by the compiled cylindrical Hole wall.
    pub fn tool_faces(&self) -> &BTreeSet<String> {
        &self.tool_faces
    }
    /// Rebuilds from the original stock and edited Hole controls, retaining history.
    pub fn rebuild(
        &self,
        model: &Model,
        body: BodyId,
        key: &str,
        span: &CylinderSpan,
        options: NumericalOptions,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        guarded(|| {
            prepare(model, body, key, span, options, dec)?
                .map(|prepared| prepared.admit(Some(&self.body), dec))
                .transpose()
        })
    }
    /// Moves the checked body into a numerical result consumer.
    pub fn into_body(self) -> NumericalBody {
        self.body
    }
}
fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("hole.numerical-assembly", detail)
}
fn q(value: i64) -> Rational {
    Rational::from_i64(value)
}
fn flip(sense: Sense) -> Sense {
    if sense == Sense::Same {
        Sense::Reversed
    } else {
        Sense::Same
    }
}
fn at(curve: &BSplineCurve, t: &Rational, dec: &Decider) -> Result<[Rational; 3], Refusal> {
    Ok(curve.evaluate_box([t, t], 128, dec)?.map(|v| v.mid()))
}
fn close(a: &[Rational; 3], b: &[Rational; 3], tol: &Rational) -> bool {
    a.iter()
        .zip(b)
        .fold(q(0), |sum, (a, b)| sum + (a - b).abs())
        <= *tol
}
fn rings(
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<BTreeMap<EdgeId, Circle>, Refusal> {
    let mut result = BTreeMap::new();
    for boundary in &model.face(face).loops {
        let Loop::Coedges(uses) = model.loop_(*boundary) else {
            return Err(bad("collar faces require complete circular boundaries"));
        };
        for use_id in uses {
            dec.charge_work(1)?;
            let id = model.coedge(*use_id).edge;
            let Curve::Circle(circle) = model.geom().curve(model.edge(id).curve) else {
                return Err(bad("collar face has a noncircular original boundary"));
            };
            result.insert(id, circle.clone());
        }
    }
    Ok(result)
}
fn height(cylinder: &Cylinder, point: &Point3, dec: &mut Decider) -> Result<Rational, Refusal> {
    let axis = cylinder
        .axis
        .scaled(&(&Real::one() / &cylinder.axis.norm_squared().sqrt()));
    Ok(dec
        .enclose_budgeted_at(&point.minus(&cylinder.origin).dot(&axis), 128)?
        .mid())
}
fn source(
    model: &Model,
    face: FaceId,
    options: &NumericalOptions,
    dec: &mut Decider,
) -> Result<NumericalCylinderSource, Refusal> {
    let Surface::Cylinder(cylinder) = model.geom().surface(model.face(face).surface) else {
        return Err(bad("expected original cylindrical wall"));
    };
    let mut heights = rings(model, face, dec)?
        .values()
        .map(|ring| height(cylinder, &ring.center, dec))
        .collect::<Result<Vec<_>, _>>()?;
    heights.sort();
    heights.dedup();
    let interval: [Rational; 2] = heights
        .try_into()
        .map_err(|_| bad("cylindrical wall needs two distinct original circular bounds"))?;
    NumericalCylinderSource::new(cylinder, interval, options.clone(), dec)
}
fn ring_curve(
    source: &NumericalCylinderSource,
    row: usize,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let source = source.source();
    BSplineCurve::new(
        source.degrees()[0],
        source.knots()[0].to_vec(),
        (0..9)
            .map(|i| source.controls()[2 * i + row].clone())
            .collect(),
        (0..9)
            .map(|i| source.weights()[2 * i + row].clone())
            .collect(),
        dec,
    )
}
fn only_use(
    part: &NumericalFacePartition,
    edge: usize,
    dec: &mut Decider,
) -> Result<Sense, Refusal> {
    let mut found = None;
    for trim in part
        .definition
        .faces
        .iter()
        .flat_map(|face| face.loops.iter().flatten())
    {
        dec.charge_work(1)?;
        if trim.edge == edge && found.replace(trim.sense).is_some() {
            return Err(bad("shared boundary has several retained cone uses"));
        }
    }
    found.ok_or_else(|| bad("shared boundary has no retained cone use"))
}
struct StockBoundaries {
    curves: Vec<CylinderWallBoundary>,
    anchors: BTreeMap<usize, Vec<[Rational; 2]>>,
}

struct Part {
    part: NumericalFacePartition,
    vertices: usize,
    edges: usize,
}
struct Assembly<'a> {
    model: &'a Model,
    span: &'a CylinderSpan,
    options: NumericalOptions,
    definition: NumericalDefinition,
    native: Vec<NumericalNativeTrim>,
    parts: Vec<Part>,
    tool: NumericalCylinderSource,
    tool_boundaries: Vec<CylinderWallBoundary>,
    source_faces: BTreeMap<String, FaceId>,
}
impl Assembly<'_> {
    fn cap(
        &mut self,
        face: FaceId,
        bore: &NumericalCylinderSource,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let key = format!("source/{face:?}");
        let Surface::Plane(plane) = self.model.geom().surface(self.model.face(face).surface) else {
            return Err(bad("expected original cap plane"));
        };
        let original = rings(self.model, face, dec)?;
        let mut radii = original
            .values()
            .map(|ring| {
                dec.enclose_budgeted_at(&ring.radius_sq, 128)
                    .map(|v| v.mid())
            })
            .collect::<Result<Vec<_>, _>>()?;
        radii.sort();
        let [inner_radius, outer_radius] = radii.as_slice() else {
            return Err(bad("cap needs an original annular pair"));
        };
        let center = &original
            .values()
            .next()
            .ok_or_else(|| bad("missing cap circle"))?
            .center;
        let mut outer = Vec::new();
        let mut known = BTreeMap::<Rational, usize>::new();
        for part in &self.parts {
            let mut vertices = BTreeSet::new();
            for link in &part.part.ring_links {
                dec.charge_work(1)?;
                if original.contains_key(&part.part.rings[link.ring].edge) {
                    outer.push((
                        part.edges + link.edge,
                        flip(only_use(&part.part, link.edge, dec)?),
                    ));
                    vertices.extend(part.part.definition.edges[link.edge].vertices);
                }
            }
            for event in &part.part.events {
                dec.charge_work(1)?;
                if let Some(vertex) = event.vertex.filter(|vertex| vertices.contains(vertex)) {
                    for hit in &event.hits {
                        let u = if hit.uv[1][0] == q(4) {
                            q(0)
                        } else {
                            hit.uv[1][0].clone()
                        };
                        if let Some(previous) = known.insert(u, vertex + part.vertices)
                            && previous != vertex + part.vertices
                        {
                            return Err(bad("cap event has conflicting original-ring owners"));
                        }
                    }
                }
            }
        }
        if outer.is_empty() {
            return Err(bad("cap has no retained original outer ring"));
        }
        let normal = plane.normal.scaled(&Real::from_i64(
            if self.model.face(face).sense == Sense::Same {
                1
            } else {
                -1
            },
        ));
        let cap_sense = if self.span.cylinder.axis.dot(&normal).to_f64() > 0. {
            Sense::Reversed
        } else {
            Sense::Same
        };
        for quarter in self.tool.quarters() {
            let Surface::BSpline(support) = quarter.surface() else {
                return Err(bad("missing retained tool quarter"));
            };
            let cap = cap_curve(
                &format!("{key}/tool-quarter-{}", quarter.id()),
                self.span,
                support.carrier(),
                plane,
                &self.options,
                dec,
            )?
            .ok_or_else(|| bad("tool axis does not meet cap plane"))?;
            let [a, b] = cap.curve().domain();
            let mut nodes: BTreeMap<_, _> = known
                .iter()
                .filter(|(u, _)| *u >= a && *u <= b)
                .map(|(u, v)| (u.clone(), Some(*v)))
                .collect();
            nodes.entry(a.clone()).or_insert(None);
            nodes.entry(b.clone()).or_insert(None);
            let nodes: Vec<_> = nodes.into_iter().collect();
            for (segment, pair) in nodes.windows(2).enumerate() {
                let midpoint = (&pair[0].0 + &pair[1].0) / q(2);
                let point = at(cap.curve(), &midpoint, dec)?;
                let [x, y, z] = point.map(Real::from);
                let radial = Point3::new(x, y, z).minus(center).norm_squared();
                let radial = dec.enclose_budgeted_at(&radial, 128)?;
                if radial.hi() <= *inner_radius || radial.lo() >= *outer_radius {
                    continue;
                }
                if radial.lo() <= *inner_radius || radial.hi() >= *outer_radius {
                    return Err(bad(
                        "cap material sample is unresolved at its original ring",
                    ));
                }
                let mut vertices = [0; 2];
                for end in 0..2 {
                    vertices[end] = if let Some(vertex) = pair[end].1 {
                        vertex
                    } else {
                        let u = if pair[end].0 == q(4) {
                            q(0)
                        } else {
                            pair[end].0.clone()
                        };
                        if let Some(vertex) = known.get(&u) {
                            *vertex
                        } else {
                            let vertex = self.definition.vertices.len();
                            self.definition.vertices.push(NumericalVertex {
                                key: format!("{key}/tool-seam-{u}"),
                                position: at(cap.curve(), &pair[end].0, dec)?,
                            });
                            known.insert(u, vertex);
                            vertex
                        }
                    };
                }
                let edge = self.definition.edges.len();
                self.definition.edges.push(NumericalEdge {
                    key: format!("{key}/tool-quarter-{}/segment-{segment}", quarter.id()),
                    curve: Curve::BSpline(cap.curve().clone()),
                    vertices,
                    parameters: [pair[0].0.clone(), pair[1].0.clone()],
                });
                outer.push((edge, cap_sense));
                self.tool_boundaries.push(CylinderWallBoundary {
                    edge,
                    quarter: quarter.id(),
                    sense: flip(cap_sense),
                    curve: cap.cylinder_uv().clone(),
                });
            }
        }
        let signed = height(bore.original(), center, dec)?;
        let row = if (&signed - &bore.axis_parameters_mm()[0]).abs() <= self.options.tolerance_mm {
            0
        } else if (&signed - &bore.axis_parameters_mm()[1]).abs() <= self.options.tolerance_mm {
            1
        } else {
            return Err(bad("cap does not meet an original bore ring"));
        };
        let circle = ring_curve(bore, row, dec)?;
        let first = self.definition.vertices.len();
        for quarter in 0..4 {
            self.definition.vertices.push(NumericalVertex {
                key: format!("{key}/bore-vertex-{quarter}"),
                position: at(&circle, &q(quarter), dec)?,
            });
        }
        let inner_sense = if bore.original().axis.dot(&normal).to_f64() > 0. {
            Sense::Reversed
        } else {
            Sense::Same
        };
        let mut inner = Vec::new();
        for quarter in 0..4 {
            let edge = self.definition.edges.len();
            self.definition.edges.push(NumericalEdge {
                key: format!("{key}/bore-edge-{quarter}"),
                curve: Curve::BSpline(circle.clone()),
                vertices: [first + quarter, first + (quarter + 1) % 4],
                parameters: [q(quarter as i64), q(quarter as i64 + 1)],
            });
            inner.push((edge, inner_sense));
        }
        append_planar_cap(
            &mut self.definition,
            &mut self.native,
            &key,
            plane,
            self.model.face(face).sense,
            &[outer, inner],
            dec,
        )?;
        self.source_faces.insert(key, face);
        Ok(())
    }
    fn stock_rings(
        &self,
        face: FaceId,
        source: &NumericalCylinderSource,
        dec: &mut Decider,
    ) -> Result<StockBoundaries, Refusal> {
        let original = rings(self.model, face, dec)?;
        let mut boundaries = Vec::new();
        let mut anchors = BTreeMap::<usize, Vec<[Rational; 2]>>::new();
        for part in &self.parts {
            let mut vertices = BTreeSet::new();
            for link in &part.part.ring_links {
                dec.charge_work(1)?;
                if !original.contains_key(&part.part.rings[link.ring].edge) {
                    continue;
                }
                let ring = &part.part.rings[link.ring].circle;
                let signed = height(source.original(), &ring.center, dec)?;
                let row = if (&signed - &source.axis_parameters_mm()[0]).abs()
                    <= self.options.tolerance_mm
                {
                    0
                } else if (&signed - &source.axis_parameters_mm()[1]).abs()
                    <= self.options.tolerance_mm
                {
                    1
                } else {
                    return Err(bad(
                        "retained stock ring does not meet either finite cylinder row",
                    ));
                };
                let edge = &part.part.definition.edges[link.edge];
                let Curve::BSpline(curve) = &edge.curve else {
                    return Err(bad("retained ring lost its original curve"));
                };
                let map = cylinder_ring_curve(
                    &edge.key,
                    curve,
                    edge.parameters.clone(),
                    source,
                    row,
                    &self.options,
                    dec,
                )?
                .ok_or_else(|| bad("original ring has no matching native cylinder quarter"))?;
                boundaries.push(CylinderWallBoundary {
                    edge: part.edges + link.edge,
                    quarter: map.quarter(),
                    sense: flip(only_use(&part.part, link.edge, dec)?),
                    curve: map.curve().clone(),
                });
                vertices.extend(edge.vertices);
            }
            for event in &part.part.events {
                dec.charge_work(1)?;
                if let Some(vertex) = event.vertex.filter(|vertex| vertices.contains(vertex)) {
                    anchors
                        .entry(vertex + part.vertices)
                        .or_default()
                        .extend(event.hits.iter().map(|hit| hit.uv[1].clone()));
                }
            }
        }
        Ok(StockBoundaries {
            curves: boundaries,
            anchors,
        })
    }
    fn stock_wall(
        &mut self,
        face: FaceId,
        stock: &NumericalCylinderSource,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let key = format!("source/{face:?}");
        let StockBoundaries {
            curves: mut boundaries,
            anchors,
        } = self.stock_rings(face, stock, dec)?;
        let incoming = numerical_native_cylinder_intersection(
            [stock.original(), &self.span.cylinder],
            [
                stock.axis_parameters_mm().clone(),
                self.span.interval_mm.clone(),
            ],
            dec,
        )?;
        if !incoming.sources()[0]
            .source()
            .same_definition(stock.source())
            || !incoming.sources()[1]
                .source()
                .same_definition(self.tool.source())
        {
            return Err(bad(
                "incoming wall query changed a retained native cylinder definition",
            ));
        }
        let mut seams = BTreeMap::<(usize, usize), Vec<(Rational, usize)>>::new();
        for (index, piece) in incoming.pieces().iter().enumerate() {
            let mut vertices = [0; 2];
            for (end, t) in piece.curve().domain().into_iter().enumerate() {
                let position = at(piece.curve(), t, dec)?;
                let uv = piece
                    .uv()
                    .each_ref()
                    .map(|curve| at(curve, t, dec))
                    .into_iter()
                    .collect::<Result<Vec<_>, _>>()?;
                let mut candidates = Vec::new();
                for (vertex, maps) in &anchors {
                    dec.charge_work(1)?;
                    if close(
                        &position,
                        &self.definition.vertices[*vertex].position,
                        &self.options.tolerance_mm,
                    ) && maps.iter().any(|map| {
                        (&map[0] - &uv[1][0]).abs() + (&map[1] - &uv[1][1]).abs()
                            <= self.options.tolerance_mm
                    }) {
                        candidates.push(*vertex)
                    }
                }
                vertices[end] = if let [vertex] = candidates.as_slice() {
                    *vertex
                } else {
                    if !candidates.is_empty() {
                        return Err(bad(
                            "incoming wall endpoint has several original-ring owners",
                        ));
                    }
                    let mut owner = None;
                    for (side, at) in uv.iter().enumerate() {
                        for seam in 0..=4 {
                            if (&at[0] - &q(seam)).abs() <= self.options.tolerance_mm {
                                owner = Some((side, seam as usize % 4, at[1].clone()));
                                break;
                            }
                        }
                        if owner.is_some() {
                            break;
                        }
                    }
                    let (side, seam, v) = owner.ok_or_else(|| {
                        bad("wall endpoint has neither original-ring nor requested-seam ownership")
                    })?;
                    let known = seams.entry((side, seam)).or_default();
                    if let Some((_, vertex)) = known
                        .iter()
                        .find(|(old, _)| (old - &v).abs() <= self.options.tolerance_mm)
                    {
                        if !close(
                            &position,
                            &self.definition.vertices[*vertex].position,
                            &self.options.tolerance_mm,
                        ) {
                            return Err(bad("owned cylinder seam endpoints exceed tolerance"));
                        }
                        *vertex
                    } else {
                        let vertex = self.definition.vertices.len();
                        self.definition.vertices.push(NumericalVertex {
                            key: format!("{key}/seam-{side}-{seam}/event-{}", known.len()),
                            position,
                        });
                        known.push((v, vertex));
                        vertex
                    }
                };
            }
            let edge = self.definition.edges.len();
            self.definition.edges.push(NumericalEdge {
                key: format!("{key}/fit-{index}"),
                curve: Curve::BSpline(piece.curve().clone()),
                vertices,
                parameters: piece.curve().domain().map(Clone::clone),
            });
            let forward = outside_left(stock.source(), &piece.uv()[0], &self.span.cylinder, dec)?;
            let direction = if forward {
                Sense::Same
            } else {
                Sense::Reversed
            };
            let sense = if self.model.face(face).sense == Sense::Same {
                direction
            } else {
                flip(direction)
            };
            boundaries.push(CylinderWallBoundary {
                edge,
                quarter: piece.source_domains()[0][0][0].to_f64_nearest() as usize,
                sense,
                curve: piece.uv()[0].clone(),
            });
            self.tool_boundaries.push(CylinderWallBoundary {
                edge,
                quarter: piece.source_domains()[1][0][0].to_f64_nearest() as usize,
                sense: flip(sense),
                curve: piece.uv()[1].clone(),
            });
        }
        let faces = append_cylinder_wall(
            &mut self.definition,
            &mut self.native,
            &key,
            stock,
            self.model.face(face).sense,
            &boundaries,
            dec,
        )?;
        for index in faces {
            self.source_faces
                .insert(self.definition.faces[index].key.clone(), face);
        }
        Ok(())
    }
    fn bore_wall(
        &mut self,
        face: FaceId,
        bore: &NumericalCylinderSource,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let key = format!("source/{face:?}");
        let mut boundaries = Vec::new();
        for cap in &self.definition.faces {
            if !matches!(cap.surface, Surface::Plane(_)) {
                continue;
            }
            for trim in &cap.loops[1] {
                let edge = &self.definition.edges[trim.edge];
                let Curve::BSpline(curve) = &edge.curve else {
                    return Err(bad("bore ring lost its retained curve"));
                };
                let mapping = (0..2)
                    .map(|row| {
                        cylinder_ring_curve(
                            &edge.key,
                            curve,
                            edge.parameters.clone(),
                            bore,
                            row,
                            &self.options,
                            dec,
                        )
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let mut matches = mapping.into_iter().flatten();
                let map = matches
                    .next()
                    .ok_or_else(|| bad("cap bore ring has no native cylinder row"))?;
                if matches.next().is_some() {
                    return Err(bad("cap bore ring matches both finite rows"));
                }
                boundaries.push(CylinderWallBoundary {
                    edge: trim.edge,
                    quarter: map.quarter(),
                    sense: flip(trim.sense),
                    curve: map.curve().clone(),
                });
            }
        }
        let faces = append_cylinder_wall(
            &mut self.definition,
            &mut self.native,
            &key,
            bore,
            self.model.face(face).sense,
            &boundaries,
            dec,
        )?;
        for index in faces {
            self.source_faces
                .insert(self.definition.faces[index].key.clone(), face);
        }
        Ok(())
    }
}

fn outside_left(
    source: &BSplineSurface,
    uv: &BSplineCurve,
    tool: &Cylinder,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let axis = tool
        .axis
        .scaled(&(&Real::one() / &tool.axis.norm_squared().sqrt()));
    let enclose = |value: &Real, dec: &mut Decider| dec.enclose_budgeted_at(value, 128);
    let origin = [
        enclose(&tool.origin.x, dec)?,
        enclose(&tool.origin.y, dec)?,
        enclose(&tool.origin.z, dec)?,
    ];
    let axis = [
        enclose(&axis.x, dec)?,
        enclose(&axis.y, dec)?,
        enclose(&axis.z, dec)?,
    ];
    let mut direction = None;
    for (_, span) in uv.spans() {
        let t = (span[0] + span[1]) / q(2);
        let at = at(uv, &t, dec)?;
        let jet = source.jet_box([&at[0], &at[0]], [&at[1], &at[1]], 128, dec)?;
        let d: [Ball; 3] = std::array::from_fn(|i| &jet.value[i] - &origin[i]);
        let h = &d[0] * &axis[0] + &d[1] * &axis[1] + &d[2] * &axis[2];
        let radial: [Ball; 3] = std::array::from_fn(|i| &d[i] - &(&h * &axis[i]));
        let gradient = jet
            .first
            .each_ref()
            .map(|d| &radial[0] * &d[0] + &radial[1] * &d[1] + &radial[2] * &d[2]);
        let tangent = uv.first_at(&Real::from(t), dec)?;
        let dx = enclose(&tangent.x, dec)?;
        let dy = enclose(&tangent.y, dec)?;
        let left = &gradient[1] * &dx - &gradient[0] * &dy;
        let forward = if left.lo().sign() == Sign::Positive {
            true
        } else if left.hi().sign() == Sign::Negative {
            false
        } else {
            return Err(bad("wall material-side sample is tangent or unresolved"));
        };
        if direction.is_some_and(|previous| previous != forward) {
            return Err(bad(
                "wall material side changes within the retained native piece",
            ));
        }
        direction = Some(forward);
    }
    direction.ok_or_else(|| bad("wall piece has no nonzero native span"))
}

/// Builds a separate numerical Hole result on a two-bevel annular cylinder.
/// The compiled tool, all reported cone/wall pieces, original rings and native
/// parameter definitions are retained. The original exact model is immutable.
/// `None` denotes another input family. An intersected bore is not handled by
/// this continuation. Supported exact Hole cases keep using the exact evaluator.
pub fn build_numerical_hole(
    model: &Model,
    body: BodyId,
    key: &str,
    span: &CylinderSpan,
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<Option<NumericalHole>, Refusal> {
    guarded(|| {
        prepare(model, body, key, span, options, dec)?
            .map(|prepared| prepared.admit(None, dec))
            .transpose()
    })
}

fn guarded<T>(run: impl FnOnce() -> Result<T, Refusal>) -> Result<T, Refusal> {
    if catallax_numbers::operation_budget_active() {
        return run();
    }
    let (result, work) = catallax_numbers::with_operation_budget(numerical::DEFAULT_WORK, run);
    if work.exhausted {
        Err(Refusal::new(
            "numbers.operation-budget",
            "numerical Hole assembly exhausted its operation allowance",
        ))
    } else {
        result
    }
}
struct Prepared {
    definition: NumericalDefinition,
    native: Vec<NumericalNativeTrim>,
    options: NumericalOptions,
    source_faces: BTreeMap<String, FaceId>,
    tool_faces: BTreeSet<String>,
}
impl Prepared {
    fn admit(
        self,
        previous: Option<&NumericalBody>,
        dec: &mut Decider,
    ) -> Result<NumericalHole, Refusal> {
        let body = if let Some(previous) = previous {
            previous.rebuild_with_native_trims(self.definition, self.native, self.options, dec)?
        } else {
            numerical::build_with_native_trims(self.definition, self.native, self.options, dec)?
        };
        Ok(NumericalHole {
            body,
            source_faces: self.source_faces,
            tool_faces: self.tool_faces,
        })
    }
}

fn prepare(
    model: &Model,
    body: BodyId,
    key: &str,
    span: &CylinderSpan,
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<Option<Prepared>, Refusal> {
    let faces = model.faces_of_body(body);
    let mut cylinders = Vec::new();
    let mut cones = Vec::new();
    let mut caps = Vec::new();
    for face in faces {
        dec.charge_work(1)?;
        match model.geom().surface(model.face(face).surface) {
            Surface::Cylinder(cylinder) => cylinders.push((
                face,
                dec.enclose_budgeted_at(&cylinder.radius_sq, 128)?.mid(),
            )),
            Surface::Cone(_) => cones.push(face),
            Surface::Plane(_) => caps.push(face),
            _ => return Ok(None),
        }
    }
    if cylinders.len() != 2 || cones.len() != 2 || caps.len() != 2 {
        return Ok(None);
    }
    cylinders.sort_by(|a, b| a.1.cmp(&b.1));
    let bore_face = cylinders[0].0;
    let stock_face = cylinders[1].0;
    let bore = source(model, bore_face, &options, dec)?;
    let stock = source(model, stock_face, &options, dec)?;
    let stock_axis = stock
        .original()
        .axis
        .scaled(&(&Real::one() / &stock.original().axis.norm_squared().sqrt()));
    let bore_axis = bore
        .original()
        .axis
        .scaled(&(&Real::one() / &bore.original().axis.norm_squared().sqrt()));
    let parallel = stock_axis.cross(&bore_axis).norm_squared();
    let centered = bore
        .original()
        .origin
        .minus(&stock.original().origin)
        .cross(&stock_axis)
        .norm_squared();
    let tolerance_squared = &options.tolerance_mm * &options.tolerance_mm;
    if dec.enclose_budgeted_at(&parallel, 128)?.hi() > tolerance_squared
        || dec.enclose_budgeted_at(&centered, 128)?.hi() > tolerance_squared
    {
        return Ok(None);
    }

    let bore_query = numerical_native_cylinder_intersection(
        [bore.original(), &span.cylinder],
        [bore.axis_parameters_mm().clone(), span.interval_mm.clone()],
        dec,
    )?;
    if !bore_query.pieces().is_empty() {
        return Err(bad(
            "this numerical continuation requires a Hole separated from the bore",
        ));
    }
    let tool = NumericalCylinderSource::new(
        &span.cylinder,
        span.interval_mm.clone(),
        options.clone(),
        dec,
    )?;
    let mut assembly = Assembly {
        model,
        span,
        options: options.clone(),
        definition: NumericalDefinition {
            name: key.into(),
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
        },
        native: Vec::new(),
        parts: Vec::new(),
        tool,
        tool_boundaries: Vec::new(),
        source_faces: BTreeMap::new(),
    };
    for face in cones {
        let part = clip_cone_cylinder_face(
            model,
            face,
            &format!("source/{face:?}"),
            &span.cylinder,
            span.interval_mm.clone(),
            options.clone(),
            dec,
        )?
        .ok_or_else(|| bad("cone has no supported finite circular band"))?;
        let vertices = assembly.definition.vertices.len();
        let edges = assembly.definition.edges.len();
        assembly
            .definition
            .vertices
            .extend(part.definition.vertices.clone());
        assembly
            .definition
            .edges
            .extend(part.definition.edges.clone());
        for edge in &mut assembly.definition.edges[edges..] {
            edge.vertices
                .iter_mut()
                .for_each(|vertex| *vertex += vertices)
        }
        for mut retained in part.definition.faces.clone() {
            retained
                .loops
                .iter_mut()
                .flatten()
                .for_each(|trim| trim.edge += edges);
            assembly.source_faces.insert(retained.key.clone(), face);
            assembly.definition.faces.push(retained);
        }
        assembly.native.extend(part.native_trims.clone());
        for trim in &part.cylinder_trims {
            assembly.tool_boundaries.push(CylinderWallBoundary {
                edge: edges + trim.edge,
                quarter: trim.quarter,
                sense: flip(only_use(&part, trim.edge, dec)?),
                curve: trim.curve.clone(),
            });
        }
        assembly.parts.push(Part {
            part,
            vertices,
            edges,
        });
    }
    for face in caps {
        assembly.cap(face, &bore, dec)?;
    }
    assembly.stock_wall(stock_face, &stock, dec)?;
    assembly.bore_wall(bore_face, &bore, dec)?;
    let tool_indices = append_cylinder_wall(
        &mut assembly.definition,
        &mut assembly.native,
        &format!("{key}/wall"),
        &assembly.tool,
        Sense::Reversed,
        &assembly.tool_boundaries,
        dec,
    )?;
    let tool_faces = tool_indices
        .into_iter()
        .map(|index| assembly.definition.faces[index].key.clone())
        .collect();
    Ok(Some(Prepared {
        definition: assembly.definition,
        native: assembly.native,
        options,
        source_faces: assembly.source_faces,
        tool_faces,
    }))
}
