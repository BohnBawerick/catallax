//! Rolling-envelope surfaces and their boundary rows. One rounded control table
//! writes the surface and every boundary curve, so a shared row cannot drift
//! from the surface it bounds. Each boundary keeps the straight parameter-line
//! pcurve that holds exactly on the rows' own parametrisation.

use super::*;
use crate::graph::Field;
use crate::junction::EnvelopeEdge;
use catallax_numbers::Ball;
use catallax_shapes::{BoundaryRow, ParameterTable, RollingEnvelope, SplineSupport};

// Inverse of a parameter table's angular tolerance. At 10^-11 rad, the table
// stays far below a reader's edge tolerance on any chart radius in the file.
const TABLE_RADIANS: i64 = 100_000_000_000;

pub(super) struct Written {
    id: Id,
    points: Vec<Id>,
    xyz: Vec<Point>,
    weights: Vec<String>,
    columns: usize,
    knots: (Field, Field),
    error: Rational,
}

fn knot_fields(knots: &[Rational], digits: u32) -> (Field, Field) {
    let mut multiplicities = Vec::new();
    let mut values = Vec::new();
    let mut start = 0;
    while start < knots.len() {
        let mut end = start + 1;
        while end < knots.len() && knots[end] == knots[start] {
            end += 1;
        }
        multiplicities.push(text((end - start).to_string()));
        values.push(text(number::rational(&knots[start], digits).text));
        start = end;
    }
    (Field::List(multiplicities), Field::List(values))
}

impl Exporter<'_> {
    fn rolling_envelope(&self, id: SurfaceId) -> Result<&RollingEnvelope, Refusal> {
        match self.model.geom().surface(id) {
            Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) => Ok(envelope),
            _ => Err(unsupported()),
        }
    }

    pub(super) fn written_envelope(&mut self, id: SurfaceId) -> Result<Rc<Written>, Refusal> {
        if let Some(written) = self.envelopes.get(&id) {
            return Ok(written.clone());
        }
        let envelope = self.rolling_envelope(id)?.clone();
        let carrier = envelope.carrier();
        let digits = self.digits();
        let mut points = Vec::with_capacity(carrier.controls().len());
        let mut xyz = Vec::with_capacity(carrier.controls().len());
        let mut weights = Vec::with_capacity(carrier.controls().len());
        let mut position = Rational::zero();
        let mut relative = Rational::zero();
        let mut bounds = None;
        for (control, weight) in envelope.enclosed_controls(self.dec)? {
            self.dec.check_deadline()?;
            let point = Point {
                xyz: control
                    .each_ref()
                    .map(|value| number::ball_number(value, digits)),
            };
            let weight = number::ball_number(&weight, digits);
            position = position.max(point.error());
            relative = relative.max(&weight.error / &(&weight.value - &weight.error));
            include_bounds(&mut bounds, &point_bounds(&point));
            points.push(self.point_entity(&point));
            xyz.push(point);
            weights.push(weight.text);
        }
        let bounds: Bounds = bounds.ok_or_else(unsupported)?;
        // Positive weights make each surface point a convex combination of
        // its controls, so a control move adds directly and a relative weight
        // change e moves the combination by at most 2e/(1-e) times the hull
        // diameter, here bounded by the box's L1 size.
        let diameter = (0..3).fold(Rational::zero(), |sum, i| {
            sum + (&bounds[i][1] - &bounds[i][0])
        });
        let error = &position + &(relative.mul_pow2(2) * diameter);
        self.deviation(error.clone());
        let [count_u, count_v] = carrier.control_count();
        let [degree_u, degree_v] = carrier.degrees();
        let [knots_u, knots_v] = carrier.knots();
        let u = knot_fields(knots_u, digits);
        let v = knot_fields(knots_v, digits);
        let grid = |values: Vec<Field>| {
            Field::List(
                values
                    .chunks(count_v)
                    .map(|row| Field::List(row.to_vec()))
                    .collect(),
            )
        };
        let surface = self.graph.add_complex(vec![
            ("BOUNDED_SURFACE", vec![]),
            (
                "B_SPLINE_SURFACE",
                vec![
                    text(degree_u.to_string()),
                    text(degree_v.to_string()),
                    grid(points.iter().map(|p| Field::Ref(*p)).collect()),
                    text(".UNSPECIFIED."),
                    boolean(false),
                    boolean(false),
                    boolean(false),
                ],
            ),
            (
                "B_SPLINE_SURFACE_WITH_KNOTS",
                vec![u.0.clone(), v.0, u.1.clone(), v.1, text(".UNSPECIFIED.")],
            ),
            ("GEOMETRIC_REPRESENTATION_ITEM", vec![]),
            (
                "RATIONAL_B_SPLINE_SURFACE",
                vec![grid(weights.iter().cloned().map(text).collect())],
            ),
            ("REPRESENTATION_ITEM", vec![text("''")]),
            ("SURFACE", vec![]),
        ]);
        debug_assert_eq!(points.len(), count_u * count_v);
        let origin = bounds.clone();
        self.surfaces.insert(
            id,
            SurfaceGeometry {
                id: surface,
                origin,
                base_error: error.clone(),
                extent_error: Rational::zero(),
            },
        );
        let written = Rc::new(Written {
            id: surface,
            points,
            xyz,
            weights,
            columns: count_u,
            knots: u,
            error,
        });
        self.envelopes.insert(id, written.clone());
        Ok(written)
    }

    // The row on its analytic neighbour, as a cubic B-spline in that chart over
    // the row's own parameter. Each piece is one cubic Hermite span of the table.
    fn neighbour_pcurve(
        &mut self,
        plan: &EnvelopeEdge,
        table: &ParameterTable,
        context: Id,
    ) -> Result<(Id, Rational), Refusal> {
        let neighbour = &plan.neighbour;
        self.surface(neighbour.surface)?;
        let surface = self.surfaces[&neighbour.surface].id;
        let digits = self.digits();
        let middle = |ball: &Ball| ball.mid();
        let value = |base: &Ball, sign: i64, data: &[Rational; 2]| {
            [
                middle(base) + Rational::from_i64(sign) * data[0].clone(),
                Rational::from_i64(sign) * data[1].clone(),
            ]
        };
        let nodes: Vec<[[Rational; 2]; 2]> = (0..table.knots.len())
            .map(|i| {
                let angle = &table.angle[i];
                match plan.row {
                    BoundaryRow::Rail => [
                        value(&neighbour.u, neighbour.u_sign, angle),
                        value(&neighbour.v, neighbour.v_sign, &table.height[i]),
                    ],
                    BoundaryRow::Mouth => [
                        value(&neighbour.u, neighbour.u_sign, angle),
                        [middle(&neighbour.v), Rational::zero()],
                    ],
                    _ => [
                        [middle(&neighbour.u), Rational::zero()],
                        value(&neighbour.v, neighbour.v_sign, angle),
                    ],
                }
            })
            .collect();
        let mut controls = Vec::with_capacity(3 * nodes.len());
        let mut rounding = &neighbour.u.rad() + &neighbour.v.rad();
        let mut point = |graph: &mut Graph, x: Rational, y: Rational| {
            let (x, y) = (number::rational(&x, digits), number::rational(&y, digits));
            rounding = rounding.clone().max(&x.error + &y.error);
            graph.add(
                "CARTESIAN_POINT",
                vec![text("''"), text(format!("({},{})", x.text, y.text))],
            )
        };
        controls.push(point(
            &mut self.graph,
            nodes[0][0][0].clone(),
            nodes[0][1][0].clone(),
        ));
        let third = Rational::from_ratio(1, 3).expect("nonzero");
        for i in 1..nodes.len() {
            let width = (&table.knots[i] - &table.knots[i - 1]) * third.clone();
            let [[u0, du0], [v0, dv0]] = &nodes[i - 1];
            let [[u1, du1], [v1, dv1]] = &nodes[i];
            controls.push(point(
                &mut self.graph,
                u0 + &(&width * du0),
                v0 + &(&width * dv0),
            ));
            controls.push(point(
                &mut self.graph,
                u1 - &(&width * du1),
                v1 - &(&width * dv1),
            ));
            controls.push(point(&mut self.graph, u1.clone(), v1.clone()));
        }
        let last = table.knots.len() - 1;
        let multiplicities: Vec<_> = (0..=last)
            .map(|i| text(if i == 0 || i == last { "4" } else { "3" }))
            .collect();
        let knots: Vec<_> = table
            .knots
            .iter()
            .map(|knot| text(number::rational(knot, digits).text))
            .collect();
        let curve = self.graph.add(
            "B_SPLINE_CURVE_WITH_KNOTS",
            vec![
                text("''"),
                text("3"),
                refs(controls),
                text(".UNSPECIFIED."),
                boolean(false),
                boolean(false),
                Field::List(multiplicities),
                Field::List(knots),
                text(".UNSPECIFIED."),
            ],
        );
        let representation = self.graph.add(
            "DEFINITIONAL_REPRESENTATION",
            vec![text("''"), refs([curve]), context.into()],
        );
        let pcurve = self.graph.add(
            "PCURVE",
            vec![text("''"), surface.into(), representation.into()],
        );
        Ok((pcurve, rounding))
    }

    // One boundary row as a rational B-spline curve over the surface's own
    // controls, with its exact line on the envelope and its fitted neighbour.
    pub(super) fn envelope_edge(
        &mut self,
        id: EdgeId,
        plan: &EnvelopeEdge,
        topology: &mut Topology,
    ) -> Result<Rc<Edge>, Refusal> {
        let written = self.written_envelope(plan.surface)?;
        let source = self.model.edge(id).clone();
        let (Some(start), Some(end)) = (source.start, source.end) else {
            return Err(unsupported());
        };
        let start = self.vertex(start, topology)?;
        let end = self.vertex(end, topology)?;
        let columns = written.columns;
        let indices: Vec<usize> = match plan.row {
            BoundaryRow::Rail => (0..columns).map(|i| i * 3).collect(),
            BoundaryRow::Mouth => (0..columns).map(|i| i * 3 + 2).collect(),
            BoundaryRow::Negative => (0..3).collect(),
            BoundaryRow::Positive => (0..3).map(|j| (columns - 1) * 3 + j).collect(),
        };
        let along_u = matches!(plan.row, BoundaryRow::Rail | BoundaryRow::Mouth);
        let (degree, multiplicities, knots) = if along_u {
            (
                self.rolling_envelope(plan.surface)?.carrier().degrees()[0],
                written.knots.0.clone(),
                written.knots.1.clone(),
            )
        } else {
            (
                2,
                Field::List(vec![text("3"), text("3")]),
                Field::List(vec![text("0."), text("1.")]),
            )
        };
        let weights: Vec<Field> = indices
            .iter()
            .map(|i| text(written.weights[*i].clone()))
            .collect();
        let curve = self.graph.add_complex(vec![
            ("BOUNDED_CURVE", vec![]),
            (
                "B_SPLINE_CURVE",
                vec![
                    text(degree.to_string()),
                    refs(indices.iter().map(|i| written.points[*i])),
                    text(".UNSPECIFIED."),
                    boolean(false),
                    boolean(false),
                ],
            ),
            (
                "B_SPLINE_CURVE_WITH_KNOTS",
                vec![multiplicities, knots, text(".UNSPECIFIED.")],
            ),
            ("CURVE", vec![]),
            ("GEOMETRIC_REPRESENTATION_ITEM", vec![]),
            ("RATIONAL_B_SPLINE_CURVE", vec![Field::List(weights)]),
            ("REPRESENTATION_ITEM", vec![text("''")]),
        ]);
        let context = self.parameter_context();
        let line = match plan.row {
            BoundaryRow::Rail => self.pcurve_line(written.id, context, ["0.", "0."], "(1.,0.)"),
            BoundaryRow::Mouth => self.pcurve_line(written.id, context, ["0.", "1."], "(1.,0.)"),
            BoundaryRow::Negative => {
                self.pcurve_line(written.id, context, ["-1.", "0."], "(0.,1.)")
            }
            BoundaryRow::Positive => self.pcurve_line(written.id, context, ["1.", "0."], "(0.,1.)"),
        };
        let tolerance = Rational::from_ratio(1, TABLE_RADIANS).expect("nonzero");
        let table = self
            .rolling_envelope(plan.surface)?
            .parameter_table(plan.row, &tolerance, self.dec)?;
        let (fitted, rounding) = self.neighbour_pcurve(plan, &table, context)?;
        let radius = number::upper(&plan.neighbour.radius, self.dec)?;
        let bounded = self.graph.add(
            "SURFACE_CURVE",
            vec![
                text("''"),
                curve.into(),
                refs([line, fitted]),
                text(".CURVE_3D."),
            ],
        );
        let gap = |point: &Point, other: &Point| {
            (0..3).fold(Rational::zero(), |sum, i| {
                sum + (&point.xyz[i].value - &other.xyz[i].value).abs()
            })
        };
        let first = &written.xyz[indices[0]];
        let last = &written.xyz[*indices.last().unwrap()];
        let same = gap(first, &start.point) + gap(last, &end.point)
            <= gap(first, &end.point) + gap(last, &start.point);
        let bounds = {
            let mut bounds = None;
            for i in &indices {
                include_bounds(&mut bounds, &point_bounds(&written.xyz[*i]));
            }
            bounds.ok_or_else(unsupported)?
        };
        self.deviation(&written.error + &(&(&table.bound + &rounding) * &radius));
        let edge = Rc::new(Edge {
            id: self.edge_entity(start.id, end.id, bounded, same),
            bounds,
        });
        topology.edges.insert(id, edge.clone());
        Ok(edge)
    }
}
