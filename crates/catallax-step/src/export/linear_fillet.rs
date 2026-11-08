//! Exact contact-section bands. The arc is `u`, the authored station is `v`.
//! Every physical boundary uses the surface's own rounded control row and its
//! straight unit-interval pcurve. A collapsed station has a constant 3D control
//! row and pcurve at the existing apex, but no additional physical model edge.

use super::*;
use crate::graph::Field;
use catallax_shapes::BoundaryRow;
use catallax_topology::LinearFilletBand;

#[derive(Clone, Copy)]
pub(super) struct Row {
    surface: SurfaceId,
    boundary: BoundaryRow,
    forward: bool,
}

fn slot(row: BoundaryRow) -> usize {
    match row {
        BoundaryRow::Negative => 0,
        BoundaryRow::Mouth => 1,
        BoundaryRow::Positive => 2,
        BoundaryRow::Rail => 3,
    }
}

impl Exporter<'_> {
    pub(super) fn plan_linear_bands(
        &mut self,
        bodies: impl Iterator<Item = BodyId>,
    ) -> Result<(), Refusal> {
        for body in bodies {
            for face in self.model.faces_of_body(body) {
                let Some(band) = LinearFilletBand::of_face(self.model, face, self.dec)? else {
                    continue;
                };
                let surface = self.model.face(face).surface;
                let mut rows = [None; 4];
                for (edge, boundary, _, forward) in band.rows() {
                    rows[slot(boundary)] = Some(edge);
                    let row = Row {
                        surface,
                        boundary,
                        forward,
                    };
                    if let Some(previous) = self.linear_rows.insert(edge, row)
                        && (previous.surface != row.surface
                            || previous.boundary != row.boundary
                            || previous.forward != row.forward)
                    {
                        return Err(Refusal::new(
                            "step.unsupported-geometry",
                            "two variable bands share a boundary outside the proved disjoint domain",
                        ));
                    }
                }
                self.linear_faces.insert(face, rows);
            }
        }
        Ok(())
    }

    // The unit-interval parameter is also the 3D control row's parameter.
    fn linear_pcurve(&mut self, surface: Id, row: BoundaryRow) -> Id {
        let coordinates = match row {
            BoundaryRow::Rail => ["(0.,0.)", "(0.,1.)"],
            BoundaryRow::Mouth => ["(1.,0.)", "(1.,1.)"],
            BoundaryRow::Negative => ["(0.,0.)", "(1.,0.)"],
            BoundaryRow::Positive => ["(0.,1.)", "(1.,1.)"],
        };
        let points = coordinates.map(|point| {
            self.graph
                .add("CARTESIAN_POINT", vec![text("''"), text(point)])
        });
        let curve = self.graph.add(
            "B_SPLINE_CURVE_WITH_KNOTS",
            vec![
                text("''"),
                text("1"),
                refs(points),
                text(".UNSPECIFIED."),
                boolean(false),
                boolean(false),
                text("(2,2)"),
                text("(0.,1.)"),
                text(".UNSPECIFIED."),
            ],
        );
        let context = self.parameter_context();
        let representation = self.graph.add(
            "DEFINITIONAL_REPRESENTATION",
            vec![text("''"), refs([curve]), context.into()],
        );
        self.graph.add(
            "PCURVE",
            vec![text("''"), surface.into(), representation.into()],
        )
    }

    pub(super) fn linear_edge(
        &mut self,
        id: EdgeId,
        topology: &mut Topology,
    ) -> Result<Rc<Edge>, Refusal> {
        let row = self.linear_rows[&id];
        let written = self.written_envelope(row.surface)?;
        let source = self.model.edge(id);
        let (Some(from), Some(to)) = (source.start, source.end) else {
            return Err(unsupported());
        };
        let start = self.vertex(from, topology)?;
        let end = self.vertex(to, topology)?;
        // The native carrier is degree (2,1), u-major/v-minor. Rails vary v;
        // positive end sections vary u. No independent circle parameter enters.
        let (indices, degree): (&[usize], usize) = match row.boundary {
            BoundaryRow::Rail => (&[0, 1], 1),
            BoundaryRow::Mouth => (&[4, 5], 1),
            BoundaryRow::Negative => (&[0, 2, 4], 2),
            BoundaryRow::Positive => (&[1, 3, 5], 2),
        };
        let curve = self.spline_row(
            &written,
            indices,
            degree,
            Field::List(vec![text((degree + 1).to_string()); 2]),
            Field::List(vec![text("0."), text("1.")]),
            false,
        );
        let pcurve = self.linear_pcurve(written.id, row.boundary);
        let bounded = self.graph.add(
            "SURFACE_CURVE",
            vec![text("''"), curve.into(), refs([pcurve]), text(".CURVE_3D.")],
        );
        let mut bounds = None;
        for index in indices {
            include_bounds(&mut bounds, &point_bounds(&written.xyz[*index]));
        }
        let edge = Rc::new(Edge {
            id: self.edge_entity(start.id, end.id, bounded, row.forward),
            bounds: bounds.ok_or_else(unsupported)?,
        });
        topology.edges.insert(id, edge.clone());
        Ok(edge)
    }

    // A pcurve-only pole leaves the raw translator without a usable 3D edge.
    // The surface's constant row closes the full parameter boundary exactly,
    // with the same [0,1] parameter as its pcurve and no positive-radius cap.
    fn linear_pole(
        &mut self,
        rail: EdgeId,
        end: usize,
        topology: &mut Topology,
    ) -> Result<(Id, Bounds), Refusal> {
        let plan = self.linear_rows[&rail];
        let written = self.written_envelope(plan.surface)?;
        let source = self.model.edge(rail);
        let vertex = if (end == 0) == plan.forward {
            source.start
        } else {
            source.end
        }
        .ok_or_else(unsupported)?;
        let vertex = self.vertex(vertex, topology)?;
        let indices = [end, 2 + end, 4 + end];
        let curve = self.spline_row(
            &written,
            &indices,
            2,
            Field::List(vec![text("3"), text("3")]),
            Field::List(vec![text("0."), text("1.")]),
            true,
        );
        let boundary = if end == 0 {
            BoundaryRow::Negative
        } else {
            BoundaryRow::Positive
        };
        let pcurve = self.linear_pcurve(written.id, boundary);
        let bounded = self.graph.add(
            "SURFACE_CURVE",
            vec![text("''"), curve.into(), refs([pcurve]), text(".CURVE_3D.")],
        );
        self.report.added_edges += 1;
        Ok((
            self.edge_entity(vertex.id, vertex.id, bounded, true),
            point_bounds(&vertex.point),
        ))
    }

    /// A counter-clockwise walk of the complete parameter boundary. A checked
    /// zero station is a constant file row sharing the native apex vertex.
    pub(super) fn linear_loop(
        &mut self,
        rows: [Option<EdgeId>; 4],
        topology: &mut Topology,
    ) -> Result<(Id, Bounds), Refusal> {
        let mut uses = [None; 4];
        let mut bounds = None;
        for (index, id) in rows.into_iter().enumerate() {
            let (edge, forward) = if let Some(id) = id {
                let row = self.linear_rows[&id];
                let edge = self.edge(id, topology)?;
                include_bounds(&mut bounds, &edge.bounds);
                (edge.id, row.forward)
            } else {
                if index != 0 && index != 2 {
                    return Err(unsupported());
                }
                let (edge, extent) =
                    self.linear_pole(rows[3].ok_or_else(unsupported)?, index / 2, topology)?;
                include_bounds(&mut bounds, &extent);
                (edge, true)
            };
            uses[index] = Some(self.oriented(edge, forward == (index < 2)));
        }
        Ok((
            self.graph.add(
                "EDGE_LOOP",
                vec![text("''"), refs(uses.into_iter().flatten())],
            ),
            bounds.ok_or_else(unsupported)?,
        ))
    }
}
