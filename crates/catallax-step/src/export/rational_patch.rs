//! Full rational charts, shared physical rows and bounded decimal knot changes.

use super::*;
use crate::graph::Field;
use catallax_numbers::Integer;
use catallax_shapes::{BSplineSurface, BoundaryRow, Curve, RationalPatch, SplineBoundary};
use catallax_topology::RationalPatchFace;
use std::collections::BTreeSet;

pub(super) struct Parameters {
    pub degrees: [usize; 2],
    pub rows: usize,
    pub v_knots: (Field, Field),
}

struct Knots {
    fields: (Field, Field),
    error: Rational,
    source_gap: Rational,
    written_gap: Rational,
}

fn knots(source: &[Rational], digits: u32, dec: &Decider) -> Result<Knots, Refusal> {
    let origin = &source[0];
    let width = source.last().expect("clamped finite axis") - origin;
    let mut multiplicities = Vec::new();
    let mut values = Vec::new();
    let mut error = Rational::zero();
    let mut previous = None::<(Rational, Rational)>;
    let mut gaps = None::<[Rational; 2]>;
    let mut start = 0;
    while start < source.len() {
        dec.check_deadline()?;
        let mut end = start + 1;
        while end < source.len() && source[end] == source[start] {
            end += 1;
        }
        let parameter = &(&source[start] - origin) / &width;
        let number = number::rational(&parameter, digits);
        if let Some((old, written)) = &previous {
            let gap = [&parameter - old, &number.value - written];
            if gap[1] <= Rational::zero() {
                return Err(number::budget(
                    "STEP decimal precision collapses a rational knot span",
                ));
            }
            if let Some(gaps) = &mut gaps {
                for i in 0..2 {
                    if gap[i] < gaps[i] {
                        gaps[i] = gap[i].clone();
                    }
                }
            } else {
                gaps = Some(gap);
            }
        }
        error = error.max(number.error);
        multiplicities.push(text((end - start).to_string()));
        values.push(text(number.text));
        previous = Some((parameter, number.value));
        start = end;
    }
    let [source_gap, written_gap] = gaps.ok_or_else(unsupported)?;
    Ok(Knots {
        fields: (Field::List(multiplicities), Field::List(values)),
        error,
        source_gap,
        written_gap,
    })
}

fn sensitivity(knots: &Knots, degree: usize) -> Rational {
    if knots.error.is_zero() {
        return Rational::zero();
    }
    let coefficient =
        &Rational::from_i64(16) / &knots.written_gap + &Rational::from_i64(4) / &knots.source_gap;
    &(&coefficient * &knots.error) * &Rational::from_integer(Integer::from_u64(degree as u64))
}

pub(super) fn parameters(
    carrier: &BSplineSurface,
    weights: &[Rational],
    diameter: &Rational,
    digits: u32,
    dec: &Decider,
) -> Result<((Field, Field), Parameters, Rational), Refusal> {
    let source = carrier.knots();
    let u = knots(source[0], digits, dec)?;
    let v = knots(source[1], digits, dec)?;
    let degrees = carrier.degrees();
    let error = if u.error.is_zero() && v.error.is_zero() {
        Rational::zero()
    } else {
        let minimum = weights.iter().min().ok_or_else(unsupported)?;
        let maximum = weights.iter().max().ok_or_else(unsupported)?;
        let ratio = maximum / minimum;
        // Convex de Boor blends cost at most 8*p*epsilon/min_written_span
        // in coefficient L1 norm. Positive-weight normalisation costs 2*ratio.
        // The parameter shift adds 4*p*ratio*diameter*epsilon/min_source_span.
        &(diameter * &ratio) * &(sensitivity(&u, degrees[0]) + sensitivity(&v, degrees[1]))
    };
    Ok((
        u.fields,
        Parameters {
            degrees,
            rows: carrier.control_count()[1],
            v_knots: v.fields,
        },
        error,
    ))
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) struct Row {
    edge: EdgeId,
    surface: SurfaceId,
    boundary: BoundaryRow,
    forward: bool,
}

// SURFACE_CURVE admits at most two associated pcurves.
pub(super) type Rows = [Option<Row>; 2];

fn bad(detail: &str) -> Refusal {
    Refusal::new("step.rational-patch-unproven", detail)
}

fn side(row: BoundaryRow) -> (usize, SplineBoundary) {
    match row {
        BoundaryRow::Rail => (0, SplineBoundary::VMin),
        BoundaryRow::Mouth => (0, SplineBoundary::VMax),
        BoundaryRow::Negative => (1, SplineBoundary::UMin),
        BoundaryRow::Positive => (1, SplineBoundary::UMax),
    }
}

fn support(model: &Model, surface: SurfaceId) -> Result<&RationalPatch, Refusal> {
    match model.geom().surface(surface) {
        Surface::BSpline(catallax_shapes::SplineSupport::RationalPatch(patch)) => Ok(patch),
        _ => Err(unsupported()),
    }
}

fn face_ids(model: &Model, body: BodyId) -> impl Iterator<Item = FaceId> + '_ {
    model
        .body(body)
        .lumps
        .iter()
        .flat_map(move |lump| model.lump(*lump).shells.iter())
        .flat_map(move |shell| model.shell(*shell).faces.iter().copied())
}

impl Exporter<'_> {
    pub(super) fn plan_rational_patches(
        &mut self,
        bodies: impl Iterator<Item = BodyId> + Clone,
    ) -> Result<(), Refusal> {
        let faces: BTreeSet<_> = bodies
            .clone()
            .flat_map(|body| face_ids(self.model, body))
            .filter(|face| {
                matches!(
                    self.model.geom().surface(self.model.face(*face).surface),
                    Surface::BSpline(catallax_shapes::SplineSupport::RationalPatch(_))
                )
            })
            .collect();
        if faces.is_empty() {
            return Ok(());
        }
        for face in faces {
            let checked =
                RationalPatchFace::of_face(self.model, face, self.dec)?.ok_or_else(|| {
                    bad("rational STEP faces require their complete checked physical domain")
                })?;
            let surface = self.model.face(face).surface;
            let rows = checked.rows().map(|(edge, boundary, _, forward)| Row {
                edge,
                surface,
                boundary,
                forward,
            });
            for row in rows {
                let uses = self.rational_rows.entry(row.edge).or_insert([None; 2]);
                if uses.iter().flatten().any(|previous| *previous == row) {
                    continue;
                }
                let slot = uses.iter_mut().find(|slot| slot.is_none()).ok_or_else(|| {
                    bad("a shared STEP curve needs more than two rational pcurve charts")
                })?;
                *slot = Some(row);
            }
            self.rational_faces.insert(face, rows);
        }
        for body in bodies {
            for face in face_ids(self.model, body) {
                let surface = self.model.geom().surface(self.model.face(face).surface);
                if matches!(
                    surface,
                    Surface::Plane(_)
                        | Surface::BSpline(catallax_shapes::SplineSupport::RationalPatch(_))
                ) {
                    continue;
                }
                for lp in &self.model.face(face).loops {
                    if let Loop::Coedges(coedges) = self.model.loop_(*lp)
                        && coedges
                            .iter()
                            .any(|c| self.rational_rows.contains_key(&self.model.coedge(*c).edge))
                    {
                        return Err(bad(
                            "a rational row on a non-planar law-specific neighbour needs a certified pcurve",
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    fn rational_pcurve(&mut self, surface: Id, row: BoundaryRow, forward: bool) -> Id {
        let mut coordinates = match row {
            BoundaryRow::Rail => ["(0.,0.)", "(1.,0.)"],
            BoundaryRow::Mouth => ["(0.,1.)", "(1.,1.)"],
            BoundaryRow::Negative => ["(0.,0.)", "(0.,1.)"],
            BoundaryRow::Positive => ["(1.,0.)", "(1.,1.)"],
        };
        if !forward {
            coordinates.reverse();
        }
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

    pub(super) fn rational_edge(
        &mut self,
        id: EdgeId,
        topology: &mut Topology,
    ) -> Result<Rc<Edge>, Refusal> {
        let plans = self.rational_rows[&id];
        let source = plans[0].expect("registered rational edge has an authority row");
        let model = self.model;
        let authority = support(model, source.surface)?;
        let written = self.written_envelope(source.surface)?;
        let parameters = written.rational.as_ref().ok_or_else(unsupported)?;
        let (axis, boundary) = side(source.boundary);
        let indices: Vec<_> = match source.boundary {
            BoundaryRow::Rail => (0..written.columns).map(|i| i * parameters.rows).collect(),
            BoundaryRow::Mouth => (0..written.columns)
                .map(|i| (i + 1) * parameters.rows - 1)
                .collect(),
            BoundaryRow::Negative => (0..parameters.rows).collect(),
            BoundaryRow::Positive => ((written.columns - 1) * parameters.rows
                ..written.columns * parameters.rows)
                .collect(),
        };
        let knots = if axis == 0 {
            &written.knots
        } else {
            &parameters.v_knots
        };
        let curve = self.spline_row(
            &written,
            &indices,
            parameters.degrees[axis],
            knots.0.clone(),
            knots.1.clone(),
            false,
        );
        let mut pcurves = [None; 2];
        for (index, row) in plans.into_iter().flatten().enumerate() {
            let direction = if row == source {
                true
            } else {
                let target = support(model, row.surface)?;
                authority
                    .boundary_orientation(
                        boundary,
                        &Curve::BSpline(target.boundary(side(row.boundary).1)),
                        self.dec,
                    )?
                    .ok_or_else(|| {
                        bad("shared rational rows need an exact common normalised parameter chart")
                    })?
            };
            if direction != (row.forward == source.forward) {
                return Err(bad(
                    "a shared rational row reverses its physical endpoint correspondence",
                ));
            }
            let target = self.written_envelope(row.surface)?;
            if row != source {
                self.deviation(&written.error + &target.error);
            }
            pcurves[index] = Some(self.rational_pcurve(target.id, row.boundary, direction));
        }
        let bounded = self.graph.add(
            "SURFACE_CURVE",
            vec![
                text("''"),
                curve.into(),
                refs(pcurves.into_iter().flatten()),
                text(".CURVE_3D."),
            ],
        );
        let source_edge = model.edge(id);
        let (Some(a), Some(b)) = (source_edge.start, source_edge.end) else {
            return Err(unsupported());
        };
        let start = self.vertex(a, topology)?;
        let end = self.vertex(b, topology)?;
        let endpoint_indices = if source.forward {
            [indices[0], *indices.last().expect("finite row")]
        } else {
            [*indices.last().expect("finite row"), indices[0]]
        };
        for (vertex, index) in [&start, &end].into_iter().zip(endpoint_indices) {
            let gap = (0..3).fold(Rational::zero(), |sum, axis| {
                sum + (&vertex.point.xyz[axis].value - &written.xyz[index].xyz[axis].value).abs()
            });
            if gap > self.options.curve_tolerance {
                return Err(number::budget(
                    "rounded rational row endpoint exceeds the file's sewing tolerance",
                ));
            }
        }
        let mut bounds = None;
        for index in indices {
            include_bounds(&mut bounds, &point_bounds(&written.xyz[index]));
        }
        let edge = Rc::new(Edge {
            id: self.edge_entity(start.id, end.id, bounded, source.forward),
            bounds: bounds.ok_or_else(unsupported)?,
        });
        topology.edges.insert(id, edge.clone());
        Ok(edge)
    }

    pub(super) fn rational_loop(
        &mut self,
        rows: [Row; 4],
        topology: &mut Topology,
    ) -> Result<(Id, Bounds), Refusal> {
        let mut uses = [None; 4];
        let mut bounds = None;
        for (walk, index) in [0, 3, 1, 2].into_iter().enumerate() {
            let row = rows[index];
            let edge = self.edge(row.edge, topology)?;
            include_bounds(&mut bounds, &edge.bounds);
            uses[walk] = Some(self.oriented(edge.id, row.forward == (walk < 2)));
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
