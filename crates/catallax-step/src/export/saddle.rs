//! A closed saddle band between two cylinders. Its rails are the surface's own
//! first and last rows, so a shared rail cannot drift from the band it bounds.
//! The one added edge is the seam at `u = -1`, the first column of controls.
//!
//! Each wall keeps its exact cylinder, and the rail is exactly on it. The rail's
//! chart curve on the wall has to run on the band's parameter, where the angle
//! is `2 atan(f)` and no polynomial matches it. So it comes from the band's
//! table of cubic Hermite spans, each raised to the degree the chart code reads.

use super::*;
use crate::fit;
use crate::graph::Field;
use crate::junction::{Neighbour, saddle_neighbour};
use crate::uv;
use catallax_numbers::{Ball, Outcome, StrictSign};
use catallax_shapes::{ParameterTable, RollingSaddle, SplineSupport};
use catallax_topology::{saddle_added, saddle_rows};

// Inverse of a rail table's angular tolerance, as for a rolling envelope.
const TABLE_RADIANS: i64 = 100_000_000_000;

pub(in crate::export) struct RailChart {
    table: ParameterTable,
    neighbour: Neighbour,
    // One whole turn of the chart.
    turn: Rational,
}

impl RailChart {
    // Chart value of a table entry, and its slope.
    fn node(&self, index: usize, turns: &Rational) -> [[Rational; 2]; 2] {
        let scaled = |sign: i64, data: &[Rational; 2]| {
            let sign = Rational::from_i64(sign);
            [&sign * &data[0], &sign * &data[1]]
        };
        let [angle, slope] = scaled(self.neighbour.u_sign, &self.table.angle[index]);
        let [height, rise] = scaled(self.neighbour.v_sign, &self.table.height[index]);
        [
            [&(&self.neighbour.u.mid() + turns) + &angle, slope],
            [&self.neighbour.v.mid() + &height, rise],
        ]
    }

    /// The chart curve, with its ends on the written vertices `ends`.
    pub(super) fn fit(&self, ends: [&[Rational; 2]; 2], digits: u32) -> uv::Fit {
        // The start vertex's chart angle fixes the turn the table lies in.
        let first = self.node(0, &Rational::zero());
        let tau = Rational::from_ratio(314_159_265_358_979, 50_000_000_000_000).expect("nonzero");
        let turns = if &ends[0][0] - &first[0][0] > tau.clone().mul_pow2(-1) {
            self.turn.clone()
        } else if &first[0][0] - &ends[0][0] > tau.mul_pow2(-1) {
            -self.turn.clone()
        } else {
            Rational::zero()
        };
        let spans = self.table.knots.len() - 1;
        let mut points = Vec::with_capacity(fit::DEGREE * spans + 1);
        let mut rounding = self.neighbour.u.rad() + self.neighbour.v.rad();
        let mut push = |point: [Rational; 2], points: &mut Vec<[number::Number; 2]>| {
            let point = point.map(|value| number::rational(&value, digits));
            rounding = rounding
                .clone()
                .max(point[0].error.clone() + point[1].error.clone());
            points.push(point);
        };
        let third = Rational::from_ratio(1, 3).expect("nonzero");
        for span in 0..spans {
            let width = (&self.table.knots[span + 1] - &self.table.knots[span]) * third.clone();
            let [[u0, du0], [v0, dv0]] = self.node(span, &turns);
            let [[u1, du1], [v1, dv1]] = self.node(span + 1, &turns);
            let cubic = [
                [u0.clone(), v0.clone()],
                [&u0 + &(&width * &du0), &v0 + &(&width * &dv0)],
                [&u1 - &(&width * &du1), &v1 - &(&width * &dv1)],
                [u1, v1],
            ];
            for index in usize::from(span != 0)..=fit::DEGREE {
                push(raised(&cubic, index), &mut points);
            }
        }
        // The stored rows meet at the written vertex; the miss joins the error.
        let mut miss = Rational::zero();
        for (at, end) in [(0, ends[0]), (points.len() - 1, ends[1])] {
            for axis in 0..2 {
                miss = miss.max((&points[at][axis].value - &end[axis]).abs());
                points[at][axis] = number::rational(&end[axis], digits);
            }
        }
        let last = self.table.knots.len() - 1;
        uv::Fit {
            points,
            knots: self.table.knots.clone(),
            multiplicities: (0..=last)
                .map(|i| {
                    if i == 0 || i == last {
                        fit::DEGREE + 1
                    } else {
                        fit::DEGREE
                    }
                })
                .collect(),
            error: &self.table.bound + &(rounding + miss),
        }
    }
}

// Control `index` of a cubic Bezier raised to `fit::DEGREE`.
fn raised(cubic: &[[Rational; 2]; 4], index: usize) -> [Rational; 2] {
    const WEIGHTS: [i64; 4] = [1, 3, 3, 1];
    let raise = fit::DEGREE - 3;
    let choose = |n: usize, k: usize| (0..k).fold(1i64, |c, i| c * (n - i) as i64 / (i as i64 + 1));
    let whole = choose(fit::DEGREE, index);
    std::array::from_fn(|axis| {
        (0..4).fold(Rational::zero(), |sum, j| {
            if j > index || index - j > raise {
                return sum;
            }
            let share = WEIGHTS[j] * choose(raise, index - j);
            let share = Rational::from_ratio(share, whole).expect("nonzero");
            sum + &cubic[j][axis] * &share
        })
    })
}

impl Exporter<'_> {
    pub(super) fn plan_saddle_rails(
        &mut self,
        bodies: impl Iterator<Item = BodyId>,
    ) -> Result<(), Refusal> {
        for body in bodies {
            for face in self.model.faces_of_body(body) {
                if self.rolling_saddle(face).is_some() {
                    let rows = saddle_rows(self.model, face).ok_or_else(unsupported)?;
                    for (row, (edge, _)) in rows.iter().enumerate() {
                        self.saddle_rails.insert(*edge, (face, row));
                    }
                }
            }
        }
        Ok(())
    }

    pub(super) fn rolling_saddle(&self, id: FaceId) -> Option<RollingSaddle> {
        match self.model.geom().surface(self.model.face(id).surface) {
            Surface::BSpline(SplineSupport::RollingSaddle(saddle)) => Some(saddle.clone()),
            _ => None,
        }
    }

    fn sign_of(&mut self, value: &Real) -> Result<bool, Refusal> {
        match self.dec.sign(value) {
            Outcome::Proven(StrictSign::Positive) => Ok(true),
            Outcome::Proven(StrictSign::Negative) => Ok(false),
            Outcome::Refused(refusal) => Err(refusal),
            _ => Err(unsupported()),
        }
    }

    // A rail as a row of the band's written controls, with its line on the band
    // and its table on the wall.
    pub(super) fn saddle_edge(
        &mut self,
        id: EdgeId,
        topology: &mut Topology,
    ) -> Result<Rc<Edge>, Refusal> {
        let (face, row) = self.saddle_rails[&id];
        let saddle = self.rolling_saddle(face).ok_or_else(unsupported)?;
        let written = self.written_envelope(self.model.face(face).surface)?;
        let source = self.model.edge(id).clone();
        let (Some(from), Some(to)) = (source.start, source.end) else {
            return Err(unsupported());
        };
        let start = self.vertex(from, topology)?;
        let end = self.vertex(to, topology)?;
        let wall = {
            let state = self.chart_state.as_ref().ok_or_else(unsupported)?;
            let uses = state.plan.uses.get(&id).ok_or_else(unsupported)?;
            let [(wall, _)] = uses.as_slice() else {
                return Err(unsupported());
            };
            *wall
        };
        let indices: Vec<usize> = (0..written.columns).map(|i| i * 3 + 2 * row).collect();
        let degree = saddle.carrier().degrees()[0];
        let curve = self.spline_row(
            &written,
            &indices,
            degree,
            written.knots.0.clone(),
            written.knots.1.clone(),
        );
        let context = self.parameter_context();
        let height = if row == 0 { "0." } else { "1." };
        let line = self.pcurve_line(written.id, context, ["0.", height], "(1.,0.)");
        let tolerance = Rational::from_ratio(1, TABLE_RADIANS).expect("nonzero");
        let table = saddle.parameter_table(row, &tolerance, self.dec)?;
        let surface = self.model.face(wall).surface;
        let neighbour = saddle_neighbour(self.model, surface, &saddle, row, self.dec)?;
        let turn = Ball::pi_checked(128, self.dec)?.mid() * Rational::from_i64(2);
        let bounds = {
            let mut bounds = None;
            for i in &indices {
                include_bounds(&mut bounds, &point_bounds(&written.xyz[*i]));
            }
            bounds.ok_or_else(unsupported)?
        };
        let same = source.sense == Sense::Same;
        let wrapper = self.chart_edge(
            id,
            charts::EdgeInput {
                curve,
                geometry: &charts::Geometry::Rail(Box::new(RailChart {
                    table,
                    neighbour,
                    turn,
                })),
                bounds: &bounds,
                error: &written.error,
                same,
                extra: Some(line),
            },
            topology,
        )?;
        let edge = Rc::new(Edge {
            id: self.edge_entity(start.id, end.id, wrapper, same),
            bounds,
        });
        topology.edges.insert(id, edge.clone());
        Ok(edge)
    }

    // The rail's vertex, which must be the chart corner `corner`.
    fn saddle_vertex(
        &mut self,
        rail: EdgeId,
        corner: &Point3,
        topology: &mut Topology,
    ) -> Result<Rc<Vertex>, Refusal> {
        let start = self.model.edge(rail).start.ok_or_else(unsupported)?;
        let at = self.model.geom().position(self.model.vertex(start).point);
        match self.dec.sign(&at.minus(corner).norm_squared()) {
            Outcome::ProvenCoincident(_) => self.vertex(start, topology),
            Outcome::Refused(refusal) => Err(refusal),
            _ => Err(Refusal::new(
                "step.unsupported-geometry",
                "a saddle rail's vertex is not at the band's seam",
            )),
        }
    }

    /// The band's one loop, walked counter-clockwise on its chart, and whether
    /// the chart's normal is the outward one.
    pub(super) fn saddle_loop(
        &mut self,
        id: FaceId,
        saddle: &RollingSaddle,
        topology: &mut Topology,
    ) -> Result<(Id, Bounds, bool), Refusal> {
        if self.chart_state.is_none() {
            return Err(Refusal::new(
                "step.unsupported-geometry",
                "a saddle's walls have no chart that cuts at the band's seam",
            ));
        }
        let surface = self.model.face(id).surface;
        let rows = saddle_rows(self.model, id).ok_or_else(unsupported)?;
        let [(wide, _), (narrow, _)] = rows;
        let rings = [self.edge(wide, topology)?, self.edge(narrow, topology)?];
        let controls = saddle.carrier().controls();
        let start = self.saddle_vertex(wide, &controls[0], topology)?;
        let end = self.saddle_vertex(narrow, &controls[2], topology)?;
        let written = self.written_envelope(surface)?;
        let curve = self.spline_row(
            &written,
            &[0, 1, 2],
            2,
            Field::List(vec![text("3"), text("3")]),
            Field::List(vec![text("0."), text("1.")]),
        );
        let context = self.parameter_context();
        let pcurves =
            ["3.", "-1."].map(|u| self.pcurve_line(written.id, context, [u, "0."], "(0.,1.)"));
        let seam_curve = self.graph.add(
            "SEAM_CURVE",
            vec![text("''"), curve.into(), refs(pcurves), text(".CURVE_3D.")],
        );
        let seam = self.edge_entity(start.id, end.id, seam_curve, true);
        self.report.added_edges += 1;
        let mut bounds = None;
        for ring in &rings {
            include_bounds(&mut bounds, &ring.bounds);
        }
        let bounds = bounds.ok_or_else(unsupported)?;
        let gap = |point: &Point, other: &Point| {
            (0..3).fold(Rational::zero(), |sum, i| {
                sum + (&point.xyz[i].value - &other.xyz[i].value).abs()
            })
        };
        let ends = gap(&start.point, &written.xyz[0]).max(gap(&end.point, &written.xyz[2]));
        self.deviation(ends + self.surface_error(surface, &bounds));
        // A rail runs up `u`, and an edge of opposite sense runs down it.
        let forward = |exporter: &Self, rail| exporter.model.edge(rail).sense == Sense::Same;
        // A reader pairs the seam's two chart curves from each use in turn and
        // keeps the last pairing, so the use up the far side comes last.
        let uses = [
            (seam, false),
            (rings[0].id, forward(self, wide)),
            (seam, true),
            (rings[1].id, !forward(self, narrow)),
        ];
        let uses: Vec<_> = uses
            .into_iter()
            .map(|(edge, same)| self.oriented(edge, same))
            .collect();
        // The chart's normal is `d/du x d/dv`, at the wide rail's first corner.
        let along = controls[3].minus(&controls[0]);
        let across = controls[1].minus(&controls[0]);
        let centre = saddle.definition().section(&Rational::from_i64(-1))?;
        let toward_ball =
            self.sign_of(&along.cross(&across).dot(&centre[0].minus(&controls[0])))?;
        // Added material lies away from the ball, so its outward side faces it.
        let added = saddle_added(self.model, id)?;
        Ok((
            self.graph.add("EDGE_LOOP", vec![text("''"), refs(uses)]),
            bounds,
            toward_ball == added,
        ))
    }
}
