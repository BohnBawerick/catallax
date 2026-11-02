//! A closed rolling-pipe band. Its rails stay the exact ellipses its neighbours
//! share. The one added edge is the seam at `u = 0`, a section arc written from
//! the surface's own first column. The wall rail's vertex sits there because
//! the wall's chart is framed to cut at `t = 0`.
//!
//! The carrier's weights grow by four over each quarter while its controls
//! close, and a reader that sees closed controls treats the surface as periodic
//! and takes the last column's weights from the first. Multiplying column `j`
//! by `2^(-j/2)` gives every knot column the same weights. That is a Mobius
//! change of parameter on each quarter, so the point set does not change. On
//! the written chart `s = sqrt(2) T / (1 + (sqrt(2) - 1) T)` for
//! `T = tan(h / 2)` and `h` the angle into the quarter.
//!
//! A rail's chart curve runs on the ellipse's own parameter, the path angle `t`
//! up to direction, and no rational curve matches `s(h)`, so it is written as
//! quintic Hermite spans. `s` has its nearest pole `3 pi / 4` from the quarter,
//! and Cauchy's bound on a disc of radius 2 gives `|s^(6)| <= 81` there, so a
//! span of width `w` misses by at most `81 w^6 / 46080`.

use super::*;
use crate::graph::Field;
use catallax_numbers::{Ball, Outcome, StrictSign};
use catallax_shapes::{RollingPipe, SplineSupport};
use catallax_topology::{pipe_outward, pipe_rows};

const BITS: u32 = 192;
// Spans per quarter turn: 32 keeps `81 w^6 / 46080` below `3e-11` in `u`.
const SPANS: i64 = 32;

impl Exporter<'_> {
    pub(super) fn plan_pipe_rails(
        &mut self,
        bodies: impl Iterator<Item = BodyId>,
    ) -> Result<(), Refusal> {
        for body in bodies {
            for face in self.model.faces_of_body(body) {
                if self.rolling_pipe(face).is_some() {
                    let rows = pipe_rows(self.model, face, self.dec)?;
                    for (row, (edge, _)) in rows.ok_or_else(unsupported)?.iter().enumerate() {
                        self.pipe_rails.insert(*edge, (face, row));
                    }
                }
            }
        }
        Ok(())
    }

    // Value, slope and curvature of `u(t)` at span ends, in increasing `t`.
    fn rail_nodes(&mut self) -> Result<Vec<[Rational; 3]>, Refusal> {
        let pi = Ball::pi_checked(BITS, self.dec)?;
        let exact = |value: i64| Ball::from_rational(&Rational::from_i64(value), BITS);
        let root = exact(2).sqrt().ok_or_else(unsupported)?;
        let gap = &root - &exact(1);
        let mut nodes = Vec::with_capacity(4 * SPANS as usize + 4);
        for quarter in 0..4 {
            for span in 0..=SPANS {
                let fraction = Rational::from_ratio(span, 4 * SPANS).expect("nonzero");
                let half = &pi * &Ball::from_rational(&fraction, BITS);
                let tangent = half.sin().checked_div(&half.cos());
                let tangent = tangent.ok_or_else(unsupported)?;
                let turn = &(&tangent.square() + &exact(1))
                    * &Ball::from_rational(&Rational::from_ratio(1, 2).expect("nonzero"), BITS);
                let bend = &tangent * &turn;
                let under = &exact(1) + &(&gap * &tangent);
                let inverse = under.recip().ok_or_else(unsupported)?;
                let inverse_sq = inverse.square();
                let value = &(&root * &tangent) * &inverse;
                let first = &root * &inverse_sq;
                let second = -&(&(&(&exact(2) * &root) * &gap) * &(&inverse_sq * &inverse));
                let slope = &first * &turn;
                let curve = &(&second * &turn.square()) + &(&first * &bend);
                let value = match span {
                    0 => Rational::zero(),
                    _ if span == SPANS => Rational::one(),
                    _ => value.mid(),
                };
                nodes.push([
                    Rational::from_i64(quarter) + value,
                    slope.mid(),
                    curve.mid(),
                ]);
            }
        }
        Ok(nodes)
    }

    // The rail on the band's chart, when `id` is a band's rail.
    pub(super) fn rail_pcurve(&mut self, id: EdgeId) -> Result<Option<Id>, Refusal> {
        let Some(&(face, row)) = self.pipe_rails.get(&id) else {
            return Ok(None);
        };
        let pipe = self.rolling_pipe(face).ok_or_else(unsupported)?;
        let written = self.written_envelope(self.model.face(face).surface)?;
        let catallax_shapes::Curve::Ellipse(ellipse) =
            self.model.geom().curve(self.model.edge(id).curve)
        else {
            return Err(unsupported());
        };
        let major = ellipse.major_sq.sqrt();
        let turn = ellipse.normal.dot(&pipe.chart_normal());
        let forward = self.positive(&turn)?;
        let digits = self.digits();
        let nodes = self.rail_nodes()?;
        let pi = Ball::pi_checked(BITS, self.dec)?;
        let width = (&pi
            * &Ball::from_rational(&Rational::from_ratio(1, 2 * SPANS).expect("nonzero"), BITS))
            .hi();
        let fifth = Rational::from_ratio(1, 5).expect("nonzero");
        let twentieth = Rational::from_ratio(1, 20).expect("nonzero");
        let rail_v = if row == 0 { "0." } else { "1." };
        let mut rounding = Rational::zero();
        let mut controls = Vec::new();
        let spans = 4 * SPANS as usize;
        for index in 0..spans {
            // Each quarter stores both ends of its spans, which join at a corner.
            let quarter = index / SPANS as usize;
            let at = index + quarter;
            let (a, b) = (&nodes[at], &nodes[at + 1]);
            let mut bezier = [
                a[0].clone(),
                &a[0] + &(&(&width * &a[1]) * &fifth),
                &(&a[0] + &(&(&width * &a[1]) * &(&fifth + &fifth)))
                    + &(&(&(&width * &width) * &a[2]) * &twentieth),
                &(&b[0] - &(&(&width * &b[1]) * &(&fifth + &fifth)))
                    + &(&(&(&width * &width) * &b[2]) * &twentieth),
                &b[0] - &(&(&width * &b[1]) * &fifth),
                b[0].clone(),
            ];
            if !forward {
                bezier.reverse();
            }
            controls.push(bezier);
        }
        if !forward {
            controls.reverse();
        }
        let mut points = Vec::with_capacity(5 * spans + 1);
        for (index, bezier) in controls.iter().enumerate() {
            for u in &bezier[usize::from(index != 0)..] {
                let u = number::rational(u, digits);
                rounding = rounding.max(u.error.clone());
                let text_u = u.text;
                points.push(self.graph.add(
                    "CARTESIAN_POINT",
                    vec![text("''"), text(format!("({text_u},{rail_v})"))],
                ));
            }
        }
        let mut multiplicities = vec!["6".to_string()];
        multiplicities.extend((1..spans).map(|_| "5".to_string()));
        multiplicities.push("6".to_string());
        let mut knots = Vec::with_capacity(spans + 1);
        for span in 0..=spans {
            let share = Rational::from_ratio(span as i64, 2 * SPANS).expect("nonzero");
            let knot = &pi * &Ball::from_rational(&share, BITS);
            knots.push(number::rational(&knot.mid(), digits).text);
        }
        let curve = self.graph.add(
            "B_SPLINE_CURVE_WITH_KNOTS",
            vec![
                text("''"),
                text("5"),
                refs(points),
                text(".UNSPECIFIED."),
                boolean(false),
                boolean(false),
                text(format!("({})", multiplicities.join(","))),
                text(format!("({})", knots.join(","))),
                text(".UNSPECIFIED."),
            ],
        );
        let context = self.parameter_context();
        let representation = self.graph.add(
            "DEFINITIONAL_REPRESENTATION",
            vec![text("''"), refs([curve]), context.into()],
        );
        let pcurve = self.graph.add(
            "PCURVE",
            vec![text("''"), written.id.into(), representation.into()],
        );
        // A miss in `u` moves the point at most `2 (major + radius)` per unit,
        // since `dt/du <= 2`, and the fit misses by at most `3e-11`.
        let reach = number::upper(&(&major + pipe.radius()), self.dec)?;
        let miss = Rational::from_ratio(3, 100_000_000_000).expect("nonzero") + rounding;
        self.deviation(miss * reach * Rational::from_i64(2) + written.error.clone());
        Ok(Some(pcurve))
    }

    pub(super) fn rolling_pipe(&self, id: FaceId) -> Option<RollingPipe> {
        match self.model.geom().surface(self.model.face(id).surface) {
            Surface::BSpline(SplineSupport::RollingPipe(pipe)) => Some(pipe.clone()),
            _ => None,
        }
    }

    fn positive(&mut self, value: &Real) -> Result<bool, Refusal> {
        match self.dec.sign(value) {
            Outcome::Proven(StrictSign::Positive) => Ok(true),
            Outcome::Proven(StrictSign::Negative) => Ok(false),
            Outcome::Refused(refusal) => Err(refusal),
            _ => Err(unsupported()),
        }
    }

    // The rail's start vertex, which must be the chart corner `corner`.
    fn rail_vertex(
        &mut self,
        rail: EdgeId,
        corner: &Point3,
        topology: &mut Topology,
    ) -> Result<Rc<Vertex>, Refusal> {
        let source = self.model.edge(rail).clone();
        let (vertex, at) = if let Some(start) = source.start {
            let at = self.model.geom().position(self.model.vertex(start).point);
            (self.vertex(start, topology)?, at)
        } else {
            let at = self
                .model
                .geom()
                .curve(source.curve)
                .parameter_origin(self.model.geom());
            (topology.ring_vertices[&rail].clone(), at)
        };
        match self.dec.sign(&at.minus(corner).norm_squared()) {
            Outcome::ProvenCoincident(_) => Ok(vertex),
            Outcome::Refused(refusal) => Err(refusal),
            _ => Err(Refusal::new(
                "step.unsupported-geometry",
                "a pipe rail's vertex is not at the band's seam",
            )),
        }
    }

    // Whether the rail runs with the chart's `u`.
    fn rail_forward(&mut self, rail: EdgeId, pipe: &RollingPipe) -> Result<bool, Refusal> {
        let source = self.model.edge(rail);
        let reversed = source.sense == Sense::Reversed;
        let catallax_shapes::Curve::Ellipse(ellipse) = self.model.geom().curve(source.curve) else {
            return Err(unsupported());
        };
        let turn = ellipse.normal.dot(&pipe.chart_normal());
        Ok(self.positive(&turn)? != reversed)
    }

    /// The band's one loop, walked counter-clockwise on its chart, and whether
    /// the chart's normal is the outward one.
    pub(super) fn pipe_loop(
        &mut self,
        id: FaceId,
        pipe: &RollingPipe,
        topology: &mut Topology,
    ) -> Result<(Id, Bounds, bool), Refusal> {
        if self.chart_state.is_none() {
            return Err(Refusal::new(
                "step.unsupported-geometry",
                "a pipe's wall has no chart that cuts at the band's seam",
            ));
        }
        let surface = self.model.face(id).surface;
        let rows = pipe_rows(self.model, id, self.dec)?.ok_or_else(unsupported)?;
        let [(low, _), (high, _)] = rows;
        let rings = [self.edge(low, topology)?, self.edge(high, topology)?];
        let zero = Rational::zero();
        let corners = [Rational::zero(), Rational::one()].map(|v| pipe.point(&zero, &v));
        let start = self.rail_vertex(low, &corners[0], topology)?;
        let end = self.rail_vertex(high, &corners[1], topology)?;
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
            ["4.", "0."].map(|u| self.pcurve_line(written.id, context, [u, "0."], "(0.,1.)"));
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
        // Each seam end is a rounded control beside a rounded vertex.
        let gap = |point: &Point, other: &Point| {
            (0..3).fold(Rational::zero(), |sum, i| {
                sum + (&point.xyz[i].value - &other.xyz[i].value).abs()
            })
        };
        let ends = gap(&start.point, &written.xyz[0]).max(gap(&end.point, &written.xyz[2]));
        self.deviation(ends + self.surface_error(surface, &bounds));
        // A reader pairs the seam's two chart curves from each use in turn and
        // keeps the last pairing, so the use up the far side comes last.
        let uses = [
            (seam, false),
            (rings[0].id, self.rail_forward(low, pipe)?),
            (seam, true),
            (rings[1].id, !self.rail_forward(high, pipe)?),
        ];
        let uses: Vec<_> = uses
            .into_iter()
            .map(|(edge, same)| self.oriented(edge, same))
            .collect();
        let controls = pipe.carrier().controls();
        let along = controls[3].minus(&controls[0]);
        let across = controls[1].minus(&controls[0]);
        let law = pipe.law();
        let chart_along_ray =
            self.positive(&(&along.cross(&across).dot(&law.normal) * &law.reach))?;
        let outward = pipe_outward(self.model, id, low, self.dec)?;
        Ok((
            self.graph.add("EDGE_LOOP", vec![text("''"), refs(uses)]),
            bounds,
            chart_along_ray == outward,
        ))
    }
}
