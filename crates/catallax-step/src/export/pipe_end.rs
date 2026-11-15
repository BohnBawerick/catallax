//! Parameter-coupled end pcurves on the finite pipe export frame.

use super::*;
use crate::{fit, uv};
use catallax_intersect::{CurveTrace, TracePiece};
use catallax_numbers::Ball;
use uv::{Jet, Parametric};

struct EndMap<'a> {
    piece: &'a TracePiece,
    source: &'a fit::TraceSpan,
    range: [Rational; 2],
    low: Rational,
    quarter: i64,
    written_quarter: i64,
}

impl Parametric for EndMap<'_> {
    fn range(&self) -> [Rational; 2] {
        self.range.clone()
    }
    fn reuse_endpoints(&self) -> bool {
        true
    }

    fn coefficients<const N: usize>(
        &self,
        t: &Ball,
        _constants: bool,
        dec: &Decider,
    ) -> Result<[[Ball; 2]; N], Refusal> {
        assert!(N > 0);
        let bits = t.precision();
        let ball = |q: &Rational| Ball::from_rational(q, bits);
        let width = &self.range[1] - &self.range[0];
        let speed = (&self.source.to - &self.source.from) / width;
        let parameter = &ball(&self.source.from) + &((t - &ball(&self.range[0])) * ball(&speed));
        let native = self
            .piece
            .section_parameters(&parameter, N - 1, dec)?
            .ok_or_else(unsupported)?;
        let jets: [Jet<N>; 2] = std::array::from_fn(|axis| {
            Jet::from_coefficients(std::array::from_fn(|i| {
                &native[i][axis] * &ball(&speed.pow(i as i32).unwrap())
            }))
        });
        let scalar = |n| Jet::constant(ball(&Rational::from_i64(n)));
        let u = jets[0]
            .scale(&ball(&Rational::from_ratio(1, 2).unwrap()))
            .add(&Jet::constant(ball(&self.low)));
        let tangent = u.sub(&scalar(self.quarter));
        let root = ball(&Rational::from_i64(2))
            .sqrt()
            .ok_or_else(unsupported)?;
        let denominator = scalar(1).add(&tangent.scale(&(&root - &ball(&Rational::one()))));
        let written = tangent
            .scale(&root)
            .div(&denominator)?
            .add(&scalar(self.written_quarter));
        Ok(std::array::from_fn(|i| {
            [
                written.coefficients[i].clone(),
                jets[1].coefficients[i].clone(),
            ]
        }))
    }
}

impl Exporter<'_> {
    pub(super) fn pipe_end_pcurve(
        &mut self,
        edge: EdgeId,
        trace: &CurveTrace,
        fit: &fit::Fit,
        spans: &[fit::TraceSpan],
    ) -> Result<Id, Refusal> {
        let face = self
            .model
            .coedges_of_edge(edge)
            .into_iter()
            .filter_map(|coedge| self.model.owner_of_coedge(coedge))
            .map(|(_, face)| face)
            .find(|face| self.rolling_pipe(*face).is_some())
            .ok_or_else(unsupported)?;
        let pipe = self.rolling_pipe(face).ok_or_else(unsupported)?;
        let surface = self.model.face(face).surface;
        let start = self.model.edge(edge).start.ok_or_else(unsupported)?;
        let point = self.model.geom().position(self.model.vertex(start).point);
        let positive = self.positive(&point.minus(&pipe.law().base).dot(&pipe.law().side))?;
        let centre = if positive { 1 } else { 3 };
        let quarter = if positive { 1 } else { 2 };
        let low = Rational::from_ratio(4 * centre - 1, 4).unwrap();
        let rotated = self.pipe_rotations.get(&surface) == Some(&true);
        let written_quarter = if rotated { (quarter + 2) % 4 } else { quarter };

        // Qualify the local chart map by spatial samples, without claiming a
        // complete-body bound for this approximate exchange path.
        for piece in trace.pieces() {
            for i in 0..=4 {
                let t = Ball::from_rational(&Rational::from_ratio(i, 4).unwrap(), 128);
                let uv = piece
                    .section_parameters(&t, 0, self.dec)?
                    .ok_or_else(unsupported)?;
                let u = &low + &uv[0][0].mid().mul_pow2(-1);
                let expected = pipe.point(&u, &uv[0][1].mid());
                let found = piece.evaluate(&t, self.dec)?;
                let mut error = Rational::zero();
                for (coordinate, found) in [&expected.x, &expected.y, &expected.z]
                    .into_iter()
                    .zip(found)
                {
                    error = error
                        + number::magnitude(&(&self.dec.enclose_at(coordinate, 128)? - &found));
                }
                if error > &self.options.curve_tolerance / &Rational::from_i64(8) {
                    return Err(Refusal::new(
                        "step.pipe-end-parameter-map",
                        "retained end chart does not match the pipe within the requested tolerance",
                    ));
                }
            }
        }
        let pieces: Vec<_> = spans
            .iter()
            .enumerate()
            .map(|(i, source)| EndMap {
                piece: &trace.pieces()[source.piece],
                source,
                range: [fit.knots[i].clone(), fit.knots[i + 1].clone()],
                low: low.clone(),
                quarter,
                written_quarter,
            })
            .collect();
        let reach = number::upper(
            &(&pipe.spine(self.dec)?.major_sq.sqrt() + pipe.radius()),
            self.dec,
        )?;
        let metric = Rational::from_i64(64) * (&Rational::one() + &reach);
        let tolerance = &self.options.curve_tolerance / &(&metric * &Rational::from_i64(4));
        let budget = uv::Budget::new(pieces.len(), self.digits(), self.dec)?;
        let mapped = crate::observe("step.pipe-end-pcurve", &edge.index(), self.dec, |dec| {
            uv::fit(
                &pieces,
                self.options.significant_digits,
                &tolerance,
                budget,
                dec,
            )
        })?;
        let written = self.written_envelope(surface)?;
        self.deviation(&(&mapped.error * &metric) + &written.error);
        let points = mapped
            .points
            .iter()
            .map(|point| {
                self.graph.add(
                    "CARTESIAN_POINT",
                    vec![
                        text("''"),
                        text(format!("({},{})", point[0].text, point[1].text)),
                    ],
                )
            })
            .collect::<Vec<_>>();
        let knots = mapped
            .knots
            .iter()
            .map(|q| number::rational(q, self.digits()).text)
            .collect::<Vec<_>>()
            .join(",");
        let multiplicities = mapped
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
                boolean(false),
                boolean(false),
                text(format!("({multiplicities})")),
                text(format!("({knots})")),
                text(".UNSPECIFIED."),
            ],
        );
        let context = self.parameter_context();
        let representation = self.graph.add(
            "DEFINITIONAL_REPRESENTATION",
            vec![text("''"), refs([curve]), context.into()],
        );
        self.approximate_pipe_ends = true;
        Ok(self.graph.add(
            "PCURVE",
            vec![text("''"), written.id.into(), representation.into()],
        ))
    }
}
