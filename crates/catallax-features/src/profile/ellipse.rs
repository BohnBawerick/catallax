//! Exact affine reduction of co-oriented, equal-aspect ellipses to circles.

use super::{P, bad, dot, point};
use catallax_numbers::{Rational, Refusal, Sign};
use catallax_recipe::{Loop, Number, Point2, Segment};

#[derive(Clone, Debug)]
pub(crate) struct Frame {
    pub axis: P,
    pub ratio: Rational,
}

impl Frame {
    pub fn new(axis: P, ratio: Rational) -> Result<Self, Refusal> {
        if dot(&axis, &axis).is_zero() || ratio.sign() != Sign::Positive {
            return Err(bad(
                "sketch.ellipse-axes",
                "ellipse axes must have positive lengths",
            ));
        }
        Ok(Self { axis, ratio })
    }

    pub fn forward(&self, p: &P) -> P {
        [
            &self.axis[0] * &p[0] - &(&self.axis[1] * &self.ratio) * &p[1],
            &self.axis[1] * &p[0] + &(&self.axis[0] * &self.ratio) * &p[1],
        ]
    }

    pub fn inverse(&self, p: &P) -> P {
        let n = dot(&self.axis, &self.axis);
        [
            &(&self.axis[0] * &p[0] + &self.axis[1] * &p[1]) / &n,
            (&self.axis[0] * &p[1] - &self.axis[1] * &p[0]) / (&n * &self.ratio),
        ]
    }

    fn mapped(&self, p: &Point2) -> Point2 {
        Point2(self.inverse(&point(p)).map(Number::new))
    }

    fn radius(&self, axis: &Point2, ratio: &Number) -> Result<Rational, Refusal> {
        let a = self.inverse(&point(axis));
        if ratio.value() != &self.ratio || !a[1].is_zero() {
            return Err(bad(
                "sketch.ellipse-frame",
                "ellipses in one profile must have parallel first axes and equal aspect ratios",
            ));
        }
        if a[0].is_zero() {
            return Err(bad(
                "sketch.ellipse-axes",
                "ellipse axes must have positive lengths",
            ));
        }
        Ok(a[0].abs())
    }
}

pub(crate) fn normalize(loops: &[Loop]) -> Result<Option<(Frame, Vec<Loop>)>, Refusal> {
    let Some((axis, ratio)) = loops.iter().flat_map(|lp| &lp.0).find_map(|s| match s {
        Segment::Ellipse { axis, ratio, .. } | Segment::EllipticArc { axis, ratio, .. } => {
            Some((point(axis), ratio.value().clone()))
        }
        _ => None,
    }) else {
        return Ok(None);
    };
    let frame = Frame::new(axis, ratio)?;
    let loops = loops
        .iter()
        .map(|lp| {
            let segments = lp
                .0
                .iter()
                .map(|s| {
                    Ok(match s {
                        Segment::Line { id, to } => Segment::Line {
                            id: id.clone(),
                            to: frame.mapped(to),
                        },
                        Segment::Ellipse {
                            id,
                            centre,
                            axis,
                            ratio,
                        } => Segment::Circle {
                            id: id.clone(),
                            centre: frame.mapped(centre),
                            radius: Number::new(frame.radius(axis, ratio)?),
                        },
                        Segment::EllipticArc {
                            id,
                            centre,
                            axis,
                            ratio,
                            to,
                            turn,
                        } => {
                            let radius = frame.radius(axis, ratio)?;
                            let c = frame.mapped(centre);
                            let end = frame.mapped(to);
                            let d = super::sub(&point(&end), &point(&c));
                            if dot(&d, &d) != &radius * &radius {
                                return Err(bad(
                                    "sketch.ellipse-endpoint",
                                    "arc endpoint must lie exactly on its ellipse",
                                ));
                            }
                            Segment::Arc {
                                id: id.clone(),
                                centre: c,
                                to: end,
                                turn: *turn,
                            }
                        }
                        Segment::Bezier { .. } | Segment::Arc { .. } | Segment::Circle { .. } => {
                            return Err(bad(
                                "sketch.ellipse-circular-mix",
                                "one profile cannot yet mix ellipses and circular curves",
                            ));
                        }
                    })
                })
                .collect::<Result<_, Refusal>>()?;
            Ok(Loop(segments))
        })
        .collect::<Result<_, Refusal>>()?;
    Ok(Some((frame, loops)))
}
