//! Whole spline charts whose boundary stations must be shared with their neighbours.

use std::collections::{BTreeMap, btree_map::Entry};

use catallax_numbers::{Ball, Decider, Interval, Outcome, Rational, Real, Refusal, Sign};
use catallax_shapes::{BSplineCurve, BoundaryRow, SplineBoundary, SplineSupport, Surface};
use catallax_topology::{EdgeId, FaceId, Model, RollingMitreBand, Sense};

use super::{BITS, Chart, SplineChart, Trim, Uv, q};
use crate::{MeshOptions, failure, junction::Seed};

pub(crate) struct FixedTrim {
    pub(crate) axes: [Vec<f64>; 2],
    pub(crate) sense: Sense,
    pub(crate) apex: bool,
    pub(super) rows: BTreeMap<EdgeId, Trim>,
}

impl FixedTrim {
    pub(super) fn new(
        model: &Model,
        face: FaceId,
        support: &SplineSupport,
        spline: &SplineChart,
        options: &MeshOptions,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        let rows = match support {
            SplineSupport::RollingMitre(_) => RollingMitreBand::checked(model, face, dec)?
                .map_err(|undecided| {
                    failure("rolling mitre has no checked trim").keeping_any(undecided)
                })?
                .rows()
                .collect(),
            SplineSupport::RollingEnvelope(envelope) => {
                if !catallax_topology::envelope_retains_chart(model, face, dec)? {
                    return Err(failure("rolling envelope has no checked complete trim"));
                }
                let low = envelope.definition().corners(Sign::Negative);
                let high = envelope.definition().corners(Sign::Positive);
                catallax_topology::envelope_rows(model, face, dec)?
                    .into_iter()
                    .map(|(edge, row, sense)| {
                        let first = match row {
                            BoundaryRow::Rail | BoundaryRow::Negative => &low[0],
                            BoundaryRow::Mouth => &low[1],
                            BoundaryRow::Positive => &high[0],
                        };
                        let start = model.edge(edge).start.expect("checked envelope endpoint");
                        let start = model.geom().position(model.vertex(start).point);
                        let forward = match dec.sign(&start.minus(first).norm_squared()) {
                            Outcome::ProvenCoincident(_) => true,
                            Outcome::Proven(_) => false,
                            Outcome::Refused(error) => return Err(error),
                            Outcome::NotExact(error) => return Err(error.into()),
                        };
                        Ok((edge, row, sense, forward))
                    })
                    .collect::<Result<Vec<_>, Refusal>>()?
            }
            _ => return Ok(None),
        };
        let axes = spline.fixed_lattice(
            &options.chord_tolerance.mul_pow2(-1),
            options.max_vertices,
            dec,
        )?;
        let bounds = axes.each_ref().map(|axis| [axis[0], axis[axis.len() - 1]]);
        let mut chart_sense = None;
        let rows = rows
            .into_iter()
            .map(|(edge, row, sense, forward)| {
                let fixed = match row {
                    BoundaryRow::Rail => (1, bounds[1][0]),
                    BoundaryRow::Mouth => (1, bounds[1][1]),
                    BoundaryRow::Negative => (0, bounds[0][0]),
                    BoundaryRow::Positive => (0, bounds[0][1]),
                };
                if row == BoundaryRow::Rail {
                    chart_sense = Some(if (sense == Sense::Same) == forward {
                        Sense::Same
                    } else {
                        Sense::Reversed
                    });
                }
                let mut ends = [bounds.map(|axis| axis[0]), bounds.map(|axis| axis[1])];
                ends.iter_mut().for_each(|uv| uv[fixed.0] = fixed.1);
                if !forward {
                    ends.reverse();
                }
                (edge, Trim { ends, fixed })
            })
            .collect();
        Ok(Some(Self {
            axes,
            sense: chart_sense.ok_or_else(|| failure("a fixed trim has no lower rail"))?,
            apex: matches!(support, SplineSupport::RollingMitre(_)),
            rows,
        }))
    }

    pub(crate) fn points(
        &self,
        trim: Trim,
    ) -> impl DoubleEndedIterator<Item = Uv> + ExactSizeIterator + '_ {
        let (axis, value) = trim.fixed;
        self.axes[1 - axis].iter().map(move |&x| {
            let mut uv = [value; 2];
            uv[1 - axis] = x;
            uv
        })
    }

    pub(crate) fn seed_edges(
        &self,
        model: &Model,
        chart: &Chart,
        seeds: &mut BTreeMap<EdgeId, Seed>,
        options: &MeshOptions,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let Surface::BSpline(support) = &chart.surface else {
            unreachable!()
        };
        let spline = chart.spline().expect("fixed spline trim");
        for (&edge, &trim) in &self.rows {
            let Entry::Vacant(slot) = seeds.entry(edge) else {
                return Err(failure("two fixed charts prescribe the same edge"));
            };
            let (axis, value) = trim.fixed;
            let points = self.points(trim);
            let mut exact = Vec::new();
            if crate::junction::meets_plane(model, edge) {
                let side = match (axis, value == self.axes[axis][0]) {
                    (0, true) => SplineBoundary::UMin,
                    (0, false) => SplineBoundary::UMax,
                    (1, true) => SplineBoundary::VMin,
                    _ => SplineBoundary::VMax,
                };
                let curve = BSplineCurve::boundary(Chart::carrier(support), side);
                exact.reserve(points.len());
                for uv in self.points(trim) {
                    let p = if let SplineSupport::RollingMitre(mitre) = support {
                        mitre.point(&q(uv[0]), &q(uv[1]))
                    } else {
                        curve.point_at(&Real::from(q(uv[1 - axis])), dec)?
                    };
                    exact.push([p.x, p.y, p.z]);
                }
            }
            let mut error = Rational::zero();
            let mut previous = None;
            let mut enclosed = Vec::with_capacity(points.len());
            for uv in points {
                if let Some(a) = previous {
                    let bound = spline.chord(&[a, a, uv], dec)?;
                    if !bound.is_finite()
                        || bound.hi() > Interval::from_rational(&options.chord_tolerance).lo()
                    {
                        return Err(failure("fixed rail has no chord certificate"));
                    }
                    error = error.max(q(bound.hi()));
                }
                previous = Some(uv);
                let xyz = spline.evaluate(uv, dec)?;
                if !xyz.iter().all(Interval::is_finite) {
                    return Err(failure("fixed rail exceeds binary64"));
                }
                enclosed.push(xyz.map(|x| {
                    Ball::from_rational(&q(x.lo()), BITS)
                        .hull(&Ball::from_rational(&q(x.hi()), BITS))
                }));
            }
            slot.insert(Seed {
                points: enclosed,
                error,
                forward: trim.ends[0][1 - axis] < trim.ends[1][1 - axis],
                exact,
            });
        }
        Ok(())
    }
}
