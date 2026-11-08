//! Exact bands over a mitre. The spine is `u`, the section is `v`. Each rail
//! is written as its exact ellipse and the end section as its exact circle,
//! which the neighbouring cylinders chart for themselves. The apex column has
//! a constant 3D control row and pcurve at the existing apex vertex, but no
//! additional model edge.
//!
//! On the band, each conic's chart curve runs on the conic's own angle `x`
//! from its foot. A rail has `u = tan(x / 2)`; the section has
//! `v = tan(x / 2) / l` for `l = across / along`, over `x` up to `2 atan l`.
//! Both are `g = tan(x / 2) / l` with `T = tan(x / 2) <= l` (`l = 1` on a
//! rail), so `|g''''| = T (1 + T^2) (2 + 3 T^2) / (2 l)`, which grows with
//! `T` and stays below `(1 + l^2) (2 + 3 l^2) / 2`. Cubic Hermite spans of
//! width `h` then miss by at most `h^4 / 384` times that.

use super::*;
use crate::graph::Field;
use catallax_numbers::{Ball, Outcome, StrictSign};
use catallax_shapes::{BoundaryRow, SplineSupport};
use catallax_topology::RollingMitreBand;

const BITS: u32 = 192;
// The chart curves' target miss, in chart units.
const CHART_MISS: f64 = 1e-10;

#[derive(Clone, Copy)]
pub(super) struct Row {
    surface: SurfaceId,
    boundary: BoundaryRow,
    forward: bool,
}

// The degree-(6, 2) carrier is `u`-major: control `(i, j)` is `3 i + j`.
const SPINE: usize = 6;

fn row_indices(boundary: BoundaryRow) -> (Vec<usize>, usize) {
    match boundary {
        BoundaryRow::Rail => ((0..=SPINE).map(|i| 3 * i).collect(), SPINE),
        BoundaryRow::Mouth => ((0..=SPINE).map(|i| 3 * i + 2).collect(), SPINE),
        BoundaryRow::Negative => (vec![0, 1, 2], 2),
        BoundaryRow::Positive => ((3 * SPINE..3 * SPINE + 3).collect(), 2),
    }
}

impl Exporter<'_> {
    pub(super) fn plan_mitre_bands(
        &mut self,
        bodies: impl Iterator<Item = BodyId>,
    ) -> Result<(), Refusal> {
        for body in bodies {
            for face in self.model.faces_of_body(body) {
                let Some(band) = RollingMitreBand::of_face(self.model, face, self.dec)? else {
                    continue;
                };
                let surface = self.model.face(face).surface;
                let mut rows = [None; 3];
                for (index, (edge, boundary, _, forward)) in band.rows().enumerate() {
                    rows[index] = Some(edge);
                    let row = Row {
                        surface,
                        boundary,
                        forward,
                    };
                    if let Some(previous) = self.mitre_rows.insert(edge, row)
                        && (previous.surface != row.surface
                            || previous.boundary != row.boundary
                            || previous.forward != row.forward)
                    {
                        return Err(Refusal::new(
                            "step.unsupported-geometry",
                            "two mitre bands share a boundary outside the proved disjoint domain",
                        ));
                    }
                }
                let [Some(a), Some(b), Some(c)] = rows else {
                    return Err(unsupported());
                };
                self.mitre_faces.insert(face, [a, b, c]);
            }
        }
        Ok(())
    }

    // The unit-interval parameter is also the 3D control row's parameter.
    fn mitre_pcurve(&mut self, surface: Id, row: BoundaryRow) -> Id {
        let coordinates = match row {
            BoundaryRow::Rail => ["(0.,0.)", "(1.,0.)"],
            BoundaryRow::Mouth => ["(0.,1.)", "(1.,1.)"],
            BoundaryRow::Negative => ["(0.,0.)", "(0.,1.)"],
            BoundaryRow::Positive => ["(1.,0.)", "(1.,1.)"],
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

    pub(super) fn mitre_edge(
        &mut self,
        id: EdgeId,
        topology: &mut Topology,
    ) -> Result<Rc<Edge>, Refusal> {
        let row = self.mitre_rows[&id];
        let Surface::BSpline(SplineSupport::RollingMitre(band)) =
            self.model.geom().surface(row.surface).clone()
        else {
            return Err(unsupported());
        };
        let written = self.written_envelope(row.surface)?;
        let source = self.model.edge(id);
        let (Some(from), Some(to)) = (source.start, source.end) else {
            return Err(unsupported());
        };
        let start = self.vertex(from, topology)?;
        let end = self.vertex(to, topology)?;
        let quarter = Ball::pi_checked(BITS, self.dec)?
            * Ball::from_rational(&Rational::from_ratio(1, 2).expect("nonzero"), BITS);
        // A charted neighbour charts the conic as written here.
        let capture = self.has_chart_edge(id, topology);
        let (conic, form, first, spread, lambda, reach) = match row.boundary {
            BoundaryRow::Rail | BoundaryRow::Mouth => {
                let rail = band.rail(usize::from(row.boundary == BoundaryRow::Mouth), self.dec)?;
                let up = band.frame()[1].clone();
                let axis = rail.major_dir.cross(&up);
                let major = number::sqrt_real(&rail.major_sq);
                let minor = number::sqrt_real(&rail.minor_sq);
                let reach = &(&major + &minor) * &Real::from_i64(2);
                let along = match self.dec.sign(&(&rail.major_sq - &rail.minor_sq)) {
                    Outcome::Proven(StrictSign::Negative) => false,
                    Outcome::Refused(refusal) => return Err(refusal),
                    _ => true,
                };
                // With `y = axis x U` along `up`, the angle is `t` itself; the
                // other way round the chart starts a quarter turn on.
                let (reference, major_sq, minor_sq, first) = if along {
                    (
                        rail.major_dir.clone(),
                        rail.major_sq.clone(),
                        rail.minor_sq.clone(),
                        Ball::from_rational(&Rational::zero(), BITS),
                    )
                } else {
                    (
                        up,
                        rail.minor_sq.clone(),
                        rail.major_sq.clone(),
                        quarter.clone(),
                    )
                };
                let major = number::real(&number::sqrt_real(&major_sq), self.digits(), self.dec)?;
                let minor = number::real(&number::sqrt_real(&minor_sq), self.digits(), self.dec)?;
                let conic =
                    self.ellipse(&rail.center, &axis, &reference, major, Some(minor), capture)?;
                let form = CurveForm::Ellipse {
                    center: rail.center.clone(),
                    normal: axis,
                    major_dir: reference,
                    major_sq,
                    minor_sq,
                };
                let one = Ball::from_rational(&Rational::one(), BITS);
                (conic, form, first, quarter, one, reach)
            }
            BoundaryRow::Negative => {
                let section = band.section(self.dec)?;
                let radius = number::sqrt_real(&section.radius_sq);
                let lambda = self.dec.enclose_at(&(band.across() / band.along()), BITS)?;
                let reach = &(&radius * &Real::from_i64(2)) * &(band.across() / band.along());
                let spread = lambda.atan_checked(self.dec)?
                    * Ball::from_rational(&Rational::from_i64(2), BITS);
                let radius = number::real(&radius, self.digits(), self.dec)?;
                let conic = self.ellipse(
                    &section.center,
                    &section.normal,
                    &section.ref_dir,
                    radius,
                    None,
                    capture,
                )?;
                let form = CurveForm::Circle {
                    center: section.center.clone(),
                    normal: section.normal.clone(),
                    radius_sq: section.radius_sq.clone(),
                    ref_dir: section.ref_dir.clone(),
                };
                (
                    conic,
                    form,
                    Ball::from_rational(&Rational::zero(), BITS),
                    spread,
                    lambda,
                    reach,
                )
            }
            BoundaryRow::Positive => return Err(unsupported()),
        };
        let (pcurve, miss) =
            self.mitre_chart(written.id, row.boundary, &first, &spread, &lambda)?;
        let reach = number::upper(&reach, self.dec)?;
        self.deviation(&(&miss * &reach) + &written.error);
        let error = arc_error(&conic.curve, &start, &end);
        let bounds = conic.curve.bounds.clone().ok_or_else(unsupported)?;
        let bounded = match conic.written.filter(|_| capture) {
            Some(written) => {
                let geometry = charts::Geometry::Analytic(Rc::new(charts::Analytic {
                    written,
                    exact: form,
                }));
                self.chart_edge(
                    id,
                    charts::EdgeInput {
                        curve: conic.curve.id,
                        geometry: &geometry,
                        bounds: &bounds,
                        error: &error,
                        same: row.forward,
                        extra: Some(pcurve),
                    },
                    topology,
                )?
            }
            None => {
                self.deviation(error);
                self.graph.add(
                    "SURFACE_CURVE",
                    vec![
                        text("''"),
                        conic.curve.id.into(),
                        refs([pcurve]),
                        text(".CURVE_3D."),
                    ],
                )
            }
        };
        let edge = Rc::new(Edge {
            id: self.edge_entity(start.id, end.id, bounded, row.forward),
            bounds,
        });
        topology.edges.insert(id, edge.clone());
        Ok(edge)
    }

    // The chart curve `g = tan(x / 2) / lambda` of a conic's angle `x`, from
    // `first` over `spread`, as cubic Hermite spans; and its chart miss.
    fn mitre_chart(
        &mut self,
        surface: Id,
        boundary: BoundaryRow,
        first: &Ball,
        spread: &Ball,
        lambda: &Ball,
    ) -> Result<(Id, Rational), Refusal> {
        let exact = |value: i64| Ball::from_rational(&Rational::from_i64(value), BITS);
        let lambda_hi = lambda.hi();
        let lambda_sq = &lambda_hi * &lambda_hi;
        let fourth = &(&(&Rational::one() + &lambda_sq)
            * &(&Rational::from_i64(2) + &(&Rational::from_i64(3) * &lambda_sq)))
            * &Rational::from_ratio(1, 2).expect("nonzero");
        let estimate = spread.hi().to_f64_nearest()
            * (fourth.to_f64_nearest() / (384.0 * CHART_MISS)).powf(0.25);
        let spans = (estimate.ceil() as i64).clamp(4, 4096);
        let width_hi = &spread.hi() * &Rational::from_ratio(1, spans).expect("nonzero");
        let width_sq = &width_hi * &width_hi;
        let mut miss =
            &(&(&width_sq * &width_sq) * &fourth) * &Rational::from_ratio(1, 384).expect("nonzero");
        let width =
            spread * &Ball::from_rational(&Rational::from_ratio(1, spans).expect("nonzero"), BITS);
        let third = Ball::from_rational(&Rational::from_ratio(1, 3).expect("nonzero"), BITS);
        let half = Ball::from_rational(&Rational::from_ratio(1, 2).expect("nonzero"), BITS);
        let inverse = lambda.recip().ok_or_else(unsupported)?;
        // Value and slope at each node; the ends are exactly 0 and 1.
        let mut nodes = Vec::with_capacity(spans as usize + 1);
        for k in 0..=spans {
            let at = spread
                * &Ball::from_rational(&Rational::from_ratio(k, spans).expect("nonzero"), BITS);
            let angle = &at * &half;
            let tangent = angle
                .sin()
                .checked_div(&angle.cos())
                .ok_or_else(unsupported)?;
            let value = &tangent * &inverse;
            let slope = &(&(&tangent.square() + &exact(1)) * &half) * &inverse;
            let value = match k {
                0 => Ball::from_rational(&Rational::zero(), BITS),
                _ if k == spans => Ball::from_rational(&Rational::one(), BITS),
                _ => value,
            };
            nodes.push((value, slope));
        }
        let digits = self.digits();
        let fixed = match boundary {
            BoundaryRow::Rail => Some("0."),
            BoundaryRow::Mouth => Some("1."),
            _ => None,
        };
        let mut points = Vec::with_capacity(3 * spans as usize + 1);
        let mut rounding = Rational::zero();
        for k in 0..spans as usize {
            let (a, da) = &nodes[k];
            let (b, db) = &nodes[k + 1];
            let step = &width * &third;
            let bezier = [a.clone(), a + &(&step * da), b - &(&step * db), b.clone()];
            for control in &bezier[usize::from(k != 0)..] {
                let value = number::rational(&control.mid(), digits);
                rounding = rounding.max(&value.error + &control.rad());
                let text_value = value.text;
                let point = match fixed {
                    Some(v) => format!("({text_value},{v})"),
                    None => format!("(0.,{text_value})"),
                };
                points.push(
                    self.graph
                        .add("CARTESIAN_POINT", vec![text("''"), text(point)]),
                );
            }
        }
        // Writing knots as decimals moves each by at most their rounding, and
        // `g` changes by at most its largest slope per unit of angle.
        let mut knots = Vec::with_capacity(spans as usize + 1);
        let mut shift = Rational::zero();
        for k in 0..=spans {
            let knot = first
                + &(spread
                    * &Ball::from_rational(
                        &Rational::from_ratio(k, spans).expect("nonzero"),
                        BITS,
                    ));
            let value = number::rational(&knot.mid(), digits);
            shift = shift.max(&value.error + &knot.rad());
            knots.push(value.text);
        }
        let steepest = &(&Rational::one() + &lambda_sq)
            * &(&Rational::from_ratio(1, 2).expect("nonzero")
                * &lambda.recip().ok_or_else(unsupported)?.hi());
        miss = &(&miss + &rounding) + &(&shift * &steepest);
        let mut multiplicities = vec!["4".to_string()];
        multiplicities.extend((1..spans).map(|_| "3".to_string()));
        multiplicities.push("4".to_string());
        let curve = self.graph.add(
            "B_SPLINE_CURVE_WITH_KNOTS",
            vec![
                text("''"),
                text("3"),
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
        Ok((
            self.graph.add(
                "PCURVE",
                vec![text("''"), surface.into(), representation.into()],
            ),
            miss,
        ))
    }

    // The apex column: a constant row closing the parameter boundary at the
    // rails' shared end, with the same unit parameter as its pcurve.
    fn mitre_pole(
        &mut self,
        rail: EdgeId,
        topology: &mut Topology,
    ) -> Result<(Id, Bounds), Refusal> {
        let plan = self.mitre_rows[&rail];
        let written = self.written_envelope(plan.surface)?;
        let source = self.model.edge(rail);
        let vertex = if plan.forward {
            source.end
        } else {
            source.start
        }
        .ok_or_else(unsupported)?;
        let vertex = self.vertex(vertex, topology)?;
        let (indices, _) = row_indices(BoundaryRow::Positive);
        let curve = self.spline_row(
            &written,
            &indices,
            2,
            Field::List(vec![text("3"), text("3")]),
            Field::List(vec![text("0."), text("1.")]),
            true,
        );
        let pcurve = self.mitre_pcurve(written.id, BoundaryRow::Positive);
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

    // One boundary row in the walk, along its row's parameter or against it.
    fn mitre_use(
        &mut self,
        id: EdgeId,
        along: bool,
        bounds: &mut Option<Bounds>,
        topology: &mut Topology,
    ) -> Result<Id, Refusal> {
        let forward = self.mitre_rows[&id].forward;
        let edge = self.edge(id, topology)?;
        include_bounds(bounds, &edge.bounds);
        Ok(self.oriented(edge.id, forward == along))
    }

    /// A counter-clockwise walk of the complete parameter boundary: the first
    /// rail, the apex column, the second rail back and the end section back.
    pub(super) fn mitre_loop(
        &mut self,
        [rail, mouth, section]: [EdgeId; 3],
        topology: &mut Topology,
    ) -> Result<(Id, Bounds), Refusal> {
        let mut bounds = None;
        let mut uses = Vec::with_capacity(4);
        uses.push(self.mitre_use(rail, true, &mut bounds, topology)?);
        let (pole, extent) = self.mitre_pole(rail, topology)?;
        include_bounds(&mut bounds, &extent);
        uses.push(self.oriented(pole, true));
        uses.push(self.mitre_use(mouth, false, &mut bounds, topology)?);
        uses.push(self.mitre_use(section, false, &mut bounds, topology)?);
        Ok((
            self.graph.add("EDGE_LOOP", vec![text("''"), refs(uses)]),
            bounds.ok_or_else(unsupported)?,
        ))
    }
}
