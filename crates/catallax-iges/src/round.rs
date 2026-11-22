use catallax_numbers::{Ball, Rational, Real};
use catallax_shapes::{Circle, Surface, Vector3};
use catallax_topology::Face;

use super::{Exporter, IgesError, Value, ex, text};

impl Exporter<'_> {
    pub(super) fn round_face(&mut self, face: &Face) -> Result<usize, IgesError> {
        if !face.loops.is_empty() {
            return Err(IgesError::UnrepresentableDefinition(
                "sphere and torus IGES export requires a whole support",
            ));
        }

        let mut loop_fields = vec![4usize.into()];
        match self.model.geom().surface(face.surface).clone() {
            Surface::Sphere(s) => {
                let axis = unit(&s.axis, self.dec)?;
                let x = unit(&s.ref_dir, self.dec)?;
                let bottom = s.center.plus(&axis.scaled(&(-&s.radius)));
                let top = s.center.plus(&axis.scaled(&s.radius));
                let bottom = self.round_vertex(&bottom)?;
                let top = self.round_vertex(&top)?;
                let meridian = self.round_edge(
                    Circle {
                        center: s.center,
                        radius_sq: &s.radius * &s.radius,
                        normal: x.cross(&axis),
                        ref_dir: axis.scaled(&Real::from_i64(-1)),
                    },
                    bottom,
                    top,
                )?;
                let half = self.circle_span(1)?;
                let vertex_range = "1";
                // IGES 196 has azimuth and latitude in degrees. OCCT's B-rep
                // reader uses the vertex-list records for the degenerate poles.
                self.round_use(
                    &mut loop_fields,
                    true,
                    bottom,
                    true,
                    [[0, -90], [360, -90]],
                    vertex_range,
                );
                self.round_use(
                    &mut loop_fields,
                    false,
                    meridian,
                    true,
                    [[360, -90], [360, 90]],
                    &half,
                );
                self.round_use(
                    &mut loop_fields,
                    true,
                    top,
                    true,
                    [[360, 90], [0, 90]],
                    vertex_range,
                );
                self.round_use(
                    &mut loop_fields,
                    false,
                    meridian,
                    false,
                    [[0, 90], [0, -90]],
                    &half,
                );
            }
            Surface::Torus(t) => {
                let axis = unit(&t.axis, self.dec)?;
                let x = unit(&t.ref_dir, self.dec)?;
                let outer = &t.major + &t.minor;
                let joint = self.round_vertex(&t.center.plus(&x.scaled(&outer)))?;
                let equator = self.round_edge(
                    Circle {
                        center: t.center.clone(),
                        radius_sq: &outer * &outer,
                        normal: axis.clone(),
                        ref_dir: x.clone(),
                    },
                    joint,
                    joint,
                )?;
                let meridian = self.round_edge(
                    Circle {
                        center: t.center.plus(&x.scaled(&t.major)),
                        radius_sq: &t.minor * &t.minor,
                        normal: x.cross(&axis),
                        ref_dir: x,
                    },
                    joint,
                    joint,
                )?;
                let full = self.circle_span(2)?;
                // IGES 198 uses the meridian angle first and the reversed
                // azimuth second (OCCT TransferSurfaceAndUVBounds).
                self.round_use(
                    &mut loop_fields,
                    false,
                    equator,
                    true,
                    [[0, 360], [0, 0]],
                    &full,
                );
                self.round_use(
                    &mut loop_fields,
                    false,
                    meridian,
                    true,
                    [[0, 0], [360, 0]],
                    &full,
                );
                self.round_use(
                    &mut loop_fields,
                    false,
                    equator,
                    false,
                    [[360, 0], [360, 360]],
                    &full,
                );
                self.round_use(
                    &mut loop_fields,
                    false,
                    meridian,
                    false,
                    [[360, 360], [0, 360]],
                    &full,
                );
            }
            _ => unreachable!("caller selects a sphere or torus"),
        }
        let surface = self.surface(face.surface)?;
        let boundary = self.graph.add(508, 1, loop_fields);
        self.faces += 1;
        Ok(self.graph.add(
            510,
            1,
            vec![
                surface.into(),
                1usize.into(),
                1usize.into(),
                boundary.into(),
            ],
        ))
    }

    fn round_vertex(&mut self, p: &catallax_shapes::Point3) -> Result<usize, IgesError> {
        let point = ex::point(p, self.dec)?;
        self.added += 1;
        Ok(self.add_vertex(point))
    }

    fn round_edge(&mut self, circle: Circle, start: usize, end: usize) -> Result<usize, IgesError> {
        let curve = self.circle(&circle, start, end, true)?;
        let index = self.edge_rows.len() / 5 + 1;
        self.edge_rows.extend([
            curve.into(),
            self.vertex_list.into(),
            start.into(),
            self.vertex_list.into(),
            end.into(),
        ]);
        self.seams += 1;
        Ok(index)
    }

    fn circle_span(&mut self, half_turns: i64) -> Result<String, IgesError> {
        // Same range enclosure as the existing circular pcurves. The degree
        // endpoints are exact. Only the shared radian curve range is rounded.
        let span = Ball::pi_checked(160, self.dec)?
            * Ball::from_rational(&Rational::from_i64(half_turns), 160);
        let decimal = span.mid().to_decimal_bounds(17).0;
        let value = Rational::parse(&decimal).expect("decimal from the number library");
        let error = (&value - &span.lo()).abs().max((&value - &span.hi()).abs());
        self.parameter_error = self.parameter_error.clone().max(&self.bound * &error);
        ex::decimal(&value)
    }

    fn round_use(
        &mut self,
        fields: &mut Vec<Value>,
        pole: bool,
        index: usize,
        forward: bool,
        [[u0, v0], [u1, v1]]: [[i32; 2]; 2],
        range: &str,
    ) {
        // IGES 126 degree-one pcurve with unit weights, parameterized on
        // [0, range]. In 508 it follows the use direction, including reversals.
        let pcurve = self.graph.add(
            126,
            0,
            vec![
                text(1),
                text(1),
                text(1),
                text(0),
                text(1),
                text(0),
                text(0),
                text(0),
                text(range),
                text(range),
                text(1),
                text(1),
                text(u0),
                text(v0),
                text(0),
                text(u1),
                text(v1),
                text(0),
                text(0),
                text(range),
                text(0),
                text(0),
                text(1),
            ],
        );
        let list = if pole {
            self.vertex_list
        } else {
            self.edge_list
        };
        fields.extend([
            usize::from(pole).into(),
            list.into(),
            index.into(),
            usize::from(forward).into(),
            1usize.into(),
            1usize.into(),
            pcurve.into(),
        ]);
    }
}

fn unit(v: &Vector3, dec: &mut catallax_numbers::Decider) -> Result<Vector3, IgesError> {
    let [x, y, z] = ex::unit(v, dec)?;
    Ok(catallax_shapes::RVector3::new(x, y, z).to_real())
}
