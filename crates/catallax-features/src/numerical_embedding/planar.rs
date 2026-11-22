//! Source face boundaries have no outer designation; numerical faces do.
//! Area only proposes their order. Admission still verifies containment.

use super::*;

fn unit(vector: &catallax_shapes::Vector3, dec: &Decider) -> Result<[f64; 3], Refusal> {
    let values = [
        scalar(&vector.x, dec)?,
        scalar(&vector.y, dec)?,
        scalar(&vector.z, dec)?,
    ]
    .map(|v| v.to_f64_nearest());
    let length = libm::hypot(libm::hypot(values[0], values[1]), values[2]);
    if !length.is_finite() || length <= 0.0 {
        return Err(fit_failure(
            "a planar boundary has no finite projection normal",
        ));
    }
    Ok(values.map(|v| v / length))
}

fn area(
    face: &NumericalFace,
    boundary: &[NumericalTrim],
    edges: &[NumericalEdge],
    normal: [f64; 3],
    dec: &mut Decider,
) -> Result<f64, Refusal> {
    let Some(first) = boundary.first() else {
        return Ok(0.0);
    };
    let origin = &first.uv[0];
    let relative = |uv: &[Rational; 2]| {
        std::array::from_fn::<_, 2, _>(|i| (&uv[i] - &origin[i]).to_f64_nearest())
    };
    let cross = |a: [f64; 2], b: [f64; 2]| a[0] * b[1] - a[1] * b[0];
    let mut area = 0.0;
    for trim in boundary {
        dec.charge_work(1)?;
        let edge = &edges[trim.edge];
        let mut parameters = edge.parameters.each_ref();
        if trim.sense == Sense::Reversed {
            parameters.reverse();
        }
        let a = relative(&trim.uv[0]);
        let b = relative(&trim.uv[1]);
        match &edge.curve {
            Curve::Circle(circle) => {
                let angle = (parameters[1] - parameters[0]).to_f64_nearest();
                let projection = dot(normal, unit(&circle.normal, dec)?);
                let radius_sq = scalar(&circle.radius_sq, dec)?.to_f64_nearest();
                // Green's integral for the projected arc, including its bulge.
                area += cross(a, b) + radius_sq * projection * (angle - libm::sin(angle));
            }
            Curve::BSpline(curve) => {
                let forward = parameters[0] < parameters[1];
                let [lo, hi] = if forward {
                    parameters
                } else {
                    [parameters[1], parameters[0]]
                };
                let mut previous = a;
                let mut visit = |span: &[Rational]| -> Result<(), Refusal> {
                    let low = (&span[0]).max(lo);
                    let high = (&span[1]).min(hi);
                    if low >= high {
                        return Ok(());
                    }
                    let [from, to] = if forward { [low, high] } else { [high, low] };
                    // Ranking uses eighths of each retained knot span, not an
                    // endpoint chord or a new geometric approval tolerance.
                    for i in 1..=8 {
                        dec.charge_work(1)?;
                        let t = from + &((to - from) * Rational::from_ratio(i, 8).unwrap());
                        let next = if &t == parameters[1] {
                            b
                        } else {
                            let point = curve.evaluate_box([&t, &t], BITS, dec)?;
                            let uv = face.surface.parameters_of(&point, dec)?.ok_or_else(|| {
                                fit_failure("a planar boundary has no finite parameter projection")
                            })?;
                            relative(&uv.map(|v| v.mid()))
                        };
                        area += cross(previous, next);
                        previous = next;
                    }
                    Ok(())
                };
                let spans = curve.knots().windows(2);
                if forward {
                    for span in spans {
                        visit(span)?;
                    }
                } else {
                    for span in spans.rev() {
                        visit(span)?;
                    }
                }
            }
            _ => area += cross(a, b),
        }
    }
    if !area.is_finite() {
        return Err(fit_failure("a planar boundary has no finite area proposal"));
    }
    Ok(area)
}

pub(super) fn outer_first(
    face: &mut NumericalFace,
    edges: &[NumericalEdge],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let Surface::Plane(plane) = &face.surface else {
        return Ok(());
    };
    if face.loops.len() < 2 {
        return Ok(());
    }
    let normal = unit(&plane.normal, dec)?;
    let mut largest = (0, -1.0);
    for (index, boundary) in face.loops.iter().enumerate() {
        let magnitude = area(face, boundary, edges, normal, dec)?.abs();
        if magnitude > largest.1 {
            largest = (index, magnitude);
        }
    }
    face.loops.swap(0, largest.0);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::{RPoint3, RVector3};

    #[test]
    fn split_arc_bulges_determine_the_outer_loop_instead_of_endpoint_chords() {
        let q = Rational::from_i64;
        let center = RPoint3::new(q(0), q(0), q(0)).to_real();
        let normal = RVector3::new(q(0), q(0), q(1)).to_real();
        let reference = RVector3::new(q(1), q(0), q(0)).to_real();
        let circle = |radius: i64| {
            Curve::Circle(catallax_shapes::Circle {
                center: center.clone(),
                normal: normal.clone(),
                ref_dir: reference.clone(),
                radius_sq: q(radius * radius).into(),
            })
        };
        let half = Ball::pi(BITS).mid();
        let turn = half.mul_pow2(1);
        let edges = vec![
            NumericalEdge {
                key: "inner".into(),
                curve: circle(1),
                parameters: [q(0), turn.clone()],
                vertices: [0, 0],
            },
            NumericalEdge {
                key: "outer-first-half".into(),
                curve: circle(3),
                parameters: [q(0), half.clone()],
                vertices: [1, 2],
            },
            NumericalEdge {
                key: "outer-second-half".into(),
                curve: circle(3),
                parameters: [half, turn],
                vertices: [2, 1],
            },
        ];
        let point = |x: i64| [q(x), q(0)];
        let mut face = NumericalFace {
            key: "unordered-annulus".into(),
            surface: Surface::Plane(Plane {
                origin: center,
                normal,
                ref_dir: reference,
            }),
            sense: Sense::Same,
            loops: vec![
                vec![NumericalTrim {
                    tag: catallax_topology::numerical::CopyTag::Independent,
                    edge: 0,
                    sense: Sense::Reversed,
                    uv: [point(1), point(1)],
                }],
                vec![
                    NumericalTrim {
                        tag: catallax_topology::numerical::CopyTag::Independent,
                        edge: 1,
                        sense: Sense::Same,
                        uv: [point(3), point(-3)],
                    },
                    NumericalTrim {
                        tag: catallax_topology::numerical::CopyTag::Independent,
                        edge: 2,
                        sense: Sense::Same,
                        uv: [point(-3), point(3)],
                    },
                ],
            ],
        };
        let mut dec = Decider::default();
        outer_first(&mut face, &edges, &mut dec).unwrap();
        let outer = area(&face, &face.loops[0], &edges, [0.0, 0.0, 1.0], &mut dec).unwrap();
        let inner = area(&face, &face.loops[1], &edges, [0.0, 0.0, 1.0], &mut dec).unwrap();
        assert!((outer - 18.0 * std::f64::consts::PI).abs() < 1e-12);
        assert!((inner + 2.0 * std::f64::consts::PI).abs() < 1e-12);
    }
}
