//! Exact rectangular trims on a Bezier extrusion graph.

use super::*;
use catallax_numbers::RealPoly;
use catallax_recipe::{ElementId, Point2, Segment};
use catallax_shapes::{BezierExtrusion, Pose, Rotation};
use catallax_topology::BezierExtrusionBand;

type UV = [Rational; 2];

fn split(controls: &[RPoint3], t: &Rational) -> [Vec<RPoint3>; 2] {
    let mut row = controls.to_vec();
    let mut left = vec![row[0].clone()];
    let mut right = vec![row.last().unwrap().clone()];
    while row.len() > 1 {
        row = row
            .windows(2)
            .map(|p| p[0].plus(&p[1].minus(&p[0]).scaled(t)))
            .collect();
        left.push(row[0].clone());
        right.push(row.last().unwrap().clone());
    }
    right.reverse();
    [left, right]
}

fn parameter(
    values: &[Rational],
    station: &Rational,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let t = RealPoly::var(1, 0);
    let one_minus = RealPoly::constant(1, Real::one()) - t.clone();
    let mut row: Vec<_> = values
        .iter()
        .map(|v| RealPoly::constant(1, (v - station).into()))
        .collect();
    while row.len() > 1 {
        row = row
            .windows(2)
            .map(|p| &p[0] * &one_minus + &p[1] * &t)
            .collect();
    }
    for parameter in row[0].real_roots(dec)? {
        if sign(&parameter, dec)? != Sign::Negative
            && sign(&(&parameter - &Real::one()), dec)? != Sign::Positive
        {
            return dec.rational_value(&parameter)?.ok_or_else(|| {
                unsupported("the Bezier trim needs an algebraic section parameter")
            });
        }
    }
    Err(bad(
        "features.limit-missed",
        "the profile extends beyond the Bezier chord",
    ))
}

pub(super) fn build(
    model: &mut Model,
    face: FaceId,
    surface: &BezierExtrusion,
    profile: &Profile,
    direction: &RVector3,
    step: &StepId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let band = BezierExtrusionBand::of_face(model, face, dec)?
        .ok_or_else(|| unsupported("the spline limit has no certified finite band"))?;
    let controls: Vec<_> = surface
        .controls()
        .iter()
        .map(Point3::as_rational)
        .collect::<Option<_>>()
        .ok_or_else(|| unsupported("the Bezier limit needs a rational section chart"))?;
    let axis = surface
        .direction()
        .as_rational()
        .ok_or_else(|| unsupported("the Bezier limit needs rational generators"))?;
    let chord = controls.last().unwrap().minus(&controls[0]);
    let across = chord.cross(&axis);
    let normal = normal(profile);
    if !direction.cross(&across).is_zero() || !normal.dot(&axis).is_zero() {
        return Err(unsupported(
            "the Bezier limit needs a rectangular footprint in its chord chart",
        ));
    }
    let [boundary] = profile.loops.as_slice() else {
        return Err(unsupported("the Bezier limit needs one rectangular loop"));
    };
    if boundary.segments.len() != 4
        || boundary
            .segments
            .iter()
            .any(|s| s.circle.is_some() || s.ids.len() != 1)
    {
        return Err(unsupported(
            "the Bezier limit needs four straight authored sides",
        ));
    }
    let base = &controls[0];
    let chart = |p: &RPoint3| -> UV {
        let p = p.minus(base);
        [
            &p.dot(&chord) / &chord.norm_squared(),
            &p.dot(&axis) / &axis.norm_squared(),
        ]
    };
    let world = |uv| profile.world(uv, &RVector3::zero());
    let points: Vec<_> = boundary
        .segments
        .iter()
        .map(|s| chart(&world(&s.end)))
        .collect();
    let [[u0, u1], [v0, v1]] = std::array::from_fn(|axis| {
        [
            points.iter().map(|p| &p[axis]).min().unwrap().clone(),
            points.iter().map(|p| &p[axis]).max().unwrap().clone(),
        ]
    });
    let mut sides = BTreeMap::new();
    for segment in &boundary.segments {
        let [a, b] = [&segment.start, &segment.end].map(|p| chart(&world(p)));
        let side = match (a[0] == b[0], a[1] == b[1]) {
            (true, false) if a[0] == u0 => "u-min",
            (true, false) if a[0] == u1 => "u-max",
            (false, true) if a[1] == v0 => "v-min",
            (false, true) if a[1] == v1 => "v-max",
            _ => {
                return Err(unsupported(
                    "the profile is not rectangular in the Bezier chord chart",
                ));
            }
        };
        sides.insert(side, segment.ids[0].clone());
    }
    if sides.len() != 4 {
        return Err(unsupported(
            "the profile is not rectangular in the Bezier chord chart",
        ));
    }
    let [min, max] = band.axial();
    if u0 < Rational::zero()
        || u1 > Rational::one()
        || sign(&(&Real::from(v0.clone()) - min), dec)? == Sign::Negative
        || sign(&(&Real::from(v1.clone()) - max), dec)? == Sign::Positive
    {
        return Err(bad(
            "features.limit-missed",
            "the profile extends beyond the selected spline trim",
        ));
    }

    let stations: Vec<_> = controls.iter().map(|p| chart(p)[0].clone()).collect();
    let lo = parameter(&stations, &u0, dec)?;
    let hi = parameter(&stations, &u1, dec)?;
    let [left, _] = split(&controls, &hi);
    let [_, span] = split(&left, &(&lo / &hi));
    let shift = axis.scaled(&v0);
    let span: Vec<_> = span.iter().map(|p| p.plus(&shift)).collect();

    let sketch_origin = world(&[Rational::zero(), Rational::zero()]);
    let denominator = direction.dot(&normal);
    let height = |point: &RPoint3| point.minus(&sketch_origin).dot(&normal) / denominator.clone();
    let rim_heights = [height(&span[0]).sign(), height(span.last().unwrap()).sign()];
    if rim_heights.contains(&Sign::Negative) {
        return Err(bad(
            "features.limit-missed",
            "a spline cap endpoint lies behind the sketch",
        ));
    }
    if rim_heights.contains(&Sign::Zero) {
        return Err(unsupported("the Bezier cap has a trim edge at zero depth"));
    }
    if span.iter().any(|p| height(p).sign() != Sign::Positive) {
        return Err(unsupported(
            "the Bezier hull does not certify positive cap height",
        ));
    }
    let on_sketch = |p: &RPoint3| p.plus(&direction.scaled(&(-height(p))));
    let rotation = Rotation::from_frame(&axis.to_real(), &chord.to_real(), dec)?;
    if rotation.rational_matrix().is_none() {
        return Err(unsupported(
            "the Bezier section frame cannot be represented with rational profile coordinates",
        ));
    }
    let frame = Pose {
        rotation,
        translation: span[0].to_real().minus(&Point3::origin()),
    };
    let inverse = frame.inverse();
    let local = |point: &RPoint3| {
        let p = crate::profile::rational_point(&inverse, point);
        Point2([p.x.into(), p.y.into()])
    };
    let id = |s| ElementId::new(s).expect("fixed section side id");
    let first = &span[0];
    let last = span.last().unwrap();
    let mut section = Profile::new(
        &catallax_recipe::Pose::identity(),
        &[catallax_recipe::Loop(vec![
            Segment::Bezier {
                id: id("cap"),
                controls: span[1..span.len() - 1].iter().map(local).collect(),
                to: local(last),
            },
            Segment::Line {
                id: id("u-max"),
                to: local(&on_sketch(last)),
            },
            Segment::Line {
                id: id("sketch"),
                to: local(&on_sketch(first)),
            },
            Segment::Line {
                id: id("u-min"),
                to: local(first),
            },
        ])],
        dec,
    )?;
    section.pose = frame;
    let column = crate::extrude::extrude_vector(
        model,
        &section,
        step,
        &axis.scaled(&(&v1 - &v0)),
        [Some(Role::Start), Some(Role::End)],
        None,
        dec,
    )?;
    let side = |s: &str| Role::Side {
        segment: sides[s].clone(),
        path: None,
    };
    Ok(transformed_body(
        model,
        column,
        &Isometry::identity(),
        |origin| {
            let role = match &origin.role {
                Role::Start => side("v-min"),
                Role::End => side("v-max"),
                Role::Side { segment, .. } => match segment.as_str() {
                    "cap" => Role::End,
                    "sketch" => Role::Start,
                    name => side(name),
                },
                _ => unreachable!("section sweep has only caps and sides"),
            };
            Origin {
                step: origin.step.clone(),
                role,
            }
        },
    ))
}
