//! Whole-disk spherical first-root certificate, including axis and closed rim.

use super::*;
use catallax_recipe::{Number, Point2};

#[derive(Clone, Debug)]
pub(super) struct Circle {
    pub radius: Rational,
    // Unit normal oriented along the raw ray, independent of its magnitude.
    pub normal: Vector3,
}

pub(super) fn compile(
    source: Scene,
    plane: &catallax_recipe::Pose,
    id: &ElementId,
    centre: &Point2,
    radius: &Number,
    direction: Option<&catallax_recipe::Vector3>,
    dec: &mut Decider,
) -> Result<AuthoredCap, Refusal> {
    dec.check_deadline()?;
    let Surface::Sphere(sphere) = &source.surface else {
        return Err(coverage(
            "circular first-hit profiles require a stock sphere",
        ));
    };
    let pose = scene::placement(plane)?;
    let u = pose.apply_vector(&axis(0));
    let v = pose.apply_vector(&axis(1));
    let normal = u.cross(&v);
    let direction = direction
        .map(|w| {
            let [x, y, z] = w.0.each_ref().map(|n| Real::from(n.value().clone()));
            Vector3::new(x, y, z)
        })
        .unwrap_or_else(|| normal.clone());
    if sign(&direction.norm_squared(), dec)? == Sign::Zero {
        return Err(Refusal::new(
            "extrude.zero-direction",
            "extrusion direction is zero",
        ));
    }
    let projection = normal.dot(&direction);
    let forward = sign(&projection, dec)?;
    if forward == Sign::Zero {
        return Err(Refusal::new(
            "extrude.zero-thickness",
            "ray lies in the sketch plane",
        ));
    }
    if sign(&normal.cross(&direction).norm_squared(), dec)? != Sign::Zero {
        return Err(coverage(
            "oblique circular sphere limits need a different whole-profile proof",
        ));
    }
    let origin = pose.apply_point(&Point3::new(
        centre.0[0].value().clone().into(),
        centre.0[1].value().clone().into(),
        Real::zero(),
    ));
    let displacement = origin.minus(&sphere.center);
    if sign(&displacement.cross(&normal).norm_squared(), dec)? != Sign::Zero {
        return Err(coverage(
            "offset circular sphere limits need a different whole-profile proof",
        ));
    }
    let normal = if forward == Sign::Negative {
        normal.scaled(&Real::from_i64(-1))
    } else {
        normal
    };
    let speed = normal.dot(&direction);
    let start = displacement.dot(&normal);
    let r = radius.value();
    if r.sign() != Sign::Positive {
        return Err(Refusal::new(
            "profile.non-positive-radius",
            "circle radius must be positive",
        ));
    }
    if sign(&(&sphere.radius - &Real::from(r.clone())), dec)? != Sign::Positive {
        return Err(coverage(
            "circular profile reaches the sphere silhouette or has missing rays",
        ));
    }
    // For every rho² in [0,r²], D = 4 speed² (R²-rho²) is strictly
    // positive. The two roots are (-start +/- sqrt(R²-rho²))/speed.
    // Strict exterior-before-entry or strict interior at the rim proves a
    // single positive first branch on the entire closed disk. Equality is
    // excluded: at the axis or rim it can change the least positive root.
    let branch = if sign(&(&start + &sphere.radius), dec)? == Sign::Negative {
        Branch::Entry
    } else if sign(
        &(&start * &start + Real::from(r * r) - &sphere.radius * &sphere.radius),
        dec,
    )? == Sign::Negative
    {
        Branch::Exit
    } else {
        return Err(coverage(
            "circular sphere start is not uniformly before entry or inside material",
        ));
    };
    // R bounds sqrt(R²-rho²) everywhere, including the pole. This bound is
    // deliberately on the raw line parameter, rather than a normalised ray.
    let upper = projected::enclose(&((&sphere.radius - &start) / speed), dec)?.hi();
    let extent = upper + Rational::one();
    let frame = Frame {
        origin,
        u,
        v,
        direction,
    };
    let bounds = [[-r, r.clone()], [-r, r.clone()]];
    let bound = source
        .radius
        .clone()
        .max(build::frame_radius(&frame, &bounds, dec)?);
    let midpoint = Plane::new(sphere.center.clone(), normal.clone(), frame.u.clone(), dec)?;
    Ok(AuthoredCap {
        definition: CurvedCapped {
            scene: source,
            frame,
            bounds,
            sides: Vec::new(),
            polygon: None,
            circle: Some(Circle {
                radius: r.clone(),
                normal,
            }),
            branch,
            midpoint,
            extent,
        },
        side_ids: vec![id.clone()],
        radius: bound,
    })
}
