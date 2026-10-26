//! Exact three-face circular cap and independent finite analytic cells.

use super::*;
use catallax_recipe::{Origin, Role};
use catallax_shapes::{Circle, Curve, Cylinder, Sphere};
use catallax_topology::{Body, BodyId, BodyKind, Edge, Lump, Model, Sense, Shell};

/// The certificate makes each rim a whole circle and every retained cap point
/// belong to the same hemisphere. No intersection sampling or polygonisation.
pub(super) fn build(
    cap: &CurvedCapped,
    model: &mut Model,
    step: &StepId,
    sides: &[ElementId],
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    dec.check_deadline()?;
    if sides.len() != 1 {
        return Err(coverage("a circular cap requires its one authored side id"));
    }
    let circle = cap.circle.as_ref().expect("circle dispatch");
    let frame = &cap.frame;
    let normal = &circle.normal;
    let Surface::Sphere(sphere) = &cap.scene.surface else {
        unreachable!("sphere certificate")
    };
    let rim_height = cap.height(&[circle.radius.clone().into(), Real::zero()], dec)?;
    let rim_center = frame.origin.plus(&frame.direction.scaled(&rim_height));
    let mut rings = Vec::new();
    for center in [frame.origin.clone(), rim_center] {
        let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
            center,
            normal: normal.clone(),
            radius_sq: (&circle.radius * &circle.radius).into(),
            ref_dir: frame.u.clone(),
        }));
        rings.push(model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: None,
            end: None,
        }));
    }
    let origin = |role| {
        vec![Origin {
            step: step.clone(),
            role,
        }]
    };
    let start = crate::extrude::face(
        model,
        Surface::Plane(Plane::new(
            frame.origin.clone(),
            normal.scaled(&Real::from_i64(-1)),
            frame.u.clone(),
            dec,
        )?),
        Sense::Same,
        vec![vec![(rings[0], Sense::Reversed)]],
        None,
        origin(Role::Start),
    );
    let wall = crate::extrude::face(
        model,
        Surface::Cylinder(Cylinder::new(
            frame.origin.clone(),
            normal.clone(),
            circle.radius.clone().into(),
            frame.u.clone(),
            dec,
        )?),
        Sense::Same,
        vec![
            vec![(rings[0], Sense::Same)],
            vec![(rings[1], Sense::Reversed)],
        ],
        None,
        origin(Role::Side {
            segment: sides[0].clone(),
            path: None,
        }),
    );
    // Coaxial rings are exact constant-latitude trims in this chart. The pole
    // is an ordinary interior point of the finite cap, not a topology vertex.
    let seam = frame.u.clone();
    let pole = normal.clone();
    let end = crate::extrude::face(
        model,
        Surface::Sphere(Sphere::new(
            sphere.center.clone(),
            sphere.radius.clone(),
            pole,
            seam,
            dec,
        )?),
        if cap.branch == Branch::Entry {
            Sense::Reversed
        } else {
            Sense::Same
        },
        vec![vec![(rings[1], Sense::Same)]],
        None,
        origin(Role::End),
    );
    let shell = model.add_shell(Shell {
        faces: vec![start, wall, end],
    });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let result = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    });
    if !model.check(result).is_empty() {
        return Err(Refusal::new(
            "features.invalid-topology",
            "circular sphere cap failed topology validation",
        ));
    }
    Ok(result)
}

/// Analytic disk column and sphere sets compiled by the field without topology.
pub(super) fn construction(
    cap: &CurvedCapped,
    sides: &[ElementId],
    dec: &mut Decider,
) -> Result<CapConstruction, Refusal> {
    dec.check_deadline()?;
    if sides.len() != 1 {
        return Err(coverage("a circular cap requires its one authored side id"));
    }
    let circle = cap.circle.as_ref().expect("circle dispatch");
    let f = &cap.frame;
    let n = &circle.normal;
    let offset = f.direction.scaled(&Real::from(cap.extent.clone()));
    let start = Plane::new(
        f.origin.clone(),
        n.scaled(&Real::from_i64(-1)),
        f.u.clone(),
        dec,
    )?;
    let end = Plane::new(f.origin.plus(&offset), n.clone(), f.u.clone(), dec)?;
    let cylinder = Surface::Cylinder(Cylinder::new(
        f.origin.clone(),
        n.clone(),
        circle.radius.clone().into(),
        f.u.clone(),
        dec,
    )?);
    // A column point is a start-disk point plus a fraction of `offset`.
    let bound = build::frame_radius(f, &cap.bounds, dec)?
        + projected::enclose(&offset.norm_squared().sqrt(), dec)?.hi();
    let column = ConvexCell::circular(
        vec![cylinder, Surface::Plane(start), Surface::Plane(end)],
        vec![
            Role::Side {
                segment: sides[0].clone(),
                path: None,
            },
            Role::Start,
            Role::End,
        ],
        bound.clone(),
    );
    let target = ConvexCell::circular(
        vec![cap.scene.surface.clone()],
        vec![Role::End],
        cap.scene.radius.clone(),
    );
    let clipper = if cap.branch == Branch::Entry {
        // A finite box around the column, with only its upper plane at the
        // sphere midpoint. Other supports are strictly outside the result.
        let margin = Real::from(&circle.radius + &Rational::one());
        let mut surfaces = Vec::new();
        for direction in [&f.u, &f.v] {
            for sign in [-1, 1] {
                let normal = direction.scaled(&Real::from_i64(sign));
                surfaces.push(Surface::Plane(Plane::new(
                    f.origin.plus(&normal.scaled(&margin)),
                    normal,
                    n.clone(),
                    dec,
                )?));
            }
        }
        surfaces.push(Surface::Plane(Plane::new(
            f.origin.plus(&n.scaled(&Real::from_i64(-1))),
            n.scaled(&Real::from_i64(-1)),
            f.u.clone(),
            dec,
        )?));
        surfaces.push(Surface::Plane(cap.midpoint.clone()));
        // The expanded start square bounds every transverse section. Axial
        // displacement is at most |centre-start| or one unit below start.
        let wide = &circle.radius + &Rational::one();
        let bounds = [[-&wide, wide.clone()], [-&wide, wide.clone()]];
        let axial = cap.midpoint.origin.minus(&f.origin);
        let bound = build::frame_radius(f, &bounds, dec)?
            + projected::enclose(&axial.norm_squared().sqrt(), dec)?.hi()
            + Rational::one();
        Some(ConvexCell::circular(surfaces, vec![Role::End; 6], bound))
    } else {
        None
    };
    Ok(CapConstruction {
        column,
        target,
        clipper,
    })
}
