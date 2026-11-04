//! Exact rational bands and shared pole topology for equal crossing cylinders.
//! Incomplete pole fans retain a corner refusal: a collapsed two-contact cone
//! is not a finite-area spherical corner.

use catallax_numbers::{Decider, Outcome, Real, Refusal, Sign, ZeroProof};
use catallax_shapes::{Curve, Cylinder, Point3, Surface};
use catallax_topology::{EdgeId, FaceId, Model, Sense};

use crate::{
    fillet::{FilletDefinition, UNSUPPORTED_CORNER},
    profile::sign,
};

fn exact_zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(value) {
        Outcome::ProvenCoincident(proof) => Ok(proof.proof != ZeroProof::SeparationBound),
        Outcome::Refused(error) => Err(error),
        _ => Ok(false),
    }
}

fn axis_foot(cylinder: &Cylinder, point: &Point3) -> Point3 {
    cylinder.origin.plus(&cylinder.axis.scaled(
        &(point.minus(&cylinder.origin).dot(&cylinder.axis) / cylinder.axis.norm_squared()),
    ))
}

fn reaches_collapsed_pole(
    model: &Model,
    edge: EdgeId,
    faces: [FaceId; 2],
    radius: &Real,
    convex: bool,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let [first, second] = faces.map(|face| model.face(face));
    let (Surface::Cylinder(a), Surface::Cylinder(b)) = (
        model.geom().surface(first.surface),
        model.geom().surface(second.surface),
    ) else {
        return Ok(false);
    };
    let edge = model.edge(edge);
    if first.sense != second.sense
        || !matches!(
            model.geom().curve(edge.curve),
            Curve::Ellipse(_) | Curve::Intersection(_)
        )
        || !exact_zero(&(&a.radius_sq - &b.radius_sq), dec)?
    {
        return Ok(false);
    }
    // Equal face senses give equal signed offsets. A zero or negative inward
    // offset has a different failure; no positive-radius pole is claimed there.
    if (first.sense == Sense::Same) == convex
        && sign(&(&a.radius_sq - &radius.square()), dec)? != Sign::Positive
    {
        return Ok(false);
    }
    let normal = a.axis.cross(&b.axis);
    if sign(&normal.norm_squared(), dec)? != Sign::Positive
        || !exact_zero(&b.origin.minus(&a.origin).dot(&normal), dec)?
    {
        return Ok(false);
    }
    if edge.start.is_none() && edge.end.is_none() {
        // A whole equal-cylinder intersection branch includes both poles.
        return Ok(true);
    }
    for vertex in [edge.start, edge.end].into_iter().flatten() {
        let point = model.geom().position(model.vertex(vertex).point);
        if exact_zero(
            &axis_foot(a, &point)
                .minus(&axis_foot(b, &point))
                .norm_squared(),
            dec,
        )? {
            return Ok(true);
        }
    }
    // A trimmed arc that stops before either pole is not this corner case.
    Ok(false)
}

pub(super) fn corner_refusal(
    model: &Model,
    edge: EdgeId,
    faces: [FaceId; 2],
    radius: &Real,
    convex: bool,
    refusal: Refusal,
    dec: &mut Decider,
) -> Refusal {
    if refusal.code == "fillet.unsupported-contact"
        && matches!(
            reaches_collapsed_pole(model, edge, faces, radius, convex, dec),
            Ok(true)
        )
    {
        Refusal::new(
            UNSUPPORTED_CORNER,
            "this equal-radius crossing has no admitted full-fan construction at its collapsed poles; a finite SphereSector patch cannot replace that boundary topology",
        )
    } else {
        refusal
    }
}

// A complete pole-to-pole arc, rather than an arbitrary trim reaching one pole.
fn is_arc(
    model: &Model,
    edge: EdgeId,
    faces: [FaceId; 2],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let [first, second] = faces.map(|face| model.face(face));
    let (Surface::Cylinder(a), Surface::Cylinder(b)) = (
        model.geom().surface(first.surface),
        model.geom().surface(second.surface),
    ) else {
        return Ok(false);
    };
    let edge = model.edge(edge);
    let (Some(start), Some(end)) = (edge.start, edge.end) else {
        return Ok(false);
    };
    if start == end
        || first.sense != second.sense
        || !matches!(model.geom().curve(edge.curve), Curve::Ellipse(_))
        || !exact_zero(&(&a.radius_sq - &b.radius_sq), dec)?
        || sign(&a.axis.cross(&b.axis).norm_squared(), dec)? != Sign::Positive
        || !exact_zero(&b.origin.minus(&a.origin).dot(&a.axis.cross(&b.axis)), dec)?
    {
        return Ok(false);
    }
    for vertex in [start, end] {
        let point = model.geom().position(model.vertex(vertex).point);
        if !exact_zero(
            &axis_foot(a, &point)
                .minus(&axis_foot(b, &point))
                .norm_squared(),
            dec,
        )? {
            return Ok(false);
        }
    }
    Ok(true)
}

fn same_cylinder(a: &Cylinder, b: &Cylinder, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(exact_zero(&(&a.radius_sq - &b.radius_sq), dec)?
        && exact_zero(&a.axis.cross(&b.axis).norm_squared(), dec)?
        && exact_zero(
            &a.origin.minus(&b.origin).cross(&a.axis).norm_squared(),
            dec,
        )?)
}

fn witness(
    model: &mut Model,
    edge: EdgeId,
    dec: &mut Decider,
) -> Result<(Point3, catallax_shapes::Vector3), Refusal> {
    let Curve::Ellipse(ellipse) = model.geom().curve(model.edge(edge).curve).clone() else {
        return Err(super::bad("a pole-to-pole strip needs its source ellipse"));
    };
    let major = ellipse
        .major_dir
        .scaled(&(&ellipse.major_sq / &ellipse.major_dir.norm_squared()).sqrt());
    for positive in [true, false] {
        let at = if positive {
            ellipse.center.plus(&major)
        } else {
            ellipse.center.plus(&major.scaled(&Real::from_i64(-1)))
        };
        let point = model
            .geom_mut()
            .add_point(catallax_shapes::PointDef::Explicit(at.clone()));
        if super::on_trimmed(model, edge, point, dec)? {
            let tangent = ellipse.normal.cross(&at.minus(&ellipse.center));
            return Ok((at, tangent));
        }
    }
    Err(super::bad(
        "neither source-ellipse half belongs to the selected arc",
    ))
}

/// Builds complete collapsed pole fans before the ordinary strip-end machinery.
/// Partial fans remain a corner refusal; a zero-area sphere is never inserted.
pub(super) fn build(
    model: &mut Model,
    body: catallax_topology::BodyId,
    selected: &[EdgeId],
    uses: &super::Uses,
    radius: &Real,
    step: &catallax_recipe::StepId,
    dec: &mut Decider,
) -> Result<Option<(catallax_topology::BodyId, Vec<FilletDefinition>)>, Refusal> {
    use catallax_numbers::OperationSpan;
    use catallax_recipe::{Origin, Role};
    use catallax_shapes::SplineSupport;
    use std::collections::BTreeSet;

    let mut arcs = BTreeSet::new();
    for edge in selected {
        let [(a, _), (b, _)] = uses[edge].as_slice() else {
            continue;
        };
        if is_arc(model, *edge, [*a, *b], dec)? {
            arcs.insert(*edge);
        }
    }
    if arcs.is_empty() {
        return Ok(None);
    }
    if selected.iter().any(|edge| !arcs.contains(edge)) {
        return Err(Refusal::new(
            UNSUPPORTED_CORNER,
            "collapsed crossing fans cannot yet share a step with other strip-end constructions",
        ));
    }
    let poles: BTreeSet<_> = arcs
        .iter()
        .flat_map(|edge| {
            let edge = model.edge(*edge);
            [edge.start, edge.end].into_iter().flatten()
        })
        .collect();
    for pole in &poles {
        let incident: Vec<_> = uses
            .keys()
            .filter(|edge| {
                let edge = model.edge(**edge);
                edge.start == Some(*pole) || edge.end == Some(*pole)
            })
            .collect();
        if incident.len() != 4 || incident.iter().any(|edge| !arcs.contains(edge)) {
            return Err(Refusal::new(
                UNSUPPORTED_CORNER,
                "a collapsed crossing pole needs all four incident arcs selected",
            ));
        }
    }
    let mut runs = Vec::with_capacity(arcs.len());
    for edge in arcs {
        let [(a, sa), (b, sb)] = uses[&edge].as_slice() else {
            return Err(super::bad("a crossing arc needs two incident faces"));
        };
        let faces = [*a, *b];
        let (at, tangent) = witness(model, edge, dec)?;
        // The other half of the same ellipse has the opposite dihedral sign;
        // classify at this trimmed arc, not at the carrier's positive axis.
        let na = super::outward(model, *a, &at, dec)?;
        let nb = super::outward(model, *b, &at, dec)?;
        let orientation = if *sa == model.edge(edge).sense { 1 } else { -1 };
        let convex = match sign(
            &(nb.dot(&na.cross(&tangent)) * Real::from_i64(orientation)),
            dec,
        )? {
            Sign::Negative => true,
            Sign::Positive => false,
            Sign::Zero => return Err(super::bad("a crossing arc's interior is tangent")),
        };
        let name = model
            .edge_name(body, edge)
            .ok_or_else(|| super::bad("a crossing arc needs an unambiguous origin-pair name"))?;
        let definition = FilletDefinition::ring(
            Origin {
                step: step.clone(),
                role: Role::Blend { edge: name },
            },
            faces.map(|face| model.geom().surface(model.face(face).surface).clone()),
            faces.map(|face| model.face(face).sense == Sense::Same),
            model.geom().curve(model.edge(edge).curve).clone(),
            &at,
            radius.clone(),
            convex,
            dec,
        )?;
        let Surface::BSpline(SplineSupport::RollingSaddle(saddle)) = &definition.surface else {
            return Err(super::bad("a crossing arc has no retained rational band"));
        };
        let [_, wide, narrow] = saddle
            .definition()
            .section(&catallax_numbers::Rational::one())?;
        for (face, point) in faces.iter().zip([wide, narrow]) {
            if !catallax_boolean::face_contains(model, *face, &point, dec)? {
                return Err(Refusal::new(
                    "fillet.unsupported-overflow",
                    "a crossing contact leaves its input face",
                ));
            }
        }
        runs.push((edge, faces, [*sa, *sb], definition));
    }
    // All non-source supports must be clear of each complete positive-weight
    // control hull. Together with one contact in each face, this proves that
    // neither the band nor either rail escapes its input trim.
    for (_, _, _, definition) in &runs {
        let Surface::BSpline(SplineSupport::RollingSaddle(saddle)) = &definition.surface else {
            unreachable!("crossing runs retain their rational band")
        };
        for face in model.faces_of_body(body) {
            let surface = model.geom().surface(model.face(face).surface);
            let mut source = false;
            if let Surface::Cylinder(other) = &surface {
                for support in &definition.supports {
                    if let Surface::Cylinder(cylinder) = support {
                        source |= same_cylinder(cylinder, other, dec)?;
                    }
                }
            }
            if source {
                continue;
            }
            let clear = if let Surface::Plane(plane) = &surface {
                saddle.clear_of(plane, dec)?
            } else {
                false
            };
            if !clear {
                return Err(Refusal::new(
                    "fillet.unsupported-overflow",
                    "a crossing band is not proved clear of another input support",
                ));
            }
        }
    }
    let trimming = OperationSpan::begin("brep", step, "fillet.trimming", dec);
    let mut edit = crate::edge_blend::Replacement::default();
    let mut definitions = Vec::with_capacity(runs.len());
    for (old, faces, senses, definition) in runs {
        let input = model.edge(old).clone();
        let mut boundary = Vec::with_capacity(2);
        for i in 0..2 {
            let Curve::BSpline(rail) = &definition.contacts[i] else {
                return Err(super::bad("a crossing band lost a rational rail"));
            };
            let start = input
                .start
                .ok_or_else(|| super::bad("a crossing arc lost its pole"))?;
            let point = model.geom().position(model.vertex(start).point);
            let sense = if exact_zero(&point.minus(rail.control(0)).norm_squared(), dec)? {
                Sense::Same
            } else {
                let end = rail.control(rail.control_count() - 1);
                if !exact_zero(&point.minus(end).norm_squared(), dec)? {
                    return Err(super::bad(
                        "a crossing rail does not end on its shared pole",
                    ));
                }
                Sense::Reversed
            };
            let curve = model.geom_mut().add_curve(definition.contacts[i].clone());
            // These exact contact laws lie on their source cylinders over the
            // whole chart; no fitted-rail incidence is being asserted.
            let support = model.face(faces[i]).surface;
            model.geom_mut().record_fact(catallax_numbers::Fact {
                kind: catallax_numbers::FactKind::Incident,
                a: catallax_shapes::GeomKey::Curve(curve),
                b: catallax_shapes::GeomKey::Surface(support),
                at: None,
                origin: catallax_numbers::FactOrigin::Construction("fillet.crossing-contact"),
            });
            let edge = model.add_edge(catallax_topology::Edge {
                curve,
                sense,
                start: input.start,
                end: input.end,
            });
            edit.edges.insert((faces[i], old), edge);
            boundary.push((edge, super::reverse(senses[i])));
        }
        let face = crate::extrude::face(
            model,
            definition.surface.clone(),
            if definition.convex {
                Sense::Reversed
            } else {
                Sense::Same
            },
            vec![boundary],
            None,
            vec![definition.origin.clone()],
        );
        edit.faces.push((faces[0], face));
        definitions.push(definition);
    }
    let result = crate::edge_blend::apply(model, body, &edit)?;
    model.record_blend_base(result, body);
    trimming.finish(dec, None, None);
    Ok(Some((result, definitions)))
}
