//! Face-bounded sweeps on the evaluated B-rep.

use super::*;
use catallax_boolean::{PatchSelection, arrange};
use catallax_numbers::MarkLog;

use catallax_shapes::{BezierExtrusion, Curve, Isometry, Line, SplineSupport, Vector3};
use catallax_topology::{
    Body, BodyKind, Edge, EdgeId, FaceId, Loop, Lump, Sense, Shell, VertexId, transformed_body,
};

pub(crate) fn may_retry(error: &Refusal) -> bool {
    matches!(
        error.code,
        "features.unsupported-limit-surface" | "features.unsupported-curved-limit"
    )
}

fn unsupported(detail: &str) -> Refusal {
    bad("features.unsupported-curved-limit", detail)
}

fn select(model: &Model, body: BodyId, reference: &FaceRef) -> Result<FaceId, Refusal> {
    let faces = model.faces_of_body(body);
    let uses: BTreeMap<_, _> = faces
        .iter()
        .map(|&face| {
            (
                face,
                paths(model, face)
                    .into_iter()
                    .flatten()
                    .map(|(edge, _)| edge)
                    .collect::<BTreeSet<_>>(),
            )
        })
        .collect();
    let matches: Vec<_> = faces
        .iter()
        .copied()
        .filter(|face| {
            reference.matches_origins(model.face_origins(*face))
                && reference.touching.iter().all(|pattern| {
                    faces.iter().any(|other| {
                        other != face
                            && !uses[face].is_disjoint(&uses[other])
                            && model
                                .face_origins(*other)
                                .iter()
                                .any(|origin| pattern.matches(origin))
                    })
                })
        })
        .collect();
    match matches.as_slice() {
        [face] => Ok(*face),
        [] => Err(bad(
            "recipe.missing-reference",
            "limit matches no trimmed target face",
        )),
        _ => Err(unsupported(
            "face placement needs one resolved trimmed face",
        )),
    }
}

fn direction(
    profile: &Profile,
    authored: Option<&catallax_recipe::Vector3>,
) -> Result<RVector3, Refusal> {
    let normal = normal(profile);
    let direction = authored.map(vector).unwrap_or_else(|| normal.clone());
    if direction.is_zero() {
        return Err(bad("extrude.zero-direction", "extrusion direction is zero"));
    }
    if direction.dot(&normal).is_zero() {
        return Err(bad(
            "extrude.zero-thickness",
            "extrusion direction lies in the sketch plane",
        ));
    }
    Ok(direction)
}

fn normal(profile: &Profile) -> RVector3 {
    crate::profile::rational_vector(
        &profile.pose,
        &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
    )
}

fn extent(profile: &Profile, target_radius: &Rational, direction: &RVector3) -> Rational {
    (&crate::evaluate::profile_radius(profile) + target_radius + Rational::one())
        / direction.dot(&normal(profile)).abs()
}

fn shifted(profile: &Profile, offset: &RVector3) -> Profile {
    let mut profile = profile.clone();
    profile.pose.translation = profile.pose.translation.plus(&offset.to_real());
    profile
}

fn move_recipe(recipe: &Recipe, profile: &StepId, offset: &RVector3) -> Recipe {
    let mut recipe = recipe.clone();
    let step = recipe
        .steps
        .iter_mut()
        .find(|step| &step.id == profile)
        .expect("checked sketch");
    let Operation::Sketch { plane, .. } = &mut step.op else {
        unreachable!("checked sketch")
    };
    for (component, delta) in plane
        .translation
        .0
        .iter_mut()
        .zip([&offset.x, &offset.y, &offset.z])
    {
        *component = (component.value() + delta).into();
    }
    recipe
}

fn origin(step: &StepId, role: Role) -> Origin {
    Origin {
        step: step.clone(),
        role,
    }
}

fn named(model: &Model, face: FaceId, role: Role) -> bool {
    model
        .face_origins(face)
        .iter()
        .any(|origin| origin.role == role)
}

fn paths(model: &Model, face: FaceId) -> Vec<Vec<(EdgeId, Sense)>> {
    model
        .face(face)
        .loops
        .iter()
        .filter_map(|id| {
            let Loop::Coedges(uses) = model.loop_(*id) else {
                return None;
            };
            Some(
                uses.iter()
                    .map(|coedge| {
                        let coedge = model.coedge(*coedge);
                        (coedge.edge, coedge.sense)
                    })
                    .collect(),
            )
        })
        .collect()
}

fn sketch_chart(profile: &Profile, direction: &RVector3, point: &Point3) -> (P, Real) {
    let inverse = profile.pose.inverse();
    let local = inverse.apply_point(point);
    let ray = inverse.apply_vector(&direction.to_real());
    let time = &local.z / &ray.z;
    (
        [&local.x - &(&time * &ray.x), &local.y - &(&time * &ray.y)],
        time,
    )
}

fn planar(
    model: &mut Model,
    selected: FaceId,
    plane: &Plane,
    profile: &Profile,
    direction: &RVector3,
    extent: &Rational,
    step: &StepId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let orientation = match sign(&plane.normal.dot(&direction.to_real()), dec)? {
        Sign::Positive => Sense::Same,
        Sign::Negative => Sense::Reversed,
        Sign::Zero => {
            return Err(bad(
                "features.limit-parallel",
                "direction lies in the limiting plane",
            ));
        }
    };
    clip(
        model,
        selected,
        profile,
        direction,
        extent,
        step,
        |point, dec| {
            let side = sign(&point.minus(&plane.origin).dot(&plane.normal), dec)?;
            Ok(match (side, orientation) {
                (Sign::Positive, Sense::Same) | (Sign::Negative, Sense::Reversed) => Sign::Positive,
                (Sign::Zero, _) => Sign::Zero,
                _ => Sign::Negative,
            })
        },
        orientation,
        dec,
    )
}

#[allow(clippy::too_many_arguments)]
fn bezier(
    model: &mut Model,
    selected: FaceId,
    surface: &BezierExtrusion,
    profile: &Profile,
    direction: &RVector3,
    extent: &Rational,
    step: &StepId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let controls = surface.controls();
    let base = &controls[0];
    let chord = controls.last().unwrap().minus(base);
    let axis = surface.direction();
    let across = chord.cross(axis);
    let direction_real = direction.to_real();
    let skew = across.cross(&direction_real);
    for component in [&skew.x, &skew.y, &skew.z] {
        if sign(component, dec)? != Sign::Zero {
            return Err(unsupported(
                "the Bezier limit needs a proved graph along the sweep direction",
            ));
        }
    }
    // The existing Bezier carrier proves monotonicity in the chord chart.
    // Its two chart coordinates are constant on each accepted sweep line.
    for segment in profile.loops.iter().flat_map(|boundary| &boundary.segments) {
        if segment.circle.is_some() {
            return Err(unsupported(
                "the Bezier limit needs a certified circular footprint",
            ));
        }
        let p = profile.world(&segment.end, &RVector3::zero()).to_real();
        let from_base = p.minus(base);
        for vector in [&chord, axis] {
            let station = from_base.dot(vector);
            if sign(&station, dec)? == Sign::Negative
                || sign(&(&station - &vector.norm_squared()), dec)? == Sign::Positive
            {
                return Err(bad(
                    "features.limit-missed",
                    "the profile extends beyond the finite Bezier chart",
                ));
            }
        }
    }
    let orientation = if sign(&across.dot(&direction_real), dec)? == Sign::Positive {
        Sense::Same
    } else {
        Sense::Reversed
    };
    clip(
        model,
        selected,
        profile,
        direction,
        extent,
        step,
        |point, dec| {
            let side = surface.side_at(point, dec)?;
            Ok(match (side, orientation) {
                (Sign::Positive, Sense::Same) | (Sign::Negative, Sense::Reversed) => Sign::Positive,
                (Sign::Zero, _) => Sign::Zero,
                _ => Sign::Negative,
            })
        },
        orientation,
        dec,
    )
}

#[allow(clippy::too_many_arguments)]
fn clip(
    model: &mut Model,
    selected: FaceId,
    profile: &Profile,
    direction: &RVector3,
    extent: &Rational,
    step: &StepId,
    mut side: impl FnMut(&Point3, &mut Decider) -> Result<Sign, Refusal>,
    orientation: Sense,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let tool = crate::extrude::extrude_vector(
        model,
        profile,
        step,
        &direction.scaled(extent),
        [Some(Role::Start), None],
        None,
        dec,
    )?;
    let target = model.face(selected).clone();
    let flip = if target.sense == orientation {
        Sense::Same
    } else {
        Sense::Reversed
    };
    let cap = model.add_named_face(target, vec![origin(step, Role::End)]);
    let answer = arrange(
        model,
        tool,
        &[cap],
        |model, owner, point, dec| {
            let keep = if owner == cap {
                let (uv, time) = sketch_chart(profile, direction, point);
                if sign(&time, dec)? == Sign::Positive
                    && profile_location(profile, &uv, dec)? == RegionLocation::Inside
                {
                    PatchSelection::Keep(flip)
                } else {
                    PatchSelection::Discard
                }
            } else {
                match side(point, dec)? {
                    Sign::Negative => PatchSelection::Keep(Sense::Same),
                    Sign::Positive => PatchSelection::Discard,
                    Sign::Zero => return Ok(None),
                }
            };
            // A partial start would create a curtain through the profile or a zero-depth ridge.
            if owner != cap
                && named(model, owner, Role::Start)
                && matches!(keep, PatchSelection::Discard)
            {
                return Err(bad(
                    "features.limit-missed",
                    "a profile line meets no positive limiting cap",
                ));
            }
            Ok(Some(keep))
        },
        dec,
    )
    .map_err(|error| error.refusal)?;
    if !answer.marks.is_empty() {
        return Err(unsupported("the trimmed limit did not yield an exact cap"));
    }
    let faces = model.faces_of_body(answer.body);
    if !faces.iter().any(|face| named(model, *face, Role::End))
        || faces
            .iter()
            .any(|face| model.face_origins(*face).is_empty())
    {
        return Err(bad(
            "features.limit-missed",
            "the trimmed face does not cover the entire profile",
        ));
    }
    Ok(answer.body)
}

// A normal extrusion of a single circle is exactly the stock-cylinder scene.
fn cylinder_recipe(
    recipe: &Recipe,
    target: &StepId,
    dec: &mut Decider,
) -> Result<Option<Recipe>, Refusal> {
    let Some(Operation::Extrude {
        profile,
        distance,
        direction,
    }) = recipe.step(target).map(|step| &step.op)
    else {
        return Ok(None);
    };
    let Some(Operation::Sketch { plane, loops }) = recipe.step(profile).map(|step| &step.op) else {
        return Ok(None);
    };
    let [boundary] = loops.as_slice() else {
        return Ok(None);
    };
    let [catallax_recipe::Segment::Circle { centre, radius, .. }] = boundary.0.as_slice() else {
        return Ok(None);
    };
    let frame = crate::profile::placement(plane, dec)?;
    let axis = crate::profile::rational_vector(
        &frame,
        &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
    );
    let direction = direction
        .as_ref()
        .map(vector)
        .unwrap_or_else(|| axis.clone());
    if !direction.cross(&axis).is_zero() {
        return Ok(None);
    }
    let backwards = direction.dot(&axis).sign() != distance.value().sign();
    let local = RPoint3::new(
        centre.0[0].value().clone(),
        centre.0[1].value().clone(),
        if backwards {
            -distance.value().abs()
        } else {
            Rational::zero()
        },
    );
    let RPoint3 { x, y, z } = crate::profile::rational_point(&frame, &local);
    let mut pose = plane.clone();
    pose.translation = catallax_recipe::Vector3([x.into(), y.into(), z.into()]);
    let mut recipe = recipe.clone();
    let step = recipe
        .steps
        .iter_mut()
        .find(|step| &step.id == target)
        .expect("checked target");
    step.op = Operation::Cylinder {
        radius: radius.clone(),
        height: distance.value().abs().into(),
        pose,
    };
    Ok(Some(recipe))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn limit(
    model: &mut Model,
    recipe: &Recipe,
    target: (&StepId, BodyId, &Rational),
    sketch: (&StepId, &Profile),
    authored: Option<&catallax_recipe::Vector3>,
    reference: Option<&FaceRef>,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(BodyId, Rational), Refusal> {
    let (target_id, target_body, bound) = target;
    let (profile_id, profile) = sketch;
    let direction = direction(profile, authored)?;
    let reach = extent(profile, bound, &direction);
    let face = reference
        .map(|reference| select(model, target_body, reference))
        .transpose()?;
    if let Some(face) = face {
        let surface = model.geom().surface(model.face(face).surface).clone();
        let body = match &surface {
            Surface::Plane(plane) => Some(planar(
                model, face, plane, profile, &direction, &reach, step, dec,
            )?),
            Surface::BSpline(SplineSupport::BezierExtrusion(support)) => Some(bezier(
                model, face, support, profile, &direction, &reach, step, dec,
            )?),
            Surface::Cylinder(_) => None,
            _ => {
                return Err(unsupported(
                    "the evaluated limiting face has no supported first-hit proof",
                ));
            }
        };
        if let Some(body) = body {
            return Ok((
                body,
                bound.clone().max(crate::evaluate::profile_radius(profile)),
            ));
        }
    }
    let cylinder = cylinder_recipe(recipe, target_id, dec)?.ok_or_else(|| {
        unsupported("the evaluated limiting target has no supported first-hit proof")
    })?;
    // The target's original role was resolved on the real body, and only its wall reaches here.
    // The equivalent cylinder has no dependency on the circle's authored segment id.
    let mut cylinder = cylinder;
    cylinder.result = target_id.clone();
    let last = cylinder
        .steps
        .iter()
        .position(|step| &step.id == profile_id)
        .expect("checked profile")
        .max(
            cylinder
                .steps
                .iter()
                .position(|step| &step.id == target_id)
                .expect("checked target"),
        );
    cylinder.steps.truncate(last + 1);
    let cap = curved::AuthoredCap::compile(&cylinder, target_id, profile_id, authored, None, dec)?;
    Ok((cap.build(model, step, dec)?, cap.radius().clone()))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn start(
    model: &mut Model,
    recipe: &Recipe,
    target: (&StepId, BodyId, &Rational),
    sketch: (&StepId, &Profile),
    authored: Option<&catallax_recipe::Vector3>,
    reference: &FaceRef,
    distance: &Rational,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(BodyId, Rational), Refusal> {
    let (target_id, target_body, bound) = target;
    let (profile_id, profile) = sketch;
    let direction = direction(profile, authored)?;
    let face = select(model, target_body, reference)?;
    if profile.is_empty() {
        return Ok((catallax_topology::empty(model), Rational::zero()));
    }
    if distance.is_zero() {
        return Err(bad(
            "extrude.zero-distance",
            "extrusion distance must not be zero",
        ));
    }
    let reach = extent(profile, bound, &direction);
    let offset = direction.scaled(&(-reach));
    let behind = move_recipe(recipe, profile_id, &offset);
    let shifted_profile = shifted(profile, &offset);
    let column = match behind.step(target_id).map(|step| &step.op) {
        Some(Operation::Torus { .. }) => {
            let cap = toroidal::AnnularCap::compile(
                &behind,
                target_id,
                profile_id,
                authored,
                Some(reference),
                dec,
            )?;
            cap.build(model, target_body, &shifted_profile, step, dec)?
        }
        Some(Operation::Sphere { .. } | Operation::Cylinder { .. } | Operation::Cone { .. })
            if !matches!(
                model.geom().surface(model.face(face).surface),
                Surface::Plane(_)
            ) =>
        {
            let cap = curved::AuthoredCap::compile(
                &behind,
                target_id,
                profile_id,
                authored,
                Some(reference),
                dec,
            )?;
            cap.build(model, step, dec)?
        }
        _ => {
            limit(
                model,
                &behind,
                target,
                (profile_id, &shifted_profile),
                authored,
                Some(reference),
                step,
                dec,
            )?
            .0
        }
    };
    let norm = direction.norm_squared();
    let length = crate::extrude::sqrt(&norm)
        .map(Real::from)
        .unwrap_or_else(|| Real::from(norm).sqrt());
    let displacement = direction
        .to_real()
        .scaled(&(&Real::from(distance.clone()) / &length));
    let body = translated_cap_sweep(model, column, displacement, distance.sign(), step, dec)?;
    Ok((body, bound + &distance.abs()))
}

fn connector(
    model: &mut Model,
    connectors: &mut BTreeMap<VertexId, EdgeId>,
    start: VertexId,
    end: VertexId,
    vector: &Vector3,
) -> EdgeId {
    *connectors.entry(start).or_insert_with(|| {
        let origin = model.geom().position(model.vertex(start).point);
        let curve = model.geom_mut().add_curve(Curve::Line(Line {
            origin,
            dir: vector.clone(),
        }));
        model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(start),
            end: Some(end),
            tolerance: None,
        })
    })
}

fn translated_cap_sweep(
    model: &mut Model,
    column: BodyId,
    vector: Vector3,
    distance: Sign,
    step: &StepId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    dec.check_deadline()?;
    let mut map = Isometry::identity();
    map.translation = vector.clone();
    let copy = transformed_body(model, column, &map, Clone::clone);
    let sources = model.faces_of_body(column);
    let copies = model.faces_of_body(copy);
    let mut sides = BTreeMap::new();
    for &face in &sources {
        if model
            .face_origins(face)
            .iter()
            .any(|origin| matches!(origin.role, Role::Side { .. }))
        {
            for (edge, _) in paths(model, face).into_iter().flatten() {
                sides.insert(edge, face);
            }
        }
    }
    let mut faces = Vec::new();
    let mut connectors = BTreeMap::new();
    let reverse = crate::extrude::reverse;
    for (old, moved) in sources.into_iter().zip(copies) {
        if !named(model, old, Role::End) {
            continue;
        }
        let old_paths = paths(model, old);
        let moved_paths = paths(model, moved);
        for (a, b) in old_paths.iter().flatten().zip(moved_paths.iter().flatten()) {
            let Some(&side) = sides.get(&a.0) else {
                continue;
            };
            let edge = model.edge(a.0).clone();
            let moved_edge = model.edge(b.0).clone();
            let use_old = if distance == Sign::Positive {
                a.1
            } else {
                reverse(a.1)
            };
            let use_moved = if distance == Sign::Positive {
                reverse(b.1)
            } else {
                b.1
            };
            let walls = if let (Some(start), Some(end), Some(new_start), Some(new_end)) =
                (edge.start, edge.end, moved_edge.start, moved_edge.end)
            {
                let start_wall = connector(model, &mut connectors, start, new_start, &vector);
                let end_wall = connector(model, &mut connectors, end, new_end, &vector);
                let (first, second) = if use_old == Sense::Same {
                    (start_wall, end_wall)
                } else {
                    (end_wall, start_wall)
                };
                vec![vec![
                    (a.0, use_old),
                    (second, Sense::Same),
                    (b.0, use_moved),
                    (first, Sense::Reversed),
                ]]
            } else {
                vec![vec![(a.0, use_old)], vec![(b.0, use_moved)]]
            };
            let surface = model.face(side).surface;
            let sense = model.face(side).sense;
            let origins = model.face_origins(side).to_vec();
            faces.push(crate::extrude::face_on(
                model, surface, sense, walls, None, origins,
            ));
        }
        for (face, mut paths, role, flip) in [
            (old, old_paths, Role::Start, distance == Sign::Positive),
            (moved, moved_paths, Role::End, distance == Sign::Negative),
        ] {
            let support = model.face(face).surface;
            let mut sense = model.face(face).sense;
            if flip {
                sense = reverse(sense);
                for path in &mut paths {
                    path.reverse();
                    for (_, sense) in path {
                        *sense = reverse(*sense);
                    }
                }
            }
            faces.push(crate::extrude::face_on(
                model,
                support,
                sense,
                paths,
                None,
                vec![origin(step, role)],
            ));
        }
    }
    if faces.is_empty() {
        return Err(bad(
            "features.limit-missed",
            "the selected start face does not cover the profile",
        ));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        wires: Vec::new(),
        marks: MarkLog::new(),
    });
    if !model.check(body).is_empty() {
        return Err(bad(
            "features.invalid-topology",
            "the translated cap and its sides do not form a closed solid",
        ));
    }
    Ok(body)
}
