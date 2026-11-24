//! Face-bounded sweeps on the evaluated B-rep.

use super::*;
use catallax_boolean::{PatchSelection, arrange};
use catallax_measure::distance::PreparedFace;
use catallax_numbers::MarkLog;
use catallax_recipe::Step;

use catallax_shapes::{Curve, Isometry, Line, SplineSupport, Vector3};

mod bezier;
use catallax_topology::{
    Body, BodyKind, Edge, EdgeId, FaceId, Loop, Lump, Sense, Shell, VertexId, transformed_body,
};

pub(crate) fn may_retry(error: &Refusal) -> bool {
    matches!(
        error.code,
        "features.unsupported-limit-surface"
            | "features.unsupported-curved-limit"
            | "features.bezier-profile"
    )
}

fn unsupported(detail: &str) -> Refusal {
    bad("features.unsupported-curved-limit", detail)
}

pub(crate) fn select(
    model: &mut Model,
    body: BodyId,
    reference: &FaceRef,
    dec: &mut Decider,
) -> Result<FaceId, Refusal> {
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
        _ => match &reference.near {
            Some(near) => nearest(model, &matches, &point(near), dec),
            None => Err(unsupported(
                "face placement needs one resolved trimmed face",
            )),
        },
    }
}

fn nearest(
    model: &mut Model,
    faces: &[FaceId],
    point: &Point3,
    dec: &mut Decider,
) -> Result<FaceId, Refusal> {
    let options = Default::default();
    let mut distances = Vec::with_capacity(faces.len());
    for &face in faces {
        let prepared = PreparedFace::new(model, face, dec)?;
        distances.push(prepared.point_to_face(model, point, &options, dec)?);
    }
    let (closest, upper) = distances
        .iter()
        .enumerate()
        .filter_map(|(i, distance)| Some((i, distance.upper.as_ref()?)))
        .min_by(|a, b| a.1.cmp(b.1))
        .ok_or_else(|| unsupported("no trimmed face has a certified near witness"))?;
    if distances
        .iter()
        .enumerate()
        .all(|(i, distance)| i == closest || upper < &distance.lower)
    {
        return Ok(faces[closest]);
    }
    if distances.iter().all(|d| upper <= &d.lower)
        && distances
            .iter()
            .filter(|d| d.upper.as_ref() == Some(upper))
            .count()
            > 1
    {
        return Err(bad(
            "recipe.ambiguous-reference",
            "near point has equal minimum distance to several trimmed faces",
        ));
    }
    Err(unsupported(
        "trimmed-face distances do not prove a unique nearest face",
    ))
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

#[allow(clippy::too_many_arguments)]
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
    let denominator = plane.normal.dot(&direction.to_real());
    let orientation = match sign(&denominator, dec)? {
        Sign::Positive => Sense::Same,
        Sign::Negative => Sense::Reversed,
        Sign::Zero => {
            return Err(bad(
                "features.limit-parallel",
                "direction lies in the limiting plane",
            ));
        }
    };
    let tool = crate::extrude::extrude_vector(
        model,
        profile,
        step,
        &direction.scaled(extent),
        [Some(Role::Start), None],
        None,
        dec,
    )?;
    for face in model.faces_of_body(tool) {
        if named(model, face, Role::Start) {
            let point = catallax_boolean::face_point(model, face, dec)?.value;
            let time = &plane.origin.minus(&point).dot(&plane.normal) / &denominator;
            if sign(&time, dec)? != Sign::Positive
                || !catallax_boolean::face_contains(
                    model,
                    selected,
                    &point.plus(&direction.to_real().scaled(&time)),
                    dec,
                )?
            {
                return Err(bad(
                    "features.limit-missed",
                    "a profile line misses the selected planar trim",
                ));
            }
        }
    }
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
                let side = sign(&point.minus(&plane.origin).dot(&plane.normal), dec)?;
                match (side, orientation) {
                    (Sign::Zero, _) => return Ok(None),
                    (Sign::Positive, Sense::Same) | (Sign::Negative, Sense::Reversed) => {
                        PatchSelection::Discard
                    }
                    _ => PatchSelection::Keep(Sense::Same),
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
fn circle_cylinder(
    recipe: &Recipe,
    target: &StepId,
    dec: &mut Decider,
) -> Result<Option<Operation>, Refusal> {
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
    Ok(Some(Operation::Cylinder {
        radius: radius.clone(),
        height: distance.value().abs().into(),
        pose,
    }))
}

fn footprint_bounds(profile: &Profile) -> [[Rational; 2]; 2] {
    std::array::from_fn(|axis| {
        let mut range: Option<[Rational; 2]> = None;
        for segment in profile.loops.iter().flat_map(|boundary| &boundary.segments) {
            let [lo, hi] = if let Some((centre, radius_sq, _)) = &segment.circle {
                let radius = Ball::from_rational(radius_sq, 64)
                    .sqrt()
                    .expect("positive profile radius")
                    .hi();
                [&centre[axis] - &radius, &centre[axis] + &radius]
            } else {
                [
                    segment.start[axis].clone().min(segment.end[axis].clone()),
                    segment.start[axis].clone().max(segment.end[axis].clone()),
                ]
            };
            range = Some(match range {
                Some([a, b]) => [a.min(lo), b.max(hi)],
                None => [lo, hi],
            });
        }
        range.expect("nonempty profile")
    })
}

#[allow(clippy::too_many_arguments)]
fn quadric(
    model: &mut Model,
    recipe: &Recipe,
    target: Step,
    sketch: (&StepId, &Profile),
    authored: Option<&catallax_recipe::Vector3>,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(BodyId, Rational), Refusal> {
    let (profile_id, profile) = sketch;
    let target_id = target.id.clone();
    let mut definition = Recipe {
        format: recipe.format.clone(),
        name: recipe.name.clone(),
        description: String::new(),
        result: target_id.clone(),
        steps: vec![
            recipe.step(profile_id).expect("checked profile").clone(),
            target,
        ],
        materials: Vec::new(),
    };
    let cap_error = match curved::AuthoredCap::compile(
        &definition,
        &target_id,
        profile_id,
        authored,
        None,
        dec,
    ) {
        Ok(cap) => return Ok((cap.build(model, step, dec)?, cap.radius().clone())),
        Err(error)
            if may_retry(&error)
                && profile.loops.iter().any(|boundary| {
                    boundary
                        .segments
                        .iter()
                        .any(|segment| segment.circle.is_some())
                }) =>
        {
            error
        }
        Err(error) => return Err(error),
    };
    let Operation::Sketch { loops, .. } = &mut definition.steps[0].op else {
        unreachable!("checked sketch")
    };
    let [[u0, u1], [v0, v1]] = footprint_bounds(profile).map(|[lo, hi]| {
        let margin = (&hi - &lo) / Rational::from_i64(128);
        [&lo - &margin, &hi + &margin]
    });
    *loops = vec![catallax_recipe::Loop(
        [
            [u1.clone(), v0.clone()],
            [u1, v1.clone()],
            [u0.clone(), v1],
            [u0, v0],
        ]
        .into_iter()
        .enumerate()
        .map(|(i, point)| catallax_recipe::Segment::Line {
            id: catallax_recipe::ElementId::new(format!("enclosure-{i}")).expect("fixed side id"),
            to: catallax_recipe::Point2(point.map(Into::into)),
        })
        .collect(),
    )];
    let cap =
        curved::AuthoredCap::compile(&definition, &target_id, profile_id, authored, None, dec)?;
    let bound = cap.radius().clone();
    let envelope = cap.build(model, step, dec)?;
    let direction = direction(profile, authored)?;
    let reach = extent(profile, &bound, &direction);
    let offset = direction.scaled(&reach);
    // Both construction caps are outside the envelope. The result inherits
    // the actual sketch rim, never the four walls of the enclosure.
    let tool = crate::extrude::extrude_vector(
        model,
        &shifted(profile, &offset.scaled(&Rational::from_i64(-1))),
        step,
        &offset.scaled(&Rational::from_i64(2)),
        [None, None],
        None,
        dec,
    )?;
    let body = catallax_boolean::boolean(
        model,
        envelope,
        tool,
        catallax_boolean::BoolOp::Intersect,
        dec,
    )
    .map_err(|error| error.refusal)?
    .body;
    if model.body(body).kind != BodyKind::Solid || !model.check(body).is_empty() {
        return Err(cap_error);
    }
    Ok((body, bound))
}

fn axial_cap(
    model: &Model,
    body: BodyId,
    profile: &Profile,
    direction: &RVector3,
) -> Option<(FaceId, Plane)> {
    let ray = direction.to_real();
    let normal = normal(profile);
    let sketch = profile.pose.apply_point(&Point3::origin());
    model
        .faces_of_body(body)
        .into_iter()
        .filter_map(|face| {
            let Surface::Plane(plane) = model.geom().surface(model.face(face).surface) else {
                return None;
            };
            let plane_normal = plane.normal.as_rational()?;
            if !normal.cross(&plane_normal).is_zero() || !direction.cross(&plane_normal).is_zero() {
                return None;
            }
            let offset = plane
                .origin
                .minus(&sketch)
                .dot(&plane.normal)
                .as_rational()?;
            let speed = ray.dot(&plane.normal).as_rational()?;
            let time = offset / speed;
            (time.sign() == Sign::Positive).then(|| (time, face, plane.clone()))
        })
        .min_by(|a, b| a.0.cmp(&b.0))
        .map(|(_, face, plane)| (face, plane))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn limit(
    model: &mut Model,
    recipe: &Recipe,
    target: (&StepId, BodyId, &Rational),
    sketch: (&StepId, &Profile),
    authored: Option<&catallax_recipe::Vector3>,
    selected: Option<FaceId>,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(BodyId, Rational), Refusal> {
    let (target_id, target_body, bound) = target;
    let (profile_id, profile) = sketch;
    let direction = direction(profile, authored)?;
    let reach = extent(profile, bound, &direction);
    if let Some(face) = selected {
        let surface = model.geom().surface(model.face(face).surface).clone();
        let body = match &surface {
            Surface::Plane(plane) => Some(planar(
                model, face, plane, profile, &direction, &reach, step, dec,
            )?),
            Surface::BSpline(SplineSupport::BezierExtrusion(support)) => Some(bezier::build(
                model, face, support, profile, &direction, step, dec,
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
    let cylinder = match recipe.step(target_id).map(|step| &step.op) {
        Some(op @ Operation::Cylinder { .. }) => op.clone(),
        _ => circle_cylinder(recipe, target_id, dec)?.ok_or_else(|| {
            unsupported("the evaluated limiting target has no supported first-hit proof")
        })?,
    };
    if selected.is_none()
        && let Some((face, plane)) = axial_cap(model, target_body, profile, &direction)
    {
        let body = planar(model, face, &plane, profile, &direction, &reach, step, dec)?;
        return Ok((
            body,
            bound.clone().max(crate::evaluate::profile_radius(profile)),
        ));
    }
    let mut target_step = recipe.step(target_id).expect("checked target").clone();
    target_step.op = cylinder;
    quadric(
        model,
        recipe,
        target_step,
        (profile_id, profile),
        authored,
        step,
        dec,
    )
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
    let face = select(model, target_body, reference, dec)?;
    if distance.is_zero() {
        return Err(bad(
            "extrude.zero-distance",
            "extrusion distance must not be zero",
        ));
    }
    if profile.is_empty() {
        return Ok((catallax_topology::empty(model), Rational::zero()));
    }
    let reach = extent(profile, bound, &direction);
    let offset = direction.scaled(&(-reach));
    let behind = move_recipe(recipe, profile_id, &offset);
    let shifted_profile = shifted(profile, &offset);
    let target_step = behind.step(target_id).expect("checked target");
    let column = match &target_step.op {
        Operation::Torus { .. } => {
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
        Operation::Sphere { .. } | Operation::Cylinder { .. } | Operation::Cone { .. }
            if !matches!(
                model.geom().surface(model.face(face).surface),
                Surface::Plane(_)
            ) =>
        {
            quadric(
                model,
                &behind,
                target_step.clone(),
                (profile_id, &shifted_profile),
                authored,
                step,
                dec,
            )?
            .0
        }
        _ => {
            limit(
                model,
                &behind,
                target,
                (profile_id, &shifted_profile),
                authored,
                Some(face),
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
    let parts = &model.body(column).lumps;
    let mut owners = BTreeMap::new();
    for (index, lump) in parts.iter().enumerate() {
        for shell in &model.lump(*lump).shells {
            owners.extend(model.shell(*shell).faces.iter().map(|face| (*face, index)));
        }
    }
    let mut shells = vec![Vec::new(); parts.len()];
    let mut connectors = BTreeMap::new();
    let reverse = crate::extrude::reverse;
    for (old, moved) in sources.into_iter().zip(copies) {
        if !named(model, old, Role::End) {
            continue;
        }
        let faces = &mut shells[owners[&old]];
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
    if shells.is_empty() || shells.iter().any(Vec::is_empty) {
        return Err(bad(
            "features.limit-missed",
            "the selected start face does not cover the profile",
        ));
    }
    let lumps = shells
        .into_iter()
        .map(|faces| {
            let shell = model.add_shell(Shell { faces });
            model.add_lump(Lump {
                shells: vec![shell],
            })
        })
        .collect();
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps,
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
