//! Feature patterns save ordinary operations with explicit bindings and order.

use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_recipe::{Operation, Recipe, Step, StepId};
use catallax_shapes::{Isometry, Point3, Pose, RVector3, Rotation, Vector3};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};

/// One authored feature placement. Persist the returned ordinary recipe, not this edit.
/// All ids are supplied by the editor so adding instances preserves earlier ids.
#[derive(Clone, Debug)]
pub struct FeatureInstance {
    /// Model-space placement of this instance.
    pub placement: Isometry,
    /// Every source graph step and its unique saved instance id.
    pub steps: BTreeMap<StepId, StepId>,
    /// Every external input and the actual earlier body/profile/path read by this instance.
    pub inputs: BTreeMap<StepId, StepId>,
    /// Every external origin creator and its explicit binding. Identity is explicit too.
    pub origins: BTreeMap<StepId, StepId>,
}

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("pattern.expansion", detail)
}
fn exact() -> Refusal {
    Refusal::new(
        "pattern.expansion-exact-parameter",
        "placement has no exact rational spelling in recipe/3; an exact construction form is required",
    )
}
fn rational(v: &Real) -> Result<catallax_recipe::Number, Refusal> {
    v.as_rational().map(Into::into).ok_or_else(exact)
}
fn vector(v: &catallax_recipe::Vector3) -> Vector3 {
    let [x, y, z] = v.0.each_ref().map(|n| Real::from(n.value().clone()));
    Vector3::new(x, y, z)
}
fn point(v: &catallax_recipe::Point3) -> Point3 {
    let [x, y, z] = v.0.each_ref().map(|n| Real::from(n.value().clone()));
    Point3::new(x, y, z)
}
fn save_vector(v: Vector3) -> Result<catallax_recipe::Vector3, Refusal> {
    Ok(catallax_recipe::Vector3([
        rational(&v.x)?,
        rational(&v.y)?,
        rational(&v.z)?,
    ]))
}
fn save_point(p: Point3) -> Result<catallax_recipe::Point3, Refusal> {
    Ok(catallax_recipe::Point3([
        rational(&p.x)?,
        rational(&p.y)?,
        rational(&p.z)?,
    ]))
}
fn pose(p: &catallax_recipe::Pose) -> Result<Pose, Refusal> {
    let [w, x, y, z] = p.rotation.each_ref().map(|n| n.value().clone());
    let [tx, ty, tz] = p.translation.0.each_ref().map(|n| n.value().clone());
    Ok(Pose::from_rational(
        Rotation::from_quaternion(w, x, y, z)?,
        RVector3::new(tx, ty, tz),
    ))
}

// Homogeneous quaternion recovered from an exact rational proper orthogonal matrix.
// No square root, normalisation, or rounded coefficient is saved.
fn save_pose(
    map: &Isometry,
    source: &catallax_recipe::Pose,
    flip: usize,
) -> Result<catallax_recipe::Pose, Refusal> {
    let composed = Isometry::from(&pose(source)?).then(map);
    let mut m = composed.matrix().clone();
    if map.reversing() {
        for row in &mut m {
            row[flip] = -&row[flip];
        }
    }
    let mut rational_matrix = std::array::from_fn(|_| std::array::from_fn(|_| Rational::zero()));
    for i in 0..3 {
        for j in 0..3 {
            rational_matrix[i][j] = m[i][j].as_rational().ok_or_else(exact)?;
        }
    }
    let m: [[Rational; 3]; 3] = rational_matrix;
    let trace = &(&m[0][0] + &m[1][1]) + &m[2][2] + Rational::one();
    let q = if !trace.is_zero() {
        [
            trace,
            &m[2][1] - &m[1][2],
            &m[0][2] - &m[2][0],
            &m[1][0] - &m[0][1],
        ]
    } else {
        let i = (0..3)
            .find(|i| m[*i][*i] != Rational::from_i64(-1))
            .expect("proper half-turn has an axis");
        let mut q = [
            Rational::zero(),
            Rational::zero(),
            Rational::zero(),
            Rational::zero(),
        ];
        for j in 0..3 {
            q[j + 1] = &m[j][i] + &Rational::from_i64(i64::from(i == j));
        }
        q
    };
    Ok(catallax_recipe::Pose {
        rotation: q.map(Into::into),
        translation: save_vector(composed.translation)?,
    })
}

fn geometry(op: &mut Operation, recipe: &Recipe, map: &Isometry) -> Result<(), Refusal> {
    use Operation::*;
    let p = |p: &mut catallax_recipe::Point3| -> Result<(), Refusal> {
        *p = save_point(map.apply_point(&point(p)))?;
        Ok(())
    };
    let v = |v: &mut catallax_recipe::Vector3| -> Result<(), Refusal> {
        *v = save_vector(map.apply_vector(&vector(v)))?;
        Ok(())
    };
    let axial = |a: &mut catallax_recipe::Vector3| -> Result<(), Refusal> {
        v(a)?;
        if map.reversing() {
            for n in &mut a.0 {
                *n = (-n.value()).into();
            }
        }
        Ok(())
    };
    match op {
        Sketch { plane, .. } => *plane = save_pose(map, plane, 2)?,
        Block { size, pose: frame } => {
            let original = pose(frame)?;
            *frame = save_pose(map, frame, 0)?;
            if map.reversing() {
                let shift = original.apply_vector(
                    &RVector3::new(size[0].value().clone(), Rational::zero(), Rational::zero())
                        .to_real(),
                );
                let shift = map.apply_vector(&shift);
                frame.translation = save_vector(vector(&frame.translation).plus(&shift))?;
            }
        }
        Cylinder { pose: frame, .. }
        | Sphere { pose: frame, .. }
        | Cone { pose: frame, .. }
        | Torus { pose: frame, .. } => *frame = save_pose(map, frame, 0)?,
        Extrude {
            profile, direction, ..
        }
        | ExtrudeThrough {
            profile, direction, ..
        }
        | ExtrudeTo {
            profile, direction, ..
        }
        | ExtrudeNext {
            profile, direction, ..
        } => {
            if let Some(direction) = direction {
                v(direction)?;
            } else if map.reversing() {
                let Some(Step {
                    op: Sketch { plane, .. },
                    ..
                }) = recipe.step(profile)
                else {
                    return Err(bad("extrusion profile is not a sketch"));
                };
                let n = pose(plane)?.apply_vector(
                    &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()).to_real(),
                );
                *direction = Some(save_vector(map.apply_vector(&n))?);
            }
        }
        Revolve { axis, .. } => {
            p(&mut axis.point)?;
            axial(&mut axis.direction)?;
        }
        Hole { at, .. } => {
            for at in at {
                p(&mut at.point)?;
                v(&mut at.direction)?;
            }
        }
        Path { start, segments } => {
            p(start)?;
            for segment in segments {
                match segment {
                    catallax_recipe::PathSegment::Line { to, .. } => p(to)?,
                    catallax_recipe::PathSegment::Arc {
                        centre, to, axis, ..
                    } => {
                        p(centre)?;
                        p(to)?;
                        axial(axis)?;
                    }
                }
            }
        }
        Pattern { layout, .. } => match layout {
            catallax_recipe::PatternLayout::Linear { vector, .. } => v(vector)?,
            catallax_recipe::PatternLayout::Circular { axis, .. } => {
                p(&mut axis.point)?;
                axial(&mut axis.direction)?;
            }
        },
        Mirror { plane, .. } => {
            p(&mut plane.point)?;
            v(&mut plane.normal)?;
        }
        Draft { pull, neutral, .. } => {
            v(pull)?;
            p(&mut neutral.point)?;
            v(&mut neutral.normal)?;
        }
        Boolean { .. }
        | Fillet { .. }
        | Chamfer { .. }
        | Shell { .. }
        | Offset { .. }
        | Thicken { .. }
        | Sweep { .. }
        | Loft { .. } => (),
        _ => {
            return Err(bad(format!(
                "feature expansion of {} is not implemented",
                op.name()
            )));
        }
    }
    Ok(())
}

fn rewrite(
    value: &mut Value,
    bindings: &BTreeMap<StepId, StepId>,
    references: &BTreeMap<StepId, StepId>,
    map: &Isometry,
    reflected_blocks: &BTreeSet<&StepId>,
) -> Result<(), Refusal> {
    match value {
        Value::Array(values) => {
            for value in values {
                rewrite(value, bindings, references, map, reflected_blocks)?;
            }
        }
        Value::Object(values) => {
            let input_object = values.contains_key("op")
                || values.contains_key("profile") && values.contains_key("start");
            if let Some(source) = values.get("step").and_then(Value::as_str)
                && reflected_blocks.iter().any(|id| id.as_str() == source)
                && let Some(segment) = values
                    .get_mut("role")
                    .and_then(|v| v.get_mut("side"))
                    .and_then(|v| v.get_mut("segment"))
            {
                match segment.as_str() {
                    Some("x_min") => *segment = Value::String("x_max".into()),
                    Some("x_max") => *segment = Value::String("x_min".into()),
                    _ => (),
                }
            }
            for (key, value) in values {
                if key == "step"
                    || input_object
                        && matches!(
                            key.as_str(),
                            "body" | "target" | "tool" | "profile" | "path" | "region"
                        )
                        && value.is_string()
                {
                    let source =
                        StepId::new(value.as_str().ok_or_else(|| bad("invalid step binding"))?)
                            .map_err(|e| bad(e.to_string()))?;
                    let table = if key == "step" { references } else { bindings };
                    let target = table.get(&source).ok_or_else(|| {
                        bad(format!(
                            "missing explicit {} binding for {source}",
                            if key == "step" { "origin" } else { "input" }
                        ))
                    })?;
                    *value = Value::String(target.as_str().into());
                } else if key == "near" && value.is_array() {
                    let source = serde_json::from_value::<catallax_recipe::Point3>(value.clone())
                        .map_err(|e| bad(e.to_string()))?;
                    *value = serde_json::to_value(save_point(map.apply_point(&point(&source)))?)
                        .expect("point serialises");
                } else {
                    rewrite(value, bindings, references, map, reflected_blocks)?;
                }
            }
        }
        _ => (),
    }
    Ok(())
}

/// Appends instances of an ordered operation graph and returns the saved recipe.
/// Each external input and reference has an explicit binding. A target-changing
/// instance can read the prior instance's saved result through `inputs`.
/// Failure is atomic and never rounds a parameter to make it serialisable.
pub fn expand_feature(
    recipe: &Recipe,
    graph: &[StepId],
    result: &StepId,
    instances: &[FeatureInstance],
    dec: &mut Decider,
) -> Result<Recipe, Refusal> {
    recipe.check().map_err(|e| bad(e.to_string()))?;
    let ids = graph.iter().collect::<BTreeSet<_>>();
    if ids.len() != graph.len() || !ids.contains(result) {
        return Err(bad("graph repeats an id or omits its result"));
    }
    let source = graph
        .iter()
        .map(|id| {
            recipe
                .step(id)
                .ok_or_else(|| bad(format!("unknown graph step {id}")))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut expanded = recipe.clone();
    for instance in instances {
        let placement = instance.placement.reduced(dec)?;
        let reflected_blocks = source
            .iter()
            .filter(|step| {
                instance.placement.reversing() && matches!(step.op, Operation::Block { .. })
            })
            .map(|step| &step.id)
            .collect();
        if instance.steps.keys().collect::<BTreeSet<_>>() != ids {
            return Err(bad("instance step ids must cover exactly the source graph"));
        }
        let mut bindings = instance.inputs.clone();
        let mut references = instance.origins.clone();
        for (a, b) in &instance.steps {
            if bindings.insert(a.clone(), b.clone()).is_some()
                || references.insert(a.clone(), b.clone()).is_some()
            {
                return Err(bad("internal step has a conflicting external binding"));
            }
        }
        for step in &source {
            let mut op = step.op.clone();
            geometry(&mut op, recipe, &placement)?;
            let mut saved = serde_json::to_value(op).expect("operation serialises");
            rewrite(
                &mut saved,
                &bindings,
                &references,
                &placement,
                &reflected_blocks,
            )?;
            let op = serde_json::from_value(saved).map_err(|e| bad(e.to_string()))?;
            expanded.steps.push(Step {
                id: instance.steps[&step.id].clone(),
                op,
            });
        }
        expanded.result = instance.steps[result].clone();
        expanded.check().map_err(|e| bad(e.to_string()))?;
    }
    Ok(expanded)
}
